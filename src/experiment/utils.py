import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import optimize, stats
from sklearn.linear_model import LinearRegression
from sympy import Expr, lambdify, latex, symbols

VOIGT = {00: 0, 11: 1, 22: 2, 12: 3, 21: 3, 2: 4, 20: 4, 1: 5, 10: 5}

def deltaij(i: int, j: int) -> int:
    return 1 if i == j else 0

def rotation_matrix(theta: float) -> np.ndarray:
    """Rotation matrix for a given angle theta in degrees."""
    theta_rad = np.deg2rad(theta)
    cos_theta, sin_theta = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([
        [cos_theta, -sin_theta, 0],
        [sin_theta, cos_theta, 0],
        [0, 0, 1]
    ])

def miller_to_cartesian(surface_plane: np.ndarray, SAW_plane: np.ndarray) -> np.ndarray:
    """Convert Miller indices to Cartesian transformation matrix."""
    surface_plane = np.array(surface_plane) / np.linalg.norm(surface_plane)
    SAW_plane = np.array(SAW_plane) / np.linalg.norm(SAW_plane)

    m1 = np.cross(surface_plane, SAW_plane)
    m2 = surface_plane
    m3 = SAW_plane
    return np.array([m1, m2, m3])

def orient_Cijkl(C_ijkl: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """Orient a 3x3x3x3 tensor according to a transformation matrix."""
    return np.einsum('ai,bj,ck,dl,ijkl->abcd', matrix, matrix, matrix, matrix, C_ijkl)

def get_Cijkl(C_ij: np.ndarray) -> np.ndarray:
    """Convert a 6x6 matrix into a 3x3x3x3 tensor according to Voigt notation."""
    return np.array([[[[C_ij[VOIGT[10*i+j]][VOIGT[10*k+l]]
                 for i in range(3)] for j in range(3)]
                 for k in range(3)] for l in range(3)])

def set_direction(C_ijkl: np.ndarray, q: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Generate a random wave vector direction if not provided.
    Returns a 3x3 Christoffel matrix.
    """
    if q is None:
        cos_theta = np.random.ranf()
        theta = np.arccos(cos_theta)
        phi = 2.0 * np.pi * np.random.ranf()
        sin_theta = np.sqrt(1 - cos_theta**2)
        q = np.array([np.cos(phi)*sin_theta, np.sin(phi)*sin_theta, cos_theta])
    return np.dot(q, np.dot(q, C_ijkl))

def read_fits(raw_path: str, fit_path: str, interval: Optional[Tuple[float, float]] = None) -> pd.DataFrame:
    """Read the fit results and return a DataFrame with the parsed data."""
    df = pd.read_csv(fit_path)
    dr = []
    for _, run in df.iterrows():
        path = os.path.join(raw_path, run['run_name'])
        data = {}
        xs, ys = [], []
        with open(path, 'r') as file:
            for i, line in enumerate(file):
                if i < 14:
                    key, value = line.strip().split('\t', 1)
                    data[key] = value
                elif i == 15:
                    xn, yn = line.strip().split('\t')
                elif i > 15:
                    x, y = map(float, line.strip().split('\t'))
                    xs.append(x)
                    ys.append(y)
                
        data[xn], data[yn] = xs, ys
        dr.append(data)
    df = pd.concat([df, pd.DataFrame(dr)], axis=1)

    df['Timestamp'] = pd.to_datetime(df['time stamp (ms)'], format='%I:%M:%S %p')
    df['time[s]'] = (df['Timestamp'] - df['Timestamp'].iloc[0]).dt.total_seconds()
    df['time[min]'] = df['time[s]'] / 60

    if interval:
        start, end = interval
        df = df[(df['time[s]'] >= start) & (df['time[s]'] <= end)]
    return df

def read_srim(raw_path: str, total_fluence: float, atomic_density: float) -> pd.DataFrame:
    """Read the SRIM output file and return a DataFrame with the parsed data."""
    with open(raw_path, 'r') as f:
        lines = f.readlines()

    start_index = next(i for i, line in enumerate(lines) if "-----------  -----------  ------------" in line)
    data_lines = lines[start_index + 1:]

    parsed_data = []
    for line in data_lines:
        if line.strip() == "" or set(line.strip()) == {'0', '.', 'E', '+', '-', ' '}:
            break
        match = re.match(r'(\S+)\s+(\S+)\s+(\S+)', line.strip())
        if match:
            parsed_data.append(list(match.groups()))

    df = pd.DataFrame(parsed_data, columns=['target_depth[A]', 'vacancies_ion', 'vacancies_recoil']).astype(float)

    df['depth[um]'] = df['target_depth[A]'] / 10000
    df['vacancies_total'] = df['vacancies_ion'] + df['vacancies_recoil']
    df['dpa'] = df['vacancies_total'] * total_fluence * 1e14 / atomic_density
    return df

def read_temperature(path: str, offset: int, duration: int) -> pd.DataFrame:
    """Read the temperature data and return a DataFrame with the parsed data."""
    df = pd.read_csv(path).iloc[offset:offset + duration]
    if 'Time' in df.columns:
        df['Time'] -= df['Time'].iloc[0]
    return df

def chi_square_filter(df: pd.DataFrame, x_col: str, y_col: str, confidence: float = 0.95) -> pd.DataFrame:
    """Filter out outliers from a DataFrame using standardized residuals."""
    x = df[x_col].values.reshape(-1, 1)
    y = df[y_col].values

    model = LinearRegression().fit(x, y)
    expected_y = model.predict(x)
    
    residuals = y - expected_y
    std_residuals = residuals / np.std(residuals)
    threshold = np.sqrt(stats.chi2.ppf(confidence, df=1))
    return df[np.abs(std_residuals) <= threshold]

def fit_performance_plot(path: str, Q1: float = 20) -> Tuple[float, float, float]:
    """Fit the performance curves of CVi Model CGR409 Cryocooler to a logarithmic function."""
    df = pd.read_csv(f'{path}/{Q1}W.csv', header=None, names=['temperature[K]', 'heat_load[W]'])
    df = df.sort_values(by='temperature[K]')

    mask = (df['temperature[K]'] > 0) & (df['heat_load[W]'] > 0)
    x, y = df['temperature[K]'][mask], df['heat_load[W]'][mask]

    def logarithmic(x, a, b, c):
        return a * np.log(b * x) + c
    
    return optimize.curve_fit(logarithmic, x, y)[0]

def solve_equation(equation: Expr, variable: str, values: Dict[str, float], initial_guesses: List[float]) -> List[float]:
    """Solve a symbolic equation for a variable using fsolve with multiple initial guesses."""
    equation_func = lambdify(symbols(variable), equation.subs(values), 'numpy')
    solutions = [optimize.fsolve(equation_func, guess)[0] for guess in initial_guesses]
    
    return solutions[0] if np.allclose(solutions, solutions[0], atol=1e-2) else solutions
    
def display_equations(equations: Dict[str, Expr]) -> None:
    """Print symbolic equations in LaTeX format."""
    for name, eq in equations.items():
        print(f"{name} = {latex(eq)}")