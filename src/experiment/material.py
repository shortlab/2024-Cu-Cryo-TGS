from typing import Any, Dict

import numpy as np

from src.experiment.utils import (
    get_Cijkl,
    miller_to_cartesian,
    orient_Cijkl,
    set_direction,
)

class Material:
    def __init__(self, config: Dict[str, Any]) -> None:
        self.name = config['name']
        self.surface_plane = config['surface_plane']
        self.atomic_mass = config['atomic_mass']
        self.mass_density = config['mass_density']
        self.atomic_density = (self.mass_density * 6.023e23 * 1e3) / self.atomic_mass
        self.yield_strength = config['yield_strength']
        self.elastic_modulus = config['elastic_modulus']
        self.poisson_ratio = config['poisson_ratio']
        self.C11 = config['C11']
        self.C12 = config['C12']
        self.C44 = config['C44']

    def thermal_conductivity(self, T: float) -> float:
        """
        Thermal conductivity of Cu in W m^-1 K^-1 as a function of temperature.
        Constants sourced from NIST Cryogenic Material Properties data for OHFC Copper RRR=50
        https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm
        """
        a, b, c, d, e, f, g, h, i = 1.8743, -0.41538, -0.6018, 0.13294, 0.26426, -0.0219, -0.051276, 0.0014871, 0.003723
        numerator = a + c * T**0.5 + e * T + g * T**1.5 + i * T**2
        denominator = 1 + b * T**0.5 + d * T + f * T**1.5 + h * T**2
        log_k = numerator / denominator
        return 10**log_k  # W m^-1 K^-1
    
    def specific_heat(self, T: float) -> float:
        """
        Specific heat of Cu in J kg^-1 K^-1 as a function of temperature.
        Constants sourced from NIST Cryogenic Material Properties data for OHFC Copper RRR=50
        https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm
        """
        a, b, c, d, e, f, g, h, i = -1.91844, -0.15973, 8.61013, -18.996, 21.9661, -12.7328, 3.54322, -0.3797, 0
        log_T = np.log10(T)
        log_y = a + b*log_T + c*log_T**2 + d*log_T**3 + e*log_T**4 + f*log_T**5 + g*log_T**6 + h*log_T**7 + i*log_T**8
        return 10**log_y  # J kg^-1 K^-1
    
    def thermal_diffusivity(self, T: float) -> float:
        """Thermal diffusivity of Cu in m^2 s^-1 as a function of temperature."""
        k = self.thermal_conductivity(T)
        Cp = self.specific_heat(T)
        return k / (self.mass_density * Cp)  # m^2 s^-1

    def get_Cij(self, T: float) -> np.ndarray:
        """
        Stiffness constants in GPa of Cu as a function of temperature.
        Constants sourced from H. M. Ledbetter, E. R. Naimon; Elastic Properties of Metals and Alloys. II. Copper. 
        J. Phys. Chem. Ref. Data 1 October 1974; 3 (4): 897–935. https://doi.org/10.1063/1.3253150
        """
        k_C11, k_C12, k_C44 = -2.1e-4, -1.2e-4, -3.5e-4  # K^-1
        C11 = self.C11 * np.exp(k_C11 * (T - 298))
        C12 = self.C12 * np.exp(k_C12 * (T - 298))
        C44 = self.C44 * np.exp(k_C44 * (T - 298))

        return np.array([
            [C11, C12, C12, 0, 0, 0],
            [C12, C11, C12, 0, 0, 0],
            [C12, C12, C11, 0, 0, 0],
            [0, 0, 0, C44, 0, 0],
            [0, 0, 0, 0, C44, 0],
            [0, 0, 0, 0, 0, C44]
        ])
    
    def get_phase_velocity(self, T: float, q: np.ndarray, SAW_plane: np.ndarray) -> np.ndarray:
        """
        Determine eigenvalues, eigenvectors of the Christoffel matrix,
        sort from low to high, then store eigens and phase velocities.
        """
        C_ij = self.get_Cij(T)
        C_ijkl = get_Cijkl(C_ij)
        matrix = miller_to_cartesian(SAW_plane)
        deg = 0
        MM = np.array([
            [np.cos(np.deg2rad(deg)), np.cos(np.deg2rad(90 - deg)), 0],
            [np.cos(np.deg2rad(90 + deg)), np.cos(np.deg2rad(deg)), 0],
            [0, 0, 1]
        ])
        matrix = matrix @ MM.T
        C_ijkl = orient_Cijkl(C_ijkl, matrix) * 1000.0 / self.mass_density
        
        christoffel = set_direction(C_ijkl, q)
        eig_val, eig_vec = np.linalg.eigh(christoffel)
        args = np.argsort(eig_val)
        eig_val = eig_val[args]
        eig_vec = eig_vec.T[args]
        return np.sign(eig_val) * np.sqrt(np.abs(eig_val))
    
    def SAW_speed(self, T: float, q: np.ndarray, SAW_plane: np.ndarray) -> float:
        """SAW speed in m s^-1 as a function of temperature, wave vector, and SAW plane."""
        v = self.get_phase_velocity(T, q, SAW_plane)
        return min(v[0], v[1])  # Return the minimum of v_T1 and v_T2