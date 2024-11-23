from typing import Any, Dict, List, Tuple

from sympy import Expr, Symbol, log, symbols

from src.experiment.utils import (
    fit_performance_plot,
    solve_equation,
)

class HeatLoadStudy:
    def __init__(self, config: Dict[str, Any]) -> None:
        self.cryohead_path = config['cryohead_path']
        self.symbols = self.define_symbols()
        self.equations = self.define_equations(self.symbols)
        self.values = config['values']

    def define_symbols(self) -> Tuple[Symbol, ...]:
        return symbols(
            'Q_cd Q_cv Q_rad Q_tgs Q_ion '
            'k_steel k_ceramic '
            'L^stage_steel A^stage_steel L^stage_ceramic A^stage_ceramic '
            'L^gimbal_steel A^gimbal_steel L^axis_ceramic A^axis_ceramic '
            'T_A T_S S sigma epsilon P_tgs E_ion I_ion a b c C T_F'
        )
    
    def define_equations(self, symbols: Tuple[Symbol, ...]) -> Dict[str, Expr]:
        (Q_cd, Q_cv, Q_rad, Q_tgs, Q_ion, k_steel, k_ceramic, 
        L_stage_steel, A_stage_steel, L_stage_ceramic, A_stage_ceramic, 
        L_gimbal_steel, A_gimbal_steel, L_axis_ceramic, A_axis_ceramic, 
        T_A, T_S, S, sigma, epsilon, P_tgs, E_ion, I_ion, a, b, c, C, T_F) = symbols

        R_T = (L_stage_steel / (k_steel * A_stage_steel) + 
                   L_stage_ceramic / (k_ceramic * A_stage_ceramic) + 
                   (2 * (k_steel * A_gimbal_steel / L_gimbal_steel + 
                         k_ceramic * A_axis_ceramic / L_axis_ceramic)) ** -1)
        Q_cd = (T_A - T_S) / R_T
        Q_rad = S * sigma * epsilon * (T_A**4 - T_S**4)
        # Q_cv = 0
        Q_tgs = P_tgs
        Q_ion = E_ion * I_ion
        Q_H = Q_cd + Q_rad + Q_tgs + Q_ion
        Q2_C = a * log(b * T_S) + c
        T_F = T_S + Q_H / C
        Q_Σ = Q_H - Q2_C
        
        return {
            'Q_cd': Q_cd,
            'Q_rad': Q_rad,
            'Q_tgs': Q_tgs,
            'Q_ion': Q_ion,
            'Q_H': Q_H,
            'Q2_C': Q2_C,
            'T_F': T_F,
            'Q_Σ': Q_Σ  
        }

    def estimate(self, Q1: float = 10, initial_guesses: List[float] = [10, 20, 30, 40, 50]) -> Tuple[float, Dict[str, Expr]]:
        popt_log = fit_performance_plot(self.cryohead_path, Q1)
        self.values['a'], self.values['b'], self.values['c'] = popt_log
        evaluated_equations = {name: eq.subs(self.values) for name, eq in self.equations.items()}

        sample_temperature = solve_equation(evaluated_equations['Q_Σ'], 'T_S', self.values, initial_guesses)
        self.values['T_S'] = sample_temperature
        heat_loads = {name: eq.subs(self.values) for name, eq in self.equations.items() if name in ['Q_cd', 'Q_rad', 'Q_tgs', 'Q_ion', 'Q_H']}
        self.values['T_F'] = float(self.equations['T_F'].subs(self.values))
        final_temperature = self.values['T_F']
        print(f'Final temperature: {final_temperature} K')
        return final_temperature
