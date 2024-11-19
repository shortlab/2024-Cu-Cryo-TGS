from typing import Any, Dict, List, Tuple

from sympy import Expr, Symbol, log, symbols

from src.experiment.utils import (
    fit_performance_plot,
    solve_equation,
    display_equations
)

class HeatLoadStudy:
    def __init__(self, config: Dict[str, Any]) -> None:
        self.cryohead_path = config['cryohead_path']
        self.symbols = self.define_symbols()
        self.equations = self.define_equations(self.symbols)
        self.values = config['values']

    def define_symbols(self) -> Tuple[Symbol, ...]:
        """Define all the necessary symbols."""
        return symbols(
            'Q_conduction Q_convection Q_radiation Q_TGS Q_ion '
            'k_steel k_ceramic '
            'L^stage_steel A^stage_steel L^stage_ceramic A^stage_ceramic '
            'L^gimbal_steel A^gimbal_steel L^axis_ceramic A^axis_ceramic '
            'T_ambient T_sample A_chamber sigma epsilon_steel P_TGS E_ion I_ion a b c'
        )
    
    def define_equations(self, symbols: Tuple[Symbol, ...]) -> Dict[str, Expr]:
        """Define heat load component equations."""
        (Q_conduction, Q_convection, Q_radiation, Q_tgs, Q_ion, k_steel, k_ceramic, 
        L_stage_steel, A_stage_steel, L_stage_ceramic, A_stage_ceramic, 
        L_gimbal_steel, A_gimbal_steel, L_axis_ceramic, A_axis_ceramic, 
        T_ambient, T_sample, A_chamber, sigma, epsilon_steel, P_tgs, E_ion, I_ion, a, b, c) = symbols

        R_total = (L_stage_steel / (k_steel * A_stage_steel) + 
                   L_stage_ceramic / (k_ceramic * A_stage_ceramic) + 
                   (2 * (k_steel * A_gimbal_steel / L_gimbal_steel + 
                         k_ceramic * A_axis_ceramic / L_axis_ceramic)) ** -1)
        Q_conduction = (T_ambient - T_sample) / R_total
        Q_radiation = A_chamber * sigma * epsilon_steel * (T_ambient**4 - T_sample**4)
        # Q_convection = 0
        Q_tgs = P_tgs
        Q_ion = E_ion * I_ion
        Q_sample = Q_conduction + Q_radiation + Q_tgs + Q_ion
        Q2_cryohead = a * log(b * T_sample) + c
        equation = Q_sample - Q2_cryohead
        
        return {
            'Q_conduction': Q_conduction,
            'Q_radiation': Q_radiation,
            'Q_tgs': Q_tgs,
            'Q_ion': Q_ion,
            'Q_sample': Q_sample,
            'Q2_cryohead': Q2_cryohead,
            'equation': equation
        }

    def estimate(self, Q1: float = 10, initial_guesses: List[float] = [10, 20, 30, 40, 50]) -> Tuple[float, Dict[str, Expr]]:
        """Estimate the sample temperature for given parameters."""
        popt_log = fit_performance_plot(self.cryohead_path, Q1)
        self.values['a'], self.values['b'], self.values['c'] = popt_log
        evaluated_equations = {name: eq.subs(self.values) for name, eq in self.equations.items()}

        sample_temperature = solve_equation(evaluated_equations['equation'], 'T_sample', self.values, initial_guesses)
        self.values['T_sample'] = sample_temperature
        heat_loads = {name: eq.subs(self.values) for name, eq in self.equations.items() 
                      if name in ['Q_conduction', 'Q_radiation', 'Q_tgs', 'Q_ion', 'Q_sample']}

        print(f"Sample temperature: {sample_temperature} K")
        display_equations(heat_loads)