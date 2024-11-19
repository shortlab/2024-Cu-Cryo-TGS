import math
import os
from typing import Any, Tuple

import pandas as pd

from src.experiment.material import Material
from src.experiment.plot import (
    plot_depth_vs_dpa,
    plot_fluence_vs_alpha,
    plot_fluence_vs_saw,
)
from src.experiment.utils import read_fits, read_srim, chi_square_filter

class Irradiation:
    def __init__(self, config: dict[str, Any], material: Material) -> None:
        self.path = config['path']
        self.srim_path = config['srim_path']

        self.probe_depth = config['grating_spacing'] / math.pi
        self.start_current = config['start_current']
        self.end_current = config['end_current']
        self.duration = config['duration']
        self.average_current = (self.start_current - self.end_current) / 2
        self.total_incident_charge = self.average_current * self.duration
        self.ion_charge = config['ion_charge']
        self.total_ions = self.total_incident_charge / (self.ion_charge * 1.602176634e-19)
        self.beam_aperture = config['beam_aperture']
        self.beam_area = math.pi * (self.beam_aperture * 0.001 / 2) * (self.beam_aperture * 0.001 / math.sqrt(2))
        self.total_fluence = self.total_ions / self.beam_area / 10000
        self.m = (self.end_current - self.start_current) / self.duration
        self.b = self.start_current
        
        self.material = material
        self.atomic_density = material.atomic_density

        self.tgs_df, self.srim_df = self.process()

    def process(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        tgs_df = read_fits(self.path, os.path.join(self.path, 'fit', 'fit.csv'), interval=(0, self.duration))
        srim_df = read_srim(self.srim_path, self.total_fluence, self.material.atomic_density)

        tgs_df['fluence[ionsm^-2]'] = tgs_df['time[s]'].apply(self.Phi)
        closest_depth = srim_df['depth[um]'].iloc[(srim_df['depth[um]'] - self.probe_depth).abs().argsort()[:1]].values[0]
        K_tgs = srim_df[srim_df['depth[um]'] == closest_depth]['vacancies_total'].values[0] * 1e10
        tgs_df['dpa'] = tgs_df['time[s]'].apply(lambda t: self.dpa(K_tgs, t))

        tgs_df = chi_square_filter(tgs_df, 'time[s]', 'f[Hz]', confidence=0.95)

        tgs_df.to_csv(os.path.join(self.path, 'process', 'irradiation.csv'), index=False)
        srim_df.to_csv(os.path.join(os.path.dirname(self.srim_path), 'srim.csv'), index=False)
        return tgs_df, srim_df
    
    def plot(self) -> None:
        plot_fluence_vs_alpha(self.tgs_df)
        plot_fluence_vs_saw(self.tgs_df)
        plot_depth_vs_dpa(self.srim_df)

    def I(self, t: float) -> float:
        """Calculate the beam current at time t."""
        return self.m * t + self.b

    def Q(self, t: float) -> float:
        """Calculate the integral of the beam current up to time t."""
        return self.m * t ** 2 / 2 + self.b * t

    def Phi(self, t: float) -> float:
        """Calculate the fluence at time t."""
        return self.Q(t) / (self.ion_charge * 1.602176634e-19 * self.beam_area)

    def dpa(self, K: float, t: float) -> float:
        """Calculate the dose at time t given the vacancies at a certain depth."""
        return K * self.Phi(t) / self.atomic_density