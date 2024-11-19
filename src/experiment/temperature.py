from typing import Any, Tuple

import numpy as np
import pandas as pd

from src.experiment.plot import (
    plot_time_vs_temperature,
    plot_time_vs_temperature_zoom,
)
from src.experiment.utils import read_temperature

class Temperature:
    def __init__(self, config: dict[str, Any]) -> None:
        self.ref_path = config['ref_path']
        self.ref_offset = config['ref_offset']
        self.tgs_path = config['tgs_path']
        self.tgs_offset = config['tgs_offset']
        self.ion_path = config['ion_path']
        self.ion_offset = config['ion_offset']
        self.ref_tgs_duration = config['ref_tgs_duration']
        self.ion_duration = config['ion_duration']
        self.zoom_duration = config['zoom_duration']

        self.ref_df, self.tgs_df, self.ion_df = self.process()

    def process(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        ref_df = read_temperature(self.ref_path, self.ref_offset, self.ref_tgs_duration)
        tgs_df = read_temperature(self.tgs_path, self.tgs_offset, self.ref_tgs_duration)
        ion_df = read_temperature(self.ion_path, self.ion_offset, self.ion_duration)

        return ref_df, tgs_df, ion_df

    def plot(self) -> None:
        ref_tgs_time = np.arange(self.ref_tgs_duration) / (60 * 60)
        ion_time = np.arange(self.ion_offset, self.ion_offset + self.ion_duration) / (60 * 60)

        plot_time_vs_temperature(self.ref_df, self.tgs_df, self.ion_df, ref_tgs_time, ion_time)
        plot_time_vs_temperature_zoom(self.ref_df, self.tgs_df, self.ion_df, ref_tgs_time, ion_time, self.zoom_duration)