import os
from typing import Any

import pandas as pd

from src.experiment.plot import plot_time_vs_alpha_saw_temperature
from src.experiment.utils import read_fits, chi_square_filter

class Cooldown:
    def __init__(self, config: dict[str, Any]) -> None:
        self.path = config['path']
        self.temperature_path = config['temperature_path']
        
        self.df = self.process()

    def process(self) -> pd.DataFrame:
        tgs_df = read_fits(self.path, os.path.join(self.path, 'fit', 'fit.csv'))
        temperature_df = pd.read_csv(self.temperature_path)

        temperature_df['Timestamp'] = pd.to_datetime(temperature_df['Timestamp']).dt.round('s')
        date = temperature_df['Timestamp'].iloc[0].date()
        tgs_df['Timestamp'] = pd.to_datetime(tgs_df['Timestamp']).apply(lambda x: x.replace(year=date.year, month=date.month, day=date.day))
        start, end = tgs_df['Timestamp'].iloc[0], tgs_df['Timestamp'].iloc[-1]
        temperature_df = temperature_df[(temperature_df['Timestamp'] >= start) & (temperature_df['Timestamp'] <= end)]

        tgs_df = chi_square_filter(tgs_df, 'time[s]', 'f[Hz]', confidence=0.95)
        
        df = pd.merge(tgs_df, temperature_df, on='Timestamp', how='inner')
        df.to_csv(os.path.join(self.path, 'process', 'cooldown.csv'), index=False)
        return df
    
    def plot(self) -> None:
        plot_time_vs_alpha_saw_temperature(self.df)