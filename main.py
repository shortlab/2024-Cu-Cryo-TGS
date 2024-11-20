import yaml

from src.core.fit import TGSAnalyzer
from src.experiment.material import Material
from src.experiment.irradiation import Irradiation
from src.experiment.study import HeatLoadStudy
from src.experiment.temperature import Temperature

if __name__ == '__main__':

    # with open('config.yaml', "r") as file: config = yaml.safe_load(file)

    # analyzer = TGSAnalyzer(config)
    # analyzer.fit()
    
    with open('experiment.yaml', "r") as file: experiment_config = yaml.safe_load(file)

    material = Material(experiment_config['material'])
    irradiation = Irradiation(experiment_config['tgs']['irradiation'], material)
    irradiation.plot()

    # heat_load_study = HeatLoadStudy(experiment_config['heat_load'])
    # heat_load_study.estimate()

    # temperature = Temperature(experiment_config['temperature'])
    # temperature.plot()
