
from datetime import datetime
from typing import Dict

from utils.config_parser import ConfigParser
from tests.wider_window_simulation_tests import simulate_the_processes


def run(config: Dict):
    config_parser = ConfigParser(config=config)
    result_savers = []
    for seed in range(200, 220):
        print(f"Seed no {seed}: {datetime.now()}")
        result_saver = simulate_the_processes(alphas_list=[0, 0.25, 0.5, 0.75, 1], seed=seed)
        result_savers.append(result_saver)
    return result_savers
