
from datetime import datetime
from typing import Dict

from executions import execute_from_config
from utils.config_parser import ConfigParser
from tests.wider_window_simulation_tests import simulate_the_processes


def run(config: Dict):
    config_parser = ConfigParser(config=config)
    result_savers = execute_from_config(config_parser=config_parser)
    return result_savers
