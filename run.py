
from datetime import datetime
from typing import Dict

from executions import execute_from_config
from utils.config_parser import ConfigParser
import utils.const as const
from tests.wider_window_simulation_tests import simulate_the_processes


def run(config: Dict):
    config_parser = ConfigParser(config=config)
    result_savers, result_savers_whole_df = execute_from_config(config_parser=config_parser)
    if const.SAVE_RESULTS_TO_CSV:
        result_savers_whole_df.to_csv(f"result_savers_whole_df_{str(datetime.now()).replace(':', '_')}.csv")
    return result_savers
