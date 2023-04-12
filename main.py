import logging
import json
import pandas as pd

from run import run
from utils.results_saver import ResultSaver


if __name__ == '__main__':
    logging.basicConfig(filename='main_logger.log', filemode='w', level=logging.INFO)
    with open("config.json", "r") as json_data:
        config_json = json.loads(json_data.read())
    result_savers, result_savers_whole_df = run(config=config_json)
    for result_saver in result_savers:
        pd.set_option('display.max_columns', 10)
        print(result_saver.results_all_df)

