import json

from run import run
from utils.results_saver import ResultSaver


if __name__ == '__main__':
    with open("config.json", "r") as json_data:
        config_json = json.loads(json_data.read())
    result_savers = run(config=config_json)
