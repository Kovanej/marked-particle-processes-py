import json

from executions import first_blood, overnight_computations
from run import run
from utils.results_saver import ResultSaver


if __name__ == '__main__':
    with open("config.json", "r") as json_data:
        config_json = json.loads(json_data.read())
    result_savers = run(config=config_json)
    # result_savers = overnight_computations()
    brkpnt = "breakpoint is here"
