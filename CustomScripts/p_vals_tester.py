
from datetime import datetime
import json
import logging
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from utils.config_parser import ConfigParser
import utils.const as const
from utils.results_saver import ResultSaver


init_seed = 69

with open("../config_p_vals_tester.json", "r") as json_data:
    config_json = json.loads(json_data.read())
envelope_count = config_json["permutation_tests_parameters"]["number_of_permutations"]
config_parser = ConfigParser(config=config_json)

for _ in range(envelope_count):
    seed = init_seed + _
    config_parser.initialize_the_processes(seed=seed)
    result_saver = config_parser.return_the_result_saver(seed=seed)
    breakpoint_var=1
