
from datetime import datetime
import json
import logging
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Dict

from utils.config_parser import ConfigParser
import utils.const as const
from utils.results_saver import ResultSaver


# set up parameters which dataset to create
model = 'ball_continuous_radius_alpha=0'
f_type = 'first_mark'
w_type = 'shared_area'


with open("../config_models.json", "r") as json_data:
    config_json = json.loads(json_data.read())
envelope_count = config_json["permutation_tests_parameters"]["number_of_permutations"]
init_seed = config_json["initial_seed"]

completed_envelopes = pd.read_csv(f"./completed_envelopes/df_envelope_values_{model}_{w_type}_{f_type}.csv", index_col = False)

seed_inside_envelope: Dict[str, bool] = {}

config_parser = ConfigParser(config=config_json)

for _ in range(envelope_count):
    seed = init_seed + _
    config_parser.initialize_the_processes(seed=seed)
    result_saver = config_parser.return_the_result_saver(seed=seed)
    df_values = result_saver.results_all_df[['Input Value', 'PWFCF Value']]
    df_to_evaluate = pd.merge(completed_envelopes, df_values, on=['Input Value'])
    df_to_evaluate = df_to_evaluate.loc[df_to_evaluate['Input Value'] > 0.04]
    df_to_evaluate.plot(x='Input Value', y=['PWFCF Value', 'PWFCF Lower', 'PWFCF Upper'])
    this_shoud_be_zero = (df_to_evaluate['PWFCF Value'] > df_to_evaluate['PWFCF Upper']).sum() + (
                df_to_evaluate['PWFCF Value'] < df_to_evaluate['PWFCF Lower']).sum()
    seed_inside_envelope[seed] = True if this_shoud_be_zero == 0 else False

print(np.mean(list(seed_inside_envelope.values())))
breakpoint_val=1
