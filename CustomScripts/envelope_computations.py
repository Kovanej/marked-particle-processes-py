
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

with open("../config_envelope.json", "r") as json_data:
    config_json = json.loads(json_data.read())
envelope_count = config_json["permutation_tests_parameters"]["number_of_permutations"]
config_parser = ConfigParser(config=config_json)

cols = ['Seed', 'Grain Type', 'Model', 'Intensity', 'f-Mark Type', 'Weight Type', 'Input Value', 'PWFCF Value']

envelopes_df = pd.DataFrame(columns=cols)

p_val_low = np.floor(envelope_count * 0.025)
p_val_high = np.ceil(envelope_count * 0.975)

for _ in range(envelope_count):
    seed = init_seed + _
    config_parser.initialize_the_processes(seed=seed)
    result_saver = config_parser.return_the_result_saver(seed=seed)
    envelopes_df = pd.concat([envelopes_df, result_saver.results_all_df])

envelopes_df['Rank'] = envelopes_df.groupby([
    'Grain Type', 'Model', 'Intensity', 'f-Mark Type', 'Weight Type', 'Input Value'])['PWFCF Value'].rank()

grouped_fw_types = envelopes_df.groupby(['Grain Type', 'Model', 'Intensity', 'f-Mark Type', 'Weight Type'])

for g_n, g_df in grouped_fw_types:
    g_df.sort_values(["Grain Type", "Model", "Intensity", "f-Mark Type", "Weight Type", "Input Value", "PWFCF Value"])
    grouped = g_df.groupby(['Seed'])
    for group_name, group_df in grouped:
        plt.plot(group_df['Input Value'], group_df['PWFCF Value'], label=str(group_name), marker=".")
    plt.show()
    plt.close()
envelopes_df.to_csv(f"./envelope_test_vals_{datetime.now().__str__().replace(':', '-')}.csv")
a=1
