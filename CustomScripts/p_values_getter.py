
from datetime import datetime
import json
import os
import pandas as pd
import numpy as np
from typing import Dict

from utils.config_parser import ConfigParser


def assign_the_lexicographic_value(g_df, envelope_count):
    # assign the "rooftop" rank
    g_df['Rank Reversed'] = envelope_count + 1 - g_df['Rank']
    g_df['Both Sided Rank'] = g_df[['Rank', 'Rank Reversed']].min(axis=1)
    g_df['Frequency'] = g_df.groupby(['Seed', 'Both Sided Rank'])['Seed'].transform('count')
    g_df_frequency = g_df[['Seed', 'Both Sided Rank', 'Frequency']].drop_duplicates().sort_values(
        ['Seed', 'Both Sided Rank', 'Frequency'])
    dfi = g_df_frequency.set_index(['Seed', 'Both Sided Rank'])[['Frequency']]
    long = dfi.unstack('Both Sided Rank', fill_value=0)
    long.columns = long.columns.droplevel(0)
    cols = list(long.columns)
    df_final = long.sort_values(by=cols, ascending=False)
    return df_final


print(datetime.now())
INIT_SEED = 69
PERMUTATIONS = 4999

model = "ball_max_shared_area_disc"
F_TYPES = ["product", "first_mark"]
W_TYPES = ["shared_area", "intersection"]
#F_TYPES = ["product"]
#W_TYPES = ["shared_area"]

results_per_w_f = {}
for f_type in F_TYPES:
    for w_type in W_TYPES:
        null_model_permutation_base_df = pd.read_csv(
            f"./results_per_model/results_splitted_mod={model}_alpha=0_w={w_type}_f={f_type}.csv")
        null_model_permutation_base_df = null_model_permutation_base_df[
            null_model_permutation_base_df['Seed'] <= INIT_SEED + PERMUTATIONS][
            ['Seed', 'Grain Type', 'Model', 'Intensity', 'f-Mark Type',
             'Weight Type', 'Input Value', 'PWFCF Value']]
        results_per_w_f[(w_type, f_type)] = null_model_permutation_base_df


with open("../config_models.json", "r") as json_data:
    config_json = json.loads(json_data.read())

print("config loaded")

envelope_count = config_json["permutation_tests_parameters"]["number_of_permutations"]
init_seed = config_json["initial_seed"]

seed_inside_envelope = {}
for f_type in F_TYPES:
    for w_type in W_TYPES:
        seed_inside_envelope[(w_type, f_type)] = {}

config_parser = ConfigParser(config=config_json)

df_list_to_save = list()

for _ in range(envelope_count):
    seed = init_seed + _
    # this is ugly, but needed now for check
    if seed in range(69, 5067):
        print("warning")
    config_parser.initialize_the_processes(seed=seed)
    result_saver = config_parser.return_the_result_saver(seed=seed)

    for f_type in F_TYPES:
        for w_type in W_TYPES:

            null_model_permutation_base_df = results_per_w_f[(w_type, f_type)]
            null_model_copy = null_model_permutation_base_df.copy()
            df_to_concat = result_saver.results_all_df
            df_to_concat = df_to_concat.loc[(df_to_concat['f-Mark Type'] == f_type) & (df_to_concat['Weight Type'] == w_type)]
            df_to_evaluate = pd.concat([null_model_copy, df_to_concat])

            df_to_evaluate['Rank'] = df_to_evaluate.groupby(['Input Value'])['PWFCF Value'].rank()
            df_final = assign_the_lexicographic_value(g_df=df_to_evaluate, envelope_count=PERMUTATIONS + 1)
            test_rank = np.where(df_final.index == seed)[0][0]
            seed_inside_envelope[(w_type, f_type)][seed] = False if test_rank < np.ceil(0.05 * (PERMUTATIONS + 1)) else True
            if seed % 100 == 0:
                print(1 - np.mean(list(seed_inside_envelope[(w_type, f_type)].values())))
    print(f"Computations done for seed: {seed}")


for w_type in W_TYPES:
    for f_type in F_TYPES:
        df_list_to_save.append(
            pd.DataFrame({
                'Model': [model], 'Alpha': [config_json['marking_parameters']['alphas'][0]],
                'Number of Permutations': envelope_count,
                'f-Mark Type': [f_type], 'Weight Type': [w_type],
                'Rejection Rate': [1 - np.mean(list(seed_inside_envelope[(w_type, f_type)].values()))]
            }))

df_to_save = pd.concat(df_list_to_save)
df_to_save.to_csv(f"./rejection_rate{config_json['marking_parameters']['alphas'][0]}_{datetime.now().__str__().replace(':', '-')}.csv", index=False)
print(datetime.now())
