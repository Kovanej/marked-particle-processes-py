
import json
import os
import pandas as pd
import numpy as np
from typing import Dict

from utils.config_parser import ConfigParser


def sort_lexicographically(df: pd.DataFrame, seed_to_be_tested: int):
    # df = pd.read_csv(csv_str)
    dfi = df.set_index(['Seed', 'Both Sided Rank'])[['Frequency']]
    long = dfi.unstack('Both Sided Rank', fill_value=0)
    long.columns = long.columns.droplevel(0)
    cols = list(long.columns)
    df_final = long.sort_values(by=cols, ascending=False)
    df_final[df_final['Seed'] == seed_to_be_tested]
    df_seeds = df_final.head(250)
    # df_seeds.to_csv(f"./pvals/new_{csv_str.split('.')[0]}_seeds_pvals.csv")


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
    long.sort_values(by=cols, ascending=False)
    long.to_csv(f"./g_df_frequency.csv")


model = "ball_bivariate_radius"
f_type = "product"
w_type = "intersection"

null_model_permutation_base_df = pd.read_csv(
    f"./results_per_model/results_splitted_mod={model}_alpha=0_w={w_type}_f={f_type}.csv")

with open("../config_models.json", "r") as json_data:
    config_json = json.loads(json_data.read())

envelope_count = config_json["permutation_tests_parameters"]["number_of_permutations"]
init_seed = config_json["initial_seed"]


seed_inside_envelope: Dict[str, bool] = {}
config_parser = ConfigParser(config=config_json)

for _ in range(envelope_count):
    seed = init_seed + _
    # this is ugly, but needed now
    if seed in range(69, 5068):
        print("warning")
    config_parser.initialize_the_processes(seed=seed)
    result_saver = config_parser.return_the_result_saver(seed=seed)
    


a=1
