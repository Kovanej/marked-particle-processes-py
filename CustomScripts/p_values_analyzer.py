
import os
import pandas as pd
import numpy as np

path_to_results = "./results/csvs/balls_69-7368"
f_types = ["product",
           #"first_mark", "square"
           ]
w_types = ["intersection",
           #"shared_area"
           ]
models = [
    #"ball_bivariate_radius_alpha=0",
    "ball_continuous_radius_alpha=0",
    #"ball_intersection_count_alpha=0",
    #"ball_max_shared_area_cont_alpha=0"
]

for model in models:
    for f_type in f_types:
        for w_type in w_types:
            df_pvals = pd.read_csv(f"./pvals/new_g_df_frequency_mod={model}_w={w_type}_f={f_type}_int=50_seeds_pvals.csv")

            p_val_clos = df_pvals[240:250]
            seeds_to_analyze = list(p_val_clos['Seed'])

            df = pd.read_csv(f"./both_sided_ranks/g_df_frequency_mod={model}_w={w_type}_f={f_type}_int=50.csv")
            df_ranks = pd.read_csv(
                    f"./ranks/g_df_pure_ranks_mod={model}_w={w_type}_f={f_type}_int=50.csv")

            for seed in reversed(seeds_to_analyze):
                df_seed = df[df['Seed'] == seed]
                results_csv_path = [f"./results/csvs/balls_69-7368/{file}" for file in os.listdir(path_to_results) if f"seed={seed}" in file][0]
                df_results = pd.read_csv(results_csv_path)

                a=1

breakpoint_var=1