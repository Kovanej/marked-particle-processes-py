
import os
import numpy as np
import pandas as pd


# path to result dfs - ball
path_to_result_dfs = "./results/csvs/balls_69-7368"

# set up parameters which dataset to create
model = 'ball_continuous_radius_alpha=0'
f_type = 'product'
w_type = 'intersection'

# according to which we have lower and upper seed
seed_lower = 4145
seed_upper = 2313

# load the full dataset
df_lower_files = [filename for filename in os.listdir(f"{path_to_result_dfs}") if f"seed={seed_lower}" in filename]
df_upper_files = [filename for filename in os.listdir(f"{path_to_result_dfs}") if f"seed={seed_upper}" in filename]

if len(df_lower_files) > 1:
    print(f"df_lower_files has more then one file {df_lower_files}")
if len(df_upper_files) > 1:
    print(f"df_lower_files has more then one file {df_upper_files}")

df_lower_file_path = f"{path_to_result_dfs}/{df_lower_files[0]}"
df_upper_file_path = f"{path_to_result_dfs}/{df_upper_files[0]}"

print(df_lower_file_path, df_upper_file_path)

df_lower = pd.read_csv(df_lower_file_path)
df_upper = pd.read_csv(df_upper_file_path)

print(df_lower.head())


