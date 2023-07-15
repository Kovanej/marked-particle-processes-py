
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

df_lower = df_lower.loc[
    (df_lower['f-Mark Type'] == f_type) & (df_lower['Weight Type'] == w_type) & (df_lower['Model'] == model)
][['Model', 'f-Mark Type', 'Weight Type', 'Input Value', 'PWFCF Value']]
df_lower.rename(columns={'PWFCF Value': 'PWFCF Lower'}, inplace=True)

df_upper = df_upper.loc[
    (df_upper['f-Mark Type'] == f_type) & (df_upper['Weight Type'] == w_type) & (df_upper['Model'] == model)
][['Model', 'f-Mark Type', 'Weight Type', 'Input Value', 'PWFCF Value']]
df_upper.rename(columns={'PWFCF Value': 'PWFCF Upper'}, inplace=True)

df_model_envelope_values = merged_df = pd.merge(
    df_lower, df_upper, on=['Model', 'f-Mark Type', 'Weight Type', 'Input Value']
)

df_model_envelope_values.to_csv(f"./completed_envelopes/df_envelope_values_{model}_{w_type}_{f_type}.csv", index=False)
breakpoint_var=1
