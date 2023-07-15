
import os
import pandas as pd
import numpy as np


csv_strs = [f"both_sided_ranks/{file}" for file in os.listdir() if "g_df_frequency_mod=" in file]

for csv_str in csv_strs:
    df = pd.read_csv(csv_str)
    dfi = df.set_index(['Seed', 'Both Sided Rank'])[['Frequency']]
    long = dfi.unstack('Both Sided Rank', fill_value=0)
    long.columns = long.columns.droplevel(0)
    cols = list(long.columns)
    df_final = long.sort_values(by=cols, ascending=False)

    df_seeds = df_final.head(250)
    df_seeds.to_csv(f"./pvals/new_{csv_str.split('.')[0]}_seeds_pvals.csv")

breakpoint_var=1