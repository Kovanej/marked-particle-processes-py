
from datetime import datetime
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from typing import List


class EnvelopeProcessor(object):

    def __init__(self, paths: List[str] = []):
        self.paths = paths
        self.number_of_paths = len(self.paths)
        self.csvs = {path: [file for file in os.listdir(f"{path}") if ".csv" in file] for path in self.paths}
        self.cols = ['Seed', 'Grain Type', 'Model', 'Intensity', 'f-Mark Type', 'Weight Type', 'Input Value', 'PWFCF Value']
        self.dfs = {path: pd.DataFrame(columns=self.cols) for path in self.paths}

    def load_the_data(self, path: str = "./", all_csv: bool = True):
        for path in self.dfs.keys():
            dfs_list = []
            for file in self.csvs[path]:
                df = pd.read_csv(f"{path}/{file}")
                dfs_list.append(df)
            self.dfs[path] = pd.concat(dfs_list)

    def load_the_whole_dataset(self):
        df = pd.read_csv("./results/csvs/balls_69-7368/df_whole_dataset_balls_new_dropped.csv", index_col=False)
        self.dfs["ball"] = df


envelope_processor = EnvelopeProcessor([
    "./results/csvs/segment_69-5112",
# "./results/csvs/balls_69-7368"
])

envelope_processor.load_the_data()

envelopes_df = envelope_processor.dfs["./results/csvs/segment_69-5112"]

# filter seeds so that we have 5000 permutations
envelopes_df = envelopes_df[envelopes_df.Seed < 5068]
# get rid of distance, it doesn't make sense
envelopes_df = envelopes_df[envelopes_df['Weight Type'] != 'distance']
envelope_count = np.unique(envelopes_df["Seed"]).shape[0]

envelopes_df['Rank'] = envelopes_df.groupby([
    'Grain Type', 'Model', 'Intensity', 'f-Mark Type', 'Weight Type', 'Input Value'])['PWFCF Value'].rank()
grouped_fw_types = envelopes_df.groupby(['Grain Type', 'Model', 'Intensity', 'f-Mark Type', 'Weight Type'])


def assign_the_lexicographic_value(g_df, envelope_count):
    # assign the "rooftop" rank
    g_df['Rank Reversed'] = envelope_count + 1 - g_df['Rank']
    g_df['Both Sided Rank'] = g_df[['Rank', 'Rank Reversed']].min(axis=1)
    g_df['Frequency'] = g_df.groupby(['Seed', 'Both Sided Rank'])['Seed'].transform('count')
    g_df_frequency = g_df[['Seed', 'Both Sided Rank', 'Frequency']].sort_values(
        ['Seed', 'Both Sided Rank', 'Frequency']).drop_duplicates()
    # dfi = g_df_frequency.set_index(['Seed', 'Both Sided Rank'])[['Frequency']]
    # long = dfi.unstack('Both Sided Rank', fill_value=0)
    # long.columns = long.columns.droplevel(0)
    # cols = list(long.columns)
    # long.sort_values(by=cols, ascending=False)
    model_str = np.unique(g_df['Model'])[0]
    intensity_str = np.unique(g_df['Intensity'])[0]
    f_str = np.unique(g_df['f-Mark Type'])[0]
    w_str = np.unique(g_df['Weight Type'])[0]
    g_pure_ranks = g_df[['Seed', 'Rank', 'Frequency']].sort_values(
        ['Seed', 'Rank', 'Frequency']).drop_duplicates()
    g_pure_ranks.to_csv(
        f"./ranks/g_df_pure_ranks_mod={model_str}_w={w_str}_f={f_str}_int={intensity_str}.csv",
        index=False
    )


for g_n, g_df in grouped_fw_types:
    g_df.sort_values(["Grain Type", "Model", "Intensity", "f-Mark Type", "Weight Type", "Input Value", "PWFCF Value"])
    model = np.unique(g_df['Model'])[0]
    f_type = np.unique(g_df['f-Mark Type'])[0]
    w_type = np.unique(g_df['Weight Type'])[0]
    g_df.to_csv(f"results_splitted_mod={model}_w={w_type}_f={f_type}.csv", index = False)
    #grouped = g_df.groupby(['Seed'])
    # for group_name, group_df in grouped:
    #     print(g_n)
    #     plt.plot(group_df['Input Value'], group_df['PWFCF Value'], label=str(group_name), marker=".")
    # plt.show()
    # plt.close() # TODO tyhle ploty by mohly byt zajimavy do DP
    #assign_the_lexicographic_value(g_df, envelope_count)
# envelopes_df.to_csv(f"./envelope_test_vals_{datetime.now().__str__().replace(':', '-')}.csv")

breakpoint_var=1