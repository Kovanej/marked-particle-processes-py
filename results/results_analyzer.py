
import matplotlib.pyplot as plt
import os
import pandas as pd


class ResultsAnalyzer(object):

    def __init__(self):
        self.dfs_list = []
        self.files_to_analyze = [f for f in os.listdir() if ".csv" in f]
        for f in self.files_to_analyze:
            df = pd.read_csv(f)
            self.dfs_list.append(df)
        self.df = pd.concat(self.dfs_list)
        self.grouped_keys = [
            (grain_type, model, f_mark_type, weight_type)
            for grain_type in self.df['Grain Type'].unique() for model in self.df['Model'].unique()
            for f_mark_type in self.df['f-Mark Type'].unique() for weight_type in self.df['Weight Type'].unique()
        ]
        self.grouped_dfs = {
            (grain_type, model, f_mark_type, weight_type): self.df[
                (self.df['Grain Type'] == grain_type) & (self.df['Model'] == model) &
                (self.df['f-Mark Type'] == f_mark_type) & (self.df['Weight Type'] == weight_type)
            ]
            for grain_type, model, f_mark_type, weight_type in self.grouped_keys
        }
        self.grouped_dfs = {k: df for k, df in self.grouped_dfs.items() if not df.empty}
        self.grouped_keys = [k for k in self.grouped_dfs.keys()]

    def plot_the_histograms(self):
        for grain_type, model, f_mark_type, weight_type in self.grouped_keys:
            ax = self.grouped_dfs[(grain_type, model, f_mark_type, weight_type)]["Both Sided p-value"].hist(
                bins=8
            )
            # TODO this doesn't work if the list of columns is passed ^^
            fig = ax.get_figure()
            ax.set_title(f"{grain_type, model, f_mark_type, weight_type}")
            plt.show()


results_analyzer = ResultsAnalyzer()
results_analyzer.plot_the_histograms()