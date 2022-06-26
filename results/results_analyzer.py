
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


results_analyzer = ResultsAnalyzer()
a=1
