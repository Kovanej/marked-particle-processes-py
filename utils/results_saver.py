import os
from datetime import datetime as dt
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from typing import Dict, Optional

import utils.const as const


class ResultSaver(object):

    def __init__(self):
        self.result_dict = {
            "Seed": [],
            "Grain Type": [],
            "Model": [],
            "f-Mark Type": [],
            "Weight Type": [],
            "Value": [],
            "Permutation Test Count": [],
            "Left p-value": [],
            "Right p-value": [],
            "Both Sided p-value": []
        }
        self.results_grouped_df_dict: Dict = {}
        self.results_all_df = None
        self.results_grouped_by_df = None

    def save_the_results(
            self, model_name: str, grain_type: str, permutations_count: int, quantile_dict: Dict, value_dict: Dict,
            seed: Optional[int] = None
    ):
        for (key, val) in quantile_dict.items():
            self.result_dict["Seed"].append(seed)
            self.result_dict["Grain Type"].append(grain_type)
            self.result_dict["Model"].append(model_name)
            self.result_dict["f-Mark Type"].append(key[0])
            self.result_dict["Weight Type"].append(key[1])
            self.result_dict["Value"].append(value_dict[key])
            self.result_dict["Permutation Test Count"].append(permutations_count)
            self.result_dict["Left p-value"].append(val)
            self.result_dict["Right p-value"].append(1 - val)
            self.result_dict["Both Sided p-value"].append(2 * min(val, 1 - val))

    def save_to_pandas(self, save_csv: bool = const.SAVE_RESULTS_TO_CSV):
        self.results_all_df = pd.DataFrame(self.result_dict)
        results_all_grouped_by = self.results_all_df.groupby([
            "Grain Type", "Model", "f-Mark Type", "Weight Type"
        ])
        self.results_grouped_df_dict = {
            k: results_all_grouped_by.get_group(k) for k in results_all_grouped_by.groups
        }
        # self.results_grouped_by_df = pd.DataFrame(self.results_grouped_df_dict)
        dtm = str(dt.now()).replace(":", "-")
        if save_csv:
            if "results" not in os.listdir():
                os.mkdir("./results")
            self.results_all_df.to_csv(f"results/results_{dtm}.csv", index=False)

    def pickle_the_result_dataframes(self):
        f_results = open(f'pickles/results_{str(dt.now()).replace(":", "-")}.txt', 'wb')
        pickle.dump(self.results_all_df, f_results)
        f_results_grouped = open(f'pickles/results_grouped_dict{str(dt.now()).replace(":", "-")}.txt', 'wb')
        pickle.dump(self.results_grouped_df_dict, f_results_grouped)

    def plot_the_p_values(
            self, grain_type: str, model: str, f_mark_type: str,
            weight_type: str
    ):
        ax = self.results_grouped_df_dict[(grain_type, model, f_mark_type, weight_type)]["Both Sided p-value"].hist(
            bins=8
        )
        # TODO this doesn't work if the list of columns is passed ^^
        fig = ax.get_figure()
        ax.set_title(f"{grain_type, model, f_mark_type, weight_type}")
        plt.show()
        plt.close(fig)

