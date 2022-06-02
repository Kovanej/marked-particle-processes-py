
from datetime import datetime as dt
from typing import Dict, Optional
import pandas as pd
import pickle

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
            "Quantile": []
        }
        self.results_df = None
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
            self.result_dict["Quantile"].append(val)

    def save_to_pandas(self, save_csv: bool = const.SAVE_RESULTS_TO_CSV):
        self.results_df = pd.DataFrame(self.result_dict)
        dtm = str(dt.now()).replace(":", "-")
        if save_csv:
            self.results_df.to_csv(f"results_{dtm}.csv", index=False)
        self.results_grouped_by_df = self.results_df.groupby(
            by=["Grain Type", "Model", "f-Mark Type", "Weight Type"]
        )['Quantile'].mean()

    def pickle_the_result_dataframes(self):
        f_results = open('results.txt', 'wb')
        pickle.dump(self.results_df, f_results)
        f_results_grouped = open('results_grouped.txt', 'wb')
        pickle.dump(self.results_grouped_by_df, f_results_grouped)


