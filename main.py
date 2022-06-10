import logging
from matplotlib import pyplot as plt
import pickle
import pandas as pd
import time

from executions import first_blood, overnight_computations
import utils.const as const
from utils.results_saver import ResultSaver

LOAD_PICKLE = False


if __name__ == '__main__':
    # if LOAD_PICKLE:
    #     result_saver = ResultSaver()
    #     file_results = open("results.txt", 'rb')
    #     result_saver.results_all_df = pickle.load(file_results)
    #     file_grouped_results = open("results_grouped.txt", 'rb')
    #     result_saver.results_grouped_by_df = pickle.load(file_grouped_results)
    # else:
    #     result_saver = first_blood(number_of_seeds=4)
    #     print("first blood done")
    # if const.PLOT_THE_P_VALUES:
    #     for key, val in result_saver.results_grouped_df_dict.items():
    #         grain_type = key[0]
    #         model = key[1]
    #         f_mark_type = key[2]
    #         weight_type = key[3]
    #         result_saver.plot_the_p_values(
    #             grain_type=grain_type, model=model, f_mark_type=f_mark_type, weight_type=weight_type
    #         )
    result_savers = overnight_computations()
    brkpnt = "breakpoint is here"
