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
    if LOAD_PICKLE:
        result_saver = ResultSaver()
        file_results = open("results.txt", 'rb')
        result_saver.results_all_df = pickle.load(file_results)
        file_grouped_results = open("results_grouped.txt", 'rb')
        result_saver.results_grouped_by_df = pickle.load(file_grouped_results)
    else:
        result_savers = overnight_computations()
    brkpnt = "breakpoint is here"
