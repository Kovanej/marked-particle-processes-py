import logging
import pickle
from executions import first_blood


LOAD_PICKLE = True


if __name__ == '__main__':
    if LOAD_PICKLE:
        file_results = open("results.txt", 'rb')
        test_results = pickle.load(file_results)
        file_grouped_results = open("results_grouped.txt", 'rb')
        test_grouped_results = pickle.load(file_grouped_results)
    else:
        result_saver = first_blood(number_of_seeds=100)
        print("first blood done")
    brkpnt = "breakpoint is here"
