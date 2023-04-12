
from datetime import datetime
import logging
import pandas as pd
import numpy as np
from skspatial.objects import Point, Vector, Circle

from Geometry.particle import Particle
from Processes.point_process import PoissonPointProcess
from Processes.particle_process import ParticleProcess
from Processes.segment_process import SegmentProcess
from Processes.ball_process import BallProcess
from Processes.markings import Mark
from Geometry.grain import Segment
from utils.config_parser import ConfigParser
import utils.const as const
from utils.results_saver import ResultSaver
from tests.wider_window_simulation_tests import simulate_the_processes


def execute_from_config(config_parser: ConfigParser):
    result_savers = []
    for seed in range(config_parser.initial_seed, config_parser.initial_seed + config_parser.number_of_realizations):
        config_parser.initialize_the_processes(seed=seed)
        result_saver = config_parser.return_the_result_saver(seed=seed)
        result_savers.append(result_saver)
        print(f"Results computed for seed={seed}.")
    result_savers_whole_df = pd.concat([r_s.results_all_df for r_s in result_savers])
    return result_savers, result_savers_whole_df


def overnight_computations():
    result_savers = []
    for seed in range(200, 220):
        print(f"Seed no {seed}: {datetime.now()}")
        result_saver = simulate_the_processes(alphas_list=[0, 0.25, 0.5, 0.75, 1], seed=seed)
        result_savers.append(result_saver)
    return result_savers
