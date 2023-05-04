import logging
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from skspatial.objects import Circle, Point

from Geometry.grain import Grain, Segment
from Processes.markings import Mark
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.segment_process import SegmentProcess
from Processes.ball_process import BallProcess
from Processes.point_process import PoissonPointProcess
import utils.const as const

LAMBDA_TEST = 10


RADIUS_TEST_MIN = 0.15
RADIUS_TEST_RANGE = 0.1

R_MIN = RADIUS_TEST_MIN
R_MAX = RADIUS_TEST_MIN + RADIUS_TEST_RANGE

P_ALT = 0.5
win_edge_start_point, win_edge_end_point = - 1 - R_MAX, 1 + R_MAX

_start = datetime.now()
np.random.seed(seed=69)

poisson_point_process = PoissonPointProcess(
    intensity=LAMBDA_TEST, window_edge_start_point=win_edge_start_point,
    window_edge_end_point=win_edge_end_point)

MARKS_TEST = [1 for _ in range(poisson_point_process.points.size)]

RADIUS_TEST = np.random.uniform(low=R_MIN, high=R_MAX)
particles = [
    Particle(
        germ=center, grain_type="ball", germ_inside_the_obs_window=True,
        grain=Circle(point=center, radius=RADIUS_TEST), mark=Mark(mark_type="discrete", mark_value=mark)
    ) for center, mark in zip(poisson_point_process.points, MARKS_TEST)]
ball_process_test = BallProcess(
    germ_intensity=LAMBDA_TEST, particles=particles, max_radius=R_MAX, min_radius=R_MIN, marked=True)

ball_process_test.plot_itself(
    show_germs=True,  win_min=win_edge_start_point, win_max=win_edge_end_point, edge_effects=R_MAX
    )

breakpoint_var=1



