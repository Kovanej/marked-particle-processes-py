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

NUMBER_OF_SEEDS = 50
LAMBDA_TEST = 30


RADIUS_TEST_MIN = 0.3
RADIUS_TEST_RANGE = 0.2

R_MIN = RADIUS_TEST_MIN
R_MAX = RADIUS_TEST_MIN + RADIUS_TEST_RANGE

P_ALT = 0.5
win_edge_start_point, win_edge_end_point = - 1 - 3 * R_MAX, 1 + 3 * R_MAX

#F_MARK_STATISTICS_TO_COMPUTE = const.F_MARK_COMBINATIONS
F_MARK_STATISTICS_TO_COMPUTE = [('product', 'intersection')]

f_mark_statistics = {k: [] for k in F_MARK_STATISTICS_TO_COMPUTE}


def simple_pwfcf_for_uniform_distribution(value: float, lam: float, p: float, a: float, b: float):
    return (lam * p * np.pi) ** 2 * (1 / 180) * (
        (210 * b ** 2 + 300 * a * b + 210 * a ** 2) * value ** 2 +
        (270 * b ** 3 + 450 * a ** 2 * b + 450 * a * b ** 2 + 210 * a ** 3) * value +
        (101 * b ** 4 + 166 * a * b ** 3 + 186 * a ** 2 * b ** 2 + 166 * a ** 3 * b + 101 * a ** 4)
    )


for _ in range(NUMBER_OF_SEEDS):
    _start = datetime.now()
    np.random.seed(seed=_)

    poisson_point_process = PoissonPointProcess(
        intensity=LAMBDA_TEST, window_edge_start_point=win_edge_start_point,
        window_edge_end_point=win_edge_end_point)

    MARKS_TEST = [np.random.binomial(1, P_ALT) for _ in range(poisson_point_process.points.size)]

    RADIUS_TEST = np.random.uniform(low=R_MIN, high=R_MAX)
    particles = [
        Particle(
            germ=center, grain_type="ball", germ_inside_the_obs_window=True,
            grain=Circle(point=center, radius=RADIUS_TEST), mark=Mark(mark_type="discrete", mark_value=mark)
        ) for center, mark in zip(poisson_point_process.points, MARKS_TEST)]
    ball_process_test = BallProcess(
        germ_intensity=LAMBDA_TEST, particles=particles, max_radius=R_MAX, min_radius=R_MIN, marked=True)
    ball_process_test.compute_the_f_mark_characteristics(set_of_f_mark_combinations=F_MARK_STATISTICS_TO_COMPUTE)
    for k in F_MARK_STATISTICS_TO_COMPUTE:
        f_mark_statistics[k].append(ball_process_test.f_mark_statistics[k])
    _end = datetime.now()
    print(f"Computation of f-mark stats for seed: {_} ended in {_end}.\n Runtime of this computation: {_end - _start}.")


f_mark_statistics_pds = {k: pd.DataFrame(v) for k, v in f_mark_statistics.items()}

# temporary solution
t_computed = list(f_mark_statistics[('product', 'intersection')][0].keys())
f_mark_statistics_max = {}
f_mark_statistics_mean = {}
f_mark_statistics_min = {}
for t in t_computed:
    val = 0
    f_mark_statistics_min[t] = np.inf
    f_mark_statistics_max[t] = 0
    for _ in f_mark_statistics.keys():
        for dct in f_mark_statistics[_]:
            one_val = dct[t]
            val += one_val
            if one_val < f_mark_statistics_min[t]:
                f_mark_statistics_min[t] = one_val
            if one_val > f_mark_statistics_max[t]:
                f_mark_statistics_max[t] = one_val
        f_mark_statistics_mean[t] = val / len(f_mark_statistics[_])

# Create a list of sorted inputs and outputs
inputs = sorted(f_mark_statistics_mean.keys())
outputs_mean = [f_mark_statistics_mean[k] for k in inputs]
outputs_min = [f_mark_statistics_min[k] for k in inputs]
outputs_max = [f_mark_statistics_max[k] for k in inputs]
outputs_theoretical = [
    simple_pwfcf_for_uniform_distribution(value=t, lam=LAMBDA_TEST, p=P_ALT, a=R_MIN, b=R_MAX) for t in inputs
]

# Create a line plot
plt.plot(inputs, outputs_mean, marker='o', label="mean")
plt.plot(inputs, outputs_min, marker='o', label="min")
plt.plot(inputs, outputs_max, marker='o', label="max")
plt.plot(inputs, outputs_theoretical, marker='o', label="theoretical")

plt.legend()

# Add axis labels and a title
plt.xlabel('Input')
plt.ylabel('Output')
plt.title('S_{w, f}(r)')

ball_process_test.plot_itself()

# Show the plot
plt.show()


break_point_var=1

