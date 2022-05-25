import logging
import os
import numpy as np
from skspatial.objects import Point, Circle

from Geometry.grain import Grain, Segment
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.ball_process import BallProcess
from Processes.markings import Mark
import utils.const as const


logging.basicConfig(filename='example.log', filemode='w', level=logging.INFO)

NO_OF_CIRCLES = 300
NO_OF_INSIDE_CIRCLES = 15

RANDOM_CENTERS = [
    (np.random.random_sample(), np.random.random_sample())
    for _ in range(NO_OF_CIRCLES)
]
# [
#     [
#         # (Point([1/4, 1/3]), 1 / dv),
#         # (Point([1/4, 2/3]), 1 / dv),
#         # (Point([2/4, 1/3]), 1 / dv),
#         # (Point([2/4, 2/3]), 1 / dv),
#         # (Point([3/4, 1/3]), 1 / dv),
#         # (Point([3/4, 2/3]), 1 / dv),
#         (Point([0.5, 0.5]), 1 / dv),
#         (Point([0.25, 0.75]), 1 / dv),
#         (Point([0.25, 0.25]), 1 / dv),
#         (Point([0.75, 0.75]), 1 / dv),
#         (Point([0.75, 0.25]), 1 / dv),
#     ]
#     for dv in range(1, NO_OF_INSIDE_CIRCLES)
# ] + \

# CENTERS_AND_RADII = [[
#     (Point([RANDOM_CENTERS[i][0], RANDOM_CENTERS[i][1]]), 1 / (dv + 1)),
#     ] for dv in range(1, NO_OF_INSIDE_CIRCLES) for i in range(NO_OF_CIRCLES)
# ]

os.chdir("../")

EPICENTER_COUNT = 10

# CENTERS_AND_RADII_LISTS = [[
#     (Point([np.random.random_sample() ** k, 1 - np.random.random_sample() ** k]), 1 / (2 + np.sqrt(np.random.random_sample() + dv))),
#     (Point([np.random.random_sample() ** k, np.random.random_sample() ** k]), 1 / (2 + np.sqrt(np.random.random_sample() + dv))),
#     (Point([1 - np.random.random_sample() ** k, 1 - np.random.random_sample() ** k]), 1 / (2 + np.sqrt(np.random.random_sample() + dv))),
#     (Point([1 - np.random.random_sample() ** k, np.random.random_sample() ** k]), 1 / (2 + np.sqrt(np.random.random_sample() + dv))),
#     (Point([np.random.random_sample(), 1 - 1/k]), 1 / (2 + np.sqrt(np.random.random_sample() + dv))),
#     (Point([np.random.random_sample(), 1/k]), 1 / (2 + np.sqrt(np.random.random_sample() + dv))),
#     (Point([1 - 1/k, np.random.random_sample()]), 1 / (2 + np.sqrt(np.random.random_sample() + dv))),
#     (Point([1/k, np.random.random_sample()]), 1 / (2 + np.sqrt(np.random.random_sample() + dv))),
# ] for dv in range(1, NO_OF_INSIDE_CIRCLES) for k in range(1, EPICENTER_COUNT)]
# CENTERS_AND_RADII = [c_r for inside_list in CENTERS_AND_RADII_LISTS for c_r in inside_list]


NO_OF_INSIDE_CIRCLES = 100
CENTERS_AND_RADII_LISTS = [[
    (Point([0, 0]), np.sqrt(2) / np.power(dv, 1/5)),
    (Point([1, 1]), np.sqrt(2) / np.power(dv, 1/5))
] for dv in range(1, NO_OF_INSIDE_CIRCLES)
]  # for k in range(1, EPICENTER_COUNT)]
CENTERS_AND_RADII = [c_r for inside_list in CENTERS_AND_RADII_LISTS for c_r in inside_list]

particles = [
    Particle(
        germ=center,
        grain_type="ball",
        grain=Circle(point=center, radius=radius),
        mark=Mark(mark_type="continuous", mark_value=radius)
    )
    for center, radius in CENTERS_AND_RADII
]

ball_process_test = BallProcess(
    particles=particles, germ_intensity=len(CENTERS_AND_RADII), marked=True
)
ball_process_test.plot_itself()

# test_process.compute_the_f_mark_characteristics()

a=1
