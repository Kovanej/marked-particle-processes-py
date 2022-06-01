import logging
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
from skspatial.objects import Circle, Point

from Geometry.grain import Grain, Segment
from Processes.markings import Mark
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.segment_process import SegmentProcess
from Processes.ball_process import BallProcess
from Processes.point_process import PoissonPointProcess
import utils.const as const


os.chdir("../")


TESTED_GRAIN_TYPE = "segment"
# TESTED_INTENSITY = 10
MAX_SEGMENT_LENGTH = 0.2
MIN_SEGMENT_LENGTH = 0.1
MAX_CIRC_RAD = 0.15
MIN_CIRC_RAD = 0.05
#
# CENTERS_AND_RADII_LISTS = [[
#     (Point([1 / 4, 1 / 2]), 1 / 4),
#     (Point([2 / 4, 1 / 2]), 1 / 4),
#     (Point([1, 1 / 2]), 1 / 4),
# ]
#
# ]  # for k in range(1, EPICENTER_COUNT)]
# CENTERS_AND_RADII = [c_r for inside_list in CENTERS_AND_RADII_LISTS for c_r in inside_list]
#
# START_AND_END_POINTS = [
#     (Point([- 1/2, 0]), Point([1, 0])),
#     (Point([-1, -1]), Point([1, 1])),
#     (Point([0, 10]), Point([0, 15])),
# ]
# SEGMENT_MARKS = [
#     1, 2, 3
# ]
# MARKS_BALLS = [
#     1, 2, 3
# ]

# TESTING POISSON PROCESS
POISSON_TEST_INTENSITY = 100
poisson_point_process = PoissonPointProcess(intensity=POISSON_TEST_INTENSITY)
segment_angles = [np.pi * np.random.random_sample() for _ in range(len(poisson_point_process.points))]
segment_lengths = [
    MIN_SEGMENT_LENGTH + (MAX_SEGMENT_LENGTH - MIN_SEGMENT_LENGTH) * np.random.random_sample()
    for _ in range(len(poisson_point_process.points))
]
ball_radii = [
    MIN_CIRC_RAD + (MAX_CIRC_RAD - MIN_CIRC_RAD) * np.random.random_sample()
    for _ in range(len(poisson_point_process.points))
]
marks_null = np.random.binomial(size=len(poisson_point_process.points), n=1, p=0.5)
marks_angles_beatles_discrete = [
    np.random.binomial(size=1, n=1, p=np.array(1/np.pi * segment_angles[k]))[0]
    for k in range(len(poisson_point_process.points))
]

logging.basicConfig(filename='example.log', filemode='w', level=logging.INFO)
# poisson_point_process = PoissonPointProcess(intensity=TESTED_INTENSITY)
if TESTED_GRAIN_TYPE == "ball":
    particles = [
        Particle(
            germ=Point(poisson_point_process.points[k]),
            grain_type="ball",
            grain=Circle(point=Point(poisson_point_process.points[k]), radius=ball_radii[k]),
            mark=Mark(mark_type="discrete", mark_value=marks_null[k])
        )
        for k in range(len(poisson_point_process.points))
    ]

    test_process_null = BallProcess(
        particles=particles, germ_intensity=poisson_point_process.intensity, marked=True
    )
elif TESTED_GRAIN_TYPE == "segment":
    particles = [
        Particle(
            grain_type="segment",
            germ=Point(poisson_point_process.points[k]),
            grain=Segment(
                start_point=poisson_point_process.points[k] - 1 / 2 * segment_lengths[k] * np.array(
                    [np.cos(segment_angles[k]), np.sin(segment_angles[k])]
                ),
                angle=segment_angles[k],
                length=segment_lengths[k]
            ),
            mark=Mark(mark_type="discrete", mark_value=marks_null[k])
        )
        for k in range(len(poisson_point_process.points))
    ]
    test_process_null = SegmentProcess(particles=particles, germ_intensity=poisson_point_process.intensity, marked=True)
    particles_angle_model_discrete = [
        Particle(
            grain_type="segment",
            germ=Point(poisson_point_process.points[k]),
            grain=Segment(
                start_point=poisson_point_process.points[k] - 1 / 2 * segment_lengths[k] * np.array(
                    [np.cos(segment_angles[k]), np.sin(segment_angles[k])]
                ),
                angle=segment_angles[k],
                length=segment_lengths[k]
            ),
            mark=Mark(mark_type="discrete", mark_value=marks_angles_beatles_discrete[k])
        )
        for k in range(len(poisson_point_process.points))
    ]
    test_process_angles_discrete = SegmentProcess(
        particles=particles_angle_model_discrete, germ_intensity=poisson_point_process.intensity, marked=True
    )
    particles_angle_model_continuous = [
        Particle(
            grain_type="segment",
            germ=Point(poisson_point_process.points[k]),
            grain=Segment(
                start_point=poisson_point_process.points[k] - 1 / 2 * segment_lengths[k] * np.array(
                    [np.cos(segment_angles[k]), np.sin(segment_angles[k])]
                ),
                angle=segment_angles[k],
                length=segment_lengths[k]
            ),
            mark=Mark(mark_type="discrete", mark_value=segment_angles[k] / np.pi)
        )
        for k in range(len(poisson_point_process.points))
    ]
    test_process_angles_continuous = SegmentProcess(
        particles=particles_angle_model_continuous, germ_intensity=poisson_point_process.intensity, marked=True
    )
else:
    raise NotImplementedError()

test_process_null.plot_itself()
test_process_null.compute_the_f_mark_characteristics()
test_process_null.perform_the_permutation_test_for_f_mark_characteristics()
print("---NULL MODEL---")
for _, __ in test_process_null.f_mark_statistics.items():
    print(f"STATISTIC: {_}, VALUE: {__}, QUANTILE: {test_process_null.f_mark_statistics_quantiles[_]}")

if TESTED_GRAIN_TYPE == "segment":
    test_process_angles_discrete.plot_itself()
    test_process_angles_discrete.compute_the_f_mark_characteristics()
    test_process_angles_discrete.perform_the_permutation_test_for_f_mark_characteristics()
    print("---ANGLE MODEL DISCRETE---")
    for _, __ in test_process_angles_discrete.f_mark_statistics.items():
        print(f"STATISTIC: {_}, VALUE: {__}, QUANTILE: {test_process_angles_discrete.f_mark_statistics_quantiles[_]}")

    test_process_angles_continuous.plot_itself()
    test_process_angles_continuous.compute_the_f_mark_characteristics()
    test_process_angles_continuous.perform_the_permutation_test_for_f_mark_characteristics()
    print("---ANGLE MODEL DISCRETE---")
    for _, __ in test_process_angles_continuous.f_mark_statistics.items():
        print(f"STATISTIC: {_}, VALUE: {__}, QUANTILE: {test_process_angles_continuous.f_mark_statistics_quantiles[_]}")


brkpnt = "breakpoint here"
