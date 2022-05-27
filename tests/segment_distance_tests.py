
import logging
import numpy as np
import os
from skspatial.objects import Point

from Geometry.grain import Grain, Segment
from Plotting.plotting import plot_the_grains
from Processes.markings import Mark
from Geometry.particle import Particle
from Processes.segment_process import SegmentProcess

logging.basicConfig(filename='example.log', filemode='w', level=logging.INFO)

os.chdir("../")
# ALL FOLLOWING CHECKED FOR CORRECT VALUES
np.random.seed(69)
# START_AND_END_POINTS = [
#     (Point([- 1/2, 0]), Point([1, 0])),
#     # --- before tested from here
#     (Point([-1, -1]), Point([1, 1])),
#     (Point([0, 10]), Point([0, 15])),
#     (Point([2, 1]), Point([5/2, 0])),
#     (Point([0, 9]), Point([1, 1])),
#     (Point([9, 0]), Point([1, 1])),
#     # --- up till here
#     # (Point([0, 2]), Point([2, 0])),
#     # (Point([1, 1]), Point([2, 2])),
#     # (Point([0, 1]), Point([-1, 0])),
#     # (Point([0, -2]), Point([0, 1])),
# ]

# ALSO TESTED FOR A COUPLE OF SEGMENTS IN R^3
# START_AND_END_POINTS = [
#     (Point([0, 0, 0]), Point([0, 0, -1])),
#     (Point([-1, 2, 0]), Point([1, 2, 0])),
#     (Point([-1, 0, -2]), Point([2, 0, -2])),
#     (Point([-1, 0, 2]), Point([2, 0, 2])),
# ]

START_AND_END_POINTS = [
    (
        Point([np.random.random_sample(), np.random.random_sample()]),
        Point([np.random.random_sample(), np.random.random_sample()])
    )
    for _ in range(1000)
]

MAX_SEGMENT_LENGTH_TEST = 0.4
MIN_SEGMENT_LENGTH_TEST = 0.1
NO_OF_SEGMENTS = 100

# particles = [
#         Particle(
#             germ=point[0] / 2 + point[1] / 2,
#             grain=Segment(
#                 start_point=point[0],
#                 end_point=point[1]
#             ),
#             grain_type="segment"
#         ) for point in START_AND_END_POINTS
#     ]

ANGLE_DEPENDENCY_WEIGHT = 0

angles = [
    np.pi * np.random.random_sample() for _ in range(NO_OF_SEGMENTS)
]
lengths = [
    MIN_SEGMENT_LENGTH_TEST + np.random.random_sample() * (MAX_SEGMENT_LENGTH_TEST - MIN_SEGMENT_LENGTH_TEST)
    for _ in range(NO_OF_SEGMENTS)
]
start_points = [
    Point(
        [np.random.random_sample() - np.cos(angles[_]) * lengths[_] / 2,
         np.random.random_sample() - np.sin(angles[_]) * lengths[_] / 2]
    )
    for _ in range(NO_OF_SEGMENTS)
]
marks = [
    ANGLE_DEPENDENCY_WEIGHT * angles[k] + (1 - ANGLE_DEPENDENCY_WEIGHT) * np.random.random_sample()
    for k in range(NO_OF_SEGMENTS)
]

particles = [
    Particle(
        grain_type="segment",
        germ=start_points[k],
        grain=Segment(
            start_point=start_points[k],
            angle=angles[k],
            length=lengths[k]
        ),
        mark=Mark(mark_type="continuous", mark_value=marks[k])
    )
    for k in range(NO_OF_SEGMENTS)
]

particle_process = SegmentProcess(
    particles=particles, germ_intensity=NO_OF_SEGMENTS, space_dimension=len(START_AND_END_POINTS[0][0]), marked=True
)
particle_process.plot_itself()
# comparison = particle_process.particles_distance_matrix == particle_process.particles_distance_matrix_vectorized
# equal_arrays = comparison.all()
# print(equal_arrays)
particle_process.compute_the_f_mark_characteristics()
particle_process.perform_the_permutation_test_for_f_mark_characteristics()
print(f"{particle_process.f_mark_statistics}")
print(f"{particle_process.f_mark_statistics_quantiles}")
