
import numpy as np
import os
from skspatial.objects import Point

from Geometry.grain import Grain, Segment
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.segment_process import SegmentProcess

os.chdir("../")
# ALL FOLLOWING CHECKED FOR CORRECT VALUES
# START_AND_END_POINTS = [
#     (Point([0, 0]), Point([1, 0])),
#     (Point([0, 1/2]), Point([0, 1/4])),
#     (Point([1, 1]), Point([2, 2])),
#     (Point([0, 1]), Point([-1, 0])),
#     (Point([0, 0]), Point([1, 1])),
#     (Point([1, 0]), Point([2, 1]))
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
    for _ in range(100)
]

particles = [
        Particle(
            germ=point[0] / 2 + point[1] / 2,
            grain=Segment(
                start_point=point[0],
                end_point=point[1]
            ),
            grain_type="segment"
        ) for point in START_AND_END_POINTS
    ]
particle_process = SegmentProcess(
    particles=particles, germ_intensity=4, space_dimension=len(START_AND_END_POINTS[0][0])
)
particle_process.plot_itself()
a=1
