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
TESTED_INTENSITY = 10
MAX_SEGMENT_LENGTH = 0.3
MIN_SEGMENT_LENGTH = 0.1
MAX_CIRC_RAD = 0.15
MIN_CIRC_RAD = 0.05

CENTERS_AND_RADII_LISTS = [[
    (Point([1 / 4, 1 / 2]), 1 / 4),
    (Point([2 / 4, 1 / 2]), 1 / 4),
    (Point([3 / 4, 1 / 2]), 1 / 4),
]

]  # for k in range(1, EPICENTER_COUNT)]
CENTERS_AND_RADII = [c_r for inside_list in CENTERS_AND_RADII_LISTS for c_r in inside_list]

START_AND_END_POINTS = [
    (Point([- 1/2, 0]), Point([1, 0])),
    (Point([-1, -1]), Point([1, 1])),
    (Point([0, 10]), Point([0, 15])),
]
SEGMENT_MARKS = [
    1, 2, 3
]


logging.basicConfig(filename='example.log', filemode='w', level=logging.INFO)
# poisson_point_process = PoissonPointProcess(intensity=TESTED_INTENSITY)
if TESTED_GRAIN_TYPE == "ball":

    particles = [
        Particle(
            germ=center,
            grain_type="ball",
            grain=Circle(point=center, radius=radius),
            mark=Mark(mark_type="continuous", mark_value=center[0])
        )
        for center, radius in CENTERS_AND_RADII
    ]

    test_process = BallProcess(
        particles=particles, germ_intensity=len(CENTERS_AND_RADII), marked=True
    )
elif TESTED_GRAIN_TYPE == "segment":
    particles = [
        Particle(
            grain_type="segment",
            germ=1/2 * (START_AND_END_POINTS[k][0] + START_AND_END_POINTS[k][1]),
            grain=Segment(
                start_point=START_AND_END_POINTS[k][0],
                end_point=START_AND_END_POINTS[k][1]
            ),
            mark=Mark(mark_type="discrete", mark_value = SEGMENT_MARKS[k])
        )
        for k in range(len(START_AND_END_POINTS))
    ]
    test_process = SegmentProcess(particles=particles, germ_intensity=len(START_AND_END_POINTS), marked=True)
else:
    raise NotImplementedError()

test_process.plot_itself()

test_process.compute_the_f_mark_characteristics()

test_process.perform_the_permutation_test_for_f_mark_characteristics()

brkpnt = "breakpoint here"
