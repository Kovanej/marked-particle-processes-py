import os
import numpy as np
from skspatial.objects import Point, Circle

from Geometry.grain import Grain, Segment
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.ball_process import BallProcess
import utils.const as const


CENTERS_AND_RADII = [
    [
        # (Point([0.5, 0.5]), 1 / dv),
        # (Point([0.25, 0.75]), 1 / dv),
        # (Point([0.25, 0.25]), 1 / dv),
        # (Point([0.75, 0.75]), 1 / dv),
        # (Point([0.75, 0.25]), 1 / dv)
        (Point([1/4, 1/3]), 1 / dv),
        (Point([1/4, 2/3]), 1 / dv),
        (Point([2/4, 1/3]), 1 / dv),
        (Point([2/4, 2/3]), 1 / dv),
        (Point([3/4, 1/3]), 1 / dv),
        (Point([3/4, 2/3]), 1 / dv),
    ] for dv in range(1, 200)
]
CENTERS_AND_RADII = [c_r for inside_list in CENTERS_AND_RADII for c_r in inside_list]
os.chdir("../")

particles = [
    Particle(
        germ=center,
        grain_type="ball",
        grain=Circle(point=center, radius=radius)
    )
    for center, radius in CENTERS_AND_RADII
]

ball_process_test = BallProcess(
    particles=particles, germ_intensity=len(CENTERS_AND_RADII)
)
ball_process_test.plot_itself()

# ball_process_test.compute_the_f_mark_characteristics()

a=1
