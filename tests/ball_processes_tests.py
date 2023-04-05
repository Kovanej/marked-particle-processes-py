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

NO_OF_CIRCLES = 10
NO_OF_INSIDE_CIRCLES = 15

RANDOM_CENTERS = [
    (np.random.random_sample(), np.random.random_sample())
    for _ in range(NO_OF_CIRCLES)
]

os.chdir("../")

EPICENTER_COUNT = 10


NO_OF_INSIDE_CIRCLES = 2
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
        germ_inside_the_obs_window=True,
        grain=Circle(point=center, radius=radius),
        mark=Mark(mark_type="continuous", mark_value=radius)
    )
    for center, radius in CENTERS_AND_RADII
]

ball_process_test = BallProcess(
    particles=particles, germ_intensity=len(CENTERS_AND_RADII), marked=True, min_radius=np.sqrt(2), max_radius=np.sqrt(2)
)
ball_process_test.plot_itself()

print(f"{ball_process_test.f_mark_statistics}")
