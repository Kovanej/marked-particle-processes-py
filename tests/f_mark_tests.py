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
TESTED_INTENSITY = 300
MAX_SEGMENT_LENGTH = 0.4
MIN_SEGMENT_LENGTH = 0.1
MAX_CIRC_RAD = 0.1
MIN_CIRC_RAD = 0.05


logging.basicConfig(filename='example.log', filemode='w', level=logging.INFO)
logging.info(f"{datetime.now()}: MAIN SCRIPT RUN STARTED & LOGGER INITIALIZED")
poisson_point_process = PoissonPointProcess(intensity=TESTED_INTENSITY)
if TESTED_GRAIN_TYPE == "segment":
    particles_null_model = []
    particles_angle_mark_model = []
    for k in range(len(poisson_point_process.points)):
        angle = np.pi * np.random.random_sample()
        length = (MAX_SEGMENT_LENGTH - MIN_SEGMENT_LENGTH) * np.random.random_sample() + MIN_SEGMENT_LENGTH
        mark_null_model = Mark(
            mark_type="discrete",
            mark_value=np.random.random_sample()
            # mark_value=np.random.binomial(n=1, p=1/2, size=1)[0]
        )
        mark_angle_mark_model = Mark(
            mark_type="discrete",
            mark_value=angle / np.pi
            # mark_value=np.random.binomial(n=1, p=angle / np.pi, size=1)[0]
        )
        particles_null_model.append(
            Particle(
                germ=Point(poisson_point_process.points[k]),
                grain=Segment(
                    start_point=Point(
                        poisson_point_process.points[k]
                    ) - [np.cos(angle) * length / 2, np.sin(angle) * length / 2],
                    end_point=Point(
                        poisson_point_process.points[k]
                    ) + [np.cos(angle) * length / 2, np.sin(angle) * length / 2]
                ),
                grain_type="segment",
                mark=mark_null_model
            )
        )
        particles_angle_mark_model.append(
            Particle(
                germ=Point(poisson_point_process.points[k]),
                grain=Segment(
                    start_point=Point(poisson_point_process.points[k]) - [np.cos(angle) * length / 2,
                                                                          np.sin(angle) * length / 2],
                    end_point=Point(poisson_point_process.points[k]) + [np.cos(angle) * length / 2,
                                                                        np.sin(angle) * length / 2]
                ),
                grain_type="segment",
                mark=mark_angle_mark_model
            )
        )

    particle_process_null_model = SegmentProcess(
        particles=particles_null_model, germ_intensity=poisson_point_process.intensity, marked=True
    )
    particle_process_non_null_model = SegmentProcess(
        particles=particles_angle_mark_model, germ_intensity=poisson_point_process.intensity, marked=True
    )

elif TESTED_GRAIN_TYPE == "ball":
    particles_null_model = []
    particles_radius_mark_model = []
    for k in range(len(poisson_point_process.points)):
        radius = np.random.random_sample() * (MAX_CIRC_RAD - MIN_CIRC_RAD) + MAX_CIRC_RAD
        ind_mark = np.random.random_sample() * (MAX_CIRC_RAD - MIN_CIRC_RAD) + MAX_CIRC_RAD
        mark_null_model = Mark(mark_type="continuous", mark_value=ind_mark)
        mark_radius_model = Mark(mark_type="continuous", mark_value=radius)
        particles_null_model.append(
            Particle(
                germ=Point(poisson_point_process.points[k]),
                grain_type="ball",
                grain=Circle(point=Point(poisson_point_process.points[k]), radius=ind_mark),
                mark=mark_null_model
            )
        )
        particles_radius_mark_model.append(
            Particle(
                germ=Point(poisson_point_process.points[k]),
                grain_type="ball",
                grain=Circle(point=Point(poisson_point_process.points[k]), radius=ind_mark),
                mark=mark_radius_model
            )
        )

    particle_process_null_model = BallProcess(
        germ_intensity=poisson_point_process.intensity, particles=particles_null_model, marked=True
    )
    particle_process_non_null_model = BallProcess(
        germ_intensity=poisson_point_process.intensity, particles=particles_radius_mark_model, marked=True
    )
else:
    raise NotImplementedError()

particle_process_null_model.compute_the_f_mark_characteristics()
particle_process_non_null_model.compute_the_f_mark_characteristics()

particle_process_null_model.plot_itself()
particle_process_non_null_model.plot_itself()

brkpnt = "breakpoint here"

print(f"NULL MODEL VALUES:")
for _, __ in particle_process_null_model.f_mark_statistics.items():
    print(_, __)
print(f"NON-NULL MODEL VALUES:")
for _, __ in particle_process_non_null_model.f_mark_statistics.items():
    print(_, __)
