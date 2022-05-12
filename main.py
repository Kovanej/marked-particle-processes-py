import logging
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from skspatial.objects import Circle, Point

from Geometry.grain import Grain, Segment
from Processes.markings import Mark
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.segment_process import SegmentProcess
from Processes.ball_process import BallProcess
from Processes.point_process import PoissonPointProcess
import utils.const as const


if __name__ == '__main__':
    # TODO resolve trouble with encoding
    # logging.basicConfig(filename='example.log', filemode='w', encoding='utf-8', level=logging.INFO)
    logging.info(f"{datetime.now()}: MAIN SCRIPT RUN STARTED & LOGGER INITIALIZED")
    poisson_point_process = PoissonPointProcess(intensity=const.POISSON_INTENSITY)
    # TODO move to different script
    if const.GRAIN_TYPE == "segment":
        particles_null_model = []
        particles_angle_mark_model = []
        for k in range(len(poisson_point_process.points)):
            angle = np.pi * np.random.random(size=1)[0]
            length = (const.MAX_SEGMENT_LENGTH - const.MIN_SEGMENT_LENGTH) * np.random.random(size=1)[0] + const.MIN_SEGMENT_LENGTH
            mark_null_model = Mark(mark_type="discrete", mark_value=np.random.binomial(n=1, p=1/2, size=1)[0])
            mark_angle_mark_model = Mark(mark_type="discrete", mark_value=np.random.binomial(n=1, p=angle / np.pi, size=1)[0])
            particles_null_model.append(
                Particle(
                    germ=Point(poisson_point_process.points[k]),
                    grain=Segment(
                        start_point=Point(poisson_point_process.points[k]) - [np.cos(angle) * length / 2, np.sin(angle) * length / 2],
                        end_point=Point(poisson_point_process.points[k]) + [np.cos(angle) * length / 2, np.sin(angle) * length / 2]
                    ),
                    grain_type="segment",
                    # mark=mark_null_model
                )
            )
            # particles_angle_mark_model.append(
            #     Particle(
            #         germ=Point(poisson_point_process.points[k]),
            #         grain=Segment(
            #             start_point=Point(poisson_point_process.points[k]) - [np.cos(angle) * length / 2,
            #                                                                   np.sin(angle) * length / 2],
            #             end_point=Point(poisson_point_process.points[k]) + [np.cos(angle) * length / 2,
            #                                                                 np.sin(angle) * length / 2]
            #         ),
            #         grain_type="segment",
            #         mark=mark_angle_mark_model
            #     )
            # )

    particle_process_null_model = SegmentProcess(
        particles=particles_null_model, germ_intensity=poisson_point_process.intensity
    )
    # particle_process_angle_mark_model = SegmentProcess(
    #     particles=particles_angle_mark_model, germ_intensity=poisson_point_process.intensity
    # )

    # particle_process_null_model.compute_the_f_mark_characteristics()
    # particle_process_angle_mark_model.compute_the_f_mark_characteristics()

    particle_process_null_model.plot_itself()
    # particle_process_angle_mark_model.plot_itself()
    logging.info(f"{datetime.now()}: MAIN SCRIPT RUN ENDED")
    brkpnt = "breakpoint here"
