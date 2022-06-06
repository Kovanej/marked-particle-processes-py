
import os
import numpy as np
from skspatial.objects import Point, Circle

from Geometry.particle import Particle, Segment
from Processes.markings import Mark
from Processes.point_process import PoissonPointProcess
from Processes.segment_process import SegmentProcess
from Processes.ball_process import BallProcess
import utils.const as const


os.chdir("../.")


poisson_point_process = PoissonPointProcess(
    intensity=const.POISSON_INTENSITY,
    window_edge_start_point=-const.MAX_CIRC_RAD,
    window_edge_end_point=1 + const.MAX_CIRC_RAD
)

# testing for balls & segments
particles_ball_null_model = []
particles_segment_null_model = []
for k in range(len(poisson_point_process.points)):
    # balls
    radius = np.random.random_sample() * (const.MAX_CIRC_RAD - const.MIN_CIRC_RAD) + const.MAX_CIRC_RAD
    ind_mark = np.random.binomial(n=1, p=1/2, size=1)[0]
    mark_null_model = Mark(mark_type="discrete", mark_value=ind_mark)
    particles_ball_null_model.append(
        Particle(
            germ=Point(poisson_point_process.points[k]),
            grain_type="ball",
            grain=Circle(point=Point(poisson_point_process.points[k]), radius=radius),
            mark=mark_null_model
        )
    )
    length = np.random.random_sample() * (const.MAX_SEGMENT_LENGTH - const.MIN_SEGMENT_LENGTH) + const.MIN_SEGMENT_LENGTH
    angle = np.random.random_sample() * (const.MAX_SEGMENT_ANGLE - const.MIN_SEGMENT_ANGLE) + const.MIN_SEGMENT_ANGLE
    particles_segment_null_model.append(
        Particle(
            germ=Point(
                poisson_point_process.points[k] - 1/2 * length * np.array(np.cos(angle), np.sin(angle))
            ),
            grain_type="segment",
            grain=Segment(start_point=poisson_point_process.points[k], length=length, angle=angle),
            mark=mark_null_model
        )
    )


ball_process_null_model = BallProcess(
    germ_intensity=poisson_point_process.intensity, particles=particles_ball_null_model, marked=True
)

ball_process_null_model.plot_itself()

segment_process_null_model = SegmentProcess(
    germ_intensity=poisson_point_process.intensity, particles=particles_segment_null_model
)
segment_process_null_model.plot_itself()
