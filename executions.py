
import numpy as np
from skspatial.objects import Point, Vector, Circle

from Geometry.particle import Particle
from Processes.point_process import PoissonPointProcess
from Processes.particle_process import ParticleProcess
from Processes.segment_process import SegmentProcess
from Processes.ball_process import BallProcess
from Geometry.grain import Segment
import utils.const as const


def execute_boolean_particle_process(
        grain_type: str, intensity: float,
        min_grain_span: float, max_grain_span: float,
        window_edge_length: float = 1,
        **kwargs
):
    # we need to set up a window length for a ground process
    window_ground_edge_start_point, window_ground_edge_end_point = set_the_window_length(
        window_edge_length=window_edge_length, max_grain_span=max_grain_span
    )
    # simulate Poisson Point Process (germ process)
    germ_process = PoissonPointProcess(
        intensity=intensity,
        window_edge_start_point=window_ground_edge_start_point, window_edge_end_point=window_ground_edge_end_point
    )


def set_the_window_length(window_edge_length: float, max_grain_span: float):
    return - max_grain_span, window_edge_length + max_grain_span
