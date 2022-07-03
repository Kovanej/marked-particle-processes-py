
import os
import numpy as np
from skspatial.objects import Point, Circle
from typing import List

from Geometry.particle import Particle, Segment
from Processes.markings import Mark
from Processes.point_process import PoissonPointProcess
from Processes.segment_process import SegmentProcess
from Processes.ball_process import BallProcess
import utils.const as const
from utils.results_saver import ResultSaver


def simulate_the_processes(alphas_list: List[float] = [0, 0.5,  1], seed: int = 1):

    poisson_point_process = PoissonPointProcess(
        intensity=const.POISSON_INTENSITY,
        window_edge_start_point=-const.MAX_CIRC_RAD,
        window_edge_end_point=1 + const.MAX_CIRC_RAD
    )

    result_saver = ResultSaver()

    for alpha in alphas_list:
        # radius & angles markings
        ball_radius_categorical_particles = []
        segment_angle_categorical_particles = []
        segment_length_categorical_particles = []
        ball_radius_continuous_particles = []
        segment_angle_continuous_particles = []
        segment_length_continuous_particles = []
        for k in range(len(poisson_point_process.points)):
            radius = np.random.random_sample() * (const.MAX_CIRC_RAD - const.MIN_CIRC_RAD) + const.MAX_CIRC_RAD
            length = np.random.random_sample() * (const.MAX_SEGMENT_LENGTH - const.MIN_SEGMENT_LENGTH) + const.MIN_SEGMENT_LENGTH
            angle = np.random.random_sample() * (const.MAX_SEGMENT_ANGLE - const.MIN_SEGMENT_ANGLE) + const.MIN_SEGMENT_ANGLE
            bernoulli_independent = np.random.binomial(n=1, p=1/2, size=1)[0]
            bernoulli_radius = np.random.binomial(n=1, p=(radius - const.MIN_CIRC_RAD)/const.MAX_CIRC_RAD, size=1)[0]
            bernoulli_length = np.random.binomial(n=1, p=(length - const.MIN_SEGMENT_LENGTH)/const.MAX_SEGMENT_LENGTH, size=1)[0]
            bernoulli_angle = np.random.binomial(n=1, p=angle/np.pi, size=1)[0]
            bernoulli_alpha = np.random.binomial(n=1, p=alpha, size=1)[0]
            mark_radius_categorical = bernoulli_alpha * bernoulli_radius + (1 - bernoulli_alpha) * bernoulli_independent
            mark_length_categorical = bernoulli_alpha * bernoulli_length + (1 - bernoulli_alpha) * bernoulli_independent
            mark_angle_categorical = bernoulli_alpha * bernoulli_angle + (1 - bernoulli_alpha) * bernoulli_independent
            ball_radius_categorical_particles.append(
                Particle(
                    germ=Point(poisson_point_process.points[k]),
                    grain_type="ball",
                    grain=Circle(point=Point(poisson_point_process.points[k]), radius=radius),
                    mark=Mark(mark_type="discrete", mark_value=mark_radius_categorical)
                )
            )
            segment_length_categorical_particles.append(
                Particle(
                    germ=Point(poisson_point_process.points[k]),
                    grain_type="segment",
                    grain=Segment(
                        start_point=poisson_point_process.points[k] - 1 / 2 * length * np.array(
                            [np.cos(angle), np.sin(angle)]),
                        length=length, angle=angle
                    ),
                    mark=Mark(mark_type="discrete", mark_value=mark_length_categorical)
                )
            )
            segment_angle_categorical_particles.append(
                Particle(
                    germ=Point(poisson_point_process.points[k]),
                    grain_type="segment",
                    grain=Segment(
                        start_point=poisson_point_process.points[k] - 1 / 2 * length * np.array(
                            [np.cos(angle), np.sin(angle)]),
                        length=length, angle=angle
                    ),
                    mark=Mark(mark_type="discrete", mark_value=mark_angle_categorical)
                )
            )
        ball_process_radius_categorical = BallProcess(
            germ_intensity=poisson_point_process.intensity, particles=ball_radius_categorical_particles, marked=True
        )
        segment_process_length_categorical = SegmentProcess(
            germ_intensity=poisson_point_process.intensity, particles=segment_length_categorical_particles, marked=True
        )
        segment_process_angle_categorical = SegmentProcess(
            germ_intensity=poisson_point_process.intensity, particles=segment_angle_categorical_particles, marked=True
        )
        processes = [
            ball_process_radius_categorical, segment_process_length_categorical, segment_process_angle_categorical
        ]
        for k in range(len(processes)):
            processes[k].compute_the_f_mark_characteristics()
            processes[k].perform_the_permutation_test_for_f_mark_characteristics()
            result_saver.save_the_results(
                model_name=f"{processes[k].grain_type}_{k}_alpha={alpha}",
                value_dict=processes[k].f_mark_statistics,
                quantile_dict=processes[k].f_mark_statistics_quantiles,
                permutations_count=const.PERMUTATION_TEST_REPEAT_COUNT,
                grain_type=f"{processes[k].grain_type}", seed=seed
            )
    result_saver.save_to_pandas(save_csv=True)
    result_saver.pickle_the_result_dataframes()
    return result_saver


def simulate_the_processes_legacy(alphas_list: List[float] = [0, 0.5,  1]):
    poisson_point_process = PoissonPointProcess(
        intensity=const.POISSON_INTENSITY,
        window_edge_start_point=-const.MAX_CIRC_RAD,
        window_edge_end_point=1 + const.MAX_CIRC_RAD
    )

    dict_of_dependent_processes = {
        "ball_processes": {alpha: [] for alpha in alphas_list},
        "segment_processes": {alpha: [] for alpha in alphas_list}
    }

    # testing for balls & segments
    particles_ball_null_model = []
    particles_segment_null_model = []
    particles_ball_unmarked_model = []
    particles_segment_unmarked_model = []
    for k in range(len(poisson_point_process.points)):
        # balls
        technically_no_mark = Mark(mark_type="discrete", mark_value=9)
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
        particles_ball_unmarked_model.append(
            Particle(
                germ=Point(poisson_point_process.points[k]),
                grain_type="ball",
                grain=Circle(point=Point(poisson_point_process.points[k]), radius=radius),
                mark=technically_no_mark
            )
        )
        length = np.random.random_sample() * (const.MAX_SEGMENT_LENGTH - const.MIN_SEGMENT_LENGTH) + const.MIN_SEGMENT_LENGTH
        angle = np.random.random_sample() * (const.MAX_SEGMENT_ANGLE - const.MIN_SEGMENT_ANGLE) + const.MIN_SEGMENT_ANGLE
        particles_segment_null_model.append(
            Particle(
                germ=Point(
                    poisson_point_process.points[k]
                ),
                grain_type="segment",
                grain=Segment(
                    start_point=poisson_point_process.points[k] - 1/2 * length * np.array([np.cos(angle), np.sin(angle)]),
                    length=length, angle=angle
                ),
                mark=mark_null_model
            )
        )
        particles_segment_unmarked_model.append(
            Particle(
                germ=Point(
                    poisson_point_process.points[k]
                ),
                grain_type="segment",
                grain=Segment(
                    start_point=poisson_point_process.points[k] - 1/2 * length * np.array([np.cos(angle), np.sin(angle)]),
                    length=length, angle=angle
                ),
                mark=technically_no_mark
            )
        )

    ball_process_unmarked_model = BallProcess(
        germ_intensity=poisson_point_process.intensity, particles=particles_ball_unmarked_model, marked=True
    )
    ball_process_unmarked_model.plot_itself(show_germs=False)

    segment_process_unmarked_model = SegmentProcess(
        germ_intensity=poisson_point_process.intensity, particles=particles_segment_unmarked_model
    )
    segment_process_unmarked_model.plot_itself(show_germs=False)

    ball_process_null_model = BallProcess(
        germ_intensity=poisson_point_process.intensity, particles=particles_ball_null_model, marked=True
    )
    ball_process_null_model.plot_itself()

    segment_process_null_model = SegmentProcess(
        germ_intensity=poisson_point_process.intensity, particles=particles_segment_null_model
    )
    segment_process_null_model.plot_itself()
