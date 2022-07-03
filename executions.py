
from datetime import datetime
import logging
import numpy as np
from skspatial.objects import Point, Vector, Circle

from Geometry.particle import Particle
from Processes.point_process import PoissonPointProcess
from Processes.particle_process import ParticleProcess
from Processes.segment_process import SegmentProcess
from Processes.ball_process import BallProcess
from Processes.markings import Mark
from Geometry.grain import Segment
from utils.config_parser import ConfigParser
import utils.const as const
from utils.results_saver import ResultSaver
from tests.wider_window_simulation_tests import simulate_the_processes


def overnight_computations():
    result_savers = []
    for seed in range(200, 220):
        print(f"Seed no {seed}: {datetime.now()}")
        result_saver = simulate_the_processes(alphas_list=[0, 0.25, 0.5, 0.75, 1], seed=seed)
        result_savers.append(result_saver)
    return result_savers


def execute_from_config(config_parser: ConfigParser):
    result_savers = []
    for seed in range(config_parser.initial_seed, config_parser.initial_seed + config_parser.number_of_realizations):
        config_parser.initialize_the_processes(seed=seed)
        result_saver = config_parser.return_the_result_saver(seed=seed)
        a=1


# testing null segment & ball processes (Bernoulli) &
# angles/radii marks (either as parameter to Bernoulli or simply its value)
def first_blood(number_of_seeds: int = 1):

    result_saver = ResultSaver()

    for seed in range(number_of_seeds):

        print(f"{datetime.now()} --- Seed no: {seed}")

        np.random.seed(seed=seed)

        # one underlying Poisson Point Process for each model
        poisson_point_process = PoissonPointProcess(intensity=const.POISSON_INTENSITY)

        segment_angles = [np.pi * np.random.random_sample() for _ in range(len(poisson_point_process.points))]

        segment_lengths = [
            const.MIN_SEGMENT_LENGTH + (const.MAX_SEGMENT_LENGTH - const.MIN_SEGMENT_LENGTH) * np.random.random_sample()
            for _ in range(len(poisson_point_process.points))
        ]

        ball_radii = [
            const.MIN_CIRC_RAD + (const.MAX_CIRC_RAD - const.MIN_CIRC_RAD) * np.random.random_sample()
            for _ in range(len(poisson_point_process.points))
        ]

        marks_null = np.random.binomial(size=len(poisson_point_process.points), n=1, p=0.5)
        marks_ball_discrete = [
            np.random.binomial(size=1, n=1, p=np.array(1 / const.MAX_CIRC_RAD * ball_radii[k]))[0]
            for k in range(len(poisson_point_process.points))
        ]
        marks_angles_beatles_discrete = [
            np.random.binomial(size=1, n=1, p=np.array(1 / np.pi * segment_angles[k]))[0]
            for k in range(len(poisson_point_process.points))
        ]

        particles_balls_null = [
            Particle(
                germ=Point(poisson_point_process.points[k]),
                grain_type="ball",
                grain=Circle(point=Point(poisson_point_process.points[k]), radius=ball_radii[k]),
                mark=Mark(mark_type="discrete", mark_value=marks_null[k])
            )
            for k in range(len(poisson_point_process.points))
        ]

        ball_process_null = BallProcess(
            particles=particles_balls_null, germ_intensity=poisson_point_process.intensity, marked=True
        )
        ball_process_null.compute_the_f_mark_characteristics()
        ball_process_null.perform_the_permutation_test_for_f_mark_characteristics()
        result_saver.save_the_results(
            model_name="null model", value_dict=ball_process_null.f_mark_statistics,
            quantile_dict=ball_process_null.f_mark_statistics_quantiles,
            permutations_count=const.PERMUTATION_TEST_REPEAT_COUNT,
            grain_type="ball", seed=seed
        )

        particles_balls_radii_continuous = [
            Particle(
                germ=Point(poisson_point_process.points[k]),
                grain_type="ball",
                grain=Circle(point=Point(poisson_point_process.points[k]), radius=ball_radii[k]),
                mark=Mark(mark_type="discrete", mark_value=ball_radii[k])
            )
            for k in range(len(poisson_point_process.points))
        ]

        ball_process_radii_continuous = BallProcess(
            particles=particles_balls_radii_continuous, germ_intensity=poisson_point_process.intensity, marked=True
        )
        ball_process_radii_continuous.compute_the_f_mark_characteristics()
        ball_process_radii_continuous.perform_the_permutation_test_for_f_mark_characteristics()
        result_saver.save_the_results(
            model_name="continuous radii", value_dict=ball_process_radii_continuous.f_mark_statistics,
            quantile_dict=ball_process_radii_continuous.f_mark_statistics_quantiles,
            permutations_count=const.PERMUTATION_TEST_REPEAT_COUNT,
            grain_type="ball", seed=seed
        )

        particles_balls_radii_discrete = [
            Particle(
                germ=Point(poisson_point_process.points[k]),
                grain_type="ball",
                grain=Circle(point=Point(poisson_point_process.points[k]), radius=ball_radii[k]),
                mark=Mark(mark_type="discrete", mark_value=marks_ball_discrete[k])
            )
            for k in range(len(poisson_point_process.points))
        ]

        balls_process_radii_discrete = BallProcess(
            particles=particles_balls_radii_discrete, germ_intensity=poisson_point_process.intensity, marked=True
        )
        balls_process_radii_discrete.compute_the_f_mark_characteristics()
        balls_process_radii_discrete.perform_the_permutation_test_for_f_mark_characteristics()
        result_saver.save_the_results(
            model_name="discrete radius", value_dict=balls_process_radii_discrete.f_mark_statistics,
            quantile_dict=balls_process_radii_discrete.f_mark_statistics_quantiles,
            permutations_count=const.PERMUTATION_TEST_REPEAT_COUNT,
            grain_type="ball", seed=seed
        )

        # SEGMENTS PART
        particles_segments_null = [
            Particle(
                grain_type="segment",
                germ=Point(poisson_point_process.points[k]),
                grain=Segment(
                    start_point=poisson_point_process.points[k] - 1 / 2 * segment_lengths[k] * np.array(
                        [np.cos(segment_angles[k]), np.sin(segment_angles[k])]
                    ),
                    angle=segment_angles[k],
                    length=segment_lengths[k]
                ),
                mark=Mark(mark_type="discrete", mark_value=marks_null[k])
            )
            for k in range(len(poisson_point_process.points))
        ]

        # SEGMENTS
        segment_process_null = SegmentProcess(
            particles=particles_segments_null, germ_intensity=poisson_point_process.intensity, marked=True
        )
        segment_process_null.compute_the_f_mark_characteristics()
        segment_process_null.perform_the_permutation_test_for_f_mark_characteristics()
        result_saver.save_the_results(
            model_name="null model", value_dict=segment_process_null.f_mark_statistics,
            quantile_dict=segment_process_null.f_mark_statistics_quantiles,
            permutations_count=const.PERMUTATION_TEST_REPEAT_COUNT,
            grain_type="segment", seed=seed
        )

        particles_segments_angle_discrete = [
            Particle(
                grain_type="segment",
                germ=Point(poisson_point_process.points[k]),
                grain=Segment(
                    start_point=poisson_point_process.points[k] - 1 / 2 * segment_lengths[k] * np.array(
                        [np.cos(segment_angles[k]), np.sin(segment_angles[k])]
                    ),
                    angle=segment_angles[k],
                    length=segment_lengths[k]
                ),
                mark=Mark(mark_type="discrete", mark_value=marks_angles_beatles_discrete[k])
            )
            for k in range(len(poisson_point_process.points))
        ]
        segment_process_angles_discrete = SegmentProcess(
            particles=particles_segments_angle_discrete, germ_intensity=poisson_point_process.intensity, marked=True
        )
        segment_process_angles_discrete.compute_the_f_mark_characteristics()
        segment_process_angles_discrete.perform_the_permutation_test_for_f_mark_characteristics()
        result_saver.save_the_results(
            model_name="discrete angle", value_dict=segment_process_angles_discrete.f_mark_statistics,
            quantile_dict=segment_process_angles_discrete.f_mark_statistics_quantiles,
            permutations_count=const.PERMUTATION_TEST_REPEAT_COUNT,
            grain_type="segment", seed=seed
        )

        particles_segments_angle_continuous = [
            Particle(
                grain_type="segment",
                germ=Point(poisson_point_process.points[k]),
                grain=Segment(
                    start_point=poisson_point_process.points[k] - 1 / 2 * segment_lengths[k] * np.array(
                        [np.cos(segment_angles[k]), np.sin(segment_angles[k])]
                    ),
                    angle=segment_angles[k],
                    length=segment_lengths[k]
                ),
                mark=Mark(mark_type="discrete", mark_value=segment_angles[k] / np.pi)
            )
            for k in range(len(poisson_point_process.points))
        ]
        segment_process_angles_continuous = SegmentProcess(
            particles=particles_segments_angle_continuous, germ_intensity=poisson_point_process.intensity, marked=True
        )
        segment_process_angles_continuous.compute_the_f_mark_characteristics()
        segment_process_angles_continuous.perform_the_permutation_test_for_f_mark_characteristics()
        result_saver.save_the_results(
            model_name="continuous angle", value_dict=segment_process_angles_continuous.f_mark_statistics,
            quantile_dict=segment_process_angles_continuous.f_mark_statistics_quantiles,
            permutations_count=const.PERMUTATION_TEST_REPEAT_COUNT,
            grain_type="segment", seed=seed
        )
    result_saver.save_to_pandas()
    if const.PICKLE_RESULTS:
        result_saver.pickle_the_result_dataframes()
    return result_saver


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
