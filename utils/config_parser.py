
from typing import Dict, List, Optional, Tuple

import copy
from datetime import datetime
import logging
import numpy as np
from skspatial.objects import Point, Vector, Circle

from Geometry.particle import Particle
from Processes.point_process import PointProcess, PoissonPointProcess
from Processes.particle_process import ParticleProcess
import Processes.segment_process as sp
import Processes.ball_process as bp
from Processes.markings import Mark
from Geometry.grain import Segment
import utils.const as const
from utils.results_saver import ResultSaver
from tests.wider_window_simulation_tests import simulate_the_processes
from utils.const import CONFIG_VALID_KEYS, CONFIG_NON_NULLABLE, CONFIG_OPTIONAL_VALUES


class ConfigParser(object):

    def __init__(self, config: Dict):
        self._set_the_attribute_as_none_values_so_that_py_charm_shows_me_their_name()
        self._config = config
        self._validate_the_config()
        self._parse_the_config()
        self.germ_processes_per_seed: Dict[int, PointProcess] = {}
        self.lists_of_processes_per_seed_and_name: Dict[Tuple[str, int], List[ParticleProcess]] = {}

    def _parse_the_config(self):
        self._parse_the_non_nullable_values()
        self._parse_the_nullable_values()

    def _parse_the_non_nullable_values(self):
        for _ in CONFIG_NON_NULLABLE:
            setattr(self, _, self._config[_])

    def _parse_the_nullable_values(self):
        for attr_name, attr_value in CONFIG_OPTIONAL_VALUES.items():
            if attr_name in self._config.keys():
                setattr(self, attr_name, self._config[attr_name])
            else:
                setattr(self, attr_name, attr_value)

    def _validate_the_config(self):
        _invalid_keys = [
            _ for _ in self._config.keys() if _ not in CONFIG_VALID_KEYS
        ]
        if len(_invalid_keys) > 0:
            self.valid_config = False
            raise ValueError(
                f"Unknown key(s) in config.json: {_invalid_keys}"
            )
        _missing = [
            _ for _ in CONFIG_NON_NULLABLE if (_ not in self._config.keys()) or (self._config[_] is None)
        ]
        if len(_missing) > 0:
            self.valid_config = False
            raise ValueError(
                f"Following non-nullable parameters either missing in config or have assigned null value: {_missing}"
            )
        self.valid_config = True

    def return_the_result_saver(self, seed: int) -> ResultSaver:
        result_saver = ResultSaver()
        processes_to_save = {key: v for key, v in self.lists_of_processes_per_seed_and_name.items() if key[1] == seed}
        for key, process in processes_to_save.items():
            if self.plot_realizations:
                process.plot_itself()
            process.compute_the_f_mark_characteristics()
            process.perform_the_permutation_test_for_f_mark_characteristics()
            for weight, fs in self.f_mark_weights_and_statistics.items():
                for f in fs:
                    result_saver.save_the_results(
                        model_name=key[0], grain_type=process.grain_type,
                        permutations_count=const.PERMUTATION_TEST_REPEAT_COUNT,
                        quantile_dict=process.f_mark_statistics_quantiles, value_dict=process.f_mark_statistics,
                        seed=seed, intensity=process.germ_intensity
                    )
        result_saver.save_to_pandas(save_csv=self.save_results)
        return result_saver

    def initialize_the_processes(self, seed: int = 23) -> None:
        np.random.seed(seed=seed)
        win_edge_start_point, win_edge_end_point = self._set_the_window_length()
        poisson_point_process = PoissonPointProcess(
            intensity=self.intensity, window_edge_start_point=win_edge_start_point,
            window_edge_end_point=win_edge_end_point
        )
        self.germ_processes_per_seed[seed] = poisson_point_process
        if self.process_type == "ball":
            processes = self._initialize_the_ball_processes(seed=seed)
        elif self.process_type == "segment":
            processes = self._initialize_the_segment_processes(seed=seed)
        else:
            raise ValueError(f"Incapable of simulating process of unknown process_type: {self.process_type}")
        for process in processes:
            self.lists_of_processes_per_seed_and_name[(process.model_name, seed)] = process

    def _set_the_window_length(self) -> Tuple[int, int]:
        if self.process_type == "ball":
            max_overlap = self.particles_parameters["ball"]["max_radius"]
        elif self.process_type == "segment":
            max_overlap = self.particles_parameters["segment"]["max_segment_length"]
        else:
            raise ValueError(f"Unknown particle process type: {self.process_type}")
        return - max_overlap, 1 + max_overlap

    def _initialize_the_ball_processes(self, seed: int) -> List[ParticleProcess]:
        max_rad = self.particles_parameters["ball"]["max_radius"]
        min_rad = self.particles_parameters["ball"]["min_radius"]
        particle_processes = []
        particles = []
        for point in self.germ_processes_per_seed[seed].points:
            radius = min_rad + (max_rad - min_rad) * np.random.random_sample()
            grain = Circle(point=point, radius=radius)
            particle = Particle(germ=point, grain=grain)
            particles.append(particle)
        for model in self.marking_type["ball"]:
            for alpha in self.marking_parameters["alphas"]:
                particles = copy.deepcopy(particles)
                if model == "radius_discrete":
                    particle_process = bp.BivariateRadiusMarksBallProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity,
                        particles=particles, alpha=alpha, max_radius=max_rad, min_radius=min_rad, seed=seed
                    )
                elif model == "radius_continuous":
                    particle_process = bp.ContinuousRadiusMarksBallProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity,
                        particles=particles, alpha=alpha, max_radius=max_rad, min_radius=min_rad, seed=seed
                    )
                elif model == "max_shared_area_discrete":
                    particle_process = bp.BivariateMaximalSharedAreaMarkBallProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity, particles=particles,
                        alpha=alpha, min_radius=min_rad, max_radius=max_rad
                    )
                elif model == "max_shared_area_continuous":
                    particle_process = bp.ContinuousMaximalSharedAreaMarkBallProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity, particles=particles,
                        alpha=alpha, min_radius=min_rad, max_radius=max_rad
                    )
                elif model == "nearest_neighbour_distance":
                    particle_process = bp.ContinuousNNDistanceMarkBallProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity, particles=particles,
                        alpha=alpha, min_radius=min_rad, max_radius=max_rad
                    )
                elif model == "intersection_counting":
                    particle_process = bp.CountingIntersectionNumberMarkBallProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity, particles=particles,
                        alpha=alpha, min_radius=min_rad, max_radius=max_rad
                    )
                else:
                    raise NotImplementedError(f"Unknown model '{model}' for a Ball Process.")
                particle_processes.append(particle_process)
        return particle_processes

    def _initialize_the_segment_processes(self, seed: int) -> List[ParticleProcess]:
        max_len = self.particles_parameters["segment"]["max_segment_length"]
        min_len = self.particles_parameters["segment"]["min_segment_length"]
        max_angle_rad = self.particles_parameters["segment"]["max_angle_in_degrees"] * np.pi / 180
        min_angle_rad = self.particles_parameters["segment"]["min_angle_in_degrees"] * np.pi / 180
        particle_processes = []
        particles = []
        for point in self.germ_processes_per_seed[seed].points:
            length = min_len + (max_len - min_len) * np.random.random_sample()
            angle = min_angle_rad + (max_angle_rad - min_angle_rad) * np.random.random_sample()
            particle = Particle(
                germ=point,
                grain_type="segment",
                grain=Segment(
                    start_point=point - 1 / 2 * length * np.array(
                        [np.cos(angle), np.sin(angle)]),
                    length=length, angle=angle
                )
            )
            particles.append(particle)
        for model in self.marking_type["segment"]:
            print(model)
            for alpha in self.marking_parameters["alphas"]:
                particles = copy.deepcopy(particles)
                if model == "angle_discrete":
                    particle_process = sp.BivariateAngleMarksSegmentProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity,
                        particles=particles, alpha=alpha, max_angle=max_angle_rad, min_angle=min_angle_rad, seed=seed,
                        max_length=max_len, min_length=min_len
                    )
                elif model == "angle_continuous":
                    particle_process = sp.ContinuousAngleMarksSegmentProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity,
                        particles=particles, alpha=alpha, max_angle=max_angle_rad, min_angle=min_angle_rad, seed=seed,
                        max_length=max_len, min_length=min_len
                    )
                elif model == "length_discrete":
                    particle_process = sp.BivariateLengthMarksSegmentProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity,
                        particles=particles, alpha=alpha, max_length=max_len, min_length=min_len, seed=seed,
                        max_angle=max_angle_rad, min_angle=min_angle_rad
                    )
                elif model == "length_continuous":
                    particle_process = sp.ContinuousLengthMarksSegmentProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity,
                        particles=particles, alpha=alpha, max_length=max_len, min_length=min_len, seed=seed,
                        max_angle=max_angle_rad, min_angle=min_angle_rad
                    )
                elif model == "nearest_neighbour_distance":
                    particle_process = sp.ContinuousNNDistanceMarkSegmentProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity,
                        particles=particles, alpha=alpha, max_length=max_len, min_length=min_len, seed=seed,
                        max_angle=max_angle_rad, min_angle=min_angle_rad
                    )
                elif model == "intersection_counting":
                    particle_process = sp.CountingIntersectionNumberMarkSegmentProcess(
                        germ_intensity=self.germ_processes_per_seed[seed].intensity,
                        particles=particles, alpha=alpha, max_length=max_len, min_length=min_len, seed=seed,
                        max_angle=max_angle_rad, min_angle=min_angle_rad
                    )
                else:
                    raise NotImplementedError(f"Unknown model '{model}' for a Segment Process.")
                particle_processes.append(particle_process)
        return particle_processes

    def _set_the_attribute_as_none_values_so_that_py_charm_shows_me_their_name(self):
        self.process_type = None
        self.intensity = None
        self.space_dimension = None
        self.marking_type = None
        self.particles_parameters = None
        self.marking_parameters = None
        self.plot_realizations = None
        self.compute_f_mark_statistics = None
        self.f_mark_weights_and_statistics = None
        self.perform_permutation_test = None
        self.permutation_tests_parameters = None
        self.initial_seed = None
        self.number_of_realizations = None
        self.save_results = None


