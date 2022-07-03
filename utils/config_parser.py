
from typing import Dict, List, Optional, Tuple

from datetime import datetime
import logging
import numpy as np
from skspatial.objects import Point, Vector, Circle

from Geometry.particle import Particle
from Processes.point_process import PointProcess, PoissonPointProcess
from Processes.particle_process import ParticleProcess
from Processes.segment_process import SegmentProcess
from Processes.ball_process import BallProcess
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
        self.list_of_processes_per_seed_and_name: Dict[Tuple[str, int], List[ParticleProcess]] = {}

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

    def initialize_the_corresponding_process(self, seed: int = 23) -> List[ParticleProcess]:
        np.random.seed(seed=seed)
        win_edge_start_point, win_edge_end_point = self._set_the_window_length()
        poisson_point_process = PoissonPointProcess(
            intensity=self.intensity, window_edge_start_point=win_edge_start_point,
            window_edge_end_point=win_edge_end_point
        )
        self.germ_processes_per_seed[seed] = poisson_point_process
        if self.process_type == "ball":
            processes = self._initialize_the_ball_processes()
        elif self.process_type == "segment":
            processes = self._initialize_the_segment_processes()
        else:
            raise ValueError(f"Incapable of simulating process of unknown process_type: {self.process_type}")
        return processes

    def _set_the_window_length(self):
        if self.process_type == "ball":
            max_overlap = self.particles_parameters["ball"]["max_radius"]
        elif self.process_type == "segment":
            max_overlap = self.particles_parameters["segment"]["max_segment_length"]
        else:
            raise ValueError(f"Unknown particle process type: {self.process_type}")
        return - max_overlap, 1 + max_overlap

    def _initialize_the_ball_processes(self):
        return []

    def _initialize_the_segment_processes(self):
        return []

    def _set_the_attribute_as_none_values_so_that_py_charm_shows_me_their_name(self):
        pass

