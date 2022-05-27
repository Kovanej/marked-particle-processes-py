
from datetime import datetime
import logging
import matplotlib.pyplot as plt
from typing import List, Optional, Union, Tuple, Dict
import numpy as np
import random
from scipy.optimize import fsolve
from sklearn.metrics import pairwise_distances
from skspatial.objects import Point, Vector

from Geometry.particle import Particle
import utils.const as const


class ParticleProcess(object):

    def __init__(
            self,
            germ_intensity: float,
            particles: List[Particle],
            grain_type: str,
            space_dimension: int,
            marked: bool,
            marks_aposteriori_type: Optional[str] = None
    ):
        logging.info(f"Particle process class initialized.")
        self.germ_intensity = germ_intensity
        self.particles = particles
        self.number_of_particles = len(self.particles)
        self.grain_type = grain_type
        self.space_dimension = space_dimension
        self.marked = marked
        if self.marked:
            self.marks = np.array([p.mark.mark_value for p in self.particles])
            self.marks_product = self.marks[..., None] * self.marks[None, ...]
            self.marks_square = (self.marks[..., None] - self.marks[None, ...]) ** 2 / 2
            self.max_mark = np.amax(self.marks)
            self.min_mark = np.amin(self.marks)
        self.particle_measure = self._compute_the_particles_measure()
        # compute the grains distance
        self.germs_distance_matrix = self._compute_the_germs_distance_matrix()
        # compute particles_null_model distance (inf {||x_i - x_j||: x_i \in \Xi_i, x_j \in \Xi_j})
        self.particles_distance_matrix = self._compute_the_particles_distance_matrix()
        # particles_null_model distance == 0 -> intersected, otherwise not intersected
        self.particles_intersection_matrix = self._compute_the_particles_intersection_matrix()
        # compute the pairwise shared corresponding Lebesgue measure
        # ("ball": shared areas, "segment": same as intersection matrix ...)
        self.pairwise_shared_measure_matrix = self._compute_the_pairwise_shared_measure_matrix()
        # if needed following attributes are computed via executing ParticleProcess.compute_f_mark_statistics
        self.f_mark_normalization_constants: Dict[str, float] = {}
        self.f_mark_statistics: Dict[Tuple[str, str], float] = {}
        self.f_mark_statistics_permutations: Dict[Tuple[str, str], np.array] = {
            (f_type, weight_type): np.array([]) for (f_type, weight_type) in const.F_MARK_COMBINATIONS
        }
        self.f_mark_statistics_quantiles: Dict[Tuple[str, str], float] = {}

    @staticmethod
    def _dot_pairwise(a, b):
        return (a[:, None, :] * b[None, ...]).sum(axis=-1)

    @staticmethod
    def _norm(a):
        return np.sqrt((a * a).sum(axis=-1))

    def _compute_the_particles_intersection_matrix(self):
        logging.info(f"{datetime.now()} :Particle intersection matrix computation start.")
        intersection_matrix = np.where(self.particles_distance_matrix == 0, 1, 0)
        logging.info(f"{datetime.now()} :Particle intersection matrix computation end.")
        return intersection_matrix

    def _compute_the_pairwise_shared_measure_matrix(self):
        # overridden for each subclass - if not specified, it cannot be computed
        raise NotImplementedError(
            f"Method called from the instance of general ParticleProcess class. "
            f"Please use the instance of some of the ParticleProcess subclasses (SegmentProcess, BallProcess)."
        )

    def _compute_the_germs_distance_matrix(self):
        logging.info(f"{datetime.now()} :Germs distance computation start.")
        grains_distance_matrix = pairwise_distances(
                [self.particles[k].germ for k in range(self.number_of_particles)]
            )
        logging.info(f"{datetime.now()} :Germs distance computation end.")
        return grains_distance_matrix

    def compute_the_f_mark_characteristics(
            self
    ):
        for f_type in const.F_MARK_TYPES:
            self._compute_the_f_mark_normalization_constants(f_type=f_type)
        for (f_type, weight_type) in const.F_MARK_COMBINATIONS:
            weight_matrix = {
                "intersection": self.particles_intersection_matrix,
                "shared_area": self.pairwise_shared_measure_matrix,
                "distance": self.particles_distance_matrix
            }.get(weight_type)
            self.f_mark_statistics[f_type, weight_type] = self._compute_the_f_mark_correlation(
                f_type=f_type, weight_matrix=weight_matrix
            )

    def perform_the_permutation_test_for_f_mark_characteristics(self):
        for _ in range(const.PERMUTATION_TEST_REPEAT_COUNT):
            permutation = np.random.permutation(self.number_of_particles)
            marks_permuted = self.marks[permutation]
            for (f_type, weight_type) in const.F_MARK_COMBINATIONS:
                weight_matrix = {
                    "intersection": self.particles_intersection_matrix,
                    "shared_area": self.pairwise_shared_measure_matrix,
                    "distance": self.particles_distance_matrix
                }.get(weight_type)
                val = self._compute_the_f_mark_correlation(
                        f_type=f_type, weight_matrix=weight_matrix, marks_vector=marks_permuted
                    )
                self.f_mark_statistics_permutations[f_type, weight_type] = np.append(
                    self.f_mark_statistics_permutations[f_type, weight_type], val
                )
        for (f_type, weight_type) in const.F_MARK_COMBINATIONS:
            self.f_mark_statistics_quantiles[f_type, weight_type] = np.where(
                self.f_mark_statistics[f_type, weight_type] <= self.f_mark_statistics_permutations[f_type, weight_type],
                0, 1
            ).mean()

    def _compute_the_f_mark_correlation(
            self, f_type: str, weight_matrix: np.array , marks_vector: Optional[np.array] = None
    ):
        # since we are working on a [0, 1]^d window -> not normed by its Lebesgue measure, since it is 1 (je to jedno)
        weight_matrix_zero_diagonal = weight_matrix.copy()
        np.fill_diagonal(weight_matrix_zero_diagonal, 0)
        if marks_vector is None:
            marks_vector = self.marks
            marks_product = self.marks_product
            marks_square = self.marks_square
        else:
            marks_product = marks_vector[..., None] * marks_vector[None, ...]
            marks_square = (marks_vector[..., None] - marks_vector[None, ...]) ** 2 / 2
        if f_type == "product":
            f_nn = (marks_product * weight_matrix_zero_diagonal).sum()
            norm_by = self.f_mark_normalization_constants[f_type] * (self.germ_intensity ** 2)
        elif f_type == "square":
            f_nn = (marks_square * weight_matrix_zero_diagonal).sum()
            norm_by = self.f_mark_normalization_constants[f_type] * (self.germ_intensity ** 2)
        elif f_type == "first_mark":
            # TODO in my opinion we should not use the np.triu-ed matrix, only the zero-diagonal one
            shared_area_zero_on_and_under_diagonal = np.triu(weight_matrix_zero_diagonal)
            f_nn = (marks_vector * weight_matrix_zero_diagonal).sum()
            norm_by = self.f_mark_normalization_constants[f_type] * (self.germ_intensity ** 2)
        else:
            raise NotImplementedError(f"f-mark shared area correlation cannot be obtained for unknown f_type={f_type}")
        # TODO for testing right now returning non-normed values
        # return f_shared_nn / norm_by
        return f_nn

    def _compute_the_f_mark_normalization_constants(self, f_type: str):
        n_times_n_minus_one = (self.number_of_particles * (self.number_of_particles - 1))
        if f_type == "product":
            marks_product_zero_diagonal = self.marks_product.copy()
            np.fill_diagonal(marks_product_zero_diagonal, 0)
            c_f_nn = marks_product_zero_diagonal.sum()
            c_f = c_f_nn / n_times_n_minus_one
        elif f_type == "square":
            # zero_diagonal by default - no need to "hardcode"
            c_f_nn = self.marks_square.sum()
            c_f = c_f_nn / n_times_n_minus_one
        elif f_type == "first_mark":
            c_f_nn = self.marks.sum()
            c_f = c_f_nn / self.number_of_particles
        else:
            raise NotImplementedError(f"f-mark normalization constant cannot be obtained for unknown f_type={f_type}")
        self.f_mark_normalization_constants[f_type] = c_f

    def _plot_particles(self, ax, fig):
        pass

    def plot_itself(self, show_germs: bool = False):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_aspect('equal', adjustable='box')
        self._plot_particles(ax=ax, fig=fig)
        if show_germs:
            for particle in self.particles:
                color, alpha = self._choose_germ_color()
                particle.germ.plot_2d(ax, c=color, alpha=alpha)
        if const.SAVE_PLOTS:
            plt.savefig(f"generated_pics/{str(datetime.now()).replace(':','-')}_plot.png", dpi=1000)
        plt.show()

    def _compute_the_particles_distance_matrix(self):
        raise NotImplementedError(
            f"Method called from the instance of general ParticleProcess class. "
            f"Please use the instance of some of the ParticleProcess subclasses (SegmentProcess, BallProcess)."
        )

    def _choose_edge_color(self):
        return "#000000"
        # return np.random.choice(const.PARTICLE_COLORS_CHOICE)

    def _choose_face_color(self, particle=None):
        if particle.mark is not None:
            alpha = (1 + particle.mark.mark_value - self.min_mark) / (1 + self.max_mark - self.min_mark)
        else:
            alpha = const.ALPHA
        col = np.random.choice(const.PARTICLE_COLORS_CHOICE)
        return col, alpha

    def _choose_germ_color(self):
        return "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]), np.random.random(1)[0]

    def _compute_the_particles_measure(self):
        pass

