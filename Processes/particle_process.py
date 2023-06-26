
from datetime import datetime
import logging
from matplotlib.colors import ListedColormap
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from typing import List, Optional, Union, Tuple, Dict
import numpy as np
import random
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
            model_name: Optional[str] = None,
            marked_aposteriori: Optional[bool] = False,
            marks_aposteriori_type: Optional[str] = None,
            seed: Optional[int] = None
    ):
        logging.info(f"Particle process class initialized.")
        self.seed = seed
        self.model_name = model_name
        self.angles_matrix = None
        self.germ_intensity = germ_intensity
        self.particles = particles
        self.number_of_particles = len(self.particles)
        self.grain_type = grain_type
        self.space_dimension = space_dimension
        self.marked = marked
        self.marked_aposteriori = marked_aposteriori
        if self.marked_aposteriori:
            if not marks_aposteriori_type:
                raise ValueError(f"Parameter marks_aposteriori_type not filled for marked_aposteriori==True!")
            self.marks_aposteriori_type = marks_aposteriori_type
        if self.marked and not self.marked_aposteriori:
            self._compute_the_marks_matrices()
        self.particle_measure = self._compute_the_particles_measure()
        # compute the grains distance
        self.germs_distance_matrix = self._compute_the_germs_distance_matrix()
        # compute particles_ball_null_model distance (inf {||x_i - x_j||: x_i \in \Xi_i, x_j \in \Xi_j})
        self.particles_distance_matrix = self._compute_the_particles_distance_matrix()
        # particles_ball_null_model distance == 0 -> intersected, otherwise not intersected
        self.particles_intersection_matrix = self._compute_the_particles_intersection_matrix()
        # compute the pairwise shared corresponding Lebesgue measure
        # ("ball": shared areas, "segment": same as intersection matrix ...)
        self.pairwise_shared_measure_matrix = self._compute_the_pairwise_shared_measure_matrix()
        # if needed following attributes are computed via executing ParticleProcess.compute_f_mark_statistics
        self.f_mark_normalization_constants: Dict[str, float] = {}
        self.f_mark_statistics: Dict[Tuple[str, str], Dict[float, float]] = {}
        self.f_mark_statistics_permutations: Dict[Tuple[str, str], Dict[float, np.array]] = {
            val: np.array([]) for val in const.F_MARK_COMBINATIONS[self.grain_type]
        }
        self.f_mark_statistics_quantiles: Dict[Tuple[str, str],  Dict[float, float]] = {}
        if self.marked_aposteriori:
            self._compute_the_marks()

    def _compute_the_marks_matrices(self):
        self.marks = np.array([p.mark.mark_value for p in self.particles])
        self.marks_product = self.marks[..., None] * self.marks[None, ...]
        self.marks_square = (self.marks[..., None] - self.marks[None, ...]) ** 2 / 2
        self.max_mark = np.amax(self.marks)
        self.min_mark = np.amin(self.marks)

    def _compute_the_marks(self):
        pass

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
        germs_distance_matrix = pairwise_distances(
                [self.particles[k].germ for k in range(self.number_of_particles)]
            )
        self.germs_distance_from_origin = np.array(
            [self._norm(self.particles[k].germ) for k in range(self.number_of_particles)])
        logging.info(f"{datetime.now()} :Germs distance computation end.")
        return germs_distance_matrix

    def compute_the_statistics(self):
        pass

    def compute_the_f_mark_characteristics(self, set_of_f_mark_combinations=None):
        for f_type in const.F_MARK_TYPES:
            self._compute_the_f_mark_normalization_constants(f_type=f_type)
        if not set_of_f_mark_combinations:
            set_of_f_mark_combinations = const.F_MARK_COMBINATIONS[self.grain_type]
        for (f_type, weight_type) in set_of_f_mark_combinations:
            weight_matrix = {
                "intersection": self.particles_intersection_matrix,
                "shared_area": self.pairwise_shared_measure_matrix,
                "distance": self.particles_distance_matrix,
                "angle": self.angles_matrix
            }.get(weight_type)
            self.f_mark_statistics[f_type, weight_type] = {}
            points_to_eval = np.round(np.arange(0.01, 1, 0.01, dtype=float), 2)
            for t in points_to_eval:
                weight_matrix_r = weight_matrix.copy()
                for k in range(len(self.particles)):
                    # TODO this works only for balls
                    x_norm = self.germs_distance_from_origin[k]
                    r = self.particles[k].grain.radius
                    if x_norm > t + r:
                        weight_matrix_r[k, :] = 0
                self.f_mark_statistics[f_type, weight_type][np.round(t, 2)] = self._compute_the_f_mark_correlation(
                    f_type=f_type, weight_matrix=weight_matrix_r
                )

    def perform_the_permutation_test_for_f_mark_characteristics(self):
        for _ in range(const.PERMUTATION_TEST_REPEAT_COUNT):
            permutation = np.random.permutation(self.number_of_particles)
            marks_permuted = self.marks[permutation]
            for (f_type, weight_type) in const.F_MARK_COMBINATIONS[self.grain_type]:
                weight_matrix = {
                    "intersection": self.particles_intersection_matrix,
                    "shared_area": self.pairwise_shared_measure_matrix,
                    "distance": self.particles_distance_matrix,
                    "angle": self.angles_matrix
                }.get(weight_type)
                val = self._compute_the_f_mark_correlation(
                        f_type=f_type, weight_matrix=weight_matrix, marks_vector=marks_permuted
                    )
                self.f_mark_statistics_permutations[f_type, weight_type] = np.append(
                    self.f_mark_statistics_permutations[f_type, weight_type], val
                )
        for (f_type, weight_type) in const.F_MARK_COMBINATIONS[self.grain_type]:
            self.f_mark_statistics_quantiles[f_type, weight_type] = np.where(
                self.f_mark_statistics[f_type, weight_type] <= self.f_mark_statistics_permutations[f_type, weight_type],
                0, 1
            ).mean()

    def _compute_the_f_mark_correlation(
            self, f_type: str, weight_matrix: np.array, marks_vector: Optional[np.array] = None
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

    def plot_itself(self, show_germs: bool = False, win_min=-1, win_max=1, edge_effects = None):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(left=win_min, right=win_max)
        ax.set_ylim(bottom=win_min, top=win_max)
        self._plot_particles(ax=ax, fig=fig)
        # levels = [str(lv) for lv in range(self.particles[0].mark.number_of_levels)]
        # cols = [const.PARTICLE_COLORS_CHOICE[level] for level in range(self.particles[0].mark.number_of_levels)]
        if show_germs:
            for particle in self.particles:
                p_array = np.array([particle.germ])
                # fig, ax = plt.subplots()
                ax.scatter(p_array[:, 0], p_array[:, 1], color='black', marker=".")
        plt.xlim(win_min, win_max)
        plt.ylim(win_min, win_max)
        if edge_effects:
            length = win_max - win_min - 2 * edge_effects
            square = patches.Rectangle(
                (win_min + edge_effects, win_min + edge_effects), length, length, linewidth=1, edgecolor='r', facecolor='none')
            ax.add_patch(square)
        if const.SAVE_PLOTS:
            plt.savefig(f"generated_pics/{str(datetime.now()).replace(':','-')}_plot.png", dpi=1000)
        plt.show()
        plt.close()

    def _compute_the_particles_distance_matrix(self):
        raise NotImplementedError(
            f"Method called from the instance of general ParticleProcess class. "
            f"Please use the instance of some of the ParticleProcess subclasses (SegmentProcess, BallProcess)."
        )

    def _add_legend(
            self,
            ax, plt, fig,
            mark_type: str,
            face_color: Optional[str] = None,
            max_alpha: Optional[float] = None,
            min_alpha: Optional[float] = None
    ):
        if mark_type == "discrete":
            handles, labels = plt.gca().get_legend_handles_labels()
            labels = [int(lb) for lb in labels]
            by_label = dict(zip(labels, handles))
            by_label_sorted = {k: by_label[k] for k in sorted(by_label)}
            plt.figlegend(by_label_sorted.values(), by_label_sorted.keys())
        else:
            alphas = np.linspace(min_alpha, max_alpha, const.GRADIENT_SMOOTHNESS)
            r, g, b = tuple(int(face_color.lstrip('#')[i:i + 2], 16) / 256 for i in (0, 2, 4))
            cols = np.array([[r, g, b, alphas[k]] for k in range(const.GRADIENT_SMOOTHNESS)])
            psm = ax.pcolormesh(
                [[0, 0], [1, 1]], cmap=ListedColormap(cols), rasterized=True, vmin=self.min_mark, vmax=self.max_mark
            )
            fig.colorbar(psm, ax=ax)

    @staticmethod
    def _choose_edge_color():
        return "#000000"
        # return np.random.choice(const.PARTICLE_COLORS_CHOICE)

    def _choose_face_color(self, particle):
        if particle.mark is not None:
            if particle.mark.mark_type == "continuous":
                col = const.PARTICLE_COLORS_CHOICE[0]
                alpha = 0.2 + 0.8 * ((particle.mark.mark_value - self.min_mark) / (self.max_mark - self.min_mark))
            else:
                col = const.PARTICLE_COLORS_CHOICE[particle.mark.mark_value]
                alpha = const.ALPHA
        else:
            col = np.random.choice(const.PARTICLE_COLORS_CHOICE)
            alpha = const.ALPHA
        return col, alpha

    @staticmethod
    def _choose_germ_color():
        return "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]), np.random.random(1)[0]

    def _compute_the_particles_measure(self):
        pass
