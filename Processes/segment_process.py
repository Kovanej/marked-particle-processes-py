import logging
from datetime import datetime
import matplotlib.pyplot as plt
from typing import List, Optional, Union
import numpy as np
from scipy.optimize import fsolve
from sklearn.metrics import pairwise_distances
from skspatial.objects import Point, Vector

from Geometry.particle import Particle
from Processes.particle_process import ParticleProcess
import utils.const as const


class SegmentProcess(ParticleProcess):

    def __init__(self, germ_intensity: float, particles: List[Particle], space_dimension: int = 2):
        logging.info(f"{datetime.now()} Segment process init start.")
        super().__init__(
            germ_intensity=germ_intensity, grain_type="segment", particles=particles, space_dimension=space_dimension
        )
        logging.info(f"{datetime.now()} Segment process angles matrix computation start.")
        self.angles_matrix = self._compute_the_angles_matrix()
        logging.info(f"{datetime.now()} Segment process angles matrix computation end.")
        logging.info(f"{datetime.now()} Segment process init end.")

    @staticmethod
    def _dot_pairwise(a, b):
        return (a[:, None, :] * b[None, ...]).sum(axis=-1)

    @staticmethod
    def _norm(a):
        return np.sqrt((a * a).sum(axis=-1))

    def _compute_the_angles_matrix(self):
        a = np.array([particle.grain.vector for particle in self.particles])
        coss = self._dot_pairwise(a, a) / (self._norm(a)[:, None] * self._norm(a))
        coss_clipped = np.clip(coss, -1, 1)
        angles_matrix = np.arccos(
            coss_clipped
        )
        return angles_matrix

    # TODO delete later, right now saved for performance comparison with vectorized version
    def _compute_the_angles_matrix_non_vectorized(self):
        angle_matrix = np.zeros(shape=self.germs_distance_matrix.shape)
        for i in range(self.number_of_particles):
            for j in range(self.number_of_particles):
                angle = self.particles[i].grain.vector.angle_between(self.particles[j].grain.vector)
                angle_matrix[i, j] = angle
                angle_matrix[j, i] = angle
        return angle_matrix

    def _compute_the_shared_corresponding_measure_matrix(self):
        # TODO this generally doesn't hold for segments on the same line - fix
        logging.info(f"{datetime.now()} :Segments shared measure matrix computation start.")
        shared_measure_matrix = self.particles_intersection_matrix
        logging.info(f"{datetime.now()} :Segments shared measure matrix computation end.")
        return shared_measure_matrix

    def _compute_the_particles_distance_matrix(self):
        logging.info(f"{datetime.now()} :Segments distance computation start.")
        timer_start = datetime.now()
        logging.info(f"{datetime.now()} DISTANCE - Setup of needed values start.")
        start_points_pre = np.array([particle.grain.start_point for particle in self.particles])
        end_points_pre = np.array([particle.grain.end_point for particle in self.particles])
        start_points_1 = start_points_pre[:, None, :]
        end_points_1 = end_points_pre[:, None, :]
        start_points_2 = start_points_pre[None, ...]
        end_points_2 = end_points_pre[None, ...]
        A_11 = ((start_points_1 - end_points_1) ** 2).sum(axis=-1)
        A_12 = ((start_points_1 - end_points_1) * (end_points_2 - start_points_2)).sum(axis=2)
        A_21 = ((start_points_1 - end_points_1) * (end_points_2 - start_points_2)).sum(axis=2)
        A_22 = ((start_points_2 - end_points_2) ** 2).sum(axis=-1)
        A_11 = np.repeat(A_11, repeats=self.number_of_particles, axis=1)
        A_22 = np.repeat(A_22, repeats=self.number_of_particles, axis=0)
        A = np.array([
            np.array([A_22[i, j], A_12[i, j], A_21[i, j], A_11[i, j]]).reshape((2, 2))
            for j in range(self.number_of_particles) for i in range(self.number_of_particles)
        ]).reshape((self.number_of_particles, self.number_of_particles, 2, 2))
        y_1 = ((end_points_2 - end_points_1) * (start_points_1 - end_points_1)).sum(axis=-1)
        y_2 = ((end_points_2 - end_points_1) * (end_points_2 - start_points_2)).sum(axis=-1)
        y = np.array([y_1, y_2]).transpose((1, 2, 0))
        idx = np.array(range(self.number_of_particles))
        # we want to skip computations for two identical segments - we can take arbitrary alpha, beta from [0, 1]
        A[idx, idx, ...] = np.eye(2)
        y[idx, idx, :] = 0
        logging.info(f"{datetime.now()} DISTANCE - Setup of needed values end.")
        logging.info(f"{datetime.now()} DISTANCE - Linalg.solve start.")
        solution = np.linalg.solve(A, y)
        logging.info(f"{datetime.now()} DISTANCE - Linalg.solve end.")
        logging.info(f"{datetime.now()} DISTANCE - Simpler solution computation start.")
        alpha_0_solution = np.array([
            np.zeros(shape=(self.number_of_particles, self.number_of_particles)), (y_2 / A_22)
        ]).T.reshape(self.number_of_particles, self.number_of_particles, 2)
        alpha_1_solution = np.array([
            np.ones(shape=(self.number_of_particles, self.number_of_particles)), ((y_2 - A_21) / A_22)
        ]).T.reshape(self.number_of_particles, self.number_of_particles, 2)
        beta_0_solution = np.array([
            (y_1 / A_11), np.zeros(shape=(self.number_of_particles, self.number_of_particles))
        ]).T.reshape(self.number_of_particles, self.number_of_particles, 2)
        beta_1_solution = np.array([
            ((y_1 - A_12) / A_11), np.ones(shape=(self.number_of_particles, self.number_of_particles))
        ]).T.reshape(self.number_of_particles, self.number_of_particles, 2)
        alpha_0_beta_0 = np.zeros(shape=(self.number_of_particles, self.number_of_particles, 2))
        alpha_1_beta_1 = np.ones(shape=(self.number_of_particles, self.number_of_particles, 2))
        alpha_0_beta_1 = np.array([
            np.zeros(shape=(self.number_of_particles, self.number_of_particles)),
            np.ones(shape=(self.number_of_particles, self.number_of_particles))
        ]).T.reshape(self.number_of_particles, self.number_of_particles, 2)
        alpha_1_beta_0 = np.array([
            np.ones(shape=(self.number_of_particles, self.number_of_particles)),
            np.zeros(shape=(self.number_of_particles, self.number_of_particles))
        ]).T.reshape(self.number_of_particles, self.number_of_particles, 2)
        logging.info(f"{datetime.now()} DISTANCE - Simpler solution computation end.")
        logging.info(f"{datetime.now()} DISTANCE - Possible distances computation start.")
        # TODO this is a slow performance part for some reason - fix this
        possible_distances = np.array([
            self._segments_vectorized_pairwise_distance(
                np.clip(solution[..., 0], 0, 1).T, np.clip(solution[..., 1], 0, 1).T, start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_0_solution[..., 0], 0, 1), np.clip(alpha_0_solution[..., 1], 0, 1), start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_1_solution[..., 0], 0, 1), np.clip(alpha_1_solution[..., 1], 0, 1), start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(beta_0_solution[..., 0], 0, 1), np.clip(beta_0_solution[..., 1], 0, 1), start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(beta_1_solution[..., 0], 0, 1), np.clip(beta_1_solution[..., 1], 0, 1), start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_0_beta_0[..., 0], 0, 1), np.clip(alpha_0_beta_0[..., 1], 0, 1), start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_0_beta_1[..., 0], 0, 1), np.clip(alpha_0_beta_1[..., 1], 0, 1), start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_1_beta_0[..., 0], 0, 1), np.clip(alpha_1_beta_0[..., 1], 0, 1), start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_1_beta_1[..., 0], 0, 1), np.clip(alpha_1_beta_1[..., 1], 0, 1), start_points_pre, end_points_pre
            ),
        ]).T
        logging.info(f"{datetime.now()} DISTANCE - Possible distances computation end.")
        distances = np.round(possible_distances.min(axis=-1), decimals=8)
        timer_end = datetime.now()
        logging.info(f"SEGMENT DISTANCE COMPUTATION LENGTH: {timer_end - timer_start}")
        logging.info(f"{datetime.now()} :Particle distance computation end.")
        return distances

    @staticmethod
    def _segments_vectorized_pairwise_distance(alpha, beta, start_points, end_points):
        alpha = alpha[..., None]
        beta = beta[..., None]
        left = alpha * start_points[:, None, None, :] + (1 - alpha) * end_points[:, None, None, :]
        right = beta * start_points[None, :, None, :] + (1 - beta) * end_points[None, :, None, :]
        res = left - right
        final = np.sqrt(np.power(res, 2).sum(axis=-1)).min(axis=-1)
        return final

    def _return_distance_vector(
            self, start_point_1: Point, end_point_1: Point, start_point_2: Point, end_point_2: Point,
            alpha: float, beta: float
    ):
        return self._return_convex_combination(
            start_point=start_point_1, end_point=end_point_1, alpha=alpha
        ) - self._return_convex_combination(
            start_point=start_point_2, end_point=end_point_2, alpha=beta
        )

    @staticmethod
    def _return_convex_combination(start_point: Point, end_point: Point, alpha: float):
        return alpha * start_point + (1 - alpha) * end_point

    def _compute_the_particles_distance_matrix_non_vectorized(self):
        timer_start = datetime.now()
        distance_matrix = np.zeros(shape=self.germs_distance_matrix.shape)
        for i in range(self.number_of_particles):
            for j in range(i, self.number_of_particles):
                if j != i:
                    initial_alpha_guess, initial_beta_guess = 1 / 2, 1 / 2
                    possible_alphas_betas = []
                    alpha_beta_opti = fsolve(self._optim_func_segments_derivative,
                                             x0=[initial_alpha_guess, initial_beta_guess],
                                             args=(self.particles[i].grain.start_point,
                                                   self.particles[i].grain.end_point,
                                                   self.particles[j].grain.start_point,
                                                   self.particles[j].grain.end_point,
                                                   )
                                             )
                    if all(0 <= _ <= 1 for _ in alpha_beta_opti):
                        possible_alphas_betas.append((alpha_beta_opti[0], alpha_beta_opti[1]))
                    alpha_opti_beta_0 = fsolve(self._optim_func_derivative_fixed_beta_0,
                                               x0=initial_alpha_guess,
                                               args=(self.particles[i].grain.start_point,
                                                     self.particles[i].grain.end_point,
                                                     self.particles[j].grain.start_point,
                                                     self.particles[j].grain.end_point,
                                                     )
                                               )[0], 0
                    if all(0 <= _ <= 1 for _ in alpha_opti_beta_0):
                        possible_alphas_betas.append(alpha_opti_beta_0)
                    alpha_opti_beta_1 = fsolve(self._optim_func_derivative_fixed_beta_1,
                                               x0=initial_alpha_guess,
                                               args=(self.particles[i].grain.start_point,
                                                     self.particles[i].grain.end_point,
                                                     self.particles[j].grain.start_point,
                                                     self.particles[j].grain.end_point,
                                                     ))[0], 1
                    if all(0 <= _ <= 1 for _ in alpha_opti_beta_1):
                        possible_alphas_betas.append(alpha_opti_beta_1)
                    alpha_0_beta_opti = 0, fsolve(self._optim_func_derivative_fixed_alpha_0,
                                                  x0=initial_beta_guess,
                                                  args=(self.particles[i].grain.start_point,
                                                        self.particles[i].grain.end_point,
                                                        self.particles[j].grain.start_point,
                                                        self.particles[j].grain.end_point,
                                                        ))[0]
                    if all(0 <= _ <= 1 for _ in alpha_0_beta_opti):
                        possible_alphas_betas.append(alpha_0_beta_opti)
                    alpha_1_beta_opti = 1, fsolve(self._optim_func_derivative_fixed_alpha_1,
                                                  x0=initial_beta_guess,
                                                  args=(self.particles[i].grain.start_point,
                                                        self.particles[i].grain.end_point,
                                                        self.particles[j].grain.start_point,
                                                        self.particles[j].grain.end_point,
                                                        ))[0]
                    if all(0 <= _ <= 1 for _ in alpha_1_beta_opti):
                        possible_alphas_betas.append(alpha_1_beta_opti)
                    alpha_0_beta_0 = 0, 0
                    possible_alphas_betas.append(alpha_0_beta_0)
                    alpha_0_beta_1 = 0, 1
                    possible_alphas_betas.append(alpha_0_beta_1)
                    alpha_1_beta_0 = 1, 0
                    possible_alphas_betas.append(alpha_1_beta_0)
                    alpha_1_beta_1 = 1, 1
                    possible_alphas_betas.append(alpha_1_beta_1)
                    possible_minimas = [
                        (
                                a * self.particles[i].grain.start_point + (1 - a) * self.particles[i].grain.end_point
                        ).distance_point(
                            b * self.particles[j].grain.start_point + (1 - b) * self.particles[j].grain.end_point
                        )
                        for a, b in possible_alphas_betas
                    ]
                    distance_matrix[i, j] = min(possible_minimas)
                    distance_matrix[j, i] = min(possible_minimas)
        distance_matrix_round = np.around(distance_matrix, decimals=8)
        timer_end = datetime.now()
        logging.info(f"SEGMENT DISTANCE PREVIOUS LENGTH: {timer_end - timer_start}")
        return distance_matrix_round

    def _compute_the_pairwise_angle_matrix(self):
        return np.zeros(shape=self.particles_distance_matrix.shape)

    @staticmethod
    def _optim_func_segments_derivative(
            params,
            start_point_1: Point, end_point_1: Point,
            start_point_2: Point, end_point_2: Point,
            alpha_fixed=None, beta_fixed=None
    ):
        if alpha_fixed is not None:
            alpha = alpha_fixed
            beta = params
        elif beta_fixed is not None:
            alpha = params
            beta = beta_fixed
        else:
            alpha, beta = params
        alpha_func = sum([
            (
                    alpha * (start_point_1[i] - end_point_1[i])
                    + beta * (end_point_2[i] - start_point_2[i])
                    + end_point_1[i] - end_point_2[i]
            ) * (
                    end_point_1[i] - start_point_1[i]
            )
            for i in range(len(start_point_1))
        ])
        beta_func = sum([
            (
                    alpha * (start_point_1[i] - end_point_1[i]) +
                    beta * (end_point_2[i] - start_point_2[i]) +
                    end_point_1[i] - end_point_2[i]
            ) * (
                    end_point_2[i] - start_point_2[i]
            )
            for i in range(len(start_point_1))
        ])
        if alpha_fixed is not None:
            return beta_func
        elif beta_fixed is not None:
            return alpha_func
        else:
            return (alpha_func, beta_func)

    def _optim_func_derivative_fixed_beta_0(
            self, params,
            start_point_1: Point, end_point_1: Point,
            start_point_2: Point, end_point_2: Point
    ):
        return self._optim_func_segments_derivative(
            params, beta_fixed=0,
            start_point_1=start_point_1, end_point_1=end_point_1, start_point_2=start_point_2, end_point_2=end_point_2
        )

    def _optim_func_derivative_fixed_alpha_0(
            self, params,
            start_point_1: Point, end_point_1: Point,
            start_point_2: Point, end_point_2: Point
    ):
        return self._optim_func_segments_derivative(
            params, alpha_fixed=0,
            start_point_1=start_point_1, end_point_1=end_point_1, start_point_2=start_point_2, end_point_2=end_point_2
        )

    def _optim_func_derivative_fixed_beta_1(
            self, params,
            start_point_1: Point, end_point_1: Point,
            start_point_2: Point, end_point_2: Point
    ):
        return self._optim_func_segments_derivative(
            params, beta_fixed=1,
            start_point_1=start_point_1, end_point_1=end_point_1, start_point_2=start_point_2, end_point_2=end_point_2
        )

    def _optim_func_derivative_fixed_alpha_1(
            self, params,
            start_point_1: Point, end_point_1: Point,
            start_point_2: Point, end_point_2: Point
    ):
        return self._optim_func_segments_derivative(
            params, alpha_fixed=1,
            start_point_1=start_point_1, end_point_1=end_point_1, start_point_2=start_point_2, end_point_2=end_point_2
        )

    def _plot_particles(self, ax, fig):
        if self.space_dimension == 3:
            ax = fig.add_subplot(111, projection='3d')
        for particle in self.particles:
            # col, alpha = self._choose_face_color()
            if particle.mark is not None:
                if particle.mark.mark_value == 0:
                    col, alpha = "#003271", 1
                elif particle.mark.mark_value == 1:
                    col, alpha = "#FEC500", 1
            else:
                col, alpha = np.random.choice(const.PARTICLE_COLORS_CHOICE), 1
            # alpha = Vector(particle.grain.start_point).norm() / np.sqrt(2)
            if self.space_dimension == 2:
                particle.grain.vector.plot_2d(
                    ax_2d=ax, point=particle.grain.start_point, head_width=0,
                    edgecolor=col, alpha=alpha
                )
            elif self.space_dimension == 3:
                particle.grain.vector.plot_3d(
                    ax_3d=ax, point=particle.grain.start_point,
                    #  edgecolor=col, alpha=alpha
                )
