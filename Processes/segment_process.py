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

    def __init__(
            self, germ_intensity: float, particles: List[Particle], marked: bool = False, space_dimension: int = 2
    ):
        logging.info(f"{datetime.now()} Segment process init start.")
        super().__init__(
            germ_intensity=germ_intensity, grain_type="segment", particles=particles, space_dimension=space_dimension,
            marked=marked
        )
        logging.info(f"{datetime.now()} Segment process angles matrix computation start.")
        self.angles_matrix = self._compute_the_angles_matrix()
        logging.info(f"{datetime.now()} Segment process angles matrix computation end.")
        logging.info(f"{datetime.now()} Segment process init end.")

    def _compute_the_angles_matrix(self):
        a = np.array([particle.grain.vector for particle in self.particles])
        coss = self._dot_pairwise(a, a) / (self._norm(a)[:, None] * self._norm(a))
        coss_clipped = np.clip(coss, -1, 1)
        angles_matrix = np.arccos(
            coss_clipped
        )
        return angles_matrix

    def _compute_the_pairwise_shared_measure_matrix(self):
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
                np.clip(solution[..., 0], 0, 1), np.clip(solution[..., 1], 0, 1), start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_0_solution[..., 0], 0, 1).T, np.clip(alpha_0_solution[..., 1], 0, 1).T, start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_1_solution[..., 0], 0, 1).T, np.clip(alpha_1_solution[..., 1], 0, 1).T, start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(beta_0_solution[..., 0], 0, 1).T, np.clip(beta_0_solution[..., 1], 0, 1).T, start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(beta_1_solution[..., 0], 0, 1).T, np.clip(beta_1_solution[..., 1], 0, 1).T, start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_0_beta_0[..., 0], 0, 1).T, np.clip(alpha_0_beta_0[..., 1], 0, 1).T, start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_0_beta_1[..., 0], 0, 1).T, np.clip(alpha_0_beta_1[..., 1], 0, 1).T, start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_1_beta_0[..., 0], 0, 1).T, np.clip(alpha_1_beta_0[..., 1], 0, 1).T, start_points_pre, end_points_pre
            ),
            self._segments_vectorized_pairwise_distance(
                np.clip(alpha_1_beta_1[..., 0], 0, 1).T, np.clip(alpha_1_beta_1[..., 1], 0, 1).T, start_points_pre, end_points_pre
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
        left = alpha * start_points[:, None, :] + (1 - alpha) * end_points[:, None, :]
        right = beta * start_points[None, :, :] + (1 - beta) * end_points[None, :, :]
        res = left - right
        final = np.sqrt(np.power(res, 2).sum(axis=-1))
        return final.T

    @staticmethod
    def _return_convex_combination(start_point: Point, end_point: Point, alpha: float):
        return alpha * start_point + (1 - alpha) * end_point

    def _plot_particles(self, ax, fig):
        if self.space_dimension == 3:
            ax = fig.add_subplot(111, projection='3d')
        for particle in self.particles:
            # col, alpha = self._choose_face_color()
            if particle.mark is not None:
                if particle.mark.mark_value == 0:
                    col, alpha = "#0075C2", 1
                elif particle.mark.mark_value == 1:
                    col, alpha = "#F87203", 1
                else:
                    col, alpha = self._choose_face_color(particle=particle)
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

    def _compute_the_particles_measure(self):
        vectors = np.array([p.grain.end_point - p.grain.start_point for p in self.particles])
        lengths = self._norm(vectors)
        return lengths
