from datetime import datetime
import matplotlib.pyplot as plt
from typing import List, Optional, Union
import numpy as np
import random
from scipy.optimize import fsolve
from sklearn.metrics import pairwise_distances
from skspatial.objects import Point, Vector

from Geometry.particle import Particle
from Processes.particle_process import ParticleProcess
import utils.const as const


class SegmentProcess(ParticleProcess):

    def __init__(self, germ_intensity: float, particles: List[Particle]):
        super().__init__(germ_intensity=germ_intensity, grain_type="segment", particles=particles)

    def _compute_the_shared_corresponding_measure_matrix(self):
        return self.particles_intersection_matrix

    def _compute_the_particles_distance_matrix(self):
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
        # distance_matrix = self.germs_distance_matrix
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
