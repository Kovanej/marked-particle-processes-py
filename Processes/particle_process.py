
from datetime import datetime
import matplotlib.pyplot as plt
from typing import List, Optional, Union
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
    ):
        self.germ_intensity = germ_intensity
        self.particles = particles
        self.number_of_particles = len(self.particles)
        self.grain_type = grain_type
        # compute the grains distance
        self.germs_distance_matrix = self._compute_the_germs_distance_matrix()
        # compute particles_null_model distance (inf {||x_i - x_j||: x_i \in \Xi_i, x_j \in \Xi_j})
        self.particles_distance_matrix = self._compute_the_particles_distance_matrix()
        # particles_null_model distance == 0 -> intersected, otherwise not intersected
        self.particles_intersection_matrix = np.where(self.particles_distance_matrix == 0, 1, 0)
        # compute the pairwise shared corresponding Lebesgue measure
        # ("ball": shared areas, "segment": same as intersection matrix ...)
        self.shared_corresponding_measure_matrix = self._compute_the_shared_corresponding_measure_matrix()
        # if needed following attributes are computed via executing ParticleProcess.compute_f_mark_statistics
        self.f_mark_normalization_constant = None
        self.f_mark_intersection_correlation = None

    def _compute_the_shared_corresponding_measure_matrix(self):
        # overridden for each subclasses - if not specified, it cannot be computed
        raise NotImplementedError(
            f"Method called from the instance of general ParticleProcess class. "
            f"Please use the instance of some of the ParticleProcess subclasses (SegmentProcess, BallProcess)."
        )

    def _compute_the_germs_distance_matrix(self):
        grains_distance_matrix = pairwise_distances(
                [self.particles[k].germ for k in range(self.number_of_particles)]
            )
        return grains_distance_matrix

    def _plot_circle_particles(self, ax):
        for particle in self.particles:
            facecolor, alpha = self._choose_face_color()
            edgecolor = self._choose_edge_color()
            particle.grain.plot_2d(
                ax, facecolor=facecolor, linestyle="-", alpha=alpha, linewidth=1, edgecolor=edgecolor,
                # alpha=0.5,
            )

    def _plot_segment_particles(self, ax):
        for particle in self.particles:
            # col, alpha = self._choose_face_color()
            if particle.mark is not None:
                if particle.mark.mark_value == 0:
                    col, alpha = "#ff0000", 1
                elif particle.mark.mark_value == 1:
                    col, alpha = "#0000ff", 1
            else:
                col, alpha = "#000000", 1
            # alpha = Vector(particle.grain.start_point).norm() / np.sqrt(2)
            particle.grain.vector.plot_2d(
                ax_2d=ax, point=particle.grain.start_point, head_width=0,
                edgecolor=col, alpha=alpha
            )

    def compute_the_f_mark_characteristics(
            self,
            f_type: str = "product"
    ):
        self.f_mark_normalization_constant = self._compute_the_f_mark_normalization_constant(f_type=f_type)
        self.f_mark_intersection_correlation = self._compute_the_f_mark_intersection_correlation(f_type=f_type)

    def _compute_the_f_mark_intersection_correlation(self, f_type: str):
        f_intersection_nn = 0
        # intersection_counter = 0
        # pairs_count = (self.number_of_particles * (self.number_of_particles - 1)) / 2
        for i in range(self.number_of_particles):
            for j in range(i + 1, self.number_of_particles):
                if self.particles_intersection_matrix[i, j] == 1:
                    # TODO temporarily hardcoded f(M_1, M_2) = M_1 * M_2
                    f_intersection_nn += self.particles[i].mark.mark_value * self.particles[j].mark.mark_value
                    # intersection_counter += 1
        # if intersection_counter == 0:
        #     intersection_counter = 1
        # multiply by 2 since we count only for (i, j) j > i
        # TODO fix this^^, since not all f(M_i, M_j) are symmetrical
        self.f_intersection_nn = f_intersection_nn
        return (2 * f_intersection_nn) / (self.f_mark_normalization_constant * (self.germ_intensity ** 2))

    def _compute_the_f_mark_normalization_constant(self, f_type: str):
        pairs_count = (self.number_of_particles * (self.number_of_particles - 1)) / 2
        c_f_nn = 0
        for i in range(self.number_of_particles):
            for j in range(i + 1, self.number_of_particles):
                # TODO temporarily hardcoded f(M_1, M_2) = M_1 * M_2
                c_f_nn += self.particles[i].mark.mark_value * self.particles[j].mark.mark_value
        c_f = c_f_nn / pairs_count
        return c_f

    def plot_itself(self, show_germs: bool = False):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_aspect('equal', adjustable='box')
        if self.grain_type == "ball":
            self._plot_circle_particles(ax=ax)
        if self.grain_type == "segment":
            self._plot_segment_particles(ax=ax)
        if show_germs:
            for particle in self.particles:
                color, alpha = self._choose_germ_color()
                particle.germ.plot_2d(ax, c=color, alpha=alpha)
        if const.SAVE_PLOTS:
            plt.savefig(f"generated_pics/{str(datetime.now()).replace(':','-')}_plot.png", dpi=600)
        plt.show()

    def _compute_the_particles_distance_matrix(self):
        raise NotImplementedError(
            f"Method called from the instance of general ParticleProcess class. "
            f"Please use the instance of some of the ParticleProcess subclasses (SegmentProcess, BallProcess)."
        )

    def _choose_edge_color(self):
        return "#" + ''.join([random.choice('ABCDEF0123456789') for i in range(6)])

    def _choose_face_color(self, particle=None):
        # TODO later adjust for marks
        alpha = np.random.random(1)[0]
        if particle is None:
            col = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
        elif self.grain_type == "segment":
            col = "#"
            r = int(particle.grain.angle / np.pi * 255)
            g = 255 - int(particle.grain.angle / np.pi * 255)
            b = 255 - int(particle.grain.angle / np.pi * 255)
            for col_part in [r, g, b]:
                hex_col_part = hex(col_part)[2:]
                if len(hex_col_part) == 1:
                    hex_col_part = "0" + hex_col_part
                col = col + hex_col_part
        else:
            # todo later
            col = "#" + ''.join([random.choice('ABCDEF0123456789') for i in range(6)])
        return col, alpha

    def _choose_germ_color(self):
        return "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]), np.random.random(1)[0]


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


class BallProcess(ParticleProcess):

    def __init__(self, germ_intensity: float, particles: List[Particle]):
        super().__init__(germ_intensity=germ_intensity, grain_type="ball", particles=particles)

    def _compute_the_shared_corresponding_measure_matrix(self):
        # TODO subtract the area of circles outside the observational window
        shared_areas_matrix = np.zeros(shape=self.particles_intersection_matrix.shape)
        for i in range(self.number_of_particles):
            for j in range(i, self.number_of_particles):
                if i == j:
                    shared_areas_matrix[i, j] = self.particles[i].lebesgue_measure
                elif self.particles_intersection_matrix[i, j] == 0:
                    shared_areas_matrix[i, j] = 0
                    shared_areas_matrix[j, i] = 0
                else:
                    dist = self.germs_distance_matrix[i, j]
                    rad_i_sq = self.particles[i].grain.radius ** 2
                    rad_j_sq = self.particles[j].grain.radius ** 2
                    x = (rad_i_sq - rad_j_sq + dist ** 2) / (2 * dist)
                    z = x ** 2

                    if rad_i_sq - z < 0 and dist > abs(self.particles[i].grain.radius - self.particles[j].grain.radius):
                        a = 1

                    if dist <= abs(self.particles[i].grain.radius - self.particles[j].grain.radius):
                        shared_areas_matrix[i, j] = np.pi * min(rad_i_sq, rad_j_sq)
                        shared_areas_matrix[j, i] = np.pi * min(rad_i_sq, rad_j_sq)
                    else:
                        y = np.sqrt(rad_i_sq - z)
                        _share = rad_i_sq * np.arcsin((y / self.particles[i].grain.radius)) + \
                                 rad_j_sq * np.arcsin((y / self.particles[j].grain.radius)) - \
                                 y * (x + np.sqrt(z + rad_j_sq - rad_i_sq))
                        shared_areas_matrix[i, j] = _share
                        shared_areas_matrix[j, i] = _share
        return shared_areas_matrix

    def _compute_the_particles_distance_matrix(self):
        radii_to_subtract = np.array([
            [
                self.particles[k].grain.radius + self.particles[j].grain.radius
                for j in range(self.number_of_particles)
            ] for k in range(self.number_of_particles)
        ])
        distance_matrix = np.where(
            self.germs_distance_matrix < radii_to_subtract, 0, self.germs_distance_matrix - radii_to_subtract
        )
        return distance_matrix

