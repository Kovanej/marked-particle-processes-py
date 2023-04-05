from datetime import datetime
import logging
import matplotlib.pyplot as plt
from typing import List, Optional, Union
import numpy as np
import random
from sklearn.metrics import pairwise_distances
from skspatial.objects import Point, Vector

from Geometry.particle import Particle
from Processes.markings import Mark
from Processes.particle_process import ParticleProcess
import utils.const as const


# TODO extend ball processes to in R^d, current state is for d = 2
class BallProcess(ParticleProcess):

    def __init__(
            self, germ_intensity: float, particles: List[Particle], max_radius: float, min_radius: float,
            marked: bool = False, model_name: Optional[str] = None, seed: Optional[int] = None,
            marked_aposteriori: Optional[bool] = False, marks_aposteriori_type: Optional[str] = None,
    ):
        self.radii_array = np.array([p.grain.radius for p in particles])
        self.max_radius = max_radius
        self.min_radius = min_radius
        super().__init__(
            germ_intensity=germ_intensity, grain_type="ball", particles=particles, space_dimension=2, marked=marked,
            model_name=model_name, seed=seed, marked_aposteriori=marked_aposteriori,
            marks_aposteriori_type=marks_aposteriori_type
        )

    def _compute_the_pairwise_shared_measure_matrix(self):
        # TODO subtract the area of circles outside the observational window
        # TODO vectorize
        logging.info(f"{datetime.now()} :Circles shared measure matrix computation start.")
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
        logging.info(f"{datetime.now()} :Circles shared measure matrix computation end.")
        return shared_areas_matrix

    def _compute_the_particles_distance_matrix(self):
        logging.info(f"{datetime.now()} :Circles distance computation start.")
        radii_to_subtract = np.array([
            [self.particles[k].grain.radius + self.particles[j].grain.radius for j in range(self.number_of_particles)]
            for k in range(self.number_of_particles)
        ])
        distance_matrix = np.where(
            self.germs_distance_matrix < radii_to_subtract, 0, self.germs_distance_matrix - radii_to_subtract
        )
        logging.info(f"{datetime.now()} :Circles distance computation end.")
        return distance_matrix

    def _plot_particles(self, ax, fig):
        min_alpha = 1
        max_alpha = 0
        for particle in self.particles:
            face_color, alpha = self._choose_face_color(particle=particle)
            if np.isnan(alpha):
                alpha = const.ALPHA
            edge_color = self._choose_edge_color()
            particle.grain.plot_2d(
                ax, facecolor=face_color, linestyle="-", alpha=alpha, linewidth=1, edgecolor=edge_color,
                label=particle.mark.mark_value
            )
            # TODO not hardcoded
            limits = plt.axis([-1, 1, -1, 1])
            min_alpha = alpha if min_alpha > alpha else min_alpha
            max_alpha = alpha if max_alpha < alpha else max_alpha
        self._add_legend(
            ax=ax, plt=plt, fig=fig, mark_type=particle.mark.mark_type, face_color=face_color, min_alpha=min_alpha,
            max_alpha=max_alpha
        )

    def _compute_the_particles_measure(self):
        # TODO currently for R^2 - improve with general formula later
        areas = (self.radii_array ** 2) * np.pi
        return areas


class RadiusMarksBallProcess(BallProcess):

    def __init__(
            self, germ_intensity: float, particles: List[Particle], model_name: Optional[str],max_radius: float,
            min_radius: float, seed: Optional[int] = None

    ):
        super().__init__(
            germ_intensity=germ_intensity, particles=particles, marked=True, model_name=model_name, seed=seed,
            min_radius=min_radius, max_radius=max_radius
        )


class BivariateRadiusMarksBallProcess(RadiusMarksBallProcess):

    def __init__(
            self, germ_intensity: float, particles: List[Particle], alpha: float, max_radius: float, min_radius: float,
            seed: Optional[int] = None
    ):
        self.alpha = alpha
        tau = np.random.binomial(n=1, p=self.alpha, size=len(particles))
        radii = np.array([particle.grain.radius for particle in particles])
        p_alt = np.where(tau == 0, 1/2, np.random.binomial(n=1, p=(radii - min_radius)/(max_radius - min_radius)))
        mark_values = np.random.binomial(n=1, p=p_alt)
        for k in range(len(particles)):
            mark_value = mark_values[k]
            mark = Mark(mark_type="discrete", mark_value=mark_value, number_of_levels=2)
            particles[k].mark = mark
        super().__init__(
            germ_intensity=germ_intensity, particles=particles, model_name=f"ball_bivariate_radius_alpha={self.alpha}",
            seed=seed, max_radius=max_radius, min_radius=min_radius
        )


class ContinuousRadiusMarksBallProcess(RadiusMarksBallProcess):

    def __init__(
            self, germ_intensity: float, particles: List[Particle], alpha: float, max_radius: float, min_radius: float,
            seed: Optional[int] = None
    ):
        self.alpha = alpha
        radii = np.array([particle.grain.radius for particle in particles])
        mark_values = self.alpha * radii + (1-self.alpha) * (min_radius + (
                max_radius - min_radius) * np.random.random(size=radii.size))
        for k in range(len(particles)):
            mark = Mark(mark_type="continuous", mark_value=mark_values[k])
            particles[k].mark = mark
        super().__init__(
            germ_intensity=germ_intensity, particles=particles, model_name=f"ball_continuous_radius_alpha={self.alpha}",
            seed=seed, max_radius=max_radius, min_radius=min_radius
        )


class BivariateMaximalSharedAreaMarkBallProcess(BallProcess):

    def __init__(
            self, germ_intensity: float, particles: List[Particle], alpha: float, min_radius: float, max_radius: float,
            seed: Optional[int] = None
    ):
        self.alpha = alpha
        super().__init__(
            germ_intensity=germ_intensity, particles=particles, marked=True, seed=seed,
            marked_aposteriori=True, marks_aposteriori_type="maximal_shared_area",
            model_name=f"ball_max_shared_area_disc_alpha={self.alpha}", max_radius=max_radius, min_radius=min_radius
        )
        self._mark_itself()

    def _mark_itself(self):
        max_shared_per_particle = self.pairwise_shared_measure_matrix.max(axis=0)
        max_possible_shared_area = (self.max_radius ** 2) * np.pi
        tau = np.random.binomial(n=1, p=self.alpha)
        p_alt = np.where(tau == 0, 1 / 2, np.random.binomial(n=1, p=max_shared_per_particle / max_possible_shared_area))
        mark_values = np.random.binomial(n=1, p=p_alt)
        for k in range(len(self.particles)):
            mark = Mark(mark_type="discrete", mark_value=mark_values[k])
            self.particles[k].mark = mark
        self._compute_the_marks_matrices()


class ContinuousMaximalSharedAreaMarkBallProcess(BallProcess):

    def __init__(
            self, germ_intensity: float, particles: List[Particle], alpha: float, min_radius: float, max_radius: float,
            seed: Optional[int] = None
    ):
        self.alpha = alpha
        super().__init__(
            germ_intensity=germ_intensity, particles=particles, marked=True, seed=seed,
            marked_aposteriori=True, marks_aposteriori_type="maximal_shared_area",
            model_name=f"ball_max_shared_area_cont_alpha={self.alpha}", max_radius=max_radius, min_radius=min_radius
        )
        self._mark_itself()

    def _mark_itself(self):
        max_shared_per_particle = self.pairwise_shared_measure_matrix.max(axis=0)
        # TODO compute correctly
        max_shared_area = self.pairwise_shared_measure_matrix.max()
        min_shared_area = self.pairwise_shared_measure_matrix.min()
        mark_values = self.alpha * max_shared_per_particle + (1 - self.alpha) * (
            min_shared_area + (max_shared_area - min_shared_area) * np.random.random(size=max_shared_per_particle.size)
        )
        for k in range(len(self.particles)):
            mark = Mark(mark_type="continuous", mark_value=mark_values[k])
            self.particles[k].mark = mark
        self._compute_the_marks_matrices()


class ContinuousNNDistanceMarkBallProcess(BallProcess):

    def __init__(
            self, germ_intensity: float, particles: List[Particle], alpha: float, min_radius: float, max_radius: float,
            seed: Optional[int] = None
    ):
        self.alpha = alpha
        super().__init__(
            germ_intensity=germ_intensity, particles=particles, marked=True, seed=seed,
            marked_aposteriori=True, marks_aposteriori_type="nn_dist",
            model_name=f"ball_N_N_dist_alpha={self.alpha}", max_radius=max_radius, min_radius=min_radius
        )
        self._mark_itself()

    def _mark_itself(self):
        part_dist_inf_diagonal = self.particles_distance_matrix.copy()
        np.fill_diagonal(part_dist_inf_diagonal, np.inf)
        nn_dist_per_particle = part_dist_inf_diagonal.min(axis=0)
        # TODO compute correctly
        max_distance = self.particles_distance_matrix.max()
        min_distance = nn_dist_per_particle.min()
        mark_values = self.alpha * nn_dist_per_particle + (1 - self.alpha) * (
                min_distance + (max_distance - min_distance) * np.random.random(size=nn_dist_per_particle.size)
        )
        for k in range(len(self.particles)):
            mark = Mark(mark_type="continuous", mark_value=mark_values[k])
            self.particles[k].mark = mark
        self._compute_the_marks_matrices()


class CountingIntersectionNumberMarkBallProcess(BallProcess):

    def __init__(
            self, germ_intensity: float, particles: List[Particle], alpha: float, min_radius: float, max_radius: float,
            seed: Optional[int] = None
    ):
        self.alpha = alpha
        super().__init__(
            germ_intensity=germ_intensity, particles=particles, marked=True, seed=seed,
            marked_aposteriori=True, marks_aposteriori_type="nn_dist",
            model_name=f"ball_intersection_count_alpha={self.alpha}", max_radius=max_radius, min_radius=min_radius
        )
        self._mark_itself()

    def _mark_itself(self):
        tau = np.random.binomial(n=1, p=self.alpha, size=len(self.particles))
        # subtracting the one, since particles intersect themselves
        intersection_count = self.particles_intersection_matrix.sum(axis=0) - 1
        # TODO compute properly
        independent_poisson = np.random.poisson(intersection_count.mean(), size=intersection_count.size)
        dependent_poisson = np.random.poisson(intersection_count)
        mark_values = np.where(tau == 0, independent_poisson, dependent_poisson)
        for k in range(len(self.particles)):
            mark = Mark(mark_type="discrete", mark_value=mark_values[k])
            self.particles[k].mark = mark
        self._compute_the_marks_matrices()
