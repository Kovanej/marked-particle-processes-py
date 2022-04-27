
import matplotlib.pyplot as plt
from typing import List
from itertools import product
import numpy as np
import random
from sklearn.metrics import pairwise_distances

from Geometry.particle import Particle


class ParticleProcess(object):

    def __init__(
            self,
            particles: List[Particle],
            grain_type: str,
    ):
        self.particles = particles
        self.number_of_particles = len(self.particles)
        self.grain_type = grain_type
        # compute the grains distance
        self.germs_distance_matrix = self._compute_the_germs_distance_matrix()
        # compute particles distance (inf {||x_i - x_j||: x_i \in \Xi_i, x_j \in \Xi_j})
        self.particles_distance_matrix = self._compute_the_particles_distance_matrix()
        # particles distance == 0 -> intersected, otherwise not intersected
        self.particles_intersection_matrix = np.where(self.particles_distance_matrix == 0, 1, 0)
        # compute the pairwise shared corresponding Lebesgue measure
        # ("circle": shared areas, "segment": same as intersection matrix ...)
        self.shared_corresponding_measure_matrix = self._compute_the_shared_corresponding_measure_matrix()

    def _compute_the_shared_corresponding_measure_matrix(self):
        if self.grain_type == "circle":
            shared_area_matrix = self._compute_the_pairwise_shared_areas_of_circles()
            return shared_area_matrix
        else:
            raise ValueError(f"Unknown value for Particle.grain_type: {self.grain_type}")

    def _compute_the_pairwise_shared_areas_of_circles(self):
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
                        a=1

                    if dist <= abs(self.particles[i].grain.radius - self.particles[j].grain.radius):
                        shared_areas_matrix[i, j] = np.pi * min(rad_i_sq, rad_j_sq)
                        shared_areas_matrix[j, i] = np.pi * min(rad_i_sq, rad_j_sq)
                    else:
                        y = np.sqrt(rad_i_sq - z)
                        # return a * asin(y / A.r) + b * asin(y / B.r) - y * (x + sqrt(z + b - a))
                        _share = rad_i_sq * np.arcsin((y / self.particles[i].grain.radius)) + \
                                 rad_j_sq * np.arcsin((y / self.particles[j].grain.radius)) - \
                                 y * (x + np.sqrt(z + rad_j_sq - rad_i_sq))
                        shared_areas_matrix[i, j] = _share
                        shared_areas_matrix[j, i] = _share
        return shared_areas_matrix

    def _compute_the_germs_distance_matrix(self):
        grains_distance_matrix = pairwise_distances(
                [self.particles[k].germ for k in range(self.number_of_particles)]
            )
        return grains_distance_matrix

    def plot_itself(self, show_germs: bool = False):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_aspect('equal', adjustable='box')
        for particle in self.particles:
            facecolor, alpha = self._choose_face_color()
            edgecolor = self._choose_edge_color()
            particle.grain.plot_2d(ax, facecolor=facecolor, linestyle="-", alpha=alpha, linewidth=1, edgecolor=edgecolor
                                   # alpha=0.5,
            )
        if show_germs:
            for particle in self.particles:
                color, alpha = self._choose_germ_color()
                particle.germ.plot_2d(ax, c=color, alpha=alpha)
        plt.show()

    def _compute_the_particles_distance_matrix(self):
        if self.grain_type == "circle":
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
        else:
            raise ValueError(f"Unknown value for Particle.grain_type: {self.grain_type}")

    def _choose_edge_color(self):
        return "#" + ''.join([random.choice('ABCDEF0123456789') for i in range(6)])

    def _choose_face_color(self):
        return "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]), np.random.random(1)[0]

    def _choose_germ_color(self):
        return "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]), np.random.random(1)[0]


