
import matplotlib.pyplot as plt
from typing import List
from itertools import product
import numpy as np
from sklearn.metrics import pairwise_distances

from Geometry.particle import Particle


class ParticleProcess(object):

    def __init__(
            self,
            particles: List[Particle],
            grain_type: str
    ):
        self.particles = particles
        self.number_of_particles = len(self.particles)
        self.grain_type = grain_type
        self.grains_distance_matrix = None
        self.particles_distance_matrix = self._compute_the_particles_distance_matrix()
        self.intersection_matrix = self._compute_the_intersection_matrix()

    def _compute_the_particles_distance_matrix(self):
        # TODO now works only for circles - otherwise return 0 matrix
        distance_matrix = np.array([
            [0 for i in range(self.number_of_particles)] for j in range(self.number_of_particles)
        ])
        if self.grain_type == "circle":
            self.grains_distance_matrix = pairwise_distances(
                [self.particles[k].germ for k in range(self.number_of_particles)]
            )
            radii_to_subtract = np.array([
                [
                    self.particles[k].grain.radius + self.particles[j].grain.radius
                    for j in range(self.number_of_particles)
                ] for k in range(self.number_of_particles)
            ])
            distance_matrix = np.where(
                self.grains_distance_matrix < radii_to_subtract, 0, self.grains_distance_matrix - radii_to_subtract
            )
        return distance_matrix


    def _compute_the_intersection_matrix(self):
        pass

    def plot_itself(self, show_germs: bool = True):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_aspect('equal', adjustable='box')
        for particle in self.particles:
            particle.grain.plot_2d(ax, facecolor="#E0AC69", linestyle="-", alpha=0.5, linewidth=1, edgecolor="black"
                                   # alpha=0.5,
            )
        if show_germs:
            for particle in self.particles:
                particle.germ.plot_2d(ax, c="red")
        plt.show()

