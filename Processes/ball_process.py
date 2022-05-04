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
