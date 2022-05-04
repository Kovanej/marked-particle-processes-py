
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
            space_dimension: int
    ):
        self.germ_intensity = germ_intensity
        self.particles = particles
        self.number_of_particles = len(self.particles)
        self.grain_type = grain_type
        self.space_dimension = space_dimension
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

    def _plot_ball_particles(self, ax):
        for particle in self.particles:
            facecolor, alpha = self._choose_face_color()
            edgecolor = self._choose_edge_color()
            particle.grain.plot_2d(
                ax, facecolor=facecolor, linestyle="-", alpha=alpha, linewidth=1, edgecolor=edgecolor,
                # alpha=0.5,
            )

    def _plot_segment_particles(self, fig, ax):
        if self.space_dimension == 3:
            ax = fig.add_subplot(111, projection='3d')
        for particle in self.particles:
            # col, alpha = self._choose_face_color()
            if particle.mark is not None:
                if particle.mark.mark_value == 0:
                    col, alpha = "#FEC500", 1
                elif particle.mark.mark_value == 1:
                    col, alpha = "#003271", 1
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
            self._plot_ball_particles(ax=ax)
        if self.grain_type == "segment":
            self._plot_segment_particles(fig=fig, ax=ax)
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
        return np.random.choice([
                "#F9EFB4", "#F4BE9A", "#BA6191", "#40008C", "#000032"
            ])

    def _choose_face_color(self, particle=None):
        # TODO later adjust for marks
        alpha = 1
        col = np.random.choice(const.PARTICLE_COLORS_CHOICE)
        return col, alpha

    def _choose_germ_color(self):
        return "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]), np.random.random(1)[0]
