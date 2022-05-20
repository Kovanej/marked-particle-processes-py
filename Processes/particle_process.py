
from datetime import datetime
import logging
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
        logging.info(f"Particle process class initialized.")
        self.germ_intensity = germ_intensity
        self.particles = particles
        self.number_of_particles = len(self.particles)
        self.grain_type = grain_type
        self.space_dimension = space_dimension
        self.marks = np.array([p.mark.mark_value for p in self.particles])
        self.marks_product = self.marks[..., None] * self.marks[None, ...]
        # compute the grains distance
        self.germs_distance_matrix = self._compute_the_germs_distance_matrix()
        # compute particles_null_model distance (inf {||x_i - x_j||: x_i \in \Xi_i, x_j \in \Xi_j})
        self.particles_distance_matrix = self._compute_the_particles_distance_matrix()
        # particles_null_model distance == 0 -> intersected, otherwise not intersected
        self.particles_intersection_matrix = self._compute_the_particles_intersection_matrix()
        # compute the pairwise shared corresponding Lebesgue measure
        # ("ball": shared areas, "segment": same as intersection matrix ...)
        self.shared_corresponding_measure_matrix = self._compute_the_shared_corresponding_measure_matrix()
        # if needed following attributes are computed via executing ParticleProcess.compute_f_mark_statistics
        self.f_mark_normalization_constant = None
        self.f_mark_intersection_correlation = None

    def _compute_the_particles_intersection_matrix(self):
        logging.info(f"{datetime.now()} :Particle intersection matrix computation start.")
        intersection_matrix = np.where(self.particles_distance_matrix == 0, 1, 0)
        logging.info(f"{datetime.now()} :Particle intersection matrix computation end.")
        return intersection_matrix

    def _compute_the_shared_corresponding_measure_matrix(self):
        # overridden for each subclass - if not specified, it cannot be computed
        raise NotImplementedError(
            f"Method called from the instance of general ParticleProcess class. "
            f"Please use the instance of some of the ParticleProcess subclasses (SegmentProcess, BallProcess)."
        )

    def _compute_the_germs_distance_matrix(self):
        logging.info(f"{datetime.now()} :Germs distance computation start.")
        grains_distance_matrix = pairwise_distances(
                [self.particles[k].germ for k in range(self.number_of_particles)]
            )
        logging.info(f"{datetime.now()} :Germs distance computation end.")
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

    def compute_the_f_mark_characteristics(
            self,
            f_type: str = "product"
    ):
        self.f_mark_normalization_constant = self._compute_the_f_mark_normalization_constant(f_type=f_type)
        self.f_mark_intersection_correlation = self._compute_the_f_mark_intersection_correlation(f_type=f_type)

    def _compute_the_f_mark_intersection_correlation(self, f_type: str):
        intersection_zero_diagonal = self.particles_intersection_matrix.copy()
        intersection_zero_diagonal[np.diag_indices(self.number_of_particles)] = 0
        f_intersection_nn = (self.marks_product * intersection_zero_diagonal).sum() / 2
        self.f_intersection_nn = f_intersection_nn
        return (2 * f_intersection_nn) / (self.f_mark_normalization_constant * (self.germ_intensity ** 2))

    def _compute_the_f_mark_normalization_constant(self, f_type: str):
        # estimate the E[f(M_i, M_j)] - evaluate it on each pair of n points and divide by (n choose 2)
        pairs_count = (self.number_of_particles * (self.number_of_particles - 1)) / 2
        marks_product_zero_diagonal = self.marks_product.copy()
        marks_product_zero_diagonal[np.diag_indices(self.number_of_particles)] = 0
        c_f_nn = marks_product_zero_diagonal.sum() / 2
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
        alpha = const.ALPHA
        col = np.random.choice(const.PARTICLE_COLORS_CHOICE)
        return col, alpha

    def _choose_germ_color(self):
        return "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]), np.random.random(1)[0]
