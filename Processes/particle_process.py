
import matplotlib.pyplot as plt
from typing import List

from Geometry.particle import Particle


class ParticleProcess(object):

    def __init__(
            self,
            particles: List[Particle]
    ):
        self.particles = particles

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

