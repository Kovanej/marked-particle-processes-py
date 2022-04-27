
import matplotlib.pyplot as plt
import numpy as np
from skspatial.objects import Circle, Point

from Geometry.grain import Grain
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.particle_process import ParticleProcess
from Processes.point_process import PoissonPointProcess

POISSON_INTENSITY = 8
MARKED = False
# everything right now for 2d
SPACE_DIMENSION = 2
GRAIN_DIMENSION = 1
GRAIN_TYPE = "circle"
MAX_CIRC_RAD = 0.1
MIN_CIRC_RAD = 0.03


if __name__ == '__main__':
    poisson_point_process = PoissonPointProcess(intensity=POISSON_INTENSITY)
    particles = [
        Particle(
            germ=Point(poisson_point_process.points[k]),
            grain=Circle(
                Point(poisson_point_process.points[k]),
                (MAX_CIRC_RAD - MIN_CIRC_RAD) * (np.random.random(1)[0]) + MIN_CIRC_RAD)
        ) for k in range(len(poisson_point_process.points))
    ]
    particle_process = ParticleProcess(particles=particles, grain_type=GRAIN_TYPE)
    particle_process.plot_itself()
    a=1

