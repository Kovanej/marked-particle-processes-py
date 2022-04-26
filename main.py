
import matplotlib.pyplot as plt
from skspatial.objects import Circle, Point

from Geometry.grain import Grain
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.particle_process import ParticleProcess
from Processes.point_process import PoissonPointProcess

MARKED = False
# everything right now for 2d
SPACE_DIMENSION = 2
GRAIN_DIMENSION = 1
GRAIN_TYPE = "circle"


if __name__ == '__main__':
    poisson_point_process = PoissonPointProcess(intensity=20)
    particles = [
        Particle(
            germ=Point(poisson_point_process.points[k]),
            grain=Circle(Point(poisson_point_process.points[k]), 0.1)
        ) for k in range(len(poisson_point_process.points))
    ]
    particle_process = ParticleProcess(particles=particles)
    particle_process.plot_itself()
    a=1

