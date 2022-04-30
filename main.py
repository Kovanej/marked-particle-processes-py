import matplotlib.pyplot as plt
import numpy as np
from skspatial.objects import Circle, Point

from Geometry.grain import Grain, Segment
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.particle_process import ParticleProcess
from Processes.point_process import PoissonPointProcess

POISSON_INTENSITY = 300
MARKED = False
MARKS_MODEL = "Lisa"
# everything right now for 2d and circle grains, soon to be added segment grains
SPACE_DIMENSION = 2
GRAIN_TYPE = "segment"
MAX_CIRC_RAD = 0.2
MIN_CIRC_RAD = 0.1

#  RECENTLY MISSING FOR SEGMENTS - ADD DISTANCE MATRIX COMPUTATION & PLOTTING

if __name__ == '__main__':
    poisson_point_process = PoissonPointProcess(intensity=POISSON_INTENSITY)
    particles = [
        Particle(
            germ=Point(poisson_point_process.points[k]),
            grain=Segment(
                start_point=Point(poisson_point_process.points[k]) - Point([0.05, 0.05]),
                end_point=Point(poisson_point_process.points[k]) + Point([0.05, 0.05])
            ),
            grain_type="segment"
            # grain=Circle(
            #     Point(poisson_point_process.points[k]),
            #     (MAX_CIRC_RAD - MIN_CIRC_RAD) * (np.random.random(1)[0]) + MIN_CIRC_RAD)
        ) for k in range(len(poisson_point_process.points))
    ]
    particle_process = ParticleProcess(particles=particles, grain_type=GRAIN_TYPE)
    particle_process.plot_itself()
