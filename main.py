import matplotlib.pyplot as plt
import numpy as np
from skspatial.objects import Circle, Point

from Geometry.grain import Grain, Segment
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.particle_process import ParticleProcess
from Processes.point_process import PoissonPointProcess
import utils.const as const

POISSON_INTENSITY = 8000
MARKED = False
MARKS_MODEL = "Lisa"
SPACE_DIMENSION = 2
GRAIN_TYPE = "segment"
MAX_CIRC_RAD = 0.15
MIN_CIRC_RAD = 0.07


if __name__ == '__main__':
    poisson_point_process = PoissonPointProcess(intensity=POISSON_INTENSITY)
    # TODO move to different script
    if GRAIN_TYPE == "segment":
        particles = []
        for k in range(len(poisson_point_process.points)):
            angle = np.pi / 2 * np.random.random(size=1)[0]
            length = (const.MAX_SEGMENT_LENGTH - const.MIN_SEGMENT_LENGTH) * np.random.random(size=1)[0]
            particles.append(
                Particle(
                    germ=Point(poisson_point_process.points[k]),
                    grain=Segment(
                        start_point=Point(poisson_point_process.points[k]) - [np.cos(angle) * length / 2, np.sin(angle) * length / 2],
                        end_point=Point(poisson_point_process.points[k]) + [np.cos(angle) * length / 2, np.sin(angle) * length / 2]
                    ),
                    grain_type="segment"
                )
            )

    # particles = [
    #     Particle(
    #         germ=Point(poisson_point_process.points[k]),
    #         grain=Segment(
    #             start_point=Point(poisson_point_process.points[k]) - Point([0.05, 0.05]),
    #             end_point=Point(poisson_point_process.points[k]) + Point([0.05, 0.05])
    #         ),
    #         grain_type="segment"
    #         # grain=Circle(
    #         #     Point(poisson_point_process.points[k]),
    #         #     (MAX_CIRC_RAD - MIN_CIRC_RAD) * (np.random.random(1)[0]) + MIN_CIRC_RAD)
    #     ) for k in range(len(poisson_point_process.points))
    # ]
    particle_process = ParticleProcess(particles=particles, grain_type=GRAIN_TYPE)
    particle_process.plot_itself()
    a=1
