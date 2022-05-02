import matplotlib.pyplot as plt
import numpy as np
from skspatial.objects import Circle, Point

from Geometry.grain import Grain, Segment
from Processes.markings import Mark
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.particle_process import ParticleProcess
from Processes.point_process import PoissonPointProcess
import utils.const as const


if __name__ == '__main__':
    poisson_point_process = PoissonPointProcess(intensity=const.POISSON_INTENSITY)
    # TODO move to different script
    if const.GRAIN_TYPE == "segment":
        particles = []
        for k in range(len(poisson_point_process.points)):
            angle = np.pi * np.random.random(size=1)[0]
            length = (const.MAX_SEGMENT_LENGTH - const.MIN_SEGMENT_LENGTH) * np.random.random(size=1)[0] + const.MIN_SEGMENT_LENGTH
            mark = Mark(mark_type="discrete", mark_value=int(np.round(np.random.random(1)[0])))
            particles.append(
                Particle(
                    germ=Point(poisson_point_process.points[k]),
                    grain=Segment(
                        start_point=Point(poisson_point_process.points[k]) - [np.cos(angle) * length / 2, np.sin(angle) * length / 2],
                        end_point=Point(poisson_point_process.points[k]) + [np.cos(angle) * length / 2, np.sin(angle) * length / 2]
                    ),
                    grain_type="segment",
                    mark=mark
                )
            )
    particle_process = ParticleProcess(particles=particles, grain_type=const.GRAIN_TYPE)
    particle_process.compute_the_f_mark_characteristics()
    particle_process.plot_itself()
    brkpnt = "breakpoint here"
