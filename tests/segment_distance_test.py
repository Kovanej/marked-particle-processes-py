
from skspatial.objects import Point

from Geometry.grain import Grain, Segment
from Plotting.plotting import plot_the_grains
from Geometry.particle import Particle
from Processes.particle_process import ParticleProcess


# ALL FOLLOWING CHECKED FOR CORRECT VALUES
# START_AND_END_POINTS = [
#     (Point([0, 0]), Point([1, 0])),
#     (Point([0, 1/2]), Point([0, 1/4])),
#     (Point([1, 1]), Point([2, 2])),
#     (Point([0, 1]), Point([-1, 0])),
#     (Point([0, 0]), Point([1, 1])),
#     (Point([1, 0]), Point([2, 1]))
# ]

# ALSO TESTED FOR A COUPLE OF SEGMENTS IN R^3
START_AND_END_POINTS = [
    (Point([0, 0, 0]), Point([0, 0, -1])),
    (Point([-1, 2, 0]), Point([1, 2, 0])),
    (Point([-1, 0, -2]), Point([2, 0, -2])),
    (Point([-1, 0, 2]), Point([2, 0, 2])),
]

particles = [
        Particle(
            germ=point[0] / 2 + point[1] / 2,
            grain=Segment(
                start_point=point[0],
                end_point=point[1]
            ),
            grain_type="segment"
        ) for point in START_AND_END_POINTS
    ]
particle_process = ParticleProcess(particles=particles, grain_type="segment")
a=1
