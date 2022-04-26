
from typing import List

from Geometry.particle import Particle


class ParticleProcess(object):

    def __init__(
            self,
            particles: List[Particle]
    ):
        self.particles = particles
