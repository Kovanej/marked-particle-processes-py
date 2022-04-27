
from skspatial.objects import Point, Circle

from Geometry.germ import Germ


class Particle(object):

    def __init__(
            self,
            germ: Point,
            # TODO extend - now testing for circles
            grain: Circle,
            grain_type: str = "circle",
            mark=None
    ):
        self.germ = germ
        self.grain = grain
        self.mark = mark
        self.grain_type = grain_type
        self.lebesgue_measure = self._compute_the_corresponding_measure()

    def _compute_the_corresponding_measure(self):
        if self.grain_type == "circle":
            return self.grain.area()
        else:
            raise ValueError(f"Unknown value for Particle.grain_type: {self.grain_type}")


