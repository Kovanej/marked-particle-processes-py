
from skspatial.objects import Point, Circle

from Geometry.germ import Germ


class Particle(object):

    def __init__(
            self,
            germ: Point,
            # TODO extend - now testing for circles
            grain: Circle
    ):
        self.germ = germ
        self.grain = grain