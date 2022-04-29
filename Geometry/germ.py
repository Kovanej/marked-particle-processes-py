
from skspatial.objects import Vector, Line, Point, Circle


class Germ(object):

    def __init__(
            self,
            space_dimension: int
    ):
        self.space_dimension = space_dimension

