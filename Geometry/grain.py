
from skspatial.objects import Vector, Line, Point


class Grain(object):

    def __init__(
            self,
            space_dimension: int,
            coordinates: list
    ):
        self.space_dimension = space_dimension
        self.coordinates = coordinates
        self.point = Point(array=coordinates)
