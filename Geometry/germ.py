
from skspatial.objects import Vector, Line, Point, Circle


class Germ(object):

    def __init__(
            self,
            space_dimension: int,
            circumcenter: Point
    ):
        self.space_dimension = space_dimension
        self.circumcenter = circumcenter


# class Circle(Germ):
#
#     def __init__(
#             self,
#             space_dimension: int,
#             circumcenter: Point,
#             radius: float
#     ):
#         super().__init__(space_dimension=space_dimension, circumcenter=circumcenter)
#         self.radius = radius
