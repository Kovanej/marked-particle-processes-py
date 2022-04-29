
import numpy as np
from typing import Optional
from skspatial.objects import Vector, Line, Point


class Grain(object):

    def __init__(
            self,
            space_dimension: int,
            circumcenter: Point,
    ):
        self.space_dimension = space_dimension
        self.circumcenter = circumcenter


class Segment(Grain):

    def __init__(
            self,
            space_dimension: int = 2,
            start_point: Optional[Point] = None,
            end_point: Optional[Point] = None,
            circumcenter: Optional[Point] = None,
            angle: Optional[float] = None,
            length: Optional[float] = None
    ):
        # right now two possible parametrizatons - a) start_point & end_point, b) start_point, angle, length
        _valid_input = self._validate_the_input(circumcenter, angle, start_point, end_point, length)
        if not _valid_input:
            raise Exception(f"Not a valid input for a class Segment.")
        super().__init__(
            space_dimension=space_dimension, circumcenter=circumcenter
        )
        self.angle = angle
        self.start_point = start_point
        self.end_point = end_point
        self.length = length
        self._compute_the_missing_parameters()

    def _compute_the_missing_parameters(self):
        if self.end_point is None:
            # assuming 2-dim for now
            x_1 = self.start_point[0] + self.length * np.cos(self.angle)
            y_1 = self.start_point[0] + self.length * np.sin(self.angle)
            self.end_point = Point(np.array([x_1, y_1]))
        if self.angle is None:
            # TODO now only for 2d
            _vector = Vector(self.end_point - self.start_point)
            if _vector[1] < 0:
                _vector = - _vector
            self.angle = Vector([1, 0]).angle_between(_vector)
        if self.length is None:
            self.length = self.start_point.distance_point(self.end_point)
        self.vector = Vector(self.end_point - self.start_point)
        if self.circumcenter is None:
            self.circumcenter = (self.start_point + self.end_point) / 2

    @staticmethod
    def _validate_the_input(circumcenter, angle, start_point, end_point, length):
        if start_point is not None and end_point is not None:
            return True
        if start_point is not None and angle is not None and length is not None:
            return True
        else:
            return False

