
from typing import Optional, List, Any
import numpy as np
from skspatial.objects import Point, Points
import matplotlib.pyplot as plt


class PointProcess(object):

    def __init__(
            self,
            intensity: float,
            process_type: str,
            space_dimension: int = 2,
            points: Optional[Points] = None,
            generate_itself: bool = True,
            window_edge_start_point: float = 0,
            window_edge_end_point: float = 1,
    ):
        self.intensity = intensity
        self.space_dimension = space_dimension
        self.process_type = process_type
        self.window_edge_start_point = window_edge_start_point
        self.window_edge_end_point = window_edge_end_point
        if generate_itself:
            self.points = self._generate_itself()
        else:
            self.points = points

    def _generate_itself(self):
        # overriden in subclasses - TODO: Raise Error if not called in subclass
        pass

    def plot_itself(self):
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_aspect('equal', adjustable='box')
        # ax.set_xlim(left=0, right=1)
        # ax.set_ylim(bottom=0, top=1)
        self.points.plot_2d(ax)
        plt.show()


class PoissonPointProcess(PointProcess):

    def __init__(
            self,
            intensity: float,
            window_edge_start_point: float = 0,
            window_edge_end_point: float = 1,
            space_dimension: int = 2,
            points: Optional[List] = None,
            generate_itself: bool = True
    ):
        super().__init__(
            intensity=intensity,
            process_type="Poisson",
            space_dimension=space_dimension,
            points=points,
            generate_itself=generate_itself,
            window_edge_start_point=window_edge_start_point,
            window_edge_end_point=window_edge_end_point
        )

    def _generate_itself(self):  # -> Any(List[Point], Points):
        # how many points to generate in square [window_edge_start_point, window_edge_end_point]^(space_dimension)
        _number_of_points = np.random.poisson(
            lam=self.intensity * (self.window_edge_end_point - self.window_edge_start_point) ** self.space_dimension
        )
        _pts = np.array([
            [np.random.random_sample() * (
                    self.window_edge_end_point - self.window_edge_start_point
            ) + self.window_edge_start_point for _ in range(self.space_dimension)] for _ in range(_number_of_points)
        ]).astype(np.float32)
        return Points(_pts)

