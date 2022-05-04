
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
            generate_itself: bool = True
    ):
        self.intensity = intensity
        self.space_dimension = space_dimension
        self.process_type = process_type
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
        self.points.plot_2d(ax)
        plt.show()


class PoissonPointProcess(PointProcess):

    def __init__(
            self,
            intensity: int,
            space_dimension: int = 2,
            points: Optional[List] = None,
            generate_itself: bool = True
    ):
        super().__init__(
            intensity=intensity,
            process_type="Poisson",
            space_dimension=space_dimension,
            points=points,
            generate_itself=generate_itself
        )

    def _generate_itself(self): # -> Any(List[Point], Points):
        _number_of_points = np.random.poisson(lam=self.intensity)
        _pts = [[np.random.random(1)[0], np.random.random(1)[0]] for _ in range(_number_of_points)]
        return Points(_pts)

