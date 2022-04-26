
from typing import List

from skspatial.objects import Vector, Line, Point
from skspatial.plotting import plot_2d

import matplotlib.pyplot as plt

from Geometry.grain import Grain


def plot_the_grains(grains: List[Grain], add=False):
    x_coords = [grain.point[0] for grain in grains]
    y_coords = [grain.point[1] for grain in grains]
    plt.scatter(x_coords, y_coords)
    plt.scatter(y_coords, x_coords)
    plt.show()


