
import matplotlib.pyplot as plt

from Geometry.grain import Grain
from Plotting.plotting import plot_the_grains
from Processes.point_process import PoissonPointProcess

MARKED = False
# everything right now for 2d
SPACE_DIMENSION = 2
GRAIN_DIMENSION = 1
GRAIN_TYPE = "segment"


if __name__ == '__main__':
    poisson_process = PoissonPointProcess(intensity=100)
    # poisson_process.plot_itself()
    a=1

