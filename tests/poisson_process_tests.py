
import matplotlib.pyplot as plt
from skspatial.objects import Point, Points

from Processes.point_process import PoissonPointProcess

fig = plt.figure()
ax = fig.add_subplot()
ax.set_aspect('equal', adjustable='box')
poisson_point_process = PoissonPointProcess(
    intensity=1, space_dimension=2, window_edge_start_point=-0.5, window_edge_end_point=1.5
)
poisson_point_process.points.plot_2d(ax_2d=ax)
plt.show()
brkpnt = "breakpoint here"
