
import matplotlib.pyplot as plt
from skspatial.objects import Point, Points

from Processes.point_process import PoissonPointProcess

poisson_point_process = PoissonPointProcess(
    intensity=100, space_dimension=2, window_edge_start_point=-0.5, window_edge_end_point=1.5
)
fig = plt.figure()
ax = fig.add_subplot()
ax.set_aspect('equal', adjustable='box')
ax.set_xlim(left=0, right=1)
ax.set_ylim(bottom=0, top=1)
poisson_point_process.points.plot_2d(ax_2d=ax, c="#000000", marker=".")
plt.show()
brkpnt = "breakpoint here"
