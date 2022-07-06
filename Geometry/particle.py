
from typing import Union, Optional
from skspatial.objects import Point, Circle

import utils.const as const
from Geometry.grain import Segment
from Processes.markings import Mark


class Particle(object):

    def __init__(
            self,
            germ: Point,
            grain: Union[Circle, Segment],
            germ_inside_the_obs_window: bool,
            grain_type: str = "ball",
            mark: Optional[Mark] = None,
            space_dimension: int = 2
    ):
        if grain_type not in const.GRAIN_VALID_TYPES:
            raise ValueError(
                f"Invalid input of grain_type=={grain_type}. Please use one of the {const.GRAIN_VALID_TYPES}."
            )
        self.germ = germ
        self.germ_inside_the_obs_window = germ_inside_the_obs_window
        self.grain = grain
        self.space_dimension = space_dimension
        self.mark = mark
        self.grain_type = grain_type
        self.lebesgue_measure = self._compute_the_corresponding_measure()

    def _compute_the_corresponding_measure(self):
        if self.grain_type == "ball":
            return self.grain.area()
        elif self.grain_type == "segment":
            return self.grain.length
        else:
            raise ValueError(f"Unknown value for Particle.grain_type: {self.grain_type}")

