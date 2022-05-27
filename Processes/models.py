

class Model(object):

    def __init__(self, grain_type: str, germ_intensity: float):
        self.grain_type = grain_type
        self.germ_intensity = germ_intensity

    def simulate_itself(self):
        raise NotImplementedError(
            f"General class Model called. Please call only the specific Model subclass."
        )


class SegmentModel(Model):

    def __init__(self, germ_intensity: float):
        super().__init__(grain_type="segment", germ_intensity=germ_intensity)

    def simulate_itself(self):
        pass


class BallModel(Model):

    def __init__(self, germ_intensity: float):
        super().__init__(grain_type="ball", germ_intensity=germ_intensity)

    def simulate_itself(self):
        pass

