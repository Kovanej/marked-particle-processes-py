
GRAIN_VALID_TYPES = ["segment", "ball"]
GRAIN_TYPES_DIMENSIONS = {
    "ball": 2,
    "segment": 1
}
SAVE_PLOTS = False

MAX_SEGMENT_LENGTH = 1
MIN_SEGMENT_LENGTH = 0

MAX_CIRC_RAD = 0.15
MIN_CIRC_RAD = 0.05

POISSON_INTENSITY = 100
MARKED = True
MARKS_MODEL = "Lisa"
SPACE_DIMENSION = 2
GRAIN_TYPE = "segment"

COMPUTE_SEGMENT_DISTANCES = False

F_MARK_TYPES = [
    "product", "square",
    # "first_mark",
]
WEIGHT_TYPES = [
    "intersection", "shared_area",
    # "distance"
]

F_MARK_COMBINATIONS = [
    (f, w) for f in F_MARK_TYPES for w in WEIGHT_TYPES
]

# COLORS

# PARTICLE_COLORS_CHOICE = [
#                 "#F9EFB4", "#F4BE9A", "#BA6191", "#40008C", "#000032"
#             ]
# PARTICLE_COLORS_CHOICE = [
#                 "#E0ADDD", "#FFF4D9", "#87D7D7", "#7488BB"
#             ]
# PARTICLE_COLORS_CHOICE = [
#                 "#0075C2", "#FDE152", "#EEB336", "#017C81", "#254479", "#F87203"
#             ]
# PARTICLE_COLORS_CHOICE = [
#                 "#22A7D9", "#2B2B2B", "#FCCF1C", "#FC6499"
#             ]
PARTICLE_COLORS_CHOICE = [
                "#FC6499", "#9FF4DF", "#F6F0A3", "#6BCEEE"
            ] + ["#22A7D9", "#2B2B2B", "#FCCF1C", "#FC6499"
                 ] + ["#0075C2", "#FDE152", "#EEB336", "#017C81", "#254479", "#F87203"
                      ] + ["#E0ADDD", "#FFF4D9", "#87D7D7", "#7488BB"]


# PARTICLE_COLORS_CHOICE = ["#FC6499", "#6BCEEE", "#87D7D7"]
ALPHA = 0.65
