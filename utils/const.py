
import numpy as np

GRAIN_VALID_TYPES = ["segment", "ball"]
GRAIN_TYPES_DIMENSIONS = {
    "ball": 2,
    "segment": 1
}
SAVE_PLOTS = False

MAX_SEGMENT_ANGLE = np.pi
MIN_SEGMENT_ANGLE = 0
MAX_SEGMENT_LENGTH = 0.2
MIN_SEGMENT_LENGTH = 0.1

MAX_CIRC_RAD = 0.1
MIN_CIRC_RAD = 0.05

POISSON_INTENSITY = 100
MARKED = True
MARKS_MODEL = "Lisa"
SPACE_DIMENSION = 2
GRAIN_TYPE = "segment"

PERMUTATION_TEST_REPEAT_COUNT = 500
SAVE_RESULTS_TO_CSV = False
PICKLE_RESULTS = False
PLOT_THE_P_VALUES = True

F_MARK_TYPES = [
    "product", "square", "first_mark",
]
WEIGHT_TYPES = {
    "ball": [
        "intersection", "shared_area", "distance"
    ],
    "segment": [
        "intersection", "distance", "angle"
    ],
}

F_MARK_COMBINATIONS = {
    k: [(f, v[i]) for f in F_MARK_TYPES for i in range(len(v))]
    for k, v in WEIGHT_TYPES.items()
}

# COLORS

# PARTICLE_COLORS_CHOICE = [
#                 "#FC6499", "#9FF4DF", "#F6F0A3"
#             ]
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
# PARTICLE_COLORS_CHOICE = [
#                 "#FC6499", "#9FF4DF", "#F6F0A3", "#6BCEEE"
#             ] + ["#22A7D9", "#2B2B2B", "#FCCF1C", "#FC6499"
#                  ] + ["#0075C2", "#FDE152", "#EEB336", "#017C81", "#254479", "#F87203"
#                       ] + ["#E0ADDD", "#FFF4D9", "#87D7D7", "#7488BB"]

# PARTICLE_COLORS_CHOICE = ["#448F58", "#90C26E", "#DE3657", "#6A559B", "#FF9443"]
PARTICLE_COLORS_CHOICE = [
    "#1B2784", "#F6AF5B", "#4DA58C", "#1B2784", "#448F58", "#90C26E", "#DE3657", "#6A559B", "#FF9443"
] + [
    "#000000"
]

# PARTICLE_COLORS_CHOICE = ["#FC6499", "#6BCEEE", "#87D7D7"]
ALPHA = 0.75

CONFIG_VALID_KEYS = [
    "process_type", "intensity", "space_dimension", "marking_type", "marking_parameters", "particles_parameters",
    "plot_realizations", "compute_f_mark_statistics", "f_mark_weights_and_statistics", "perform_permutation_test",
    "permutation_tests_parameters", "initial_seed", "number_of_realizations", "save_results"
]
CONFIG_NON_NULLABLE = ["process_type", "intensity"]
CONFIG_OPTIONAL_VALUES = {
    "space_dimension": 2,
    "marking_type": {
        "ball": ["radius_discrete", "radius_continuous"],
        "segment": ["angle_discrete", "angle_continuous", "length_discrete", "length_continuous"]
    },
    "particles_parameters": {
        "ball": {
          "max_radius": 0.15,
          "min_radius": 0.05
        },
        "segment": {
          "max_segment_length": 0.2,
          "min_segment_length": 0.1,
          "max_angle_in_degrees": 180,
          "min_angle_in_degrees": 0
        }
      },
    "marking_parameters": {
        "alphas": [0, 1]
      },
    "plot_realizations": True,
    "compute_f_mark_statistics": True,
    "f_mark_weights_and_statistics": {
        "intersection": ["product"]
    },
    "perform_permutation_test": False,
    "permutation_tests_parameters": {
        "number_of_permutations": 5000,
        "plot_histograms": False
    },
    "initial_seed": 0,
    "number_of_realizations": 1,
    "save_results": True
}
