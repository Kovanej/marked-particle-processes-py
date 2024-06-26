
import numpy as np

GRAIN_VALID_TYPES = ["segment", "ball"]
GRAIN_TYPES_DIMENSIONS = {
    "ball": 2,
    "segment": 1
}
SAVE_PLOTS = True

MAX_SEGMENT_ANGLE = np.pi
MIN_SEGMENT_ANGLE = 0
MAX_SEGMENT_LENGTH = 0.2
MIN_SEGMENT_LENGTH = 0.1

MAX_CIRC_RAD = 0.1
MIN_CIRC_RAD = 0.05

POISSON_INTENSITY = 50
MARKED = True
SPACE_DIMENSION = 2
GRAIN_TYPE = "segment"

PERMUTATION_TEST_REPEAT_COUNT = 5000
SAVE_RESULTS_TO_CSV = False
PICKLE_RESULTS = False
PLOT_THE_P_VALUES = True

RADIUS_DISCRETE = "radius_discrete"
RADIUS_CONTINUOUS = "radius_continuous"
NND_DIST = "nearest_neighbour_distance"
MAX_SHARED_AREA_DISCRETE = "max_shared_area_discrete"
MAX_SHARED_AREA_CONTINUOUS = "max_shared_area_continuous"
INTERSECTION_COUNTING = "intersection_counting"

F_MARK_TYPES = [
    "product", "square", "first_mark",
]
WEIGHT_TYPES = {
    "ball": [
        "intersection", "shared_area"
    ],
    "segment": [
        "intersection", "angle"
    ],
}

F_MARK_COMBINATIONS = {
    k: [(f, v[i]) for f in F_MARK_TYPES for i in range(len(v))]
    for k, v in WEIGHT_TYPES.items()
}

PARTICLE_COLORS_CHOICE = [
    "#0089FF", "#FF0000", "#F6AF5B", "#4DA58C", "#90C26E", "#1B2784", "#000000",
    "#4DA58C", "#DE3657", "#FBD0C0",
    "#E16071", "#86007D", "#008018",
    # "#000000", "#FBD0C0",
    "#448F58", "#90C26E", "#DE3657", "#6A559B",
    "#613A6B", "#77ABB5", "#FBD0C0",
    "#A6BC3E", "#92A43B",
    "#A2E0DB", "#F55E55", "#FF857A", "#61A8E8", "#4FDDC3",
    "#FF6D74", "#FFC562", "#476088",
    "#1B8781", "#F61F51", "#4DF584", "#1CE784", "#04FF58", "#000000",
    "#448F58", "#90C26E", "#DE3657", "#6A559B",
    "#613A6B", "#77ABB5", "#FBD0C0",
    "#A6BC3E", "#92A43B",
    "#A2E0DB", "#F55E55", "#FF857A", "#61A8E8", "#4FDDC3",
    "#FF6D74", "#FFC562", "#476088",
    "#448F58", "#90C26E", "#DE3657", "#6A559B",
    "#613A6B", "#77ABB5", "#FBD0C0",
    "#A6BC3E", "#92A43B",
    "#A2E0DB", "#F55E55", "#FF857A", "#61A8E8", "#4FDDC3",
    "#FF6D74", "#FFC562", "#476088",
    "#448F58", "#90C26E", "#DE3657", "#6A559B",
    "#613A6B", "#77ABB5", "#FBD0C0",
    "#A6BC3E", "#92A43B",
    "#A2E0DB", "#F55E55", "#FF857A", "#61A8E8", "#4FDDC3",
    "#FF6D74", "#FFC562", "#476088",
]

ALPHA = 0.85
GRADIENT_SMOOTHNESS = 100

CONFIG_VALID_KEYS = [
    "process_type", "intensity", "space_dimension", "marking_type", "marking_parameters", "particles_parameters",
    "plot_realizations", "compute_f_mark_statistics", "f_mark_weights_and_statistics", "perform_permutation_test",
    "permutation_tests_parameters", "initial_seed", "number_of_realizations", "save_results"
]
CONFIG_NON_NULLABLE = ["process_type", "intensity"]
CONFIG_OPTIONAL_VALUES = {
    "space_dimension": 2,
    "marking_type": {
        "ball": [
            "radius_discrete", "radius_continuous", "nearest_neighbour_distance", "max_shared_area_continuous",
            "intersection_counting"
        ],
        "segment": [
            "angle_discrete", "angle_continuous", "length_discrete", "length_continuous", "nearest_neighbour_distance",
            "intersection_counting"
        ]
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
        #"intersection": ["product"] #,  "square", "first_mark"]
        #"shared_area": ["product",  "square", "first_mark"]

        "shared_area": ["first_mark"],
        #"shared_area": ["square"],
        #"shared_area": ["product"],
        #"intersection": ["square"],
        #"intersection": ["first_mark"]
        # "intersection": ["product"]
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
