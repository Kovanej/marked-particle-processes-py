
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils.config_parser import ConfigParser
import Processes.ball_process as bp


def s_wA_fm_ball_radius_discrete(t: float, gamma: float, alpha: float, a: float, b: float, p: float = 1 / 2):
    s_wA_fm_alpha_0 = gamma ** 2 * p * (np.pi ** 3) * (
            (b ** 2 + a * b + a ** 2) / (3 * (b - a))
    ) * (
            t ** 2 * (b ** 3 - a ** 3) / 3 +
            t * (b ** 4 - a ** 4) / 2 +
            (b ** 5 - a ** 5) / 5
    )
    s_wA_fm_alpha_1 = gamma ** 2 * (np.pi ** 3) * (
            b ** 2 + a * b + a ** 2
    ) / (3 * ((b - a) ** 2)) * (
            t ** 2 * ((b ** 4 - a ** 4) / 4 - a * (b ** 3 - a ** 3) / 3) +
            t * (2 * (b ** 5 - a ** 5) / 5 - 2 * a * (b ** 4 - a ** 4) / 4) +
            (b ** 6 - a ** 6) / 6 - a * (b ** 5 - a ** 5) / 5
    )
    return alpha * s_wA_fm_alpha_1 + (1 - alpha) * s_wA_fm_alpha_0


def s_WI_fp_ball_radius_discrete(t: float, gamma: float, p: float, a: float, b: float):
    return (gamma * p * np.pi) ** 2 * (1 / 180) * (
            (210 * b ** 2 + 300 * a * b + 210 * a ** 2) * t ** 2 +
            (270 * b ** 3 + 450 * a ** 2 * b + 450 * a * b ** 2 + 210 * a ** 3) * t +
            (101 * b ** 4 + 166 * a * b ** 3 + 186 * a ** 2 * b ** 2 + 166 * a ** 3 * b + 101 * a ** 4)
    )


points_to_eval = sorted(np.round(np.arange(0.01, 1, 0.01, dtype=float), 2))
gamma = 50
a = 0.05
b = 0.15
p = 1/2
alphas_to_plot = [0, 0.25, 0.5, 0.75, 1]
alphas_to_plot = [0]

outputs_to_plot = {
    alpha: [s_wA_fm_ball_radius_discrete(t=t, gamma=gamma, alpha=alpha, a=a, b=b, p=p) for t in points_to_eval]
    for alpha in alphas_to_plot
}


f_type = 'product'
w_type = 'shared_area'

with open("../config.json", "r") as json_data:
    config_json = json.loads(json_data.read())

config_parser = ConfigParser(config=config_json)

envelope_count = config_json["permutation_tests_parameters"]["number_of_permutations"]
init_seed = config_json["initial_seed"]

dfs = {}

for _ in range(100, 131):
    seed = init_seed + _
    config_parser.initialize_the_processes(seed=seed)
    result_saver = config_parser.return_the_result_saver(seed=seed)
    dfs[_] = result_saver.results_all_df

for k, v in dfs.items():
    plot_v = v.loc[(v['f-Mark Type'] == f_type) & (v['Weight Type'] == w_type)][['Input Value', 'PWFCF Value']].sort_values('Input Value')
    plt.plot(plot_v['Input Value'], plot_v['PWFCF Value'], marker=',')

# Show the plot


for k, v in outputs_to_plot.items():
    plt.plot(points_to_eval, v, marker='o', label=f"alpha={k}")

plt.legend()

# Add axis labels and a title
plt.xlabel('t')
plt.ylabel('S_{w_A, f_m}(t)')
plt.title('S_{w_A, f_m}(t) for Model 1 with respect to alpha.')

#plt.grid()

plt.show()



