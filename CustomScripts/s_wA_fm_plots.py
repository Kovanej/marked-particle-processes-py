
import json
import matplotlib.pyplot as plt
import numpy as np

import utils.const as const
from utils.config_parser import ConfigParser


def s_WI_f_uniform_segment(t: float, gamma: float, alpha: float, a: float, b: float, p: float = 1 / 2,
                                 f_type="first_mark"):
    if f_type == "product":
        c_f = p ** 2
    elif f_type == "square":
        c_f = p * (1 - p)
    else:
        c_f = p

    val = gamma ** 2 * c_f * (b + a) * (
        t ** 2 * (b + a) / 2 +
        t * 2 / (3 * np.pi) * (b ** 2 + a * b + a ** 2)
    )
    return val


def s_wI_fp_angle_discrete(t: float, gamma: float, alpha: float, a: float, b: float):

    if alpha == 1:
        val = gamma ** 2 * (2 * b ** 2 - a * b - a ** 2) / (18 * (b - a) ** 2) * (
            t ** 2 * (2 * b ** 2 - a * b - a ** 2) +
            t * (3 * b ** 3 - a * b ** 2 - a ** 2 * b - a ** 3) / np.pi
        )
    else:
        val = s_WI_f_uniform_segment(t=t, gamma=gamma, alpha=0, a=a, b=b, f_type="product")
    return val


def s_wA_f_ball_radius_discrete(t: float, gamma: float, alpha: float, a: float, b: float, p: float = 1 / 2,
                                 f_type="first_mark"):
    if f_type == "product":
        c_f = p ** 2
    elif f_type == "square":
        c_f = p * (1 - p)
    else:
        c_f = p
    s_wA_fm_alpha_0 = gamma ** 2 * c_f * (np.pi ** 3) * (
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


with open("../config.json", "r") as json_data:
    config_json = json.loads(json_data.read())

config_parser = ConfigParser(config=config_json)

gamma = config_parser.intensity
alphas_to_plot = config_parser.marking_parameters['alphas']
if config_parser.process_type == 'ball':
    a = config_parser.particles_parameters['ball']['min_radius']
    b = config_parser.particles_parameters['ball']['max_radius']
else:
    a = config_parser.particles_parameters['segment']['min_segment_length']
    b = config_parser.particles_parameters['segment']['max_segment_length']

p = 1/2
col_ind = 0

# DOUBLE CHECK, WHAT I ACTUALLY WANT
plot_simulations = True
f_type = 'product'
w_type = 'intersection'

outputs_to_plot = {
    alpha: [s_wA_f_ball_radius_discrete(t=t, gamma=gamma, alpha=alpha, a=a, b=b, p=p, f_type=f_type) for t in points_to_eval]
    for alpha in alphas_to_plot
}

outputs_to_plot = {
    alpha: [s_WI_f_uniform_segment(t=t, gamma=gamma, alpha=alpha, a=a, b=b, p=p, f_type=f_type) for t in points_to_eval]
    for alpha in alphas_to_plot
}

outputs_to_plot = {
    alpha: [s_wI_fp_angle_discrete(t=t, gamma=gamma, alpha=alpha, a=a, b=b) for t in points_to_eval]
    for alpha in alphas_to_plot
}

fig = plt.figure()

axes = fig.add_axes([0,0,1,1])

if plot_simulations:
    envelope_count = config_json["permutation_tests_parameters"]["number_of_permutations"]
    init_seed = config_json["initial_seed"]

    dfs = {}

    for _ in range(69, 110):
        seed = init_seed + _
        config_parser.initialize_the_processes(seed=seed)
        result_saver = config_parser.return_the_result_saver(seed=seed)
        dfs[_] = result_saver.results_all_df

    cols = {}

    for k, v in dfs.items():
        v_gb = v.groupby(['Model'])
        for g_n, g_df in v_gb:
            if g_n not in cols.keys():
                cols[g_n] = const.PARTICLE_COLORS_CHOICE[col_ind]
                col_ind += 1
            color = const.PARTICLE_COLORS_CHOICE[col_ind]
            plot_v = g_df.loc[(g_df['f-Mark Type'] == f_type) & (g_df['Weight Type'] == w_type)][
                ['Input Value', 'PWFCF Value']].sort_values('Input Value')
            axes.plot(plot_v['Input Value'], plot_v['PWFCF Value'], marker=',',
                     color=cols[g_n], alpha=0.3)

# Show the plot

# fuj
col_ind = 1

for k, v in outputs_to_plot.items():
    plt.plot(points_to_eval, v, marker='o', label=f"alpha={k}", color=const.PARTICLE_COLORS_CHOICE[col_ind])
    col_ind -= 1

plt.legend()

# Add axis labels and a title
plt.xlabel('t')
plt.ylabel('S_{w, f}(t)')
plt.title('S_{w, f}(t) for Model 1.')

plt.grid()

plt.show()



