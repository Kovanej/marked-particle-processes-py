
import matplotlib.pyplot as plt
import numpy as np


def simple_pwfcf_for_uniform_distribution(value: float, lam: float, p: float, a: float, b: float):
    return (lam * p * np.pi) ** 2 * (1 / 180) * (
        (210 * b ** 2 + 300 * a * b + 210 * a ** 2) * value ** 2 +
        (270 * b ** 3 + 450 * a ** 2 * b + 450 * a * b ** 2 + 210 * a ** 3) * value +
        (101 * b ** 4 + 166 * a * b ** 3 + 186 * a ** 2 * b ** 2 + 166 * a ** 3 * b + 101 * a ** 4)
    )


points_to_eval = np.round(np.arange(0.02, 1, 0.02, dtype=float), 2)

# Create a list of sorted inputs and outputs
inputs = sorted(points_to_eval)

outputs_theoretical = {
    "lamE40_pE0D5_aE0D2_bE0D2": [simple_pwfcf_for_uniform_distribution(value=t, lam=40, p=0.5, a=0.2, b=0.2) for t in inputs],
    "lamE40_pE0D5_aE0D2_bE0D3": [simple_pwfcf_for_uniform_distribution(value=t, lam=40, p=0.5, a=0.2, b=0.3) for t in inputs],
    "lamE40_pE0D5_aE0D2_bE0D5": [simple_pwfcf_for_uniform_distribution(value=t, lam=40, p=0.5, a=0.2, b=0.5) for t in inputs],
    "lamE40_pE0D3_aE0D2_bE0D5": [simple_pwfcf_for_uniform_distribution(value=t, lam=40, p=0.3, a=0.2, b=0.5) for t in inputs],
    "lamE40_pE0D7_aE0D2_bE0D5": [simple_pwfcf_for_uniform_distribution(value=t, lam=40, p=0.7, a=0.2, b=0.5) for t in inputs],
    "lamE40_pE0D5_aE0D5_bE0D8": [simple_pwfcf_for_uniform_distribution(value=t, lam=40, p=0.5, a=0.2, b=0.7) for t in inputs]
}

# Create a line plot
for k, v in outputs_theoretical.items():
    plt.plot(inputs, v, marker='', label=f"{k.replace('E', '=').replace('D', '.').replace('_', ', ')}")

plt.legend()

# Add axis labels and a title
plt.xlabel('t')
plt.ylabel('S_{w, f}(t)')
plt.title('PWFCF for unif. dist. Boolean ball process')

# Show the plot
plt.show()


break_point_var=1

