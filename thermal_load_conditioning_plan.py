import numpy as np
import matplotlib.pyplot as plt


""" Define Functions """

def average_power_in_gun(CavFwdPwr, pulse_length_s, tau, rep_rate):
    """

    :param CavFwdPwr:
    :param tau:
    :return:
    """
    p_peak = CavFwdPwr * (1. - np.exp(-pulse_length_s/tau))

    p_avg = p_peak * rep_rate * pulse_length_s

    return p_avg

def cavity_forward_power_from_P_avg_in_cavity(p_avg, pulse_length_s, tau, rep_rate):
    """

    :param p_avg:
    :param pulse_length_s:
    :param tau:
    :param rep_rate:
    :return:
    """
    p_peak = p_avg / (rep_rate * pulse_length_s)

    CavFwdPwr = p_peak / (1. - np.exp(-pulse_length_s/tau))

    return CavFwdPwr

""" Define Variables """

pulse_length_vector = [
    250.e-9,
    500.e-9,
    750.e-9,
    1000.e-9,
    1250.e-9,
    1500.e-9,
    1750.e-9,
    2000.e-9,
    2250.e-9,
    2500.e-9,
    2750.e-9,
    3000.e-9,
    ]

pulse_length_vector_ns = [i*1.e9 for i in pulse_length_vector]

starting_pulse_length_peak_power = 1.e6
tau = 640.e-6
rep_rate = 100.
target_CavFwdPwr = 7.e6
target_CavFwdPwr_MW = target_CavFwdPwr/1.e6
ramp_rate_per_watt = 5000./120000.
current_pulse_length = 750.e-9
current_CavFwdPwr = 3.093e6

""" Calculate starting Cavity Forward Powers for each pulse length change """

starting_CavFwdPwrs = [starting_pulse_length_peak_power]
for idx in range(len(pulse_length_vector)-1):
    p_avg = average_power_in_gun(target_CavFwdPwr, pulse_length_vector[idx], tau, rep_rate)
    p_avg_90_percent = .9 * p_avg
    start_CavFwdPwr = cavity_forward_power_from_P_avg_in_cavity(p_avg_90_percent, pulse_length_vector[idx+1], tau, rep_rate)
    starting_CavFwdPwrs.append(start_CavFwdPwr)

for idx in range(len(starting_CavFwdPwrs)-1):
    print(f"pulse length of {pulse_length_vector[idx+1]*1.e6} ns, start at a cavity forward power of {starting_CavFwdPwrs[idx]/1.e6:1.3f} MW.")

starting_CavFwdPwrs_MW = [i/1.e6 for i in starting_CavFwdPwrs]
pulses_to_7MW = [(target_CavFwdPwr - i) / ramp_rate_per_watt for i in starting_CavFwdPwrs]

print("\n")
for idx in range(len(starting_CavFwdPwrs)-1):
    print(f"pulse length of {pulse_length_vector[idx+1]*1.e6} ns, {pulses_to_7MW[idx]/1.e6:1.3f} million pulses will have to be delivered.")

cumulative_pulses_vector = [pulses_to_7MW[i] + np.sum(pulses_to_7MW[:i]) for i in range(len(pulses_to_7MW))]
cumulative_pulses_vector.insert(0, 0.)
print(f"{cumulative_pulses_vector = }")
cumulative_pulses_vector_millions = [i/1.e6 for i in cumulative_pulses_vector]

""" calculate current position on plot """

current_pulse_length_index = pulse_length_vector.index(current_pulse_length)
current_cumulative_pulses = cumulative_pulses_vector[current_pulse_length_index] + (current_CavFwdPwr - starting_CavFwdPwrs[current_pulse_length_index]) / ramp_rate_per_watt
print(f"{current_cumulative_pulses = }")

print("\n")
for idx in range(len(cumulative_pulses_vector)-1):
    print(f"pulse length of {pulse_length_vector[idx]*1.e6} ns, the cumulative number of pulses delivered will be {cumulative_pulses_vector[idx]/1.e6:1.3f} million pulses.")

""" Calculate pulses and days until completion """

N_pulses_to_finish = cumulative_pulses_vector[-1] - current_cumulative_pulses
N_days_to_finish = N_pulses_to_finish/100./3600./24.
N_months_to_finish = N_days_to_finish/30.5
print(f"{N_pulses_to_finish = }")
print(f"{N_days_to_finish = }")
print(f"{N_months_to_finish = }")

""" Plot data """

for idx in range(len(cumulative_pulses_vector)-1):
    plt.plot([cumulative_pulses_vector_millions[idx],cumulative_pulses_vector_millions[idx+1]], [starting_CavFwdPwrs_MW[idx], target_CavFwdPwr_MW], label=r"$t_{pulse}$"f" = {int(pulse_length_vector_ns[idx])} ns. "r"$N_{pulses}$"f" = {cumulative_pulses_vector_millions[idx+1]:1.1f}")

plt.scatter(current_cumulative_pulses/1.e6, current_CavFwdPwr/1.e6, marker='x', color='r', s=40, label="Current position")
plt.ylim(-2., 7.5)
plt.legend(loc="lower right",fontsize=8, ncol=2,handleheight=2.4, labelspacing=0.05)
plt.xlabel("Millions of Pulses Delivered")
plt.ylabel("Cavity Forward Power [MW]")
plt.show()
