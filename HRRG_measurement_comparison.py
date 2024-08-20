import numpy as np
import argparse
import sys, os, csv, pickle
from matplotlib import pyplot as plt
import pandas as pd
from snps import SNPs
pmm_addr = r'C:\Users\zup98752\PycharmProjects\PhD'
sys.path.insert(1,pmm_addr)
import HRRG_Master_Module as hmm
import skrf as rf

read_in_pickle = True
C16_APR23_colour = "b"
C20_APR23_colour = "g"
C8_JAN24_colour = "r"

savepath_APR23 = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20230418\analysis'
snp_addr_APR23 = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20230418\data'
pkl_addr_APR23 = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20230418'
fnames_APR23 = os.listdir(snp_addr_APR23)
fname_APR23 = fnames_APR23[0]
fname_and_addr_APR23 = f'{snp_addr_APR23}\\{fname_APR23}'
print(f'{fnames_APR23 = }')

savepath_JAN24 = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20240123\analysis'
snp_addr_JAN24 = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20240123\data'
pkl_addr_JAN24 = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20240123'
fnames_JAN24 = os.listdir(snp_addr_JAN24)
fname_JAN24 = fnames_JAN24[0]
fname_and_addr_JAN24 = f'{snp_addr_JAN24}\\{fname_JAN24}'
print(f'{fnames_JAN24 = }')


savepath_COMP = r"C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_comparison\analysis"

""" Read in data from previously saved pickles """

with open(f'{pkl_addr_APR23}\\data_dict.pkl', 'rb') as handle:
    data_dict_APR23 = pickle.load(handle)

with open(f'{pkl_addr_JAN24}\\data_dict.pkl', 'rb') as handle:
    data_dict_JAN24 = pickle.load(handle)

freq_C16_APR23 = []
max_gamma_C16_APR23 = []
Q_0_C16_APR23 = []
beta_C16_APR23 = []
bandwidth_C16_APR23 = []

freq_C20_APR23 = []
max_gamma_C20_APR23 = []
Q_0_C20_APR23 = []
beta_C20_APR23 = []
bandwidth_C20_APR23 = []
print(f'\nAPR23 data')
for key, val in data_dict_APR23.items():
    print(f'\n{key}')
    for key_2, val_2 in val.items():
        try:
            print(f'{key_2}: length = {len(val_2)}')
        except:
            print(f'{key_2}: {val_2}')
        if "C16" in str(data_dict_APR23[key]["name"]): # and "removed" in str(data_dict_APR23[key]["name"]):
            if str(key_2) == "res_freq":
                freq_C16_APR23.append(val_2)
            elif str(key_2) == "max_gamma_s21":
                max_gamma_C16_APR23.append(val_2)
            elif str(key_2) == "Q_0":
                Q_0_C16_APR23.append(val_2)
            elif str(key_2) == "beta":
                beta_C16_APR23.append(val_2)
            elif str(key_2) == "bandwidth":
                bandwidth_C16_APR23.append(val_2)
            else:
                pass

        if "C20" in str(data_dict_APR23[key]["name"]): # and "removed" in str(data_dict_APR23[key]["name"]):
            if str(key_2) == "res_freq":
                freq_C20_APR23.append(val_2)
            elif str(key_2) == "max_gamma_s21":
                max_gamma_C20_APR23.append(val_2)
            elif str(key_2) == "Q_0":
                Q_0_C20_APR23.append(val_2)
            elif str(key_2) == "beta":
                beta_C20_APR23.append(val_2)
            elif str(key_2) == "bandwidth":
                bandwidth_C20_APR23.append(val_2)
            else:
                pass
            
        

freq_JAN24 = []
max_gamma_JAN24 = []
Q_0_JAN24 = []
beta_JAN24 = []
bandwidth_C8_JAN24 = []
print(f'\nJAN24 data')
for key, val in data_dict_JAN24.items():
    print(f'\n{key}')
    for key_2, val_2 in val.items():
        try:
            print(f'{key_2}: length = {len(val_2)}')
        except:
            print(f'{key_2}: {val_2}')

        if str(key_2) == "res_freq":
            freq_JAN24.append(val_2)
        elif str(key_2) == "max_gamma_s21":
            max_gamma_JAN24.append(val_2)
        elif str(key_2) == "Q_0":
            Q_0_JAN24.append(val_2)
        elif str(key_2) == "beta":
            beta_JAN24.append(val_2)
        elif str(key_2) == "bandwidth":
            bandwidth_C8_JAN24.append(val_2)
        else:
            pass

# design frequency of 2.9985 GHz as measured under experimental conditions
# shift frequencies to design and centre on zero
measured_freq_Hz = 2999982788.
design_freq_Hz = 2998500000.
delta_freq = measured_freq_Hz - design_freq_Hz
print(f'{delta_freq = }')
freq_C16_APR23 = [((i-delta_freq)-design_freq_Hz)/1.e3 for i in freq_C16_APR23]
freq_C20_APR23 = [((i-delta_freq)-design_freq_Hz)/1.e3 for i in freq_C20_APR23]
freq_JAN24 = [(i-design_freq_Hz)/1.e3 for i in freq_JAN24]


# remove data taken without cathode in
freq_JAN24 = freq_JAN24[:-1]
max_gamma_JAN24 = max_gamma_JAN24[:-1]
Q_0_JAN24 = Q_0_JAN24[:-1]
beta_JAN24 = beta_JAN24[:-1]
# freq max gamma

plt.scatter(freq_C16_APR23, max_gamma_C16_APR23, color=C16_APR23_colour, label="C16  APR23")
plt.scatter(freq_C20_APR23, max_gamma_C20_APR23, color=C20_APR23_colour, label="C20  APR23")
plt.scatter(freq_JAN24, max_gamma_JAN24, color=C8_JAN24_colour, label="C8    JAN24")
plt.xlabel(r'$\Delta$''f [kHz]')
plt.ylabel('Max 'r'$\Gamma_{S21}$'' [dB]')
plt.legend(loc='lower left')
plt.savefig(f'{savepath_COMP}\\freq_maxGamma_all.png')
plt.close('all')

plt.scatter(freq_JAN24, max_gamma_JAN24, color=C8_JAN24_colour, label="C8    JAN24")
plt.xlabel(r'$\Delta$''f [kHz]')
plt.ylabel('Max 'r'$\Gamma_{S21}$'' [dB]')
plt.legend(loc='lower left')
plt.savefig(f'{savepath_COMP}\\freq_maxGamma_JAN2024.png')
plt.close('all')


# freq Q0

plt.scatter(freq_C16_APR23, Q_0_C16_APR23, color=C16_APR23_colour, label="C16  APR23")
plt.scatter(freq_C20_APR23, Q_0_C20_APR23, color=C20_APR23_colour, label="C20  APR23")
plt.scatter(freq_JAN24, Q_0_JAN24, color=C8_JAN24_colour, label="C8    JAN24")
plt.xlabel(r'$\Delta$''f [kHz]')
plt.ylabel('$Q_0$')
plt.legend(loc='lower left')
plt.savefig(f'{savepath_COMP}\\freq_Q0_all.png')
plt.close('all')

plt.scatter(freq_JAN24, Q_0_JAN24, color=C8_JAN24_colour, label="C8    JAN24")
plt.xlabel(r'$\Delta$''f [kHz]')
plt.ylabel('$Q_0$')
plt.legend(loc='lower left')
plt.savefig(f'{savepath_COMP}\\freq_Q0_JAN2024.png')
plt.close('all')

""" Compare Distributions of freqs and max gammas """
# frequency distribution comparison

nbins = 4
# C16 APR23
mean_freq_C16_APR23 = np.mean(freq_C16_APR23)
std_freqC16_APR23 = np.std(freq_C16_APR23)
freq_C16_APR23_centred = [i-mean_freq_C16_APR23 for i in freq_C16_APR23]
n, bedges, patches = plt.hist(freq_C16_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
plt.text(min(freq_C16_APR23_centred), max(n), f'mean = {std_freqC16_APR23:1.3f} [kHz]\n'
                                              f's. dev. = {std_freqC16_APR23:1.3f} [kHz]\n')

plt.xlabel(r'frequency [kHz]')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\freq_C16_APR23_centred_hist.png')
plt.close('all')

# C20 APR23
mean_freq_C20_APR23 = np.mean(freq_C20_APR23)
std_freqC20_APR23 = np.std(freq_C20_APR23)
freq_C20_APR23_centred = [i-mean_freq_C20_APR23 for i in freq_C20_APR23]
n, bedges, patches = plt.hist(freq_C20_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C20_APR23_colour)
plt.text(min(freq_C20_APR23_centred), max(n), f'mean = {mean_freq_C20_APR23:1.3f} [kHz]\n'
                                              f's. dev. = {std_freqC20_APR23:1.3f} [kHz]\n')
plt.xlabel(r'frequency [kHz]')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\freq_C20_APR23_centred_hist.png')
plt.close('all')


C16_C20_freq_APR23 = freq_C16_APR23 + freq_C20_APR23
mean_freq_C16_C20_APR23 = np.mean(C16_C20_freq_APR23)
std_freqC16_C20_APR23 = np.std(C16_C20_freq_APR23)
C16_C20_freq_APR23_centred = [i-mean_freq_C16_C20_APR23 for i in C16_C20_freq_APR23]
n, bedges, patches = plt.hist(C16_C20_freq_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
plt.text(min(C16_C20_freq_APR23_centred), max(n),  f'mean = {mean_freq_C16_C20_APR23:1.3f} [kHz]\n'
                                              f's. dev. = {std_freqC16_C20_APR23:1.3f} [kHz]\n')
plt.xlabel(r'frequency [kHz]')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\freq_C16_C20_APR23_centred_hist.png')
plt.close('all')

# C8 JAN24
mean_freq_JAN24 = np.mean(freq_JAN24)
std_freq_JAN24 = np.std(freq_JAN24)
freq_JAN24_centred = [i-mean_freq_JAN24 for i in freq_JAN24]
n, bedges, patches = plt.hist(freq_JAN24_centred, bins = nbins, histtype='step', lw=0.8, color=C8_JAN24_colour)
plt.text(min(freq_JAN24_centred), max(n), f'mean = {mean_freq_JAN24:1.3f} [kHz]\n'
                                              f's. dev. = {std_freq_JAN24:1.3f} [kHz]\n')
plt.xlabel(r'frequency [kHz]')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\freq_JAN24_centred_hist.png')
plt.close('all')


Q_0_C16_APR23_mean = np.mean(Q_0_C16_APR23)
Q_0_C16_APR23_std = np.std(Q_0_C16_APR23)
beta_C16_APR23_mean = np.mean(beta_C16_APR23)
beta_C16_APR23_std = np.std(beta_C16_APR23)
Q_0_C20_APR23_mean = np.mean(Q_0_C20_APR23)
Q_0_C20_APR23_std = np.std(Q_0_C20_APR23)
beta_C20_APR23_mean = np.mean(beta_C20_APR23)
beta_C20_APR23_std = np.std(beta_C20_APR23)


# max_gamma distribution comparison

nbins = 4
# C16 APR23
mean_max_gamma_C16_APR23 = np.mean(max_gamma_C16_APR23)
std_C16_APR23 = np.std(max_gamma_C16_APR23)
max_gamma_C16_APR23_centred = [i-mean_max_gamma_C16_APR23 for i in max_gamma_C16_APR23]
n, bedges, patches = plt.hist(max_gamma_C16_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
plt.text(min(max_gamma_C16_APR23_centred), max(n), f's. dev. = {std_C16_APR23:1.3f} [dB]')
plt.xlabel(r'$\Delta$'r'Max 'r'$\Gamma_{S21}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\max_gamma_C16_APR23_centred_hist.png')
plt.close('all')

# C20 APR23
mean_max_gamma_C20_APR23 = np.mean(max_gamma_C20_APR23)
std_C20_APR23 = np.std(max_gamma_C20_APR23)
max_gamma_C20_APR23_centred = [i-mean_max_gamma_C20_APR23 for i in max_gamma_C20_APR23]
n, bedges, patches = plt.hist(max_gamma_C20_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C20_APR23_colour)
plt.text(min(max_gamma_C20_APR23_centred), max(n), f's. dev. = {std_C20_APR23:1.3f} [dB]')
plt.xlabel(r'$\Delta$'r'Max 'r'$\Gamma_{S21}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\max_gamma_C20_APR23_centred_hist.png')
plt.close('all')

# C16 C20 combined
C16_C20_max_gamma_APR23 = max_gamma_C16_APR23 + max_gamma_C20_APR23
mean_max_gamma_C16_C20_APR23 = np.mean(C16_C20_max_gamma_APR23)
std_max_gamma_C16_C20_APR23 = np.std(C16_C20_max_gamma_APR23)
C16_C20_max_gamma_APR23_centred = [i-mean_max_gamma_C16_C20_APR23 for i in C16_C20_max_gamma_APR23]
n, bedges, patches = plt.hist(C16_C20_max_gamma_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
plt.text(min(C16_C20_max_gamma_APR23_centred), max(n), f's. dev. = {std_max_gamma_C16_C20_APR23:1.3f} [dB]')
plt.xlabel(r'$\Delta$'r'Max 'r'$\Gamma_{S21}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\max_gamma_C16_C20_APR23_centred_hist.png')
plt.close('all')

# C8 JAN24
mean_max_gamma_JAN24 = np.mean(max_gamma_JAN24)
std_max_gamma_JAN24 = np.std(max_gamma_JAN24)
max_gamma_JAN24_centred = [i-mean_max_gamma_JAN24 for i in max_gamma_JAN24]
n, bedges, patches = plt.hist(max_gamma_JAN24_centred, bins = nbins, histtype='step', lw=0.8, color=C8_JAN24_colour)
plt.text(min(max_gamma_JAN24_centred), max(n), f's. dev. = {std_max_gamma_JAN24:1.3f} [dB]')
plt.xlabel(r'$\Delta$'r'Max 'r'$\Gamma_{S21}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\max_gamma_JAN24_centred_hist.png')
plt.close('all')



# All stats

Q_0_APR23 = Q_0_C16_APR23 + Q_0_C20_APR23
beta_APR23 = beta_C16_APR23 + beta_C20_APR23
Q_0_APR23_mean = np.mean(Q_0_APR23)
Q_0_APR23_std = np.std(Q_0_APR23)
beta_APR23_mean = np.mean(beta_APR23)
beta_APR23_std = np.std(beta_APR23)

Q_0_JAN24_mean = np.mean(Q_0_JAN24)
Q_0_JAN24_std = np.std(Q_0_JAN24)
beta_JAN24_mean = np.mean(beta_JAN24)
beta_JAN24_std = np.std(beta_JAN24)


print('\n2023')
print(f'Q_0 mean = {Q_0_APR23_mean:1.3f}\n'
      f'Q_0 SDev = {Q_0_APR23_std:1.3f}')
print(f'B mean = {beta_APR23_mean:1.3f}\n'
      f'B SDev = {beta_APR23_std:1.6f}')
print(f'max gamma mean = {mean_max_gamma_C16_C20_APR23:1.3f}\n'
      f'max gamma SDev = {std_max_gamma_C16_C20_APR23:1.6f}')
print(f'freq mean = {mean_freq_C16_C20_APR23:1.3f}\n'
      f'freq SDev = {std_freqC16_C20_APR23:1.3f}')

print('\n2024')
print(f'Q_0 mean = {Q_0_JAN24_mean:1.3f}\n'
      f'Q_0 SDev = {Q_0_JAN24_std:1.3f}')
print(f'B mean = {beta_JAN24_mean:1.3f}\n'
      f'B SDev = {beta_JAN24_std:1.6f}')
print(f'max gamma mean = {mean_max_gamma_JAN24:1.3f}\n'
      f'max gamma SDev = {std_max_gamma_JAN24:1.6f}')
print(f'freq mean = {mean_freq_JAN24:1.3f}\n'
      f'freq SDev = {std_freq_JAN24:1.3f}')


plt.scatter(Q_0_C16_APR23, max_gamma_C16_APR23, color=C16_APR23_colour, label="C16  APR23")
plt.scatter(Q_0_C20_APR23, max_gamma_C20_APR23, color=C20_APR23_colour, label="C20  APR23")
plt.scatter(Q_0_JAN24, max_gamma_JAN24, color=C8_JAN24_colour, label="C8    JAN24")
plt.xlabel(r'$Q_0$')
plt.ylabel('Max 'r'$\Gamma_{S21}$'' [dB]')
plt.legend(loc='lower right')
plt.savefig(f'{savepath_COMP}\\Q0_maxGamma_all.png')
plt.close('all')

# combined 23 24 histograms



# C8 JAN24
mean_max_gamma_JAN24 = np.mean(max_gamma_JAN24)
std_JAN24 = np.std(max_gamma_JAN24)
max_gamma_JAN24_centred = [i-mean_max_gamma_JAN24 for i in max_gamma_JAN24]
n, bedges, patches = plt.hist(max_gamma_JAN24_centred, bins = nbins, histtype='step', lw=0.8, color=C8_JAN24_colour)
n, bedges, patches = plt.hist(C16_C20_max_gamma_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
plt.text(min(max_gamma_JAN24_centred), max(n), f's. dev. APR23 = {std_max_gamma_C16_C20_APR23:1.6f} [dB]\n'
                                               f's. dev. JAN24 = {std_max_gamma_JAN24:1.6f} [dB]')
plt.xlabel(r'$\Delta$'r'Max 'r'$\Gamma_{S21}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\max_gamma_APR23_JAN24_centred_hist.png')
plt.close('all')



mean_freq_JAN24 = np.mean(freq_JAN24)
std_JAN24 = np.std(freq_JAN24)
freq_JAN24_centred = [i-mean_freq_JAN24 for i in freq_JAN24]
n, bedges, patches = plt.hist(freq_JAN24_centred, bins = nbins, histtype='step', lw=0.8, color=C8_JAN24_colour)
n, bedges, patches = plt.hist(C16_C20_freq_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
plt.text(75., 4., f's. dev. APR23 = {std_freqC16_C20_APR23:1.3f} [kHz]\n'
                                          f's. dev. JAN24 = {std_JAN24:1.3f} [kHz]')
plt.xlabel(r'frequency [kHz]')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\freq_APR23_JAN24_centred_hist.png')
plt.close('all')


mean_max_gamma_JAN24 = np.mean(max_gamma_JAN24)
std_JAN24 = np.std(max_gamma_JAN24)
max_gamma_JAN24_centred = [i-mean_max_gamma_JAN24 for i in max_gamma_JAN24]
n, bedges, patches = plt.hist(max_gamma_JAN24_centred, bins = nbins, histtype='step', lw=0.8, color=C8_JAN24_colour)
n, bedges, patches = plt.hist(C16_C20_max_gamma_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
plt.text(min(max_gamma_JAN24_centred), max(n), f's. dev. APR23 = {std_max_gamma_C16_C20_APR23:1.6f} [dB]\n'
                                               f's. dev. JAN24 = {std_max_gamma_JAN24:1.6f} [dB]')
plt.xlabel(r'$\Delta$'r'Max 'r'$\Gamma_{S21}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\max_gamma_APR23_JAN24_centred_hist.png')
plt.close('all')



freq_JAN24_centred = [i-mean_freq_JAN24 for i in freq_JAN24]
n, bedges, patches = plt.hist(freq_JAN24_centred, bins = nbins, histtype='step', lw=0.8, color=C8_JAN24_colour)
n, bedges, patches = plt.hist(freq_C16_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
n, bedges, patches = plt.hist(freq_C20_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C20_APR23_colour)

plt.text(75., 3., f's. dev. C16 APR23 = {std_freqC16_APR23:1.3f} [kHz]\n'
                  f's. dev. C20 APR23 = {std_freqC20_APR23:1.3f} [kHz]\n'
                  f's. dev. C8   JAN24 = {std_JAN24:1.3f} [kHz]\n'
                  )
plt.xlabel(r'$\Delta$''f [kHz]')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\freq_C8_C16_C20_centred_hist.png')
plt.close('all')


max_gamma_JAN24_centred = [i-mean_max_gamma_JAN24 for i in max_gamma_JAN24]
n, bedges, patches = plt.hist(max_gamma_JAN24_centred, bins = nbins, histtype='step', lw=0.8, color=C8_JAN24_colour)
n, bedges, patches = plt.hist(max_gamma_C16_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
n, bedges, patches = plt.hist(max_gamma_C20_APR23_centred, bins = nbins, histtype='step', lw=0.8, color=C20_APR23_colour)

plt.text(-0.1, 3., f's. dev. C16 APR23 = {std_C16_APR23:1.3f} [dB]\n'
                  f's. dev. C20 APR23 = {std_C20_APR23:1.3f} [dB]\n'
                  f's. dev. C8   JAN24 = {std_JAN24:1.3f} [dB]\n'
                  )
plt.xlabel(r'$\Delta$'r'Max 'r'$\Gamma_{S21}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\max_gamma_C8_C16_C20_centred_hist.png')
plt.close('all')


# noise on power

temp_operational_C = 50.
temp_operational_K = 273.15 + temp_operational_C
temp_experiment_C = 20.8
temp_experiment_K = 273.15 + temp_experiment_C

# C16_APR23

power_noise_C16_APR23 = [hmm.get_noise_power(temp_experiment_K, bandwidth_C16_APR23[i]) for i in range(len(bandwidth_C16_APR23))]
P_out_C16_APR23 = [hmm.get_P_out_from_dB_and_P_in(max_gamma_C16_APR23[i], 10.e-2) for i in range(len(max_gamma_C16_APR23))]

print(f'\nC16_APR23\nBandwidths [Hz]          P_out [W]                    P_noise [W]')
for i in range(len(bandwidth_C16_APR23)):
    print(f'{bandwidth_C16_APR23[i]}        {P_out_C16_APR23[i]}       {power_noise_C16_APR23[i]}')

n, bedges, patches = plt.hist(power_noise_C16_APR23, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)

plt.xlabel(r'$P_{noise}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\P_noise_hist.png')
plt.close('all')

n, bedges, patches = plt.hist(P_out_C16_APR23, bins = nbins, histtype='step', lw=0.8, color=C16_APR23_colour)
plt.xlabel(r'$P_{out}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\P_out_hist.png')
plt.close('all')


# C8_JAN24

bandwidth_C8_JAN24 = bandwidth_C8_JAN24[:-1]  #  remove cathode out outlier

power_noise_C8_JAN24 = [hmm.get_noise_power(temp_operational_K, bandwidth_C8_JAN24[i]) for i in range(len(bandwidth_C8_JAN24))]
P_out_C8_JAN24 = [hmm.get_P_out_from_dB_and_P_in(max_gamma_JAN24[i], 10.e-2) for i in range(len(max_gamma_JAN24))]

print(f'\nC8_JAN24\nBandwidths [Hz]          P_out [W]                    P_noise [W]')
for i in range(len(bandwidth_C8_JAN24)):
    print(f'{bandwidth_C8_JAN24[i]}        {P_out_C8_JAN24[i]}       {power_noise_C8_JAN24[i]}')

n, bedges, patches = plt.hist(power_noise_C8_JAN24, bins = nbins, histtype='step', lw=0.8, color=C8_JAN24_colour)

plt.xlabel(r'$P_{noise}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\P_noise_hist.png')
plt.close('all')

n, bedges, patches = plt.hist(P_out_C8_JAN24, bins = nbins, histtype='step', lw=0.8, color=C8_JAN24_colour)
plt.xlabel(r'$P_{out}$')
plt.ylabel('N')
plt.savefig(f'{savepath_COMP}\\P_out_C8_JAN24_hist.png')
plt.close('all')