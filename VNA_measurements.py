import numpy as np
import argparse
import sys, os, csv
from matplotlib import pyplot as plt
import pandas as pd
from snps import SNPs
import PhD_Master_Module as pmm

savepath = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20230418\analysis'
snp_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20230418\data'
fnames = os.listdir(snp_addr)
fname = fnames[0]
fname_and_addr = f'{snp_addr}\\{fname}'
print(f'{fnames = }')



C16_with_stop_fnames = ['C16_20.7degc_insert_1.s2p', 'C16_20.7degc_insert_2.s2p', 'C16_20.7degc_insert_3.s2p']
C16_without_stop_fnames = ['C16_20.8degc_insert_1_hard_stop_totally_removed.s2p', 'C16_20.8degc_insert_2hard_stop_totally_removed.s2p']
C20_with_stop_fnames = ['C20_20.7degc_insert_1.s2p', 'C20_20.7degc_insert_2.s2p', 'C20_20.7degc_insert_3.s2p',]
C20_without_stop_fnames = ['C20_20.7degc_insert_4_hard_stop_totally_removed.s2p',
                           'C20_20.7degc_insert_5_hard_stop_totally_removed.s2p',
                           'C20_20.7degc_insert_1_hard_stop_removed.s2p',
                           'C20_20.7degc_insert_2_hard_stop_removed.s2p',
                           'C20_20.7degc_insert_3_hard_stop_removed.s2p',
                           ]

print(f'{C16_with_stop_fnames = }')
print(f'{C16_without_stop_fnames = }')
print(f'{C20_with_stop_fnames = }')
print(f'{C20_without_stop_fnames = }')

# from CONVERTF CERN app at 20.8 deg C
design_res_freq_MHz = 2999.982788

def get_response_freq(freq_list, S11, print_plot=False):
    min_gamma, min_idx = pmm.get_min_val_idx_from_list(S11)
    res_freq = freq_list[min_idx]

    if print_plot:
        plt.scatter(freq_list, s11, marker='.', s=3, color='k')
        plt.scatter(res_freq, min_gamma, marker='x', s=60, color='r')
        plt.show()

    return res_freq, min_gamma

def get_Q_L(freq_list, s21, savename, print_plot=False):
    '''
    for a given s21 trace it will find the peak resonant frequency
    and find the -3dB bandwidth
    The loaded Q is given by: res_freq/-3dB bandwidth
    :param s21:
    :return:
    '''
    # plt.scatter(freq_list, s21)
    # plt.show()
    # db = 10.*np.log10(y/max_gamma)
    # m3db = 10.**(-3./10.)*max_gamma
    clip_idx = 100
    freq_list = freq_list[clip_idx:-clip_idx]
    s21_smooth = pmm.smooth(s21, 100)[clip_idx:-clip_idx]
    s21 = s21[clip_idx:-clip_idx]
    max_gamma, max_idx = pmm.get_max_val_idx_from_list(s21_smooth)
    max_gamma_raw, max_idx_raw = pmm.get_max_val_idx_from_list(s21)
    res_freq = freq_list[max_idx]
    minus3dB_gamma = max_gamma -3.

    min_idx_above_minus3dB = [i for i in range(len(s21_smooth)) if s21_smooth[i] >= minus3dB_gamma][0]
    max_idx_above_minus3dB = [i for i in range(len(s21_smooth)) if s21_smooth[i] >= minus3dB_gamma][-1]
    min_feq_minus3dB = freq_list[min_idx_above_minus3dB]
    max_feq_minus3dB = freq_list[max_idx_above_minus3dB]
    bandwidth = max_feq_minus3dB - min_feq_minus3dB

    Q_L = res_freq/bandwidth

    if print_plot:
        print(f'{max_gamma = }\n{minus3dB_gamma = }\n{res_freq = }')
        print(f'{min_idx_above_minus3dB = }\n{max_idx_above_minus3dB = }\n{bandwidth = }')

        plt.scatter(freq_list, s21, marker='.', s=0.5, color='k')
        plt.plot(freq_list, s21_smooth, ls='-', lw=1., color='darkorange')
        plt.scatter(res_freq, max_gamma, marker='x', s=70, color='r', zorder=4)
        plt.hlines(minus3dB_gamma, freq_list[min_idx_above_minus3dB], freq_list[max_idx_above_minus3dB], ls='--', lw=0.8, color='r')
        plt.savefig(f'{savepath}\\{savename}.png')
        plt.close('all')
        print(f'{Q_L = }')

    return Q_L, s21_smooth, freq_list, res_freq, max_gamma, max_gamma_raw, min_feq_minus3dB, max_feq_minus3dB, bandwidth

def get_beta(s11, overcoupled=False, print_plot=False):
    min_gamma, min_idx = pmm.get_min_val_idx_from_list(s11)
    max_gamma, max_idx = pmm.get_max_val_idx_from_list(s11)

    if overcoupled:
        beta = (max_gamma + min_gamma) / (max_gamma - min_gamma)
    else:
        beta = (max_gamma - min_gamma) / (max_gamma + min_gamma)

    if beta < 0.:
        beta = -beta
    else:
        pass

    if print_plot:
        print(f'{min_gamma = }')
        print(f'{max_gamma = }')
        print(f'{beta = }')

    return beta

def get_Q_0(Q_L, beta, print_plot=False):
    Q_0 = (1. + beta)*Q_L

    if print_plot:
        print(f'{Q_0 = }')

    return Q_0

def get_zcoord_from_freq_delta(freq_delta):
    gradient = -1437. # Hz/um
    z = freq_delta*(1./gradient) # um

    return z

def is_overcoupled(s11, s11_ph, name, print_plot=False):
    s11_ph_unwrapped = s11_ph
    # s11_ph_unwrapped = [i if i>0. else -i for i in s11_ph]
    s11_ph_radians = [i/(180.*np.pi) for i in s11_ph_unwrapped]
    re_s11 = [s11[i]*np.cos(s11_ph_radians[i]) for i in range(len(s11))]
    im_s11 = [s11[i]*np.sin(s11_ph_radians[i]) for i in range(len(s11))]
    quadrant_1 = [i for i in range(len(s11)) if re_s11[i] > 0. and im_s11[i] > 0.]
    quadrant_2 = [i for i in range(len(s11)) if re_s11[i] > 0. and im_s11[i] < 0.]
    quadrant_3 = [i for i in range(len(s11)) if re_s11[i] < 0. and im_s11[i] > 0.]
    quadrant_4 = [i for i in range(len(s11)) if re_s11[i] < 0. and im_s11[i] < 0.]

    if print_plot:
        print(f'{len(quadrant_1) = }')
        print(f'{len(quadrant_2) = }')
        print(f'{len(quadrant_3) = }')
        print(f'{len(quadrant_4) = }')

        plt.scatter(re_s11, im_s11, marker='o', s=5, color='k')
        plt.show()

def read_s2p(fname_and_addr):
    freq = []
    S11_a = []
    S11_b = []
    S21_a = []
    S21_b = []
    S12_a = []
    S12_b = []
    S22_a = []
    S22_b = []
    with open(fname_and_addr, 'r') as f:
        for idx, row in enumerate(f):
            if '!' in row or '#' in row:
                pass
            else:
                row = row.split(' ')
                # if idx == 10:
                #     print(type(row))
                #     print(row[1])
                #     print(float(row[1]))

                #if float(row[0]) > 2.99:
                freq.append(float(row[0]))
                S11_a.append(float(row[1]))
                S11_b.append(float(row[2]))
                S21_a.append(float(row[3]))
                S21_b.append(float(row[4]))
                S12_a.append(float(row[5]))
                S12_b.append(float(row[6]))
                S22_a.append(float(row[7]))
                S22_b.append(float(row[8]))
                # else:
                #     pass

    return freq, S11_a, S11_b, S21_a, S21_b, S12_a, S12_b, S22_a, S22_b



data_dict = {}
name_keys = []
for f in fnames:

    if 'remove' in f:
        name = f'{f[:3]}_{f[20]}_removed'
        data_dict[name] = {}
        data_dict[name]['stop_present?'] = False
    else:
        name = f'{f[:3]}_{f[20]}'
        data_dict[name] = {}
        data_dict[name]['stop_present?'] = True

    data_dict[name]['fname'] = f
    name_keys.append(name)

    freq, s11, s11_ph, s21, s21_ph, s12, s12_ph, s22, s22_ph = read_s2p(f'{snp_addr}\\{data_dict[name]["fname"]}')
    data_dict[name]['freq_list'] = freq
    data_dict[name]['s11'] = s11
    data_dict[name]['s21'] = s21
    res_freq, min_gamma_s11 = get_response_freq(freq, s11)
    data_dict[name]['res_freq'] = res_freq
    data_dict[name]['min_gamma_s11'] = min_gamma_s11
    Q_L, s21_smooth, freq_list_s21_smooth, res_freq_s21, max_gamma_s21, max_gamma_raw_s21, min_feq_minus3dB_s21, max_feq_minus3dB_s21, bandwidth = get_Q_L(freq, s21, name, print_plot=False)
    data_dict[name]['Q_L'] = Q_L
    data_dict[name]['s21_smooth'] = s21_smooth
    data_dict[name]['freq_list_s21_smooth'] = freq_list_s21_smooth
    data_dict[name]['res_freq_s21'] = res_freq_s21
    data_dict[name]['max_gamma_s21'] = max_gamma_s21
    data_dict[name]['max_gamma_raw_s21'] = max_gamma_raw_s21
    data_dict[name]['min_feq_minus3dB_s21'] = min_feq_minus3dB_s21
    data_dict[name]['max_feq_minus3dB_s21'] = max_feq_minus3dB_s21
    data_dict[name]['bandwidth'] = bandwidth

    data_dict[name]['overcoupled'] = is_overcoupled(s11, s11_ph, name, print_plot=False)
    beta = get_beta(s11, overcoupled=True)
    data_dict[name]['beta'] = beta
    Q_0 = get_Q_0(Q_L, beta)
    data_dict[name]['Q_0'] = Q_0
    design_freq_Hz = design_res_freq_MHz*1e6
    delta_freq = res_freq-design_freq_Hz
    z_coord = get_zcoord_from_freq_delta(delta_freq)
    data_dict[name]['z_coord'] = z_coord

    print(f'\n{name = }')
    print(f'{res_freq = }')
    print(f'{Q_L = }')
    print(f'{Q_0 = }')
    print(f'{beta = }')

'''Plots'''

# s21
for name in name_keys:
    s11= data_dict[name]['s11']
    s21_smooth = data_dict[name]['s21_smooth']
    freq_list_s21_smooth = data_dict[name]['freq_list_s21_smooth']
    freq_list = data_dict[name]['freq_list']
    freq_list_GHz = [i/1e9 for i in freq_list]
    s21 = data_dict[name]['s21']
    res_freq = data_dict[name]['res_freq']
    res_freq_GHz = res_freq/1e9
    res_freq_s21 = data_dict[name]['res_freq_s21']
    max_gamma_s21 = data_dict[name]['max_gamma_s21']
    max_gamma_raw_s21 = data_dict[name]['max_gamma_raw_s21']
    min_feq_minus3dB_s21 = data_dict[name]['min_feq_minus3dB_s21']
    max_feq_minus3dB_s21 = data_dict[name]['max_feq_minus3dB_s21']
    bandwidth = data_dict[name]['bandwidth']
    Q_0 = data_dict[name]['Q_0']
    Q_L = data_dict[name]['Q_L']
    beta = data_dict[name]['beta']
    z_coord = data_dict[name]['z_coord']

    plt.scatter(freq_list, s21, marker='.', s=0.5, color='k')
    plt.plot(freq_list_s21_smooth, s21_smooth, ls='-', lw=1., color='darkorange')
    plt.scatter(res_freq, max_gamma_s21, marker='x', s=70, color='r', zorder=4)
    plt.hlines(max_gamma_s21-3., max_feq_minus3dB_s21, min_feq_minus3dB_s21, ls='--', lw=0.8, color='r')
    plt.text(3.002e9, -85., f'f = {res_freq_GHz:1.6f} GHz\n'r'$\beta$'f' = {beta:1.4f}\n'r'$Q_L$'f' = {Q_L:1.0f}\n'r'$Q_0$'f' = {Q_0:1.0f}\n')
    plt.savefig(f'{savepath}\\{name}_s21.png')
    plt.close('all')

# beta
for name in name_keys:
    s11= data_dict[name]['s11']
    s21_smooth = data_dict[name]['s21_smooth']
    freq_list_s21_smooth = data_dict[name]['freq_list_s21_smooth']
    freq_list = data_dict[name]['freq_list']
    freq_list_GHz = [i/1e9 for i in freq_list]
    s21 = data_dict[name]['s21']
    res_freq = data_dict[name]['res_freq']
    res_freq_GHz = res_freq/1e9
    res_freq_s21 = data_dict[name]['res_freq_s21']
    max_gamma_s21 = data_dict[name]['max_gamma_s21']
    max_gamma_raw_s21 = data_dict[name]['max_gamma_raw_s21']
    min_feq_minus3dB_s21 = data_dict[name]['min_feq_minus3dB_s21']
    max_feq_minus3dB_s21 = data_dict[name]['max_feq_minus3dB_s21']
    bandwidth = data_dict[name]['bandwidth']
    Q_0 = data_dict[name]['Q_0']
    Q_L = data_dict[name]['Q_L']
    beta = data_dict[name]['beta']
    z_coord = data_dict[name]['z_coord']

    plt.scatter(z_coord, beta, marker='o', s=10, color='k')
    # plt.plot(freq_list_s21_smooth, s21_smooth, ls='-', lw=1., color='darkorange')
    # plt.scatter(res_freq, max_gamma_s21, marker='x', s=70, color='r', zorder=4)
    # plt.hlines(max_gamma_s21-3., max_feq_minus3dB_s21, min_feq_minus3dB_s21, ls='--', lw=0.8, color='r')
    # plt.text(3.002e9, -85., f'f = {res_freq_GHz:1.6f} GHz\n'r'$\beta$'f' = {beta:1.4f}\n'r'$Q_L$'f' = {Q_L:1.0f}\n'r'$Q_0$'f' = {Q_0:1.0f}\n')
plt.savefig(f'{savepath}\\zcoord_beta.png')
plt.close('all')

# Q_0
for name in name_keys:
    s11= data_dict[name]['s11']
    s21_smooth = data_dict[name]['s21_smooth']
    freq_list_s21_smooth = data_dict[name]['freq_list_s21_smooth']
    freq_list = data_dict[name]['freq_list']
    freq_list_GHz = [i/1e9 for i in freq_list]
    s21 = data_dict[name]['s21']
    res_freq = data_dict[name]['res_freq']
    res_freq_GHz = res_freq/1e9
    res_freq_s21 = data_dict[name]['res_freq_s21']
    max_gamma_s21 = data_dict[name]['max_gamma_s21']
    max_gamma_raw_s21 = data_dict[name]['max_gamma_raw_s21']
    min_feq_minus3dB_s21 = data_dict[name]['min_feq_minus3dB_s21']
    max_feq_minus3dB_s21 = data_dict[name]['max_feq_minus3dB_s21']
    bandwidth = data_dict[name]['bandwidth']
    Q_0 = data_dict[name]['Q_0']
    Q_L = data_dict[name]['Q_L']
    beta = data_dict[name]['beta']
    z_coord = data_dict[name]['z_coord']

    plt.scatter(z_coord, Q_0, marker='o', s=10, color='k')
    # plt.plot(freq_list_s21_smooth, s21_smooth, ls='-', lw=1., color='darkorange')
    # plt.scatter(res_freq, max_gamma_s21, marker='x', s=70, color='r', zorder=4)
    # plt.hlines(max_gamma_s21-3., max_feq_minus3dB_s21, min_feq_minus3dB_s21, ls='--', lw=0.8, color='r')
    # plt.text(3.002e9, -85., f'f = {res_freq_GHz:1.6f} GHz\n'r'$\beta$'f' = {beta:1.4f}\n'r'$Q_L$'f' = {Q_L:1.0f}\n'r'$Q_0$'f' = {Q_0:1.0f}\n')
plt.savefig(f'{savepath}\\zcoord_Q_0.png')
plt.close('all')

# Q_L
for name in name_keys:
    s11= data_dict[name]['s11']
    s21_smooth = data_dict[name]['s21_smooth']
    freq_list_s21_smooth = data_dict[name]['freq_list_s21_smooth']
    freq_list = data_dict[name]['freq_list']
    freq_list_GHz = [i/1e9 for i in freq_list]
    s21 = data_dict[name]['s21']
    res_freq = data_dict[name]['res_freq']
    res_freq_GHz = res_freq/1e9
    res_freq_s21 = data_dict[name]['res_freq_s21']
    max_gamma_s21 = data_dict[name]['max_gamma_s21']
    max_gamma_raw_s21 = data_dict[name]['max_gamma_raw_s21']
    min_feq_minus3dB_s21 = data_dict[name]['min_feq_minus3dB_s21']
    max_feq_minus3dB_s21 = data_dict[name]['max_feq_minus3dB_s21']
    bandwidth = data_dict[name]['bandwidth']
    Q_0 = data_dict[name]['Q_0']
    Q_L = data_dict[name]['Q_L']
    beta = data_dict[name]['beta']
    z_coord = data_dict[name]['z_coord']

    plt.scatter(z_coord, Q_L, marker='o', s=10, color='k')
    # plt.plot(freq_list_s21_smooth, s21_smooth, ls='-', lw=1., color='darkorange')
    # plt.scatter(res_freq, max_gamma_s21, marker='x', s=70, color='r', zorder=4)
    # plt.hlines(max_gamma_s21-3., max_feq_minus3dB_s21, min_feq_minus3dB_s21, ls='--', lw=0.8, color='r')
    # plt.text(3.002e9, -85., f'f = {res_freq_GHz:1.6f} GHz\n'r'$\beta$'f' = {beta:1.4f}\n'r'$Q_L$'f' = {Q_L:1.0f}\n'r'$Q_0$'f' = {Q_0:1.0f}\n')
plt.savefig(f'{savepath}\\zcoord_Q_L.png')
plt.close('all')

# max_gamma_s21
for name in name_keys:
    s11= data_dict[name]['s11']
    s21_smooth = data_dict[name]['s21_smooth']
    freq_list_s21_smooth = data_dict[name]['freq_list_s21_smooth']
    freq_list = data_dict[name]['freq_list']
    freq_list_GHz = [i/1e9 for i in freq_list]
    s21 = data_dict[name]['s21']
    res_freq = data_dict[name]['res_freq']
    res_freq_GHz = res_freq/1e9
    res_freq_s21 = data_dict[name]['res_freq_s21']
    max_gamma_s21 = data_dict[name]['max_gamma_s21']
    max_gamma_raw_s21 = data_dict[name]['max_gamma_raw_s21']
    min_feq_minus3dB_s21 = data_dict[name]['min_feq_minus3dB_s21']
    max_feq_minus3dB_s21 = data_dict[name]['max_feq_minus3dB_s21']
    bandwidth = data_dict[name]['bandwidth']
    Q_0 = data_dict[name]['Q_0']
    Q_L = data_dict[name]['Q_L']
    beta = data_dict[name]['beta']
    z_coord = data_dict[name]['z_coord']

    plt.scatter(z_coord, max_gamma_s21, marker='o', s=10, color='k')
    # plt.scatter(z_coord, max_gamma_raw_s21, marker='o', s=10, color='r')
    # plt.plot(freq_list_s21_smooth, s21_smooth, ls='-', lw=1., color='darkorange')
    # plt.scatter(res_freq, max_gamma_s21, marker='x', s=70, color='r', zorder=4)
    # plt.hlines(max_gamma_s21-3., max_feq_minus3dB_s21, min_feq_minus3dB_s21, ls='--', lw=0.8, color='r')
    # plt.text(3.002e9, -85., f'f = {res_freq_GHz:1.6f} GHz\n'r'$\beta$'f' = {beta:1.4f}\n'r'$Q_L$'f' = {Q_L:1.0f}\n'r'$Q_0$'f' = {Q_0:1.0f}\n')
plt.savefig(f'{savepath}\\zcoord_max_gamma_s21.png')
plt.close('all')

# max_gamma_raw_s21
for name in name_keys:
    s11= data_dict[name]['s11']
    s21_smooth = data_dict[name]['s21_smooth']
    freq_list_s21_smooth = data_dict[name]['freq_list_s21_smooth']
    freq_list = data_dict[name]['freq_list']
    freq_list_GHz = [i/1e9 for i in freq_list]
    s21 = data_dict[name]['s21']
    res_freq = data_dict[name]['res_freq']
    res_freq_GHz = res_freq/1e9
    res_freq_s21 = data_dict[name]['res_freq_s21']
    max_gamma_s21 = data_dict[name]['max_gamma_s21']
    max_gamma_raw_s21 = data_dict[name]['max_gamma_raw_s21']
    min_feq_minus3dB_s21 = data_dict[name]['min_feq_minus3dB_s21']
    max_feq_minus3dB_s21 = data_dict[name]['max_feq_minus3dB_s21']
    bandwidth = data_dict[name]['bandwidth']
    Q_0 = data_dict[name]['Q_0']
    Q_L = data_dict[name]['Q_L']
    beta = data_dict[name]['beta']
    z_coord = data_dict[name]['z_coord']

    # plt.scatter(z_coord, max_gamma_s21, marker='o', s=10, color='k')
    plt.scatter(z_coord, max_gamma_raw_s21, marker='o', s=10, color='r')
    # plt.plot(freq_list_s21_smooth, s21_smooth, ls='-', lw=1., color='darkorange')
    # plt.scatter(res_freq, max_gamma_s21, marker='x', s=70, color='r', zorder=4)
    # plt.hlines(max_gamma_s21-3., max_feq_minus3dB_s21, min_feq_minus3dB_s21, ls='--', lw=0.8, color='r')
    # plt.text(3.002e9, -85., f'f = {res_freq_GHz:1.6f} GHz\n'r'$\beta$'f' = {beta:1.4f}\n'r'$Q_L$'f' = {Q_L:1.0f}\n'r'$Q_0$'f' = {Q_0:1.0f}\n')
plt.savefig(f'{savepath}\\zcoord_max_gamma_raw_s21.png')
plt.close('all')