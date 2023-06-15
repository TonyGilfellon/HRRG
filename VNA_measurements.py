import numpy as np
import argparse
import sys, os, csv, pickle
from matplotlib import pyplot as plt
import pandas as pd
from snps import SNPs
import PhD_Master_Module as pmm
import skrf as rf

# correlate conditioning data with experimental

# I need the same cavity forward power readings
#   1. straight line fit for both C20 & C16
#   2. delta intesect, c = abs(C20_c - C16_c)
#   3. for C20: 2.37 (C16 probe power at a cavity for pow of 4.45) - delta_c = answer needed

# read in raw data or load pickle?
read_data = False
# exit()

savepath = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20230418\analysis'
snp_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20230418\data'
pkl_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20230418'
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
design_freq_Hz = 2999982788.



convert_freq_Hz = design_freq_Hz - 2998500000.
convert_freq_GHz = convert_freq_Hz/1e9
print(f'{convert_freq_Hz = }')
print(f'{convert_freq_GHz = }')
# def get_gamma(s11_list):
#     return [(2.*50.*i)/(1.-i) for i in s11_list]

def get_FF_from_zcoord(zcoord):
    A_poly_hrrg = 8.725565888064103e-06
    B_poly_hrrg = 0.015100094187422454
    C_poly_hrrg = 96.58430786759465
    gradient = 0.051
    FF = A_poly_hrrg*zcoord**2. + B_poly_hrrg*zcoord + C_poly_hrrg

    return FF

def get_response_freq(freq_list, S11, print_plot=False):
    min_gamma, min_idx = pmm.get_min_val_idx_from_list(S11)
    res_freq = freq_list[min_idx]

    if print_plot:
        plt.scatter(freq_list, s11, marker='.', s=3, color='k')
        plt.scatter(res_freq, min_gamma, marker='x', s=60, color='r')
        plt.show()

    return res_freq, min_gamma

def get_s21_at_res_freq(freq_list, s21, res_freq):
    idx = pmm.find_nearest_value_index_in_list(freq_list, res_freq)
    s21_at_res = s21[idx] 
    
    return s21_at_res

def get_crude_noise(freq_list, s21):
    noise = [s21[i+1] - s21[i] for i in range(len(s21)-1)]
    freq_list_noise = freq_list[:-1]

    return freq_list_noise, noise


def get_s21_data(freq_list, s21, delta_Db, savename, print_plot=False):
    '''
    for a given s21 trace it will find the peak resonant frequency
    and find the -3dB bandwidth
    The loaded Q is given by: res_freq/-3dB bandwidth
    :param s21:
    :return:
    '''
    delta_Db = float(delta_Db)
    # plt.scatter(freq_list, s21)
    # plt.show()
    # db = 10.*np.log10(y/max_gamma)
    # m3db = 10.**(-3./10.)*max_gamma
    clip_idx = 100
    # freq_list_for_smoothing = freq_list[clip_idx:-clip_idx]
    freq_list_smooth, s21_smooth = pmm.smooth(freq_list, s21, 100, clip_idx)
    s21 = s21
    max_gamma, res_freq_idx = pmm.get_max_val_idx_from_list(s21_smooth)
    max_gamma_raw, max_idx_raw = pmm.get_max_val_idx_from_list(s21)
    res_freq = freq_list_smooth[res_freq_idx]

    res_freq_idx_smooth = res_freq_idx #+ clip_idx



    minus3dB_gamma = max_gamma -3.
    min_idx_above_minus3dB = [i for i in range(len(s21_smooth)) if s21_smooth[i] >= minus3dB_gamma][0]
    max_idx_above_minus3dB = [i for i in range(len(s21_smooth)) if s21_smooth[i] >= minus3dB_gamma][-1]
    min_feq_minus3dB = freq_list_smooth[min_idx_above_minus3dB]
    max_feq_minus3dB = freq_list_smooth[max_idx_above_minus3dB]
    bandwidth = max_feq_minus3dB - min_feq_minus3dB

    minus_delta_dB_gamma = max_gamma - delta_Db
    min_idx_above_minus_delta_dB = [i for i in range(len(s21_smooth)) if s21_smooth[i] >= minus_delta_dB_gamma][0]
    max_idx_above_minus_delta_dB = [i for i in range(len(s21_smooth)) if s21_smooth[i] >= minus_delta_dB_gamma][-1]
    min_feq_minus_delta_dB = freq_list_smooth[min_idx_above_minus_delta_dB]
    max_feq_minus_delta_dB = freq_list_smooth[max_idx_above_minus_delta_dB]

    freq_list_noise, noise = get_crude_noise(freq_list, s21)
    freq_list_noise_centred = [i-res_freq for i in freq_list_noise]
    freq_list_noise_centred_zoom = freq_list_noise_centred[min_idx_above_minus_delta_dB+clip_idx:max_idx_above_minus_delta_dB+clip_idx]
    noise_zoom = noise[min_idx_above_minus_delta_dB:max_idx_above_minus_delta_dB]

    plt.scatter(freq_list_noise_centred, noise, marker='o', s=1, color='k')
    plt.xlabel(r'$\Delta f$')
    plt.ylabel(r'$\Delta Db$')
    plt.savefig(f'{savepath}\\{savename}_noise.png')
    plt.close('all')

    plt.scatter(freq_list_noise_centred_zoom, noise_zoom, marker='o', s=4, color='k')
    plt.xlabel(r'$\Delta f$')
    plt.ylabel(r'$\Delta Db$')
    plt.savefig(f'{savepath}\\{savename}_noise_zoom.png')
    plt.close('all')

    plt.hist(noise, bins=int(np.sqrt(len(noise))), histtype='step', lw=0.8, color='k')
    plt.ylabel('N')
    plt.xlabel(r'$\Delta Db$')
    plt.savefig(f'{savepath}\\{savename}_noise_hist.png')
    plt.close('all')

    plt.hist(noise_zoom, bins=int(np.sqrt(len(noise_zoom))), histtype='step', lw=0.8, color='k')
    plt.ylabel('N')
    plt.xlabel(r'$\Delta Db$')
    plt.savefig(f'{savepath}\\{savename}_noise_zoom_hist.png')
    plt.close('all')

    std_noise_zoom = np.std(noise_zoom)
    mean_noise_zoom = np.mean(noise_zoom)



    half_delta_f_delta_Db = (max_feq_minus_delta_dB - min_feq_minus_delta_dB)/2.

    print(f'{delta_Db = }')
    print(f'{min_feq_minus_delta_dB = }')
    print(f'{max_feq_minus_delta_dB = }')
    print(f'{half_delta_f_delta_Db = }')
    print(f'{min_feq_minus3dB = }')
    print(f'{max_feq_minus3dB = }')
    half_bandwidth = bandwidth/2.
    print(f'{half_bandwidth = }')
    print(f'{len(freq_list_smooth) = }')
    print(f'{len(s21_smooth) = }')


    Q_L = res_freq/bandwidth

    if print_plot:
        # print(f'{max_gamma = }\n{minus3dB_gamma = }\n{res_freq = }')
        # print(f'{min_idx_above_minus3dB = }\n{max_idx_above_minus3dB = }\n{bandwidth = }')

        plt.scatter(freq_list, s21, marker='.', s=0.5, color='k')
        plt.plot(freq_list_smooth, s21_smooth, ls='-', lw=1., color='darkorange')
        plt.scatter(res_freq, max_gamma, marker='x', s=70, color='r', zorder=4)
        plt.hlines(minus3dB_gamma, min_feq_minus3dB, max_feq_minus3dB, ls='--', lw=0.8, color='r')
        plt.hlines(minus_delta_dB_gamma, min_feq_minus_delta_dB, max_feq_minus_delta_dB, ls='--', lw=0.8, color='r')
        plt.savefig(f'{savepath}\\{savename}.png')
        plt.close('all')
        #print(f'{Q_L = }')

    return Q_L, s21_smooth, freq_list_smooth, freq_list , res_freq, max_gamma, max_gamma_raw, min_feq_minus3dB, max_feq_minus3dB, bandwidth, half_delta_f_delta_Db, res_freq_idx_smooth, std_noise_zoom

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
  
    #intercept = -5.0000
 
    z = freq_delta/gradient # um

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

def get_snp_freq(snp):
    res_freq = snp.frequency.f[np.argmin(snp.s_mag[:, 0, 0])]
    return res_freq

def read_s2p(fname_and_addr):
    snp = rf.Network(fname_and_addr)
    snp.plot_s_db()
    plt.legend(['S11', 'S21', 'S12', 'S22'])
    plt.savefig(f'{savepath}\\{name}_scikit_rf_db.png')
    plt.close('all')

    snp.plot_s_deg()
    plt.legend(['S11', 'S21', 'S12', 'S22'])
    plt.savefig(f'{savepath}\\{name}_scikit_rf_deg.png')
    plt.close('all')

    snp.plot_s_smith()
    plt.legend(['S11', 'S21', 'S12', 'S22'])
    plt.savefig(f'{savepath}\\{name}_scikit_rf_smith.png')
    plt.close('all')

    snp.plot_s_complex()
    plt.legend(['S11', 'S21', 'S12', 'S22'])
    plt.savefig(f'{savepath}\\{name}_scikit_rf_complex.png')
    plt.close('all')

    snp_freq =  get_snp_freq(snp)
    print(f'{snp_freq = }')

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

def get_delta_Db(power_1, power_2):
    power_1 = float(power_1)
    power_2 = float(power_2)
    return abs(10.*np.log10(power_1/power_2))

def get_delta_Db_ratio(power_ratio):
    power_ratio = float(power_ratio)
    return abs(10.*np.log10(power_ratio))

def get_delta_Db_from_MCR_data():
    # C20_CFP = [4.43, 4.7]
    # C20_PP = [1.97, 2.08]
    # C16_CFP = [4.45, 4.75]
    # C16_PP = [2.37, 2.5]

    C20_CFP = [1.648, 2.329]
    C20_PP = [3.144, 4.481]
    C16_CFP = [1.676, 2.302]
    C16_PP = [3.021, 4.206]

    C20_c, C20_m = pmm.best_fit(C20_CFP, C20_PP)
    C16_c, C16_m = pmm.best_fit(C16_CFP, C16_PP)
    
    xs = np.linspace(0.3, 2.4, 1000, endpoint=True)
    C20_fit = [x*C20_m + C20_c for x in xs]
    C16_fit = [x*C16_m + C16_c for x in xs]

    plt.plot(xs, C20_fit, ls='-', lw=0.8, color='g', label='C20')
    plt.plot(xs, C16_fit, ls='-', lw=0.8, color='b', label='C16')
    plt.scatter(C20_CFP, C20_PP, marker='x', s=40, c='k')
    plt.scatter(C16_CFP, C16_PP, marker='x', s=40, c='k')
    plt.xlabel('Cavity Forward Power [MW]')
    plt.ylabel('Cavity Probe Power [MW]')
    plt.ylabel('Cavity Probe Power [MW]')
    plt.legend(loc='best')
    plt.savefig(f'{savepath}\\MCR_CFP_CPP_C20_C16.png')
    plt.close('all')

    delta_c = abs(C20_c - C16_c)
    C20_4p45 = 2.37 - delta_c

    C20_pp = C20_4p45
    C16_pp = 2.37
    delta_Db = get_delta_Db(C16_pp, C20_pp)

    C20_4p75 = 2.5 - delta_c
    C20_pp_2 = C20_4p75
    C16_pp_2 = 2.5
    delta_Db_2 = get_delta_Db(C16_pp_2, C20_pp_2)

    # Db = 10.*np.log10(C16_pp/C20_pp)
    LP_ratio = C20_pp / C16_pp
    HP_ratio = C20_pp_2 / C16_pp_2

    print(f'{delta_c = }')
    print(f'{C20_4p45 = }')
    print(f'{delta_Db = }')
    print(f'{delta_Db_2 = }')
    print(f'{LP_ratio = }')
    print(f'{HP_ratio = }')

    average_delta_Db = np.mean([delta_Db, delta_Db_2])

    return abs(average_delta_Db)


######################
# Db_test = get_delta_Db(237., 197.)
# input(f'{Db_test = }')
#####################

delta_Db = get_delta_Db_from_MCR_data()

print(f'{delta_Db = }')

if read_data:
    data_dict = {}
    name_keys = []
    C16_no_stop_freqs = []
    C16_stop_freqs = []
    C20_no_stop_freqs = []
    C20_stop_freqs = []
    Q_0_list = []
    beta_list = []
    for f in fnames:

        if 'remove' in f:
            print(f'{len(f) = }')
            name = f'{f[:3]}_{f[20]}_removed'
            data_dict[name] = {}
            data_dict[name]['stop_present?'] = False
            data_dict[name]['name'] = name

        else:
            print(f'{len(f) = }')
            name = f'{f[:3]}_{f[20]}'
            data_dict[name] = {}
            data_dict[name]['stop_present?'] = True
            data_dict[name]['name'] = name

        data_dict[name]['fname'] = f
        name_keys.append(name)

        freq, s11, s11_ph, s21, s21_ph, s12, s12_ph, s22, s22_ph = read_s2p(f'{snp_addr}\\{data_dict[name]["fname"]}')
        data_dict[name]['freq_list'] = freq
        data_dict[name]['s11'] = s11
        data_dict[name]['s21'] = s21
        res_freq, min_gamma_s11 = get_response_freq(freq, s11)
        data_dict[name]['res_freq'] = res_freq

        if '16' in name:
            if data_dict[name]['stop_present?']:
                C16_stop_freqs.append(res_freq)
            else:
                C16_no_stop_freqs.append(res_freq)
        elif '20' in name:
            if data_dict[name]['stop_present?']:
                C20_stop_freqs.append(res_freq)
            else:
                C20_no_stop_freqs.append(res_freq)

        data_dict[name]['min_gamma_s11'] = min_gamma_s11
        Q_L, s21_smooth, freq_list_s21_smooth, freq_list, res_freq_s21, max_gamma_s21, max_gamma_raw_s21, min_feq_minus3dB_s21, max_feq_minus3dB_s21, bandwidth, half_delta_f_delta_Db, res_freq_idx_smooth, std_noise_zoom = get_s21_data(freq, s21, delta_Db, name, print_plot=True)
        data_dict[name]['Q_L'] = Q_L
        data_dict[name]['s21_smooth'] = s21_smooth
        data_dict[name]['freq_list_s21_smooth'] = freq_list_s21_smooth
        data_dict[name]['res_freq_s21'] = res_freq_s21
        data_dict[name]['res_freq_idx_smooth'] = res_freq_idx_smooth
        data_dict[name]['max_gamma_s21'] = max_gamma_s21
        data_dict[name]['max_gamma_raw_s21'] = max_gamma_raw_s21
        data_dict[name]['min_feq_minus3dB_s21'] = min_feq_minus3dB_s21
        data_dict[name]['max_feq_minus3dB_s21'] = max_feq_minus3dB_s21
        data_dict[name]['bandwidth'] = bandwidth
        data_dict[name]['delta_f_tune'] = half_delta_f_delta_Db
        temp_freq_gradient = -50000. # -50kHz per deg C
        delta_temp = 1./temp_freq_gradient * half_delta_f_delta_Db
        data_dict[name]['delta_temp_tune'] = delta_temp
        print(f'{half_delta_f_delta_Db = }')
        print(f'{delta_temp = }')
        data_dict[name]['std_noise_zoom'] = std_noise_zoom

        response_at_res_freq_s21 = get_s21_at_res_freq(freq, s21, design_freq_Hz)
        data_dict[name]['response_at_res_freq_s21'] = response_at_res_freq_s21
        data_dict[name]['overcoupled'] = is_overcoupled(s11, s11_ph, name, print_plot=False)
        beta = get_beta(s11, overcoupled=True)
        data_dict[name]['beta'] = beta
        Q_0 = get_Q_0(Q_L, beta)
        data_dict[name]['Q_0'] = Q_0
        design_freq_Hz = design_res_freq_MHz*1e6
        delta_freq = res_freq-design_freq_Hz
        z_coord = get_zcoord_from_freq_delta(delta_freq)
        data_dict[name]['z_coord'] = z_coord

        Q_0_list.append(Q_0)

        beta_list.append(beta)

        print(f'\n{name = }')
        print(f'{res_freq = }')
        print(f'{Q_L = }')
        print(f'{Q_0 = }')
        print(f'{beta = }')

    # pickle dump data for ease of reading in again

    # with open(f'{snp_addr}', 'wb') as f:
    #     pickle.dump(data_dict, f, pickle.HIGHEST_PROTOCOL)

    with open(f'{pkl_addr}\\data_dict.pkl', 'wb') as handle:
        pickle.dump(data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f'{name_keys = }')
    print(f'{C16_no_stop_freqs = }')
    print(f'{C16_stop_freqs = }')
    print(f'{C20_no_stop_freqs = }')
    print(f'{C20_stop_freqs = }')
    print(f'{Q_0_list = }')
    print(f'{beta_list = }')





else:

    # read in data

    # with open(f'{savepath}\\', 'rb') as f:
    #     data_dict =  pickle.load(f)

    with open(f'{pkl_addr}\\data_dict.pkl', 'rb') as handle:
        data_dict = pickle.load(handle)

    name_keys = ['C16_1', 'C16_2', 'C16_3', 'C16_1_removed', 'C16_2_removed', 'C20_1', 'C20_1_removed', 'C20_2',
                 'C20_2_removed', 'C20_3', 'C20_3_removed', 'C20_4_removed', 'C20_5_removed', 'C20_6_removed']
    C16_no_stop_freqs = [2999896249.0623, 2999958764.6912]
    C16_stop_freqs = [2999983770.9427, 2999951262.8157, 3000058789.6974]
    C20_no_stop_freqs = [3000081295.3238, 2999833733.4334, 2999861240.3101, 2999848737.1843, 2999928757.1893,
                         2999833733.4334]
    C20_stop_freqs = [3000128807.2018, 2999866241.5604, 2999881245.3113]
    Q_0_list = [12547.97994807139, 12617.657297656851, 12537.330040632458, 12550.004752512048, 12546.808700265308,
                12583.56071336736, 12458.871458176514, 12545.365317140619, 12419.932268171007, 12478.660109259128,
                12481.326229441232, 12484.358002160503, 12540.965629097118, 12485.258396662344]
    beta_list = [1.0186490215817914, 1.0193628061135411, 1.016888651120073, 1.0190336643521531, 1.0184757397847044,
                 1.0137900760029726, 1.0146366560141389, 1.018310829336793, 1.018857251884908, 1.0179677372548406,
                 1.0184123448471052, 1.0189127242335867, 1.0175609581821787, 1.019068430258647]


'''Plots'''

circle_size = 15
cross_size = 50
plus_size = 80
square_size = 20
C16_colour = 'b'
C20_colour = 'r'

# hitogram of peak max smoothed s21

all_max_s21 = []
for name in name_keys:
    max_gamma_s21 = data_dict[name]['max_gamma_s21']
    all_max_s21.append(max_gamma_s21)
mean_peak_s21 = np.mean(all_max_s21)
std_peak_s21 = np.std(all_max_s21)
n, bedges, patches = plt.hist(all_max_s21, bins=9, histtype='step', lw=0.8, color='k')
plt.vlines(mean_peak_s21, 0., max(n), ls='--', lw=0.8, color='r')
plt.text(-69.68, 3.5, r'$\mu$'f' = {mean_peak_s21:1.4f}\n'r'$\sigma$'f' = {std_peak_s21:1.4f}')
plt.ylabel('N')
plt.xlabel('Peak S21 [Db]')
plt.savefig(f'{savepath}\\all_max_s21.png')
plt.close('all')

# s21_zoom_overplot
legend_done_16_guard_present = False
legend_done_16_guard_removed = False
legend_done_20_guard_present = False
legend_done_20_guard_removed = False

STD_noise = []
for name in name_keys:

    res_freq_s21 = data_dict[name]['res_freq_s21']
    freq_list_s21_smooth = data_dict[name]['freq_list_s21_smooth']
    s21_smooth = data_dict[name]['s21_smooth']
    res_freq_idx_smooth = data_dict[name]['res_freq_idx_smooth']
    buffer_idx = 25
    min_idx_zoom = res_freq_idx_smooth - buffer_idx
    max_idx_zoom = res_freq_idx_smooth + buffer_idx
    freq_list_s21_smooth_zoom = freq_list_s21_smooth[min_idx_zoom:max_idx_zoom]
    freq_list_s21_smooth_zoom_centred = [f - res_freq_s21 for f in freq_list_s21_smooth_zoom]
    s21_smooth_zoom = s21_smooth[min_idx_zoom:max_idx_zoom]

    std_noise_zoom = data_dict[name]['std_noise_zoom']
    STD_noise.append(std_noise_zoom)
    #
    # freq_list = data_dict[name]['freq_list']
    # s11= data_dict[name]['s11']
    # freq_list_GHz = [i/1e9 for i in freq_list]
    # s21 = data_dict[name]['s21']
    # res_freq = data_dict[name]['res_freq']
    # res_freq_GHz = res_freq/1e9
    #
    # max_gamma_s21 = data_dict[name]['max_gamma_s21']
    # max_gamma_raw_s21 = data_dict[name]['max_gamma_raw_s21']
    # min_feq_minus3dB_s21 = data_dict[name]['min_feq_minus3dB_s21']
    # max_feq_minus3dB_s21 = data_dict[name]['max_feq_minus3dB_s21']
    # bandwidth = data_dict[name]['bandwidth']
    # Q_0 = data_dict[name]['Q_0']
    # Q_L = data_dict[name]['Q_L']
    # beta = data_dict[name]['beta']
    # z_coord = data_dict[name]['z_coord']
    #
    # FF = get_FF_from_zcoord(z_coord)

    size = 1

    if '16' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_16_guard_present:
                plt.plot(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, ls='-', lw=size, color='b')
                # plt.scatter(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, marker='.', s=size, color='b')
            else:
                plt.plot(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, ls='-', lw=size, color='b', label='C16 guard present')
                # plt.scatter(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, marker='.', s=size, color='b',
                #             label='C16 guard present')
                legend_done_16_guard_present = True
        else:
            if legend_done_16_guard_removed:
                plt.plot(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, ls='--', lw=size, color='b')
                # plt.scatter(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, marker='.', s=size, color='b')
            else:
                plt.plot(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, ls='--', lw=size, color='b',  label='C16 guard removed')
                # plt.scatter(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, marker='.', s=size, color='b',
                #             label='C16 guard removed')
                legend_done_16_guard_removed = True

    elif '20' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_20_guard_present:
                plt.plot(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, ls='-', lw=size, color='r')
                #plt.scatter(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, marker='.', s=size, color='r')
            else:
                plt.plot(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, ls='-', lw=size, color='r',  label='C20 guard present')
                # plt.scatter(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, marker='.', s=size, color='r',
                #             label='C20 guard present')
                legend_done_20_guard_present = True
        else:
            if legend_done_20_guard_removed:
                plt.plot(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, ls='--', lw=size, color='r')
                #plt.scatter(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, marker='.', s=size, color='r')
            else:
                plt.plot(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, ls='--', lw=size, color='r', label='C20 guard removed')
                # plt.scatter(freq_list_s21_smooth_zoom_centred, s21_smooth_zoom, marker='.', s=size, color='r',
                #             label='C20 guard removed')
                legend_done_20_guard_removed = True

mean_STD_noise = np.mean(STD_noise)
plt.text(-60000., -69.6, f'{mean_STD_noise = :1.3f} Db')
plt.ylim(-70., -69.55)
plt.xlabel(r'$\Delta$''f [Hz]')
plt.ylabel(r'S21 [Db]')
plt.legend(loc='lower center')
plt.savefig(f'{savepath}\\S21_overplot.png')
plt.close('all')

# histogram of the noise at the peak s21 region
nbins = 9
# plt.hist(STD_noise, bins=int(np.sqrt(len(STD_noise))), histtype='step', lw=0.8, color='k')
n, bedges, patches = plt.hist(STD_noise, bins=nbins, histtype='step', lw=1., color='k')
plt.vlines(mean_STD_noise, 0, max(n), ls='-', lw=0.8, color='r', label=f'{mean_STD_noise = :1.3f}')
plt.ylabel('N')
plt.xlabel(r'$\sigma$'' [Db]')
plt.legend(loc='upper right')
plt.savefig(f'{savepath}\\STD_noise_hist.png')
plt.close('all')

# beta & Q_0 histograms

nbins = 5

n, bedges, patches = plt.hist(beta_list, bins = nbins, histtype='step', lw=0.8, color='k')
plt.vlines(np.mean(beta_list), 0., max(n), ls='--', lw=0.8, color='r')
plt.text(min(bedges), max(n), f'mean = {np.mean(beta_list):1.3f}')
plt.xlabel(r'$\beta$')
plt.ylabel('N')

plt.savefig(f'{savepath}\\Beta_histpng')
plt.close('all')

n, bins, patches = plt.hist(Q_0_list, bins = nbins, histtype='step', lw=0.8, color='k')
plt.vlines(np.mean(Q_0_list), 0., max(n), ls='--', lw=0.8, color='r')
plt.text(min(bins), max(n), f'mean = {np.mean(Q_0_list):1.3f}')
plt.xlabel(r'$Q_0$')
plt.ylabel('N')
plt.savefig(f'{savepath}\\Q_0_hist.png')
plt.close('all')

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
legend_done_16_guard_present = False
legend_done_16_guard_removed = False
legend_done_20_guard_present = False
legend_done_20_guard_removed = False
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

    # print(f'{str(data_dict[name]["name"]) = }')
    if '16' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_16_guard_present:
                plt.scatter(z_coord, beta, marker='x', s=cross_size, color='b')
            else:
                plt.scatter(z_coord, beta, marker='x', s=cross_size, color='b', label='C16 guard present')
                legend_done_16_guard_present = True
        else:
            if legend_done_16_guard_removed:
                plt.scatter(z_coord, beta, marker='o', s=circle_size, color='b')
            else:
                plt.scatter(z_coord, beta, marker='o', s=circle_size, color='b', label='C16 guard removed')
                legend_done_16_guard_removed = True

    elif '20' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_20_guard_present:
                plt.scatter(z_coord, beta, marker='x', s=cross_size, color='r')
            else:
                plt.scatter(z_coord, beta, marker='x', s=cross_size, color='r', label='C20 guard present')
                legend_done_20_guard_present = True
        else:
            if legend_done_20_guard_removed:
                plt.scatter(z_coord, beta, marker='o', s=circle_size, color='r')
            else:
                plt.scatter(z_coord, beta, marker='o', s=circle_size, color='r', label='C20 guard removed')
                legend_done_20_guard_removed = True


plt.xlabel('Z-coordinate ['r'$\mu$''m]')
plt.ylabel(r'$\beta$')
plt.legend(loc='upper left')
plt.savefig(f'{savepath}\\zcoord_beta.png')
plt.close('all')

# Q_0
legend_done_16_guard_present = False
legend_done_16_guard_removed = False
legend_done_20_guard_present = False
legend_done_20_guard_removed = False
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

    # print(f'{str(data_dict[name]["name"]) = }')
    if '16' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_16_guard_present:
                plt.scatter(z_coord, Q_0, marker='x', s=cross_size, color='b')
            else:
                plt.scatter(z_coord, Q_0, marker='x', s=cross_size, color='b', label='C16 guard present')
                legend_done_16_guard_present = True
        else:
            if legend_done_16_guard_removed:
                plt.scatter(z_coord, Q_0, marker='o', s=circle_size, color='b')
            else:
                plt.scatter(z_coord, Q_0, marker='o', s=circle_size, color='b', label='C16 guard removed')
                legend_done_16_guard_removed = True

    elif '20' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_20_guard_present:
                plt.scatter(z_coord, Q_0, marker='x', s=cross_size, color='r')
            else:
                plt.scatter(z_coord, Q_0, marker='x', s=cross_size, color='r', label='C20 guard present')
                legend_done_20_guard_present = True
        else:
            if legend_done_20_guard_removed:
                plt.scatter(z_coord, Q_0, marker='o', s=circle_size, color='r')
            else:
                plt.scatter(z_coord, Q_0, marker='o', s=circle_size, color='r', label='C20 guard removed')
                legend_done_20_guard_removed = True

plt.xlabel('Z-coordinate ['r'$\mu$''m]')
plt.ylabel(r'$Q_0$')
plt.legend(loc='lower center')
plt.savefig(f'{savepath}\\zcoord_Q_0.png')
plt.close('all')

# Q_L

legend_done_16_guard_present = False
legend_done_16_guard_removed = False
legend_done_20_guard_present = False
legend_done_20_guard_removed = False
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

    # print(f'{str(data_dict[name]["name"]) = }')
    if '16' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_16_guard_present:
                plt.scatter(z_coord, Q_L, marker='x', s=cross_size, color='b')
            else:
                plt.scatter(z_coord, Q_L, marker='x', s=cross_size, color='b', label='C16 guard present')
                legend_done_16_guard_present = True
        else:
            if legend_done_16_guard_removed:
                plt.scatter(z_coord, Q_L, marker='o', s=circle_size, color='b')
            else:
                plt.scatter(z_coord, Q_L, marker='o', s=circle_size, color='b', label='C16 guard removed')
                legend_done_16_guard_removed = True

    elif '20' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_20_guard_present:
                plt.scatter(z_coord, Q_L, marker='x', s=cross_size, color='r')
            else:
                plt.scatter(z_coord, Q_L, marker='x', s=cross_size, color='r', label='C20 guard present')
                legend_done_20_guard_present = True
        else:
            if legend_done_20_guard_removed:
                plt.scatter(z_coord, Q_L, marker='o', s=circle_size, color='r')
            else:
                plt.scatter(z_coord, Q_L, marker='o', s=circle_size, color='r', label='C20 guard removed')
                legend_done_20_guard_removed = True

plt.xlabel('Z-coordinate ['r'$\mu$''m]')
plt.ylabel(r'$Q_L$')
plt.legend(loc='lower center')
plt.savefig(f'{savepath}\\zcoord_Q_L.png')
plt.close('all')



legend_done_16_guard_present = False
legend_done_16_guard_removed = False
legend_done_20_guard_present = False
legend_done_20_guard_removed = False

all_res_freq_data_z_coords = []
all_res_freq_data_s21 = []


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
    response_at_res_freq_s21 = data_dict[name]['response_at_res_freq_s21']
    all_res_freq_data_z_coords.append(z_coord)
    all_res_freq_data_s21.append(response_at_res_freq_s21)



    # print(f'{str(data_dict[name]["name"]) = }')
    if '16' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_16_guard_present:
                plt.scatter(z_coord, response_at_res_freq_s21, marker='x', s=cross_size, color='b')
            else:
                plt.scatter(z_coord, response_at_res_freq_s21, marker='x', s=cross_size, color='b', label='C16 guard present')
                legend_done_16_guard_present = True
        else:
            if legend_done_16_guard_removed:
                plt.scatter(z_coord, response_at_res_freq_s21, marker='o', s=circle_size, color='b')
            else:
                plt.scatter(z_coord, response_at_res_freq_s21, marker='o', s=circle_size, color='b', label='C16 guard removed')
                legend_done_16_guard_removed = True

    elif '20' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_20_guard_present:
                plt.scatter(z_coord, response_at_res_freq_s21, marker='x', s=cross_size, color='r')
            else:
                plt.scatter(z_coord, response_at_res_freq_s21, marker='x', s=cross_size, color='r', label='C20 guard present')
                legend_done_20_guard_present = True
        else:
            if legend_done_20_guard_removed:
                plt.scatter(z_coord, response_at_res_freq_s21, marker='o', s=circle_size, color='r')
            else:
                plt.scatter(z_coord, response_at_res_freq_s21, marker='o', s=circle_size, color='r', label='C20 guard removed')
                legend_done_20_guard_removed = True

all_res_freq_data_z_coords_2, all_res_freq_data_s21_2 = zip(*sorted(zip(all_res_freq_data_z_coords, all_res_freq_data_s21)))
for idx in range(len(all_res_freq_data_z_coords)):
    print(all_res_freq_data_z_coords[idx],    all_res_freq_data_z_coords_2[idx])


# interp_x, interp_y = pmm.cubic_spline_interpolation(all_res_freq_data_z_coords_2[:-1], all_res_freq_data_s21_2[:-1], 10000)
smooth_x, smooth_y = pmm.smooth(all_res_freq_data_z_coords_2, all_res_freq_data_s21_2, box_pts=100, clip_index=100)

plt.plot(smooth_x, smooth_y, ls='-', lw=0.8, color='k', label='smoothed data')
plt.xlabel('Z-coordinate ['r'$\mu$''m]')
plt.ylabel(r'$\Gamma_{MAX}$''[dB]')
plt.legend(loc=8)
plt.savefig(f'{savepath}\\zcoord_response_at_res_freq_s21.png')
plt.close('all')

# max_gamma_raw_s21
legend_done_16_guard_present = False
legend_done_16_guard_removed = False
legend_done_20_guard_present = False
legend_done_20_guard_removed = False
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

    if '16' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_16_guard_present:
                plt.scatter(z_coord, max_gamma_raw_s21, marker='x', s=cross_size, color='b')
            else:
                plt.scatter(z_coord, max_gamma_raw_s21, marker='x', s=cross_size, color='b', label='C16 guard present')
                legend_done_16_guard_present = True
        else:
            if legend_done_16_guard_removed:
                plt.scatter(z_coord, max_gamma_raw_s21, marker='o', s=circle_size, color='b')
            else:
                plt.scatter(z_coord, max_gamma_raw_s21, marker='o', s=circle_size, color='b', label='C16 guard removed')
                legend_done_16_guard_removed = True

    elif '20' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_20_guard_present:
                plt.scatter(z_coord, max_gamma_raw_s21, marker='x', s=cross_size, color='r')
            else:
                plt.scatter(z_coord, max_gamma_raw_s21, marker='x', s=cross_size, color='r', label='C20 guard present')
                legend_done_20_guard_present = True
        else:
            if legend_done_20_guard_removed:
                plt.scatter(z_coord, max_gamma_raw_s21, marker='o', s=circle_size, color='r')
            else:
                plt.scatter(z_coord, max_gamma_raw_s21, marker='o', s=circle_size, color='r', label='C20 guard removed')
                legend_done_20_guard_removed = True

plt.xlabel('Z-coordinate ['r'$\mu$''m]')
plt.ylabel('raw ''$\Gamma_{MAX}$''[dB]')
# plt.legend(loc='upper left')
plt.legend(loc=8)
plt.savefig(f'{savepath}\\zcoord_max_gamma_raw_s21.png')
plt.close('all')


# max_gamma_s21
legend_done_16_guard_present = False
legend_done_16_guard_removed = False
legend_done_20_guard_present = False
legend_done_20_guard_removed = False
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

    if '16' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_16_guard_present:
                plt.scatter(z_coord, max_gamma_s21, marker='x', s=cross_size, color='b')
            else:
                plt.scatter(z_coord, max_gamma_s21, marker='x', s=cross_size, color='b', label='C16 guard present')
                legend_done_16_guard_present = True
        else:
            if legend_done_16_guard_removed:
                plt.scatter(z_coord, max_gamma_s21, marker='o', s=circle_size, color='b')
            else:
                plt.scatter(z_coord, max_gamma_s21, marker='o', s=circle_size, color='b', label='C16 guard removed')
                legend_done_16_guard_removed = True

    elif '20' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_20_guard_present:
                plt.scatter(z_coord, max_gamma_s21, marker='x', s=cross_size, color='r')
            else:
                plt.scatter(z_coord, max_gamma_s21, marker='x', s=cross_size, color='r', label='C20 guard present')
                legend_done_20_guard_present = True
        else:
            if legend_done_20_guard_removed:
                plt.scatter(z_coord, max_gamma_s21, marker='o', s=circle_size, color='r')
            else:
                plt.scatter(z_coord, max_gamma_s21, marker='o', s=circle_size, color='r', label='C20 guard removed')
                legend_done_20_guard_removed = True

plt.xlabel('Z-coordinate ['r'$\mu$''m]')
plt.ylabel('$\Gamma_{MAX}$''[dB]')
plt.legend(loc='upper left')
# plt.legend(loc=8)
plt.savefig(f'{savepath}\\zcoord_max_gamma_s21.png')
plt.close('all')

# insertion frequencies#
for name in name_keys:
    print(name)
    print(type(name))
    print(data_dict[name]['stop_present?'])
C20_guard_removed_freqs = [data_dict[str(n)]['res_freq'] for n in name_keys if 'C20' in str(data_dict[n]['name']) and str(data_dict[n]['stop_present?']) == 'False']
C20_guard_present_freqs = [data_dict[str(n)]['res_freq'] for n in name_keys if 'C20' in str(data_dict[n]['name']) and str(data_dict[n]['stop_present?']) == 'True']
C16_guard_removed_freqs = [data_dict[str(n)]['res_freq'] for n in name_keys if 'C16' in str(data_dict[n]['name']) and str(data_dict[n]['stop_present?']) == 'False']
C16_guard_present_freqs = [data_dict[str(n)]['res_freq'] for n in name_keys if 'C16' in str(data_dict[n]['name']) and str(data_dict[n]['stop_present?']) == 'True']
C20_guard_removed_freqs = [(design_freq_Hz - i)/1e3 for i in C20_guard_removed_freqs]
C20_guard_present_freqs = [(design_freq_Hz - i)/1e3 for i in C20_guard_present_freqs]
C16_guard_removed_freqs = [(design_freq_Hz - i)/1e3 for i in C16_guard_removed_freqs]
C16_guard_present_freqs = [(design_freq_Hz - i)/1e3 for i in C16_guard_present_freqs]


plt.scatter([1]*len(C16_guard_present_freqs), C16_guard_present_freqs, marker='x', s=cross_size, color='b', label='C16 guard present')
plt.scatter([2]*len(C16_guard_removed_freqs), C16_guard_removed_freqs, marker='o', s=cross_size, color='b', label='C16 guard removed')
plt.scatter([3]*len(C20_guard_present_freqs), C20_guard_present_freqs, marker='x', s=cross_size, color='r', label='C20 guard present')
plt.scatter([4]*len(C20_guard_removed_freqs), C20_guard_removed_freqs, marker='o', s=cross_size, color='r', label='C20 guard removed')

plt.ylabel(r'$\Delta$''f [kHz]')
plt.ylim(-200, 175)
plt.legend(loc='lower left')
plt.savefig(f'{savepath}\\Insertion_frequencies.png')
plt.close('all')


# z-coord vs freq

C16_freq_delta_kHz = 44.999999
C20_freq_delta_kHz = -44.999999
print(f'{C20_freq_delta_kHz = }')
print(f'{C16_freq_delta_kHz = }')
C16_zcoord_um = (1./-1.437) * C16_freq_delta_kHz # gradient = Hz/um
C20_zcoord_um = (1./-1.437) * C20_freq_delta_kHz # gradient = Hz/um

sim_z_coords = np.linspace(-200, 200, 1000, endpoint=True)
gradient = -1.4375
intercept = -5.
sim_freqs = [z * gradient for z in sim_z_coords]

plt.plot(sim_z_coords, sim_freqs, ls='--', lw=0.8, color='k', label='CST Sim')
plt.scatter(C16_zcoord_um, C16_freq_delta_kHz, marker='s', s=square_size, color='magenta', label='C16 from run')
plt.scatter(C20_zcoord_um, C20_freq_delta_kHz, marker='s', s=square_size, color='cyan', label='C20 from run')

legend_done_16_guard_present = False
legend_done_16_guard_removed = False
legend_done_20_guard_present = False
legend_done_20_guard_removed = False
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

    res_freq = (res_freq - design_freq_Hz) / 1e3 # um kHz

    if '16' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_16_guard_present:
                plt.scatter(z_coord, res_freq, marker='+', s=plus_size, color='b')
            else:
                plt.scatter(z_coord, res_freq, marker='+', s=plus_size, color='b', label='C16 guard present')
                legend_done_16_guard_present = True
        else:
            if legend_done_16_guard_removed:
                plt.scatter(z_coord, res_freq, marker='o', s=circle_size, color='b')
            else:
                plt.scatter(z_coord, res_freq, marker='o', s=circle_size, color='b', label='C16 guard removed')
                legend_done_16_guard_removed = True

    elif '20' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_20_guard_present:
                plt.scatter(z_coord, res_freq, marker='+', s=plus_size, color='r')
            else:
                plt.scatter(z_coord, res_freq, marker='+', s=plus_size, color='r', label='C20 guard present')
                legend_done_20_guard_present = True
        else:
            if legend_done_20_guard_removed:
                plt.scatter(z_coord, res_freq, marker='o', s=circle_size, color='r')
            else:
                plt.scatter(z_coord, res_freq, marker='o', s=circle_size, color='r', label='C20 guard removed')
                legend_done_20_guard_removed = True

plt.xlabel('Z-coordinate ['r'$\mu$''m]')
plt.ylabel(r'$\Delta$''f [kHz]')
plt.legend(loc='upper right')
plt.savefig(f'{savepath}\\zcoord_freqs_sim_line.png')
plt.close('all')

bad_insertion_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\HRRG_measurements_20230418\raw_data'
bad_insertion_fnames = os.listdir(bad_insertion_addr)

bad_insertion_fnames_addr = [f'{bad_insertion_addr}\\{i}' for i in bad_insertion_fnames if 'out' in i or 'bad' in i]
cathode_out_fnames_addr = [f'{bad_insertion_addr}\\{i}' for i in bad_insertion_fnames if 'out' in i]
cathode_out_res_freq = res_freq = get_snp_freq(rf.Network(cathode_out_fnames_addr[0]))
print(f'{cathode_out_res_freq = }')
print(bad_insertion_fnames)

bad_insertion_frequencies = []
for snp_file in bad_insertion_fnames_addr:
    snp = rf.Network(snp_file)
    res_freq = get_snp_freq(snp)
    print(res_freq)
    bad_insertion_frequencies.append(res_freq)

bad_insertion_freqs_GHz = [i/1e9 for i in bad_insertion_frequencies]

# C16_no_stop_freqs = []
# C16_stop_freqs = []
# C20_no_stop_freqs = []
# C20_stop_freqs = []


plt.scatter([0]*len(bad_insertion_frequencies), bad_insertion_frequencies, marker= 's', s = square_size, color=C16_colour, label='C16 bad insertion')
plt.scatter(0, cathode_out_res_freq, marker= 's', s = square_size, color='purple', label='Cathode out')
plt.scatter([2]*len(C16_no_stop_freqs), C16_no_stop_freqs, marker= 'o', s = circle_size, color=C16_colour, label='C16 guard removed')
plt.scatter([1]*len(C16_stop_freqs), C16_stop_freqs, marker= 'x', s = cross_size, color=C16_colour, label='C16 guard present')
plt.scatter([4]*len(C20_no_stop_freqs), C20_no_stop_freqs, marker= 'o', s = circle_size, color=C20_colour, label='C20 guard removed')
plt.scatter([3]*len(C20_stop_freqs), C20_stop_freqs, marker= 'x', s = cross_size, color=C20_colour, label='C20 guard present')

plt.ylabel('Frequency [GHz]')
plt.legend(loc='upper right')
plt.savefig(f'{savepath}\\Freqs_incl_bad_insertions.png')
plt.close('all')


# field flatness

HRRG_Field_Flatness = [96.60317211254477, 95.16144487001026, 93.89552062661363, 93.33856511314178, 92.38771204756279, 98.16246497556632, 99.96578947538544, 100.92562161509852, 102.91949593981819]

legend_done_16_guard_present = False
legend_done_16_guard_removed = False
legend_done_20_guard_present = False
legend_done_20_guard_removed = False
zcoords = [0., -100., -200., -250., -350., 100., 200., 250., 350.]
ff_c_hrrg, ff_m_hrrg = pmm.best_fit(zcoords, HRRG_Field_Flatness)
ff_x_fit_hrrg = np.linspace(min(zcoords), max(zcoords), 1000)


y_poly_hrrg = np.polyfit(zcoords, HRRG_Field_Flatness, 2)
print(f'{y_poly_hrrg = }')
A_poly_hrrg = float(y_poly_hrrg[0])
B_poly_hrrg = float(y_poly_hrrg[1])
C_poly_hrrg = float(y_poly_hrrg[2])

print(f'{A_poly_hrrg = }')
print(f'{B_poly_hrrg = }')
print(f'{C_poly_hrrg = }')

y_fit_poly_hrrg = np.poly1d(y_poly_hrrg)

C16_FF = get_FF_from_zcoord(C16_zcoord_um)
C20_FF = get_FF_from_zcoord(C20_zcoord_um)

plt.scatter(C16_zcoord_um, C16_FF, marker='s', s=square_size, color='magenta', label='C16 from run')
plt.scatter(C20_zcoord_um, C20_FF, marker='s', s=square_size, color='cyan', label='C20 from run')
plt.plot(ff_x_fit_hrrg, y_fit_poly_hrrg(ff_x_fit_hrrg), ls='--', lw=0.8, color='k', label='fit')
plt.scatter(zcoords, HRRG_Field_Flatness, marker='o', s=15, color='k', label='sim')

# zcoord_field_flatness_sim_line

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

    FF = get_FF_from_zcoord(z_coord)



    if '16' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_16_guard_present:
                plt.scatter(z_coord, FF, marker='+', s=plus_size, color='b')
            else:
                plt.scatter(z_coord, FF, marker='+', s=plus_size, color='b', label='C16 guard present')
                legend_done_16_guard_present = True
        else:
            if legend_done_16_guard_removed:
                plt.scatter(z_coord, FF, marker='o', s=circle_size, color='b')
            else:
                plt.scatter(z_coord, FF, marker='o', s=circle_size, color='b', label='C16 guard removed')
                legend_done_16_guard_removed = True

    elif '20' in str(data_dict[name]['name']):
        if data_dict[name]['stop_present?']:
            if legend_done_20_guard_present:
                plt.scatter(z_coord, FF, marker='+', s=plus_size, color='r')
            else:
                plt.scatter(z_coord, FF, marker='+', s=plus_size, color='r', label='C20 guard present')
                legend_done_20_guard_present = True
        else:
            if legend_done_20_guard_removed:
                plt.scatter(z_coord, FF, marker='o', s=circle_size, color='r')
            else:
                plt.scatter(z_coord, FF, marker='o', s=circle_size, color='r', label='C20 guard removed')
                legend_done_20_guard_removed = True

plt.xlabel('Z-coordinate ['r'$\mu$''m]')
plt.ylabel(r'Field Flatness [%]')
plt.legend(loc='upper left')
plt.savefig(f'{savepath}\\zcoord_field_flatness_sim_line.png')
plt.close('all')





