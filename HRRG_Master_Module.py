import numpy as np
from matplotlib import pyplot as plt
import sys
pmm_addr = r'C:\Users\zup98752\PycharmProjects\PhD'
sys.path.insert(1,pmm_addr)
import PhD_Master_Module as pmm
import skrf as rf

f = 2998.538e6
mu0 = 4. * np.pi * 1.e-7
w = 2. * np.pi * f
sigma = 5.8e7
Rs = np.sqrt(w * mu0 / (2. * sigma))
Q = 13199.
G = Q * Rs
print(f'{G = }')

# exit()

def get_FF_from_zcoord(zcoord):
    A_poly_hrrg = 8.725565888064103e-06
    B_poly_hrrg = 0.015100094187422454
    C_poly_hrrg = 96.58430786759465
    gradient = 0.051
    FF = A_poly_hrrg * zcoord ** 2. + B_poly_hrrg * zcoord + C_poly_hrrg

    return FF


def get_response_freq(freq_list, s11, print_plot=False):
    min_gamma, min_idx = pmm.get_min_val_idx_from_list(s11)
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
    noise = [s21[i + 1] - s21[i] for i in range(len(s21) - 1)]
    freq_list_noise = freq_list[:-1]

    return freq_list_noise, noise


def get_s21_data(freq_list, s21, delta_Db, savepath, savename, print_plot=False):
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

    res_freq_idx_smooth = res_freq_idx  # + clip_idx

    minus3dB_gamma = max_gamma - 3.
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
    freq_list_noise_centred = [i - res_freq for i in freq_list_noise]
    freq_list_noise_centred_zoom = freq_list_noise_centred[
                                   min_idx_above_minus_delta_dB + clip_idx:max_idx_above_minus_delta_dB + clip_idx]
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

    # plt.hist(noise_zoom, bins=int(np.sqrt(len(noise_zoom))), histtype='step', lw=0.8, color='k')
    # plt.ylabel('N')
    # plt.xlabel(r'$\Delta Db$')
    # plt.savefig(f'{savepath}\\{savename}_noise_zoom_hist.png')
    # plt.close('all')

    std_noise_zoom = np.std(noise_zoom)
    mean_noise_zoom = np.mean(noise_zoom)

    half_delta_f_delta_Db = (max_feq_minus_delta_dB - min_feq_minus_delta_dB) / 2.

    print(f'{delta_Db = }')
    print(f'{min_feq_minus_delta_dB = }')
    print(f'{max_feq_minus_delta_dB = }')
    print(f'{half_delta_f_delta_Db = }')
    print(f'{min_feq_minus3dB = }')
    print(f'{max_feq_minus3dB = }')
    half_bandwidth = bandwidth / 2.
    print(f'{half_bandwidth = }')
    print(f'{len(freq_list_smooth) = }')
    print(f'{len(s21_smooth) = }')

    Q_L = res_freq / bandwidth

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
        # print(f'{Q_L = }')

    return Q_L, s21_smooth, freq_list_smooth, freq_list, res_freq, max_gamma, max_gamma_raw, min_feq_minus3dB, max_feq_minus3dB, bandwidth, half_delta_f_delta_Db, res_freq_idx_smooth, std_noise_zoom


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
    Q_0 = (1. + beta) * Q_L

    if print_plot:
        print(f'{Q_0 = }')

    return Q_0


def get_zcoord_from_freq_delta(freq_delta):
    gradient = -1437.  # Hz/um

    # intercept = -5.0000

    z = freq_delta / gradient  # um

    return z


def is_overcoupled(s11, s11_ph, name, print_plot=False):
    s11_ph_unwrapped = s11_ph
    # s11_ph_unwrapped = [i if i>0. else -i for i in s11_ph]
    s11_ph_radians = [i / (180. * np.pi) for i in s11_ph_unwrapped]
    re_s11 = [s11[i] * np.cos(s11_ph_radians[i]) for i in range(len(s11))]
    im_s11 = [s11[i] * np.sin(s11_ph_radians[i]) for i in range(len(s11))]
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


def read_s2p(fname_and_addr, savepath, name):
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

    snp_freq = get_snp_freq(snp)
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

                # if float(row[0]) > 2.99:
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
    return abs(10. * np.log10(power_1 / power_2))


def get_delta_Db_ratio(power_ratio):
    power_ratio = float(power_ratio)
    return abs(10. * np.log10(power_ratio))


def get_delta_Db_from_MCR_data(savepath):
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
    C20_fit = [x * C20_m + C20_c for x in xs]
    C16_fit = [x * C16_m + C16_c for x in xs]

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

