import numpy as np
import sys, os, csv, pickle
from matplotlib import pyplot as plt
pmm_addr = r'C:\Users\zup98752\PycharmProjects\PhD'
sys.path.insert(1,pmm_addr)
import HRRG_Master_Module as hmm


start_freq_Hz = 2988500000.
end_freq_Hz = 3008500000.
span = end_freq_Hz - start_freq_Hz
print(f"{start_freq_Hz = }")
print(f"{end_freq_Hz = }")
print(f"Span = {span/1e6} MHz")
fout = 2.9996816
fin = 2.9993384
fdiff = fout-fin
print(f"{fdiff = }")

read_data = True
delta_Db = 0.
design_freq_Hz = 2998500000.
design_temp = 50.
design_freq_MHz = design_freq_Hz/1.e6
design_freq_GHz = design_freq_Hz/1.e9
print(f'{design_freq_MHz = } MHz')
print(f'{design_freq_GHz = } GHz')

savepath = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\Measurements\HRRG_measurements_20240517\analysis'
snp_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\Measurements\HRRG_measurements_20240517\data'
pkl_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\Measurements\HRRG_measurements_20240517'
fnames = os.listdir(snp_addr)
###################
# just look at good insertions
fnames = [f for f in fnames if "in" in f][1:]
###################
fname = fnames[0]
fname_and_addr = f'{snp_addr}\\{fname}'
print(f'{fnames = }')

if read_data:
    data_dict = {}
    name_keys = []
    Q_0_list = []
    Q_L_list = []
    beta_list = []
    res_freq_list = []
    min_gamma_s11_list = []
    res_freq_s21_list = []
    max_gamma_s21_list = []
    max_gamma_raw_s21_list = []
    min_feq_minus3dB_s21_list = []
    max_feq_minus3dB_s21_list = []
    bandwidth_list = []
    response_at_res_freq_s21_list = []
    overcoupled_list = []
    z_coord_list = []
    for f in fnames:

        print(f'{len(f) = }')
        name = f'{f[:-4]}'
        print(f'{name = }')
        data_dict[name] = {}
        data_dict[name]['name'] = name

        data_dict[name]['fname'] = f
        name_keys.append(name)

        freq, s11, s11_ph, s21, s21_ph, s12, s12_ph, s22, s22_ph = hmm.read_s2p(f'{snp_addr}\\{data_dict[name]["fname"]}', savepath, name)
        data_dict[name]['freq_list'] = freq
        data_dict[name]['s11'] = s11
        data_dict[name]['s21'] = s21
        res_freq, min_gamma_s11 = hmm.get_response_freq(freq, s11)
        data_dict[name]['res_freq'] = res_freq
        res_freq_list.append(res_freq)
        data_dict[name]['min_gamma_s11'] = min_gamma_s11
        min_gamma_s11_list.append(min_gamma_s11)
        Q_L, s21_smooth, freq_list_s21_smooth, freq_list, res_freq_s21, max_gamma_s21, max_gamma_raw_s21, min_feq_minus3dB_s21, max_feq_minus3dB_s21, bandwidth, half_delta_f_delta_Db, res_freq_idx_smooth, std_noise_zoom = hmm.get_s21_data(
            freq, s21, delta_Db, savepath, name, print_plot=True)
        data_dict[name]['Q_L'] = Q_L
        Q_L_list.append(Q_L)
        data_dict[name]['s21_smooth'] = s21_smooth
        data_dict[name]['freq_list_s21_smooth'] = freq_list_s21_smooth
        data_dict[name]['res_freq_s21'] = res_freq_s21
        res_freq_s21_list.append(res_freq_s21)
        data_dict[name]['res_freq_idx_smooth'] = res_freq_idx_smooth
        data_dict[name]['max_gamma_s21'] = max_gamma_s21
        max_gamma_s21_list.append(max_gamma_s21)
        data_dict[name]['max_gamma_raw_s21'] = max_gamma_raw_s21
        max_gamma_raw_s21_list.append(max_gamma_raw_s21)
        data_dict[name]['min_feq_minus3dB_s21'] = min_feq_minus3dB_s21
        min_feq_minus3dB_s21_list.append(min_feq_minus3dB_s21)
        data_dict[name]['max_feq_minus3dB_s21'] = max_feq_minus3dB_s21
        max_feq_minus3dB_s21_list.append(max_feq_minus3dB_s21)
        data_dict[name]['bandwidth'] = bandwidth
        bandwidth_list.append(bandwidth)
        data_dict[name]['delta_f_tune'] = half_delta_f_delta_Db
        temp_freq_gradient = -50000.  # -50kHz per deg C
        delta_temp = 1. / temp_freq_gradient * half_delta_f_delta_Db
        data_dict[name]['delta_temp_tune'] = delta_temp
        print(f'{half_delta_f_delta_Db = }')
        print(f'{delta_temp = }')
        data_dict[name]['std_noise_zoom'] = std_noise_zoom

        response_at_res_freq_s21 = hmm.get_s21_at_res_freq(freq, s21, design_freq_Hz)
        data_dict[name]['response_at_res_freq_s21'] = response_at_res_freq_s21
        response_at_res_freq_s21_list.append(response_at_res_freq_s21)
        overcoupled = hmm.is_overcoupled(s11, s11_ph, name, print_plot=False)
        data_dict[name]['overcoupled'] = overcoupled
        overcoupled_list.append(overcoupled)
        beta = hmm.get_beta(s11, overcoupled=overcoupled)
        data_dict[name]['beta'] = beta
        Q_0 = hmm.get_Q_0(Q_L, beta)
        data_dict[name]['Q_0'] = Q_0

        delta_freq = res_freq - design_freq_Hz
        z_coord = hmm.get_zcoord_from_freq_delta(delta_freq)
        data_dict[name]['z_coord'] = z_coord
        z_coord_list.append(z_coord)

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
    print(f'{Q_0_list = }')
    print(f'{Q_L_list = }')
    print(f'{beta_list = }')
    print(f'{res_freq_list = }')
    print(f'{min_gamma_s11_list = }')
    print(f'{res_freq_s21_list = }')
    print(f'{max_gamma_s21_list = }')
    print(f'{max_gamma_raw_s21_list = }')
    print(f'{min_feq_minus3dB_s21_list = }')
    print(f'{max_feq_minus3dB_s21_list = }')
    print(f'{bandwidth_list = }')
    print(f'{response_at_res_freq_s21_list = }')
    print(f'{overcoupled_list = }')
    print(f'{z_coord_list = }')


else:

    # read in data pickle

    with open(f'{pkl_addr}\\data_dict.pkl', 'rb') as handle:
        data_dict = pickle.load(handle)

    # for key, item in data_dict.items():
    #     input(f'{key}:    {item}')

    name_keys = ['ins00_no_cathode_50degs', 'ins01a_cat8_50deg', 'ins01b_cat8_50p1deg', 'ins01_cath8_50degs',
                 'ins02_cath8_49p8degs', 'ins03a_cath8_50degs', 'ins03_cath8_50degs', 'ins04_cath8_50degs',
                 'ins05_cath8_50degs', 'ins06_cath8_50p2degs']
    Q_0_list = [11837.396995158368, 11461.022995475754, 11543.083181090495, 11878.054235832227, 11881.753530339533,
                11463.302293046796, 11375.561481992734, 11300.980722974575, 11849.197933937508, 13386.803260686087]
    Q_L_list = [5833.793774319066, 5636.345864661655, 5679.25, 5856.62109375, 5856.5078125, 5636.3646616541355,
                5594.451492537313, 5553.051851851852, 5833.708171206225, 6579.916666666667]
    beta_list = [1.029107893266257, 1.033413717091641, 1.0325013304733013, 1.028141149255517, 1.0288120345335119,
                 1.03381132719021, 1.0333649325884944, 1.0350936790201017, 1.0311605562345898, 1.034494346790525]
    res_freq_list = [2998580000.0, 2998548000.0, 2998644000.0, 2998594000.0, 2998538000.0, 2998546000.0, 2998628000.0,
                     2998652000.0, 2998526000.0, 3000444000.0]
    min_gamma_s11_list = [-20.932287, -18.52713, -19.131866, -21.665035, -21.38899, -18.132273, -17.761797, -18.523979,
                          -20.99095, -18.099392]
    res_freq_s21_list = [2998570000.0, 2998536000.0, 2998644000.0, 2998590000.0, 2998532000.0, 2998546000.0,
                         2998626000.0, 2998648000.0, 2998526000.0, 3000442000.0]
    max_gamma_s21_list = [-70.04373457000001, -70.34012320000001, -70.23740544, -69.9672446, -69.99944525, -70.33473394,
                          -70.39605728, -70.32988387, -70.05706568, -69.77073468]
    max_gamma_raw_s21_list = [-69.536636, -69.681618, -69.590088, -69.368622, -69.1455, -69.740921, -69.912674,
                              -69.823639, -69.488625, -69.008713]
    min_feq_minus3dB_s21_list = [2998316000.0, 2998272000.0, 2998372000.0, 2998328000.0, 2998274000.0, 2998272000.0,
                                 2998354000.0, 2998374000.0, 2998262000.0, 3000210000.0]
    max_feq_minus3dB_s21_list = [2998830000.0, 2998804000.0, 2998900000.0, 2998840000.0, 2998786000.0, 2998804000.0,
                                 2998890000.0, 2998914000.0, 2998776000.0, 3000666000.0]
    bandwidth_list = [514000.0, 532000.0, 528000.0, 512000.0, 512000.0, 532000.0, 536000.0, 540000.0, 514000.0,
                      456000.0]
    response_at_res_freq_s21_list = [-70.30513, -70.214592, -71.014854, -70.394852, -69.79213, -70.575706, -71.352661,
                                     -71.077904, -69.860779, -86.739891]
    overcoupled_list = [None, None, None, None, None, None, None, None, None, None]
    z_coord_list = [-55.67153792623521, -33.40292275574113, -100.20876826722338, -65.41405706332637,
                    -26.443980514961726, -32.011134307585245, -89.07446068197633, -105.7759220598469,
                    -18.093249826026444, -1352.8183716075157]


temerature_list = [50., 50., 50., 50.1, 49.8, 50., 50., 50., 50., 50.2]
temp_50_freq = 3000444000.0
temp_52_freq = 3000372000.0
delta_freq_per_degree_Hz = (temp_52_freq - temp_50_freq) / 2.
delta_freq_per_degree_kHz = delta_freq_per_degree_Hz / 1000.
print(f'{delta_freq_per_degree_kHz = } kHz')

temp_delta_list = [i-design_temp for i in temerature_list]
# [f(x) if x is not None else '' for x in xs]
print(f'{temp_delta_list = }')
temp_corrected_freq_list = [res_freq_list[i] for i in range(len(res_freq_list)) ]

# plt.scatter(res_freq_list[:-1], max_gamma_s21_list[:-1])
# plt.scatter(temp_corrected_freq_list[:-1], max_gamma_s21_list[:-1], marker='x', s=80)

freqs_centred_kHz = [(i-design_freq_Hz)/1.e3 for i in res_freq_list[:-1]]
freqs_GHz_for_plots = [i/1.e9 for i in res_freq_list[:-1]]
# plt.scatter(freqs_centred_kHz, max_gamma_raw_s21_list[:-1], color = 'r')
plt.scatter(freqs_centred_kHz, max_gamma_s21_list[:-1], color = 'r')
# plt.scatter(freqs_GHz_for_plots, max_gamma_raw_s21_list[:-1], marker='x', s=80)
plt.xlabel(r'$\Delta$''f [kHz]')
plt.ylabel('$\Gamma_{MAX}$''[dB]')
plt.savefig(f'{savepath}\\freq_centred_vs_max_s21_smoothed.png')
plt.close('all')


# histogram of the noise at the peak s21 region
# nbins = 9
# # plt.hist(STD_noise, bins=int(np.sqrt(len(STD_noise))), histtype='step', lw=0.8, color='k')
# n, bedges, patches = plt.hist(STD_noise, bins=nbins, histtype='step', lw=1., color='k')
# plt.vlines(mean_STD_noise, 0, max(n), ls='-', lw=0.8, color='r', label=f'{mean_STD_noise = :1.3f}')
# plt.ylabel('N')
# plt.xlabel(r'$\sigma$'' [Db]')
# plt.legend(loc='upper right')
# plt.savefig(f'{savepath}\\STD_noise_hist.png')
# plt.close('all')

# beta & Q_0 histograms

nbins = 4
# res_freq_list = res_freq_list
# beta_list = beta_list
# Q_0_list = Q_0_list

res_freq_list_centred = [(i-design_freq_Hz)/1.e3 for i in res_freq_list]
n, bedges, patches = plt.hist(res_freq_list_centred, bins = nbins, histtype='step', lw=0.8, color='r')
plt.vlines(np.mean(res_freq_list_centred), 0., max(n), ls='--', lw=0.8, color='r')
plt.text(min(res_freq_list_centred), max(n), f'mean = {np.mean(res_freq_list_centred):1.3f}\n'
                                             f'SDev = {np.std(res_freq_list_centred):1.3f})')
plt.xlabel(r'$\Delta$''f [kHz]')
plt.ylabel('N')

plt.savefig(f'{savepath}\\freq_histpng')
plt.close('all')

n, bedges, patches = plt.hist(beta_list, bins = nbins, histtype='step', lw=0.8, color='r')
plt.vlines(np.mean(beta_list), 0., max(n), ls='--', lw=0.8, color='r')
plt.text(min(bedges), max(n), f'mean = {np.mean(beta_list):1.3f}\n'
                                             f'SDev = {np.std(beta_list):1.3f})')
plt.xlabel(r'$\beta$')
plt.ylabel('N')

plt.savefig(f'{savepath}\\Beta_histpng')
plt.close('all')

n, bins, patches = plt.hist(Q_0_list, bins = nbins, histtype='step', lw=0.8, color='r')
plt.vlines(np.mean(Q_0_list), 0., max(n), ls='--', lw=0.8, color='r')
plt.text(min(bins), max(n), f'mean = {np.mean(Q_0_list):1.3f}\n'
                                             f'SDev = {np.std(Q_0_list):1.3f})')
plt.xlabel(r'$Q_0$')
plt.ylabel('N')
plt.savefig(f'{savepath}\\Q_0_hist.png')
plt.close('all')

n, bins, patches = plt.hist(max_gamma_s21_list, bins = nbins, histtype='step', lw=0.8, color='r')
plt.vlines(np.mean(max_gamma_s21_list), 0., max(n), ls='--', lw=0.8, color='r')
plt.text(min(bins), max(n), f'mean = {np.mean(max_gamma_s21_list):1.3f}\n'
                                             f'SDev = {np.std(max_gamma_s21_list):1.3f})')
plt.xlabel(r'$Q_0$')
plt.ylabel('N')
plt.savefig(f'{savepath}\\max_gamma_s21.png')
plt.close('all')

print(f'Q_0 mean = {np.mean(Q_0_list):1.3f}\n'
      f'Q_0 SDev = {np.std(Q_0_list):1.3f}')
print(f'B mean = {np.mean(beta_list):1.3f}\n'
      f'B SDev = {np.std(beta_list):1.3f}')
print(f'max gamma mean = {np.mean(max_gamma_s21_list):1.3f}\n'
      f'max gamma SDev = {np.std(max_gamma_s21_list):1.3f}')
print(f'freq mean = {np.mean(res_freq_list):1.3f}\n'
      f'freq SDev = {np.std(res_freq_list):1.3f}')

