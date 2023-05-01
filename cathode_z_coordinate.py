print('imported czc')
import csv, sys, os
import numpy as np
import matplotlib.pyplot as plt



import PhD_Master_Module as pmm

# r_0 = (2.4048*2.998e8)/(2.*np.pi*2.9985e9)
# print(f'{r_0 = }')

'''
C1Radius -> 4.0765,
BlendRadius -> 0.4,
IrisRadius -> 2.6 - 1.2,
MajorRadius -> 1.2,
MinorRadius -> 0.8,
C2Radius -> 4.0865
'''

radius = [42.0,
          42.5,
          43.0,
          43.5,
          44.0]

freq_GHz = [3.1270790206506,
            3.0853956613406,
            3.0446676125191,
            3.0048592213129,
            2.9659554987067]

c, m = pmm.best_fit(radius, freq_GHz)
r_pi = (2.9985-c)/m
print(f'{r_pi = }')
# exit()

save_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\CST_cathode_coordinate_simulation\analysis'
data_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\CST_cathode_coordinate_simulation\data'
fnames = os.listdir(data_addr)
print(fnames)
frq_re_im_fnames = [i for i in fnames if 'um.csv' in i]
freqs = []
for f in frq_re_im_fnames:

    freq, re, im = pmm.read_csv_three_columns_delimiter(f'{data_addr}\\{f}', delimiter=',', fname='unknown')
    abs = [np.sqrt(re[i]**2. + im[i]**2.) for i in range(len(freq))]
    min_val, min_idx = pmm.get_min_val_idx_from_list(abs)
    #print(f'{min_idx = }')
    res_freq = freq[min_idx]
    freqs.append(res_freq)
    if f == '0um.csv':
        design_freq_GHz = res_freq
    plt.plot(freq, re)
    plt.plot(freq, im)
    plt.plot(freq, abs)
    plt.savefig(f'{save_addr}\\{f}reimabs.png')
    plt.close('all')

z_coords = [0., 100., 200., -100., -200.]
c, m = pmm.best_fit(z_coords, freqs)
print(f'Gradient = {m} GHz/um')
m_Hz_um = m*1e9
m_Hz_m = m_Hz_um * 1e6
m_Hz_mm = m_Hz_um*1e3
m_kHz_mm = m_Hz_um/1e3
m_kHz_um = m_Hz_um/1e3
m_um_kHz = 1./m_kHz_um





for idx, i in enumerate(z_coords):
    print(f'{i}     {freqs[idx]}')

# assume -50kHz per deg C
# freq_temp_gradient_Hz_degC = -5e4
# temp_freq_gradient_Hz_degC = 1./freq_temp_gradient_Hz_degC
# freq_temp_gradient_kHz_degC = freq_temp_gradient_Hz_degC/1e3
# temp_freq_gradient_degC_kHz = 1./freq_temp_gradient_kHz_degC
# print(f'{freq_temp_gradient_kHz_degC = }\n{temp_freq_gradient_degC_kHz = }')
#
# um_degC = m_kHz_um * temp_freq_gradient_degC_kHz
# print(f'{um_degC = }')
#
# m_degC = m_Hz_m * temp_freq_gradient_Hz_degC
# print(f'{m_degC = }')


plot_freqs_kHz = [(i-design_freq_GHz)*1e6 for i in freqs]
delta_f_plot = max(plot_freqs_kHz) - min(plot_freqs_kHz)
x_fit = np.linspace(min(z_coords), max(z_coords), 1000)
c_plot, m_plot = pmm.best_fit(z_coords, plot_freqs_kHz)
print(f'gradient = {m_Hz_mm:1.0f} Hz/mm')
print(f'gradient = {m_Hz_um:1.0f} Hz/um')
print(f'c_plot = {c_plot:1.4f} Hz')
y_fit = [c_plot+m_plot*x for x in x_fit]

temp_vector = [i*(1./-50.) for i in plot_freqs_kHz]
x_fit_degC_um = np.linspace(min(temp_vector), max(temp_vector), 1000)
c_z_temp, m_z_temp = pmm.best_fit(temp_vector, z_coords)
y_fit_degC_um = [c_z_temp+x*m_z_temp for x in x_fit_degC_um]

design_temp = 50.
C20_temp = 50.9
C20_temp_delta = design_temp - C20_temp
C20_zcoord_um = m_z_temp * C20_temp_delta
C20_freq_kHz_delta = -50.*C20_temp_delta
C20_freq_GHz = (C20_freq_kHz_delta/1e3)+design_freq_GHz
print(f'{design_freq_GHz = }')
print(f'{C20_freq_GHz = }')
print(f'{C20_zcoord_um = }')

C16_temp = 49.1
C16_temp_delta = design_temp - C16_temp
C16_zcoord_um = m_z_temp * C16_temp_delta
C16_freq_kHz_delta = -50.*C16_temp_delta
C16_freq_GHz = (C16_freq_kHz_delta/1e6)+design_freq_GHz
print(f'{design_freq_GHz = }')
print(f'{C16_freq_GHz = }')
print(f'{C16_zcoord_um = }')

plt.scatter(z_coords, plot_freqs_kHz, marker='o', s=15, color='k', label='sims')
plt.scatter(C16_zcoord_um, C16_freq_kHz_delta, marker='o', s=20, color='b', label='C16')
plt.scatter(C20_zcoord_um, C20_freq_kHz_delta, marker='o', s=20, color='g', label='C20')
plt.plot(x_fit, y_fit, ls='--', lw=0.6, color='r', label='fit')
plt.text(-180, min(plot_freqs_kHz) + (delta_f_plot * 0.1), f'design frequency = {design_freq_GHz:1.4f} GHz\ngradient = {m_kHz_mm:1.4f} kHz/'r'$\mu$m')
# plt.text(-180, min(freqs) + (delta_f * 0.1), f'gradient = {m_Hz_mm:1.0f} Hz/mm')
plt.xlabel('Z-coordinate ('r'$\mu$''m)')
plt.ylabel(r'$\Delta$''f (kHz)')
plt.legend(loc='upper right')
plt.savefig(f'{save_addr}\\zcoord_freqs.png')
plt.close('all')

print(f'{m_z_temp = }')
plt.scatter(temp_vector, z_coords, marker='o', s=15, color='k', label='sims')
plt.scatter(C16_temp_delta, C16_zcoord_um, marker='o', s=20, color='b', label='C16')
plt.scatter(C20_temp_delta, C20_zcoord_um, marker='o', s=20, color='g', label='C20')
plt.plot(x_fit_degC_um, y_fit_degC_um, ls='--', lw=0.6, color='r', label='fit')
plt.text(0., min(z_coords) , f'gradient = {m_z_temp:1.4f} 'r'$\mu$''m / 'r'$^{\circ}$''C')
# plt.text(-180, min(freqs) + (delta_f * 0.1), f'gradient = {m_Hz_mm:1.0f} Hz/mm')
plt.xlabel(r'$\Delta$''T ('r'$^{\circ}$''C')
plt.ylabel(r'$\Delta$''z ('r'$\mu$''m)')
plt.legend(loc='upper left')
plt.savefig(f'{save_addr}\\temp_zcoord.png')
plt.close('all')

c_degC_kHz, m_degC_kHz = pmm.best_fit(temp_vector, plot_freqs_kHz)
x_fit_degC_kHz = x_fit_degC_um
y_fit_degC_kHz = [c_degC_kHz+x*m_degC_kHz for x in x_fit_degC_kHz]
plt.scatter(temp_vector, plot_freqs_kHz, marker='o', s=15, color='k', label='sims')
plt.scatter(C16_temp_delta, C16_freq_kHz_delta, marker='o', s=20, color='b', label='C16')
plt.scatter(C20_temp_delta, C20_freq_kHz_delta, marker='o', s=20, color='g', label='C20')
plt.plot(x_fit_degC_kHz, y_fit_degC_kHz, ls='--', lw=0.6, color='r', label='fit')
plt.text(-5.5, min(z_coords) , f'gradient = {m_degC_kHz:1.4f} kHz / 'r'$^{\circ}$''C')
# plt.text(-180, min(freqs) + (delta_f * 0.1), f'gradient = {m_Hz_mm:1.0f} Hz/mm')
plt.xlabel(r'$\Delta$''T ('r'$^{\circ}$''C')
plt.ylabel(r'$\Delta$''f (kHz)')
plt.legend(loc='upper right')
plt.savefig(f'{save_addr}\\temp_freq.png')
plt.close('all')
