import sys, os, csv
import numpy as np
from matplotlib import pyplot as plt
import PhD_Master_Module as pmm

save_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\CST_cathode_coordinate_simulation\analysis'
data_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\CST_cathode_coordinate_simulation\data'
fnames = os.listdir(data_addr)

ff_fnames = [i for i in fnames if 'FieldFlatness.csv' in i]
print(ff_fnames)

def field_flatness(z, Ez, Ncells=1.5, flip=False, plot=False):
    '''
    given a z vs Ez field profile and the number of cells, it returns the field flatness as a percentage
    :param z:
    :param Ez:
    :param Ncells: either 1.5 or 2.5
    :return: %FF
    '''
    if flip:
        Ez = Ez[::-1]

    if Ncells == 1.5:
        iris_1_index = 1680

        half_cell_max = max(Ez[:iris_1_index])
        full_cell_max = max(Ez[iris_1_index:])

        # global_max = max([half_cell_max, full_cell_max])
        # global_2nd_max = min([half_cell_max, full_cell_max])

        FF = full_cell_max / half_cell_max* 100.

    elif Ncells == 2.5:

        iris_1_index = 1150
        iris_2_index = 3450

        half_cell_max = max(Ez[:iris_1_index])
        full_1_cell_max = max(Ez[iris_1_index:iris_2_index])
        full_2_cell_max = max(Ez[iris_2_index:])
        min_full_cells = min([full_1_cell_max, full_2_cell_max])

        if plot:
            plt.plot(z, Ez)
            plt.hlines(half_cell_max, 0., max(z), ls='--', lw=0.8, color='r')
            plt.hlines(min_full_cells, 0., max(z), ls='--', lw=0.8, color='r')
            plt.vlines(z[iris_1_index], 0., max(Ez), ls='--', lw=0.8, color='g')
            plt.vlines(z[iris_2_index], 0., max(Ez), ls='--', lw=0.8, color='g')
            plt.show()
        # global_min = min([half_cell_min, full_cell_min])
        # global_2nd_min = min([half_cell_min, full_cell_min])

        FF = min_full_cells / half_cell_max * 100.

    else:
        print(f'\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              f'\nNcells must equal 1.5 or 2.5 only, we have Ncells = {Ncells}'
              f'\naborting analysis.....')
        exit()

    return FF

HRRG_Field_Flatness = []
for f in ff_fnames:
    print(f)
    z, Ez = pmm.read_csv_two_columns_delimiter(f'{data_addr}\\{f}', delimiter='\t')
    '''
    plt.plot(z, Ez)
    plt.show()
    print(f'{len(Ez) = }')
    min_val, min_idx = pmm.get_min_val_idx_from_list(Ez[:5333])
    iris_index = pmm.find_nearest_value_index_in_list(z, z[min_idx])
    print(z[min_idx])
    print(f'{iris_index = }')
    '''
    FF = field_flatness(z, Ez)
    print(f'{FF = }')
    HRRG_Field_Flatness.append(FF)


zcoords = [0., -100., -200., -250., -350., 100., 200., 250., 350.]
ff_c, ff_m = pmm.best_fit(zcoords, HRRG_Field_Flatness)
ff_x_fit = np.linspace(min(zcoords), max(zcoords), 1000)
hrrg_A = 9e-6
hrrg_B = 0.0151
hrrg_C = 96.584
ff_y_fit = [hrrg_A*x**2. * hrrg_B*x + hrrg_C for x in ff_x_fit]
plt.plot(ff_x_fit, ff_y_fit, ls='-', lw=0.8, color='r', label='fit')
plt.scatter(zcoords, HRRG_Field_Flatness, marker='o', s=15, color='k', label='sim')
plt.text(-300, 98, f'FF = 9e-6'r'$z^2$'f'x + 1.15e-2z + {96.584:1.3f}')
plt.legend()
plt.savefig(f'{save_addr}\\HRRG_z_FF.png')
plt.close('all')

print(f'\nHRRG z vs FF\n')
for idx, z in enumerate(zcoords):
    print(f'{z}, {HRRG_Field_Flatness[idx]}')


# LRRG Field Flatness

lrrg_data_addr = r'C:\Users\zup98752\OneDrive - Science and Technology Facilities Council\HRRG\CST_cathode_coordinate_simulation\data\lrrg_data'
fnames = os.listdir(lrrg_data_addr)

lrrg_ff_fnames = [i for i in fnames if 'Amp_Pos_Neg_Pen_' in i]
print(lrrg_ff_fnames)

lrrg_zcoords = [-100., -200., 0., 100., 200.]

LRRG_Field_Flatness = []
for f in lrrg_ff_fnames:
    print(f)
    z, Ez = pmm.read_csv_two_columns_delimiter(f'{lrrg_data_addr}\\{f}', delimiter=',')


    # print(f'{len(Ez) = }')
    # min_val, min_idx = pmm.get_min_val_idx_from_list(Ez[:5333])
    # iris_index = pmm.find_nearest_value_index_in_list(z, z[min_idx])
    # print(z[min_idx])
    # print(f'{iris_index = }')

    FF = field_flatness(z, Ez, flip=True, Ncells=2.5, plot=False)
    print(f'{FF = }')
    LRRG_Field_Flatness.append(FF)

ff_c, ff_m = pmm.best_fit(lrrg_zcoords, LRRG_Field_Flatness)
ff_x_fit = np.linspace(min(lrrg_zcoords), max(lrrg_zcoords), 1000)
ff_y_fit = [ff_c + x*ff_m for x in ff_x_fit]
plt.plot(ff_x_fit, ff_y_fit, ls='-', lw=0.8, color='r', label='fit')
plt.scatter(lrrg_zcoords, LRRG_Field_Flatness, marker='o', s=15, color='k', label='sim')
plt.text(-175, 110, f'FF = {ff_m:1.3}z + {ff_c:1.3f}')
plt.legend()
plt.savefig(f'{save_addr}\\LRRG_z_FF.png')
plt.close('all')

print(f'\nLRRG z vs FF\n')
for idx, z in enumerate(lrrg_zcoords):
    print(f'{z}, {LRRG_Field_Flatness[idx]}')