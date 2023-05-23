import sys, os, csv
import numpy as np
from matplotlib import pyplot as plt
from decimal import Decimal
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

print(f'{HRRG_Field_Flatness = }')
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

plt.plot(ff_x_fit_hrrg, y_fit_poly_hrrg(ff_x_fit_hrrg), ls='-', lw=0.8, color='r', label='fit')
plt.scatter(zcoords, HRRG_Field_Flatness, marker='o', s=15, color='k', label='sim')
plt.text(-300, 101, f'FFfit = {A_poly_hrrg:.2e}'r'$z^2$'f'x + {B_poly_hrrg:.2e}z + {C_poly_hrrg:.2e}')
plt.text(-300, 100, f'SLF grad = {ff_m_hrrg:.4f}')
plt.legend()
plt.xlabel(f'z-coordinate ('r'$\mu$''m)')
plt.ylabel(f'Field Flatness (%)')
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

ff_c_lrrg, ff_m_lrrg = pmm.best_fit(lrrg_zcoords, LRRG_Field_Flatness)
ff_x_fit_lrrg = np.linspace(min(lrrg_zcoords), max(lrrg_zcoords), 1000)



y_poly_lrrg = np.polyfit(lrrg_zcoords, LRRG_Field_Flatness, 2)
print(f'{y_poly_lrrg = }')
A_poly_lrrg = float(y_poly_lrrg[0])
B_poly_lrrg = float(y_poly_lrrg[1])
C_poly_lrrg = float(y_poly_lrrg[2])

y_fit_poly_lrrg = np.poly1d(y_poly_lrrg)

plt.plot(ff_x_fit_lrrg, y_fit_poly_lrrg(ff_x_fit_lrrg), ls='-', lw=0.8, color='r', label='fit')
plt.scatter(lrrg_zcoords, LRRG_Field_Flatness, marker='o', s=15, color='k', label='sim')
plt.text(-195, 105, f'FFfit = {A_poly_lrrg:.2e}'r'$z^2$'f'x + {B_poly_lrrg:.2e}z + {C_poly_lrrg:.2e}')
plt.text(-195, 101, f'SLF grad = {ff_m_lrrg:.4f}')
plt.legend()
plt.xlabel(f'z-coordinate ('r'$\mu$''m)')
plt.ylabel(f'Field Flatness (%)')
plt.savefig(f'{save_addr}\\LRRG_z_FF.png')
plt.close('all')

print(f'\nLRRG z vs FF\n')
for idx, z in enumerate(lrrg_zcoords):
    print(f'{z}, {LRRG_Field_Flatness[idx]}')

# combined HRRG LRRG plot

plt.plot(ff_x_fit_hrrg, y_fit_poly_hrrg(ff_x_fit_hrrg), ls='-', lw=0.8, color='r', label='HRRG fit')
plt.scatter(zcoords, HRRG_Field_Flatness, marker='o', s=15, color='darkorange', label='HRRG sim')
# plt.text(-300, 101, f'FFfit = {A_poly_hrrg:.2e}'r'$z^2$'f'x + {B_poly_hrrg:.2e}z + {C_poly_hrrg:.2e}')
# plt.text(-300, 100, f'SLF grad = {ff_m_hrrg:.4f}')
plt.plot(ff_x_fit_lrrg, y_fit_poly_lrrg(ff_x_fit_lrrg), ls='-', lw=0.8, color='g', label='LRRG fit')
plt.scatter(lrrg_zcoords, LRRG_Field_Flatness, marker='o', s=15, color='b', label='LRRG sim')
# plt.text(-195, 105, f'FF = {A_poly_lrrg:.2e}'r'$z^2$'f'x + {B_poly_lrrg:.2e}z + {C_poly_lrrg:.2e}')
plt.text(200, 111, f'LRRG')
plt.text(200, 95, f'HRRG')
plt.xlabel(f'z-coordinate ('r'$\mu$''m)')
plt.ylabel(f'Field Flatness (%)')
plt.legend()
plt.savefig(f'{save_addr}\\HRRG_LRRG_z_FF.png')
plt.close('all')