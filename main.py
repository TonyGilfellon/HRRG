# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import sys, os, csv
pmm_addr = r'C:\Users\zup98752\PycharmProjects\PhD'
sys.path.insert(1,pmm_addr)
# import cathode_z_coordinate as czc


def print_bye(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'\nBye, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #import field_flatness as ff
    #import cathode_z_coordinate
    import VNA_measurements as vm

    print_bye('HRRG')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
