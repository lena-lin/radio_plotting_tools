import numpy as np
import os
import click

from ..coordinates import get_pixel_coordinates
from ..gaussian_components import get_gaussian_components

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from datetime import datetime
from astropy.io import fits
import astropy.units as u
import glob
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit


def get_obs_dates(list_files):
    dates = sorted([fits.open(file)['PRIMARY'].header['DATE-OBS'] for file in list_files])
    return dates



@click.command()
@click.argument('input_path', type=click.Path(file_okay=False, dir_okay=True))
@click.argument('file_name', type=click.STRING)
def main(
    input_path,
    file_name,
):
    files = sorted(glob.glob('{}/*/{}.fits'.format(input_path, file_name)))
    dates = get_obs_dates(files)

    fig = plt.figure(figsize=(12,20))
    ax = fig.add_subplot(111, aspect='equal')

    for num_epoch, f in enumerate(files):
        epoch = fits.open(f)
        x,y = get_pixel_coordinates(epoch['PRIMARY'].header)
        clean_map = epoch['PRIMARY'].data[0][0]
        df_components = get_gaussian_components(epoch['AIPS CC'].data)

        time_diff = (([(np.datetime64(i) - np.datetime64(dates[0])) for i in dates]) / np.timedelta64(1, 'D') / 60) * 2
        y_pos_components = df_components.y_positions.values - time_diff[num_epoch]
        y = y - time_diff[num_epoch]

        ax.contour(x, y, clean_map, norm=LogNorm())#, levels=np.logspace(np.log10(loglevs_min), np.log10(flux_max), 10))
        ax.plot(df_components.x_positions, y_pos_components, color='black',  marker='.', markersize= 2, linestyle='none')
    plt.show()












if __name__ == "__main__":
    main()









# datapaths = ([])
# for data in glob.glob('../../../../radiokram/modeling/1638+118/*/final_mod.fits'):
#     datapaths = np.append(datapaths, data)
#
# # Array With All Dates
# dates = sorted([fits.open(file)['PRIMARY'].header['DATE-OBS'] for file in datapaths])
#
# def time_difference(dates, shortest_difference, scaling):
#     time_diff = ( ([(np.datetime64(i) - np.datetime64(dates[0])) for i in dates]) / np.timedelta64(1, 'D') / shortest_difference ) * scaling
#     return time_diff
#
# fig = plt.figure(figsize=(12,20))
# ax = fig.add_subplot(111, aspect='equal')
#
# # Define Empty Arrays
# num_epoch = 0
# y_ticks = ([])
# x_core = ([])
# y_core = ([])
# x_c1 = ([])
# y_c1 = ([])
#
# # Define Linear Function
# def linear_fit(f, m, b):
#     return m * f + b
#
# # Start Iteration Through Different Epochs
#
# for epoch in sorted(datapaths):
#
#     # Loading Data
#     input_file = epoch
#     difmap_data = fits.open(input_file)
#
#     x_n_pixel = difmap_data['PRIMARY'].header['NAXIS1']
#     x_ref_pixel = difmap_data['PRIMARY'].header['CRPIX1']
#     x_inc = (difmap_data['PRIMARY'].header['CDELT1'] * u.degree).to(u.mas)
#     x_ref_value = (difmap_data['PRIMARY'].header['CRVAL1'] * u.degree).to(u.mas)
#     y_n_pixel = difmap_data['PRIMARY'].header['NAXIS2']
#     y_ref_pixel = difmap_data['PRIMARY'].header['CRPIX2']
#     y_inc = (difmap_data['PRIMARY'].header['CDELT2'] * u.degree).to(u.mas)
#     y_ref_value = (difmap_data['PRIMARY'].header['CRVAL2'] * u.degree).to(u.mas)
#
#     x = np.linspace(x_ref_pixel * x_inc, -x_ref_pixel * x_inc, x_n_pixel)
#     y = np.linspace(-y_ref_pixel * y_inc, y_ref_pixel * y_inc, y_n_pixel)
#
#     clean_map = difmap_data['PRIMARY'].data[0][0]
#     flux_min = clean_map.min()
#     flux_max = clean_map.max()
#     loglevs_min = flux_min + (flux_max - flux_min) * 0.007
#
#     delta_x = (difmap_data['AIPS CC'].data['DELTAX'] * u.degree).to(u.mas)
#     delta_y = (difmap_data['AIPS CC'].data['DELTAY'] * u.degree).to(u.mas)
#     major_ax = (difmap_data['AIPS CC'].data['MAJOR AX'] * u.degree).to(u.mas)
#     minor_ax = (difmap_data['AIPS CC'].data['MINOR AX'] * u.degree).to(u.mas)
#     phi_ellipse = (difmap_data['AIPS CC'].data['POSANGLE'])
#
#     # Shift Core to (0,0)
#     delta_x = delta_x - delta_x[0]
#     delta_y = delta_y - delta_y[0]
#
#     x = - (x - delta_x[0])
#     y = y - delta_y[0]
#
#     # Calculate Timedifferences
#     time_diff = time_difference(dates, 60, 2)
#
#     # y data has to be shifted for the different epochs
#     delta_y = delta_y.value - time_diff[num_epoch]
#     y = y.value - time_diff[num_epoch]
#
#     # Plot Data
#     plt.plot(delta_x, delta_y, color='black',  marker='.', markersize= 2, linestyle='none')
#
#     ax.set_aspect(1)
#     ax.contour(x, y, clean_map, norm=LogNorm(), levels=np.logspace(np.log10(loglevs_min), np.log10(flux_max), 10))
#
#     # Plot Ellipses
#     for i in range(len(delta_x)):
#         ellipse = Ellipse((delta_x[i].value, delta_y[i]), height=major_ax[i].value, width=minor_ax[i].value, angle=phi_ellipse[i], edgecolor='black', facecolor='none', linewidth=0.5)
#         ax.add_artist(ellipse)
#
#     # Generate Arrays For Same Corecomponents
#     x_core = np.append(x_core, delta_x[0].value)
#     y_core = np.append(y_core, delta_y[0])
#     x_c1 = np.append(x_c1, delta_x[1].value)
#     y_c1 = np.append(y_c1, delta_y[1])
#
#
#     # Generate Y_Ticks
#     y_ticks = np.append(y_ticks ,delta_y[0])
#     num_epoch += 1
#
# x = np.linspace(-40, 20, 10)
# # Line Core
# plt.plot(x_core, y_core, color='red', marker='+', label='Core', linestyle='none', markersize=5)
# plt.axvline(0, ymin=0, ymax=1, color='red', linewidth=1)
#
# # Linear Fit C1
# print(x_c1)
# print(y_c1)
# params, covariance = curve_fit(linear_fit, x_c1, y_c1, p0=[-100000, 10000000])
# errors = np.sqrt(np.diag(covariance))
# plt.plot(x_c1, y_c1, color='blue', marker='+', label='C1', linestyle='none', markersize=5)
# plt.plot(x, linear_fit(x,*params), color='blue', linewidth=1)
# print(params[0], params[1])
#
# plt.yticks(y_ticks, dates)
# plt.xlabel('Relative RA / mas')
# plt.legend()
#
# plt.xlim(5, -15)
# plt.ylim(-35, 3)
# plt.tight_layout()
# #plt.show()
# plt.savefig('dynamic_plot.pdf', bbox_inches='tight', pad_inches=0.1)
