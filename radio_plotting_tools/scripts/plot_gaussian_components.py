import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from datetime import datetime
from astropy.io import fits
import astropy.units as u
import glob
from matplotlib.colors import LogNorm

datapaths = ([])
for data in glob.glob('../../../../radiokram/modeling/1638+118/*/final_mod.fits'):
    datapaths = np.append(datapaths, data)

# Array With All Dates
dates = sorted([fits.open(file)['PRIMARY'].header['DATE-OBS'] for file in datapaths])

def time_difference(dates, shortest_difference, scaling):
    time_diff = ( ([(np.datetime64(i) - np.datetime64(dates[0])) for i in dates]) / np.timedelta64(1, 'D') / shortest_difference ) * scaling
    return time_diff

plt.figure(figsize=(12,8))
fig, ax = plt.subplots()
num_epoch = 0
y_ticks = ([])

for epoch in sorted(datapaths):

    # Loading Data
    input_file = epoch
    difmap_data = fits.open(input_file)

    x_n_pixel = difmap_data['PRIMARY'].header['NAXIS1']
    x_ref_pixel = difmap_data['PRIMARY'].header['CRPIX1']
    x_inc = (difmap_data['PRIMARY'].header['CDELT1'] * u.degree).to(u.mas)
    x_ref_value = (difmap_data['PRIMARY'].header['CRVAL1'] * u.degree).to(u.mas)
    y_n_pixel = difmap_data['PRIMARY'].header['NAXIS2']
    y_ref_pixel = difmap_data['PRIMARY'].header['CRPIX2']
    y_inc = (difmap_data['PRIMARY'].header['CDELT2'] * u.degree).to(u.mas)
    y_ref_value = (difmap_data['PRIMARY'].header['CRVAL2'] * u.degree).to(u.mas)

    x = np.linspace(x_ref_pixel * x_inc, -x_ref_pixel * x_inc, x_n_pixel)
    y = np.linspace(-y_ref_pixel * y_inc, y_ref_pixel * y_inc, y_n_pixel)

    clean_map = difmap_data['PRIMARY'].data[0][0]
    flux_max = clean_map.max()

    delta_x = (difmap_data['AIPS CC'].data['DELTAX'] * u.degree).to(u.mas)
    delta_y = (difmap_data['AIPS CC'].data['DELTAY'] * u.degree).to(u.mas)
    major_ax = (difmap_data['AIPS CC'].data['MAJOR AX'] * u.degree).to(u.mas)
    minor_ax = (difmap_data['AIPS CC'].data['MINOR AX'] * u.degree).to(u.mas)
    phi_ellipse = (difmap_data['AIPS CC'].data['POSANGLE'])

    # Shift Core to (0,0)
    delta_x = delta_x - delta_x[0]
    delta_y = delta_y - delta_y[0]

    x = - (x - delta_x[0])
    y = y - delta_y[0]

    # Calculate Timedifferences
    time_diff = time_difference(dates, 60, 1.5)

    # y data has to be shifted for the different epochs
    delta_y = delta_y.value - time_diff[num_epoch]
    y = y.value - time_diff[num_epoch]

    # Plot Data
    plt.plot(delta_x, delta_y,  marker='.', markersize= 3, linestyle='none', label='Gaussian Components')

    ax.set_aspect(1)
    ax.contour(x, y, clean_map, norm=LogNorm(), levels=np.logspace(np.log10(flux_max*1e-2), np.log10(flux_max), 10))

    # Plot Ellipses
    for i in range(len(delta_x)):
        ellipse = Ellipse((delta_x[i].value, delta_y[i]), height=major_ax[i].value, width=minor_ax[i].value, angle=phi_ellipse[i], edgecolor='black', facecolor='none')
        ax.axis('equal')
        ax.add_artist(ellipse)

    # Generate Y_Ticks
    y_ticks = np.append(y_ticks ,delta_y[0])
    num_epoch += 1

plt.yticks(y_ticks, dates)
plt.xlabel('Relative RA / mas')
#plt.legend()

plt.xlim(5, -15)
plt.ylim(-30, 5)
plt.tight_layout()
plt.show()

#plt.savefig('dynamic_plot.pdf', bbox_inches='tight', pad_inches=0.1)
