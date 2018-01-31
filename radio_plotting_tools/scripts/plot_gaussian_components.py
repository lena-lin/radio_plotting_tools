import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from datetime import datetime
from astropy.io import fits
import astropy.units as u

# Loading Data
input_file = "../../../modeling/1638+118/2013-12-15/final_mod.fits"
difmap_data = fits.open(input_file)

delta_x = (difmap_data['AIPS CC'].data['DELTAX'] * u.degree).to(u.mas)
delta_y = (difmap_data['AIPS CC'].data['DELTAY'] * u.degree).to(u.mas)
major_ax = (difmap_data['AIPS CC'].data['MAJOR AX'] * u.degree).to(u.mas)
minor_ax = (difmap_data['AIPS CC'].data['MINOR AX'] * u.degree).to(u.mas)
phi_ellipse = (difmap_data['AIPS CC'].data['POSANGLE'])

# Shift Core to (0,0)
delta_x = delta_x - delta_x[0]
delta_y = delta_y - delta_y[0]


# Array With Dates
dates = difmap_data['PRIMARY'].header['DATE-OBS']

# Calculate Timedifferences
def time_difference(dates, shortest_difference, scaling):
    time_diff = ( ([(np.datetime64(i) - np.datetime64(dates[0])) for i in dates]) / np.timedelta64(1, 'D') / shortest_difference ) * scaling
    return time_diff

# y data has to be shifted for the different epochs
#y = y - time_diff[i]

# Plot Data
plt.figure(figsize=(16,12))
fig, ax = plt.subplots()
plt.plot(delta_x, delta_y,  marker='.', markersize= 3, linestyle='none', label='Gaussian Components')

# Plot Ellipses
for i in range(len(delta_x)):
    ellipse = Ellipse((delta_x[i].value, delta_y[i].value), height=major_ax[i].value, width=minor_ax[i].value, angle=phi_ellipse[i], edgecolor='black', facecolor='none')
    ax.axis('equal')
    ax.add_artist(ellipse)

# Generate Y_Ticks
y_ticks = delta_y[i].value
plt.yticks([y_ticks], [str(dates)])

plt.xlim(5, -15)
plt.ylim(-10, 10)
plt.xlabel('Relative RA / mas')
plt.legend()
plt.tight_layout()
plt.savefig('dynamic_plot.pdf', bbox_inches='tight', pad_inches=0.1)
