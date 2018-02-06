import numpy as np
import matplotlib as mpl
mpl.use('pgf')
mpl.rcParams.update(
{'font.size': 20,
'font.family': 'serif',
'text.usetex': True,
'pgf.rcfonts': False,
'pgf.texsystem': 'lualatex',
'text.latex.unicode': True,
})

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from datetime import datetime
from astropy.io import fits
import astropy.units as u
import glob
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit

# Array with datapaths for all epochs
datapaths = ([])
for data in glob.glob('../../../../radiokram/modeling/1638+118/*/final_mod.fits'):
    datapaths = sorted(np.append(datapaths, data))

# Array with all sorted dates
dates = sorted([fits.open(file)['PRIMARY'].header['DATE-OBS'] for file in datapaths])

# Function to calculate time differences
def time_difference(dates, shortest_difference, scaling):
    time_diff = ( ([(np.datetime64(i) - np.datetime64(dates[0])) for i in dates]) / np.timedelta64(1, 'D') / shortest_difference ) * scaling
    return time_diff

# Define empty arrays
y_ticks = ([])
delta_x = ([])
delta_y = ([])
major_ax = ([])
minor_ax = ([])
phi_ellipse = ([])
x_n_pixel = ([])
x_ref_pixel = ([])
x_inc = ([])
x_ref_value = ([])
y_n_pixel = ([])
y_ref_pixel = ([])
y_inc = ([])
y_ref_value = ([])
clean_map = ([])
flux_min = ([])
flux_max = ([])
loglevs_min = ([])
c_x = ([])
c_y = ([])

# x_values for linear regression
x_values = np.linspace(-40, 100, 10)

# Colors for different c_i
colors = (['red', 'blue', 'magenta', 'orange', 'green', 'brown', 'cyan', 'darkorange', 'grey', 'olive'])

# Function for linear fit
def linear_fit(f, m, b):
    return m * f + b

# Function to rotate positions of components for linear fit
def rotate(x, y, angle, length):
    ox = np.zeros(length)
    oy = np.zeros(length)
    angle = np.zeros(length) + angle

    qx = ox + np.cos(angle) * (x - ox) - np.sin(angle) * (y - oy)
    qy = oy + np.sin(angle) * (x - ox) + np.cos(angle) * (y - oy)
    return qx, qy

# Function to find nearest points
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

# Load data for different epochs

for epoch in sorted(datapaths):

    input_file = epoch
    difmap_data = fits.open(input_file)

    x_n_pixel.append(difmap_data['PRIMARY'].header['NAXIS1'])
    x_ref_pixel.append(difmap_data['PRIMARY'].header['CRPIX1'])
    x_inc.append((difmap_data['PRIMARY'].header['CDELT1'] * u.degree).to(u.mas))
    x_ref_value.append((difmap_data['PRIMARY'].header['CRVAL1'] * u.degree).to(u.mas))
    y_n_pixel.append(difmap_data['PRIMARY'].header['NAXIS2'])
    y_ref_pixel.append(difmap_data['PRIMARY'].header['CRPIX2'])
    y_inc.append((difmap_data['PRIMARY'].header['CDELT2'] * u.degree).to(u.mas))
    y_ref_value.append((difmap_data['PRIMARY'].header['CRVAL2'] * u.degree).to(u.mas))


    clean_map.append(difmap_data['PRIMARY'].data[0][0])

    delta_x.append(((difmap_data['AIPS CC'].data['DELTAX'] * u.degree).to(u.mas)).value)
    delta_y.append(((difmap_data['AIPS CC'].data['DELTAY'] * u.degree).to(u.mas)).value)
    major_ax.append(((difmap_data['AIPS CC'].data['MAJOR AX'] * u.degree).to(u.mas)).value)
    minor_ax.append(((difmap_data['AIPS CC'].data['MINOR AX'] * u.degree).to(u.mas)).value)
    phi_ellipse.append((difmap_data['AIPS CC'].data['POSANGLE']))

# Shift core to (0,0) and calculate timedifference
time_diff = np.asarray(time_difference(dates, 60, 2))

for i in range(len(dates)):
    delta_x[i] = delta_x[i] - delta_x[i][0]
    delta_y[i] = delta_y[i] - delta_y[i][0] - time_diff[i]

# Make contour plot and plot all components
fig = plt.figure(figsize=(12,20))
ax = fig.add_subplot(111, aspect='equal')

for i in range(len(dates)):
    c_x = np.append(c_x, delta_x[i])
    c_y = np.append(c_y, delta_y[i])

    x = np.linspace(x_ref_pixel[i] * x_inc[i], -x_ref_pixel[i] * x_inc[i], x_n_pixel[i])
    y = np.linspace(-y_ref_pixel[i] * y_inc[i], y_ref_pixel[i] * y_inc[i], y_n_pixel[i]).value - time_diff[i]
    flux_min.append(clean_map[i].min())
    flux_max.append(clean_map[i].max())
    loglevs_min.append(flux_min[i] + (flux_max[i] - flux_min[i]) * 0.007)



    plt.plot(delta_x[i], delta_y[i], marker='.', color='black', linestyle='none', markersize=4)

    ax.set_aspect(1)
    ax.contour(-x, y, clean_map[i], norm=LogNorm(), levels=np.logspace(np.log10(loglevs_min[i]), np.log10(flux_max[i]), 10))

    # Plot ellipses
    for j in range(len(delta_x)):
        ellipse = Ellipse((delta_x[i][j], delta_y[i][j]), height=major_ax[i][j], width=minor_ax[i][j], angle=phi_ellipse[i][j], edgecolor='black', facecolor='none', linewidth=1.5)
        ax.add_artist(ellipse)

    x_c = ([])
    y_c = ([])
    params = ([])

# Search for epoch with most components / search for suitable components in other epochs / make arrays with positions for component c_i

number_c = ([])
([number_c.append(len(delta_x[i])) for i in range(len(delta_x))])
#ref_pos = delta_x[np.array(number_c).argmax(axis=0)]
ref_pos = delta_x[0]

for i in range(len(ref_pos)):
    if i == 0:
        x_c.append(c_x[np.abs(c_x-ref_pos[i]) == 0])
        z = ([])
        for j in range(len(x_c[i])):
            z.append(c_y[c_x == x_c[i][j]][0] - time_diff[j])
        y_c.append(z)
    else:
        x_c.append(c_x[(np.abs(c_x) - np.abs(ref_pos[i]) > -0.1) & (np.abs(c_x) - np.abs(ref_pos[i]) < 0.4)])
        z = ([])
        for j in range(len(x_c[i])):
            z.append(c_y[c_x == x_c[i][j]][0])
        y_c.append(z)

# Plot identified Cs and make linear fit with shifted positions / colored ellipses still missing

for c in range(len(ref_pos)):
    if len(x_c[c]) == 1:
        plt.plot(x_c[c], y_c[c], marker='.', color=colors[c], linestyle='none', markersize=8, label='C_'+str(c))
    else:
        x_c_neu, y_c_neu  = rotate(x_c[c], y_c[c], -np.deg2rad(90), len(x_c[c]))    # Shift positions for 90degree to avoid infinite slopes

        params, covariance = curve_fit(linear_fit, x_c_neu, y_c_neu)
        errors = np.sqrt(np.diag(covariance))
        plt.plot(x_c[c], y_c[c], marker='.', color=colors[c], linestyle='none', markersize=8, label='C_'+str(c))
        #for j in range(len(x_c[c])):
        #    ellipse = Ellipse((x_c[c][j], y_c[c][j]), height=1, width=1, angle=1, edgecolor=colors[c], facecolor='none', linewidth=1.5)
        #    ax.add_artist(ellipse)
        plt.plot(-linear_fit(x_values, *params), x_values, color=colors[c], linewidth=1)

# Change axes and save plot

plt.yticks(-time_diff, dates)
plt.xlabel('Relative RA / mas')
plt.legend()

plt.xlim(5, -15)
plt.ylim(-35, 3)
plt.tight_layout()
plt.savefig('dynamic_plot.pdf', bbox_inches='tight', pad_inches=0.1)
