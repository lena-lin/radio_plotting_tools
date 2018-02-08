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
import pandas as pd

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

# Shift core to (0,0) and calculate timedifference
time_diff = np.asarray(time_difference(dates, 60, 2))


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

########################################################################################################################################

df_components = pd.DataFrame(columns=['date', 'x_positions', 'y_positions', 'major_axes', 'minor_axes',
                                      'phi_ellipse', 'c_i', 'y_no_time', 'radial_dist'])

i = 0
for epoch in sorted(datapaths):
    input_file = epoch
    date = fits.open(epoch)['PRIMARY'].header['DATE-OBS']
    difmap_data = fits.open(input_file)

    x_positions = ((difmap_data['AIPS CC'].data['DELTAX'] * u.degree).to(u.mas)).value
    y_positions = ((difmap_data['AIPS CC'].data['DELTAY'] * u.degree).to(u.mas)).value
    major_axes = ((difmap_data['AIPS CC'].data['MAJOR AX'] * u.degree).to(u.mas)).value
    minor_axes = ((difmap_data['AIPS CC'].data['MINOR AX'] * u.degree).to(u.mas)).value
    phi_ellipse = (difmap_data['AIPS CC'].data['POSANGLE'])

    # Shift Core to (0,0)
    x_positions_corected = x_positions - x_positions[0] # in case first component fits core!!
    y_positions_corected = y_positions - y_positions[0]
    radial_dist = np.sqrt(x_positions**2 + y_positions**2)
    y_positions_time = y_positions - y_positions[0] - time_diff[i]

    df_components_i = pd.DataFrame({'date': date,
                                  'x_positions': x_positions_corected,
                                  'y_positions': y_positions_time,
                                  'major_axes': major_axes,
                                  'minor_axes': minor_axes,
                                  'phi_ellipse': phi_ellipse,
                                  'c_i': '',
                                  'y_no_time': y_positions_corected,
                                  'radial_dist': radial_dist
                                  })
    df_components = pd.concat([df_components, df_components_i], ignore_index=True)
    i+=1

# Choose reference epoch
ref_pos = df_components['x_positions'][df_components['date'] == '2013-12-15']

# Find components
c_i = np.arange(len(ref_pos))
for i in range(len(ref_pos)):
    index = df_components[(np.abs(df_components['x_positions']) - np.abs(ref_pos[i]) > -0.1 ) & (np.abs(df_components['x_positions']) - np.abs(ref_pos[i]) < 0.25)].index
    df_components.loc[index, 'c_i'] = c_i[i]

# Make contour plot and plot all components
fig = plt.figure(figsize=(12,20))
ax = fig.add_subplot(111, aspect='equal')

for i in range(len(dates)):

    x = np.linspace(x_ref_pixel[i] * x_inc[i], -x_ref_pixel[i] * x_inc[i], x_n_pixel[i])
    y = np.linspace(-y_ref_pixel[i] * y_inc[i], y_ref_pixel[i] * y_inc[i], y_n_pixel[i]).value - time_diff[i]
    flux_min.append(clean_map[i].min())
    flux_max.append(clean_map[i].max())
    loglevs_min.append(flux_min[i] + (flux_max[i] - flux_min[i]) * 0.007)



    plt.plot(df_components.loc[df_components['date'] == dates[i], 'x_positions'], df_components.loc[df_components['date'] == dates[i], 'y_positions'],
             marker='.', color='black', linestyle='none', markersize=4, label='')

    ax.set_aspect(1)
    ax.contour(-x, y, clean_map[i], norm=LogNorm(), levels=np.logspace(np.log10(loglevs_min[i]), np.log10(flux_max[i]), 10))

    # Plot ellipses
    for j in df_components['x_positions'][df_components['date'] == dates[i]].index:
        ellipse = Ellipse((df_components['x_positions'][df_components['date'] == dates[i]][j],
                               df_components['y_positions'][df_components['date'] == dates[i]][j]),
                               height = df_components['minor_axes'][df_components['date'] == dates[i]][j],
                               width = df_components['major_axes'][df_components['date'] == dates[i]][j],
                               angle = df_components['phi_ellipse'][df_components['date'] == dates[i]][j],
                               edgecolor='black',
                               facecolor='none',
                               linewidth=1)
        ax.add_artist(ellipse)


# Plot identified Cs and make linear fit with shifted positions / colored ellipses still missing

for c in range(len(ref_pos)):
    if len(df_components[df_components['c_i'] == c]) == 1:
        plt.plot(df_components['x_positions'][df_components['c_i'] == c], df_components['y_positions'][df_components['c_i'] == c], marker='.', color=colors[c], linestyle='none', markersize=8, label='C_'+str(c))
    else:
        x_c_neu, y_c_neu  = rotate(df_components['x_positions'][df_components['c_i'] == c],     # Shift positions for 90degree to avoid infinite slopes
                            df_components['y_positions'][df_components['c_i'] == c],
                            -np.deg2rad(90),
                            len(df_components['x_positions'][df_components['c_i'] == c]))

        params, covariance = curve_fit(linear_fit, x_c_neu, y_c_neu)
        errors = np.sqrt(np.diag(covariance))
        plt.plot(df_components['x_positions'][df_components['c_i'] == c],
                 df_components['y_positions'][df_components['c_i'] == c],
                 marker='.',
                 color=colors[c],
                 linestyle='none',
                 markersize=8,
                 label='C_'+str(c))
        for j in df_components['x_positions'][df_components['c_i'] == c].index:
            ellipse = Ellipse((df_components['x_positions'][df_components['c_i'] == c][j],
                               df_components['y_positions'][df_components['c_i'] == c][j]),
                               height = df_components['minor_axes'][df_components['c_i'] == c][j],
                               width = df_components['major_axes'][df_components['c_i'] == c][j],
                               angle = df_components['phi_ellipse'][df_components['c_i'] == c][j],
                               edgecolor=colors[c],
                               facecolor='none',
                               linewidth=1.5)
            ax.add_artist(ellipse)
        plt.plot(-linear_fit(x_values, *params), x_values, color=colors[c], linewidth=1)

# Change axes and save plot

plt.yticks(-time_diff, dates)
plt.xlabel('Relative RA / mas')
plt.legend()

plt.xlim(5, -15)
plt.ylim(-35, 3)
plt.tight_layout()
plt.savefig('dynamic_plot.pdf', bbox_inches='tight', pad_inches=0.1)

df_components.to_csv('components.csv')
