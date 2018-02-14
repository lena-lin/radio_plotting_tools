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
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit

import pandas as pd

import glob
import plot_clean_map


################################################################################################################################################

# Function to rotate positions of components for linear fit
def rotate(x, y, angle, length):
    ox = np.zeros(length)
    oy = np.zeros(length)
    angle = np.zeros(length) + angle

    qx = ox + np.cos(angle) * (x - ox) - np.sin(angle) * (y - oy)
    qy = oy + np.sin(angle) * (x - ox) + np.cos(angle) * (y - oy)
    return qx, qy

# Function to calculate time differences
def time_difference(dates, shortest_difference, scaling):
    time_diff = ( ([(np.datetime64(i) - np.datetime64(dates[0])) for i in dates]) / np.timedelta64(1, 'D') / shortest_difference ) * scaling
    return time_diff

# Function for linear fit
def linear_fit(f, m, b):
    return m * f + b

# x_values for linear regression
x_values = np.linspace(-40, 100, 10)

# Colors for different c_i
colors = (['red', 'blue', 'magenta', 'orange', 'green', 'brown', 'cyan', 'darkorange', 'grey', 'olive'])



################################################################################################################################################



def make_dynamic_plot(
    input_file,
    ref_epoch,
    data_path
):
    # Array with datapaths for all epochs
    datapaths = ([])
    for data in glob.glob(data_path + '/*/final_mod.fits'):
        datapaths = sorted(np.append(datapaths, data))

    # Load df
    df_components = pd.read_csv(input_file, index_col=0)
    dates = sorted(np.unique(df_components['date']))
    time_diff = np.asarray(time_difference(dates, 60, 2))

    # Choose reference epoch
    ref_pos = df_components['x_positions'][df_components['date'] == ref_epoch]


    # Find components
    c_i = np.arange(len(ref_pos))
    for i in range(len(ref_pos)):
        index = df_components[(np.abs(df_components['x_positions']) - np.abs(ref_pos[i]) > -0.1 )
                               &
                               (np.abs(df_components['x_positions']) - np.abs(ref_pos[i]) < 0.25)].index
        df_components.loc[index, 'c_i'] = c_i[i]

    # Subtract time differences from y_positions
    for i in range(len(dates)):
        df_components.loc[df_components['date'] == dates[i], 'y_positions'] = df_components['y_positions'] - time_diff[i]

    ################################################################################################################################################

    # Plot all components and contours
    fig = plt.figure(figsize=(12,20))
    ax = fig.add_subplot(111, aspect='equal')

    for i in range(len(dates)):
        # Plot components
        plt.plot(df_components.loc[df_components['date'] == dates[i], 'x_positions'],
                 df_components.loc[df_components['date'] == dates[i], 'y_positions'],
                 marker='.',
                 color='black',
                 linestyle='none',
                 markersize=4,
                 label=''
                 )

        # Plot clean map
        plot_clean_map.main(str(datapaths[i]) ,0.007, time_diff[i])

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


    # Plot identified Components and make linear fit with shifted positions
    for c in range(len(ref_pos)):
        if len(df_components[df_components['c_i'] == c]) == 1:
            plt.plot(df_components.loc[df_components['c_i'] == c, 'x_positions'],
                     df_components.loc[df_components['c_i'] == c, 'y_positions'],
                     marker='.',
                     color=colors[c],
                     linestyle='none',
                     markersize=8,
                     label='C_'+str(c)
                     )
        else:
            x_c_neu, y_c_neu  = rotate(df_components['x_positions'][df_components['c_i'] == c],     # Shift positions for 90degree to avoid infinite slopes
                                df_components['y_positions'][df_components['c_i'] == c],
                                -np.deg2rad(90),
                                len(df_components['x_positions'][df_components['c_i'] == c]))

            params, covariance = curve_fit(linear_fit, x_c_neu, y_c_neu)
            errors = np.sqrt(np.diag(covariance))
            plt.plot(df_components.loc[df_components['c_i'] == c, 'x_positions'],
                     df_components.loc[df_components['c_i'] == c, 'y_positions'],
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

    # Save identified components to df
    df_components.to_csv('components.csv')

if __name__ == '__main__':
  make_dynamic_plot()
