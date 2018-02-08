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

import pandas as pd
import glob
from astropy.io import fits
from scipy.optimize import curve_fit

'''
Change structure after saving r in df

Infos the programm needs when we write it as a function:

- Dataframe to use
- Reference epoch
- Output filenames

'''

#################################################################################################################################################

# Define functions

def linear_fit(f, m, b):
    return m * f + b

# Colors for different c_i
colors = (['red', 'blue', 'magenta', 'orange', 'green', 'brown', 'cyan', 'darkorange', 'grey', 'olive'])

#################################################################################################################################################

# Load data / array with all observation dates

df_components = pd.read_csv('components.csv')
dates = np.unique(df_components['date'])

#################################################################################################################################################

# Plot major_axis vs distance without corrections
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111)

for i in range(len(df_components.loc[df_components['date'] == '2013-12-15', 'x_positions'])):
    plt.plot(df_components.loc[df_components['c_i'] == i, 'radial_dist'],
             df_components.loc[df_components['c_i'] == i, 'major_axes'],
             linestyle='none',
             marker='.',
             markersize=8,
             label='C_'+str(i),
             color=colors[i]
             )

plt.ylabel('Majoraxis / mas')
plt.xlabel('Distance / mas')
plt.xlim(-1, 12)
plt.ylim(0.0000001,10)
plt.yscale('log')
plt.legend(loc=4)
plt.tight_layout()
plt.savefig('maj_vs_dist_no_correction.pdf', bbox_inches='tight', pad_inches=0.1)
plt.close()

#####################################################################################################################################################

# Linear regression with not diverging points

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111)

r = np.sqrt(df_components.loc[(df_components['c_i'] >= 0) & (df_components['major_axes'] > 1e-3), 'x_positions']**2 +
            df_components.loc[(df_components['c_i'] >= 0) & (df_components['major_axes'] > 1e-3), 'y_no_time']**2)

params, covariance = curve_fit(linear_fit,
                     df_components.loc[(df_components['c_i'] >= 0) & (df_components['major_axes'] > 1e-3), 'radial_dist'],
                     df_components.loc[(df_components['c_i'] >= 0) & (df_components['major_axes'] > 1e-3) , 'major_axes']
                     )

errors = np.sqrt(np.diag(covariance))

x_values = np.linspace(-5, 15, 1000)
for i in range(len(df_components.loc[df_components['date'] == '2013-12-15', 'x_positions'])):
    plt.plot(df_components.loc[(df_components['c_i'] == i) & (df_components['major_axes'] > 1e-3), 'radial_dist'],
             df_components.loc[(df_components['c_i'] == i) & (df_components['major_axes'] > 1e-3), 'major_axes'],
             linestyle='none',
             marker='.',
             markersize=8,
             label='C_'+str(i),
             color=colors[i]
             )

plt.plot(x_values, linear_fit(x_values, *params), linewidth=1.5, color='black', linestyle='--')

plt.xlim(-1, 12)
plt.ylim(0.01,10)
plt.yscale('log')
plt.ylabel('Majoraxis / mas')
plt.xlabel('Distance / mas')
plt.legend(loc=4)
plt.tight_layout()
plt.savefig('maj_vs_dist_regression.pdf', bbox_inches='tight', pad_inches=0.1)
plt.close()


#####################################################################################################################################################

# Plot major_axis vs distance with corrections

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111)

for i in range(len(df_components.loc[df_components['date'] == '2013-12-15', 'x_positions'])):
    plt.plot(df_components.loc[(df_components['c_i'] == i) & (df_components['major_axes'] > 1e-3), 'radial_dist'],
             df_components.loc[(df_components['c_i'] == i) & (df_components['major_axes'] > 1e-3), 'major_axes'],
             linestyle='none',
             marker='.',
             markersize=8,
             label='C_'+str(i),
             color=colors[i]
             )

    r_2 = np.sqrt(df_components.loc[(df_components['c_i'] == i) & (df_components['major_axes'] < 1e-3), 'x_positions']**2 +
                df_components.loc[(df_components['c_i'] == i) & (df_components['major_axes'] < 1e-3), 'y_no_time']**2)
    plt.plot(r_2,
         linear_fit(df_components.loc[(df_components['c_i'] == i) & (df_components['major_axes'] < 1e-3), 'major_axes'], *params),
         linestyle='none',
         marker='.',
         markersize=8,
         color=colors[i],
         label=''
         )

plt.plot(x_values, linear_fit(x_values, *params), linewidth=1.5, color='black', linestyle='--')

plt.xlim(-1, 12)
plt.ylim(0.01,10)
plt.yscale('log')
plt.ylabel('Majoraxis / mas')
plt.xlabel('Distance / mas')
plt.legend(loc=4)
plt.tight_layout()
plt.savefig('maj_vs_dist_corrected.pdf', bbox_inches='tight', pad_inches=0.1)
plt.close()
