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
from scipy.optimize import curve_fit
import matplotlib.dates as mdates
from datetime import datetime

from pandas.tseries import converter
converter.register()


#################################################################################################################################################

# Load data / array with all observation dates

df_components = pd.read_csv('components.csv')
dates = np.unique(df_components['date'])

# Add delta_days to df

date = (pd.to_datetime(df_components['date'], format='%Y-%m-%d'))
delta_days = ((date - date.min()) / np.timedelta64(1,'D'))
delta_days = pd.DataFrame({'delta_days': delta_days})
df_components = pd.concat([df_components, delta_days], axis=1)

# Define linear function
def linear_fit(f, m, b):
    return m * f + b

x_values = np.linspace(16000, 17000, 10)

# Colors for different c_i
colors = (['red', 'blue', 'magenta', 'orange', 'green', 'brown', 'cyan', 'darkorange', 'grey', 'olive'])

#################################################################################################################################################

# Plot radial_dist vs dates
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111)

for i in range(len(df_components.loc[df_components['date'] == '2013-12-15', 'x_positions'])):
    plt.plot(np.array(df_components.loc[df_components['c_i'] == i, 'date'], dtype='datetime64[D]'),
             df_components.loc[df_components['c_i'] == i, 'radial_dist'],
             linestyle='none',
             marker='.',
             markersize=8,
             label='C_'+str(i),
             color=colors[i]
             )
    print(df_components.loc[df_components['c_i'] == i, 'delta_days'])
    print(df_components.loc[df_components['c_i'] == i, 'radial_dist'])
    params, covariance = curve_fit(linear_fit,
                     df_components.loc[df_components['c_i'] == i, 'delta_days'],
                     df_components.loc[df_components['c_i'] == i, 'radial_dist']
                     )
    errors = np.sqrt(np.diag(covariance))
    print(params)
    print(errors)

    plt.plot(np.array(x_values, dtype='datetime64[D]'), linear_fit(x_values, *params), linewidth=0.5, color=colors[i], linestyle='--')
    plt.show()


plt.ylabel('Distance / mas')
plt.xlabel('Date')
plt.xlim(np.array(x_values, dtype='datetime64[D]')[0],np.array(x_values, dtype='datetime64[D]')[-1])
plt.ylim(-1,13)
plt.ylim()
plt.legend(loc=9, ncol=5)
plt.tight_layout()
plt.savefig('jet_motion.pdf', bbox_inches='tight', pad_inches=0.1)
