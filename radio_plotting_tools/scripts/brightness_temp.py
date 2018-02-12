import numpy as np

import pandas as pd
import uncertainties.unumpy as unp
import scipy.constants as const


df_components = pd.read_csv('components.csv')

lamb = const.c / 15e9
kB = unp.uarray(const.physical_constants['Boltzmann constant'][0], const.physical_constants['Boltzmann constant'][2])
S = df_components.loc[df_components['c_i'] == 0, 'flux']
z = 0.078
a_maj = df_components.loc[df_components['c_i'] == 0, 'major_axes']
a_min = df_components.loc[df_components['c_i'] == 0, 'minor_axes']

bright_temp = (2 * np.log(2) / np.pi * kB) * ( (S * lamb**2 * (1 + z)) / (a_maj * a_min) )

print(bright_temp)

lamb = const.c / 8.4e9
print((2 * np.log(2) / np.pi * kB) * ( (1.29 * lamb**2 * (1 + 0.056)) / (0.095/2 * 0.066/2) ))
