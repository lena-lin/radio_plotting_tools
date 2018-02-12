import numpy as np


import uncertainties.unumpy as unp
import scipy.constants as const


df_components = pd.read_csv('components.csv')

lamb = const.c / 15e9
kB = unp.uarray(const.physical_constants['Boltzmann constant'][0], const.physical_constants['Boltzmann constant'][2])
S = df_components['']
z = 0.078
a_maj = df_components['major_axes']
a_min = df_components['minor_axes']


bright_temp = (2 * np.log(2) / np.pi * kB) * (S * lamb**2 (1 + z)) / (a_maj * a_min)

print(bright_temp)
