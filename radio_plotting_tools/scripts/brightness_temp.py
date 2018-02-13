import numpy as np

import pandas as pd
import uncertainties.unumpy as unp
import scipy.constants as const
import astropy.units as u



def brightness_temp(
    input_file,
    redshift,
    frequency=15e9
):
    df_components = pd.read_csv(input_file)

    lamb = const.c / frequency
    kB = unp.uarray(const.physical_constants['Boltzmann constant'][0], const.physical_constants['Boltzmann constant'][2])
    S = (df_components.loc[df_components['c_i'] == 0, 'flux'])
    z = redshift
    a_maj = u.mas.to(u.rad, df_components.loc[df_components['c_i'] == 0, 'major_axes'])
    a_min = u.mas.to(u.rad, df_components.loc[df_components['c_i'] == 0, 'minor_axes'])


    bright_temp = (2 * np.log(2) / (np.pi * kB)) * ( (S*1e-26 * lamb**2 * (1 + z)) / ( a_maj * a_min ))

    df_bright_temp = pd.DataFrame({'bright_temp': bright_temp})

    df_bright_temp.to_csv('bright_temp.csv')



if __name__ == '__main__':
    brightness_temp()
