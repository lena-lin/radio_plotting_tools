import astropy.units as u
import numpy as np
import pandas as pd

def get_gaussian_components(aips_cc_data):
    x_positions = (aips_cc_data['DELTAX'] * u.degree).to(u.mas)
    y_positions = (aips_cc_data['DELTAY'] * u.degree).to(u.mas)
    major_axes = (aips_cc_data['MAJOR AX'] * u.degree).to(u.mas)
    minor_axes = (aips_cc_data['MINOR AX'] * u.degree).to(u.mas)
    phi_ellipse = (aips_cc_data['POSANGLE'])

    # Shift Core to (0,0)
    x_positions = x_positions - x_positions[0] # in case first component fits core!!
    y_positions = y_positions - y_positions[0]

    df_components = pd.DataFrame({'x_positions': x_positions,
                                  'y_positions': y_positions,
                                  'major_axes': major_axes,
                                  'minor_axes': minor_axes,
                                  'phi_ellipse': phi_ellipse
                                  })

    return df_components
