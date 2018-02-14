import numpy as np
from astropy.io import fits
import astropy.units as u
import glob
import pandas as pd

'''

Get DataFrame with parameters of Gaussian components from all epochs of a source, sorted by date.
inputs:
input_path: path to folder containing fits files of all epochs in subfolders
file_name: fits files containig clean map and model components must have the same name for every epoch!

'''


def get_gaussian_components_from_fits(
    input_path,
    file_name
):

    datapaths = sorted(glob.glob('{}/*/{}'.format(input_path, file_name)))
    df_components = pd.DataFrame()

    for epoch in datapaths:

        difmap_data = fits.open(epoch)
        gaussian_components_data = difmap_data['AIPS CC'].data

        date = difmap_data['PRIMARY'].header['DATE-OBS']
        x_positions = ((gaussian_components_data['DELTAX'] * u.degree).to(u.mas)).value
        y_positions = ((gaussian_components_data['DELTAY'] * u.degree).to(u.mas)).value
        major_axes = ((gaussian_components_data['MAJOR AX'] * u.degree).to(u.mas)).value
        minor_axes = ((gaussian_components_data['MINOR AX'] * u.degree).to(u.mas)).value
        phi_ellipse = gaussian_components_data['POSANGLE']
        flux = gaussian_components_data['FLUX']

        # Shift Core to (0,0)
        x_positions_corrected = x_positions - x_positions[0]  # in case first component fits core!!
        y_positions_corrected = y_positions - y_positions[0]
        radial_dist = np.sqrt(x_positions_corrected**2 + y_positions_corrected**2)

        df_components_i = pd.DataFrame({'date': date,
                                        'x_positions': x_positions_corrected,
                                        'y_positions': y_positions_corrected,
                                        'major_axes': major_axes,
                                        'minor_axes': minor_axes,
                                        'phi_ellipse': phi_ellipse,
                                        'c_i': '',
                                        'radial_dist': radial_dist,
                                        'flux': flux
                                        })

        df_components = pd.concat([df_components, df_components_i], ignore_index=True)

    return df_components
