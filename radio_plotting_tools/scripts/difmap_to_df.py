import numpy as np
from datetime import datetime
from astropy.io import fits
import astropy.units as u
import glob
import pandas as pd

'''

Ich überleg noch, ob es Sinn macht die CleanMap Daten zwischen zu speichern.
Da das array clean_map nicht so einfach ins df passt, scheint es auf dem ersten Blick keine gute Idee zu sein.

'''

def load_data(
    data_path
):
    # Array with datapaths for all epochs
    datapaths = ([])
    for data in glob.glob(data_path + '/*/final_mod.fits'):
        datapaths = sorted(np.append(datapaths, data))


    # df_clean_map = pd.DataFrame(columns=['x_n_pixel', 'x_ref_pixel', 'x_inc', 'x_ref_value', 'y_n_pixel', 'y_ref_pixel',
    #                                      'y_inc', 'y_ref_value', 'clean_map'])
    #
    # i=0
    # for epoch in sorted(datapaths):
    #     difmap_data = fits.open(epoch)
    #     date = fits.open(epoch)['PRIMARY'].header['DATE-OBS']
    #
    #     x_n_pixel = difmap_data['PRIMARY'].header['NAXIS1']
    #     x_ref_pixel = difmap_data['PRIMARY'].header['CRPIX1']
    #     x_inc = ((difmap_data['PRIMARY'].header['CDELT1'] * u.degree).to(u.mas)).value
    #     x_ref_value = ((difmap_data['PRIMARY'].header['CRVAL1'] * u.degree).to(u.mas)).value
    #     y_n_pixel = difmap_data['PRIMARY'].header['NAXIS2']
    #     y_ref_pixel = difmap_data['PRIMARY'].header['CRPIX2']
    #     y_inc = ((difmap_data['PRIMARY'].header['CDELT2'] * u.degree).to(u.mas)).value
    #     y_ref_value = ((difmap_data['PRIMARY'].header['CRVAL2'] * u.degree).to(u.mas)).value
    #     clean_map = difmap_data['PRIMARY'].data[0][0]
    #     print(x_n_pixel)
    #
    #
    #     df_clean_map_i = pd.DataFrame({'date': date,
    #                               'x_n_pixel': x_n_pixel,
    #                               'x_ref_pixel': x_ref_pixel,
    #                               'x_inc': x_inc,
    #                               'x_ref_value': x_ref_value,
    #                               'y_n_pixel': y_n_pixel,
    #                               'y_ref_pixel': y_ref_pixel,
    #                               'y_inc': y_inc,
    #                               'y_ref_value': y_ref_value
    #                               #'clean_map': clean_map
    #                               }, index=[0])
    #
    #     df_clean_map = pd.concat([df_clean_map, df_clean_map_i], ignore_index=True)
    #     i+=1
    #
    # df_clean_map.to_csv('clean_map.csv')


    df_components = pd.DataFrame(columns=['date', 'x_positions', 'y_positions', 'major_axes', 'minor_axes',
                                          'phi_ellipse', 'c_i', 'radial_dist', 'flux'])

    i = 0
    for epoch in sorted(datapaths):
        date = fits.open(epoch)['PRIMARY'].header['DATE-OBS']
        difmap_data = fits.open(epoch)

        x_positions = ((difmap_data['AIPS CC'].data['DELTAX'] * u.degree).to(u.mas)).value
        y_positions = ((difmap_data['AIPS CC'].data['DELTAY'] * u.degree).to(u.mas)).value
        major_axes = ((difmap_data['AIPS CC'].data['MAJOR AX'] * u.degree).to(u.mas)).value
        minor_axes = ((difmap_data['AIPS CC'].data['MINOR AX'] * u.degree).to(u.mas)).value
        phi_ellipse = difmap_data['AIPS CC'].data['POSANGLE']
        flux = difmap_data['AIPS CC'].data['FLUX']

        # Shift Core to (0,0)
        x_positions_corrected = x_positions - x_positions[0] # in case first component fits core!!
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
        i+=1


    df_components.to_csv('components.csv')

if __name__ == '__main__':
  load_data()
