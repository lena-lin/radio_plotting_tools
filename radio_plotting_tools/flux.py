from astropy.io import fits
import astropy.units as u
from astropy.table import Table
import numpy as np
from datetime import datetime
import pandas as pd


def get_flux_from_clean_map(file, peak_flux=False):

    difmap_data = fits.open(file)
    clean_map = difmap_data[0].data[0][0]
    beam_maj = (difmap_data[0].header['BMAJ']*u.deg).to(u.mas) / 2
    beam_min = (difmap_data[0].header['BMIN']*u.deg).to(u.mas) / 2
    beamsize = beam_maj * beam_min * np.pi
    pixel_inc_ra = (difmap_data[0].header['CDELT1'] * u.deg).to(u.mas)
    pixel_inc_dec = (difmap_data[0].header['CDELT2'] * u.deg).to(u.mas)
    pixelsize = abs(pixel_inc_ra * pixel_inc_dec)
    if peak_flux is True:
        flux = clean_map.max() * (pixelsize/beamsize)
    else:
        flux = clean_map.sum() * (pixelsize/beamsize)
    date = datetime.strptime(difmap_data[0].header['DATE-OBS'].strip(), '%Y-%m-%d')

    return date, flux


def get_flux_from_mod_file(f, peak_flux=False, ignore_negative_components=False):

    df = pd.read_csv(
                f,
                sep=" ",
                skipinitialspace=True,
                skiprows=4,
                names=['Flux', 'Radius', 'Theta', 'nan'],
                engine='python'
            )
    df = df.drop(['nan'], axis=1)

    if ignore_negative_components is True:
        df['Flux'][df['Flux'] < 0] = 0

    if peak_flux is True:
        flux = df['Flux'].max()
    else:
        flux = df['Flux'].sum()

    return flux


def get_flux_components_from_fits(file, abs_ra_max=5, abs_dec_max=5):
    header = fits.open(file)[0].header
    components = fits.open(file)[1].data
    comp_array = Table(
                    [components['DELTAX'], components['DELTAY'], components['FLUX']],
                    names=['ra', 'dec', 'flux']
                )

    comp_array['ra'] = comp_array['ra'] * u.deg.to(u.mas)
    comp_array['dec'] = comp_array['dec'] * u.deg.to(u.mas)

    comp_masked = comp_array[(abs(comp_array['ra']) < abs_ra_max) & (abs(comp_array['dec']) < abs_dec_max)]

    x_inc = header['CDELT1'] * u.deg.to(u.mas)
    y_inc = header['CDELT2'] * u.deg.to(u.mas)

    x_ref_pixel = header['CRPIX1']
    y_ref_pixel = header['CRPIX2']

    comp_masked['col'] = (comp_masked['ra'].data/x_inc + x_ref_pixel).astype(int)
    comp_masked['row'] = (comp_masked['dec'].data/y_inc + y_ref_pixel).astype(int)

    return comp_masked
