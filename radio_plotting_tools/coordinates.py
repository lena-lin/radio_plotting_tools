import astropy.units as u
import numpy as np


def get_pixel_coordinates(header, x_ref_pixel=None, y_ref_pixel=None, relative=True):
    x_n_pixel = header['NAXIS1']

    x_inc = (header['CDELT1'] * u.degree).to(u.mas)
    x_ref_value = (header['CRVAL1'] * u.degree).to(u.mas)
    y_n_pixel = header['NAXIS2']

    y_inc = (header['CDELT2'] * u.degree).to(u.mas)
    y_ref_value = (header['CRVAL2'] * u.degree).to(u.mas)

    if x_ref_pixel == None:
        x_ref_pixel = header['CRPIX1']
    if y_ref_pixel == None:
        y_ref_pixel = header['CRPIX2']

    x = np.linspace(x_ref_pixel * -x_inc, (x_n_pixel - x_ref_pixel) * x_inc, x_n_pixel)
    y = np.linspace(y_ref_pixel * -y_inc, (y_n_pixel - y_ref_pixel) * y_inc, y_n_pixel)

    if not relative:
        x += x_ref_value
        y += y_ref_value

    return x, y


def relative_position(ra, dec, ra_0, dec_0):
    ra_rel = (ra - ra_0)*u.deg.to(u.mas)
    dec_rel = (dec - dec_0)*u.deg.to(u.mas)
    return ra_rel, dec_rel


def get_mask_from_mas(header, ra_min, ra_max, dec_min, dec_max):
    x_ref_pixel = header['CRPIX1']
    x_inc = (header['CDELT1'] * u.degree).to(u.mas)
    y_ref_pixel = header['CRPIX2']
    y_inc = (header['CDELT2'] * u.degree).to(u.mas)

    mask_row_min = int(y_ref_pixel - abs(dec_min/y_inc.value))
    mask_row_max = int(y_ref_pixel + abs(dec_max/y_inc.value))
    mask_col_min = int(x_ref_pixel - abs(ra_min/x_inc.value))
    mask_col_max = int(x_ref_pixel + abs(ra_max/x_inc.value))

    return mask_row_min, mask_row_max, mask_col_min, mask_col_max


def get_point_coordinates(header, row, col, x_ref_pixel=None, y_ref_pixel=None,):
    x_inc = (header['CDELT1'] * u.degree).to(u.mas)
    y_inc = (header['CDELT2'] * u.degree).to(u.mas)

    if x_ref_pixel == None:
        x_ref_pixel = header['CRPIX1']
    if y_ref_pixel == None:
        y_ref_pixel = header['CRPIX2']

    ra = (col - x_ref_pixel)*x_inc.value
    dec = (row - y_ref_pixel)*y_inc.value

    return ra, dec
