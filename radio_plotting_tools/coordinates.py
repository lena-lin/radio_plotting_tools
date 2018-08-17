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

    print(x_ref_pixel)

    x = np.linspace(x_ref_pixel * -x_inc, (x_n_pixel - x_ref_pixel) * x_inc, x_n_pixel)
    y = np.linspace(y_ref_pixel * -y_inc, (y_n_pixel - y_ref_pixel) * y_inc, y_n_pixel)

    if not relative:
        x += x_ref_value
        y += y_ref_value

    return x, y
