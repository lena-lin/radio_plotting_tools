import numpy as np
from astropy.io import fits
from scipy import interpolate
from .coordinates import get_pixel_coordinates

def interpolated_map(file, ra_min=4, ra_max=-4, dec_min=-4, dec_max=4):
    x, y = get_pixel_coordinates(fits.open(file)[0].header)
    cm = fits.open(file)[0].data[0][0]
    func_interpol_cm = interpolate.interp2d(x, y, cm, kind='cubic')

    x = np.linspace(-5, 5, 101)
    y = np.linspace(-5, 5, 101)

    return x, y, func_interpol_cm(x,y)
