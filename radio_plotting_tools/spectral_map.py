import numpy as np
from astropy.io import fits
from scipy import interpolate
from .coordinates import get_pixel_coordinates
from skimage.feature import peak_local_max


def interpolated_map(clean_map, header, x_ref=None, y_ref=None, ra=[5, -5], dec=[-5, 5]):

    x, y = get_pixel_coordinates(header, x_ref_pixel=x_ref, y_ref_pixel=y_ref)
    func_interpol_cm = interpolate.interp2d(x, y, clean_map, kind='cubic')

    ra_linspace = np.linspace(ra[0], ra[1], 501)
    dec_linspace = np.linspace(dec[0], dec[1], 501)

    return ra_linspace, dec_linspace, func_interpol_cm(ra_linspace, dec_linspace)


def find_peaks(file, c1_x_ref=1054, c1_y_ref=1026, sigma=3000):
    header = fits.open(file)[0].header
    clean_map = fits.open(file)[0].data[0][0]

    threshold = sigma * header['NOISE']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']

    peaks = peak_local_max(clean_map, min_distance=10, threshold_abs=threshold)

    dist = 1000
    for peak in peaks:
        if abs(np.subtract(peak, [c1_x_ref, c1_y_ref])).sum() < dist:
            x_ref, y_ref = peak
            dist = abs(np.subtract(peak, [c1_x_ref, c1_y_ref])).sum()

    return y_ref, x_ref
