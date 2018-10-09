import astropy.units as u
import numpy as np
from skimage.feature import peak_local_max
from skimage.feature import blob_log



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


def get_mask_from_mas(header, ra_min, ra_max, dec_min, dec_max, x_ref_pixel=None, y_ref_pixel=None):
    if x_ref_pixel == None:
        x_ref_pixel = header['CRPIX1']
    if y_ref_pixel == None:
        y_ref_pixel = header['CRPIX2']

    x_inc = header['CDELT1'] * u.deg.to(u.mas)
    y_inc = header['CDELT2'] * u.deg.to(u.mas)

    mask_row_min = int(y_ref_pixel - abs(dec_min/y_inc))
    mask_row_max = int(y_ref_pixel + abs(dec_max/y_inc))
    mask_col_min = int(x_ref_pixel - abs(ra_min/x_inc))
    mask_col_max = int(x_ref_pixel + abs(ra_max/x_inc))

    return mask_row_min, mask_row_max, mask_col_min, mask_col_max


def get_point_coordinates(header, row, col, x_ref_pixel=None, y_ref_pixel=None):
    if x_ref_pixel == None:
        x_ref_pixel = header['CRPIX1']
    if y_ref_pixel == None:
        y_ref_pixel = header['CRPIX2']

    x_inc = (header['CDELT1'] * u.degree).to(u.mas)
    y_inc = (header['CDELT2'] * u.degree).to(u.mas)

    ra = (col - x_ref_pixel)*x_inc.value
    dec = (row - y_ref_pixel)*y_inc.value

    return ra, dec


def find_peaks_max_y(header, clean_map, max_y_offset, sigma=150, max_sigma_=50, num_sigma_=1, threshold_=0.001):
    #threshold = sigma * header['NOISE']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']

    x_header = x_ref
    y_header = y_ref
    #peaks = peak_local_max(clean_map, min_distance=1, threshold_abs=threshold)

    peaks = blob_log(
                    clean_map/abs(clean_map).max(),
                    max_sigma=max_sigma_,
                    num_sigma=num_sigma_,
                    threshold=threshold_,
                    )

    max_y = y_ref
    for peak in peaks:
        if (peak[0] > max_y) and (abs(peak[0] - y_header) < max_y_offset) and (abs(peak[1] - x_header) < 5):
            max_y = peak[0]
            y_ref = peak[0]
            x_ref = peak[1]

    return x_ref, y_ref
