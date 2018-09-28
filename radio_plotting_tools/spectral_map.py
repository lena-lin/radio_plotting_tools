import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
from scipy import interpolate
from .coordinates import get_pixel_coordinates
from skimage.feature import peak_local_max
from radio_plotting_tools.noise import noise_level


def interpolated_map(header, clean_map, offset, ra=[5, -5], dec=[-5, 5], convolution_beam_header=None):
    x_ref, y_ref = find_peaks_max_y(header, clean_map, max_y_offset=offset)

    noise = noise_level(clean_map)
    clean_map = (clean_map > noise).astype(int) * clean_map

    if convolution_beam_header is not None:
        bmaj = convolution_beam_header['BMAJ']/abs(convolution_beam_header['CDELT1'])
        bmin = convolution_beam_header['BMIN']/abs(convolution_beam_header['CDELT2'])
        bpa = convolution_beam_header['BPA']
        clean_map = convolve_cleanmap(clean_map, bmin, bmaj, bpa)

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


def find_peaks_max_y(header, clean_map, max_y_offset, sigma=150):
    threshold = sigma * header['NOISE']
    x_ref = header['CRPIX1']
    y_ref = header['CRPIX2']

    x_header = x_ref
    y_header = y_ref
    peaks = peak_local_max(clean_map, min_distance=5, threshold_abs=threshold)
    max_y = y_ref
    for peak in peaks:
        if (peak[0] > max_y) and (abs(peak[0] - y_header) < max_y_offset) and (abs(peak[1] - x_header) < 10):
            max_y = peak[0]
            y_ref, x_ref = peak

    return int(x_ref), int(y_ref)


def beamsize(header):
    major = (header['BMAJ']/2.) * u.deg.to(u.mas)
    minor = (header['BMIN']/2.) * u.deg.to(u.mas)

    return major * minor * np.pi


def spectral_indices(header1, data1, header2, data2, offset1, offset2, epsilon=1e-3):
    _, _, interp_map_1 = interpolated_map(header1, data1, offset=offset1)
    _, _, interp_map_2 = interpolated_map(header2, data2, offset=offset2, convolution_beam_header=header1)

    interp_map_1[interp_map_1 < epsilon] = np.nan
    interp_map_2[interp_map_2 < epsilon] = np.nan

    beam1 = beamsize(header1)
    beam2 = beamsize(header2)

    freq1 = header1['CRVAL3']
    freq2 = header2['CRVAL3']

    spectral_indices = np.log((interp_map_1 / interp_map_2) * (beam2 / beam1)) / np.log(freq2/freq1)

    return spectral_indices


def convolve_cleanmap(clean_map, beam_min_fhwm, beam_maj_fwhm, beam_pos_angle):
    kernel = Gaussian2DKernel(x_stddev=beam_min_fhwm/2.3548, y_stddev=beam_maj_fwhm/2.3548, theta=beam_pos_angle)

    return convolve(clean_map, kernel)
