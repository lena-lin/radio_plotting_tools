import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
from scipy import interpolate
from .coordinates import get_pixel_grid_coordinates
from radio_plotting_tools.coordinates import find_peaks_max_y

from radio_plotting_tools.noise import noise_level


def interpolated_map(header, clean_map, ra=[5, -5], dec=[-5, 5], convolution_beam_header=None):
    x_ref, y_ref = find_peaks_max_y(header, clean_map, max_x_offset=None, max_y_offset=None)

    noise = noise_level(clean_map)
    clean_map = (clean_map > noise).astype(int) * clean_map

    if convolution_beam_header is not None:
        bmaj = convolution_beam_header['BMAJ']/abs(convolution_beam_header['CDELT1'])
        bmin = convolution_beam_header['BMIN']/abs(convolution_beam_header['CDELT2'])
        bpa = convolution_beam_header['BPA']
        clean_map = convolve_cleanmap(clean_map, bmin, bmaj, bpa)

    x, y = get_pixel_grid_coordinates(header, x_ref_pixel=x_ref, y_ref_pixel=y_ref)
    func_interpol_cm = interpolate.interp2d(x, y, clean_map, kind='cubic')

    ra_linspace = np.linspace(ra[0], ra[1], 501)
    dec_linspace = np.linspace(dec[0], dec[1], 501)

    return ra_linspace, dec_linspace, func_interpol_cm(ra_linspace, dec_linspace)





def beamsize(header):
    major = (header['BMAJ']/2.) * u.deg.to(u.mas)
    minor = (header['BMIN']/2.) * u.deg.to(u.mas)

    return major * minor * np.pi


def spectral_indices(header1, data1, header2, data2, epsilon=1e-3):
    _, _, interp_map_1 = interpolated_map(header1, data1)
    _, _, interp_map_2 = interpolated_map(header2, data2, convolution_beam_header=header1)

    interp_map_1[interp_map_1 < epsilon] = np.nan
    interp_map_2[interp_map_2 < epsilon] = np.nan

    beam1 = beamsize(header1)
    beam2 = beamsize(header2)

    freq1 = header1['CRVAL3']
    freq2 = header2['CRVAL3']

    spectral_indices = np.log((interp_map_2 / interp_map_1) * (beam1 / beam2)) / np.log(freq2/freq1)

    return spectral_indices


def convolve_cleanmap(clean_map, beam_min_fhwm, beam_maj_fwhm, beam_pos_angle):
    kernel = Gaussian2DKernel(x_stddev=beam_min_fhwm/2.3548, y_stddev=beam_maj_fwhm/2.3548, theta=beam_pos_angle)

    return convolve(clean_map, kernel)
