import numpy as np
from astropy.io import fits
from radio_plotting_tools.coordinates import get_mask_from_mas
from radio_plotting_tools.noise import noise_from_header
from radio_plotting_tools.coordinates import find_peaks_max_y


def load_data(file, mask_ra_min=5, mask_ra_max=-5, mask_dec_min=-5, mask_dec_max=5,  n_sigma=3, threshold_=0.009, crop_noise=True, shift_ref=False, max_x_offset=45, max_y_offset=300):
    difmap_data = fits.open(file)
    header = difmap_data[0].header
    clean_map = difmap_data['PRIMARY'].data[0][0]

    if crop_noise is True:
        noise = noise_from_header(header, n_sigma=n_sigma)
        clean_map = (clean_map > noise).astype(int) * clean_map

    if shift_ref is True:

            x_ref, y_ref = find_peaks_max_y(
                            header,
                            clean_map,
                            max_x_offset,
                            max_y_offset,
                            min_sigma=1,
                            max_sigma=50,
                            num_sigma=1,
                            threshold=threshold_,
                            overlap=0.5)

    else:
        x_ref = header['CRPIX1']
        y_ref = header['CRPIX2']

    mask_row_min, mask_row_max, mask_col_min, mask_col_max = get_mask_from_mas(
                                                                                difmap_data[0].header,
                                                                                mask_ra_min,
                                                                                mask_ra_max,
                                                                                mask_dec_min,
                                                                                mask_dec_max,
                                                                                x_ref_pixel=x_ref,
                                                                                y_ref_pixel=y_ref
                                                                            )

    mask_array = np.zeros([clean_map.shape[0], clean_map.shape[1]], dtype=bool) 
    mask_array[mask_row_min:mask_row_max, mask_col_min:mask_col_max] = True
    masked_clean_map = clean_map*mask_array

    return masked_clean_map, difmap_data, x_ref, y_ref
