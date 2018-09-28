
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from glob import glob
from astropy.io import fits
from radio_plotting_tools.noise import noise_level
from radio_plotting_tools.blob_detection import find_blobs
from radio_plotting_tools.spectral_map import find_peaks_max_y
from radio_plotting_tools.coordinates import get_pixel_coordinates, get_point_coordinates
from tqdm import tqdm


files_15 = sorted(glob('/home/lena/Dokumente/Radio/NGC1275/VLBA_15GHz/Data/FITS/201*/*.fits'))

for f in tqdm(files_15[-10:-9]):
    clean_map = fits.open(f)[0].data[0][0]
    header = fits.open(f)[0].header
    print(header['NAXIS1'])
    x_ref, y_ref = find_peaks_max_y(header, clean_map, max_y_offset=40)
    noise = noise_level(clean_map)
    clean_map = (clean_map > noise).astype(int) * clean_map
    mask_array = np.zeros([clean_map.shape[0], clean_map.shape[1]], dtype=bool)
    mask_center_dist = int(abs(5./(header['CDELT1'] * u.deg.to(u.mas))))
    row_min = y_ref - mask_center_dist
    row_max = y_ref + mask_center_dist
    col_min = x_ref - mask_center_dist
    col_max = x_ref + mask_center_dist
    mask_array[row_min:row_max, col_min:col_max] = True
    masked_clean_map = clean_map*mask_array

    x_grid, y_grid = get_pixel_coordinates(header, x_ref_pixel=x_ref, y_ref_pixel=y_ref)
    blobs = find_blobs(masked_clean_map)
    from IPython import embed; embed()
    fig, ax = plt.subplots()
    ax.pcolorfast(x_grid, y_grid, clean_map)
    for blob in blobs:
        print(blob)
        row, col, r = blob
        x, y = get_point_coordinates(header, row, col, x_ref_pixel=x_ref, y_ref_pixel=y_ref)
        ax.scatter(x, y, s=r)
    ax.set_aspect(1)
    ax.set_xlim(5, -5)
    ax.set_ylim(-5, 5)

    plt.savefig('/home/lena/Dokumente/Radio/NGC1275/VLBA_15GHz/Plots/Blobs/{}_blobs.png'.format(header['DATE-OBS']))
    plt.close()
