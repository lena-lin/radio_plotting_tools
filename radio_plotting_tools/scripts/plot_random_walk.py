
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from glob import glob
from astropy.io import fits
from skimage.segmentation import random_walker, watershed
from radio_plotting_tools.noise import noise_level
from radio_plotting_tools.blob_detection import find_blobs
from radio_plotting_tools.spectral_map import find_peaks_max_y
from radio_plotting_tools.coordinates import get_pixel_coordinates, get_point_coordinates
from tqdm import tqdm


#files_15 = sorted(glob('/home/lena/Dokumente/Radio/NGC1275/VLBA_15GHz/Data/FITS/*/*.fits'))
files_43 = sorted(glob('/home/lena/Dokumente/Radio/NGC1275/VLBA-BU/Data_Web/FITS/20*-*/0316*.IMAP'))

for f in tqdm(files_43):
    clean_map = fits.open(f)[0].data[0][0]
    header = fits.open(f)[0].header
    #print(header['NAXIS1'])
    x_ref, y_ref = find_peaks_max_y(header, clean_map, max_y_offset=200)
    noise = noise_level(clean_map, n_sigma=5)
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

    markers = masked_clean_map.copy()
    markers[markers == 0] = -1
    markers[markers > 0] = 0

    for i, [row, col, r] in enumerate(blobs):
        markers[int(row), int(col)] = i+1
    #print(noise, masked_clean_map.min(), masked_clean_map.max())
    #from IPython import embed; embed()
    labels_rw = random_walker(masked_clean_map, markers)

    fig = plt.figure(figsize=(16, 9))
    ax1 = fig.add_subplot(121)
    ax1.pcolorfast(x_grid, y_grid, masked_clean_map)
    ax1.set_aspect(1)
    ax1.set_xlim(5, -5)
    ax1.set_ylim(-5, 5)

    ax2 = fig.add_subplot(122)
    cax = ax2.pcolorfast(x_grid, y_grid, labels_rw)
    ax2.contour(
                    x_grid,
                    y_grid,
                    masked_clean_map,
                    levels=np.logspace(
                       np.log10(noise),
                       np.log10(masked_clean_map.max()),
                       10
                    ),
                    cmap='Pastel2'
                    )
    for blob in blobs:
        #print(blob)
        row, col, r = blob
        x, y = get_point_coordinates(header, row, col, x_ref_pixel=x_ref, y_ref_pixel=y_ref)
        ax2.scatter(x, y, s=2, c='k')
    ax2.set_aspect(1)
    ax2.set_xlim(5, -5)
    ax2.set_ylim(-5, 5)
    #cbar = fig.colorbar(cax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    #plt.show()
    plt.savefig('/home/lena/Dokumente/Radio/NGC1275/VLBA-BU/Plots/RandomWalker/{}_random_walker.png'.format(header['DATE-OBS']))
    plt.close()
