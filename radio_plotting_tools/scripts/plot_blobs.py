import matplotlib.pyplot as plt
from glob import glob
from astropy.io import fits
from radio_plotting_tools.noise import noise_level
from radio_plotting_tools.blob_detection import find_blobs
from radio_plotting_tools.spectral_map import find_peaks_max_y
from radio_plotting_tools.coordinates import get_pixel_coordinates, get_point_coordinates
from tqdm import tqdm


files_15 = sorted(glob('/home/lena/Dokumente/Radio/NGC1275/VLBA_15GHz/Data/FITS/201*/*.fits'))

for f in tqdm(files_15[:3]):
    clean_map = fits.open(f)[0].data[0][0]
    header = fits.open(f)[0].header
    x_ref, y_ref = find_peaks_max_y(f, max_y_offset=40)
    noise = noise_level(clean_map)
    clean_map = (clean_map > noise).astype(int) * clean_map

    x, y = get_pixel_coordinates(header, x_ref_pixel=x_ref, y_ref_pixel=y_ref)

    blobs = find_blobs(clean_map)
    fig, ax = plt.subplots()
    ax.pcolorfast(-x, y, clean_map)
    for blob in blobs:
        col, row, r = blob
        x, y = get_point_coordinates(header, row, col)
        ax.scatter(row, col, s=r)
        ax.set_aspect(1)
        ax.set_xlim(5, -5)
        ax.set_ylim(-5, 5)

        plt.savefig('/home/lena/Dokumente/Radio/NGC1275/VLBA_15GHz/Plots/Blobs/{}_blobs.png'.format(header['DATE-OBS']))
