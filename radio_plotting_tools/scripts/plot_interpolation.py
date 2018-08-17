import numpy as np
from radio_plotting_tools.spectral_map import interpolated_map
import click
import matplotlib.pyplot as plt
from astropy.io import fits
from radio_plotting_tools.coordinates import get_pixel_coordinates
from radio_plotting_tools.spectral_map import find_peaks
from radio_plotting_tools.noise import noise_level


@click.command()
@click.argument('file1', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('file2', type=click.Path(file_okay=True, dir_okay=False))
def main(
    file1,
    file2,
    epsilon=1e-3,
):
    x_ref, y_ref = find_peaks(file2)

    fits_file1 = fits.open(file1)[0]
    clean_map1 = fits_file1.data[0][0]
    header1 = fits_file1.header
    noise1 = noise_level(clean_map1)

    fits_file2 = fits.open(file2)[0]
    clean_map2 = fits_file2.data[0][0]
    header2 = fits_file2.header
    noise2 = noise_level(clean_map2, window_width=300)

    clean_map1 = (clean_map1 > noise1).astype(int) * clean_map1
    clean_map2 = (clean_map2 > noise2).astype(int) * clean_map2

    x1, y1, interp_map_1 = interpolated_map(clean_map1, header1)
    x2, y2, interp_map_2 = interpolated_map(clean_map2, header2, x_ref=x_ref, y_ref=y_ref)

    interp_map_1[interp_map_1 < epsilon] = np.nan
    interp_map_2[interp_map_2 < epsilon] = np.nan

    spectral_indices = np.log(interp_map_1 / interp_map_2) / np.log(43/15)

    #spectral_indices = spectral_indices * (clean_map1 > noise1).astype(int) * (clean_map2 > noise2).astype(int)

    from IPython import embed; embed()
    #x1_cm, y1_cm = get_pixel_coordinates(fits.open(file1)[0].header)
    #x2_cm, y2_cm = get_pixel_coordinates(fits.open(file2)[0].header, x_ref_pixel=x_ref, y_ref_pixel=y_ref)

    #from IPython import embed; embed()

    fig = plt.figure()

    ax1 = plt.subplot(131)
    ax1.pcolorfast(-x1, y1, interp_map_1)
    ax1.invert_xaxis()
    ax1.set_aspect(1)

    ax2 = plt.subplot(132)
    ax2.pcolorfast(-x2, y2, interp_map_2)
    ax2.invert_xaxis()
    ax2.set_aspect(1)

    ax3 = plt.subplot(133)
    cax = ax3.pcolorfast(-x2, y2, spectral_indices)
    ax3.invert_xaxis()
    ax3.set_aspect(1)
    plt.tight_layout()
    plt.colorbar(cax, fraction=0.046, pad=0.04)

    # ax3 = plt.subplot(223)
    # ax3.pcolorfast(x1_cm, y1_cm, fits.open(file1)[0].data[0][0])
    # ax3.invert_xaxis()
    # ax3.set_aspect(1)
    # ax4 = plt.subplot(224)
    # ax4.pcolorfast(x2_cm, y2_cm, fits.open(file2)[0].data[0][0])
    # ax4.invert_xaxis()
    # ax4.set_aspect(1)

    plt.show()


if __name__ == '__main__':
    main()
