import click
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.table import Table
from matplotlib.patches import Ellipse
import matplotlib as mpl
from radio_plotting_tools.io import load_data
from radio_plotting_tools.coordinates import get_point_coordinates, get_pixel_grid_coordinates
from radio_plotting_tools.noise import noise_from_header
from glob import glob
from tqdm import tqdm

import matplotlib.pyplot as plt


def plot_clean_map(input_file, output_path, ra_min, ra_max, dec_min, dec_max):



    masked_clean_map, difmap_data, x_ref, y_ref = load_data(input_file, ra_min, ra_max, dec_min, dec_max, threshold_=0.005, n_sigma=5, shift_ref=True)
    header = difmap_data[0].header
    date = header['DATE-OBS']
    source = header['OBJECT']

    x_grid, y_grid = get_pixel_grid_coordinates(header, x_ref_pixel=x_ref, y_ref_pixel=y_ref)
    noise = noise_from_header(header, n_sigma=5)

    fig, ax1 = plt.subplots(1, 1, figsize=(10, 10))

    cax = ax1.pcolorfast(x_grid, y_grid, masked_clean_map, cmap='afmhot')
    ax1.set_aspect(1)
    ax1.set_xlim(ra_min, ra_max)
    ax1.set_ylim(dec_min, dec_max)
    ax1.contour(
                x_grid,
                y_grid,
                masked_clean_map,
                levels=np.logspace(
                    np.log10(noise),
                    np.log10(masked_clean_map.max()),
                    5
                ),
                cmap='YlOrBr',
                linewidths=1,
                )


    ax1.set_xlabel('Relative Right Ascension in mas')
    ax1.set_ylabel('Relative Declination in mas')
    # ax1.text(
    #     -2.7,
    #     1.1,
    #     date,
    #     color='lightyellow',
    #     fontsize=16)
    # ax1.text(
    #     2,
    #     1.1,
    #     source,
    #     color='lightyellow',
    #     fontsize=16)


    cbar = fig.colorbar(cax, fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel('Flux in Jy/beam', fontsize=16)
    ax1.set_xlabel('Relative Right Ascension in mas', fontsize=16)
    ax1.set_ylabel('Relative Declination in mas', fontsize=16)
    ax1.set_aspect(1)
    plt.tight_layout()
    # plt.show()
    # plt.savefig('{}/{}_{}_cleanmap.pdf'.format(output_path, date, source))
    plt.savefig('{}/{}_{}_cleanmap.png'.format(output_path, date, source), dpi=150)
    plt.close()


@click.command()
@click.argument('input_wildcard', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('output_path', type=click.Path(file_okay=False, dir_okay=True))
@click.option(
    '--ra_min',
    type=click.INT,
    help='relative RA min for plotting',
    default='5'
)
@click.option(
    '--ra_max',
    type=click.INT,
    help='relative RA max for plotting',
    default='-5'
)
@click.option(
    '--dec_min',
    type=click.INT,
    help='relative DEC min for plotting',
    default='-5'
)
@click.option(
    '--dec_max',
    type=click.INT,
    help='relative DEC max for plotting',
    default='5'
)
def main(
    input_wildcard,
    output_path,
    ra_min,
    ra_max,
    dec_min,
    dec_max
):
    files = sorted(glob(input_wildcard))
    for f in tqdm(files):
        plot_clean_map(f, output_path, ra_min, ra_max, dec_min, dec_max)


if __name__ == '__main__':
    main()
