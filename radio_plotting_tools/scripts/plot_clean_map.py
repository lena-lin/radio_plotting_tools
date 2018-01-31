import click
import numpy as np
from astropy.io import fits


import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from ..coordinates import get_pixel_coordinates


@click.command()
@click.argument('input_file', type=click.Path(file_okay=True, dir_okay=False))
@click.option(
    '-loglevs',
    type=click.FLOAT,
    default=0.02,
)
def main(
    input_file,
    loglevs,
):
    difmap_data = fits.open(input_file)
    clean_map = difmap_data['PRIMARY'].data[0][0]

    flux_min = clean_map.min()
    flux_max = clean_map.max()
    loglevs_min = flux_min + (flux_max - flux_min) * loglevs
    x, y = get_pixel_coordinates(difmap_data['PRIMARY'].header)

    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.contour(x, y, clean_map, norm=LogNorm(), levels=np.logspace(np.log10(loglevs_min), np.log10(flux_max), 10))
    ax.set_xlabel('Relative Right Ascension in mas')
    ax.set_ylabel('Relative Declination in mas')
    plt.show()


if __name__ == '__main__':
    main()
