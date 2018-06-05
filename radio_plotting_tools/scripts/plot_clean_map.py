import click
import numpy as np
from astropy.io import fits
import astropy.units as u


import matplotlib.pyplot as plt

from ..coordinates import get_pixel_coordinates


@click.command()
@click.argument('input_file', type=click.Path(file_okay=True, dir_okay=False))
def main(
    input_file,
):
    difmap_data = fits.open(input_file)
    clean_map = difmap_data['PRIMARY'].data[0][0]

    date = difmap_data['PRIMARY'].header['DATE-OBS']
    source = difmap_data['PRIMARY'].header['OBJECT']
    x, y = get_pixel_coordinates(difmap_data['PRIMARY'].header)

    window_width = 100
    y_max = clean_map.shape[0]
    x_max = clean_map.shape[1]
    x_1 = 0 + window_width
    x_2 = x_max - window_width
    y_1 = 0 + window_width
    y_2 = y_max - window_width

    n_sigma = 3
    n_contours = 10

    rms = np.array(
           [clean_map[0:y_max, 0:x_1].std(),
            clean_map[0:y_max, x_2:x_max].std(),
            clean_map[0:y_1, 0:x_max].std(),
            clean_map[y_2:y_max, 0:x_max].std()]
        ).mean()
    print(rms)
    fig, ax = plt.subplots(figsize=(10, 10))
    cax = ax.imshow(clean_map[150:300, 180:330],
                    origin='lower'
                    )
    ax.contour(clean_map[150:300, 180:330],
               origin='lower',
               levels=np.logspace(
                   np.log10(rms*n_sigma),
                   np.log10(clean_map.max()),
                   n_contours
               ),
               cmap='tab20c'
               )
    ax.set_xlabel('Relative Right Ascension in mas')
    ax.set_ylabel('Relative Declination in mas')
    ax.text(
        5,
        140,
        date,
        color='lightblue',
        fontsize=16)
    ax.text(
        120,
        140,
        source,
        color='lightblue',
        fontsize=16)

    cbar = fig.colorbar(cax, fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel('Flux in Jy/beam')
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
