from radio_plotting_tools.spectral_map import interpolated_map
import click
import matplotlib.pyplot as plt
from astropy.io import fits
from radio_plotting_tools.coordinates import get_pixel_coordinates


@click.command()
@click.argument('file1', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('file2', type=click.Path(file_okay=True, dir_okay=False))
def main(
    file1,
    file2,
):
    x1, y1, interp_map_1 = interpolated_map(file1)
    x2, y2, interp_map_2 = interpolated_map(file2)

    x1_cm, y1_cm = get_pixel_coordinates(fits.open(file1)[0].header)
    x2_cm, y2_cm = get_pixel_coordinates(fits.open(file2)[0].header)




    fig = plt.figure()
    ax1 = plt.subplot(221)
    ax1.pcolorfast(x1,y1, interp_map_1)
    ax1.set_aspect(1)
    ax2 = plt.subplot(222)
    ax2.pcolorfast(x2,y2, interp_map_2)
    ax2.set_aspect(1)
    ax3 = plt.subplot(223)
    ax3.pcolorfast(x1_cm, y1_cm, fits.open(file1)[0].data[0][0])
    ax3.set_aspect(1)
    ax4 = plt.subplot(224)
    ax4.pcolorfast(x2_cm, y2_cm, fits.open(file2)[0].data[0][0])
    ax4.set_aspect(1)

    plt.show()

if __name__ == '__main__':
    main()
