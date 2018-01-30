import click
import numpy as np
from astropy.io import fits
import astropy.units as u


import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


@click.command()
@click.argument('input_file', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('source', type=click.STRING)
def main(
    input_file,
    source,
):
    difmap_data = fits.open(input_file)
    clean_map = difmap_data['PRIMARY'].data[0][0]
    x_n_pixel = difmap_data['PRIMARY'].header['NAXIS1']
    x_ref_pixel = difmap_data['PRIMARY'].header['CRPIX1']
    x_inc = (difmap_data['PRIMARY'].header['CDELT1'] * u.degree).to(u.mas)
    x_ref_value = (difmap_data['PRIMARY'].header['CRVAL1'] * u.degree).to(u.mas)
    y_n_pixel = difmap_data['PRIMARY'].header['NAXIS2']
    y_ref_pixel = difmap_data['PRIMARY'].header['CRPIX2']
    y_inc = (difmap_data['PRIMARY'].header['CDELT2'] * u.degree).to(u.mas)
    y_ref_value = (difmap_data['PRIMARY'].header['CRVAL2'] * u.degree).to(u.mas)
    flux_min = difmap_data['PRIMARY'].header['DATAMIN']
    flux_max = difmap_data['PRIMARY'].header['DATAMAX']



    X = np.linspace(x_ref_pixel * x_inc, -x_ref_pixel * x_inc, x_n_pixel)
    Y = np.linspace(-y_ref_pixel * y_inc, y_ref_pixel * y_inc, y_n_pixel)


    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.contour(X, Y, clean_map, norm=LogNorm(), levels=np.logspace(np.log10(flux_max*1e-2), np.log10(flux_max), 10))
    # plt.savefig('build/{}_clean_map.png'.format(source))
    ax.set_xlabel('Relative Right Ascension in mas')
    ax.set_ylabel('Relative Declination in mas')
    plt.show()


if __name__ == '__main__':
    main()
