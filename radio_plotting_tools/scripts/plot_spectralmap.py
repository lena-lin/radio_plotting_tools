import numpy as np
from radio_plotting_tools.spectral_map import interpolated_map, spectral_indices
import click
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits


@click.command()
@click.argument('file1', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('file2', type=click.Path(file_okay=True, dir_okay=False))
def main(
    file1,
    file2,
    epsilon=1e-3,
):
    header1 = fits.open(file1)[0].header
    header2 = fits.open(file2)[0].header

    x1, y1, interp_map_1 = interpolated_map(file1)
    x2, y2, interp_map_2 = interpolated_map(file2)

    interp_map_1[interp_map_1 < epsilon] = np.nan
    interp_map_2[interp_map_2 < epsilon] = np.nan

    freq1 = header1['CRVAL3']
    freq2 = header2['CRVAL3']

    spec_ind = spectral_indices(file1, file2)

    fig = plt.figure(figsize=(18, 6))

    ax1 = plt.subplot(131)
    ax1.pcolorfast(-x1, y1, interp_map_1)
    ax1.invert_xaxis()
    ax1.set_aspect(1)
    ax1.set_xlabel('Relative RA in mas')
    ax1.set_ylabel('Relative DEC in mas')
    ax1.set_title('Clean Map 43 GHz')

    ax2 = plt.subplot(132)
    ax2.pcolorfast(-x2, y2, interp_map_2)
    ax2.invert_xaxis()
    ax2.set_aspect(1)
    ax2.set_xlabel('Relative RA in mas')
    ax2.set_ylabel('Relative DEC in mas')
    ax2.set_title('Clean Map 15 GHz')

    ax3 = plt.subplot(133)
    cax = ax3.pcolorfast(-x2, y2, spec_ind, cmap='Spectral')
    ax3.invert_xaxis()
    ax3.set_aspect(1)
    ax3.set_xlim(5, -3.8)
    ax3.text(4.5, 4, 'Source: {}'.format(header1['OBJECT']), size=10)
    ax3.text(4.5, 3.5, 'Freq.: {: .2f} GHz, Date: {}'.format((freq1 * u.Hz).to(u.GHz), header1['DATE-OBS']), size=10)
    ax3.text(4.5, 3, 'Freq.: {: .2f} GHz, Date: {}'.format((freq2 * u.Hz).to(u.GHz), header2['DATE-OBS']), size=10)
    ax3.set_xlabel('Relative RA in mas')
    ax3.set_ylabel('Relative DEC in mas')
    ax3.set_title('Spectral Map')

    plt.colorbar(cax, fraction=0.052, pad=0.05)

    plt.savefig('/home/lena/Dokumente/Radio/NGC1275/VLBA-BU/Plots/Spectral_maps/{}_spectral_map_15_43.png'.format(header1['DATE-OBS']))
    plt.close()


if __name__ == '__main__':

    main()
