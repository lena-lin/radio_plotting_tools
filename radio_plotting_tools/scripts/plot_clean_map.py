import click
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.table import Table
from matplotlib.patches import Ellipse
import matplotlib as mpl


import matplotlib.pyplot as plt

from ..coordinates import get_pixel_coordinates


@click.command()
@click.argument('input_file', type=click.Path(file_okay=True, dir_okay=False))
@click.argument('output_path', type=click.Path(file_okay=False, dir_okay=True))
@click.option(
    '--ra_min',
    type=click.INT,
    help='RA min for plotting',
    default='180'
)
@click.option(
    '--ra_max',
    type=click.INT,
    help='RA max for plotting',
    default='330'
)
@click.option(
    '--dec_min',
    type=click.INT,
    help='DEC min for plotting',
    default='150'
)
@click.option(
    '--dec_max',
    type=click.INT,
    help='DEC max for plotting',
    default='300'
)
def main(
    input_file,
    output_path,
    ra_min,
    ra_max,
    dec_min,
    dec_max
):
    difmap_data = fits.open(input_file)
    bdsm_list = fits.open('/home/lena/Dokumente/Radio/NGC1275/VLBA-BU/Data/2015-02/2015-02_gauss_pybdsf', ignore_missing_end=True)
    t = Table(bdsm_list[1].data)
    clean_map = difmap_data['PRIMARY'].data[0][0]

    date = difmap_data['PRIMARY'].header['DATE-OBS']
    source = difmap_data['PRIMARY'].header['OBJECT']
    x, y = get_pixel_coordinates(difmap_data['PRIMARY'].header)
    X, Y = np.meshgrid(x.value[ra_min:ra_max], y.value[dec_min:dec_max])
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

    ra_rel = ((t['RA']-t['RA'].mean()) * u.deg).to(u.mas).value
    dec_rel = ((t['DEC']-t['DEC'].mean()) * u.deg).to(u.mas).value

    ra_cor = ((t['RA']-t['RA'][4]) * u.deg).to(u.mas).value
    dec_cor = ((t['DEC']-t['DEC'][4]) * u.deg).to(u.mas).value

    ra_rel = ra_cor
    dec_rel = dec_cor

    cmap = plt.cm.get_cmap('viridis')
    norm = mpl.colors.Normalize(vmin=t['Source_id'].min(), vmax=t['Source_id'].max())

    ells = [Ellipse(
        (-ra_rel[i], dec_rel[i]),
        (t['Maj'][i]*u.deg).to(u.mas).value,
        (t['Min'][i]*u.deg).to(u.mas).value,
        t['PA'][i] + 90)
            for i in range(len(t))]






    fig, ax = plt.subplots(figsize=(10, 10))
    ax.pcolormesh(
                    X,
                    Y,
                    clean_map[dec_min:dec_max, ra_min:ra_max],
                    cmap='magma'
                )
    cax = ax.contour(
                    X,
                    Y,
                    clean_map[dec_min:dec_max, ra_min:ra_max],
                    origin='lower',
                    levels=np.logspace(
                       np.log10(rms*n_sigma),
                       np.log10(clean_map.max()),
                       n_contours
                    ),
                    cmap='Pastel2'
                    )
    ax.set_xlabel('Relative Right Ascension in mas')
    ax.set_ylabel('Relative Declination in mas')
    ax.text(
        -2.7,
        1.1,
        date,
        color='lightyellow',
        fontsize=16)
    ax.text(
        2,
        1.1,
        source,
        color='lightyellow',
        fontsize=16)

    for e, c in zip(ells, t['Source_id']):
        ax.add_artist(e)
        e.set_edgecolor(cmap(norm(c))[:3])
        e.set_alpha(0.7)
        print(e)

    cbar = fig.colorbar(cax, fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel('Flux in Jy/beam', fontsize=16)
    ax.set_xlabel('Relative Right Ascension in mas', fontsize=16)
    ax.set_ylabel('Relative Declination in mas', fontsize=16)
    ax.set_aspect(1)
    plt.tight_layout()
    plt.show()
    # plt.savefig('{}/{}_{}_cleanmap.pdf'.format(output_path, date, source))
    # plt.savefig('{}/{}_{}_cleanmap.png'.format(output_path, date, source))


if __name__ == '__main__':
    main()
