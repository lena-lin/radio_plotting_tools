from datetime import datetime, timedelta
from astropy.io import fits
from glob import glob
import radio_plotting_tools.spectral_map as sm
import matplotlib.pyplot as plt
import astropy.units as u
from tqdm import tqdm
import numpy as np


files_15 = sorted(glob('/home/lena/Documents/Radio/NGC1275/VLBA_15GHz/Data/FITS/201*/*.fits'))
files_43 = sorted(glob('/home/lena/Documents/Radio/NGC1275/VLBA-BU/Data_Web/FITS/20*-*/0316*.IMAP'))


epoch_partners = []

for file15 in files_15:
    date15 = datetime.strptime(fits.open(file15)[0].header['DATE-OBS'].strip(), '%Y-%m-%d')
    min_delta = timedelta(1000)
    fitting_epoch = None
    for file43 in files_43:
        date43 = datetime.strptime(fits.open(file43)[0].header['DATE-OBS'].strip(), '%Y-%m-%d')
        if abs(date15 - date43) < min_delta:
            min_delta = abs(date15 - date43)
            fitting_epoch = file43
    if min_delta <= timedelta(20):
        epoch_partners.append([file15, fitting_epoch, min_delta])


for f15, f43, td in tqdm(epoch_partners):
    header1 = fits.open(f15)[0].header
    data1 = fits.open(f15)[0].data[0][0]
    header2 = fits.open(f43)[0].header
    data2 = fits.open(f43)[0].data[0][0]

    x1, y1, ip_15 = sm.interpolated_map(header1, data1, offset=40)
    x2, y2, ip_43 = sm.interpolated_map(header2, data2, offset=200, convolution_beam_header=header1)

    header1 = fits.open(f15)[0].header
    header2 = fits.open(f43)[0].header

    ip_15[ip_15 < 1e-3] = np.nan
    ip_43[ip_43 < 1e-3] = np.nan

    freq1 = header1['CRVAL3']
    freq2 = header2['CRVAL3']

    spec_ind = sm.spectral_indices(header1, data1, header2, data2, offset1=40, offset2=200)

    fig = plt.figure(figsize=(18, 6))

    ax1 = plt.subplot(131)
    ax1.pcolorfast(-x1, y1, ip_15)
    ax1.invert_xaxis()
    ax1.set_aspect(1)
    ax1.set_xlabel('Relative RA in mas')
    ax1.set_ylabel('Relative DEC in mas')
    ax1.set_title('Clean Map {: .2f} GHz'.format((freq1 * u.Hz).to(u.GHz)))

    ax2 = plt.subplot(132)
    ax2.pcolorfast(-x2, y2, ip_43)
    ax2.invert_xaxis()
    ax2.set_aspect(1)
    ax2.set_xlabel('Relative RA in mas')
    ax2.set_ylabel('Relative DEC in mas')
    ax2.set_title('Clean Map {: .2f} GHz'.format((freq2 * u.Hz).to(u.GHz)))

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
    plt.savefig('/home/lena/Documents/Radio/NGC1275/VLBA-BU/Plots/Spectral_maps/{}_spectral_map_15_43_conv.png'.format(header1['DATE-OBS']))
    plt.close()
