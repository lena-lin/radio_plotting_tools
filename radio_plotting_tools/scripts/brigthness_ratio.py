import numpy as np
import matplotlib as mpl
mpl.use('pgf')
mpl.rcParams.update(
{'font.size': 20,
'font.family': 'serif',
'text.usetex': True,
'pgf.rcfonts': False,
'pgf.texsystem': 'lualatex',
'text.latex.unicode': True,
})
import matplotlib.pyplot as plt

import pandas as pd
import glob
from astropy.io import fits

df_components = pd.read_csv('components.csv')
datapaths = ([])
for data in glob.glob('../../../../radiokram/modeling/1638+118/*/final_mod.fits'):
    datapaths = sorted(np.append(datapaths, data))
dates = sorted([fits.open(file)['PRIMARY'].header['DATE-OBS'] for file in datapaths])

# Plot major_axis vs distance
fig = plt.figure(figsize=(12,20))
ax = fig.add_subplot(111, aspect='equal')

r = np.sqrt(df_components.loc[df_components['c_i'] == 4, 'x_positions']**2 + df_components.loc[df_components['c_i'] == 4, 'y_positions']**2)
print(df_components.loc[df_components['c_i'] == 4, 'y_shifted'])

plt.plot(r,
         df_components.loc[df_components['c_i'] == 4, 'major_axes'],
         linestyle='none',
         marker='.',
         markersize=8
         )


plt.ylabel('Majoraxis / mas')
plt.xlabel('Distance / minor_axes')
#plt.ylim(0,18)


plt.savefig('brightness_ratio.pdf', bbox_inches='tight', pad_inches=0.1)
