from astropy.coordinates import Angle, SkyOffsetFrame, ICRS
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from skimage.measure import block_reduce

import img_scale as IS

plt.style.use('paper')

# grab WCS
hdr = fits.getheader('Downloads/SAGE_LMC_IRAC8.0_2_mosaic.fits')
hdr['CD1_1'] *= 2
hdr['CD1_2'] *= 2
hdr['CD2_1'] *= 2
hdr['CD2_2'] *= 2
hdr['CRPIX1'] /= 2
hdr['CRPIX2'] /= 2
hdr['PIXSCAL1'] *= 2
hdr['PIXSCAL2'] *= 2
wcs = WCS(hdr)

center = ICRS(Angle(str(hdr['CRVAL1'])+'deg'), Angle(str(hdr['CRVAL2'])+'deg'))

# grab image data
data = list()
for wl in ['8.0', '4.5', '3.6']:
    data.append(fits.getdata('Downloads/SAGE_LMC_IRAC'+wl+'_2_mosaic.fits'))
data = np.array(data)
data = np.swapaxes(data, 0, 2)
data = np.swapaxes(data, 0, 1)

data = block_reduce(data, block_size=(2, 2, 1))

for i in range(3):
    tmp = np.nanmin(data[:, :, i])
    data[:, :, i] += np.abs(tmp)
for i in range(3):
    tmp = np.nanmax(data[:, :, i])
    data[:, :, i] *= 1.0/tmp

data[:, :, 0] = IS.linear(data[:, :, 0], scale_min=0.0027, scale_max=0.006)
data[:, :, 1] = IS.linear(data[:, :, 1], scale_min=0.00025, scale_max=0.0016)
data[:, :, 2] = IS.linear(data[:, :, 2], scale_min=0.0005, scale_max=0.0021)

'''
fig = plt.figure(dpi=100)
ax = fig.add_subplot(111)
col = ['r', 'g', 'b']
for i in range(3):
    tmp = data[:, :, i]
    ax.hist(tmp[~np.isnan(tmp)], bins='auto', color=col[i], histtype='step')
#ax.set_xscale('log')
#ax.set_yscale('log')
plt.show()
'''
fig = plt.figure(figsize=(4.983, 4.983))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection=wcs)

ax.imshow(data)

# add field annotations
fieldPos = [[Angle('05h38m48s'), Angle('-69d04m48s')],
            [Angle('05h39m37s'), Angle('-69d45m48s')],
            [Angle('05h40m09s'), Angle('-69d44m44s')],
            [Angle('05h13m18s'), Angle('-69d22m25s')],
            [Angle('05h47m09s'), Angle('-70d40m16s')],
            [Angle('05h44m29s'), Angle('-69d25m43s')],
            [Angle('05h24m09s'), Angle('-71d53m37s')]]
fieldNames = ['30 Dor-10', 'N159W', 'N159E', 'N113', 'GMC 225', 'N166',
              'PCC 11546']
xytext = [(-89, 25), (-75, -18), (-72, 10), (10, -25), (10, -25),
          (-60, 15), (-57, 20)]
annConns = ['angle,angleA=0,angleB=90',
            'angle,angleA=0,angleB=270',
            'angle,angleA=0,angleB=90',
            'angle,angleA=0,angleB=270',
            'angle,angleA=0,angleB=270',
            'angle,angleA=0,angleB=90',
            'angle,angleA=0,angleB=90',]
c = (200.0/256, 198.0/256, 54.0/256)
for i in range(len(fieldPos)):
    pix = wcs.all_world2pix(fieldPos[i][0].deg, fieldPos[i][1].deg, 1)
    ax.annotate(fieldNames[i], xy=pix,
                xytext=xytext[i], textcoords='offset pixels', color=c,
                arrowprops=dict(arrowstyle='-', color=c,
                                connectionstyle=annConns[i]),
                usetex=False, family='sans-serif')
    ax.plot(fieldPos[i][0].deg, fieldPos[i][1].deg, color=c, marker='o',
            fillstyle='none', markersize=4, transform=ax.get_transform('world'))

# add N/E direction arrows
ax.arrow(69.5, -69.85, 2, 0, color=c, head_width=0.15, head_length=0.27,
         transform=ax.get_transform('world'))
ax.text(72.2, -69.95, 'E', color=c, transform=ax.get_transform('world'),
        usetex=False, family='sans-serif')
ax.arrow(69.5, -69.85, 0, 0.68, color=c, head_width=0.4, head_length=0.1,
         transform=ax.get_transform('world'))
ax.text(69.4, -69.03, 'N', color=c, transform=ax.get_transform('world'),
        usetex=False, family='sans-serif')

# turn off ticks and tick labels
ax.coords[0].set_ticks_visible(False)
#ax.coords[0].set_ticklabel_visible(False)
ax.coords[1].set_ticks_visible(False)
#ax.coords[1].set_ticklabel_visible(False)

# set field of view
ax.set_xlim([250, 6700])
ax.set_ylim([1200, 6800])

ax.set_xlabel('Right Ascension (J2000)')
ax.set_ylabel('Declination (J2000)')

lon = ax.coords[0]
lat = ax.coords[1]
lon.set_major_formatter('hh:mm')

#plt.show()
fig.savefig('lmc.pdf', pad_inches=0)
plt.close()
