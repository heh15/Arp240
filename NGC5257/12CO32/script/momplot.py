import cube
from astropy.utils import data
from spectral_cube import SpectralCube, Projection
from radio_beam import Beam
from astropy import units as u
import time
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from regions import read_ds9
from photutils import SkyEllipticalAperture
from photutils import EllipticalAperture
import matplotlib.colors as colors
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'
from regions import read_ds9

############################################################
# directories and files

Dir = '/1/home/heh15/workingspace/Arp240/NGC5257/12CO32/'
imageDir = Dir + 'image/'
picDir = Dir + 'picture/'
regionDir = Dir + 'region/'

mom0file = imageDir + 'NGC5257co32_all_map40r_2rms_mom0_shift.fits'

############################################################
# basic information

galaxy = 'NGC5257'
line = '12CO32'

position = SkyCoord(dec=50.4167 * u.arcmin,
                    ra=204.9706 * u.degree,
                    frame='icrs')
beamra = 204 * u.degree + 58 * u.arcmin + 30 * u.arcsec
beamdec = 50 * u.arcmin + 3 * u.arcsec
beamposition = SkyCoord(dec=beamdec, ra=beamra, frame='icrs')

beammajor = 3.789 * u.arcsec/2.0
beamminor = 2.989 * u.arcsec/2.0
pa = -17.357 * u.degree

############################################################
# basic settings

testra = 204.97609228
testdec = 0.84611111

linelabel = '$^{12}$CO 3-2'
vmax = 120
vmin = 0

############################################################
# function

def fits_import(fitsimage, item=0):
    hdr = fits.open(fitsimage)[item].header
    wcs = WCS(hdr).celestial
    data = fits.open(fitsimage)[item].data
    data = np.squeeze(data)
    data_masked = np.ma.masked_invalid(data)

    return wcs, data_masked


def cut_2d(data, position, size, wcs):
    cut = Cutout2D(data=data, position=position, size=size, wcs=wcs)
    data_cut = cut.data
    wcs_cut = cut.wcs

    return wcs_cut, data_cut


############################################################
# main program

### moment 0 plot
fitsimage = mom0file
wcs = fits_import(fitsimage)[0]
mom0 = fits_import(fitsimage)[1]

size = u.Quantity((54, 42), u.arcsec)
wcs_cut = cut_2d(mom0, position, size, wcs)[0]
mom0_cut = cut_2d(mom0, position, size, wcs)[1]

# plot the beam size
beamellipse = SkyEllipticalAperture(positions=beamposition,
                                    a=beammajor,
                                    b=beamminor,
                                    theta=pa)
beamellipse_pix = beamellipse.to_pixel(wcs_cut)

try:
    gama
except:
    gama = 1.0

fig = plt.figure()
ax = plt.subplot('111', projection=wcs_cut)

x, y = wcs_cut.wcs_world2pix(testra, testdec, 1)
props = dict(facecolor='white', alpha=0.5, edgecolor='None')
plt.text(x, y, galaxy + ' ' + linelabel, fontsize=15, bbox=props)

ax.tick_params(labelsize=8, direction='in')
    
im = plt.imshow(mom0_cut,
                norm=colors.PowerNorm(gamma=gama),
		vmax=vmax,
                vmin=vmin, 
                origin='lower',
                cmap='gist_ncar_r')
plt.xlabel('J2000 Right Ascension')
plt.ylabel('J2000 Declination')


beamellipse_pix.plot(color='black')

center_sky = read_ds9(regionDir+'center_SFR_co32.reg')[0]
center_pix = center_sky.to_pixel(wcs_cut)
center_pix.plot(color='black', linewidth=2.0)

co32_sky = read_ds9(regionDir+'co32_concentration.reg')[0]
co32_pix = co32_sky.to_pixel(wcs_cut)
co32_pix.plot(color='black', linewidth=2.0)

cbar = plt.colorbar(im)
cbar.set_label('$Jy\ beam^{-1}\ km\ s^{-1} $', fontsize=15)
plt.savefig(picDir + galaxy + '_' + line + '_2rms_mom0.png',
            bbox_inches='tight')

