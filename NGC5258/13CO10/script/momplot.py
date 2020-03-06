# import cube
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
rcParams['mathtext.default']='regular'

############################################################
# directories and files

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/13CO10/'
imageDir=Dir+'casa5.4/'
picDir=Dir+'picture/'
regionDir=Dir+'region/'

mom0file=imageDir+'NGC5258_13CO10_contsub_pbcor_mom0.fits'

############################################################
# basic information
galaxy='NGC5258'
line='13CO10'

ra=15*(13*u.degree+39*u.arcmin+57.7*u.arcsec)
dec=49*u.arcmin+53*u.arcsec
position=SkyCoord(dec=dec,ra=ra,frame='icrs')

beamra=204*u.degree+59*u.arcmin+41*u.arcsec
beamdec=49*u.arcmin+31*u.arcsec
beamposition=SkyCoord(dec=beamdec,ra=beamra,frame='icrs')

beammajor=2.057*u.arcsec/2.0
beamminor=1.645*u.arcsec/2.0
pa=-88.6*u.degree

rms_mom0=6.4e-4*20*math.sqrt(25)

############################################################
# function 

def fits_import(fitsimage, item=0):
    hdr = fits.open(fitsimage)[item].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[item].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def cut_2d(data,position,size,wcs):
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

############################################################
# basic settings

testra = 204.99581309
testdec = 0.83741222

linelabel = '$^{13}$CO 1-0'

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

fig = plt.figure()
ax = plt.subplot('111', projection=wcs_cut)

x, y = wcs_cut.wcs_world2pix(testra, testdec, 1)
props = dict(facecolor='white', alpha=0.5, edgecolor='None')
plt.text(x, y, galaxy + ' ' + linelabel, fontsize=15, bbox=props)

ax.tick_params(labelsize=8, direction='in')
im = plt.imshow(mom0_cut,
########## private region
#                 norm=colors.PowerNorm(gamma=0.7),
#                 vmax=15,
####################
                vmin=0, 
                origin='lower',
                cmap='gist_ncar_r')

beamellipse_pix.plot(color='black')
cbar = plt.colorbar(im)
cbar.set_label('$Jy\ beam^{-1}\ km\ s^{-1} $', fontsize=15)
plt.savefig(picDir + galaxy + '_' + line + '_pbcor_cube_mom0.png',
            bbox_inches='tight',
            pad_inches=0.2)

