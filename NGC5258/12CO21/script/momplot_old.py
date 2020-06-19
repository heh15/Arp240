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

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# directories and files

Dir = '/1/home/heh15/workingspace/Arp240/NGC5257/12CO21/'
imageDir = Dir + 'casa5.4/'
picDir = Dir + 'picture/'
regionDir = Dir + 'region/'

mom0file = imageDir + 'NGC5257_12CO21_combine_pbcor_mom0.fits'
mom1file = imageDir + 'NGC5257_12CO21_pbcor_cube_mom1.fits'
mom2file = imageDir + 'NGC5257_12CO21_pbcor_cube_mom2.fits'

############################################################
# basic information

galaxy = 'NGC5257'
line = '12CO21'

ra=15*(13*u.degree+39*u.arcmin+57.7*u.arcsec)
dec=49*u.arcmin+53*u.arcsec
position=SkyCoord(dec=dec,ra=ra,frame='icrs')

beamra=204*u.degree+59*u.arcmin+41*u.arcsec
beamdec=49*u.arcmin+31*u.arcsec
beamposition=SkyCoord(dec=beamdec,ra=beamra,frame='icrs')

beammajor=1.008*u.arcsec
beamminor=0.513*u.arcsec
pa=-64.574*u.degree

############################################################
# function 

def cut_2d(data,position,size,wcs):
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

def fits_import(fitsimage, item=0):

    hdr = fits.open(fitsimage)[item].header
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

############################################################
# main program

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/12CO21/'
imageDir=Dir+'casa5.4/'
picDir=Dir+'picture/'
regionDir=Dir+'region/'

fitsimage=imageDir+'NGC5258_12CO21_combine_pbcor_mom0.fits'
wcs=fits_import(fitsimage)[0]
mom0=fits_import(fitsimage)[1]

size=u.Quantity((54,42),u.arcsec)
wcs_cut=cut_2d(mom0,position,size,wcs)[0]
mom0_cut=cut_2d(mom0,position,size,wcs)[1]

# plot the beam size
beamellipse=SkyEllipticalAperture(positions=beamposition, a=beammajor, b=beamminor, theta=pa)
beamellipse_pix=beamellipse.to_pixel(wcs_cut)

fig=plt.figure()
ax=plt.subplot('111',projection=wcs_cut)
im=plt.imshow(mom0_cut,origin='lower',cmap='gist_ncar_r')
beamellipse_pix.plot(color='black')
cbar=plt.colorbar(im)
cbar.set_label('$Jy km \dot s^{-1} beam^{-1}$')
plt.savefig(picDir+'NGC5258_12CO21_pbcor_cube_mom0.png',bbox_inches='tight',pad_inches=0.2)

### moment 1 plot

fitsimage=imageDir+'NGC5258_12CO21_pbcor_cube_mom1.fits'
wcs=fits_import(fitsimage)[0]
mom1=fits_import(fitsimage)[1]

size=u.Quantity((54,42),u.arcsec)
wcs_cut=cut_2d(mom1,position,size,wcs)[0]
mom1_cut=cut_2d(mom1,position,size,wcs)[1]

fig=plt.figure()
ax=plt.subplot('111',projection=wcs_cut)
im=plt.imshow(mom1_cut,origin='lower',cmap='rainbow')
beamellipse_pix.plot(color='black')
cbar=plt.colorbar(im)
cbar.set_label('km/s')
plt.savefig(picDir+'NGC5258_12CO21_pbcor_cube_mom1.png', bbox_inches='tight',pad_inches=0.2)


### moment 2 plot

fitsimage=imageDir+'NGC5258_12CO21_pbcor_cube_mom2.fits'
wcs=fits_import(fitsimage)[0]
mom2=fits_import(fitsimage)[1]

size=u.Quantity((54,42),u.arcsec)
wcs_cut=cut_2d(mom2,position,size,wcs)[0]
mom2_cut=cut_2d(mom2,position,size,wcs)[1]

fig=plt.figure()
ax=plt.subplot('111',projection=wcs_cut)
im=plt.imshow(mom2_cut,origin='lower',cmap='rainbow')
beamellipse_pix.plot(color='black')
cbar=plt.colorbar(im)
cbar.set_label('km/s')
plt.savefig(picDir+'NGC5258_12CO21_pbcor_cube_mom2.png',bbox_inches='tight',pad_inches=0.2)
