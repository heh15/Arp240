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
rcParams['mathtext.default']='regular'

############################################################
# basic information

position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
beamra=204*u.degree+58*u.arcmin+30*u.arcsec
beamdec=50*u.arcmin+3*u.arcsec
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

############################################################
# function.

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
# main program

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/12CO21/'
imageDir=Dir+'casa5.4/'
picDir=Dir+'picture/'
regionDir=Dir+'region/'

### moment 0 plot
fitsimage=imageDir+'NGC5257_12CO21_combine_pbcor_mom0.fits'
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
im=plt.imshow(mom0_cut,origin='lower',cmap='gist_ncar_r',norm=colors.PowerNorm(gamma=0.7), vmax=15)
beamellipse_pix.plot(color='black')
cbar=plt.colorbar(im)
cbar.set_label('$Jy km \dot s^{-1} beam^{-1}$')
plt.savefig(picDir+'NGC5257_12CO21_pbcor_cube_mom0.png', bbox_inches='tight',pad_inches=0.2)

### moment 1 plot

fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube_mom1.fits'
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
plt.savefig(picDir+'NGC5257_12CO21_pbcor_cube_mom1.png', bbox_inches='tight',pad_inches=0.2)

### moment 2 plot

fitsimage=imageDir+'NGC5257_12CO21_pbcor_cube_mom2.fits'
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
plt.savefig(picDir+'NGC5257_12CO21_pbcor_cube_mom2.png', bbox_inches='tight',pad_inches=0.2)

