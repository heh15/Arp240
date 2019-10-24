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

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# basic information

ra=15*(13*u.degree+39*u.arcmin+57.7*u.arcsec)
dec=49*u.arcmin+53*u.arcsec
position=SkyCoord(dec=dec,ra=ra,frame='icrs')

beamra=204*u.degree+59*u.arcmin+41*u.arcsec
beamdec=49*u.arcmin+31*u.arcsec
beamposition=SkyCoord(dec=beamdec,ra=beamra,frame='icrs')

beammajor=2.057*u.arcsec
beamminor=1.645*u.arcsec
pa=-88.6*u.degree

rms_mom0=6.4e-4*20*math.sqrt(25)

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

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/13CO10/'
imageDir=Dir+'casa5.4/'
picDir=Dir+'picture/'
regionDir=Dir+'region/'

fitsimage=imageDir+'NGC5258_13CO10_contsub_pbcor_mom0.fits'
wcs=fits_import(fitsimage)[0]
mom0=fits_import(fitsimage)[1].data

threshold=mom0<rms_mom0
mom0[threshold]='nan'

size=u.Quantity((54,42),u.arcsec)
wcs_cut=cut_2d(mom0,position,size,wcs)[0]
mom0_cut=cut_2d(mom0,position,size,wcs)[1]

# draw the beam
beamellipse=SkyEllipticalAperture(positions=beamposition, a=beammajor, b=beamminor, theta=pa)
beamellipse_pix=beamellipse.to_pixel(wcs_cut)

fig=plt.figure()
ax=plt.subplot('111',projection=wcs_cut)
im=plt.imshow(mom0_cut,origin='lower',cmap='gist_ncar_r')
beamellipse_pix.plot(color='black')
cbar=plt.colorbar(im)
cbar.set_label('$Jy km \dot s^{-1} beam^{-1}$')
plt.savefig(picDir+'NGC5258_13CO10_2rms_mom0.png',bbox_inches='tight',pad_inches=0.2)
