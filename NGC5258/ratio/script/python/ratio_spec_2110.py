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

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# directory

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
imageDir=Dir+'image/'
picDir=Dir+'picture/'
regionDir=Dir+'region/'
ratioDir=imageDir+'ratio/2110/spec/'

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
### 12CO 1-0

fitsfile=imageDir+'12CO10/NGC5258_12CO10_uvrange_pbcor_cube_mom0.fits'
co10=fits_import(fitsfile)[1].data

hdul=fits.open(fitsfile)
hdr=hdul[0].header
hdul.close()

### 12CO 2-1
fitsfile=imageDir+'12CO21/NGC5258_12CO21_uvtaper_pbcor_cube_mom0.fits'
co21=fits_import(fitsfile)[1].data

ratio=co21/co10*112.71**2/225.46**2

### ratio fitsfile
ratiofits=ratioDir+'NGC5258_12CO21_uvtaper_pbcor_cube_ratio.fits'
hdu=fits.PrimaryHDU(ratio)
hdu.header=hdr
hdu.writeto(ratiofits, overwrite=True)

null=np.isnan(ratio)
co10[null]='nan'
co21[null]='nan'

co10fits=imageDir+'12CO10/NGC5258_12CO10_uvrange_pbcor_cube_masked.fits'
hdu=fits.PrimaryHDU(co10)
hdu.header=hdr
hdu.writeto(co10fits, overwrite=True)

co21fits=imageDir+'12CO21/NGC5258_12CO21_uvtaper_pbcor_cube_masked.fits'
hdu=fits.PrimaryHDU(co21)
hdu.header=hdr
hdu.writeto(co21fits, overwrite=True)
