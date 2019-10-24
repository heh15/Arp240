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

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# basic information

ra=15*(13*u.degree+39*u.arcmin+57.7*u.arcsec)
dec=49*u.arcmin+53*u.arcsec
position=SkyCoord(dec=dec,ra=ra,frame='icrs')

############################################################
# function 

def cut_2d(data,position,size,wcs):
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

############################################################
# main program

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/12CO10/'
imageDir=Dir+'casa5.4/'
picDir=Dir+'picture/'
regionDir=Dir+'region/'

name=imageDir+'NGC5258_12CO10_combine_contsub_pbcor.fits'
imagecube=SpectralCube.read(name)
# imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=112.73*10**9*u.Hz)
common_beam = imagecube.beams.common_beam(tolerance=1e-5)
Imcube = imagecube.convolve_to(common_beam)

## create rms cube
rmscube=cube.calc_noise_in_cube(Imcube)

# # mask the the low value of rmscube. 
# mask=rmscube<3.0e-3*u.Jy/u.beam
# lowrms=rmscube.with_mask(~mask)
# newrms=lowrms.with_fill_value(3.0e-3)

## find the signal of the cube. 
outcube=cube.find_signal_in_cube(Imcube,rmscube,snr_hi=5)
outcube.write(imageDir+'NGC5258_12CO10_pbcor_cube.fits',overwrite=True)
# kcube=outcube.to(u.K)
# kcube.write(imageDir+'NGC5258_12CO10_kcube.fits')

fitsfile=imageDir+'NGC5258_12CO10_pbcor_cube.fits'
outcube=SpectralCube.read(fitsfile)

### moment 0 map ###
imcube=outcube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=112.73*10**9*u.Hz)
moment0=imcube.moment(order=0)
moment0.write(imageDir+'NGC5258_12CO10_pbcor_cube_mom0.fits',overwrite=True)
mom0_data=moment0.hdu.data

wcs=WCS(moment0.hdu.header)
size=u.Quantity((54,42),u.arcsec)
wcs_cut=cut_2d(mom0_data,position,size,wcs)[0]
mom0_cut=cut_2d(mom0_data,position,size,wcs)[1]

### moment 1 map ### 
imcube=outcube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=112.73*10**9*u.Hz)
moment1=imcube.moment(order=1)
moment1.write(imageDir+'NGC5258_12CO10_pbcor_cube_mom1.fits',overwrite=True)

wcs=WCS(moment1.hdu.header)
mom1_data=moment1.hdu.data

# size=u.Quantity((54,42),u.arcsec)
# wcs_cut=cut_2d(mom1_data,position,size,wcs)[0]
# mom1_cut=cut_2d(mom1_data,position,size,wcs)[1]

fig=plt.figure()
ax=plt.subplot('111',projection=wcs_cut)
im=plt.imshow(mom0_cut,origin='lower',cmap='rainbow')
cbar=plt.colorbar(im)
# cbar.set_label('km/s')
