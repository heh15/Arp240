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
import radio_beam

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# directory

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/uvtaper/'
imageDir=Dir+'image/'
picDir=Dir+'picture/'
regionDir=Dir+'region/'

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

### 12CO 1-0 cube processing. ### 
name=imageDir+'12CO10/NGC5257_12CO10_combine_contsub_uvrange_pbcor.fits'
imagecube=SpectralCube.read(name)
common_beam = imagecube.beams.common_beam(tolerance=1e-5)
# Imcube = imagecube.convolve_to(common_beam)
smooth_beam=radio_beam.Beam(major=2.186*u.arcsec, minor=1.896*u.arcsec, pa=-87.6*u.deg)
Imcube = imagecube.convolve_to(smooth_beam)

## create rms cube
rmscube=cube.calc_noise_in_cube(Imcube, spatial_average_nbeam=1.0, spectral_average_nchan=2)

## find the signal of the cube. 
outcube=cube.find_signal_in_cube(Imcube,rmscube,snr_hi=4,nchan_hi=1,nchan_lo=2)
outcube.write(imageDir+'12CO10/NGC5257_12CO10_uvrange_pbcor_cube.fits',overwrite=True)

imcube=outcube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=112.73*10**9*u.Hz)
mom0=imcube.moment(order=0)
mom0.write(imageDir+'12CO10/NGC5257_12CO10_uvrange_pbcor_cube_mom0.fits',overwrite=True)

co10mom0=mom0

### 12CO 2-1 cube processing. ###
name=imageDir+'12CO21/NGC5257_12CO21_combine_contsub_uvtaper_pbcor.fits'
imagecube=SpectralCube.read(name)
common_beam = imagecube.beams.common_beam(tolerance=1e-5)

smooth_beam=radio_beam.Beam(major=2.186*u.arcsec, minor=1.896*u.arcsec, pa=-87.6*u.deg)
Imcube = imagecube.convolve_to(smooth_beam)

## create rms cube
rmscube=cube.calc_noise_in_cube(Imcube, spatial_average_nbeam=1.0, spectral_average_nchan=2)

## find the signal of the cube. 
outcube=cube.find_signal_in_cube(Imcube,rmscube,snr_hi=4,nchan_hi=1,nchan_lo=2)
outcube.write(imageDir+'12CO21/NGC5257_12CO21_uvtaper_pbcor_cube.fits',overwrite=True)

imcube=outcube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)
mom0=imcube.moment(order=0)
mom0.write(imageDir+'12CO21/NGC5257_12CO21_uvtaper_pbcor_cube_mom0.fits',overwrite=True)

co21mom0=mom0

ratio=co21mom0/co10mom0*112.73**2/225.46**2


ratio.write(imageDir+'NGC5257_2110_ratio.fits', overwrite=True)
