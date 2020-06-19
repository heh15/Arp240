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


#### smooth the 12CO 2-1 image to 1.1 x 0.8 arcsec and make moment 2 map.
 
# fitsfile='NGC5257_12CO21_combine_pbcor.fits'
# imagecube=SpectralCube.read(fitsfile)
# imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)
# beam=Beam(major=1.1*u.arcsec, minor=0.8*u.arcsec, pa=-64.5*u.degree)
# imcube_smooth=imcube.convolve_to(beam)

# rmscube=cube.calc_noise_in_cube(imcube_smooth)
# outcube=cube.find_signal_in_cube(imcube_smooth,rmscube,snr_hi=5)
# outcube.write('NGC5257_12CO21_pbcor_smooth_cube_signal.fits')

fitsfile='NGC5257_12CO21_pbcor_smooth_cube_signal.fits'
outcube=SpectralCube.read(fitsfile)
_12CO21_mom2=outcube.linewidth_sigma()
_12CO21_mom2.write('NGC5257_12CO21_pbcor_smooth_cube_signal_mom2.fits')

#### imregrid the 33 GHz image 

imregrid(imagename='NGC5257_33GHz_pbcor_smooth_co21.image', template='NGC5257_12CO21_combine_pbcor.mom0', output='NGC5257_33GHz_pbcor_smooth_co21_regrid.image')
exportfits(imagename='NGC5257_33GHz_pbcor_smooth_co21_regrid.image', fitsimage='NGC5257_33GHz_pbcor_smooth_co21_regrid.fits')
