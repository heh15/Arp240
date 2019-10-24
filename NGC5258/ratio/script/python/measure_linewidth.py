#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re
import itertools
from shutil import copytree
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import make_lupton_rgb
import aplpy
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from matplotlib.patches import Arrow
from photutils import SkyCircularAnnulus
from photutils import SkyCircularAperture
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import SkyEllipticalAperture
from photutils import SkyEllipticalAnnulus
from photutils import EllipticalAperture
from photutils import EllipticalAnnulus
from photutils import aperture_photometry
import numpy.ma as ma
import math
import shutil
import pandas as pd
from astropy.utils import data
from spectral_cube import SpectralCube, Projection
from radio_beam import Beam
from astropy import units as u
import cube

############################################################
# directory

Dir='/home/heh15/workingspace/Arp240/NGC5258/ratio/'
scriptDir=Dir+'script/measure_1213/'
workDir=Dir+'mask/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
ratioDir=Dir+'1213/contsub/'
maskDir=Dir+'mask/'

############################################################
# basic settings
beam_area=3.789*2.989*1.1331
beam_area_pix=beam_area/0.01
regions=['south','rest','all']
values=['FWMH']
linewidth=pd.DataFrame(index=regions,columns=values)
rms=0.0042*u.Jy/u.beam

############################################################
# function

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    wcs_cut=cut.wcs
    data_masked=np.ma.masked_invalid(data)
    return wcs, data_masked


def masked_convert(data_masked,region_masked):
    data_mask=data_masked.mask
    region_mask=np.ma.make_mask(region_masked==0)
    region_mask=np.ma.mask_or(data_mask,region_mask)
    data_region=np.ma.masked_where(region_mask,data_masked)
    return data_region


############################################################
# main program

# fitsimage=imageDir+'12CO21/NGC5258_12CO21_combine_uvtaper.fits'
# imagecube=SpectralCube.read(fitsimage)
# imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)
# common_beam=imcube.beams.common_beam()
# imcube2=imcube.convolve_to(common_beam)
# threshold=4*rms
# mask=imcube2>threshold
# Imcube=imcube2.with_mask(mask)

# fitsimage=imageDir+'12CO21/NGC5258_12CO21_combine_uvtaper_smooth_co32.fits'
# imagecube=SpectralCube.read(fitsimage)
# imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)
# rmscube=cube.calc_noise_in_cube(imcube, spatial_average_nbeam=1.0, spectral_average_nchan=2)
# Imcube=cube.find_signal_in_cube(imcube,rmscube,snr_hi=5)
# Imcube.write(imageDir+'12CO21/NGC5258_12CO21_combine_signal_smooth_co32.fits', overwrite=True)

fitsimage=imageDir+'12CO21/NGC5258_12CO21_combine_signal_smooth_co32.fits'
Imcube=SpectralCube.read(fitsimage)
wcs=WCS(Imcube.hdu.header)

'''Try different regions'''
## south arms
mask=maskDir+'south_SFR_mask.fits'
region_masked=fits.open(mask)[0].data[0][0]
region_mask=np.ma.make_mask(region_masked==1)
region=Imcube.with_mask(region_mask)
# spectrum=region.mean(axis=(1,2))

fwhm_map=region.linewidth_fwhm()
array=np.array(fwhm_map)
FWHM=np.nanmedian(array)
linewidth['FWMH']['south']=FWHM
print(linewidth)

width=fwhm_map[135,187]
print(width)



