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
from regions import read_ds9
import cube

############################################################
# directory

Dir='/home/heh15/workingspace/Arp240/NGC5257/ratio/'
scriptDir=Dir+'script/measure_1213/'
workDir=Dir+'mask/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
ratioDir=Dir+'1213/contsub/'
maskDir=Dir+'mask/'
regionDir=Dir+'region/'

############################################################
# basic settings
beam_area=1.019*0.522*1.1331
beam_area_pix=beam_area/0.01
regions=['center','hinge','south','rest','all']
values=['FWMH']
linewidth=pd.DataFrame(index=regions,columns=values)
rms=0.003*u.Jy/u.beam

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

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

def Regmask3d(data,region_pix,lowchan,highchan):
    region_mask=region_pix.to_mask()
    shape=np.shape(data)
    mask=region_mask.to_image(shape=((shape[1],shape[2])))
    mask3d=np.zeros((shape[0],shape[1],shape[2]))
    mask3d[lowchan:highchan]=mask
    maskTF=mask3d==1

    data_masked=np.copy(data)
    data_masked[maskTF]='nan'

    return data_masked

############################################################
# main program

# fitsimage=imageDir+'12CO21/NGC5257_12CO21_combine.fits'
# imagecube=SpectralCube.read(fitsimage)
# imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)
# common_beam=imcube.beams.common_beam()
# imcube2=imcube.convolve_to(common_beam)
# threshold=4*rms
# mask=imcube2>threshold
# Imcube=imcube2.with_mask(mask)
# Imcube.write(imageDir+'12CO21/NGC5257_12CO21_combine_signal.fits', overwrite=True)

## smooth the image first and then pick the signal. 
# fitsimage=imageDir+'12CO21/NGC5257_12CO21_combine_contsub_uvtaper.fits'
# imagecube=SpectralCube.read(fitsimage)
# imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)

# beam=Beam(major=3.789*u.arcsec, minor=2.989*u.arcsec, pa=-17.357*u.deg)
# imcube_smooth=imcube.convolve_to(beam)
# rmscube=cube.calc_noise_in_cube(imcube_smooth, spatial_average_nbeam=1.0, spectral_average_nchan=2)
# Imcube_smooth=cube.find_signal_in_cube(imcube_smooth,rmscube,snr_hi=5)

# Imcube_smooth.write(imageDir+'12CO21/NGC5257_12CO21_combine_signal_smooth_co32.fits', overwrite=True)

fitsimage=imageDir+'12CO21/NGC5257_12CO21_combine_signal_smooth_co32.fits'
Imcube=SpectralCube.read(fitsimage)
wcs=WCS(Imcube.hdu.header)

'''Try different regions'''
## center
# mask=maskDir+'center_SFR_regrid_mask.fits'
mask=maskDir+'center_SFR_mask.fits'
region_masked=fits.open(mask)[0].data[0][0]
region_mask=np.ma.make_mask(region_masked==1)
center=Imcube.with_mask(region_mask)
# spectrum=center.mean(axis=(1,2))

fwhm_map=center.linewidth_fwhm()
array=np.array(fwhm_map)
FWHM=np.nanmedian(array)
linewidth['FWMH']['center']=FWHM

fwhm_map=center.linewidth_sigma()
array=np.array(fwhm_map)
sigma=np.nanmedian(array)

width=Imcube.linewidth_fwhm()[160,160]
print('width')

## south arms
# mask=maskDir+'hinge_SFR_regrid_mask.fits'
mask=maskDir+'hinge_SFR_mask.fits'
region_masked=fits.open(mask)[0].data[0][0]
region_mask=np.ma.make_mask(region_masked==1)
region=Imcube.with_mask(region_mask)
# spectrum=region.mean(axis=(1,2))

fwhm_map=region.linewidth_fwhm()
array=np.array(fwhm_map)
FWHM=np.nanmedian(array)
linewidth['FWMH']['hinge']=FWHM


## south source
# mask=maskDir+'south_SFR_regrid_mask.fits'
mask=maskDir+'south_SFR_mask.fits'
region_masked=fits.open(mask)[0].data[0][0]
region_mask=np.ma.make_mask(region_masked==1)
region=Imcube.with_mask(region_mask)
# spectrum=region.mean(axis=(1,2))

fwhm_map=region.linewidth_fwhm()
array=np.array(fwhm_map)
FWHM=np.nanmedian(array)
linewidth['FWMH']['south']=FWHM


# ## rest disk
# # mask=maskDir+'rest_regrid_mask.fits'
# mask=maskDir+'rest_SFR_mask.fits'
# region_masked=fits.open(mask)[0].data[0][0]
# region_mask=np.ma.make_mask(region_masked==1)
# region=Imcube.with_mask(region_mask)
# # spectrum=region.mean(axis=(1,2))

# fwhm_map=region.linewidth_fwhm()
# array=np.array(fwhm_map)
# FWHM=np.nanmedian(array)
# linewidth['FWMH']['rest']=FWHM

# ## all the regions
# mask=maskDir+'whole_init_mask.fits'
# region_masked=fits.open(mask)[0].data[0][0]
# region_mask=np.ma.make_mask(region_masked==1)
# region=Imcube.with_mask(region_mask)
# # spectrum=region.mean(axis=(1,2))

# fwhm_map=region.linewidth_fwhm()
# array=np.array(fwhm_map)
# FWHM=np.nanmedian(array)
# linewidth['FWMH']['all']=FWHM

## co 3-2 concentration. 
file=regionDir+'co32_concentration.reg'
co32_sky=read_ds9(file,errors='warn')[0]
co32_pix=co32_sky.to_pixel(wcs)

data=np.array(Imcube)
highchan=np.shape(data)[0]-1;lowchan=0
co32_masked=Regmask3d(data,co32_pix,lowchan,highchan)
co32_mask=np.isnan(co32_masked)

region=Imcube.with_mask(co32_mask)
fwhm_map=region.linewidth_fwhm()
fwhm_map.quicklook()

array=np.array(fwhm_map)
FWHM=np.nanmean(array)
