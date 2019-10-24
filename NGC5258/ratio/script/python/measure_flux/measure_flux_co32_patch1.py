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
from regions import read_ds9

Dir='/home/heh15/workingspace/Arp240/NGC5258/ratio/'
scriptDir=Dir+'script/measure_1213/'
workDir=Dir+'mask/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
ratioDir=Dir+'1213/contsub/'
maskDir=Dir+'mask/'
regionDir=Dir+'region/'

############################################################
# basic settings

beam_area=3.789*2.989*1.1331

############################################################
# function
def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def masked_convert(data_masked,region_masked):
    data_mask=data_masked.mask
    region_mask=np.ma.make_mask(region_masked==0)
    region_mask=np.ma.mask_or(data_mask,region_mask)
    data_region=np.ma.masked_where(region_mask,data_masked)
    return data_region

def flux_mask_get(data_region,rms,chans,chan_width):
    flux=np.ma.sum(data_region)/beam_area_pix
    chans_tmp=chans+np.zeros((np.shape(data_region)[0],np.shape(data_region)[1]))
    error=np.sqrt(chans_tmp)*rms*chan_width/math.sqrt(beam_area_pix)
    error_masked=np.ma.masked_where(data_region.mask,error)
    uncertainty=math.sqrt(np.ma.sum(np.power(error_masked,2)))
    return flux, uncertainty

def ratio_cal(flux1,flux2,uncertainty1,uncertainty2):
    ratio=flux1/flux2
    uncertainty=math.sqrt((uncertainty1/flux1)**2+(uncertainty2/flux2)**2+2*0.05**2)*ratio
    return ratio, uncertainty

def Regmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

############################################################
# main program

regions=['southarm']
values=['12CO32 intensity','12CO32 uncertainty']
# values2=['12CO32 intensity','12CO32 uncertainty', '12CO10 intensity', '12CO10 uncertainty', '12CO21 intensity', '12CO21 uncertainty', '13CO10 intensity', '13CO10 uncertainty']

# CO13=['13CO10 flux', '13CO10 uncertainty']
ratio_table=pd.DataFrame(index=values,columns=regions)

intensities=pd.DataFrame(columns=regions,index=values)
peaks=pd.DataFrame(columns=regions,index=values)

### 12CO 3-2
## south

fitsimage=imageDir+'12CO32/NGC5258co32_all_map40r_2rms_mom0.fits'
wcs=fits_import(fitsimage)[0]
data_masked=fits_import(fitsimage)[1]

pixel_area=np.abs(wcs.wcs.cdelt[0]*wcs.wcs.cdelt[1]*3600**2)
beam_area_pix=beam_area/pixel_area

#import the south arms.
region=regionDir+'south_SFR.reg'
southarm_sky=read_ds9(region)[0]
southarm_pix=southarm_sky.to_pixel(wcs)

southarm_masked=Regmask_convert(southarm_pix, data_masked)
peaks['southarm']['12CO32 intensity']=np.nanmax(southarm_masked)
intensities['southarm']['12CO32 intensity']=np.nanmean(southarm_masked)

imagename=imageDir+'12CO32/NGC5258co32_all_map40r_nchan.fits'
chans=fits_import(imagename)[1].data
rms=0.1
chan_width=40

flux=flux_mask_get(southarm_masked, rms, chans, chan_width)[0]
uncertainty=flux_mask_get(southarm_masked, rms, chans, chan_width)[1]
peaks['southarm']['12CO32 uncertainty']= peaks['southarm']['12CO32 intensity']*uncertainty/flux
intensities['southarm']['12CO32 uncertainty']= intensities['southarm']['12CO32 intensity']*uncertainty/flux
