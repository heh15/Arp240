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

##################################################
# basic settings
alpha=4.6

incl=0.45
PA=110*u.degree
ra=(13*u.degree+39*u.arcmin+57.133*u.arcsec)*15
dec=49*u.arcmin+43.943*u.arcsec
position=SkyCoord(dec=dec,ra=ra,frame='icrs')
beammaj=2.1
beammin=1.6
d=99

rms_mom0=3.1e-3*10*math.sqrt(50)
############################################################
# function

# import and cut the file
# Nov 16th,2018. If the data is not in the first item of the header, try importfits and exportfits to convert image. 
def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

############################################################
# main program
Dir='/1/home/heh15/workingspace/Arp240/NGC5258/12CO10/'
imageDir=Dir+'casa5.4/'
regionDir=Dir+'region/'

fitsimage=imageDir+'NGC5258_12CO10_combine_contsub_pbcor_mom0.fits'
wcs=fits_import(fitsimage)[0]
data_masked=fits_import(fitsimage)[1]
data=data_masked.data

# regionfile=regionDir+'southarm.reg'
# southarm=read_ds9(regionfile)[0]
r=0.5/0.485*u.arcsec
southarm=SkyCircularAperture(positions=position,r=r)
southarm_pix=southarm.to_pixel(wcs)
southarm_masked=Apmask_convert(southarm_pix,data)

intensity_southarm=np.ma.mean(southarm_masked)*incl

regionfile=regionDir+'whole.reg'
whole_pix=read_ds9(regionfile)[0]
# whole_pix=whole.to_pixel(wcs)

fig=plt.figure()
ax=plt.subplot('111',projection=wcs)
plt.imshow(data,origin='lower')
southarm_pix.plot()

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

whole_masked=Apmask_convert(whole_pix,data)
# major=53.3*u.arcsec;minor=major*incl;PA=95*u.degree
# r25sky=SkyEllipticalAperture(positions=position,a=major,b=minor,theta=PA)
# r25pix=r25sky.to_pixel(wcs)
# r25_masked=Apmask_convert(r25pix,data)

flux_whole=np.ma.sum(whole_masked)
intensity_whole=flux_whole/(math.pi*49.8**2/0.3**2)
flux_addnoise=flux_whole+(49.8**2*incl/(0.3**2)-np.ma.count(flux_whole))*rms_mom0
intensity_addnoise=flux_addnoise/(math.pi*53.3**2/0.3**2)

ratio_upper=intensity_southarm/intensity_whole
ratio_lower=intensity_southarm/intensity_addnoise

# fig=plt.figure()
# ax=plt.subplot('111',projection=wcs)
# plt.imshow(data,origin='lower')
# r25pix.plot()

flux=flux_whole/(1.1331*beammaj*beammin)*0.3**2
L=2453*flux*d**2
mass=alpha*L

#### center surface density

ra=15*(13*u.degree+39*u.arcmin+57.7*u.arcsec)
dec=49*u.arcmin+51.1*u.arcsec
position=SkyCoord(dec=dec,ra=ra,frame='icrs')
r=0.5/0.485*u.arcsec
center=SkyCircularAperture(positions=position,r=r)
center_pix=center.to_pixel(wcs)

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked


center_masked=Apmask_convert(center_pix,data)

intensity_center=np.ma.mean(center_masked)*incl
ratio_upper=intensity_center/intensity_whole
ratio_lower=intensity_center/intensity_addnoise
