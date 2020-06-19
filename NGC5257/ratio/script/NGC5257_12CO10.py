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


############################################################
# basic setting

beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.09

############################################################
# function
def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data[0][0]
    position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
    size=u.Quantity((54,42),u.arcsec)
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs
    data_masked=np.ma.masked_invalid(data_cut)
    return wcs_cut, data_masked

def masked_convert(data_masked,region_masked):
    data_mask=data_masked.mask
    region_mask=np.ma.make_mask(region_masked==0)
    region_mask=np.ma.mask_or(data_mask,region_mask)
    data_region=np.ma.masked_where(region_mask,data_masked)
    return data_region

############################################################

# calculating 12CO10 flux
fitsimage='NGC5257_12CO10_combine_smooth_masked_mom0.fits'
wcs=fits_import(fitsimage)[0]
data12_masked=fits_import(fitsimage)[1]
data=data12_masked.data

maskimage='anomaly_mod_mask.fits'
anomaly_masked=fits_import(maskimage)[1]
anomaly=masked_convert(data12_masked,anomaly_masked)

maskimage='spiralarm_mask.fits'
spiral_masked=fits_import(maskimage)[1]
spiral= masked_convert(data12_masked,spiral_masked)

# draw the image of 12CO10 map.  

levels=[0,0.5]
colors='red'


fig=plt.figure()
ax=plt.subplot(projection=wcs)
im=ax.imshow(data,cmap='rainbow',origin='lower',vmax=7)
ax.contour(spiral_masked.data,levels=levels,transform=ax.get_transform(wcs),colors='red') 
ax.contour(anomaly_masked.data,levels=levels,transform=ax.get_transform(wcs),colors='black')
ax.xaxis.set_ticks_position('none')
# lon = ax.coords[0]
# lat = ax.coords[1]
# lon.set_major_formatter('hh:mm:ss')
plt.savefig('../../picture/NGC5257_12CO10_region.png')
