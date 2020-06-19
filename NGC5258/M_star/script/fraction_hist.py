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
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import pandas as pd

##################################################
# basic settings
alpha=1.12

incl=0.45
PA=110*u.degree
ra=204*u.degree+58*u.arcmin+13.8*u.arcsec
dec=50*u.arcmin+24.5*u.arcsec
position=SkyCoord(dec=dec,ra=ra,frame='icrs')
beammaj=2.1
beammin=1.6
d=99

rms_mom0=1.6e-3*10*math.sqrt(50)

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
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

############################################################
# main program
subnum=64

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/SM/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'
picDir=Dir+'picture/'
fitsDir=Dir+'Fits/'

fitsimage=imageDir+'NGC5258_12CO10_pbcor_cube_mom0.fits'
wcs=fits_import(fitsimage)[0]
gas_mom0=fits_import(fitsimage)[1]
sd=gas_mom0/(0.0109*beammaj*beammin)*alpha
sd_binned=np.nanmean(np.nanmean(sd.reshape(subnum,5,subnum,5),axis=-1),axis=1)


fitsimage=fitsDir+'mass_map.fits'
star_mass=fits_import(fitsimage)[1]
star_mass_binned=np.nanmean(np.nanmean(star_mass.reshape(subnum,5,subnum,5),axis=-1),axis=1)


fraction=sd_binned/star_mass_binned

fraction_log=np.log10(fraction)
bins=np.linspace(-3,0,31)

fig=plt.figure()
ax=fig.add_subplot(111)
plt.hist(fraction_log.flatten(), bins=11, range=(-3,0),weights=sd_binned.flatten(), rwidth=1, density=True)
ax.text(-3.0, 1.3, 'NGC 5258', fontsize=20)
plt.xlabel('$f_{gas}$', fontsize=25)
# plt.ylabel('Mass weighted number', fontsize=20)
fig.tight_layout()
plt.savefig(picDir+'NGC5258_fraction_hist_LIRG.png')

# fig=plt.figure()
# plt.imshow(fraction, origin='lower')
