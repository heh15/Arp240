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
from regions import read_ds9

############################################################
# directory

Dir='/home/heh15/workingspace/Arp240/NGC5257/SM/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
pictureDir=Dir+'picture/'
logDir=Dir+'log/'
regionDir=Dir+'region/'


############################################################
# basic settings
D=99
sr_arcsec=(180/math.pi*60**2)**2
arcsec_pc=480


############################################################
# functions

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial

    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def aperture_ring(radius_arcsec,wcs):
    a_in=radius_arcsec-1.5
    a_out=radius_arcsec
    b_out=a_out*incl
    ring_sky=SkyEllipticalAnnulus(position,a_in=a_in*u.arcsec,a_out=a_out*u.arcsec,b_out=b_out*u.arcsec,theta=PA*u.degree)
    ring_pix=ring_sky.to_pixel(wcs=wcs)
    return ring_pix

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

def mass_calc(flux_36,flux_45):
    mass=math.pow(10,5.65)*np.power(flux_36,2.85)*np.power(flux_45,-1.85)*(D/0.05)**2*0.7

    return mass

############################################################
# main program

Image=imageDir+'SPITZER_I1_39933184_0000_2_E11349850_maic.fits'
wcs=fits_import(Image)[0]
data_36=fits_import(Image)[1].data
data_36=data_36*1e6

Image=imageDir+'SPITZER_I2_39933184_0000_2_E11352627_maic.fits'
data_45=fits_import(Image)[1].data
data_45=data_45*1e6

mass_map=mass_calc(data_36,data_45)
mass_map=mass_map/sr_arcsec/arcsec_pc**2

fig=plt.figure()
plt.imshow(mass_map,origin='lower', vmax=2000, vmin=0)

# fitsimage=imageDir+'spitzer_45um_regrid.fits'
# hdr = fits.open(fitsimage)[0].header

# outfits=imageDir+'mass_map.fits'
# hdu=fits.PrimaryHDU(mass_map)
# hdu.header=hdr
# hdu.writeto(outfits, overwrite=True)

#### calculate the gas fraction in different regions. 

