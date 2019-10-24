import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re
import itertools
from shutil import copytree
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import make_lupton_rgb
# import aplpy
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.wcs import WCS
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

############################################################
def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def TIR(f24, f70, pixelsize=2.45): 
    flux=f24*10**(-17)/(4.25*10**10)*pixelsize**2
    L24=4*math.pi*(100*10**6*3.086*10**18)**2*flux
    nuL24=L24*1.2*10**13
    flux_erg=f70*10**(-23)
    L70=4*math.pi*(100*10**6*3.086*10**18)**2*flux_erg
    nuL70=4.283*10**12*L70
    LTIR=3.98*nuL24+1.553*nuL70

    return LTIR

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

############################################################
Dir='/home/heh15/workingspace/Arp240/NGC5257/SFR/'
workDir=Dir+'SFR_map/'
logDir=Dir+'log/'
picDir=Dir+'picture/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'

fitsimage=imageDir+'spitzer_24um.fits'
wcs_spi=fits_import(fitsimage)[0]
data_24um=fits_import(fitsimage)[1]

file=regionDir+'stellar_total.reg'
whole_sky=read_ds9(file)[0]
whole_pix=whole_sky.to_pixel(wcs_spi)
whole24um_masked=Apmask_convert(whole_pix,data_24um)
f24=np.ma.sum(whole24um_masked)

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[1].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[1].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

fitsimage=imageDir+'herschel_70um.fits'
wcs_her=fits_import(fitsimage)[0]
data_70um=fits_import(fitsimage)[1]

file=regionDir+'stellar_total.reg'
whole_sky=read_ds9(file)[0]
whole_pix=whole_sky.to_pixel(wcs_her)
whole70um_masked=Apmask_convert(whole_pix,data_70um)
f24=np.ma.sum(whole70um_masked)

whole70um_masked=Apmask_convert(whole_pix,data_70um)
f70=np.ma.sum(whole70um_masked)


L=TIR(f24,f70, pixelsize=2.45)
LTIR=L/(3.328e33)
