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

D=99

############################################################
# function  
def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
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

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/SM/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'

from regions import read_ds9
fitsimage=imageDir+'mass_map.fits'
hdr = fits.open(fitsimage)[0].header
wcs=fits_import(fitsimage)[0]
mass=fits_import(fitsimage)[1]

file=regionDir+'stellar_total.reg'
whole_star=read_ds9(file)[0]
wholestar_pix=whole_star.to_pixel(wcs)
whole_masked=Apmask_convert(wholestar_pix, mass)
SM=np.ma.sum(whole_masked)*(480*0.3)**2

#### check total stellar mass

flux_36=19448.31*1e6/(4.25e10)*0.3**2
flux_45=13959*1e6/(4.25e10)*0.3**2

Mstar_total=mass_calc(flux_36, flux_45)
