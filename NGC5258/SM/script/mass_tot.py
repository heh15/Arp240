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

############################################################
# main program

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/SM/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'
fitsDir=Dir+'Fits/'

from regions import read_ds9
fitsimage=fitsDir+'mass_map.fits'
hdr = fits.open(fitsimage)[0].header
wcs=fits_import(fitsimage)[0]
mass=fits_import(fitsimage)[1]

file=regionDir+'stellar_total.reg'
whole_star=read_ds9(file)[0]
wholestar_pix=whole_star.to_pixel(wcs)
whole_masked=Apmask_convert(wholestar_pix, mass)
SM=np.ma.sum(whole_masked)*(480*0.3)**2


