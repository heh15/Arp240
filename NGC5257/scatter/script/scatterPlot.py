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

Dir='/1/home/heh15/workingspace/Arp240/ratio/'
images={}
images['12CO10']='NGC5258_12CO10_combine_uvrange_smooth_regrid21_masked_rebin.fits'
images['12CO21']='NGC5258_12CO21_combine_uvtaper_smooth_masked_rebin.fits'
images['ratio']='NGC5258_2110_ratio_uvtaper_rebin.fits'

wcs={}
data={}


for image in images.keys(): 
    hdr = fits.open(images[image])[0].header
    wcs[image] = WCS(hdr).celestial
    data[image]=fits.open(images[image])[0].data[0][0]

y=data['ratio'].ravel()
x=data['12CO10'].ravel()
plt.plot(x,y,'.')
