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

fitsimage='NGC5257_12CO10_vel_9.res.fits'
hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
data=fits.open(fitsimage)[0].data

fitsimage='../NGC5257_12CO10_combine_contsub_mom1.fits'
hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
data_CO=fits.open(fitsimage)[0].data
data_CO=data_CO[0][0]

fig=plt.figure()
ax=plt.subplot()
ax.imshow(data,origin='lower')

fig=plt.figure()
ax=plt.subplot()
ax.imshow(data_CO,origin='lower')
