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
ratio=0.77
XCO=0.5
D=99
z=0.024
beammaj=1.004
beammin=0.556
freq=230.54
alpha=4.3

############################################################
# function

##################################################
# main program

image_mom0='NGC5257_12CO21_combine_noise40_rebin_mom0.fits'
image_mom2='NGC5257_12CO21_combine_noise40_rebin_mom2.fits'
nchan='test_2chan/NGC5257_12CO21_combine_rebin_nchan_mask.fits'

hdr=fits.open(image_mom0)[0].header
wcs=WCS(hdr).celestial
data_mom0=fits.open(image_mom0)[0].data[0][0]
data_K=data_mom0/(0.0109*beammaj*beammin*(freq/115.27)**2)
sd=alpha/ratio*data_K

hdr=fits.open(image_mom2)[0].header
wcs=WCS(hdr).celestial
data_mom2=fits.open(image_mom2)[0].data[0][0]
position=np.isnan(data_mom2)
data_mom0[position]='nan'

hdr=fits.open(nchan)[0].header
wcs=WCS(hdr).celestial
data_nchan=fits.open(nchan)[0].data[0][0]
chan1_value_flags=data_nchan<2

y=np.copy(data_mom2)
low_value_flags=data_mom2<0.1
y[chan1_value_flags]='nan'
# y[low_value_flags]='nan'
x=sd
fig=plt.figure()
ax=plt.axes(xscale='log',yscale='log')
plt.plot(x,y,'.',color='blue')
# plt.ylim(top=300,bottom=8)
# plt.xlim(left=40)
plt.ylabel('velocity dispersion (km/s)')
plt.xlabel('surface density (solar mass/arcsec**2)')
plt.title('conversion factor 0.86')
plt.savefig('NGC5257_sd_disper_086.png')

P=np.log10(x*np.power(y,2))
H=np.log10(np.power(y,2)/data_mom0)

# fig=plt.figure()
# plt.plot(x,P,'.',color='black')

# fig=plt.figure()
# plt.plot(x,H,'.',color='black')

# load the data from the table


filename='Sigma_sigv.txt'
skiprows=28
data=pd.read_csv(filename,header=None,sep=r"\s*",skiprows=skiprows)
data_120pc=data.loc[(data[1]==120)]
object=data_120pc[0].unique()
special=np.array(['M31','M33','Antennae'])
object_main=np.setdiff1d(object,special)


data_M31=data.loc[(data[0]=='M31')]
data_M33=data.loc[(data[0]=='M33')]
data_ant=data.loc[(data[0]=='Antennae')]
data_main=data.loc[data[0].isin(object_main)]

fig=plt.figure()
ax=plt.axes(xscale='log',yscale='log')
plt.scatter(data_M31[3],data_M31[4],marker='.',color='green')
plt.scatter(data_M33[3],data_M33[4],marker='.',color='green')
plt.scatter(data_main[3],data_main[4],marker='.',color='blue')
plt.scatter(data_ant[3],data_ant[4],marker='.',color='yellow')
plt.plot(x,y,'.',color='black')
plt.savefig('all.png')
