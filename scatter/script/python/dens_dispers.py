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
import seaborn as sns
import cube

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
rms=0.0024
rms_mom0=rms*10*math.sqrt(50) 
incl=0.45

H_coeff=10**6/(2e30*3.14*6.67e-11/(3.1e16))

############################################################
# function
def Apmask_convert(aperture,data_masked):
    data_cut=data_masked.data
    data_mask=data_masked.mask
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    mask_tmp=mask==0
    ap_mask=np.ma.mask_or(mask_tmp,data_mask)
 #   ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_mask


##################################################
# main program
Dir='/1/home/heh15/workingspace/Arp240/scatter/'
picDir=Dir+'picture/'
imageDir=Dir+'image/'
logDir=Dir+'log/'

############################################################ 
galaxies=['NGC5257','NGC5258']
sd_all=dict.fromkeys(galaxies)
vd_all=dict.fromkeys(galaxies)
H_all=dict.fromkeys(galaxies)

# load the data for NGC5257

# regbin the image with 5 pixels together.
image_mom0=imageDir+'NGC5257/NGC5257_12CO21_combine_5bin_mom0.fits'
image_mom2=imageDir+'NGC5257/NGC5257_12CO21_combine_5bin_mom2.fits'
nchan=imageDir+'test_2chan/NGC5257_12CO21_combine_rebin_nchan_mask.fits'

incl=0.45
hdr=fits.open(image_mom0)[0].header
wcs=WCS(hdr).celestial
data_mom0=fits.open(image_mom0)[0].data[0][0]*incl

## selection for the data
threshold=data_mom0<5.0*rms_mom0
data_mom0[threshold]='nan'

data_K=data_mom0/(0.0109*beammaj*beammin*(freq/115.27)**2)
sd=alpha/ratio*data_K
x=sd

hdr=fits.open(image_mom2)[0].header
wcs=WCS(hdr).celestial
data_mom2=fits.open(image_mom2)[0].data[0][0]
position=np.isnan(data_mom2)

y=np.copy(data_mom2) 

## Do the channel selection
# hdr=fits.open(nchan)[0].header
# wcs=WCS(hdr).celestial
# data_nchan=fits.open(nchan)[0].data[0][0]
# chan1_value_flags=data_nchan<5
low_value_flags=data_mom2<1.0
# y[chan1_value_flags]='nan'
y[low_value_flags]='nan'

## exclude the point within 1.5 arcsec.

# create mask with NAN data
mask1=np.isnan(data_mom2)
mask2=np.isnan(data_mom0)
mask=np.ma.mask_or(mask1,mask2)
data_masked=np.ma.masked_where(mask,y)

a=1.5;incl=0.45
b=1.5*incl;PA=110
ra=204.97066052
dec=50.406569999999995
position=SkyCoord(dec=dec*u.arcmin,ra=ra*u.degree,frame='icrs')
core_sky=SkyEllipticalAperture(position,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
core_pix=core_sky.to_pixel(wcs=wcs)
core_mask=Apmask_convert(core_pix,data_masked)
core_masked=np.ma.masked_where(core_mask,y)
disk_masked=np.ma.masked_where(~core_mask,y)

# ## plot the scattering plot
# x=sd
# fig=plt.figure()
# ax=plt.axes(xscale='log',yscale='log')
# plt.plot(x.flatten(),y.flatten(),'.',color='red',label='disk')
# plt.plot(x.flatten(),core_masked.flatten(),'.',color='blue',label='core')
# # plt.ylim(top=300,bottom=8)
# # plt.xlim(left=40)
# plt.ylabel('velocity dispersion (km/s)')
# plt.xlabel('surface density (solar mass/pc**2)')
# plt.title('conversion factor 4.3')
# plt.legend()
# plt.savefig(picDir+'NGC5257_sd_disper_43.png')


P=x*np.power(core_masked,2)
alpha1=5.77*np.power(disk_masked,2)/x*40/165
alpha2=5.77*np.power(core_masked,2)/x*40/165

H=H_coeff*y**2/x
H_disk=H_coeff*np.power(disk_masked,2)/x
H_core=H_coeff*np.power(core_masked,2)/x

# fig=plt.figure()
# plt.plot(x,P,'.',color='black')


# fig=plt.figure()
# ax=plt.axes(xscale='log',yscale='log')
# plt.plot(x.flatten(),alpha1.flatten(),'.',color='black',label='disk')
# plt.plot(x.flatten(),alpha2.flatten(),'.',color='blue',label='core')
# plt.xlabel(r'$\Sigma(M_{\odot}/pc^2)$')
# plt.ylabel(r'Virial parameter')
# plt.legend()
# plt.title('conversion factor 4.3')
# plt.savefig(picDir+'NGC5257_alpha_sd.png')

############################################################
# load the data for NGC5258

## moment 0 map. 

# regbin the image with 5 pixels together.
image_mom0=imageDir+'NGC5258_12CO21_combine_5bin_mom0.fits'
image_mom2=imageDir+'NGC5258_12CO21_combine_5bin_mom2.fits'

hdr=fits.open(image_mom0)[0].header
wcs=WCS(hdr).celestial
data_mom0=fits.open(image_mom0)[0].data[0][0]*incl
threshold=data_mom0<5.0*rms_mom0
data_mom0[threshold]='nan'
data_K=data_mom0/(0.0109*beammaj*beammin*(freq/115.27)**2)
sd=alpha/ratio*data_K

## moment 2 map

hdr=fits.open(image_mom2)[0].header
wcs=WCS(hdr).celestial
data_mom2=fits.open(image_mom2)[0].data[0][0]
position=np.isnan(data_mom2)

vd=np.copy(data_mom2) 
low_value_flags=data_mom2<1.0
vd[low_value_flags]='nan'

## plot the scattering point

fig=plt.figure()
ax=plt.axes(xscale='log',yscale='log')
plt.scatter(sd,vd,marker='.',color='red')
plt.scatter(x,y,marker='.',color='blue')
plt.xlabel('surface density')
plt.ylabel('dispersion (km/s)')
plt.savefig(picDir+'initial.png')

## the virial parameter
alpha8=5.77*np.power(vd,2)/sd*40/165
alpha7=5.77*np.power(y,2)/x*40/165

## The scale height
H_8=H_coeff*vd**2/sd

fig=plt.figure()
ax=plt.axes(xscale='log',yscale='log')
plt.scatter(sd,alpha8,marker='.',color='red',label='NGC5258')
plt.scatter(x,alpha7,marker='.',color='blue',label='NGC5257')
plt.legend()
plt.xlabel('surface density(solar mass/pc^2)')
plt.ylabel('virial parameter')
plt.savefig(picDir+'Virial.png')

fig=plt.figure()
ax=plt.axes(xscale='log',yscale='log')
plt.scatter(x.flatten(),H_disk.flatten(),marker='.',label='NGC 5257 disk')
plt.scatter(x.flatten(),H_core.flatten(),marker='.',label='NGC 5257 core')
plt.scatter(sd.flatten(),H_8.flatten(),marker='.',label='NGC 5258')
plt.ylabel(r'H(pc)')
plt.xlabel(r'$\Sigma_{M_{\odot}}$')
plt.legend()
plt.savefig(picDir+'Height.png')

############################################################
# ##  load the data from the table

# filename=logDir+'Sigma_sigv.txt'
# skiprows=28
# data=pd.read_csv(filename,header=None,sep=r"\s*",skiprows=skiprows)
# data_120pc=data.loc[(data[1]==120)]
# object=data_120pc[0].unique()
# special=np.array(['M31','M33','Antennae'])
# object_main=np.setdiff1d(object,special)


# data_M31=data.loc[(data[0]=='M31')]
# data_M33=data.loc[(data[0]=='M33')]
# frames=[data_M31,data_M33]
# data_M3=pd.concat(frames)
# data_ant=data.loc[(data[0]=='Antennae')]
# data_main=data.loc[data[0].isin(object_main)]

# fig=plt.figure()
# ax=plt.axes(xscale='log',yscale='log')
# plt.scatter(data_M31[3],data_M31[4],marker='.',color='green')
# plt.scatter(data_M33[3],data_M33[4],marker='.',color='green')
# plt.scatter(data_main[3],data_main[4],marker='.',color='blue',alpha=0.2)
# plt.scatter(data_ant[3],data_ant[4],marker='.',color='yellow',alpha=0.2)
# x1=plt.plot(x.flatten(),disk_masked.flatten(),'.',color='red',label='conversion factor 4.3 disk')
# x2=plt.plot(x.flatten(),core_masked.flatten(),'.',color='black',label='conversion factor 4.3 core')
# z=x/5.0
# x3=plt.plot(z.flatten(),disk_masked.flatten(),'.',color='cyan',label='conversion factor 0.86')
# x4=plt.plot(z.flatten(),core_masked.flatten(),'.',color='gray',label='conversion factor 0.86')
# #sns.kdeplot(data_main[3],data_main[4],n_levels=6,cmap='Blues')
# # plt.plot(x,y,'.',color='black')
# plt.legend(fontsize=8)
# plt.savefig(picDir+'all.png')
