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

PA= 110
incl=0.45
xcenter=159.71
ycenter=162.11
ra=204.97066052
dec=50.406569999999995
radius_arcsec=np.array([0.75, 2.25, 3.75, 5.25, 6.75,8.25,9.75,11.25,12.75])
radius_kpc=radius_arcsec*0.48
size=radius_arcsec.shape[0]-1
position=SkyCoord(dec=dec*u.arcmin,ra=ra*u.degree,frame='icrs')
rings=dict.fromkeys((range(size)))
rings_mask=dict.fromkeys((range(size)))
pixel_area=0.3*0.3
pixel_sr=pixel_area/(60**2*180/math.pi)**2
D=99

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

def aperture_ring(radius_arcsec,wcs):
    a_in=radius_arcsec-1.5
    a_out=radius_arcsec
    b_out=a_out*incl
    ring_sky=SkyEllipticalAnnulus(position,a_in=a_in*u.arcsec,a_out=a_out*u.arcsec,b_out=b_out*u.arcsec,theta=PA*u.degree)
    ring_pix=ring_sky.to_pixel(wcs=wcs)
    return ring_pix

def Apmask_convert(aperture,data_masked):
    data_cut=data_masked.data
    data_mask=data_masked.mask
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    mask_tmp=mask==0
    ap_mask=np.ma.mask_or(mask_tmp,data_mask)
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

def cut_3d(data,position,size,wcs):
    for i in range(data_3d.shape[0]):
        cut=Cutout2D(data=data[i],position=position,size=size,wcs=wcs)
        if i==0:
            data_cut=cut.data
        elif i==1:
            data_cut=np.stack((data_cut,cut.data))
        else:
            temp=np.expand_dims(cut.data,axis=0)
            data_cut=np.concatenate((data_cut,temp))
    wcs_cut=cut.wcs
    return data_cut, wcs_cut

############################################################
# Directory

Dir='/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/'
imageDir=Dir+'image/'
picDir=Dir+'picture/'


## comparision for 12CO10 data
fitsimage=imageDir+'NGC5257_12CO10_combine_contsub_mom2.fits'
wcs=fits_import(fitsimage)[0]
data_masked=fits_import(fitsimage)[1]
data=data_masked.data

a=radius_arcsec[0]
b=radius_arcsec[0]*incl
center_sky=SkyEllipticalAperture(position,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
center_pix=center_sky.to_pixel(wcs=wcs)
center_mask=Apmask_convert(center_pix,data_masked)

for i in range(8):
    rings[i]=aperture_ring(radius_arcsec[i+1],wcs)
    rings_mask[i]=Apmask_convert(rings[i],data_masked)

dispersions=np.empty((0,0))
dispersions=np.append(dispersions,np.ma.median(center_mask))
for i in range(8):
    dispersion=np.ma.median(rings_mask[i])
    dispersions=np.append(dispersions,dispersion)

# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# im=ax.imshow(data,cmap='rainbow',origin='lower',vmax=150)
# rings[6].plot(color='red')
# cbar=fig.colorbar(im)
# lon = ax.coords[0]
# lat = ax.coords[1]
# lon.set_major_formatter('hh:mm:ss')
# plt.savefig('../picture/CO10_dispersion.png')

## load the data from Bbarolo
dir='/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/Bbarolo/NGC5257_12CO10/'
dataB=np.transpose(np.loadtxt(dir+'ringlog1.txt'))
disperB=dataB[3]
radius=np.delete(radius_arcsec,-1)

# fig=plt.figure()
# mom2=plt.plot(radius_arcsec,dispersions,label='mom2 dispersion')
# barolo=plt.plot(radius,disperB,label='Bbarolo')
# plt.legend()
# plt.xlabel('velocity dispersion (km/s)')
# plt.ylabel('radius (arcsec)')
# plt.savefig('dispersion_radius.png')

############################################################
## comparison for 12CO21 data 

radius_arcsec=np.linspace(0.5,11.5,12)

fitsimage=imageDir+'NGC5257_12CO21_combine_noise40_mom2.fits'
wcs=fits_import(fitsimage)[0]
data_masked=fits_import(fitsimage)[1]
data=data_masked.data

a=radius_arcsec[0]
b=radius_arcsec[0]*incl
center_sky=SkyEllipticalAperture(position,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
center_pix=center_sky.to_pixel(wcs=wcs)
center_mask=Apmask_convert(center_pix,data_masked)

for i in range(11):
    rings[i]=aperture_ring(radius_arcsec[i+1],wcs)
    rings_mask[i]=Apmask_convert(rings[i],data_masked)

dispersions=np.empty((0,0))
dispersions=np.append(dispersions,np.ma.median(center_mask))
for i in range(11):
    dispersion=np.ma.median(rings_mask[i])
    dispersions=np.append(dispersions,dispersion)

# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# im=ax.imshow(data,cmap='rainbow',origin='lower',vmax=150)
# rings[6].plot(color='red')
# cbar=fig.colorbar(im)
# lon = ax.coords[0]
# lat = ax.coords[1]
# lon.set_major_formatter('hh:mm:ss')
# plt.savefig('../picture/CO10_dispersion.png')

## load the data from Bbarolo
dir='/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/Bbarolo/NGC5257_12CO21/'
dataB=np.transpose(np.loadtxt(dir+'ringlog1.txt'))
disperB=dataB[3]
radius=np.delete(radius_arcsec,-1)

fig=plt.figure()
mom2=plt.plot(radius_arcsec,dispersions,label='mom2 dispersion')
barolo=plt.plot(radius_arcsec,disperB,label='Bbarolo')
plt.legend()
plt.ylabel('velocity dispersion (km/s)')
plt.xlabel('radius (arcsec)')

