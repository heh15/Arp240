'''
Oct.9th

Apmask_convert doesn't include the nan value in the data. 
'''

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

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
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
fitsimage='NGC5257_12CO10_combine_contsub_mom2.fits'
wcs=fits_import(fitsimage)[0]
data_masked=fits_import(fitsimage)[1]
data=data_masked.data

a=radius_arcsec[0]
b=radius_arcsec[0]*incl
center_sky=SkyEllipticalAperture(position,a=a*u.arcsec,b=b*u.arcsec,theta=PA*u.degree)
center_pix=center_sky.to_pixel(wcs=wcs)
center_mask=Apmask_convert(center_pix,data)

# draw the ring of the image
for i in range(8):
    rings[i]=aperture_ring(radius_arcsec[i+1],wcs)
    rings_mask[i]=Apmask_convert(rings[i],data)


# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# im=ax.imshow(data,cmap='rainbow',origin='lower',vmax=150)
# rings[5].plot(color='red')
# cbar=fig.colorbar(im)
# lon = ax.coords[0]
# lat = ax.coords[1]
# lon.set_major_formatter('hh:mm:ss')
# plt.savefig('../picture/CO10_dispersion.png')


# import the line of the aperture ring. 
fitsimage='NGC5257_12CO10_combine_contsub_image.fits'

hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
data_3d=fits.open(fitsimage)[0].data[0]
position=SkyCoord(dec=50.4167*u.arcmin,ra=204.9706*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)
data_cut,wcs_cut=cut_3d(data_3d,position,size,wcs)

i=7
mask=rings_mask[3].mask
mask_3d = np.repeat(mask[np.newaxis,:, :], 70, axis=0)
data_3dma=np.ma.masked_where(mask_3d,data_cut)

spectrum=np.empty((0,0))

# for i in range(data_3dma.shape[0]):
#     # temp=np.ma.sum(data_3dma[i]) 
#     temp=np.ma.masked_where(mask,data_cut[i])
#     value=np.ma.sum(temp)
#     spectrum=np.append(spectrum,value)

for i in range(data_3dma.shape[0]):
    temp=np.ma.sum(data_3dma[i])
    spectrum=np.append(spectrum,temp)

velocity=np.linspace(-300,390,70)
plt.plot(velocity,spectrum)
