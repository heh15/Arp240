'''
Oct. 25th, 2018

make the star formation rate map

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

############################################################
Dir='/home/heh15/workingspace/Arp240/NGC5258/SFR/'
scriptDir=Dir+'script/'
picDir=Dir+'picture/'
logDir=Dir+'log/'
workDir=Dir+'SFR_map/'


# basic setting (NGC5258)

ra= 204*u.degree+58*u.arcmin+50*u.arcsec
dec=50*u.arcmin+7*u.arcsec
center_point= SkyCoord(dec=dec,ra=ra,frame='icrs')
size=u.Quantity((125,125),u.arcsec)

############################################################
# function

def fits_import(fitsimage, item=0):
    hdr = fits.open(fitsimage)[item].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[item].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def cut_2d(data_masked,position,size,wcs):
    cut=Cutout2D(data=data_masked,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

############################################################
# main program

os.chdir(workDir)

# ## add the header to the image. 
# filename='herschel_70um.fits'
# hdul=fits.open('herschel_70um.fits')
# hdr=hdul[1].header

# hdr['BMAJ']=0.0015166666666666666
# hdr['BMIN']=0.0016000000000000000
# hdr['BPA']=0.00

# hdul.writeto(filename.strip('.fits')+'_header.fits')
# hdul.close()

# filename='galex_FUV.fits'
# hdul=fits.open(filename)
# hdr=hdul[0].header

# hdr['BMAJ']=4.3/3600
# hdr['BMIN']=4.3/3600
# hdr['BPA']=0.00

# hdul.writeto(filename.strip('.fits')+'_header.fits')
# hdul.close()

# filename='spitzer_24um.fits'
# hdul=fits.open(filename)
# hdr=hdul[0].header

# hdr['BMAJ']=6.0/3600
# hdr['BMIN']=6.0/3600
# hdr['BPA']=0.00

# hdul.writeto(filename.strip('.fits')+'_header.fits')
# hdul.close()

## execute the imsmooth.py

## add the data together. 

filename='herschel_70um_back.fits'
wcs=fits_import(filename)[0]
data_70um=fits_import(filename)[1]

def cut_2d(data_masked,position,size,wcs):
    data=data_masked.data
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

herschel_wcs_cut=cut_2d(data_70um, center_point, size, wcs)[0]
herschel_cut=cut_2d(data_70um, center_point,size, wcs)[1]

filename='spitzer_24um_regrid.fits'
wcs=fits_import(filename)[0]
# hdul=fits.open(filename)
# hdr=hdul[0].header
# wcs=WCS(hdr)
# data_24um=hdul[0].data
data_24um=fits_import(filename)[1]
data_24um=data_24um-48
# hdul.close()

Spitzer_wcs_cut=cut_2d(data_24um, center_point, size, wcs)[0]
Spitzer_cut=cut_2d(data_24um, center_point,size, wcs)[1]

filename='galex_FUV_header_smooth_regrid.fits'
# hdul=fits.open(filename)
# hdr=hdul[0].header
# wcs=WCS(hdr)
# data_FUV=hdul[0].data
# hdul.close()
data_FUV=fits_import(filename)[1]
data_FUV_tmp=data_FUV*1.4*10**(-15)*1528**2/(3*10**18)
data_FUV_tmp=data_FUV_tmp*10**17/(1.5**2)*4.25*10**10


SFR=8.1*10**(-2)*data_FUV_tmp+3.2*10**-3*data_24um
SFR_ob=3.2*10**-3*data_24um
SFR_uob=8.1*10**(-2)*data_FUV_tmp

fig=plt.figure()
ax=plt.subplot('111',projection=Spitzer_wcs_cut)
ax.title.set_text('spitzer 24um')
im=ax.imshow(Spitzer_cut,origin='lower')

fig=plt.figure()
ax1=plt.subplot('111',projection=herschel_wcs_cut)
ax1.title.set_text('herschel 70um')
im=ax1.imshow(herschel_cut,origin='lower')

# fig=plt.figure()
# ax=plt.subplot('121',projection=wcs)
# ax.title.set_text('spitzer 24um')
# im=ax.imshow(SFR_ob,origin='lower',vmax=0.5,vmin=0)
# ax2=plt.subplot('122',projection=wcs)
# ax2.imshow(SFR_uob,origin='lower',vmax=0.5,vmin=0)
# ax2.title.set_text('galex FUV')
# plt.draw()
# p0 = ax.get_position().get_points().flatten()
# p1 = ax2.get_position().get_points().flatten()
# ax_cbar = fig.add_axes([p0[0],0.05, p1[2]-p0[0], 0.05])
# ax_cbar.set_label('$M_{\solar}/(kpc^2) $'  )
# plt.colorbar(im, cax=ax_cbar, orientation='horizontal')

# fig=plt.figure()
# ax=plt.subplot('111',projection=wcs)
# im=ax.imshow(SFR_uob,origin='lower',vmax=0.02)
# plt.colorbar(im)

counts=33
flux=counts*1.4*10**(-15)*1528**2/(3*10**18)*0.68*10**(23)
flux_erg=flux/10**23
L_uv=4*math.pi*(100*10**6*3.086*10**18)**2*flux_erg
test=0.68*10**(-28)*L_uv
print(test)

pixel_area=2.45*2.45
pixel_sr=pixel_area/(60**2*180/math.pi)**2

f_mjsr=9503
flux=f_mjsr*10**(-17)/(4.25*10**10)*2.45**2
flux=f_mjsr*pixel_sr*10**6
#flux=1.34*10**(-23)
L_nu=4*math.pi*(100*10**6*3.086*10**18)**2*flux
L=L_nu*1.2*10**13
SFR_spitzer=10**(-42.69)*L

# center_sky=SkyCircularAperture(position,a=a*u.arcsec)
# center_pix=center_sky.to_pixel(wcs=wcs_cut)
# center_mask=Apmask_convert(center_pix,data_cut)


L_uvnu=L_uv*1.96*10**15+3.89*L
SFR=10**(-43.35)*L_uvnu
