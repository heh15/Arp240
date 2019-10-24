'''
Oct. 25th, 2018

Compare the SFR measured by different tracers

'''
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

############################################################
Dir='/home/heh15/workingspace/Arp240/NGC5258/SFR/'
workDir=Dir+'SFR_map/'
logDir=Dir+'log/'
picDir=Dir+'picture/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
# os.chdir(workDir)


############################################################
# basic setting (NGC5257)

ra=204.9706
dec=50.4167
center_point= SkyCoord(dec=dec*u.arcmin,ra=ra*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)
d=99

# common pixelsize while regriding 
pix_size=1.6

ratio=2.0


############################################################
# function

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    # cut=Cutout2D(data=data,position=center_point,size=size,wcs=wcs)
    # data_cut=cut.data
    # wcs_cut=cut.wcs
    return wcs, data

def sfr_24um(f_mjsr,pixelsize=2.45):
    flux=f_mjsr*10**(-17)/(4.25*10**10)*pixelsize**2
    L_nu=4*math.pi*(100*10**6*3.086*10**18)**2*flux
    L=L_nu*1.2*10**13
    SFR_tot=10**(-42.69)*L
    
    return SFR_tot

def sfr_70um(f_her):
    flux_erg=f_her*10**(-23)
    L_nu=4*math.pi*(100*10**6*3.086*10**18)**2*flux_erg
    L=4.283*10**12*L_nu
    SFR_tot=10**(-43.23)*L

    return SFR_tot

def sfr_radio(flux,freq,d):
    '''
    flux in Jy/beam
    freq in Hz
    d in Mpc
    '''
    f_erg=flux*1e-23
    L=4*math.pi*(d*1e6*3.086e18)**2*f_erg
    SFR_tot=1e-27*(2.18*freq**-0.1+15.1*freq**-0.7)**(-1)*L

    return SFR_tot

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask(method='center')[0]
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)

array_round=np.vectorize(round_sig)

def flux_aperture_get(data_cut,aperture):
    flux=aperture_photometry(data_cut,apertures=aperture)['aperture_sum'][0]

    return flux

def sfr_24and70(f24,f70,pixelsize=2.45):
    '''
    from Galametz et al. 2013
    f24 in MJy/sr
    f70 in Jy
    pixelsize of spitzer in arcsec
    '''
    flux=f24*10**(-17)/(4.25*10**10)*pixelsize**2
    L24=4*math.pi*(100*10**6*3.086*10**18)**2*flux
    nuL24=L24*1.2*10**13
    flux_erg=f70*10**(-23)
    L70=4*math.pi*(100*10**6*3.086*10**18)**2*flux_erg
    nuL70=4.283*10**12*L70
    LTIR=3.98*nuL24+1.553*nuL70
    error=0.283/3.98+0.058/1.553
    SFR=10**-43.41*LTIR
    
    return SFR, error


############################################################
# main program

wavelength=['24um','70um','24nd70', '33GHz']
value=['region','flux','uncertainty','SFR']
df=pd.DataFrame(index=wavelength,columns=value)
df_south=pd.DataFrame(index=wavelength,columns=value)

### south arm peak

## 33 GHz image
pixsize=((5.81776417e-07*u.rad)**2).to(u.arcsec**2)  # from the imhead in casa
beamarea=6*u.arcsec*6*u.arcsec*1.1331

# # measure the flux
# fitsimage='ngc5258_Ka_c_r0.5_ms.fits'
# fitsimagepb= 'ngc5258_Ka_c_r0.5_ms.pbcor.fits'

# # find the shape to measure the flux
# wcs=fits_import(fitsimage)[0]
# data=fits_import(fitsimage)[1]

# flux measured from original image via casa.
flux_33GHz=2.16e-3;error=1.1e-5/3.4e-4*flux_33GHz
# compare to the smoothed image
filename=imageDir+'NGC5258_33GHz_pbcor_smooth.fits'
wcs_33GHz=fits_import(filename)[0]
data_33GHz=fits_import(filename)[1]

ra=(13*u.degree+39*u.arcmin+57.14*u.arcsec)*15 ;dec=49*u.arcmin+44.309*u.arcsec
position=SkyCoord(dec=dec,ra=ra,frame='icrs')
southarm_sky=SkyCircularAperture(positions=position,r=3*ratio*u.arcsec)
southarm_pix=southarm_sky.to_pixel(wcs_33GHz)
flux=aperture_photometry(data_33GHz,apertures=southarm_pix)['aperture_sum'][0]
flux=float(flux/(beamarea/pixsize))

freq=33
SFR_33GHz=sfr_radio(flux_33GHz,freq,d) 

df['region']['33GHz']='south'
df['flux']['33GHz']=flux_33GHz
df['SFR']['33GHz']=SFR_33GHz
df['uncertainty']['33GHz']=error/flux_33GHz*SFR_33GHz

fig=plt.figure()
plt.imshow(data_33GHz,origin='lower')
southarm_pix.plot(color='red')

## herschel flux 

fitsimage=imageDir+'herschel_70um.fits'
header=fits.open(fitsimage)[1].header
wcs_her=WCS(header).celestial
data_70um=fits.open(fitsimage)[1].data

southarm_sky=SkyCircularAperture(positions=position,r=3*ratio*u.arcsec)
southarm_pix=southarm_sky.to_pixel(wcs_her)


fig=plt.figure()
ax=plt.subplot('111',projection=wcs_her)
ax.imshow(data_70um,origin='lower',vmax=0.05)
southarm_pix.plot()

flux=aperture_photometry(data_70um,apertures=southarm_pix)['aperture_sum'][0]
SFR_70um=sfr_70um(flux)

df['region']['70um']='south'
df['flux']['70um']=flux
df['SFR']['70um']=SFR_70um


## spitzer flux
fitsimage=imageDir+'spitzer_24um.fits'
header=fits.open(fitsimage)[0].header
wcs_spi=WCS(header).celestial
data_24um=fits.open(fitsimage)[0].data-48 # subtract the noise

southarm_sky=SkyCircularAperture(positions=position,r=3*ratio*u.arcsec)
southarm_pix=southarm_sky.to_pixel(wcs_spi)

# fig=plt.figure()
# ax=plt.subplot('111',projection=wcs_spi)
# ax.imshow(data_24um,origin='lower')
# southarm_pix.plot()

flux_mjsr=aperture_photometry(data_24um,apertures=southarm_pix)['aperture_sum'][0]
pixelsize=2.45
flux=flux_mjsr*10**6/(4.25*10**10)*pixelsize**2
SFR_24um=sfr_24um(flux_mjsr)

df['region']['24um']='south'
df['flux']['24um']=flux
df['SFR']['24um']=SFR_24um


## combine 24um and 70um flux

SFR=sfr_24and70(df['flux']['24um'],df['flux']['70um'])
df['SFR']['24nd70']=SFR[0]
df['uncertainty']['24nd70']=SFR[0]*SFR[1]
