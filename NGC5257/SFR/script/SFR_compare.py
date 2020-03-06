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
from regions import read_ds9

############################################################
Dir='/home/heh15/workingspace/Arp240/NGC5257/SFR/'
workDir=Dir+'SFR_map/'
logDir=Dir+'log/'
picDir=Dir+'picture/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'
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

ratio=1.24

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
    LTIR=10**1.336*(L/3.83e33)**0.954
    
    return SFR_tot, LTIR

def sfr_70um(f_her):
    flux_erg=f_her*10**(-23)
    L_nu=4*math.pi*(100*10**6*3.086*10**18)**2*flux_erg
    L=4.283*10**12*L_nu
    SFR_tot=10**(-43.23)*L
    LTIR=10**0.567*(L/3.83e33)**0.973

    return SFR_tot, LTIR

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

# convert region to the mask. 
def Regmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
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

def flux_mask_get(data_region,rms,chans,chan_width):
    flux=np.ma.sum(data_region)/beam_area_pix
    chans_tmp=chans+np.zeros((np.shape(data_region)[0],np.shape(data_region)[1]))
    error=np.sqrt(chans_tmp)*rms*chan_width/sqrt(beam_area_pix)
    error_masked=np.ma.masked_where(data_region.mask,error)
    uncertainty=math.sqrt(np.ma.sum(np.power(error_masked,2)))
    return flux, uncertainty


############################################################
# main program

wavelength=['24um','70um','33GHz','24nd70']
value=['region','flux','uncertainty','SFR', 'SFR_TIR']
df=pd.DataFrame(index=wavelength,columns=value)
df_center=pd.DataFrame(index=wavelength,columns=value)
df_arm=pd.DataFrame(index=wavelength,columns=value)

### south region 

## 33 GHz image
rms_33=5.64e-5
pixsize=((5.81776417e-07*u.rad)**2).to(u.arcsec**2)  # from the imhead in casa
beamarea=6*u.arcsec*6*u.arcsec*1.1331

# flux measured from original image
flux_33GHz=6.8e-4;error=1.4e-5
flux_unpb=6.4e-4
# compare to the smoothed image
filename=imageDir+'NGC5257_33GHz_pbcor_smooth_6arcsec.fits'
wcs_33GHz=fits_import(filename)[0]
data_33GHz=fits_import(filename)[1]

ra=-2.7057736283005824;dec=0.014605178692092591
position=SkyCoord(dec=dec*u.rad,ra=ra*u.rad,frame='icrs')
south_sky=SkyCircularAperture(positions=position,r=3*ratio*u.arcsec)
south_pix=south_sky.to_pixel(wcs_33GHz)
flux=aperture_photometry(data_33GHz,apertures=south_pix)['aperture_sum'][0]
flux=float(flux/(beamarea/pixsize))

freq=33
SFR_33GHz=sfr_radio(flux_33GHz,freq,d) 

error=np.full(np.shape(data_33GHz), rms_33)
flux_error=aperture_photometry(data_33GHz,apertures=south_pix,error=error)['aperture_sum_err'][0]/(np.sqrt(1.1331*6*6/0.12**2))
# error=np.sqrt(np.ma.count(arm_masked)/(1.1331*6*6/0.12**2))*rms_33

df['region']['33GHz']='south'
df['flux']['33GHz']=flux_33GHz
df['SFR']['33GHz']=SFR_33GHz
df['uncertainty']['33GHz']=flux_error/flux_33GHz*SFR_33GHz

# Modify the Murphy equation. 

# def sfr_radio(flux,freq,d):
#     '''
#     flux in Jy/beam
#     freq in Hz
#     d in Mpc
#     '''
#     f_erg=flux*1e-23
#     L=4*math.pi*(d*1e6*3.086e18)**2*f_erg
#     SFR_tot=1e-27*(2.18*freq**-0.1+15.1*freq**-0.85)**(-1)*L

#     return SFR_tot

# test=sfr_radio(flux_33GHz,freq,d) 


## herschel flux 

fitsimage=imageDir+'herschel_70um.fits'
header=fits.open(fitsimage)[1].header
wcs_her=WCS(header).celestial
data_70um=fits.open(fitsimage)[1].data

ra=204*u.degree+58*u.arcmin+15*u.arcsec
dec=50*u.arcmin+13*u.arcsec
position=SkyCoord(ra=ra,dec=dec,frame='icrs')
south_sky=SkyCircularAperture(positions=position,r=3*ratio*u.arcsec)
south_pix=south_sky.to_pixel(wcs_her)


# fig=plt.figure()
# ax=plt.subplot('111',projection=wcs_her)
# ax.imshow(data_70um,origin='lower',vmax=0.05)
# south_pix.plot()

flux=aperture_photometry(data_70um,apertures=south_pix)['aperture_sum'][0]
SFR_70um=sfr_70um(flux)[0]
LTIR_70um=sfr_70um(flux)[1]
df['SFR_TIR']['70um']=LTIR_70um*3.83e33*10**(-43.41)

df['region']['70um']='south'
df['flux']['70um']=flux
df['SFR']['70um']=SFR_70um


## spitzer flux
fitsimage=imageDir+'spitzer_24um.fits'
header=fits.open(fitsimage)[0].header
wcs_spi=WCS(header).celestial
data_24um=fits.open(fitsimage)[0].data-48 # subtract the noise

ra=204*u.degree+58*u.arcmin+13*u.arcsec
dec=50*u.arcmin+13*u.arcsec
position=SkyCoord(ra=ra,dec=dec,frame='icrs')
south_sky=SkyCircularAperture(positions=position,r=3*ratio*u.arcsec)
south_pix=south_sky.to_pixel(wcs_spi)

# fig=plt.figure()
# ax=plt.subplot('111',projection=wcs_spi)
# ax.imshow(data_24um,origin='lower')
# south_pix.plot()

flux_mjsr=aperture_photometry(data_24um,apertures=south_pix)['aperture_sum'][0]
pixelsize=2.45
flux=flux_mjsr*10**6/(4.25*10**10)*pixelsize**2
SFR_24um=sfr_24um(flux_mjsr)[0]
LTIR_24um=sfr_24um(flux_mjsr)[1]
df['SFR_TIR']['24um']=LTIR_24um*3.83e33*10**(-43.41)

df['region']['24um']='south'
df['flux']['24um']=flux
df['SFR']['24um']=SFR_24um

SFR=sfr_24and70(flux_mjsr, df['flux']['70um'])
df['SFR']['24nd70']=SFR[0]
df['uncertainty']['24nd70']=SFR[0]*SFR[1]

### center

## 33 GHz image
pixsize=((5.81776417e-07*u.rad)**2).to(u.arcsec**2)  # from the imhead in casa
beamarea=6*u.arcsec*6*u.arcsec*1.1331

# flux measured from original image
flux_33GHz=5.25e-4;error=1.06e-5
# compare to the smoothed image
filename=imageDir+'NGC5257_33GHz_pbcor_smooth_6arcsec.fits'
wcs_33GHz=fits_import(filename)[0]
data_33GHz=fits_import(filename)[1]

ra=-2.7057752478718142,;dec=0.014663001991737822
position=SkyCoord(dec=dec*u.rad,ra=ra*u.rad,frame='icrs')
center_sky=SkyCircularAperture(positions=position,r=3*ratio*u.arcsec)
center_pix=center_sky.to_pixel(wcs_33GHz)
flux=aperture_photometry(data_33GHz,apertures=center_pix)['aperture_sum'][0]
flux=float(flux/(beamarea/pixsize))

freq=33
SFR_33GHz=sfr_radio(flux_33GHz,freq,d) 

error=np.full(np.shape(data_33GHz), rms_33)
flux_error=aperture_photometry(data_33GHz,apertures=center_pix,error=error)['aperture_sum_err'][0]/(np.sqrt(1.1331*6*6/0.12**2))

df_center['region']['33GHz']='center'
df_center['flux']['33GHz']=flux_33GHz
df_center['SFR']['33GHz']=SFR_33GHz
df_center['uncertainty']['33GHz']=flux_error/flux_33GHz*SFR_33GHz

## herschel flux 

fitsimage=imageDir+'herschel_70um.fits'
header=fits.open(fitsimage)[1].header
wcs_her=WCS(header).celestial
data_70um=fits.open(fitsimage)[1].data

ra=204*u.degree+58*u.arcmin+15*u.arcsec
dec=50*u.arcmin+25*u.arcsec
position=SkyCoord(ra=ra,dec=dec,frame='icrs')
center_sky=SkyCircularAperture(positions=position,r=3*ratio*u.arcsec)
center_pix=center_sky.to_pixel(wcs_her)


# fig=plt.figure()
# ax=plt.subplot('111',projection=wcs_her)
# ax.imshow(data_70um,origin='lower',vmax=0.05)
# center_pix.plot)

flux=aperture_photometry(data_70um,apertures=center_pix)['aperture_sum'][0]
SFR_70um=sfr_70um(flux)[0]

LTIR_70um=sfr_70um(flux)[1]
df_center['SFR_TIR']['70um']=LTIR_70um*3.83e33*10**(-43.41)

df_center['region']['70um']='center'
df_center['flux']['70um']=flux
df_center['SFR']['70um']=SFR_70um

## spitzer flux
fitsimage=imageDir+'spitzer_24um.fits'
header=fits.open(fitsimage)[0].header
wcs_spi=WCS(header).celestial
data_24um=fits.open(fitsimage)[0].data-48 # subtract the noise

ra=204*u.degree+58*u.arcmin+14*u.arcsec
dec=50*u.arcmin+23*u.arcsec
position=SkyCoord(ra=ra,dec=dec,frame='icrs')
center_sky=SkyCircularAperture(positions=position,r=3*ratio*u.arcsec)
center_pix=center_sky.to_pixel(wcs_spi)

# fig=plt.figure()
# ax=plt.subplot('111',projection=wcs_spi)
# ax.imshow(data_24um,origin='lower',vmax=95)
# center_pix.plot()

flux_mjsr=aperture_photometry(data_24um,apertures=center_pix)['aperture_sum'][0]
pixelsize=2.45
flux=flux_mjsr*10**6/(4.25*10**10)*pixelsize**2
SFR_24um=sfr_24um(flux_mjsr)[0]

LTIR_24um=sfr_24um(flux_mjsr)[1]
df_center['SFR_TIR']['24um']=LTIR_24um*3.83e33*10**(-43.41)

df_center['region']['24um']='center'
df_center['flux']['24um']=flux
df_center['SFR']['24um']=SFR_24um

SFR=sfr_24and70(flux_mjsr, df_center['flux']['70um'])
df_center['SFR']['24nd70']=SFR[0]
df_center['uncertainty']['24nd70']=SFR[0]*SFR[1]

### south west arm

## 33 GHz
rms_33=5.64e-5
# flux=7.56e-4 for 33 GHz continuum

# fig=plt.figure() 
# ax=plt.subplot('111',projection=wcs_33GHz)
# ax.imshow(data_33GHz,origin='lower')
# ax.contour(data_33GHz, levels=[3e-4])

from regions import read_ds9
file=regionDir+'NGC5257_arm.reg'
arm_sky=read_ds9(file)[0]
arm_pix=arm_sky.to_pixel(wcs_33GHz)
# arm_pix.plot(color='red')
arm_masked=Regmask_convert(arm_pix, data_33GHz)
flux=np.ma.sum(arm_masked)/(1.1331*6*6/0.12**2)
error=np.sqrt(np.ma.count(arm_masked)/(1.1331*6*6/0.12**2))*rms_33
SFR=sfr_radio(flux, 33, 100)

df_arm['region']['33GHz']='arm'
df_arm['flux']['33GHz']=flux
df_arm['SFR']['33GHz']=SFR
df_arm['uncertainty']['33GHz']=error/flux*SFR


## spitzer

filename=regionDir+'NGC5257_arm.reg'
arm_sky=read_ds9(filename)[0]
arm_pix=arm_sky.to_pixel(wcs_spi)
arm_masked=Regmask_convert(arm_pix, data_24um)

fig=plt.figure()
ax=plt.subplot('111', projection=wcs_spi)
plt.imshow(data_24um, origin='lower')
arm_pix.plot(color='red')

flux_24um=np.ma.sum(arm_masked)
pixelsize=2.45
flux=flux_24um*10**6/(4.25*10**10)*pixelsize**2
SFR_24um=sfr_24um(flux_24um,pixelsize=2.45)[0]

LTIR_24um=sfr_24um(flux_24um, pixelsize=2.45)[1]
df_arm['SFR_TIR']['24um']=LTIR_24um*3.83e33*10**(-43.41)

df_arm['region']['24um']='arm'
df_arm['flux']['24um']=flux
df_arm['SFR']['24um']=SFR_24um

## herschel
arm_sky=read_ds9(filename)[0]
arm_pix=arm_sky.to_pixel(wcs_her)
arm_masked=Regmask_convert(arm_pix, data_70um)

fig=plt.figure()
ax=plt.subplot('111', projection=wcs_her)
ax.imshow(data_70um, origin='lower')
arm_pix.plot(color='red')

flux_70um=np.ma.sum(arm_masked)
SFR_70um=sfr_70um(flux_70um)[0]

LTIR_70um=sfr_70um(flux_70um)[1]
df_arm['SFR_TIR']['70um']=LTIR_70um*3.83e33*10**(-43.41)

df_arm['region']['70um']='arm'
df_arm['flux']['70um']=flux_70um
df_arm['SFR']['70um']=SFR_70um

# spitzer and herschel combined
result=sfr_24and70(flux_24um,flux_70um,pixelsize=2.45)
SFR_combine=result[0]
error_combine=SFR_combine*result[1]

df_arm['region']['24nd70']='arm'
df_arm['flux']['24nd70']=flux_70um
df_arm['SFR']['24nd70']=SFR_combine
df_arm['uncertainty']['24nd70']=error_combine

### save the data frame

# with open(logDir+'SFR.txt') as out:
#     df.to_string(out)



############################################################
## test region
