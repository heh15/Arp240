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

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# directories and files

Dir='/home/heh15/workingspace/Arp240/NGC5257/SFR/'
logDir=Dir+'log/'
picDir=Dir+'picture/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'

wavelengths=['33 GHz', '12CO21']
regions=['center', 'south', 'arm', 'west']
quantities=['dispersion','dispersion_err', 'SFR', 'SFR_err', 'SD_err']

fitsfiles=dict.fromkeys(wavelengths)
fitsfiles['33 GHz']=imageDir+'NGC5257_33GHz_pbcor_smooth_co21_regrid.fits'
fitsfiles['12CO21']=imageDir+'NGC5257_12CO21_pbcor_smooth_cube_signal_mom2.fits'
fitsfiles['12CO21 mon0'] = imageDir+'NGC5257_12CO21_combine_pbcor_mom0.fits'

regionfiles=dict.fromkeys(regions)
regionfiles['center']=regionDir+'center_co2-1.reg'
regionfiles['south']=regionDir+'south_co2-1.reg'
regionfiles['arm']=regionDir+'arm_co2-1.reg'
regionfiles['west']=regionDir+'west_co2-1.reg'

regionobjects=dict.fromkeys(regions)

regionproperties=pd.DataFrame(columns=quantities, index=regions)
# regionproperties=dict.fromkeys(regions)
# for key in regionproperties.keys():
#     regionproperties[key]=dict.fromkeys(quantities)

############################################################
# basic information
galaxy='NGC5257'
linelabel='33GHz'
alpha=1.1
ratio=0.85

d=99
ra=204.9706
dec=50.4167
center_point= SkyCoord(dec=dec*u.arcmin,ra=ra*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)

beamra=204*u.degree+58*u.arcmin+30*u.arcsec
beamdec=50*u.arcmin+3*u.arcsec
beamposition=SkyCoord(dec=beamdec,ra=beamra,frame='icrs')
beammajor=1.1*u.arcsec/2.0
beamminor=0.8*u.arcsec/2.0
pa=-64.5*u.degree

rms_33GHz=1.0e-5

rms_CO=3.1e-3/(0.0109*1.1*0.8*(225.46/115.27)**2)
rms_mom0=rms_CO*10*math.sqrt(50)
SD_err_pix=rms_mom0/ratio*alpha

############################################################
# basic settings. 

testra = 204.97609228
testdec = 0.84611111

continuum = '33 GHz'

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

def sfr_24um(f_mjsr):
    flux=f_mjsr*10**(-17)/(4.25*10**10)*2.45**2
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

def Regmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

def sfr_radio(flux,freq,d):
    '''
    flux in Jy
    freq in Hz
    d in Mpc
    '''
    f_erg=flux*1e-23
    L=4*math.pi*(d*1e6*3.086e18)**2*f_erg
    SFR_tot=1e-27*(2.18*freq**-0.1+15.1*freq**-0.7)**(-1)*L

    return SFR_tot

def ssfr_radio(intensity,beammajor,beamminor,freq,d):
    '''
    flux in Jy
    beammajor and beamminor in arcsec
    freq in Hz
    d in Mpc
    sSFR in msun/kpc^2/yr
    '''
    SFR_tot=sfr_radio(intensity,freq,d)
    sSFR=SFR_tot/(1.1331*beammajor*beamminor)/(4.85*d/1000)**2

    return sSFR

############################################################
# main program

#### Import the 12CO2-1 mom0 maps
wcs, data=fits_import(fitsfiles['12CO21 mon0'])
wcs_cut, contour_cut=cut_2d(data, center_point, size, wcs)

#### Draw the map of 33 GHz coninuum. 
wcs, data=fits_import(fitsfiles['33 GHz'])
wcs_cut, data_cut=cut_2d(data, center_point, size, wcs)

for key in regionfiles.keys():
    regionobject=read_ds9(regionfiles[key])[0]    
    regionobjects[key]=regionobject

for key in regionobjects.keys():
    region_pix=regionobjects[key].to_pixel(wcs_cut)
    region_mask=Regmask_convert(region_pix, data_cut)
    intensity=np.ma.mean(region_mask)
    number=np.ma.count(region_mask)
    regionproperties['SFR'][key]=ssfr_radio(intensity, 1.1, 0.8, 33, d)
    regionproperties['SFR_err']=rms_33GHz/np.sqrt(number)/intensity*regionproperties['SFR'][key]
    print(str(number)+'\n')


# import the star cluster coordinates. 
filename=logDir+'NGC5257_cluster.txt'
colname = ['RA', 'Dec','MB', 'MB_err', 'MI', 'MI_err', 'MFUV', 'MFUV_err'] 
info=pd.read_csv(filename, names = colname, sep=r"\s*", skiprows=5, index_col=0)

filename = logDir+'NGC5257_derived.txt'
colname = ['Age', 'Age_err', 'Mass', 'Mass_err', 'Av', 'Av_err']
derived = pd.read_csv(filename, names=colname, sep=r"\s*", skiprows=5, index_col=0)

clusters = pd.merge(info, derived, left_index=True, right_index=True)
clusters_young = clusters.loc[(clusters['Age'] < 7.0) & (clusters['Mass'] > 6.0)]

coordinates = clusters_young[['RA', 'Dec']]

positions=list();apertures=list()
for i in coordinates.index:
    position=SkyCoord(ra=coordinates['RA'][i]*u.degree,dec=coordinates['Dec'][i]*u.degree, frame='fk5')
    position_icrs=position.transform_to('icrs')
    positions.append(position_icrs)
    aperture=SkyCircularAperture(positions=position_icrs, r=0.3*u.arcsec)
    apertures.append(aperture)


### Draw figures. 
fig=plt.figure()
ax=plt.subplot(projection=wcs_cut)
im=plt.imshow(data_cut*1000, 
		origin='lower', 
		cmap='viridis_r',
		vmax=0.4,
		vmin=0.0)

# coordinates. 
plt.xlabel('J2000 Right Ascension')
plt.ylabel('J2000 Declination')

ax.tick_params(labelsize=8, direction='in')
cbar=plt.colorbar(im)
cbar.set_label('$mJy \ beam^{-1} $', fontsize=15)

apertures_pix=list()
for aperture in apertures:
    aperture_pix=aperture.to_pixel(wcs_cut)
    apertures_pix.append(aperture_pix)
    aperture_pix.plot(color='magenta', linewidth=1.0)


ax.contour(contour_cut, levels=[1.1, 2.2], colors='black',  linewidths=0.6)

x, y = wcs_cut.wcs_world2pix(testra, testdec, 1)
plt.text(x, y, galaxy + ' ' + continuum, fontsize=15)

for key in regionfiles.keys():
    regionobject=regionobjects[key].to_pixel(wcs_cut)
    regionobject.plot(color='red', linewidth=1.5)

beamellipse = SkyEllipticalAperture(positions=beamposition,
                                    a=beammajor,
                                    b=beamminor,
                                    theta=pa)
beamellipse_pix = beamellipse.to_pixel(wcs_cut)
beamellipse_pix.plot(color='black')
# fig.tight_layout()
plt.savefig(picDir+galaxy+'_'+linelabel+'_'+'_pbcor.png', bbox_inches = 'tight')

#### measure the median velocity dispersion in each region. 
wcs, data=fits_import(fitsfiles['12CO21'])
wcs_cut, data_cut=cut_2d(data, center_point, size, wcs)

for key in regionobjects.keys():
    region_pix=regionobjects[key].to_pixel(wcs_cut)
    region_mask=Regmask_convert(region_pix, data_cut)
    regionproperties['dispersion'][key]=np.ma.median(region_mask)
    regionproperties['dispersion_err'][key]=np.ma.std(region_mask)
    regionproperties['SD_err'][key]=SD_err_pix/np.sqrt(np.ma.count(region_mask))

filename=logDir+galaxy+'_aperture_output.csv'
regionproperties.to_csv(filename)
