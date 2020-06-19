'''
Apr. 24th
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
# basic settings

Dir='/home/heh15/workingspace/Arp240/NGC5258/ratio/'
scriptDir=Dir+'script/'
ratioDir=Dir+'image/ratio/2110/uvtaper_mask/'
regionDir=Dir+'region/'
imageDir=Dir+'image/'

beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.01
ra= 204.9908
dec= 0.8319
center_point=SkyCoord(dec=0.8319*u.degree,ra=204.9908*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)
freq10=112.73
freq21= 225.46

regions=['northarm','southarm','center','ring','southsfr', 'total']
regions_ap=['northarm','southarm','center','ring']
values=['12CO10 flux','12CO10 uncertainty','12CO21 flux', '12CO21 uncertainty', 'ratio','uncertainty']
ratio_table=pd.DataFrame(index=regions,columns=values)
apertures=dict.fromkeys(regions_ap)


############################################################
# function 

def fits_import(fitsimage, item=0):
    hdr = fits.open(fitsimage)[item].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[item].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def masked_convert(data_masked,region_masked):
    data_mask=data_masked.mask
    region_mask=np.ma.make_mask(region_masked==0)
    region_mask=np.ma.mask_or(data_mask,region_mask)
    data_region=np.ma.masked_where(region_mask,data_masked)
    return data_region

def flux_mask_get(data_region,rms,chans,chan_width):
    flux=np.ma.sum(data_region)/beam_area_pix
    error_value=rms*chan_width*sqrt(50/beam_area_pix)
    error = np.full((180,140),error_value)
    uncertainty=error_value*math.sqrt(np.ma.count(data_region))
    return flux, uncertainty

def flux_aperture_get(data_masked,aperture,rms,chans,chan_width):
    data_cut=data_masked.data
    mask=data_masked.mask
    flux=aperture_photometry(data_cut,apertures=aperture,mask=mask)['aperture_sum'][0]
    if np.shape(chans) == ():
        chans = np.full(np.shape(data_cut), chans)
    error=np.sqrt(chans)*rms*chan_width/sqrt(beam_area_pix)
    uncertainty=aperture_photometry(data_cut,apertures=aperture,mask=mask,error=error)['aperture_sum_err'][0]

    return flux, uncertainty

def ratio_cal(flux1,flux2,uncertainty1,uncertainty2):
    ratio=flux1/flux2*freq10**2/freq21**2
    uncertainty=math.sqrt((uncertainty1/flux1)**2+(uncertainty2/flux2)**2+2*0.05**2)*ratio
    return ratio, uncertainty

def Regmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

##########################################################
# main program

##############################
# calculate the flux of 12CO10 data

# import the 12CO10 data
# fitsimage=imageDir+'12CO10/NGC5258_12CO10_uvrange_pbcor_cube_masked.fits'
fitsimage=ratioDir+'NGC5258_12CO10_combine_contsub_uvrange_smooth_masked_pbcor_mom0.fits'
wcs=fits_import(fitsimage)[0]
data10_masked=fits_import(fitsimage)[1]

# import the channel used map. 
# fitsimage='NGC5258_12CO10_combine_uvrange_smooth_regrid21_nchan.fits'
# chans=fits_import(fitsimage)[1]
chans=50
chans_10=chans

# define the aperture
position=SkyCoord(dec=0.8309*u.degree,ra=204.9906*u.degree,frame='icrs')
center_sky =SkyCircularAperture(position,r=3*u.arcsec)
center_pix = center_sky.to_pixel(wcs=wcs)
apertures['center']=center_pix

ring_sky=SkyCircularAnnulus(position,r_in=3*u.arcsec,r_out=7*u.arcsec)
ring_pix=ring_sky.to_pixel(wcs=wcs)
apertures['ring']=ring_pix

position=SkyCoord(dec=0.8340*u.degree,ra=204.9935*u.degree,frame='icrs')
northarm_sky=SkyEllipticalAperture(position,a=13*u.arcsec,b=4*u.arcsec,theta=185*u.degree)
northarm_pix=northarm_sky.to_pixel(wcs=wcs)
apertures['northarm']=northarm_pix

position=SkyCoord(dec=0.8283*u.degree,ra=204.9882*u.degree,frame='icrs')
southarm_sky = SkyEllipticalAperture(position,a=9*u.arcsec,b=4*u.arcsec,theta=360*u.degree)
southarm_pix = southarm_sky.to_pixel(wcs=wcs)
apertures['southarm'] = southarm_pix


# measure the flux 
rms_12CO10=0.00158
chan_width=10
data10_masked=data10_masked/beam_area_pix

for region in regions_ap:
    flux, uncertainty=flux_aperture_get(data10_masked,apertures[region],rms_12CO10,chans,chan_width)
    ratio_table['12CO10 flux'][region]=flux
    ratio_table['12CO10 uncertainty'][region]=uncertainty


## south arm sfr. 
from regions import read_ds9

file=regionDir+'south_SFR.reg'

southsfr_sky=read_ds9(file,errors='warn')[0]
southsfr_pix=southsfr_sky.to_pixel(wcs)
southsfr_masked=Regmask_convert(southsfr_pix, data10_masked*beam_area_pix)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(southsfr_masked,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(southsfr_masked,rms,chans, chan_width)[1]

ratio_table['12CO10 flux']['southsfr']=flux_12CO10
ratio_table['12CO10 uncertainty']['southsfr']=uncertainty_12CO10

## total region
rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(data10_masked*beam_area_pix,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(data10_masked*beam_area_pix,rms,chans, chan_width)[1]

ratio_table['12CO10 flux']['total']=flux_12CO10
ratio_table['12CO10 uncertainty']['total']=uncertainty_12CO10

############################################################
# calculate the flux of 12CO21 data
#position=SkyCoord(dec=0.8319*u.degree,ra=204.9908*u.degree,frame='icrs')

# import the 12CO21 data. 
# fitsimage=imageDir+'12CO21/NGC5258_12CO21_uvtaper_pbcor_cube_masked.fits'
fitsimage=ratioDir+'NGC5258_12CO21_combine_uvtaper_smooth_masked_pbcor_mom0.fits'
wcs=fits_import(fitsimage)[0]
data21_masked=fits_import(fitsimage)[1]

# import the channel used map. 
# fitsimage='NGC5258_12CO21_combine_uvtaper_smooth_nchan.fits'
# chans=fits_import(fitsimage)[1]
chans=50
chans_21=chans

# measure the flux 
rms_12CO21= 0.0047
chan_width=10
data21_masked=data21_masked/beam_area_pix

for region in regions_ap:
    flux, uncertainty=flux_aperture_get(data21_masked,apertures[region],rms_12CO21,chans,chan_width)
    ratio_table['12CO21 flux'][region]=flux
    ratio_table['12CO21 uncertainty'][region]=uncertainty

## south arm sfr. 
from regions import read_ds9

file=regionDir+'south_SFR.reg'

southsfr_sky=read_ds9(file,errors='warn')[0]
southsfr_pix=southsfr_sky.to_pixel(wcs)
southsfr_masked=Regmask_convert(southsfr_pix, data21_masked*beam_area_pix)

rms=0.0047
chan_width=10
flux_12CO21=flux_mask_get(southsfr_masked,rms,chans,chan_width)[0]
uncertainty_12CO21=flux_mask_get(southsfr_masked,rms,chans, chan_width)[1]

ratio_table['12CO21 flux']['southsfr']=flux_12CO21
ratio_table['12CO21 uncertainty']['southsfr']=uncertainty_12CO21

## total regions
rms=0.0047
chan_width=10
flux_12CO21=flux_mask_get(data21_masked*beam_area_pix,rms,chans,chan_width)[0]
uncertainty_12CO21=flux_mask_get(data21_masked*beam_area_pix,rms,chans, chan_width)[1]

ratio_table['12CO21 flux']['total']=flux_12CO21
ratio_table['12CO21 uncertainty']['total']=uncertainty_12CO21

############################################################
# calculate the ratio

for region in regions:
    ratio_tmp=ratio_cal(ratio_table['12CO21 flux'][region],
                    ratio_table['12CO10 flux'][region],
                    ratio_table['12CO21 uncertainty'][region],
                    ratio_table['12CO10 uncertainty'][region])
    ratio_table['ratio'][region]=ratio_tmp[0]
    ratio_table['uncertainty'][region]=ratio_tmp[1]

############################################################
# output the data frame result

os.chdir(scriptDir)
ratio_table.to_csv(path_or_buf='../log/NGC5258_2110_ratio_error_cal.csv')

# draw apertures in the image
fig=plt.figure()
ax=plt.subplot(projection=cut.wcs)
ax.imshow(data_cut,cmap='gray',origin='lower')
apertures['center']['pix'].plot(color='red')
apertures['ring around center']['pix'].plot(color='red')
apertures['northarm']['pix'].plot(color='red')
apertures['southarm']['pix'].plot(color='red')
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_major_formatter('hh:mm:ss')
plt.show()
