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
ratioDir=Dir+'image/ratio/1213/contsub_pbcor/'
regionDir=Dir+'region/'
imageDir=Dir+'image/'

beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.09
ra= 204.9908
dec= 0.8319
center_point=SkyCoord(dec=0.8319*u.degree,ra=204.9908*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)

regions=['northarm','southarm','center','ring', 'southsfr','total']
regions_ap=['northarm','southarm','center','ring']
values=['12CO10 flux','12CO10 uncertainty','13CO10 flux', '13CO10 uncertainty', 'ratio','uncertainty']
ratio_table=pd.DataFrame(index=regions,columns=values)
apertures=dict.fromkeys(regions_ap)


############################################################
# function 

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data[0][0] 
    cut=Cutout2D(data=data,position=center_point,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs
    data_masked=np.ma.masked_invalid(data_cut)
    return wcs_cut, data_masked

def masked_convert(data_masked,region_masked):
    data_mask=data_masked.mask
    region_mask=np.ma.make_mask(region_masked==0)
    region_mask=np.ma.mask_or(data_mask,region_mask)
    data_region=np.ma.masked_where(region_mask,data_masked)
    return data_region

def flux_mask_get(data_region,rms,chans,chan_width):
    flux=np.ma.sum(data_region)/beam_area_pix
    chans_tmp=chans+np.zeros((np.shape(data_region)[0],np.shape(data_region)[1]))
    error=np.sqrt(chans_tmp)*rms*chan_width/sqrt(beam_area_pix)
    error_masked=np.ma.masked_where(data_region.mask,error)
    uncertainty=math.sqrt(np.ma.sum(np.power(error_masked,2)))
    return flux, uncertainty

def flux_aperture_get(data_masked,aperture,rms,chans,chan_width):
    data_cut=data_masked.data
    mask=data_masked.mask
    flux=aperture_photometry(data_cut,apertures=aperture,mask=mask)['aperture_sum'][0]
    error=np.sqrt(chans)*rms*chan_width/sqrt(beam_area_pix)
    uncertainty=aperture_photometry(data_cut,apertures=aperture,mask=mask,error=error)['aperture_sum_err'][0]

    return flux, uncertainty

def ratio_cal(flux1,flux2,uncertainty1,uncertainty2):
    ratio=flux1/flux2*107**2/112**2
    uncertainty=math.sqrt((uncertainty1/flux1)**2+(uncertainty2/flux2)**2+2*0.05**2)*ratio
    return ratio, uncertainty

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

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
fitsimage=ratioDir+'NGC5258_12CO10_combine_contsub_uvrange_smooth_masked_pbcor_mom0.fits'
wcs=fits_import(fitsimage)[0]
data12_masked=fits_import(fitsimage)[1]

# import the channel used map. 
# fitsimage='NGC5258_12CO10_combine_smooth_nchan.fits'
# chans=fits_import(fitsimage)[1]
chans=50
chans_12=chans

# define the aperture
position=SkyCoord(dec=0.8309*u.degree,ra=204.9906*u.degree,frame='icrs')
center_sky =SkyCircularAperture(position,r=3*u.arcsec)
center_pix = center_sky.to_pixel(wcs=wcs)
apertures['center']=center_pix

ring_sky=SkyCircularAnnulus(position,r_in=3*u.arcsec,r_out=7*u.arcsec)
ring_pix=ring_sky.to_pixel(wcs=wcs)
apertures['ring']=ring_pix

position=SkyCoord(dec=0.8349*u.degree,ra=204.9930*u.degree,frame='icrs')
northarm_sky=SkyEllipticalAperture(position,a=15*u.arcsec,b=7*u.arcsec,theta=340*u.degree)
northarm_pix=northarm_sky.to_pixel(wcs=wcs)
apertures['northarm']=northarm_pix

position=SkyCoord(dec=0.8275*u.degree,ra=204.9884*u.degree,frame='icrs')
southarm_sky=SkyEllipticalAperture(position,a=10*u.arcsec,b=5*u.arcsec,theta=340*u.degree)
southarm_pix=southarm_sky.to_pixel(wcs=wcs)
apertures['southarm']=southarm_pix

data12_masked=data12_masked/beam_area_pix
rms_12CO10=1.6e-3
chan_width=10

for region in regions_ap:
    flux, uncertainty=flux_aperture_get(data12_masked,apertures[region],rms_12CO10,chans,chan_width)
    ratio_table['12CO10 flux'][region]=flux
    ratio_table['12CO10 uncertainty'][region]=uncertainty


### south spiral arm regions. 
from regions import read_ds9

file=regionDir+'south_SFR.reg'

southsfr_sky=read_ds9(file,errors='warn')[0]
southsfr_pix=southsfr_sky.to_pixel(wcs)
southsfr_masked=Regmask_convert(southsfr_pix, data12_masked*beam_area_pix)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(southsfr_masked,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(southsfr_masked,rms,chans, chan_width)[1]

ratio_table['12CO10 flux']['southsfr']=flux_12CO10
ratio_table['12CO10 uncertainty']['southsfr']=uncertainty_12CO10

### total 
rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(data12_masked*beam_area_pix,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(data12_masked*beam_area_pix,rms,chans, chan_width)[1]

ratio_table['12CO10 flux']['total']=flux_12CO10
ratio_table['12CO10 uncertainty']['total']=uncertainty_12CO10

############################################################
# calculate the flux of 13CO10 data
#position=SkyCoord(dec=0.8319*u.degree,ra=204.9908*u.degree,frame='icrs')

# import the 13CO10 data. 
fitsimage=ratioDir+'NGC5258_13CO10_12m_uvrange_smooth_masked_pbcor_mom0.fits'
wcs=fits_import(fitsimage)[0]
data13_masked=fits_import(fitsimage)[1]

# import the channel used map. 
fitsimage=imageDir+'13CO10/NGC5258_13CO10_12m_smooth_nchan.fits'
chans=fits_import(fitsimage)[1]
chans_13=chans

# measure the flux 
rms_13CO10=0.00062
chan_width=20
data13_masked=data13_masked/beam_area_pix

for region in regions_ap:
    flux, uncertainty=flux_aperture_get(data13_masked,apertures[region],rms_13CO10,chans,chan_width)
    ratio_table['13CO10 flux'][region]=flux
    ratio_table['13CO10 uncertainty'][region]=uncertainty

### south spiral arm regions. 
from regions import read_ds9

file=regionDir+'south_SFR.reg'

southsfr_sky=read_ds9(file,errors='warn')[0]
southsfr_pix=southsfr_sky.to_pixel(wcs)
southsfr_masked=Regmask_convert(southsfr_pix, data13_masked*beam_area_pix)

rms=0.0006
chan_width=20
flux_13CO10=flux_mask_get(southsfr_masked,rms,chans,chan_width)[0]
uncertainty_13CO10=flux_mask_get(southsfr_masked,rms,chans, chan_width)[1]

ratio_table['13CO10 flux']['southsfr']=flux_13CO10
ratio_table['13CO10 uncertainty']['southsfr']=uncertainty_13CO10

### total region. 
rms=0.0006
chan_width=20
flux_13CO10=flux_mask_get(data13_masked*beam_area_pix,rms,chans,chan_width)[0]
uncertainty_13CO10=flux_mask_get(data13_masked*beam_area_pix,rms,chans, chan_width)[1]

ratio_table['13CO10 flux']['total']=flux_13CO10
ratio_table['13CO10 uncertainty']['total']=uncertainty_13CO10

############################################################
# calculate the ratio

for region in regions:
    ratio_tmp=ratio_cal(ratio_table['12CO10 flux'][region],
                    ratio_table['13CO10 flux'][region],
                    ratio_table['12CO10 uncertainty'][region],
                    ratio_table['13CO10 uncertainty'][region])
    ratio_table['ratio'][region]=ratio_tmp[0]
    ratio_table['uncertainty'][region]=ratio_tmp[1]

############################################################
# output the data frame result

os.chdir(scriptDir)
ratio_table.to_csv(path_or_buf='../log/NGC5258_1213_ratio_error_cal.csv')


# # draw the image
# data=data12_masked.data*beam_area_pix
# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# im=ax.imshow(data,origin='lower')
# for region in regions:
#     apertures[region].plot(color='red')
# fig.colorbar(im)
# plt.savefig('../picture/NGC5258_12CO10_aperture.png')
