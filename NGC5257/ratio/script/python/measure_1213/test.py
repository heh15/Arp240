'''
Nov. 29th, 2018

Change the input file from the 
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

Dir='/home/heh15/workingspace/Arp240/NGC5257/ratio/'
scriptDir=Dir+'script/python/measure_1213/'
workDir=Dir+'mask/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
ratioDir=Dir+'1213/contsub/'
maskDir=Dir+'mask/'

beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.09
regions=['spiral','anomaly','restdisk','center']
values=['12CO10 flux','12CO10 uncertainty','13CO10 flux', '13CO10 uncertainty', 'ratio','uncertainty']
CO13=['13CO10 flux', '13CO10 uncertainty']
ratio_table=pd.DataFrame(index=regions,columns=values)

chans=50

Im12CO10=imageDir+'12CO10/NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor_mom0.fits'
Im13CO10=imageDir+'13CO10/NGC5257_13CO10_12m_contsub_smooth_pbcor_mom0.fits'

# Im12CO10=ratioDir+'NGC5257_12CO10_combine_contsub_smooth_masked_mom0.fits'
# Im13CO10=ratioDir+'NGC5257_13CO10_12m_contsub_smooth_masked_mom0.fits'

Im12CO10=imageDir+'12CO10/NGC5257_12CO10_combine_contsub_smooth_13COcut_pbcor_mom0.fits'
Im13CO10=imageDir+'13CO10/NGC5257_13CO10_12m_contsub_smooth_pbcor_mom0.fits'
outfile='test.txt'
comment='Image from no cut image \n'

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

def ratio_cal(flux1,flux2,uncertainty1,uncertainty2):
    ratio=flux1/flux2
    uncertainty=math.sqrt((uncertainty1/flux1)**2+(uncertainty2/flux2)**2+2*0.05**2)*ratio
    return ratio, uncertainty


############################################################
os.chdir(maskDir)

# calculating 12CO10 flux
fitsimage=Im12CO10
wcs=fits_import(fitsimage)[0]
data12_masked=fits_import(fitsimage)[1]

# # import the channel used map. 
# fitsimage='NGC5257_12CO10_combine_smooth_nchan.fits'
# wcs=fits_import(fitsimage)[0]
# chans=fits_import(fitsimage)[1]

# spiral arm
maskimage='spiralarm_mask.fits'
spiral_masked=fits_import(maskimage)[1]
spiral= masked_convert(data12_masked,spiral_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(spiral,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(spiral,rms,chans, chan_width)[1]
ratio_table['12CO10 flux']['spiral']=flux_12CO10
ratio_table['12CO10 uncertainty']['spiral']=uncertainty_12CO10


# anomaly 
maskimage='anomaly_mask.fits'
anomaly_masked=fits_import(maskimage)[1]
anomaly=masked_convert(data12_masked,anomaly_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(anomaly,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(anomaly,rms,chans, chan_width)[1]

ratio_table['12CO10 flux']['anomaly']=flux_12CO10
ratio_table['12CO10 uncertainty']['anomaly']=uncertainty_12CO10

# rest disk

maskimage='disk_mask.fits'
restdisk_masked=fits_import(maskimage)[1]
restdisk=masked_convert(data12_masked,restdisk_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(restdisk,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(restdisk,rms,chans,chan_width)[1]

ratio_table['12CO10 flux']['restdisk']=flux_12CO10
ratio_table['12CO10 uncertainty']['restdisk']=uncertainty_12CO10

# # nonarm 
# maskimage='nonarm_mask.fits'
# nonarm_masked=fits_import(maskimage)[1]
# nonarm=masked_convert(data12_masked,nonarm_masked)

# rms=0.0016
# chan_width=10
# flux_12CO10=flux_mask_get(nonarm,rms,chans,chan_width)[0]
# uncertainty_12CO10=flux_mask_get(nonarm,rms,chans,chan_width)[1]

# ratio_table['12CO10 flux']['nonarm']=flux_12CO10
# ratio_table['12CO10 uncertainty']['nonarm']=uncertainty_12CO10

# center
maskimage='center_mask.fits'
center_masked=fits_import(maskimage)[1]
center=masked_convert(data12_masked,center_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(center,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(center,rms,chans, chan_width)[1]

ratio_table['12CO10 flux']['center']=flux_12CO10
ratio_table['12CO10 uncertainty']['center']=uncertainty_12CO10

############################################################
# calculating the flux of 13CO data 
chans=25

fitsimage=Im13CO10
wcs=fits_import(fitsimage)[0]
data13_masked=fits_import(fitsimage)[1]

# # import the channel used map. 
# fitsimage='NGC5257_13CO10_12m_smooth_nchan.fits'
# wcs=fits_import(fitsimage)[0]
# chans=fits_import(fitsimage)[1]

# spiral arm 
maskimage='spiralarm_mask.fits'
spiral_masked=fits_import(maskimage)[1]
spiral= masked_convert(data13_masked,spiral_masked)

rms=0.0006
chan_width=20
flux_13CO10=flux_mask_get(spiral,rms,chans,chan_width)[0]
uncertainty_13CO10=flux_mask_get(spiral,rms,chans,chan_width)[1] 

ratio_table['13CO10 flux']['spiral']=flux_13CO10
ratio_table['13CO10 uncertainty']['spiral']=uncertainty_13CO10

# anomaly 
maskimage='anomaly_mask.fits'
anomaly_masked=fits_import(maskimage)[1]
anomaly=masked_convert(data13_masked,anomaly_masked)

rms=0.0006
chan_width=20
flux_13CO10=flux_mask_get(anomaly,rms,chans,chan_width)[0]
uncertainty_13CO10=flux_mask_get(anomaly,rms,chans,chan_width)[1]

ratio_table['13CO10 flux']['anomaly']=flux_13CO10
ratio_table['13CO10 uncertainty']['anomaly']=uncertainty_13CO10

# rest disk 
maskimage='disk_mask.fits'
restdisk_masked=fits_import(maskimage)[1]
restdisk=masked_convert(data13_masked,restdisk_masked)

rms=0.0006
chan_width=20
flux_13CO10=flux_mask_get(restdisk,rms,chans,chan_width)[0]
uncertainty_13CO10=flux_mask_get(restdisk,rms,chans, chan_width)[1]

ratio_table['13CO10 flux']['restdisk']=flux_13CO10
ratio_table['13CO10 uncertainty']['restdisk']=uncertainty_13CO10

# # nonarm 
# maskimage='nonarm_mask.fits'
# nonarm_masked=fits_import(maskimage)[1]
# nonarm=masked_convert(data13_masked,nonarm_masked)

# rms=0.0006
# chan_width=20
# flux_13CO10=flux_mask_get(nonarm,rms,chans,chan_width)[0]
# uncertainty_13CO10=flux_mask_get(nonarm,rms,chans,chan_width)[1]

# ratio_table['13CO10 flux']['nonarm']=flux_13CO10
# ratio_table['13CO10 uncertainty']['nonarm']=uncertainty_13CO10

# center
maskimage='center_mask.fits'
center_masked=fits_import(maskimage)[1]
center=masked_convert(data13_masked,center_masked)

rms=0.0006
chan_width=20
flux_13CO10=flux_mask_get(center,rms,chans,chan_width)[0]
uncertainty_13CO10=flux_mask_get(center,rms,chans, chan_width)[1]

ratio_table['13CO10 flux']['center']=flux_13CO10
ratio_table['13CO10 uncertainty']['center']=uncertainty_13CO10

############################################################
# calculating the ratio

for region in regions:
    ratio_tmp=ratio_cal(ratio_table['12CO10 flux'][region],
                    ratio_table['13CO10 flux'][region],
                    ratio_table['12CO10 uncertainty'][region],
                    ratio_table['13CO10 uncertainty'][region])
    ratio_table['ratio'][region]=ratio_tmp[0]
    ratio_table['uncertainty'][region]=ratio_tmp[1]

############################################################
# draw the region on ratio map

# # draw the image of 12CO/13CO map. 
# fitsimage='NGC5257_1213_ratio.fits'
# wcs=fits_import(fitsimage)[0]
# data=fits_import(fitsimage)[1]

# levels=[0,0.5]
# colors='red'

# fig=plt.figure()
# ax=plt.subplot(projection=wcs)
# im=ax.imshow(data,cmap='rainbow',origin='lower',vmin=0,vmax=40)
# ax.contour(spiral_masked.data,levels=levels,transform=ax.get_transform(wcs),colors='red') 
# ax.contour(anomaly_masked.data,levels=levels,transform=ax.get_transform(wcs),colors='black')
# cbar=fig.colorbar(im)
# lon = ax.coords[0]
# lat = ax.coords[1]
# lon.set_major_formatter('hh:mm:ss')
# lon.set_ticks_visible(False)
# lon.set_ticklabel_visible(False)
# lat.set_ticks_visible(False)
# lat.set_ticklabel_visible(False)
# plt.savefig('../../picture/NGC5257_1213_region.png')

############################################################
# output the data frame result
os.chdir(scriptDir)

filename=outfile
with open(filename,'w') as out:
    out.write('from script NGC5257_measure_mod_error.txt \n')
    out.write('\n')
    out.write(comment)
    ratio_table.to_string(out)

# with open(filename,'a') as out:
#     out.write('\n')
#     out.write('\n')
#     out.write(comment)
#     ratio_table.to_string(out)
