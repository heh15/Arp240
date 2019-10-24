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
# measure the flux ratio of 12CO2-1/1-0

beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.09
regions=['spiral','anomaly','total','center', 'south']
values=['12CO10 flux','12CO10 uncertainty','12CO21 flux', '12CO21 uncertainty', 'ratio','uncertainty']
ratio_table=pd.DataFrame(index=regions,columns=values)

Dir='/home/heh15/workingspace/Arp240/NGC5257/ratio/'
scriptDir=Dir+'script/python/measure_2110/'
workDir=Dir+'mask/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
ratioDir=imageDir+'ratio/2110/uvtaper_mask/'
regionDir=Dir+'region/'

chans=50

Im12CO10=ratioDir+'NGC5257_12CO10_combine_contsub_uvrange_smooth_masked_pbcor_mom0.fits'
Im12CO21=ratioDir+'NGC5257_12CO21_combine_contsub_uvtaper_smooth_masked_pbcor_mom0.fits'
outfile='test.txt'

# Im12CO10=imageDir+'12CO10/NGC5257_12CO10_combine_contsub_uvrange_smooth_pbcor_mom0.fits'
# Im12CO21=imageDir+'12CO21/NGC5257_12CO21_combine_contsub_uvtaper_smooth_pbcor_mom0.fits'

# Im12CO10=imageDir+'NGC5257_12CO10_combine_contsub_smooth_2rms_mom0.fits'
# Im12CO21=imageDir+'NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_2rms_mom0.fits'
# outfile='test.txt'
# comment='2rms cut moment map \n'

# Im12CO10=imageDir+'NGC5257_12CO10_combine_contsub_smooth_mom0.fits'
# # Im12CO21=imageDir+'NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_mom0.fits'
# outfile='test.txt'
comment='no cut moment map \n'

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
    chans_tmp=chans+np.zeros((np.shape(data_region)[0],np.shape(data_region)[1]))
    error=np.sqrt(chans_tmp)*rms*chan_width/sqrt(beam_area_pix)
    error_masked=np.ma.masked_where(data_region.mask,error)
    uncertainty=math.sqrt(np.ma.sum(np.power(error_masked,2)))
    return flux, uncertainty

def ratio_cal(flux1,flux2,uncertainty1,uncertainty2):
    ratio=flux2/flux1*112.73**2/225.46**2
    uncertainty=math.sqrt((uncertainty1/flux1)**2+(uncertainty2/flux2)**2+2*0.05**2)*ratio
    return ratio, uncertainty

def Apmask_convert(aperture,data_cut):
    apmask=aperture.to_mask()
    shape=data_cut.shape
    mask=apmask.to_image(shape=((shape[0],shape[1])))
    ap_mask=mask==0
    ap_masked=np.ma.masked_where(ap_mask,data_cut)

    return ap_masked

##############################
# calculate the flux and uncertainty of 12CO1-0.
os.chdir(workDir)

# calculating 12CO10 flux
fitsimage=Im12CO10
# fitsimage=testDir+'NGC5257_12CO10_combine_contsub_uvrange_smooth_masked.image_mom0.fits'
wcs=fits_import(fitsimage)[0]
data10_masked=fits_import(fitsimage)[1]

# # import the channel used map. 
# fitsimage='NGC5257_12CO10_combine_smooth_nchan.fits'
# wcs=fits_import(fitsimage)[0]
# chans=fits_import(fitsimage)[1]
 
# spiral arm
maskimage='spiralarm_mask.fits'
spiral_masked=fits_import(maskimage)[1]
spiral= masked_convert(data10_masked,spiral_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(spiral,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(spiral,rms,chans, chan_width)[1]
ratio_table['12CO10 flux']['spiral']=flux_12CO10
ratio_table['12CO10 uncertainty']['spiral']=uncertainty_12CO10


# anomaly 
maskimage='anomaly_mask.fits'
anomaly_masked=fits_import(maskimage)[1]
anomaly=masked_convert(data10_masked,anomaly_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(anomaly,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(anomaly,rms,chans, chan_width)[1]

ratio_table['12CO10 flux']['anomaly']=flux_12CO10
ratio_table['12CO10 uncertainty']['anomaly']=uncertainty_12CO10


# # nonarm 
# maskimage='nonarm_mod_mask.fits'
# nonarm_masked=fits_import(maskimage)[1]
# nonarm=masked_convert(data21_masked,nonarm_masked)

# rms=0.00301
# chan_width=10
# flux_12CO10=flux_mask_get(nonarm,rms,chans,chan_width)[0]
# uncertainty_12CO10=flux_mask_get(nonarm,rms,chans, chan_width)[1]

# ratio_table['12CO10 flux']['nonarm']=flux_12CO10
# ratio_table['12CO10 uncertainty']['nonarm']=uncertainty_12CO10

# center
maskimage='center_mask.fits'
center_masked=fits_import(maskimage)[1]
center=masked_convert(data10_masked,center_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(center,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(center,rms,chans, chan_width)[1]

ratio_table['12CO10 flux']['center']=flux_12CO10
ratio_table['12CO10 uncertainty']['center']=uncertainty_12CO10

# south continuum source. 
from regions import read_ds9

file=regionDir+'south.reg'

south_sky=read_ds9(file,errors='warn')[0]
south_pix=south_sky.to_pixel(wcs)
south_masked=Apmask_convert(south_pix, data10_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(south_masked,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(south_masked,rms,chans, chan_width)[1]

ratio_table['12CO10 flux']['south']=flux_12CO10
ratio_table['12CO10 uncertainty']['south']=uncertainty_12CO10

# total

file=regionDir+'whole.reg'

whole_sky=read_ds9(file,errors='warn')[0]
whole_pix=whole_sky.to_pixel(wcs)
whole_masked=Apmask_convert(whole_pix, data10_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_mask_get(whole_masked,rms,chans,chan_width)[0]
uncertainty_12CO10=flux_mask_get(whole_masked,rms,chans,chan_width)[1]

ratio_table['12CO10 flux']['total']=flux_12CO10
ratio_table['12CO10 uncertainty']['total']=uncertainty_12CO10


##############################
# Calculate the flux and uncertainty of 12CO2-1

# calculating 12CO21 flux
fitsimage=Im12CO21

wcs=fits_import(fitsimage)[0]
data21_masked=fits_import(fitsimage)[1]

# # import the channel used map. 
# fitsimage='NGC5257_12CO21_combine_smooth_nchan.fits'
# wcs=fits_import(fitsimage)[0]
# chans=fits_import(fitsimage)[1]

# spiral arm
maskimage='spiralarm_mask.fits'
spiral_masked=fits_import(maskimage)[1]
spiral= masked_convert(data21_masked,spiral_masked)

rms=0.0046
chan_width=10
flux_12CO21=flux_mask_get(spiral,rms,chans,chan_width)[0]
uncertainty_12CO21=flux_mask_get(spiral,rms,chans, chan_width)[1]
ratio_table['12CO21 flux']['spiral']=flux_12CO21
ratio_table['12CO21 uncertainty']['spiral']=uncertainty_12CO21


# anomaly 
maskimage='anomaly_mask.fits'
anomaly_masked=fits_import(maskimage)[1]
anomaly=masked_convert(data21_masked,anomaly_masked)

rms=0.0046
chan_width=10
flux_12CO21=flux_mask_get(anomaly,rms,chans,chan_width)[0]
uncertainty_12CO21=flux_mask_get(anomaly,rms,chans, chan_width)[1]

ratio_table['12CO21 flux']['anomaly']=flux_12CO21
ratio_table['12CO21 uncertainty']['anomaly']=uncertainty_12CO21

# # nonarm 
# maskimage='nonarm_mod_mask.fits'
# nonarm_masked=fits_import(maskimage)[1]
# nonarm=masked_convert(data21_masked,nonarm_masked)

# rms=0.0046
# chan_width=10
# flux_12CO21=flux_mask_get(nonarm,rms,chans,chan_width)[0]
# uncertainty_12CO21=flux_mask_get(nonarm,rms,chans, chan_width)[1]

# ratio_table['12CO21 flux']['nonarm']=flux_12CO21
# ratio_table['12CO21 uncertainty']['nonarm']=uncertainty_12CO21

# center
maskimage='center_mask.fits'
center_masked=fits_import(maskimage)[1]
center=masked_convert(data21_masked,center_masked)

rms=0.0046
chan_width=10
flux_12CO21=flux_mask_get(center,rms,chans,chan_width)[0]
uncertainty_12CO21=flux_mask_get(center,rms,chans, chan_width)[1]

ratio_table['12CO21 flux']['center']=flux_12CO21
ratio_table['12CO21 uncertainty']['center']=uncertainty_12CO21

# south continuum source 

file=regionDir+'south.reg'

south_sky=read_ds9(file,errors='warn')[0]
south_pix=south_sky.to_pixel(wcs)
south_masked=Apmask_convert(south_pix, data21_masked)

rms=0.0046
chan_width=10
flux_12CO21=flux_mask_get(south_masked,rms,chans,chan_width)[0]
uncertainty_12CO21=flux_mask_get(south_masked,rms,chans, chan_width)[1]

ratio_table['12CO21 flux']['south']=flux_12CO21
ratio_table['12CO21 uncertainty']['south']=uncertainty_12CO21

# total
file=regionDir+'whole.reg'

whole_sky=read_ds9(file,errors='warn')[0]
whole_pix=whole_sky.to_pixel(wcs)
whole_masked=Apmask_convert(whole_pix, data21_masked)

rms=0.0046
chan_width=10
flux_12CO21=flux_mask_get(whole_masked,rms,chans,chan_width)[0]
uncertainty_12CO21=flux_mask_get(whole_masked,rms,chans,chan_width)[1]

ratio_table['12CO21 flux']['total']=flux_12CO21
ratio_table['12CO21 uncertainty']['total']=uncertainty_12CO21



############################################################
# calculating the ratio

for region in regions:
    ratio_tmp=ratio_cal(ratio_table['12CO10 flux'][region],
                    ratio_table['12CO21 flux'][region],
                    ratio_table['12CO10 uncertainty'][region],
                    ratio_table['12CO21 uncertainty'][region])
    ratio_table['ratio'][region]=ratio_tmp[0]
    ratio_table['uncertainty'][region]=ratio_tmp[1]

##############################
# record the value into the file.

os.chdir(scriptDir)

# with open(outfile,'w') as out:
#     out.write('from script NGC5257_measure_2110_mod_error.txt \n')
#     out.write('\n')
#     ratio_table.to_string(out)

with open(outfile,'a') as out:
    out.write('\n')
    out.write('\n')
    out.write(comment)
    ratio_table.to_string(out)