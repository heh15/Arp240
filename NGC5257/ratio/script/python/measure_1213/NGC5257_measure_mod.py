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
# basic setting

beam_area=2.186*1.896*1.1331
beam_area_pix=beam_area/0.09
regions=['spiral','anomaly','restdisk','nonarm']
values=['12CO10 flux','12CO10 uncertainty','13CO10 flux', '13CO10 uncertainty', 'ratio','uncertainty']
ratio_table=pd.DataFrame(index=regions,columns=values)

Dir='/home/heh15/workingspace/Arp240/ratio/'
scriptDir=Dir+'script/'
workDir=Dir+'NGC5257/mask/'

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

def flux_get(data_region,rms,chan_width):
    flux=np.ma.sum(data_region)/beam_area_pix
    error_value=rms*chan_width*sqrt(50/beam_area_pix)
    error = np.full((180,140),error_value)
    uncertainty=error_value*math.sqrt(np.ma.count(data_region))
    return flux, uncertainty

def ratio_cal(flux1,flux2,uncertainty1,uncertainty2):
    ratio=flux1/flux2
    uncertainty=math.sqrt((uncertainty1/flux1)**2+(uncertainty2/flux2)**2+2*0.05**2)*ratio
    return ratio, uncertainty


############################################################
os.chdir(workDir)

# calculating 12CO10 flux
fitsimage='NGC5257_12CO10_combine_smooth_masked_mom0.fits'
wcs=fits_import(fitsimage)[0]
data12_masked=fits_import(fitsimage)[1]

# spiral arm
maskimage='spiralarm_mask.fits'
spiral_masked=fits_import(maskimage)[1]
spiral= masked_convert(data12_masked,spiral_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_get(spiral,rms,chan_width)[0]
uncertainty_12CO10=flux_get(spiral,rms, chan_width)[1]
ratio_table['12CO10 flux']['spiral']=flux_12CO10
ratio_table['12CO10 uncertainty']['spiral']=uncertainty_12CO10


# anomaly 
maskimage='anomaly_mod_mask.fits'
anomaly_masked=fits_import(maskimage)[1]
anomaly=masked_convert(data12_masked,anomaly_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_get(anomaly,rms,chan_width)[0]
uncertainty_12CO10=flux_get(anomaly,rms, chan_width)[1]

ratio_table['12CO10 flux']['anomaly']=flux_12CO10
ratio_table['12CO10 uncertainty']['anomaly']=uncertainty_12CO10

# rest disk

maskimage='disk_mask.fits'
restdisk_masked=fits_import(maskimage)[1]
restdisk=masked_convert(data12_masked,restdisk_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_get(restdisk,rms,chan_width)[0]
uncertainty_12CO10=flux_get(restdisk,rms, chan_width)[1]

ratio_table['12CO10 flux']['restdisk']=flux_12CO10
ratio_table['12CO10 uncertainty']['restdisk']=uncertainty_12CO10

# nonarm 
maskimage='nonarm_mod_mask.fits'
nonarm_masked=fits_import(maskimage)[1]
nonarm=masked_convert(data12_masked,nonarm_masked)

rms=0.0016
chan_width=10
flux_12CO10=flux_get(nonarm,rms,chan_width)[0]
uncertainty_12CO10=flux_get(nonarm,rms, chan_width)[1]

ratio_table['12CO10 flux']['nonarm']=flux_12CO10
ratio_table['12CO10 uncertainty']['nonarm']=uncertainty_12CO10

############################################################
# calculating the flux of 13CO data 

fitsimage='NGC5257_13CO10_12m_smooth_masked_mom0.fits'
wcs=fits_import(fitsimage)[0]
data13_masked=fits_import(fitsimage)[1]

# spiral arm 
maskimage='spiralarm_mask.fits'
spiral_masked=fits_import(maskimage)[1]
spiral= masked_convert(data13_masked,spiral_masked)

rms=0.0006
chan_width=20
flux_13CO10=flux_get(spiral,rms,chan_width)[0]
uncertainty_13CO10=flux_get(spiral,rms, chan_width)[1] 

ratio_table['13CO10 flux']['spiral']=flux_13CO10
ratio_table['13CO10 uncertainty']['spiral']=uncertainty_13CO10

# anomaly 
maskimage='anomaly_mod_mask.fits'
anomaly_masked=fits_import(maskimage)[1]
anomaly=masked_convert(data13_masked,anomaly_masked)

rms=0.0006
chan_width=20
flux_13CO10=flux_get(anomaly,rms,chan_width)[0]
uncertainty_13CO10=flux_get(anomaly,rms, chan_width)[1]

ratio_table['13CO10 flux']['anomaly']=flux_13CO10
ratio_table['13CO10 uncertainty']['anomaly']=uncertainty_13CO10

# rest disk 
maskimage='disk_mask.fits'
restdisk_masked=fits_import(maskimage)[1]
restdisk=masked_convert(data13_masked,restdisk_masked)

rms=0.0006
chan_width=20
flux_13CO10=flux_get(restdisk,rms,chan_width)[0]
uncertainty_13CO10=flux_get(restdisk,rms, chan_width)[1]

ratio_table['13CO10 flux']['restdisk']=flux_13CO10
ratio_table['13CO10 uncertainty']['restdisk']=uncertainty_13CO10

# nonarm 
maskimage='nonarm_mod_mask.fits'
nonarm_masked=fits_import(maskimage)[1]
nonarm=masked_convert(data13_masked,nonarm_masked)

rms=0.0006
chan_width=20
flux_13CO10=flux_get(nonarm,rms,chan_width)[0]
uncertainty_13CO10=flux_get(nonarm,rms, chan_width)[1]

ratio_table['13CO10 flux']['nonarm']=flux_13CO10
ratio_table['13CO10 uncertainty']['nonarm']=uncertainty_13CO10


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
# output the data frame result
os.chdir(scriptDir)

ratio_table.to_csv(path_or_buf='../log/NGC5257_1213_ratio_mod1.csv')

