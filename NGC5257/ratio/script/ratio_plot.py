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
regions=['spiral','anomaly','restdisk','center']
values=['12CO10 flux','12CO10 uncertainty','12CO21 flux', '12CO21 uncertainty', 'ratio','uncertainty']
ratio_table=pd.DataFrame(index=regions,columns=values)

Dir='/home/heh15/workingspace/Arp240/NGC5257/ratio/'
scriptDir=Dir+'script/measure_2110/'
workDir=Dir+'mask/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
ratioDir=Dir+'2110/contsub/'

chans=50

Im12CO10=ratioDir+'NGC5257_12CO10_combine_contsub_uvrange_smooth_masked_mom0.fits'
Im12CO21=ratioDir+'NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid10_masked_mom0.fits'
outfile='test.txt'

# Im12CO10=imageDir+'NGC5257_12CO10_combine_contsub_smooth_2rms_mom0.fits'
# Im12CO21=imageDir+'NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_2rms_mom0.fits'
# outfile='test.txt'
# comment='2rms cut moment map \n'

# Im12CO10=imageDir+'NGC5257_12CO10_combine_contsub_smooth_mom0.fits'
# # Im12CO21=imageDir+'NGC5257_12CO21_combine_contsub_uvtaper_smooth_regrid_mom0.fits'
# outfile='test.txt'
# comment='no cut moment map \n'

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
    ratio=flux2/flux1*115.27**2/230.54**2
    uncertainty=math.sqrt((uncertainty1/flux1)**2+(uncertainty2/flux2)**2+2*0.05**2)*ratio
    return ratio, uncertainty


