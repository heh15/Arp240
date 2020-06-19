from subprocess import call
import time
from datetime import datetime
from astropy.io import fits
import glob
import os
from astropy.wcs import WCS

Dir='/1/home/heh15/workingspace/Arp240/NGC5257/ratio/'
imageDir=Dir+'image/'

fitsfile=imageDir+'12CO32/NGC5257co32_all.map40r.mom0.fits' 
outputfits=imageDir+'12CO32/NGC5257co32_all_map40r_shift.fits'
hdul=fits.open(fitsfile)
hdr=hdul[0].header
hdr['CRVAL1']=204.970621791
hdr['CRVAL2']=0.8399
hdul.writeto(outputfits, overwrite=True)
hdul.close()

# fitsfile=imageDir+'12CO32/NGC5257co32_all.map40r.mom0.fits'
# hdr=fits.open(fitsfile)[0].header
# wcs=WCS(hdr)

#### shift the coordinates for moment 0 map made by Hao. 

fitsfile=imageDir+'12CO32/NGC5257co32_all_map40r_2rms_mom0.fits' 
outputfits=imageDir+'12CO32/NGC5257co32_all_map40r_2rms_mom0_shift.fits'
hdul=fits.open(fitsfile)
hdr=hdul[0].header
hdr['CRVAL1']=204.970621791
hdr['CRVAL2']=0.8399
hdul.writeto(outputfits, overwrite=True)
hdul.close()

# fitsfile=imageDir+'12CO32/NGC5257co32_all.map40r.mom0.fits'
# hdr=fits.open(fitsfile)[0].header
# wcs=WCS(hdr)

#### shift the coordinateds for nchan map

fitsfile=imageDir+'12CO32/NGC5257co32_all_map40r_nchan.fits' 
outputfits=imageDir+'12CO32/NGC5257co32_all_map40r_nchan_shift.fits'
hdul=fits.open(fitsfile)
hdr=hdul[0].header
hdr['CRVAL1']=204.970621791
hdr['CRVAL2']=0.8399
hdul.writeto(outputfits, overwrite=True)
hdul.close()
