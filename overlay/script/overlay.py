from subprocess import call
import time
from datetime import datetime
from astropy.io import fits
import glob
import os

Dir='/1/home/heh15/workingspace/Arp240/overlay/'
scriptDir=Dir+'script/'
imageDir=Dir+'image/'

fitsfile=imageDir+'NGC5258_12CO21_combine_noise45_mom0.fits'
hdul=fits.open(fitsfile)
hdr=hdul[0].header
ra=hdr['CRVAL1']
dec=hdr['CRVAL2']
hdul.close()

fitsfile=imageDir+'NGC5257_12CO21_combine_noise40_mom0.fits'
outputfits=imageDir+'NGC5257_12CO21_combine_noise40_mom0_header.fits'
hdul=fits.open(fitsfile)
hdr=hdul[0].header
hdr['CRVAL1']=ra
hdr['CRVAL2']=dec
hdul.writeto(outputfits)
hdul.close()
