import cube
from astropy.utils import data
from spectral_cube import SpectralCube, Projection
from radio_beam import Beam
from astropy import units as u
import time
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits
from astropy.wcs import WCS


############################################################
# directory

Dir='/1/home/heh15/workingspace/Arp240/scatter/'
imageDir=Dir+'image/'
logDir=Dir+'log/'


############################################################

# basic setting
G=6.67e-11
pc=3.1e16
k=1.38e-23

alpha=4.3
ratio=0.77

incl=0.45
beammaj=1.004
beammin=0.556
freq=230.54


############################################################
# main program

name=imageDir+'NGC5257/NGC5257_12CO21_combine_sinbeam_cube.fits'
imagecube=SpectralCube.read(name)
Imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)


## create rms cube
rmscube=cube.calc_noise_in_cube(Imcube)

# mask the the low value of rmscube. 
mask=rmscube<3.0e-3*u.Jy/u.beam
lowrms=rmscube.with_mask(~mask)
newrms=lowrms.with_fill_value(3.0e-3)

## find the signal of the cube. 
outcube=cube.find_signal_in_cube(Imcube,newrms,snr_hi=5)
kcube=outcube.to(u.K)

## surface density
mom0=kcube.moment(order=0)
array=np.array(mom0)
array_inan=np.nan_to_num(array)
mom0_binned=array_inan.reshape(192,5,192,5).mean(-1).mean(1)
SD=alpha/ratio*mom0_binned*incl
SD_tmp=alpha/ratio*np.array(mom0)*incl

rms_mom0=5*0.17*incl
threshold=alpha/ratio*rms_mom0/(0.0109*beammaj*beammin*(freq/115.27)**2)
lowsd=SD_tmp<threshold
SD_tmp[lowsd]='nan'

## velocity dispersion
mom2=kcube.linewidth_sigma()
mom2_binned=mom2.reshape(192,5,192,5).mean(-1).mean(1)
disp=np.array(mom2_binned)

## scale height
height=disp**2/(SD*math.pi*G)*(1e6/(2e30/pc**2))
Height=height/pc

## gravitational pressure
P_grav=0.5*math.pi*G*(2e30/pc**2*SD)**2
P_grav=P_grav/(k*10**6)

Comparison between spectral cube processed data and casa processed data.
data=np.transpose(np.loadtxt(logDir+'Sd_H.txt'))
sd=data[0].astype(np.float)
height=data[1]

## surface density and scale height
fitsfile=imageDir+'NGC5257/sd_casa.fits'
sd=fits.open(fitsfile)[0].data

fitsfile=imageDir+'NGC5257/height_casa.fits'
height=fits.open(fitsfile)[0].data


fig=plt.figure()
plt.title('H-sd')
ax=plt.axes(xscale='log',yscale='log')
plt.ylim(10,10000)
plt.scatter(sd,height,marker='.',color='blue')
plt.scatter(SD,Height,marker='.',color='black')


# check the surface densiy less than 5*rms_mom0
rms_mom0=5*0.17*incl
threshold=alpha/ratio*rms_mom0/(0.0109*beammaj*beammin*(freq/115.27)**2)
lowsd=SD<threshold
SD[lowsd]='nan'
Height[lowsd]='nan'

fig=plt.figure()
ax=plt.axes(xscale='log',yscale='log')
plt.ylim(10,10000)
plt.scatter(SD,Height,marker='.',color='black')
plt.scatter(sd,height,marker='.')

fig=plt.figure()
plt.scatter(sd,SD)

fig=plt.figure()
# ax=plt.axes(xscale='log',yscale='log')
plt.scatter(height,Height,marker='.')
plt.xlabel('casa height')
plt.ylabel('script height')
plt.xlim(0,2000)
plt.ylim(0,2000)

fig=plt.figure()
plt.imshow(height,origin='lower')

## store the calculated scale height. 
outputfits=imageDir+'NGC5257/height_script.fits'
hdul=fits.open(imageDir+'NGC5257/NGC5257_12CO21_pbcor_rebin_mom0.fits')
hdul[0].data=Height
hdul.writeto(outputfits)
hdul.close()

outputfits=imageDir+'NGC5257/height_casa2.fits'
hdul=fits.open(imageDir+'NGC5257/NGC5257_12CO21_pbcor_rebin_mom0.fits')
hdul[0].data=height
hdul.writeto(outputfits)
hdul.close()
