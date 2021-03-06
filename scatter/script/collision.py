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
import pandas as pd

from matplotlib import rcParams
rcParams['mathtext.default']='regular'

############################################################
# directory

Dir='/1/home/heh15/workingspace/Arp240/scatter/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
picDir=Dir+'picture/'
mapDir=Dir+'map/'

############################################################

# basic setting
G=6.67e-11
pc=3.1e16
k=1.38e-23

alpha=1.1
ratio=0.77
# incl=0.45
incl=1.0
beammaj=1.1
beammin=0.8
freq=225.46
freqSFR=33
rms_mom0=3.1e-3*10*math.sqrt(50)
threshold_CO=5*rms_mom0
# rms_mom0=5*0.17*incl
D=99

bin=10
threshold_33GHz=3*1.0e-5
subnum=53


############################################################
# function

def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)

    return wcs, data_masked

def sfr_radio(flux, beamajor,beaminor,freq,d):
    '''
    flux in Jy
    freq in Hz
    d in Mpc
    '''
    f_erg=flux*1e-23
    L=4*math.pi*(d*1e6*3.086e18)**2*f_erg
    SFR_tot=1e-27*(2.18*freq**-0.1+15.1*freq**-0.7)**(-1)*L

    return SFR_tot


############################################################
# main program

##### Data Extraction ##### 

Galaxies=['NGC5257','NGC5258']
quantities=['SD','vd','pressure','SFR','tdep','SD_err','vd_err']
Quantities=dict.fromkeys(Galaxies)
for key in Quantities.keys():
    Quantities[key]=dict.fromkeys(quantities)

#### NGC 5257

# name=imageDir+'NGC5257/NGC5257_12CO21_combine_smooth.fits'
# imagecube=SpectralCube.read(name)
# Imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)


# ## create rms cube
# rmscube=cube.calc_noise_in_cube(Imcube)

# # mask the the low value of rmscube. 
# mask=rmscube<3.0e-3*u.Jy/u.beam
# lowrms=rmscube.with_mask(~mask)
# newrms=lowrms.with_fill_value(3.0e-3)

# ## find the signal of the cube. 
# outcube=cube.find_signal_in_cube(Imcube,newrms,snr_hi=5)
# kcube=outcube.to(u.K)
# kcube.write(imageDir+'NGC5257/NGC5257_kcube.fits')

kcube=SpectralCube.read(imageDir+'NGC5257/12CO21/NGC5257_kcube.fits')

## linewidth in each pixel
cube_array=np.array(kcube)
count=(~np.isnan(cube_array)).sum(0)
linewidth=count*10

xpos=480;ypos=483 # center position
low=[int(xpos-(subnum-1)/2*bin-bin/2), int(ypos-(subnum-1)/2*bin-bin/2)]
high=[int(xpos+(subnum-1)/2*bin+bin/2),int(ypos+(subnum-1)/2*bin+bin/2)]
linewidth_sub=linewidth[low[1]:high[1],low[0]:high[0]]

linewidth_binned=linewidth_sub.reshape(int(subnum),bin,int(subnum),bin).mean(-1).mean(1)

## surface density
mom0=kcube.moment(order=0)

# exctract the subimage
xpos=480;ypos=483 # center position
low=[int(xpos-(subnum-1)/2*bin-bin/2), int(ypos-(subnum-1)/2*bin-bin/2)]
high=[int(xpos+(subnum-1)/2*bin+bin/2),int(ypos+(subnum-1)/2*bin+bin/2)]
mom0sub=mom0[low[1]:high[1],low[0]:high[0]]

# rebin the moment 0 map. 
array=np.array(mom0sub)
array_inan=np.nan_to_num(array)
mom0_binned=array_inan.reshape(int(subnum),bin,int(subnum),bin).mean(-1).mean(1)
# threshold=threshold_CO/(0.0109*beammaj*beammin*(freq/115.27)**2)
# lowvalue=mom0_binned<threshold
# mom0_binned[lowvalue]='nan'
# mom0_err=rms*10*np.sqrt(linewidth/10)
SD=alpha/ratio*mom0_binned*incl
Quantities['NGC5257']['SD']=SD


## velocity dispersion
mom2=kcube.linewidth_sigma()

# exctract the subimage
xpos=480;ypos=483 # center position
mom2sub=mom2[low[0]:high[0],low[1]:high[1]]

# bin the image. 
mom2_binned=mom2sub.reshape(int(subnum),bin,int(subnum),bin).mean(-1).mean(1)
disp=np.array(mom2_binned)
# disp_err=mom0_err/mom0_binned*linewidth_binned**2/(4*math.sqrt(10))/(2*disp)
Quantities['NGC5257']['vd']=disp
# Quantities['NGC5257']['vd_err']=disp_err

## scale height
height=disp**2/(SD*math.pi*G)*(1e6/(2e30/pc**2))
Height=height/pc

## gravitational pressure
P_grav=0.5*math.pi*G*(2e30/pc**2*SD)**2
P_grav=P_grav/(k*10**6)
Quantities['NGC5257']['pressure']=P_grav

### import the SFR in different regions.

fitsimage=imageDir+'NGC5257/33GHz/NGC5257_33GHz_pbcor_regrid_smooth.fits'
data=fits_import(fitsimage)[1].data
datasub=data[low[1]:high[1],low[0]:high[0]]
datasub_binned=np.nanmean(np.nanmean(datasub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

## filter out the low SFR pixel value

# calculate the rms noise in pbcor image. 
# pbimage=imageDir+'NGC5257/33GHz/NGC5257_33GHz_regrid_pb.fits'
# pbdata=fits_import(fitsimage)[1].data
# pbdatasub=pbdata[low[1]:high[1],low[0]:high[0]]
# pbdatasub_binned=np.nanmean(np.nanmean(pbdatasub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

threshold=threshold_33GHz
filter=datasub_binned<threshold
datasub_binned[filter]='nan'

I_erg=datasub_binned*1e-23/(beammaj*beammin*1.1331)*(3600*180/math.pi)**2
sig_SFR=1e-27*(2.18*freqSFR**-0.1+15.1*freqSFR**-0.7)**(-1)*4*math.pi*I_erg
Sig_SFR=sig_SFR*(3.1e18*1000)**2*incl

Quantities['NGC5257']['SFR']=Sig_SFR


## depletion time
tau=SD/Sig_SFR*1000**2
# tau_binned=np.nanmean(np.nanmean(tau.reshape(int(subnum),bin,int(subnum),bin),axis=-1),axis=1)
tau_binned=np.copy(tau)
Quantities['NGC5257']['tdep']=tau

##  calculate the orbital time efficiency. 
filename=mapDir+'NGC5257_Omega_map.fits'
Omega=fits_import(filename)[1]
Omega_sub=Omega[low[1]:high[1],low[0]:high[0]]
Omegasub_binned=np.nanmean(np.nanmean(Omega_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)
torb1=2*math.pi/Omegasub_binned*10**6
epsorb1=torb1/Quantities['NGC5257']['tdep']

##  import the Toomre factor and beta. 
filename=mapDir+'NGC5257_Toomre_map.fits'
Qtot=fits_import(filename)[1] 
Qtot_sub=Qtot[low[1]:high[1],low[0]:high[0]]
Qtotsub_binned1=np.nanmean(np.nanmean(Qtot_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

filename=mapDir+'NGC5257_beta_map.fits'
beta=fits_import(filename)[1] 
beta_sub=beta[low[1]:high[1],low[0]:high[0]]
betasub_binned1=np.nanmean(np.nanmean(beta_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

## import the the data. 

#### NGC 5258 (repeat:L145)

# name=imageDir+'NGC5258/NGC5258_12CO21_combine_noise45_smooth.fits'
# imagecube=SpectralCube.read(name)
# Imcube=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)

# ## create rms cube
# rmscube=cube.calc_noise_in_cube(Imcube)

# # mask the the low value of rmscube. 
# mask=rmscube<3.0e-3*u.Jy/u.beam
# lowrms=rmscube.with_mask(~mask)
# newrms=lowrms.with_fill_value(3.0e-3)

# ## find the signal of the cube. 
# outcube=cube.find_signal_in_cube(Imcube,newrms,snr_hi=5)
# kcube=outcube.to(u.K)
# kcube.write(imageDir+'NGC5258/NGC5258_kcube.fits')

## import the cube
kcube=SpectralCube.read(imageDir+'NGC5258/NGC5258_kcube.fits')
# mom0_Jy=outcube.moment(order=0)
# mom0_Jy.write(imageDir+'NGC5258/NGC5258_12CO21_cube_mom0.fits')

# ## surface density
# mom0=kcube.moment(order=0)
# # mom0.write(imageDir+'NGC5258/NGC5258_kcube_mom0.fits')

mom0=kcube.moment(order=0)

xpos=560;ypos=406 # center position
right=15;left=subnum-right-1
lower=15;upper=subnum-lower-1
low=[int(xpos-left*bin-bin/2), int(ypos-lower*bin-bin/2)]
high=[int(xpos+right*bin+bin/2),int(ypos+upper*bin+bin/2)]
mom0sub=mom0[low[1]:high[1],low[0]:high[0]] # xpos and ypos different in spectral cube? 

array=np.array(mom0sub)
array_inan=np.nan_to_num(array)
mom0_binned=array_inan.reshape(int(subnum),bin,int(subnum),bin).mean(-1).mean(1)
threshold=threshold_CO/(0.0109*beammaj*beammin*(freq/115.27)**2)
lowvalue=mom0_binned<threshold
mom0_binned[lowvalue]='nan'
SD=alpha/ratio*mom0_binned*incl
Quantities['NGC5258']['SD']=SD

## velocity dispersion
mom2=kcube.linewidth_sigma()
mom2sub=mom2[low[1]:high[1],low[0]:high[0]]

# create the subimage

mom2_binned=mom2sub.reshape(subnum,10,subnum,10).mean(-1).mean(1)
disp=np.array(mom2_binned)
Quantities['NGC5258']['vd']=disp

## scale height
height=disp**2/(SD*math.pi*G)*(1e6/(2e30/pc**2))
Height=height/pc

## gravitational pressure
P_grav=0.5*math.pi*G*(2e30/pc**2*SD)**2
P_grav=P_grav/(k*10**6)


### import the SFR in different regions.

fitsimage=imageDir+'NGC5258/NGC5258_33GHz_pbcor_regrid_smooth.fits'

## filter out the low SFR pixel value
data=fits_import(fitsimage)[1].data
datasub=data[low[1]:high[1],low[0]:high[0]]
datasub_binned=np.nanmean(np.nanmean(datasub.reshape(subnum,10,subnum,10),axis=-1),axis=1)
threshold=threshold_33GHz
filter=datasub_binned<threshold
datasub_binned[filter]='nan'

I_erg=datasub_binned*1e-23/(beammaj*beammin*1.1331)
L_erg=I_erg*4*math.pi*(D*1e6*3.1e18)**2
sig_SFR=1e-27*(2.18*freqSFR**-0.1+15.1*freqSFR**-0.7)**(-1)*L_erg
Sig_SFR=sig_SFR/(4.85*D/1000)**2*incl

Quantities['NGC5258']['SFR']=Sig_SFR

## depletion time
tau=SD/Sig_SFR*1000**2
tau_binned=np.copy(tau)
Quantities['NGC5258']['tdep']=tau

## orbital time efficiency
filename=mapDir+'NGC5258_Omega_map.fits'
Omega=fits_import(filename)[1]
Omega_sub=Omega[low[1]:high[1],low[0]:high[0]]
Omegasub_binned=np.nanmean(np.nanmean(Omega_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)
torb2=2*math.pi/Omegasub_binned*10**6
epsorb2=torb2/Quantities['NGC5258']['tdep']

## Toomre factor
filename=mapDir+'NGC5258_Toomre_map.fits'
Qtot=fits_import(filename)[1] 
Qtot_sub=Qtot[low[1]:high[1],low[0]:high[0]]
Qtotsub_binned2=np.nanmean(np.nanmean(Qtot_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

## beta
filename=mapDir+'NGC5258_beta_map.fits'
beta=fits_import(filename)[1] 
beta_sub=beta[low[1]:high[1],low[0]:high[0]]
betasub_binned2=np.nanmean(np.nanmean(beta_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

##### Analysis #####

# ## plot the efficiency per orbital time vs Q 
# fig=plt.figure()
# plt.yscale('log')

# epsorb1[np.isnan(Qtotsub_binned1)]='nan'
# epsorb2[np.isnan(Qtotsub_binned2)]='nan'

# sc1=plt.scatter(Qtotsub_binned1, epsorb1, c=betasub_binned1)
# sc2=plt.scatter(Qtotsub_binned2, epsorb2, c=betasub_binned2)
# cbar=plt.colorbar(sc1)

color1='red';color2='orange'

fig=plt.figure()
ax=fig.add_subplot(111)
plt.xscale('log')
plt.yscale('log')

tff2=math.sqrt(3)/(4*G)*(Quantities['NGC5258']['vd']*1000)/(Quantities['NGC5258']['SD']*2e30/pc**2)/(3600*24*365)
tff1=math.sqrt(3)/(4*G)*(Quantities['NGC5257']['vd']*1000)/(Quantities['NGC5257']['SD']*2e30/pc**2)/(3600*24*365)

tff2[np.isnan(Quantities['NGC5258']['tdep'])]='nan'
tff1[np.isnan(Quantities['NGC5257']['tdep'])]='nan'


plt.scatter(tff2, torb2, color=color2, marker='.', label='NGC5258')
plt.scatter(tff1, torb1, color=color1, marker='.', label='NGC5257')
plt.legend(loc='lower left')

lower=max(ax.set_xlim()[0], ax.set_ylim()[0])
upper=min(ax.set_xlim()[1], ax.set_ylim()[1])
ax.plot([lower, upper],[lower,upper],ls='--', color='black')
ax.text(3e7, 6e7, r'$t_{orb}=t_{ff}$', fontsize=20)
ax.set_xlabel(r'$t_{ff}$ (years)', fontsize=20)
ax.set_ylabel(r'$t_{orb} (years)$', fontsize=20)
fig.tight_layout()
plt.savefig(picDir+'torb_tff.png')
