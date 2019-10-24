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
from scipy.stats.stats import pearsonr

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

## import Toomre factor and beta. 
filename=mapDir+'NGC5257_Toomre_map.fits'
Qtot=fits_import(filename)[1] 
Qtot_sub=Qtot[low[1]:high[1],low[0]:high[0]]
Qtotsub_binned1=np.nanmean(np.nanmean(Qtot_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

filename=mapDir+'NGC5257_beta_map.fits'
beta=fits_import(filename)[1] 
beta_sub=beta[low[1]:high[1],low[0]:high[0]]
betasub_binned1=np.nanmean(np.nanmean(beta_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

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


## import Toomre factor
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

#### Add the result from Chris k-s law table

filename=logDir+'datafile_final.txt'
data_lit=pd.read_csv(filename,header=0,sep=r"\s*",skiprows=29)

### Add results from aperture measurements

# # 12CO 1-0 measurement
# SFR_center=0.58;SD_center=163;t_center=SD_center/SFR_center*1000**2
# SFR_south=1.1; SD_south=28.6; t_south=SD_south/SFR_south*1000**2
# SFR_belt=0.38; SD_belt=51.0; t_belt=SD_belt/SFR_belt*1000**2
# SFR_southarm=0.79; SD_southarm=133 ; t_southarm=SD_southarm/SFR_southarm*1000**2

# 12CO 2-1 measurement

SFR_center=2.06;SD_center=607.574;t_center=SD_center/SFR_center*1000**2
SFR_south=2.41; SD_south=86.86; t_south=SD_south/SFR_south*1000**2
SFR_arm=0.91; SD_arm=94.947; t_arm=SD_arm/SFR_arm*1000**2
SFR_west=0.66; SD_west=50.081; t_west=SD_west/SFR_west*1000**2
SFR_southarm=1.94; SD_southarm= 336; t_southarm=SD_southarm/SFR_southarm*1000**2

#### Results for Arp 240. 

## free-fall time
tff2=math.sqrt(3)/(4*G)*(Quantities['NGC5258']['vd']*1000)/(Quantities['NGC5258']['SD']*2e30/pc**2)/(3600*24*365)
tff2=math.sqrt(3)/(4*math.sqrt(G))*math.sqrt((10/0.7)**2*pc)/np.sqrt(Quantities['NGC5258']['SD']/0.1*2e30/pc**2)/(3600*24*365)
epsiff2=tff2/tau_binned

tff1=math.sqrt(3)/(4*G)*(Quantities['NGC5257']['vd']*1000)/(Quantities['NGC5257']['SD']*2e30/pc**2)/(3600*24*365)
tff1=math.sqrt(3)/(4*math.sqrt(G))*math.sqrt((10/0.7)**2*pc)/np.sqrt(Quantities['NGC5257']['SD']/0.1*2e30/pc**2)/(3600*24*365)
epsiff1=tff1/Quantities['NGC5257']['tdep']

mask=np.where(epsiff2=='nan')
Quantities['NGC5257']['vd'][mask]='nan'
# fig=plt.figure()
# plt.imshow(Quantities['NGC5257']['vd'],origin='lower')


color1='red';color2='orange'

fig=plt.figure()
ax=plt.subplot()
# plt.title('ULIRG conversion factor, original equation')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\Sigma_{mol}$ ($M_{\odot}/pc^2$)', fontsize=20)
plt.ylabel('SFR ($M_{\odot}/pc^2/year$)', fontsize=20)
# plt.scatter(Quantities['NGC5257']['SD'],Quantities['NGC5257']['SFR'],marker='.',label='NGC 5257',color=color1)
# plt.scatter(Quantities['NGC5258']['SD'],Quantities['NGC5258']['SFR'],marker='.', label='NGC 5258',color=color2)
# plt.legend()
# plt.savefig(picDir+'ks_test.png')

# ks_x=np.linspace(100, 1000, 9)
# ks_y=ks_x*8e-4
# plt.plot(ks_x, ks_y)

molsd=10**data_lit['log-MMol']
SFRsd=10**data_lit['log-SFR']
plt.scatter(molsd,SFRsd,label='Wilson et al. 2019',marker='.')
plt.scatter(Quantities['NGC5257']['SD'],Quantities['NGC5257']['SFR'],marker='.',label='NGC 5257',color=color1)
plt.scatter(Quantities['NGC5258']['SD'],Quantities['NGC5258']['SFR'],marker='.', label='NGC 5258',color=color2)
plt.scatter(SD_center, SFR_center, marker='o',color='red', label='NGC5257 center')
plt.scatter(SD_south, SFR_south, marker='o', color='darkgreen', label='NGC5257 south')
plt.scatter(SD_west, SFR_west, marker='o', color='navy', label='NGC5257 west')
plt.scatter(SD_arm, SFR_arm, marker='o', color='maroon', label='NGC5257 arm')
plt.scatter(SD_southarm, SFR_southarm, marker='o', color='black', label='NGC5258 southarm')
plt.legend()
# plt.title('ULIRG conversion factor, alpha=1.1')
ax.tick_params(labelsize=18)
fig.tight_layout()
plt.savefig(picDir+'ks-test_compare_binned_poster.png', bbox_inches='tight', pad_inches=0.5)


fig=plt.figure()
ax=plt.subplot()
plt.ylabel('Depletion time (years)',fontsize=20)
plt.xlabel('$\Sigma_{mol} (M_{\odot}/pc^2)$', fontsize=20)
plt.xscale('log')
plt.yscale('log')
# plt.title('ULIRG conversion factor, alpha=1.1')

molsd=10**data_lit['log-MMol']
tdep=10**data_lit['log-tdep']
plt.scatter(molsd, tdep, marker='.',label='Wilson et al. 2019')

plt.scatter(SD,tau_binned,marker='.',color=color2,label='NGC 5258')
plt.scatter(Quantities['NGC5257']['SD'],Quantities['NGC5257']['tdep'],label='NGC 5257',marker='.',color=color1)
# plt.ylim(0,1e9)
plt.scatter(SD_center, t_center, marker='o',color='red',label='NGC 5257 center')
plt.scatter(SD_south, t_south, marker='o', color='darkgreen', label='NGC5257 south')
plt.scatter(SD_west, t_west, marker='o', color='navy', label='NGC5257 west')
plt.scatter(SD_arm, t_arm, marker='o', color='maroon', label='NGC5257 arm')
plt.scatter(SD_southarm, t_southarm, marker='o', color='black', label='NGC5258 southarm')
plt.legend()
ax.tick_params(labelsize=18)
fig.tight_layout()
plt.savefig(picDir+'tdep_SD_binned_poster.png')


fig=plt.figure()
ax=plt.subplot()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('$\epsilon_{ff}$', fontsize=20)
plt.xlabel('$\Sigma_{mol} (M_{\odot}/pc^2)$', fontsize=20)

eff=10**data_lit['log-eff']
plt.scatter(molsd,eff,marker='.',label='Wilson et al. 2019')
plt.scatter(Quantities['NGC5258']['SD'],epsiff2,marker='.',label='NGC5258',color=color2)
plt.scatter(Quantities['NGC5257']['SD'],epsiff1,marker='.',label='NGC5257',color=color1)
plt.legend()
ax.tick_params(labelsize=18)
# plt.title('ULIRG conversion factor, filling factor=0.1')
fig.tight_layout()
plt.savefig(picDir+'SFE_perfreefall_binned_poster_GMC.png')


### Qtotsub_binned1, epsiff1; Qtotsub_binned2, epsiff2. correlation. 

Qwhole=np.hstack((Qtotsub_binned1.flatten(), Qtotsub_binned2.flatten()))
epsiff_whole=np.hstack((epsiff1.flatten(), epsiff2.flatten()))

bad=~np.logical_or(np.isnan(Qwhole), np.isnan(epsiff_whole))
Qwhole_filter=np.compress(bad, Qwhole)
epsiff_whole_filter=np.compress(bad, epsiff_whole)

pearsonr(Qwhole_filter, epsiff_whole_filter)

fig=plt.figure()
plt.yscale('log')
plt.scatter(Qwhole_filter, epsiff_whole_filter)

fig=plt.figure()
plt.yscale('log')
plt.ylabel('$\epsilon_{ff}$')
plt.xlabel('Toomre factor Q')
sc=plt.scatter(Qtotsub_binned1, epsiff1,c=betasub_binned1,  marker='.', label='NGC5257', cmap='brg_r')
plt.scatter(Qtotsub_binned2, epsiff2,c=betasub_binned2,  marker='.', label='NGC5258', cmap='brg_r')
# plt.legend()
cbar=plt.colorbar(sc)
cbar.set_label(r'$ \beta $', rotation='vertical')
plt.savefig(picDir+'perfreefall_Q_beta.png')

