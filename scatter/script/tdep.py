import cube
from astropy.utils import data
from spectral_cube import SpectralCube, Projection
from radio_beam import Beam
from astropy import units as u
import time
import numpy as npl
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
# basic information

G=6.67e-11
pc=3.1e16
k=1.38e-23

alpha=1.1
ratio=0.85
# incl=0.45
incl=1.0
beammaj=1.1
beammin=0.8
freq=225.46
freqSFR=33

rms_CO=3.1e-3/(0.0109*beammaj*beammin*(225.46/115.27)**2)
rms_mom0=rms_CO*10*math.sqrt(50)
threshold_CO=5*rms_mom0
# rms_mom0=5*0.17*incl
D=99

chan_width=10
chans=50

bin=10
rms_33GHz=1.0e-5
threshold_33GHz=3*1.0e-5
subnum=53

############################################################
# basic settings. 

Galaxies=['NGC5257','NGC5258']
quantities=['SD','vd','pressure','SFR','tdep','tff']
quantities2=['value', 'error']

quantities3=['galaxy','SD','SD_err','vd','vd_err','SFR', 'SFR_err', 'pressure','tdep']
Arp240=pd.DataFrame(columns=quantities3)
NGC5257=pd.DataFrame(columns=quantities3)
NGC5258=pd.DataFrame(columns=quantities3)

Quantities=dict.fromkeys(Galaxies)
for key in Quantities.keys():
    Quantities[key]=dict.fromkeys(quantities)
    for key2 in Quantities[key]:
        Quantities[key][key2]=dict.fromkeys(quantities2)

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

def mom0_uncertainty(rms, chans, chan_width):

    sigmaI=rms*chan_width*np.sqrt(chans)

    return sigmaI

def mom2_uncertainty(sigmaI, mom0, mom2, linewidth):
    SigmaD=sigmaI/mom0*linewidth**2/(4*math.sqrt(10))/(2*mom2)

    return SigmaD

# def mom2_uncertainty(imcube, rms, chans=0,corr=None):
#     if chans==0:
#         chans=np.shape(imcube[:,0,0])[0]
#     mom0=imcube.moment(order=0)
#     mom2=imcube.linewidth_sigma()
#     cubenp=np.array(imcube)
#     cubenp[cubenp<4*rms]='nan'
#     chan_width=(imcube.spectral_axis[1]-imcube.spectral_axis[0]).value
#     chan_width=round(chan_width,1)
#     sigmaI=rms*deltaV*np.sqrt(chans)
#     linewidth=chan_width*chans
#     sigma_dis1=sigmaI/np.array(mom0)*linewidth**2/(4*math.sqrt(10))/(2*np.array(mom2))
#     err1=sigmaI/mom0_array
#     err2=sigma_dis1/mom2_array
#     # calculate the correlation coefficient
#     if corr==None:
#         corr=np.ndarray(shape=(np.shape(cubenp)[1],np.shape(cubenp)[2]))
#         for i in range(np.shape(cubenp)[1]):
#             for j in range(np.shape(cubenp)[2]):
#                 dispers=(vel- mom1_array[i][j])**2*cubenp[:,i,j]
#                 df=pd.DataFrame(np.transpose(np.vstack((dispers, cubenp[:,i,j]))))
#                 corr[i][j]=df.corr()[0][1]
    
#     sigma_dis2=np.array(mom2)*np.sqrt(err1**2+err2**2-2*corr*err1*err2)

#     return sigma_dis1, sigma_dis2, corr


def reproj_binning(data, wcs, bin_num):
    
    map_in_shape=np.shape(data)
    nx_in, ny_in=map_in_shape
    nx_out=math.trunc(nx_in/bin_num);ny_out=math.trunc(ny_in/bin_num)
    xs,ys=np.meshgrid(np.arange(nx_out), np.arange(ny_out))
    wcs_out=wcs.deepcopy()
    wcs_out.wcs.crpix =[math.trunc(nx_out/2), math.trunc(ny_out/2)]
    wcs_out.wcs.cdelt=wcs.wcs.cdelt*bin_num
    wcs_out.wcs.ctype = ['RA---SIN', 'DEC--SIN']
    coords_out=pixel_to_skycoord(xs, ys, wcs_out)
    coords_out_flat=coords_out.flatten()
    pixel_labels_out = np.arange(xs.size)
    data_binned=np.zeros((nx_out, ny_out)).flatten()
    map_out_shape=(nx_out, ny_out)
    
    xs_in, ys_in = np.meshgrid(np.arange(nx_in), np.arange(ny_in))
    coords_in = pixel_to_skycoord(xs_in, ys_in, wcs)
    pixel_map_arr = np.full((nx_in, ny_in), np.nan).flatten()

    i_in=0
    npix_in = coords_in.flatten().size
    dra, ddec = np.zeros(npix_in), np.zeros(npix_in)
    i_out, d2d, d3d = match_coordinates_sky(coords_in.flatten(), coords_out_flat)
    dra, ddec = (coords_in.flatten()).spherical_offsets_to(
        coords_out_flat[i_out])
    dra = dra.arcsec
    ddec = ddec.arcsec

    good = (-0.5001 <= dra) & (dra < 0.5001) & (-0.5001 <= ddec) & (ddec < 0.5001)
    pixel_map_arr[good]=pixel_labels_out[i_out[good]]
    data_labeled=np.stack((data.flatten(),pixel_map_arr), axis=1)
    nan_index=np.where(np.isnan(data_labeled[:,1]))
    data_labeled=np.delete(data_labeled, nan_index,axis=0)
    data_labeled=data_labeled[np.argsort(data_labeled[:,1])]
    data_group=np.split(data_labeled[:,0], np.cumsum(np.unique(data_labeled[:,1], return_counts=True)[1])[:-1])
    for i in pixel_labels_out:
        data_binned[i]=np.nanmean(data_group[i])

    data_binned=data_binned.reshape(map_out_shape)
        
    return wcs_out, data_binned


############################################################
# main program

##### Data Extraction ##### 

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

# test
fig=plt.figure()
plt.imshow(linewidth_binned, origin='lower')

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

chan_width=10
chans_binned=linewidth_binned/chan_width
mom0_error=mom0_uncertainty(rms_CO, chans_binned, chan_width)

SD=alpha/ratio*mom0_binned*incl
Quantities['NGC5257']['SD']['value']=SD
Quantities['NGC5257']['SD']['error']=mom0_error/mom0_binned*SD

## velocity dispersion
mom2=kcube.linewidth_sigma()

# exctract the subimage
xpos=480;ypos=483 # center position
mom2sub=mom2[low[0]:high[0],low[1]:high[1]]

# bin the image. 
mom2_binned=mom2sub.reshape(int(subnum),bin,int(subnum),bin).mean(-1).mean(1)
disp=np.array(mom2_binned)
# disp_err=mom0_err/mom0_binned*linewidth_binned**2/(4*math.sqrt(10))/(2*disp)
Quantities['NGC5257']['vd']['value']=disp
# Quantities['NGC5257']['vd_err']['value']=disp_err
# linewidth_binned=np.full((53,53),500)
Quantities['NGC5257']['vd']['error']=mom2_uncertainty(mom0_error, np.array(mom0_binned), np.array(mom2_binned), linewidth_binned)

## scale height
height=disp**2/(SD*math.pi*G)*(1e6/(2e30/pc**2))
Height=height/pc

## gravitational pressure
P_grav=0.5*math.pi*G*(2e30/pc**2*SD)**2
P_grav=P_grav/(k*10**6)
Quantities['NGC5257']['pressure']['value']=P_grav

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

Quantities['NGC5257']['SFR']['value']=Sig_SFR
Quantities['NGC5257']['SFR']['error']=rms_33GHz/datasub_binned*Sig_SFR


## depletion time
tau=SD/Sig_SFR*1000**2
# tau_binned=np.nanmean(np.nanmean(tau.reshape(int(subnum),bin,int(subnum),bin),axis=-1),axis=1)
tau_binned=np.copy(tau)
Quantities['NGC5257']['tdep']['value']=tau
Quantities['NGC5257']['tdep']['error']=tau*np.sqrt((Quantities['NGC5257']['SFR']['error']/Sig_SFR)**2+\
                                                   (Quantities['NGC5257']['SD']['error']/SD)**2)


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

## linewidth in each pixel
cube_array=np.array(kcube)
count=(~np.isnan(cube_array)).sum(0)
linewidth=count*10

xpos=560;ypos=406 # center position
right=15;left=subnum-right-1
lower=15;upper=subnum-lower-1
low=[int(xpos-left*bin-bin/2), int(ypos-lower*bin-bin/2)]
high=[int(xpos+right*bin+bin/2),int(ypos+upper*bin+bin/2)]
linewidth_sub=linewidth[low[1]:high[1],low[0]:high[0]]

linewidth_binned=linewidth_sub.reshape(int(subnum),bin,int(subnum),bin).mean(-1).mean(1)

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


chan_width=10
chans_binned=linewidth_binned/chan_width
mom0_error=mom0_uncertainty(rms_CO, chans_binned, chan_width)

SD=alpha/ratio*mom0_binned*incl
Quantities['NGC5258']['SD']['value']=SD
Quantities['NGC5258']['SD']['error']=mom0_error/mom0_binned*SD

## velocity dispersion
mom2=kcube.linewidth_sigma()
mom2sub=mom2[low[1]:high[1],low[0]:high[0]]

# create the subimage

mom2_binned=mom2sub.reshape(subnum,10,subnum,10).mean(-1).mean(1)
disp=np.array(mom2_binned)
Quantities['NGC5258']['vd']['value']=disp
Quantities['NGC5258']['vd']['error']=mom2_uncertainty(mom0_error, np.array(mom0_binned), np.array(mom2_binned), linewidth_binned)

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

Quantities['NGC5258']['SFR']['value']=Sig_SFR
Quantities['NGC5258']['SFR']['error']=rms_33GHz/datasub_binned*Sig_SFR

## depletion time
tau=SD/Sig_SFR*1000**2
tau_binned=np.copy(tau)
Quantities['NGC5258']['tdep']['value']=tau

Quantities['NGC5258']['tdep']['error']=tau*np.sqrt((Quantities['NGC5258']['SFR']['error']/Sig_SFR)**2+\
                                                   (Quantities['NGC5258']['SD']['error']/SD)**2)

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
regions=['center', 'south', 'arm', 'west', 'southarm']
properties=['SFR', 'SFR_err', 'SD', 'SD_err', 'tdep','tdep_err', 'vd','vd_err', 'tff','tff_err', 'eff', 'eff_err', 'fraction']
table=pd.DataFrame(columns=properties, index=regions)

# 12CO 2-1 measurement

table['SFR']['center']=2.06;table['SD']['center']=607.574;table['tdep']['center']=table['SD']['center']/table['SFR']['center']*1000**2
table['SFR']['south']=2.41; table['SD']['south']=86.86; table['tdep']['south']=table['SD']['south']/table['SFR']['south']*1000**2
table['SFR']['arm']=0.91; table['SD']['arm']=94.947; table['tdep']['arm']=table['SD']['arm']/table['SFR']['arm']*1000**2
table['SFR']['west']=0.66; table['SD']['west']=50.081; table['tdep']['west']=table['SD']['west']/table['SFR']['west']*1000**2
table['SFR']['southarm']=1.94; table['SD']['southarm']= 336; table['tdep']['southarm']=table['SD']['southarm']/table['SFR']['southarm']*1000**2

filename=logDir+'NGC5257_aperture_output.csv'
Aperture1_property=pd.read_csv(filename, index_col=0)
filename=logDir+'NGC5258_aperture_output.csv'
Aperture2_property=pd.read_csv(filename, index_col=0)
Aperture_property=pd.concat([Aperture1_property, Aperture2_property])

table['SFR_err']=Aperture_property['SFR_err']
table['SD_err']=Aperture_property['SD_err']
table['vd']=Aperture_property['dispersion']
table['vd_err']=Aperture_property['dispersion_err']
table['tff']=math.sqrt(3)/(4*G)*(table['vd']*1000)/(table['SD']*2e30/pc**2)/(3600*24*365)
table['tff_err']=table['vd_err']/table['vd']*table['tff']

table['tdep_err']=table['tdep']*((table['SFR_err']/table['SFR'])**2+(table['SD_err']/table['SD'])**2)**0.5
table['eff']=table['tff']/table['tdep']
table['eff_err']=table['eff']*((table['tdep_err']/table['tdep'])**2+(table['tff_err']/table['tff'])**2)**0.5

#### Results for Arp 240. 

## free-fall time
tff2=math.sqrt(3)/(4*G)*(Quantities['NGC5258']['vd']['value']*1000)/(Quantities['NGC5258']['SD']['value']*2e30/pc**2)/(3600*24*365)
Quantities['NGC5258']['tff']['value']=tff2
Quantities['NGC5258']['tff']['error']=tff2*np.sqrt((Quantities['NGC5258']['vd']['error']/Quantities['NGC5258']['vd']['value'])**2+
                                              (Quantities['NGC5258']['SD']['error']/Quantities['NGC5258']['SD']['value'])**2)

# tff2=math.sqrt(3)/(4*math.sqrt(G))*math.sqrt((10/0.7)**2*pc)/np.sqrt(Quantities['NGC5258']['SD']['value']/0.1*2e30/pc**2)/(3600*24*365)
epsiff2=tff2/tau_binned
epsiff2_err=epsiff2*np.sqrt((Quantities['NGC5258']['tff']['error']/Quantities['NGC5258']['tff']['value'])**2+
             (Quantities['NGC5258']['tdep']['error']/Quantities['NGC5258']['tdep']['value'])**2)



tff1=math.sqrt(3)/(4*G)*(Quantities['NGC5257']['vd']['value']*1000)/(Quantities['NGC5257']['SD']['value']*2e30/pc**2)/(3600*24*365)
Quantities['NGC5257']['tff']['value']=tff1
Quantities['NGC5257']['tff']['error']=tff1*np.sqrt((Quantities['NGC5257']['vd']['error']/Quantities['NGC5257']['vd']['value'])**2+
                                              (Quantities['NGC5257']['SD']['error']/Quantities['NGC5257']['SD']['value'])**2)
# tff1=math.sqrt(3)/(4*math.sqrt(G))*math.sqrt((10/0.7)**2*pc)/np.sqrt(Quan tities['NGC5257']['SD']/0.1*2e30/pc**2)/(3600*24*365)
epsiff1=tff1/Quantities['NGC5257']['tdep']['value']
epsiff1_err=epsiff1*np.sqrt((Quantities['NGC5257']['tff']['error']/Quantities['NGC5257']['tff']['value'])**2+
             (Quantities['NGC5257']['tdep']['error']/Quantities['NGC5257']['tdep']['value'])**2)

mask=np.where(epsiff2=='nan')
Quantities['NGC5257']['vd']['value'][mask]='nan'
# fig=plt.figure()
# plt.imshow(Quantities['NGC5257']['vd']['value'],origin='lower')


color1='red';color2='orange'

fig=plt.figure()
ax=plt.subplot()
# plt.title('ULIRG conversion factor, original equation')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\Sigma_{mol}$ ($M_{\odot}/pc^2$)', fontsize=20)
plt.ylabel('$\Sigma_{SFR}$ ($M_{\odot}/pc^2/year$)', fontsize=20)
# plt.scatter(Quantities['NGC5257']['SD']['value'],Quantities['NGC5257']['SFR']['value'],marker='.',label='NGC 5257',color=color1)
# plt.scatter(Quantities['NGC5258']['SD']['value'],Quantities['NGC5258']['SFR']['value'],marker='.', label='NGC 5258',color=color2)
# plt.legend()
# plt.savefig(picDir+'ks_test.png')

# ks_x=np.linspace(100, 1000, 9)
# ks_y=ks_x*8e-4
# plt.plot(ks_x, ks_y)

molsd=10**data_lit['log-MMol']
SFRsd=10**data_lit['log-SFR']
SFRerr=SFRsd*np.log(10)*data_lit['u_log-SFR']
plt.errorbar(molsd,SFRsd, SFRerr, label='Wilson et al. 2019',marker='.', linestyle='none')
plt.errorbar(Quantities['NGC5257']['SD']['value'].flatten(),Quantities['NGC5257']['SFR']['value'].flatten(),\
             xerr=Quantities['NGC5257']['SD']['error'].flatten(), yerr=Quantities['NGC5257']['SFR']['error'].flatten(), \
             marker='.',label='NGC 5257',color=color1, linestyle='none')
plt.errorbar(Quantities['NGC5258']['SD']['value'].flatten(),Quantities['NGC5258']['SFR']['value'].flatten(),\
             xerr=Quantities['NGC5258']['SD']['error'].flatten(), yerr=Quantities['NGC5258']['SFR']['error'].flatten(), \
             marker='.',label='NGC 5258',color=color2, linestyle='none')
plt.errorbar(table['SD']['center'], table['SFR']['center'],yerr=table['SFR_err']['center'],  marker='o',color='red', label='NGC5257 center')
plt.errorbar(table['SD']['south'], table['SFR']['south'], yerr=table['SFR_err']['south'], marker='o', color='darkgreen', label='NGC5257 south')
plt.errorbar(table['SD']['west'], table['SFR']['west'], yerr=table['SFR_err']['west'],  marker='o', color='navy', label='NGC5257 west')
plt.errorbar(table['SD']['arm'], table['SFR']['arm'], yerr=table['SFR_err']['arm'], marker='o', color='maroon', label='NGC5257 arm')
plt.errorbar(table['SD']['southarm'], table['SFR']['southarm'], yerr=table['SFR_err']['southarm'], marker='o', color='black', label='NGC5258 southarm')
plt.legend()
# plt.title('ULIRG conversion factor, alpha=1.1')
ax.tick_params(labelsize=18)
fig.tight_layout()
plt.savefig(picDir+'ks-test_compare_binned_poster.png', bbox_inches='tight', pad_inches=0.5)
 
tdep=10**data_lit['log-tdep']
tdeperr=tdep*np.log(10)*data_lit['u_log-tdep']

fig=plt.figure()
ax=plt.subplot()
plt.ylabel('Depletion time (years)',fontsize=20)
plt.xlabel('$\Sigma_{mol} (M_{\odot}/pc^2)$', fontsize=20)
plt.xscale('log')
plt.yscale('log')
# plt.title('ULIRG conversion factor, alpha=1.1')

molsd=10**data_lit['log-MMol']
tdep=10**data_lit['log-tdep']
plt.errorbar(molsd, tdep,yerr=tdeperr, marker='.', linestyle='none', label='Wilson et al. 2019')

plt.errorbar(SD.flatten(),tau_binned.flatten(),Quantities['NGC5258']['tdep']['error'].flatten(), linestyle='none', marker='.',color=color2,label='NGC 5258')
plt.errorbar(Quantities['NGC5257']['SD']['value'].flatten(),Quantities['NGC5257']['tdep']['value'].flatten(),
            Quantities['NGC5257']['tdep']['error'].flatten(), label='NGC 5257',marker='.',linestyle='none', color=color1)
# plt.ylim(0,1e9)
plt.errorbar(table['SD']['center'], table['tdep']['center'],yerr=table['tdep_err']['center'],  marker='o',color='red', label='NGC5257 center')
plt.errorbar(table['SD']['south'], table['tdep']['south'], yerr=table['tdep_err']['south'], marker='o', color='darkgreen', label='NGC5257 south')
plt.errorbar(table['SD']['west'], table['tdep']['west'], yerr=table['tdep_err']['west'],  marker='o', color='navy', label='NGC5257 west')
plt.errorbar(table['SD']['arm'], table['tdep']['arm'], yerr=table['tdep_err']['arm'], marker='o', color='maroon', label='NGC5257 arm')
plt.errorbar(table['SD']['southarm'], table['tdep']['southarm'], yerr=table['tdep_err']['southarm'], marker='o', color='black', label='NGC5258 southarm')
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
efferr=eff*np.log(10)*data_lit['u_log-eff']
plt.errorbar(molsd,eff, efferr, marker='.', linestyle='none', label='Wilson et al. 2019')
plt.errorbar(Quantities['NGC5258']['SD']['value'].flatten(),epsiff2.flatten(), epsiff2_err.flatten(), linestyle='none', marker='.',label='NGC5258',color=color2)
plt.errorbar(Quantities['NGC5257']['SD']['value'].flatten(),epsiff1.flatten(), epsiff1_err.flatten(), linestyle='none', marker='.',label='NGC5257',color=color1)
ax.tick_params(labelsize=18)
plt.ylim(top=2.0)

plt.errorbar(table['SD']['center'], table['eff']['center'],yerr=table['eff_err']['center'],  marker='o',color='red', label='NGC5257 center')
plt.errorbar(table['SD']['south'], table['eff']['south'], yerr=table['eff_err']['south'], marker='o', color='darkgreen', label='NGC5257 south')
plt.errorbar(table['SD']['west'], table['eff']['west'], yerr=table['eff_err']['west'],  marker='o', color='navy', label='NGC5257 west')
plt.errorbar(table['SD']['arm'], table['eff']['arm'], yerr=table['eff_err']['arm'], marker='o', color='maroon', label='NGC5257 arm')
plt.errorbar(table['SD']['southarm'], table['eff']['southarm'], yerr=table['eff_err']['southarm'], marker='o', color='black', label='NGC5258 southarm')

plt.legend()

# plt.title('ULIRG conversion factor, filling factor=0.1')
fig.tight_layout()
plt.savefig(picDir+'SFE_perfreefall_binned_poster_Gal.png')


############################################################
# counting gas fraction component

### input and bin the fraction image. 
## NGC 5257
fitsfile=mapDir+'NGC5257_gas_vol_fraction_regrid.fits'
fraction_5257=fits_import(fitsfile)[1].data

fitsfile=mapDir+'NGC5257_mstar.fits'
mstar_5257=fits_import(fitsfile)[1].data
vstar_5257=fits_import(fitsfile)[1].data
Phi_5257=1+mstar_5257/Quantities['NGC5257']['SD']['value']*\
                Quantities['NGC5257']['vd']['value']/vstar_5257

xpos=480;ypos=483 # center position
low=[int(xpos-(subnum-1)/2*bin-bin/2), int(ypos-(subnum-1)/2*bin-bin/2)]
high=[int(xpos+(subnum-1)/2*bin+bin/2),int(ypos+(subnum-1)/2*bin+bin/2)]
fraction5257_sub=fraction_5257[low[1]:high[1],low[0]:high[0]]
fraction5257_binned=np.nanmean(np.nanmean(fraction5257_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

Phi5257_sub=Phi_5257[low[1]:high[1],low[0]:high[0]]
Phi5257_binned=np.nanmean(np.nanmean(Phi5257_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

## NGC 5258
fitsfile=mapDir+'NGC5258_gas_vol_fraction_regrid.fits'
fraction_5258=fits_import(fitsfile)[1]

fitsfile=mapDir+'NGC5258_mstar.fits'
mstar_5258=fits_import(fitsfile)[1].data
vstar_5258=fits_import(fitsfile)[1].data
Phi_5258=1+mstar_5258/Quantities['NGC5258']['SD']['value']*\
                Quantities['NGC5258']['vd']['value']/vstar_5258

fitsfile=mapDir+'NGC5258_gas_vol_fraction.fits'
fraction_5258=fits_import(fitsfile)[1].data

xpos=560;ypos=406 # center position
right=15;left=subnum-right-1
lower=15;upper=subnum-lower-1
low=[int(xpos-left*bin-bin/2), int(ypos-lower*bin-bin/2)]
high=[int(xpos+right*bin+bin/2),int(ypos+upper*bin+bin/2)]
fraction5258_sub=fraction_5258[low[1]:high[1],low[0]:high[0]] # xpos and ypos different in spectral cube?
fraction5258_binned=np.nanmean(np.nanmean(fraction5258_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

Phi5258_sub=Phi_5258[low[1]:high[1],low[0]:high[0]]
Phi5258_binned=np.nanmean(np.nanmean(Phi5258_sub.reshape(subnum,10,subnum,10),axis=-1),axis=1)

### Calculate the efficiency per freefall time
epsiff2_star=epsiff2*np.sqrt(fraction5258_binned/0.5)
epsiff2star_err=epsiff2_err*np.sqrt(fraction5258_binned/0.5)
epsiff1_star=epsiff1*np.sqrt(fraction5257_binned/0.5)
epsiff1star_err=epsiff1_err*np.sqrt(fraction5257_binned/0.5)

filename=logDir+'NGC5257_aperture_Mstar.csv'
fraction_aperture1=pd.read_csv(filename, index_col=0)
filename=logDir+'NGC5258_aperture_Mstar.csv'
fraction_aperture2=pd.read_csv(filename, index_col=0)
fraction_aperture=pd.concat([fraction_aperture1, fraction_aperture2])
table['fraction']=fraction_aperture['Mstar']

table['eff']=table['eff']*(table['fraction']/0.5)**0.5
table['eff_err']=table['eff_err']*(table['fraction']/0.5)**0.5

fig=plt.figure()
ax=plt.subplot()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('$\epsilon_{ff}$', fontsize=20)
plt.xlabel('$\Sigma_{mol} (M_{\odot}/pc^2)$', fontsize=20)

eff=10**data_lit['log-eff']
efferr=eff*np.log(10)*data_lit['u_log-eff']
plt.errorbar(molsd,eff, efferr, marker='.', linestyle='none', label='Wilson et al. 2019')
plt.errorbar(Quantities['NGC5258']['SD']['value'].flatten(),epsiff2_star.flatten(), epsiff2star_err.flatten(), linestyle='none', marker='.',label='NGC5258',color=color2)
plt.errorbar(Quantities['NGC5257']['SD']['value'].flatten(),epsiff1_star.flatten(), epsiff1star_err.flatten(), linestyle='none', marker='.',label='NGC5257',color=color1)

plt.errorbar(table['SD']['center'], table['eff']['center'],yerr=table['eff_err']['center'],  marker='o',color='red', label='NGC5257 center')
plt.errorbar(table['SD']['south'], table['eff']['south'], yerr=table['eff_err']['south'], marker='o', color='darkgreen', label='NGC5257 south')
plt.errorbar(table['SD']['west'], table['eff']['west'], yerr=table['eff_err']['west'],  marker='o', color='navy', label='NGC5257 west')
plt.errorbar(table['SD']['arm'], table['eff']['arm'], yerr=table['eff_err']['arm'], marker='o', color='maroon', label='NGC5257 arm')
plt.errorbar(table['SD']['southarm'], table['eff']['southarm'], yerr=table['eff_err']['southarm'], marker='o', color='black', label='NGC5258 southarm')

plt.legend()
ax.tick_params(labelsize=18)
plt.ylim(top=2.0)
# plt.title('ULIRG conversion factor, filling factor=0.1')
fig.tight_layout()
plt.savefig(picDir+'SFE_perfreefall_binned_poster_Gal_star.png')



############################################################
# relationship with Toomre factor

### Qtotsub_binned1, epsiff1; Qtotsub_binned2, epsiff2. correlation. 

Qwhole=np.hstack((Qtotsub_binned1.flatten(), Qtotsub_binned2.flatten()))
epsiff_whole=np.hstack((epsiff1.flatten(), epsiff2.flatten()))
tdep_whole=np.hstack((Quantities['NGC5257']['tdep']['value'].flatten(), Quantities['NGC5258']['tdep']['value'].flatten()))

bad=~np.logical_or(np.isnan(Qwhole), np.isnan(epsiff_whole))
Qwhole_filter=np.compress(bad, Qwhole)
epsiff_whole_filter=np.compress(bad, epsiff_whole)
tdep_whole_filter=np.compress(bad, tdep_whole)

pearsonr(Qwhole_filter, epsiff_whole_filter)

SFE_NGC5257=1/Quantities['NGC5257']['tdep']['value']
SFEerr_NGC5257=SFE_NGC5257*Quantities['NGC5257']['tdep']['error']/Quantities['NGC5257']['tdep']['value']

SFE_NGC5258=1/Quantities['NGC5258']['tdep']['value']
SFEerr_NGC5258=SFE_NGC5258*Quantities['NGC5258']['tdep']['error']/Quantities['NGC5258']['tdep']['value']

# import Leroy et al. 2008 data. 
SFE_ref=10**np.loadtxt(logDir+'SFE_Q_value.csv', delimiter=',')[:,1]
SFEup_ref=10**np.loadtxt(logDir+'SFE_Q_top.csv', delimiter=',')[:,1]-SFE_ref
SFElow_ref=SFE_ref-10**np.loadtxt(logDir+'SFE_Q_bottom.csv', delimiter=',')[:,1]
SFEerr_ref=np.vstack([SFElow_ref, SFEup_ref])

Q_ref=10**np.loadtxt(logDir+'SFE_Q_value.csv', delimiter=',')[:,0]
Qup_ref=10**np.loadtxt(logDir+'SFE_Q_right.csv', delimiter=',')[:,0]-Q_ref
Qlow_ref=Q_ref-10**np.loadtxt(logDir+'SFE_Q_left.csv', delimiter=',')[:,0]
Qerr_ref=np.vstack([Qlow_ref, Qup_ref])

fig=plt.figure()
plt.yscale('log')
plt.xscale('log')
plt.ylabel('$1/t_{dep}\ (year^{-1})$', fontsize=15)
plt.xlabel('Toomre factor $Q_{tot}$', fontsize=15)
plt.errorbar(Qtotsub_binned1.flatten(), 1/Quantities['NGC5257']['tdep']['value'].flatten(), SFEerr_NGC5257.flatten(), marker='.', label='NGC5257', color=color1, linestyle='none')
plt.errorbar(Qtotsub_binned2.flatten(), 1/Quantities['NGC5258']['tdep']['value'].flatten(),  SFEerr_NGC5258.flatten(), marker='.', label='NGC5258', color=color2, linestyle='none')
plt.errorbar(Q_ref, SFE_ref, xerr=Qerr_ref, yerr=SFEerr_ref, linestyle='none', label='Leroy et al. 2008')
plt.legend(fontsize=15)
# cbar=plt.colorbar(sc)
# cbar.set_label(r'$ \beta $', rotation='vertical')
fig.tight_layout()
plt.savefig(picDir+'tdep_Q.png')


fig=plt.figure()
plt.yscale('log')
plt.xscale('log')
plt.ylabel('$\epsilon_{ff}$', fontsize=15)
plt.xlabel('Toomre factor $Q_{tot}$', fontsize=15)
plt.errorbar(Qtotsub_binned1.flatten(), epsiff1.flatten(),yerr=epsiff1_err.flatten(), marker='.', label='NGC5257', color=color1, linestyle='none')
plt.errorbar(Qtotsub_binned2.flatten(), epsiff2.flatten(),yerr=epsiff2_err.flatten(), marker='.', label='NGC5258', color=color2, linestyle='none')
plt.legend(fontsize=15)
fig.tight_layout()
plt.savefig(picDir+'perfreefall_Q.png')

# # considering the shear beta. 
# fig=plt.figure()
# plt.yscale('log')
# plt.ylabel('$\epsilon_{ff}$')
# plt.xlabel('Toomre factor Q')
# sc=plt.scatter(Qtotsub_binned1, epsiff1,c=betasub_binned1,  marker='.', label='NGC5257', cmap='brg_r')
# plt.scatter(Qtotsub_binned2, epsiff2,c=betasub_binned2,  marker='.', label='NGC5258', cmap='brg_r')
# # plt.legend()
# cbar=plt.colorbar(sc)
# cbar.set_label(r'$ \beta $', rotation='vertical')
# plt.savefig(picDir+'perfreefall_Q_beta.png')

