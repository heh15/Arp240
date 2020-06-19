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

##################################################
# main program
Dir='/1/home/heh15/workingspace/Arp240/scatter/'
picDir=Dir+'picture/'
imageDir=Dir+'image/'
logDir=Dir+'log/'


############################################################
# basic setting
alpha=4.3
ratio=0.77
G=6.67e-11
pc=3.1e16
incl=0.45
beammaj=1.004
beammin=0.556
freq=230.54
deltaV=10

Galaxy=['NGC5257','NGC5258']
sigmas=['analytic','python','dispersion']
rms_disper=dict.fromkeys(Galaxy)
for key in Galaxy:
    rms_disper[key]=dict.fromkeys(sigmas)

############################################################
# function

def mom2_uncertainty(imcube, rms, chans=np.shape(imcube[:,0,0])[0],corr=None):
    mom0=imcube.moment(order=0)
    mom2=imcube.linewidth_sigma()
    cubenp=np.array(imcube)
    cubenp[cubenp<4*rms]='nan'
    chan_width=(imcube.spectral_axis[1]-imcube.spectral_axis[0]).value
    chan_width=round(chan_width,1)
    sigmaI=rms*deltaV*np.sqrt(chans)
    linewidth=chan_width*chans
    sigma_dis1=sigmaI/np.array(mom0)*linewidth**2/(4*math.sqrt(10))/(2*np.array(mom2))
    err1=sigmaI/mom0_array
    err2=sigma_dis1/mom2_array
    # calculate the correlation coefficient
    if corr==None:
        corr=np.ndarray(shape=(np.shape(cubenp)[1],np.shape(cubenp)[2]))
        for i in range(np.shape(cubenp)[1]):
            for j in range(np.shape(cubenp)[2]):
                dispers=(vel- mom1_array[i][j])**2*cubenp[:,i,j]
                df=pd.DataFrame(np.transpose(np.vstack((dispers, cubenp[:,i,j]))))
                corr[i][j]=df.corr()[0][1]
    
    sigma_dis2=np.array(mom2)*np.sqrt(err1**2+err2**2-2*corr*err1*err2)

    return sigma_dis1, sigma_dis2, corr

############################################################
#### NGC 5257

# name=imageDir+'NGC5257/NGC5257_12CO21_combine_noise40_image.fits'
# imagecube=SpectralCube.read(name)
# Imcube_tmp=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)
# Imcube_tmp2=Imcube_tmp.to(u.K)

# ### common beam selection
# common_beam = Imcube_tmp.beams.common_beam()
# Imcube = Imcube_tmp.convolve_to(common_beam)

# rmscube=cube.calc_noise_in_cube(Imcube)
# outcube=cube.find_signal_in_cube(Imcube,rmscube,snr_hi=5)
# outcube.write(imageDir+'NGC5257_12CO21_cube.fits')

outcube=SpectralCube.read(imageDir+'NGC5257_12CO21_cube.fits')

mom0=outcube.moment(order=0)
# mom0.quicklook()
mom0_array=np.array(mom0)
threshold=mom0_array<5.0
mom0_array[threshold]='nan'

mom1=outcube.moment(order=1)
mom1_array=np.array(mom1)

mom2=outcube.linewidth_sigma()
# mom2.quicklook()
mom2_array=np.array(mom2)
mom2_array[threshold]='nan'
filter2=mom2_array<5.0
mom2_array[filter2]='nan'

#####################################################################
### error calculation

N=70
sigmaI=3e-3*deltaV*np.sqrt(N)
rms2=sigmaI/np.array(mom0)*700**2/(4*math.sqrt(10))/(2*np.array(mom2))

cubenp=np.array(outcube)

vel=np.linspace(-300,390,70)
corr=np.ndarray(shape=(np.shape(cubenp)[1],np.shape(cubenp)[2]))


# for i in range(np.shape(cubenp)[1]):
#     for j in range(np.shape(cubenp)[2]):
#         dispers=(vel- mom1_array[i][j])**2*cubenp[:,i,j]
#         df=pd.DataFrame(np.transpose(np.vstack((dispers, cubenp[:,i,j]))))
#         corr[i][j]=df.corr()[0][1]

# outfits=logDir+'NGC5257_correlation.fits'
# hdu=fits.PrimaryHDU(corr)
# hdu.writeto(outfits,overwrite=True)

fitsimage=logDir+'NGC5257_correlation.fits'
corr=fits.open(fitsimage)[0].data

err1=sigmaI/mom0_array
err2=rms2/mom2_array
rms_tot=mom2_array*np.sqrt(err1**2+err2**2-2*corr*err1*err2)
filter3=np.isnan(rms_tot)
rms2[filter3]='nan'

# figure=plt.figure()
# plt.title('uncertainty of velocity dispersion')
# plt.xlabel('error(nominator only) km/s')
# plt.ylabel('error(correlation considered) km/s')
# # plt.scatter(rms2.flatten(),rms2.flatten())
# # plt.fill_between(rms2.flatten(),upper.flatten(),lower.flatten(),color='gray')
# plt.plot(rms2.flatten(),rms2.flatten())
# sc=plt.scatter(rms2.flatten(),rms_tot.flatten(),c=mom2_array.flatten(),marker='.',cmap='viridis_r')
# cbar=plt.colorbar(sc)
# cbar.set_label('velocity dispersion (km/s)')
# # plt.plot(rms2.flatten(),upper.flatten(),color='blue',label='50 km/s uncertainty')
# # plt.plot(rms2.flatten(),lower.flatten(),color='blue')
# plt.savefig(picDir+'error.png')

rms_disper['NGC5257']['analytic']= rms2
rms_disper['NGC5257']['python']=rms_tot
rms_disper['NGC5257']['dispersion']=mom2_array



rms=3e-3
# results=mom2_uncertainty(outcube, rms)
sigma_dis1=results[0]
sigma_dis2=results[1]
corr=results[2]

rms_disper['NGC5257']['analytic']= sigma_dis1
rms_disper['NGC5257']['python']=sigma_dis2
rms_disper['NGC5257']['dispersion']=mom2_array

# # fig=plt.figure()
# # plt.errorbar(mom0_array.flatten(),mom2_array.flatten(),rmsf.flatten(),linestyle='none',marker='o')

# # rmscube.median(axis=2)
# # fig=plt.figure()
# # err1=sigmaI/mom0_array
# # err2=rms2/mom2_array
# # ratio=err1/err2
# # plt.plot(mom2.flatten(),ratio.flatten(),linestyle='none',marker='o')

# ############################################################
# #### ng 5258

# # name=imageDir+'NGC5258/NGC5258_12CO21_combine_noise45.fits'
# # imagecube=SpectralCube.read(name)
# # Imcube_tmp=imagecube.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=225.46*10**9*u.Hz)
# # Imcube_tmp2=Imcube_tmp.to(u.K)

# # ### common beam selection
# # common_beam = Imcube_tmp.beams.common_beam(tolerance=1e-5)
# # Imcube = Imcube_tmp2.convolve_to(common_beam)

# # rmscube=cube.calc_noise_in_cube(Imcube)
# # outcube=cube.find_signal_in_cube(Imcube,rmscube,snr_hi=5)
# # outcube.write(imageDir+'NGC5258/NGC5258_12CO21_cube.fits')

# outcube=SpectralCube.read(imageDir+'NGC5258/NGC5258_12CO21_cube.fits')

# mom0=outcube.moment(order=0)
# # mom0.quicklook()
# mom0_array=np.array(mom0)
# threshold=mom0_array<5.0
# mom0_array[threshold]='nan'

# mom1=outcube.moment(order=1)
# mom1_array=np.array(mom1)

# mom2=outcube.linewidth_sigma()
# # mom2.quicklook()
# mom2_array=np.array(mom2)
# mom2_array[threshold]='nan'
# filter2=mom2_array<5.0
# mom2_array[filter2]='nan'

# #####################################################################
# ### error calculation

# N=70
# sigmaI=3e-3*deltaV*np.sqrt(N)
# rms2=sigmaI/np.array(mom0)*700**2/(4*math.sqrt(10))/(2*np.array(mom2))


# cubenp=np.array(outcube)

# vel=np.linspace(-300,390,70)
# corr=np.ndarray(shape=(np.shape(cubenp)[1],np.shape(cubenp)[2]))


# for i in range(np.shape(cubenp)[1]):
#     for j in range(np.shape(cubenp)[2]):
#         dispers=(vel- mom1_array[i][j])**2*cubenp[:,i,j]
#         df=pd.DataFrame(np.transpose(np.vstack((dispers, cubenp[:,i,j]))))
#         corr[i][j]=df.corr()[0][1]

# outfits=logDir+'NGC5258_correlation.fits'
# hdu=fits.PrimaryHDU(corr)
# hdu.writeto(outfits,overwrite=True)

# err1=sigmaI/mom0_array
# err2=rms2/mom2_array
# rms_tot=mom2_array*np.sqrt(err1**2+err2**2-2*corr*err1*err2)
# filter3=np.isnan(rms_tot)
# rms2[filter3]='nan'

# rms_disper['NGC5258']['analytic']= rms2
# rms_disper['NGC5258']['python']=rms_tot
# rms_disper['NGC5258']['dispersion']=mom2_array


############################################################
#### plot the figure

figure=plt.figure()
plt.title('uncertainty of velocity dispersion')
plt.xlabel('error(nominator only) km/s')
plt.ylabel('error(correlation considered) km/s')
# plt.scatter(rms2.flatten(),rms2.flatten())
# plt.fill_between(rms2.flatten(),upper.flatten(),lower.flatten(),color='gray')
# plt.plot(rms_disper['NGC5258']['analytic'].flatten(),rms_disper['NGC5258']['analytic'].flatten(),color='gray')
# sc=plt.scatter(rms_disper['NGC5258']['analytic'].flatten(),rms_disper['NGC5258']['python'].flatten(),c=rms_disper['NGC5258']['dispersion'].flatten(),marker='.',cmap='viridis_r')
plt.plot(rms_disper['NGC5257']['analytic'].flatten(),rms_disper['NGC5257']['analytic'].flatten())
sc=plt.scatter(rms_disper['NGC5257']['analytic'].flatten(),rms_disper['NGC5257']['python'].flatten(),c=rms_disper['NGC5257']['dispersion'].flatten(),marker='.',cmap='viridis_r')
cbar=plt.colorbar(sc)
cbar.set_label('velocity dispersion (km/s)')
# plt.plot(rms2.flatten(),upper.flatten(),color='blue',label='50 km/s uncertainty')
# plt.plot(rms2.flatten(),lower.flatten(),color='blue')
plt.savefig(picDir+'error.png')




