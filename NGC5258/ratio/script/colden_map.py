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

import matplotlib.colors as colors
############################################################
# basic setting
bmaj=2.186
bmin=1.896
pa=-87.314
freq=110
fillings=[0.1, 1.0]

subnum=int(320/5)
############################################################
# function
def fits_import(fitsimage):
    hdr = fits.open(fitsimage)[0].header
    wcs = WCS(hdr).celestial
    data=fits.open(fitsimage)[0].data
    data=np.squeeze(data)
    data_masked=np.ma.masked_invalid(data)
    return wcs, data_masked

def cut_2d(data_masked,position,size,wcs):
    cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

def Jy2K(jansky,beammaj,beammin,freq):
    kelvin=jansky/(0.0109*beammaj*beammin*(freq/115.271)**2)

    return kelvin

############################################################
# directory

Dir='/1/home/heh15/workingspace/Arp240/NGC5258/ratio/'
scriptDir=Dir+'script/'
logDir=Dir+'log/'
picDir=Dir+'picture/'
imageDir=Dir+'image/'
ratioDir=imageDir+'ratio/1213/contsub_pbcor/'

ratioimage=ratioDir+'NGC5258_1213_ratio_pbcor.fits'
_13COmom0=ratioDir+'NGC5258_13CO10_12m_uvrange_smooth_masked_pbcor_mom0.fits'
_12COmom0=ratioDir+'NGC5258_12CO10_combine_contsub_uvrange_smooth_masked_pbcor_mom0.fits'
_12COmom8=imageDir+'12CO10/NGC5258_12CO10_combine_contsub_uvrange_smooth_pbcor_mom8.fits'
_12COmom2=imageDir+'12CO10/NGC5258_12CO10_combine_contsub_pbcor_smooth_mom2.fits'

############################################################
# main program

for filling in fillings:

    ### import ratio map
    fitsimage=ratioimage
    wcs=fits_import(fitsimage)[0]
    data_masked=fits_import(fitsimage)[1]
    ratio_binned=np.nanmean(np.nanmean(data_masked.reshape(subnum,5,subnum,5),axis=-1),axis=1)
    tau=1.0/ratio_binned

    ### import intensity map
    fitsimage=_13COmom0
    wcs=fits_import(fitsimage)[0]
    data_masked=fits_import(fitsimage)[1]

    I13=np.copy(data_masked)
    I13_binned=np.nanmean(np.nanmean(I13.reshape(subnum,5,subnum,5),axis=-1),axis=1)
    I13_binned[I13_binned<6.4e-4*20*5]='nan'
    I13K=I13_binned/(0.0109*bmaj*bmin*(107/115.27)**2)

    ## 12CO10 map
    fitsimage=_12COmom0
    wcs=fits_import(fitsimage)[0]
    data_masked=fits_import(fitsimage)[1]

    I12=np.copy(data_masked)
    I12_binned=np.nanmean(np.nanmean(I12.reshape(subnum,5,subnum,5),axis=-1),axis=1)
    # I12_binned[I12_binned<5*1.6e-3*10*np.sqrt(50)]='nan'
    I12K=I12_binned/(0.0109*bmaj*bmin*(112/115.27)**2)

    ## 12CO10 moment 8 map
    fitsimage=_12COmom8
    data_masked=fits_import(fitsimage)[1]

    I12_peak=np.copy(data_masked)
    I12peak_binned=np.nanmean(np.nanmean(I12_peak.reshape(subnum,5,subnum,5),axis=-1),axis=1)
    I12peak_K=I12peak_binned/(0.0109*bmaj*bmin*(112/115.27)**2)

    ## velocity dispersion
    fitsimage=_12COmom2
    data_masked=fits_import(fitsimage)[1]

    linewidth=np.copy(data_masked*2*np.sqrt(2*np.log(2)))
    linewidth_binned=np.nanmean(np.nanmean(linewidth.reshape(subnum,5,subnum,5),axis=-1),axis=1)

    ### exctitation temperature
    T_ex=(I13K/linewidth_binned)/filling/(1-np.exp(-tau))
    # T_ex=20
    T_ex=5.5/(np.log(1+5.5/(I12peak_K/filling+0.82)))

    ### column density
    N_13=3e14/(1-np.exp(-5.29/T_ex))*tau/(1-np.exp(-tau))*I13K
    N_12=N_13*50
    N_H2=N_12*10000
    mass=N_H2*3.32e-24*(3.1e18)**2/2e33*1.4
    alpha=mass/I12K
    # alpha[alpha>2]='nan'
    alpha[alpha<0]='nan'

    alpha_mean=np.sum(alpha*mass)/np.sum(mass)
    print(str(alpha_mean)+'\n')

# fig=plt.figure()
# plt.scatter(I12K, alpha, marker='.')

# South_I12=9.52;South_alpha=1.33
# beammaj=3.8;beammin=3.0;freq=112
# South_I12K=Jy2K(South_I12,beammaj,beammin,freq)

# South_I13=0.54
# South_I13K=Jy2K(South_I13,beammaj,beammin, 108)
# linewidth_tmp=87
# tau_tmp=South_I13K/South_I12K
# Tex_tmp=(South_I13K/linewidth_tmp)/filling/(1-np.exp(-tau_tmp))
# South_alpha_LTE=3e14/(1-np.exp(-5.29/Tex_tmp))*tau_tmp/(1-np.exp(-tau_tmp))*South_I13K*50*10000*3.32e-24*(3.1e18)**2/2e33*1.4/South_I12K

# plt.scatter([South_I12K], [South_alpha], color='green',label='NGC 5258 south arm')
# # plt.scatter([South_I12K], [South_alpha_LTE], facecolor='none', edgecolor='green', label='NGC 5258 south arm LTE')
# plt.legend(fontsize=15)
# # plt.title('abun=50, filling_factor=0.1')
# plt.xlabel('12CO 1-0 intensity (K km s $^{-1}$)', fontsize=20)
# plt.ylabel(r'$\alpha_{CO}$ ($M_{\odot} pc^{-2} K km s^{-1}$)', fontsize=20)
# fig.tight_layout()
# plt.savefig(picDir+'NGC5258_conversion_factor.png')


# fig=plt.figure()
# ax=plt.subplot()

# position=(32,32); size=(30,30)
# alpha_cut=Cutout2D(data=alpha,position=position,size=size).data

# im=ax.imshow(alpha_cut,origin='lower', norm=colors.PowerNorm(gamma=0.5))
# ax.set_xticks([])
# ax.set_yticks([])
# cbar=plt.colorbar(im)
# cbar.set_label(r'$\alpha_{CO}$', fontsize=25)
# cbar.ax.tick_params(labelsize=20)
# plt.savefig(picDir+'NGC5258_col.png', bbox_inches='tight', pad_inches=0.5)

# # fig=plt.figure()
# # plt.imshow(I12)
 

# fig=plt.figure()
# ax=plt.subplot()

# position=(32,32); size=(30,30)
# alpha_cut=Cutout2D(data=T_ex,position=position,size=size).data

# im=ax.imshow(alpha_cut,origin='lower', norm=colors.PowerNorm(gamma=0.5))
# ax.set_xticks([])
# ax.set_yticks([])
# cbar=plt.colorbar(im)
# cbar.ax.tick_params(labelsize=20)
# cbar.set_label(r'$T_{ex} (K)$', fontsize=25)
# fig.tight_layout()
# plt.savefig(picDir+'NGC5258_Tex.png')

# fig=plt.figure()
# plt.imshow(Tk, origin='lower')
