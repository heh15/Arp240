'''
Oct. 23rd, 2018

mask the pixel with only one channel selected

Dir=/home/heh15/workingspace/Arp240/scatter/test_2chan/
'''

from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D

############################################################
# basic settings

ra=204.97066052
dec=50.406569999999995
position=SkyCoord(dec=dec*u.arcmin,ra=ra*u.degree,frame='icrs')
size=u.Quantity((54,42),u.arcsec)

############################################################
# function

def cut_3d(data,position,size,wcs):
    for i in range(data.shape[0]):
        cut=Cutout2D(data=data[i],position=position,size=size,wcs=wcs)
        if i==0:
            data_cut=cut.data
        elif i==1:
            data_cut=np.stack((data_cut,cut.data))
        else:
            temp=np.expand_dims(cut.data,axis=0)
            data_cut=np.concatenate((data_cut,temp))
    wcs_cut=cut.wcs
    return data_cut, wcs_cut


############################################################
# main program

fitsimage='NGC5257_12CO21_combine_mask_4sig.fits'
hdr = fits.open(fitsimage)[0].header
wcs = WCS(hdr).celestial
mask_temp=fits.open(fitsimage)[0].data[0]


mask_temp[isnan(mask_temp)]=0
mask_new=np.empty((mask_temp.shape[0],mask_temp.shape[1],mask_temp.shape[2]))
for i in range(mask_new.shape[0]):
    if i==0:
        mask_new[i,:,:]=mask_new[i,:,:]+mask_new[i+1,:,:]
    elif i==(mask_new.shape[0]-1):
        mask_new[i,:,:]=mask_temp[i,:,:]+mask_temp[i-1,:,:]
    else:
        mask_new[i,:,:]=mask_temp[i,:,:]+mask_temp[i+1,:,:]+mask_temp[i-1,:,:]

# use this mask to make the moment 2 map. 
