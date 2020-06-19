from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import make_lupton_rgb
# import img_scale as IS
import aplpy
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from matplotlib.patches import Arrow
# import montage_wrapper as montage
import matplotlib as mpl 
mpl.rcParams['lines.linewidth']=3



############################################################
# Coordinate parameters 

ra = 15*(13*u.deg + 39*u.arcmin + 55*u.arcsec)
dec = 50*u.arcmin + 10*u.arcsec
center = SkyCoord(ra = ra, dec = dec)
size = u.Quantity((1.5, 2.5),u.arcmin)

############################################################
# functions 

def cut_2d(data_masked,position,size,wcs):
    cut=Cutout2D(data=data_masked,position=position,size=size,wcs=wcs)
    data_cut=cut.data
    wcs_cut=cut.wcs

    return wcs_cut, data_cut

############################################################
# main program

# Make grey scale image with aplpy and overlay the contour.

# contour='/home/heh15/workingspace/Arp240/'\
#     '12CO10/NGC5257/finalImage/'\
#     'NGC5257_12CO10_noise42_mom0.fits'
# contour1='/home/heh15/workingspace/Arp240/'\
#     '12CO10/NGC5258/finalImage/'\
#     'NGC5258_12CO10_noise45_mom0.fits'
# levels=[0.3,2.3,4.3,6.3,8.3]
# gc = aplpy.FITSFigure('j9cv49010_drz.fits',north=True)
# gc.show_grayscale(invert=False)
# gc.add_scalebar(0.000562,color='white')
# gc.scalebar.set_corner('top left')
# gc.scalebar.set_label('1kpc')
# gc.show_contour(contour,levels=levels,colors=['white','green','yellow','orange','red'])
# gc.show_contour(contour1,levels=levels,colors=['white','green','yellow','orange','red'])
# gc.axis_labels.hide()
# gc.tick_labels.hide()
# gc.save('test.png')


# Make grey scale image with asptropy and overlay the contour. not changing the direction of the image.
init=0.565
levels = [2*0.565, 4*0.565, 6*0.565, 8*0.565]
# levels=np.linspace(init,8+init,4)
colors=['white','green','yellow','orange','red']

contour_dir='/home/heh15/workingspace/Arp240/'\
    'NGC5257/12CO10/casa5.4/'\
    'NGC5257_12CO10_combine_contsub_pbcor_mom0.fits'
contour1_dir='/home/heh15/workingspace/Arp240/'\
    'NGC5258/12CO10/casa5.4/'\
    'NGC5258_12CO10_combine_contsub_pbcor_mom0.fits'

hdr = fits.getheader('test_uv_cube.fits')
wcs = WCS(hdr)
data=fits.getdata('test_uv_cube.fits')
contour=fits.open(contour_dir)[0]
contour_data = contour.data[0][0]
contour_wcs=WCS(contour.header).celestial
contour1=fits.open(contour1_dir)[0]
contour1_data = contour1.data[0][0]
contour1_wcs=WCS(contour1.header).celestial

origin=[205.006889,0.8391667]
length=0.0041667
font = {'color':  'white',
        'weight': 'normal',
        }
rotation=117

wcs=wcs.celestial
r=data[0]
g=data[1]
b=data[2]
rgb = make_lupton_rgb(r,g, b,Q=1,stretch=1)
rgb_cut = []

for i in range(np.shape(rgb)[2]):
    wcs_cut, rgb_cut_temp = cut_2d(rgb[:,:,i], center, size, wcs)
    rgb_cut.append(rgb_cut_temp)
rgb_cut = np.dstack(tuple(rgb_cut))

fig=plt.figure()
ax=plt.subplot(projection=wcs_cut)
ax.imshow(rgb_cut,vmin=0,vmax=1.5,origin='lower')
ax.grid(color='white',ls='dotted')
ax.contour(contour_data,levels=levels,transform=ax.get_transform(contour_wcs),colors=colors, linewidths=2.0) 
ax.contour(contour1_data,levels=levels,transform=ax.get_transform(contour1_wcs),colors=colors, linewidths=2.0)
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_major_formatter('hh:mm:ss')
length=0.000583;origin[1]=origin[1]+0.006;origin[0]=origin[0]-0.0125
plt.hlines(y=origin[1],xmin=origin[0],xmax=origin[0]+length,color='w',transform=ax.get_transform('fk5'))
# plt.text(origin[0]+length/2,origin[1]+0.0003,'1kpc',ha='center',va='bottom',color='w',transform=ax.get_transform('fk5'))
# arrow=ax.arrow(origin[0],origin[1],length,0,width=0.0005,transform=ax.get_transform('fk5'),facecolor='white',head_width=0.0015,length_includes_head=True,head_length=0.00075)
# arrow1=ax.arrow(origin[0],origin[1],0,length,width=0.0005,facecolor='white',transform=ax.get_transform('fk5'),head_width=0.0015,length_includes_head=True,head_length=0.0007)
# ax.text(origin[0]+length+0.001,origin[1]+0.001,'E',fontsize=12,transform=ax.get_transform('fk5'),fontdict=font,rotation=rotation)
# ax.text(origin[0],origin[1]+length+0.002,'N',fontsize=12,transform=ax.get_transform('fk5'),fontdict=font,rotation=rotation)
plt.show()
fig.tight_layout()
ax.tick_params(direction='in')
plt.savefig('astropy.png')

'''
# Try to rotate the image with Cutout2D, Cutout2D can't rotate the image.The image could be rotated by changing the wcs direction and the data array direction. 

position=SkyCoord(dec=0.8336*u.degree,ra=204.9833*u.degree,frame='fk5')
size=u.Quantity((2.5,2.5),u.arcmin)
cut=Cutout2D(data=data,position=position,size=size,wcs=wcs)
data_cut=cut.data


w=cut.wcs
w.wcs.crpix=[1367.5, 1447.0]     
w.wcs.cd[0]=[-1.2474e-05, 6.1065e-06]
w.wcs.cd[1]=[-6.1065e-06,-1.2474e-05]

# rotate the image by changing the wcs vector and the data array.

#data_cut=np.rot90(data_cut)

fig=plt.figure()
ax=plt.subplot(projection=w)
ax.imshow(data_cut,vmin=0,vmax=1.5,cmap='gray')
ax.grid(color='white',ls='dotted')
lon=ax.coords[0]
lat=ax.coords[1]
plt.show()
'''




