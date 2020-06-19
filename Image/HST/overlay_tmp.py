from astropy.coordinates import Angle, SkyOffsetFrame, ICRS
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import make_lupton_rgb
import img_scale as IS
import aplpy
# from reproject import reproject_interp


# from skimage.measure import block_reduce
# hdulist=fits.open('j9cv49010_drc.fits');hdr= hdulist[1];hdulist.close()
'''
hdr = fits.getheader('color_hst_09124_62_wfpc2_f814w_f300w_wf_sci.fits')
wcs = WCS(hdr)
data=fits.getdata('color_hst_09124_62_wfpc2_f814w_f300w_wf_sci.fits')
# data_extract=np.percentile(data,[0,100])

image_r=data[0]
image_g=data[1]
image_b=data[2]
rgb_default = make_lupton_rgb(image_r, image_g, image_b,Q=10,stretch=0.5)

data = np.swapaxes(data, 0, 2)
data = np.swapaxes(data, 0, 1)
w=wcs.celestial
# print data.shape
fig=plt.figure()
ax=fig.add_subplot(111,projection=w)
plt.imshow(rgb_default,vmin=0,vmax=10)
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_major_formatter('hh:mm:ss')
plt.show()
'''

# ax.contour(hdu.data, levels=np.logspace(-4.7, -3., 10), colors='white', alpha=0.5)
# image_hist = plt.hist(final_image.flatten(), 1000)

'''
contour='/home/heh15/workingspace/Arp240/'\
    '12CO10/NGC5257/finalImage/'\
    'NGC5257_12CO10_noise42_mom0.fits'
contour1='/home/heh15/workingspace/Arp240/'\
    '12CO10/NGC5258/finalImage/'\
    'NGC5258_12CO10_noise45_mom0.fits'
levels=[0.3,2.3,4.3,6.3,8.3]
gc = aplpy.FITSFigure('j9cv49010_drz.fits',north=True)
gc.show_grayscale(invert=False)
gc.add_scalebar(0.000562,color='white')
gc.scalebar.set_corner('top left')
gc.scalebar.set_label('1kpc')
gc.show_contour(contour,levels=levels,colors=['white','green','yellow','orange','red'])
gc.show_contour(contour1,levels=levels,colors=['white','green','yellow','orange','red'])
gc.axis_labels.hide()
gc.tick_labels.hide()
gc.save('test.png')
'''
# files1=['icom15030_drz.fits','j9cv49020_drz.fits','j9cv49010_drz.fits']
# aplpy.make_rgb_cube(files1,'hst_combine.fits')
# aplpy.make_rgb_image('test2_cube.fits','test_cube.png')

'''
files1=['icom15030_drz.fits','j9cv49020_drz.fits','j9cv49010_drz.fits']
files=['j9cv49020_drz.fits','j9cv49010_drz.fits','u6dw6201r_drz.fits']

def isolate_image_extension(fits_file, extension):
    header = fits.getheader(fits_file, extension)
    data = fits.getdata(fits_file, extension)
    fits.writeto('%s_image.fits' % fits_file.rstrip('.fits'), data, header)

for i in range(len(files)):
    files[i]=files[i].replace('.fits','_image.fits')

aplpy.make_rgb_cube(files,output='test_uv'+'_cube.fits',north=True)
aplpy.make_rgb_image('test_uv'+'_cube.fits','test_uv_cube.png')


isolate_image_extension('u6dw6201r_drz.fits',1)
'''


'''
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
hdu1 = fits.open(get_pkg_data_filename('j9cv49020_drz.fits'))[1]
hdu2 = fits.open(get_pkg_data_filename('j9cv49010_drz.fits'))[1]

array, footprint = reproject_interp(hdu2, hdu1.header)

wcs = WCS(hdu1.header)
r=hdu1.data
b=hdu2.data
avg_g=0.5*(r+b)
rgb = make_lupton_rgb(r,avg_g, b,Q=0,stretch=0)

fig=plt.figure()
ax=fig.add_subplot(111)
plt.imshow(rgb,vmin=0,vmax=10)
plt.show()
'''

# hdr = fits.getheader('test2_cube.fits')
# wcs = WCS(hdr)
# data=fits.getdata('test2_cube.fits')
# # data_extract=np.percentile(data,[0,100])


# max_value=[50,6,1]

# for i in range(3):
#     tmp = max_value[i]
#     data[i, :, :] *= 1.0/tmp


# image_r=data[0]
# image_b=data[2]
# image_g=data[1]
# rgb_default = make_lupton_rgb(image_r, image_g, image_b,stretch=1,Q=10)

# w=wcs.celestial

# fig=plt.figure()
# ax=fig.add_subplot(111,projection=w)
# plt.imshow(rgb_default,vmin=0,vmax=1)
# lon = ax.coords[0]
# lat = ax.coords[1]
# lon.set_major_formatter('hh:mm:ss')
# plt.show()

# aplpy.make_rgb_image('test_uv'+'_cube.fits','test_uv_cube.png',vmin_r=0,vmin_g=0, vmin_b=0,stretch_r='arcsinh',stretch_g='arcsinh',stretch_b='arcsinh',vmid_g=0.1,vmid_b=0.003)


files1=['icom15030_drz.fits','j9cv49020_drz.fits','j9cv49010_drz.fits']
files=['j9cv49020_drz.fits','j9cv49010_drz.fits','u6dw6201r_drz.fits']

def isolate_image_extension(fits_file, extension):
    header = fits.getheader(fits_file, extension)
    data = fits.getdata(fits_file, extension)
    fits.writeto('%s_image.fits' % fits_file.rstrip('.fits'), data, header)

for i in range(len(files)):
    files[i]=files[i].replace('.fits','_image.fits')
'''
aplpy.make_rgb_cube(files,output='test_uv'+'_cube1.fits')
aplpy.make_rgb_image('test_uv'+'_cube.fits','test_uv_cube.png',vmin_r=0,vmin_g=0, vmin_b=0,stretch_r='arcsinh',stretch_g='arcsinh',stretch_b='arcsinh')
'''

'''
isolate_image_extension('u6dw6201r_drz.fits',1)
'''

aplpy.make_rgb_image('test_uv_cube.fits','test_uv_cube.png',vmin_r=0,vmin_g=0, vmin_b=0,stretch_r='arcsinh',stretch_g='arcsinh',stretch_b='arcsinh')
