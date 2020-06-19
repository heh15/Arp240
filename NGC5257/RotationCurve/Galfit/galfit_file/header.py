'''
Sept.18th, 2018
'''

from astropy.io import fits

filename='NGC5257_12CO21_combine_noise40_nocut_mom0.fits'
hdul=fits.open(filename)
hdr=hdul[0].header

hdr['GAIN']= 1.0
hdr['RDNOISE']=0.217

hdul.writeto(filename.replace('.fits', '_header.fits'))
hdul.close()
