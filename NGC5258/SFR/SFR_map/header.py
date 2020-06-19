from astropy.io import fits


filename='test_regrid_image.fits'
hdul=fits.open('test_regrid_image.fits')
hdr=hdul[0].header

# check the existence of the keyword 'EPOCH'
'EPOCH' in hdr

hdr['EPOCH']=hdr['EQUINOX']

hdul.writeto(filename.strip('.fits')+'_header.fits')
hdul.close()

filename='NGC5257_12CO21_combine_region_regrid.fits'
hdul=fits.open(filename)
hdr=hdul[0].header
hdr['EPOCH']=2000.0
hdr['CTYPE3']='VELO'
hdul.writeto(filename.replace('.fits','_header.fits'))
hdul.close()

filename='NGC5257_12CO10_test2.fits'
hdul=fits.open(filename)
hdr=hdul[0].header
hdr['RESTFREQ']=1.1273E+11
hdul.writeto(filename.replace('.fits','_header.fits'))
hdul.close()
