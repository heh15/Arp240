'''
Oct 26th, 2018

smooth the image to 6 arcsec.
regrid the image to herschel image.
'''

# smooth the herschel and galex to 6 arcsec resolution
imsmooth(imagename='herschel_70um_header.fits',
         kernel='gauss',
         major='6.0arcsec',
         minor='6.0arcsec',
         pa='0deg',
         targetres=True,
         outfile='herschel_70um_header_smooth.image')

imsmooth(imagename='galex_FUV_header.fits',
         kernel='gauss',
         major='6.0arcsec',
         minor='6.0arcsec',
         pa='0deg',
         targetres=True,
         outfile='galex_FUV_header_smooth.fits')

# regrid the spitzer and galex to herschel image. 
imregrid(imagename='galex_FUV_header_smooth.image/',template='herschel_70um_header_smooth.image/',output='galex_FUV_header_smooth_regrid.image/')

imregrid(imagename='spitzer_24um.fits',template='herschel_70um_header_smooth.image/',output='spitzer_24um_regrid.image')

# export the image to fits file

exportfits(imagename='galex_FUV_header_smooth_regrid.image/',fitsimage='galex_FUV_header_smooth_regrid.fits')

exportfits(imagename='herschel_70um_header_smooth.image/',fitsimage='herschel_70um_header_smooth.fits')

exportfits(imagename='spitzer_24um_regrid.image/',fitsimage='spitzer_24um_regrid.fits')
