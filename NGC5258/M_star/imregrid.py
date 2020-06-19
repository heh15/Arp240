'''Jan 16th, 2019

Directory in the $workingdir+image/
'''
exportfits(image='NGC5258_12CO10_combine_contsub.mom0',fitsimage='NGC5258_12CO10_contsub_mom0.fits')

imregrid(imagename='SPITZER_I1_39933184_0000_2_E11349850_maic.fits',template='NGC5258_12CO10_combine_contsub.mom0',output='spitzer_36um_regrid.image')
exportfits(imagename='spitzer_36um_regrid.image',fitsimage='spitzer_36um_regrid.fits')


imregrid(imagename='SPITZER_I2_39933184_0000_2_E11352627_maic.fits',template='NGC5258_12CO10_combine_contsub.mom0',output='spitzer_45um_regrid.image')
exportfits(imagename='spitzer_45um_regrid.image',fitsimage='spitzer_45um_regrid.fits')
