# regrid the 33GHz image to the 12CO 2-1 image. 

imregrid(imagename='ngc5258_Ka_c_r0.5_ms.pbcor.fits',template='NGC5258_12CO21_combine_noise45.mom0',output='NGC5258_33GHz_regrid.image')
exportfits(imagename='NGC5258_33GHz_regrid.image',fitsimage='NGC5258_33GHz_regrid.fits')

