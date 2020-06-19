exportfits(imagename='NGC5257_12CO21_combine_noise40.image.mom0',fitsimage='NGC5257_12CO21_combine_noise40_mom0.fits') 

## 33 GHz image
imregrid(imagename='ngc5257_Ka_c_r0.5_ms.pbcor.fits',template='NGC5257_12CO21_combine_noise40.image.mom0',output='NGC5257_33GHz_pbcor_regrid.image')
exportfits(imagename='NGC5257_33GHz_pbcor_regrid.image',fitsimage='NGC5257_33GHz_pbcor_regrid.fits')

# smooth the pbcor image
imsmooth(imagename='NGC5257_33GHz_pbcor_regrid.image',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5257_33GHz_pbcor_regrid_smooth.image')
exportfits(imagename='NGC5257_33GHz_pbcor_regrid_smooth.image', fitsimage='NGC5257_33GHz_pbcor_regrid_smooth.fits')

# smooth the original image and measure the rms
imsmooth(imagename='ngc5257_Ka_c_r0.5_ms.fits',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5257_33GHz_smooth.image')

# stellar mass map. 

