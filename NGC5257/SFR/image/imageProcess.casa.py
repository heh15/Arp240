### 33 GHz image

# import the 33 GHz image
importfits(fitsimage='ngc5257_Ka_c_r0.5_ms.fits',imagename='NGC5257_33GHz.image')

# smooth the un pbcor image to 2.186X1.896 arcsec
imsmooth(imagename='NGC5257_33GHz.image',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.314deg',
         targetres=True,
         outfile='NGC5257_33GHz_smooth.image')

# smooth the un pbcor image to 1.1X0.8 arcsec
imsmooth(imagename='NGC5257_33GHz.image',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5257_33GHz_smooth_co21.image')

# smooth the 33 GHz image to 6 arcsec
imsmooth(imagename='NGC5257_33GHz.image',major='6arcsec',minor='6arcsec',pa='0deg',targetres=True, outfile='NGC5257_33GHz_smooth_6arcsec.image')
exportfits(imagename='NGC5257_33GHz_smooth_6arcsec.image',fitsimage='NGC5257_33GHz_smooth_6arcsec.fits')

# smooth the primary beam image to 6 arcsec
importfits(fitsimage='ngc5257_Ka_c_r0.5_ms.pbcor.fits',imagename='NGC5257_33GHz_pbcor.image')

imsmooth(imagename='NGC5257_33GHz_pbcor.image',major='6arcsec',minor='6arcsec',pa='0deg',targetres=True, outfile='NGC5257_33GHz_pbcor_smooth_6arcsec.image')
exportfits(imagename='NGC5257_33GHz_pbcor_smooth_6arcsec.image',fitsimage='NGC5257_33GHz_pbcor_smooth_6arcsec.fits')

# create the pbimage. 
immath(imagename=['NGC5257_33GHz.image','NGC5257_33GHz_pbcor.image'],
       expr='IM0/IM1',
       outfile='NGC5257_33GHz.pb')
exportfits(imagename='NGC5257_33GHz.pb',fitsimage='NGC5257_33GHz_pb.fits')

# smooth the 33 GHz image to 2.186 X 1.896 arcsec. 

imsmooth(imagename='NGC5257_33GHz_pbcor.image',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.314deg',
         targetres=True,
         outfile='NGC5257_33GHz_pbcor_smooth.image')

# smooth the 33 GHz to 1.1X0.8 arcsec. 
imsmooth(imagename='NGC5257_33GHz_pbcor.image',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5257_33GHz_pbcor_smooth_co21.image')
exportfits(imagename='NGC5257_33GHz_pbcor_smooth_co21.image',fitsimage='NGC5257_33GHz_pbcor_smooth_co21.fits')

### spitzer

## imregrid the spitzer image and herschel image to the same frame
imregrid(imagename='spitzer_24um.fits',template='NGC5257_12CO10_combine_pbcor.mom0',output='spitzer_24um_regrid.image')
exportfits(imagename='spitzer_24um_regrid.image',fitsimage='spitzer_24um_regrid.fits')

imregrid(imagename='spitzer_24um.fits',template='herschel_70um.image',output='spitzer_24um_regrid70.image')
exportfits(imagename='spitzer_24um_regrid70.image',fitsimage='spitzer_24um_regrid70.fits')


### herschel
imregrid(imagename='herschel_70um.fits',template='NGC5257_12CO10_combine_pbcor.mom0',output='herschel_70um_regrid.image')
exportfits(imagename='herschel_70um_regrid.image',fitsimage='herschel_70um_regrid.fits')

importfits(fitsimage='herschel_70um.fits', imagename='herschel_70um.image')


### 12CO 2-1 image. 
imsmooth(imagename='NGC5257_12CO21_combine_pbcor.mom0/',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5257_12CO21_combine_pbcor_smooth_co21.mom0')

# output the 12CO 2-1 moment 0 map to fits format
exportfits(imagename='NGC5257_12CO21_combine_pbcor.mom0', fitsimage='NGC5257_12CO21_combine_pbcor_mom0.fits')

## smooth the spectral cube moment 0 map. 
imsmooth(imagename='NGC5257_12CO21_pbcor_cube_mom0.fits',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5257_12CO21_pbcor_cube_smooth_co21.mom0')
exportfits(imagename='NGC5257_12CO21_pbcor_cube_smooth_co21.mom0', fitsimage='NGC5257_12CO21_pbcor_cube_smooth_co21_mom0.fits')


### spitzer 3.6um image
exportfits(imagename='spitzer_regrid.image/', fitsimage='spitzer_regrid_36um.fits')
