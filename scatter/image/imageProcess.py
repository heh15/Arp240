############################################################
# NGC 5257
os.chdir('NGC5257')

# smooth the image to 1.02,0.55,-64.5. A little bit small. 
imsmooth(imagename='NGC5257_12CO21_combine_noise40.image/',
         major='1.1arcsec',
         minor='0.6arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5257_12CO21_combine_sinbeam.image')

exportfits(imagename='NGC5257_12CO21_combine_sinbeam.image',
           fitsimage='NGC5257_12CO21_combine_sinbeam_cube.fits')

imsmooth(imagename='NGC5257_12CO21_combine_noise40.image/',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5257_12CO21_combine_smooth.image')

exportfits(imagename='NGC5257_12CO21_combine_smooth.image',
           fitsimage='NGC5257_12CO21_combine_smooth.fits')

# export the pb corrected moment 0 map into fits file
exportfits(imagename = 'NGC5257_12CO21_combine_noise40.pbcor.mom0', fitsimage = 'NGC5257_12CO21_combine_noise40_pbcor_mom0.fits')

# pbcor the cube. 
impbcor(imagename='NGC5257_12CO21_combine_smooth.image/', pbimage='NGC5257_12CO21_combine_noise40.pb/', outfile='NGC5257_12CO21_combine_smooth.image.pbcor')
exportfits(imagename = 'NGC5257_12CO21_combine_smooth.image.pbcor', fitsimage='NGC5257_12CO21_combine_smooth_pbcor.fits')


### 33GHz
os.chdir('33GHz')

# regrid the 33GHz image 12CO 2-1 image. 
imregrid(imagename='ngc5257_Ka_c_r0.5_ms.pbcor.fits',template='NGC5257_12CO21_combine_noise40.image.mom0/',output='NGC5257_33GHz_pbcor_regrid.image')

exportfits(imagename='NGC5257_33GHz_pbcor_regrid.image',fitsimage='NGC5257_33GHz_pbcor_regrid.fits')

imsmooth(imagename='NGC5257_33GHz_pbcor_regrid.image',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5257_33GHz_pbcor_regrid_smooth.image')

exportfits(imagename='NGC5257_33GHz_pbcor_regrid_smooth.image',fitsimage='NGC5257_33GHz_pbcor_regrid_smooth.fits')

# get the pbimage

imregrid(imagename='ngc5257_Ka_c_r0.5_ms.fits',template='NGC5257_12CO21_combine_noise40.image/',output='NGC5257_33GHz_regrid.image')

immath(imagename=['NGC5257_33GHz_regrid.image','NGC5257_33GHz_pbcor_regrid.image'],expr='IM0/IM1',outfile='NGC5257_33GHz_regrid.pb')

exportfits(imagename='NGC5257_33GHz_regrid.pb',fitsimage='NGC5257_33GHz_regrid_pb.fits')

############################################################
# NGC 5258
os.chdir('NGC5258')

## Bin the image

imrebin(imagename='NGC5258_12CO21_combine_noise45.image/',factor=[5,5],outfile='NGC5258_12CO21_combine_5bin.image')

# rms=0.0023
imstat(imagename='NGC5258_12CO21_combine_5bin.image/',region='rmsMeasure.crtf')['rms']

immoments(imagename='NGC5258_12CO21_combine_5bin.image/',moments=0,chans='10~60',includepix=[0.0024*2,100],outfile='NGC5258_12CO21_combine_5bin.image.mom0')
exportfits(imagename='NGC5258_12CO21_combine_5bin.image.mom0/',fitsimage='NGC5258_12CO21_combine_5bin_mom0.fits')

immoments(imagename='NGC5258_12CO21_combine_5bin.image/',moments=2,chans='10~60',includepix=[0.0024*4,100],outfile='NGC5258_12CO21_combine_5bin.image.mom2')
exportfits(imagename='NGC5258_12CO21_combine_5bin.image.mom2/',fitsimage='NGC5258_12CO21_combine_5bin_mom2.fits')

## smooth the image to 1.02,0.55,-64.5. A little bit small. 
imsmooth(imagename='NGC5258_12CO21_combine_noise45.image/',
         major='1.1arcsec',
         minor='0.6arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5258_12CO21_combine_sinbeam.image')

exportfits(imagename='NGC5258_12CO21_combine_sinbeam.image',fitsimage='NGC5258_12CO21_combine_sinbeam.fits')

imsmooth(imagename='NGC5258_12CO21_combine_noise45.image/',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5258_12CO21_combine_noise45_smooth.image')

exportfits(imagename='NGC5258_12CO21_combine_noise45_smooth.image',fitsimage='NGC5258_12CO21_combine_noise45_smooth.fits')

# export pb corrected mom0 map into fits files
exportfits(imagename='NGC5258_12CO21_combine.pbcor.mom0/', fitsimage='NGC5258_12CO21_combine_pbcor_mom0.fits')

# pbcor the cube. 
impbcor(imagename='NGC5258_12CO21_combine_noise45_smooth.image/', pbimage='NGC5258_12CO21_combine_noise45.pb/', outfile='NGC5258_12CO21_combine_smooth.image.pbcor')
exportfits(imagename='NGC5258_12CO21_combine_smooth.image.pbcor', fitsimage='NGC5258_12CO21_combine_smooth_pbcor.fits')

### 33GHz image

## regrid the 33 GHz image to the 12CO 2-1 image

imregrid(imagename='ngc5258_Ka_c_r0.5_ms.pbcor.fits',template='NGC5258_12CO21_combine_noise45.image/',output='NGC5258_33GHz_pbcor_regrid.image')

exportfits(imagename='NGC5258_33GHz_pbcor_regrid.image',fitsimage='NGC5258_33GHz_pbcor_regrid.fits')

## smooth the 33 GHz image to 1.1 arcsec x 0.7 arcesec

imsmooth(imagename='NGC5258_33GHz_pbcor_regrid.image',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5258_33GHz_pbcor_regrid_smooth.image')

exportfits(imagename='NGC5258_33GHz_pbcor_regrid_smooth.image',fitsimage='NGC5258_33GHz_pbcor_regrid_smooth.fits')
