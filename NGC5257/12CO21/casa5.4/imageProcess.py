imcollapse(imagename='NGC5257_12CO21_combine_noise40.pb/',
           axes=3,
           function="mean",
           outfile='NGC5257_12CO21_combine_mom0.pb')

impbcor(imagename='NGC5257_12CO21_combine_noise40.image.mom0/',
        pbimage='NGC5257_12CO21_combine_mom0.pb',
        outfile='NGC5257_12CO21_combine_pbcor.mom0')
exportfits(imagename='NGC5257_12CO21_combine_pbcor.mom0',fitsimage='NGC5257_12CO21_combine_pbcor_mom0.fits')
        
impbcor(imagename='NGC5257_12CO21_combine_noise40.image/',
        pbimage='NGC5257_12CO21_combine_noise40.pb/',
        outfile='NGC5257_12CO21_combine_noise40_pbcor.image/')
exportfits(imagename='NGC5257_12CO21_combine_noise40_pbcor.image/',fitsimage='NGC5257_12CO21_combine_pbcor.fits')

# test the pbcor and original map
immath(imagename=['NGC5257_12CO21_combine_noise40.image.mom0','NGC5257_12CO21_combine_pbcor.mom0'],expr='IM0/IM1', outfile='Test_NGC5257_12CO21.pb')

## smooth the 12CO21 to sma resolution. 

imsmooth(imagename='NGC5257_12CO21_combine_noise40.image.mom0/',
         major='3.516arcsec',
         minor='2.83arcsec',
         pa='-7.73deg',
         targetres=True,
         outfile='NGC5257_12CO21_combine_co32.mom0')

## make the moment 8 map
rmtables('NGC5257_12CO21_combine.mom8')
immoments(imagename='NGC5257_12CO21_combine_noise40.image',moments=8,chans='10~60',outfile='NGC5257_12CO21_combine.mom8',includepix=[5*0.0031,100])

## exportfits the cube
exportfits(imagename='NGC5257_12CO21_combine_noise40.image',fitsimage='NGC5257_12CO21_combine.fits')
