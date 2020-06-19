impbcor(imagename='NGC5258_12CO21_combine_noise45.image/',
        pbimage='NGC5258_12CO21_combine_noise45.pb/',
        outfile='NGC5258_12CO21_combine_pbcor.image')
exportfits(imagename='NGC5258_12CO21_combine_pbcor.image', fitsimage='NGC5258_12CO21_combine_pbcor.fits')

imcollapse(imagename='NGC5258_12CO21_combine_noise45.pb/',
           axes=3,
           function="mean",
           outfile='NGC5258_12CO21_combine_mom0.pb')

impbcor(imagename='NGC5258_12CO21_combine_noise45.image.mom0/',
        pbimage='NGC5258_12CO21_combine_mom0.pb',
        outfile='NGC5258_12CO21_combine.pbcor.mom0')
exportfits(imagename='NGC5258_12CO21_combine.pbcor.mom0', fitsimage='NGC5258_12CO21_combine_pbcor_mom0.fits')

# moment 8 map
immoments(imagename='NGC5258_12CO21_combine_noise45.image',moments=8,outfile='NGC5258_12CO21_combine.mom8')

# exportfits the cube
exportfits(imagename='NGC5258_12CO21_combine_noise45.image',fitsimage='NGC5258_12CO21_combine.fits')
