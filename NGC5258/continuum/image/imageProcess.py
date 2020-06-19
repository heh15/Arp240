importfits(fitsimage='ngc5258_Ka_c_r0.5_ms.fits', imagename='NGC5258_33GHz.image')

importfits(fitsimage='ngc5258_Ka_c_r0.5_ms.pbcor.fits', imagename='NGC5258_33GHz_pbcor.image')

# create the primary beam correction image
immath(imagename=['NGC5258_33GHz.image', 'NGC5258_33GHz_pbcor.image'],
       expr='IM0/IM1',
       outfile='NGC5258_33GHz.pb')
