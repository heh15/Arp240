importfits(fitsimage='ngc5257_Ka_c_r0.5_ms.fits',imagename='NGC5257_33GHz.image')

importfits(fitsimage='ngc5257_Ka_c_r0.5_ms.pbcor.fits',imagename='NGC5257_33GHz_pbcor.image')


immath(imagename=['NGC5257_33GHz.image','NGC5257_33GHz_pbcor.image'],
       expr='IM0/IM1',
       outfile='NGC5257_33GHz.pb')
