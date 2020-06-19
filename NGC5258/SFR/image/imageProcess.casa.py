# smooth the 33GHz image to the beam size of 6 arcsec. 

imsmooth(imagename='NGC5258_33GHz_pbcor.image',
         major='6arcsec',
         minor='6arcsec',
         pa='0deg',
         targetres=True,
         outfile='NGC5258_33GHz_pbcor_smooth.image')
exportfits(imagename='NGC5258_33GHz_pbcor_smooth.image',fitsimage='NGC5258_33GHz_pbcor_smooth.fits')

# smooth the unpbcor 33 GHz to the 12CO 1-0 beam size. 
imsmooth(imagename='NGC5258_33GHz.image',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.314deg',
         targetres=True,
         outfile='NGC5258_33GHz_smooth_co10.image')

# smooth the unpbcor 33 GHz to the 12CO 2-1 beam size.
imsmooth(imagename='NGC5258_33GHz.image',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5258_33GHz_smooth_co21.image')


# smooth the pbcor 33GHz to the 12CO1-0 beam size
imsmooth(imagename='NGC5258_33GHz_pbcor.image',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.314deg',
         targetres=True,
         outfile='NGC5258_33GHz_pbcor_smooth_co10.image')

# smooth the pbcor 33 GHz to the 12CO 2-1 beam size.
imsmooth(imagename='NGC5258_33GHz_pbcor.image',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5258_33GHz_pbcor_smooth_co21.image')

exportfits(imagename='NGC5258_33GHz_pbcor_smooth_co21.image', fitsimage='NGC5258_33GHz_pbcor_smooth_co21.fits')

# smooth the 12CO 1-0 moment 0 map.
imsmooth(imagename='NGC5258_12CO10_combine_contsub.mom0/',
         major='2.186arcsec',
         minor='1.896arcsec',
         pa='-87.314deg',
         targetres=True,
         outfile='NGC5258_12CO10_combine_contsub_smooth.mom0/')

# smooth the 12CO 2-1 moment 0 map 
imsmooth(imagename='NGC5258_12CO21_combine.pbcor.mom0',
         major='1.1arcsec',
         minor='0.8arcsec',
         pa='-64.5deg',
         targetres=True,
         outfile='NGC5258_12CO21_combine_pbcor_smooth_co21.mom0')
