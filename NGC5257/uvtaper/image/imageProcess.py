import os

os.chdir('12CO10')

imcollapse(imagename='NGC5257_12CO10_combine_contsub_uvrange.pb/',
           axes=3,
           function="mean",
           outfile='NGC5257_12CO10_combine_uvrange_mom0.pb')

## primary beam correction 
impbcor(imagename='NGC5257_12CO10_combine_contsub_uvrange.image',pbimage='NGC5257_12CO10_combine_contsub_uvrange.pb/',outfile='NGC5257_12CO10_combine_contsub_uvrange_pbcor.image')
exportfits(imagename='NGC5257_12CO10_combine_contsub_uvrange_pbcor.image',fitsimage='NGC5257_12CO10_combine_contsub_uvrange_pbcor.fits')
        
os.chdir('../')

# 12CO 21 processing

os.chdir('12CO21')
imcollapse(imagename='NGC5257_12CO21_combine_contsub_uvtaper.pb/',
           axes=3,
           function="mean",
           outfile='NGC5257_12CO21_combine_uvtaper_mom0.pb')

exportfits(imagename='NGC5257_12CO21_combine_contsub_uvtaper.image',fitsimage='NGC5257_12CO21_combine_contsub_uvtaper.fits')

## primary beam correction. 

impbcor(imagename='NGC5257_12CO21_combine_contsub_uvtaper.image',pbimage='NGC5257_12CO21_combine_contsub_uvtaper.pb/',outfile='NGC5257_12CO21_combine_contsub_uvtaper_pbcor.image')
exportfits(imagename='NGC5257_12CO21_combine_contsub_uvtaper_pbcor.image',fitsimage='NGC5257_12CO21_combine_contsub_uvtaper_pbcor.fits')


## 13CO 1-0 process. 
os.chdir('13CO10')

imcollapse(imagename='NGC5257_13CO10_12m_uvrange.pb',
           axes=3,
           function="mean",
           outfile='NGC5257_13CO10_12m_uvrange_mom0.pb')
