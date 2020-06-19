### 12CO10 ###

# exportfits the cube

exportfits(imagename='NGC5258_12CO21_combine_uvtaper.image',fitsimage='NGC5258_12CO21_combine_uvtaper.fits')

# do the primary beam correction for the cube
impbcor(imagename='NGC5258_12CO10_combine_contsub_uvrange.image',pbimage='NGC5258_12CO10_combine_contsub_uvrange.pb',outfile='NGC5258_12CO10_combine_contsub_uvrange_pbcor.image/')

# moment 0 map pb. 

imcollapse(imagename='NGC5258_12CO10_combine_contsub_uvrange.pb/',
           axes=3,
           function="mean",
           outfile='NGC5258_12CO10_combine_contsub_uvrange_mom0.pb/')

### 13CO10 ###

# create moment 0 pb image. 
imcollapse(imagename='NGC5258_13CO10_12m_uvrange.pb',
           axes=3,
           function="mean",
           outfile='NGC5258_13CO10_12m_uvrange_mom0.pb')

### 12CO21 ###

