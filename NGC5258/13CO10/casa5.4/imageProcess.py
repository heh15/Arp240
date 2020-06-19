# make the moment0 map of 13CO1-0
immoments(imagename='NGC5258_CO13_10_12m_mthresh_split_contsub_noise45.image/',moments=0,includepix=[2*0.00064,100],chans='5~30',outfile='NGC5258_13CO10_contsub.mom0')

# primary beam correct 13CO10 data
impbcor(imagename='NGC5258_CO13_10_12m_mthresh_split_contsub_noise45.image/',
        pbimage='NGC5258_CO13_10_12m_mthresh_split_contsub_noise45.pb/',
        outfile='NGC5258_13CO10_contsub_pbcor.image')

imcollapse(imagename='NGC5258_CO13_10_12m_mthresh_split_contsub_noise45.pb/',
           axes=3,
           function="mean",
           outfile='NGC5258_13CO10_contsub_mom0.pb')

impbcor(imagename='NGC5258_13CO10_contsub.mom0/',
        pbimage='NGC5258_13CO10_contsub_mom0.pb',
        outfile='NGC5258_13CO10_contsub_pbcor.mom0')
exportfits(imagename='NGC5258_13CO10_contsub_pbcor.mom0',fitsimage='NGC5258_13CO10_contsub_pbcor_mom0.fits')


# make the moment 8 map

immoments(imagename='NGC5258_CO13_10_12m_mthresh_split_contsub_noise45.image',moments=8,outfile='NGC5258_13CO10_12m_contsub.mom8')

exportfits(imagename='NGC5258_13CO10_12m_contsub.mom8',fitsimage='NGC5258_13CO10_12m_contsub_mom8.fits')

# export fits of the cube
exportfits(imagename='NGC5258_CO13_10_12m_mthresh_split_contsub_noise45.image',fitsimage='NGC5258_13CO10_12m_contsub.fits')
