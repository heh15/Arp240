imcollapse(imagename='NGC5257_CO13_10_12m_mthresh_split_contsub_noise45.pb/',
           axes=3,
           function="mean",
           outfile='NGC5257_13CO10_combine_mom0.pb')

impbcor(imagename='NGC5257_CO13_10_12m_mthresh_split_contsub_noise45.mom0/',
        pbimage='NGC5257_13CO10_combine_mom0.pb',
        outfile='NGC5257_13CO10_combine_pbcor.mom0')
exportfits(imagename='NGC5257_13CO10_combine_pbcor.mom0',fitsimage='NGC5257_13CO10_12m_pbcor_contsub_mom0.fits')

# make the moment 8 map

immoments(imagename='NGC5257_CO13_10_12m_mthresh_split_contsub_noise45.image/',moments=8,outfile='NGC5257_13CO10_12m_contsub.mom8')

exportfits(imagename='NGC5257_13CO10_12m_contsub.mom8',fitsimage='NGC5257_13CO10_12m_contsub_mom8.fits')

# exportfits the cube

exportfits(imagename='NGC5257_CO13_10_12m_mthresh_split_contsub_noise45.image',fitsimage='NGC5257_13CO10_12m_contsub.fits')


# correct for the primary beam for the cube. 
impbcor(imagename='NGC5257_CO13_10_12m_mthresh_split_contsub_noise45.image',pbimage='NGC5257_CO13_10_12m_mthresh_split_contsub_noise45.pb',outfile='NGC5257_13CO10_12m_contsub_pbcor.image')
exportfits(imagename='NGC5257_13CO10_12m_contsub_pbcor.image',fitsimage='NGC5257_13CO10_12m_contsub_pbcor.fits')
