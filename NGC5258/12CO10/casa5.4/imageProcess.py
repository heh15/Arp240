impbcor(imagename='NGC5258_12CO10_combine_contsub.image/',
        pbimage='NGC5258_12CO10_combine_contsub.pb/',
        outfile='NGC5258_12CO10_combine_contsub_pbcor.image')
exportfits(imagename='NGC5258_12CO10_combine_contsub_pbcor.image',fitsimage='NGC5258_12CO10_combine_contsub_pbcor.fits')

imcollapse(imagename='NGC5258_12CO10_combine_contsub.pb/',
           axes=3,
           function="mean",
           outfile='NGC5258_12CO10_combine_contsub_mom0.pb')

impbcor(imagename='NGC5258_12CO10_combine_contsub.mom0/',
        pbimage='NGC5258_12CO10_combine_contsub_mom0.pb',
        outfile='NGC5258_12CO10_combine_contsub_pbcor.mom0')

exportfits(imagename= 'NGC5258_12CO10_combine_contsub_pbcor.mom0',fitsimage='NGC5258_12CO10_combine_contsub_pbcor_mom0.fits')

# make the moment 8 map of the NGC5258.
immoments(imagename='NGC5258_12CO10_combine_contsub.image',moments=8,outfile='NGC5258_12CO10_combine_contsub.mom8')

exportfits(imagename='NGC5258_12CO10_combine_contsub.mom8',fitsimage='NGC5258_12CO10_combine_contsub_mom8.fits')

# exportfits the cube
exportfits(imagename='NGC5258_12CO10_combine_contsub.image',fitsimage='NGC5258_12CO10_combine_contsub.fits')


## 2rms cut for moment 0 map
immoments(imagename='NGC5258_12CO10_combine_contsub.image', moments=0, chans='10~60', includepix=[2*1.6e-3, 100], outfile='NGC5258_12CO10_combine_contsub_2rms.mom0')

impbcor(imagename='NGC5258_12CO10_combine_contsub_2rms.mom0', 
        pbimage='NGC5258_12CO10_combine_contsub_mom0.pb/', 
        outfile='NGC5258_12CO10_combine_contsub_2rms_pbcor.mom0')

exportfits(imagename='NGC5258_12CO10_combine_contsub_2rms_pbcor.mom0', fitsimage='NGC5258_12CO10_combine_contsub_2rms_pbcor_mom0.fits')
