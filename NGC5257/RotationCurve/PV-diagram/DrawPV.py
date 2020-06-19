
imagename='NGC5257_12CO21_combine_noise40.image'
outfile='../PV-diagram/NGC5257_12CO21_pv.image'
mask='NGC5257_12CO21_combine_noise40.mask'
mode='length'
center=["13h39m52.921s","00d50m24.332s"]
length = {"value": 28, "unit": "arcsec"}
pa={"value": 110, "unit": "deg"}

impv(imagename=imagename,outfile=outfile,mode=mode,center=center,length=length,pa=pa,mask=mask)
exportfits(imagename=outfile,fitsimage='NGC5257_12CO21_pv.fits',velocity=True)
