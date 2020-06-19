Dir='/1/home/heh15/workingspace/Arp240/NGC5257/RotationCurve/'
imageDir=Dir+'image/'

imregrid(imagename='mass_map.fits',template=imageDir+'NGC5257_12CO21_combine_noise40.image.mom0/',output='mass_map_regrid.image')
exportfits(imagename='mass_map_regrid.image',fitsimage='mass_map_regrid.fits')
