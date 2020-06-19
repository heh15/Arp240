Dir='/home/heh15/workingspace/Arp240/NGC5257/ratio/'
scriptDir=Dir+'script/measure_1213/'
workDir=Dir+'mask/'
imageDir=Dir+'image/'
logDir=Dir+'log/'
ratioDir=Dir+'1213/contsub/'
maskDir=Dir+'mask/'

# regrid the mask with template of 12CO 2-1 image. 

imregrid(imagename='center_SFR.mask',template=imageDir+'12CO21/NGC5257_12CO21_combine_noise40.image',output='center_SFR_regrid.mask')
exportfits(imagename='center_SFR_regrid.mask',fitsimage='center_SFR_regrid_mask.fits')

imregrid(imagename='hinge_SFR.mask',template=imageDir+'12CO21/NGC5257_12CO21_combine_noise40.image',output='hinge_SFR_regrid.mask')
exportfits(imagename='hinge_SFR_regrid.mask',fitsimage='hinge_SFR_regrid_mask.fits')

imregrid(imagename='south_SFR.mask',template=imageDir+'12CO21/NGC5257_12CO21_combine_noise40.image',output='south_SFR_regrid.mask')
exportfits(imagename='south_SFR_regrid.mask',fitsimage='south_SFR_regrid_mask.fits')

imregrid(imagename='rest.mask',template=imageDir+'12CO21/NGC5257_12CO21_combine_noise40.image',output='rest_regrid.mask')
exportfits(imagename='rest_regrid.mask',fitsimage='rest_regrid_mask.fits')

imregrid(imagename='whole_init.mask',template=imageDir+'12CO21/NGC5257_12CO21_combine_noise40.image',output='whole_init_regrid.mask')
exportfits(imagename='whole_init_regrid.mask',fitsimage='whole_init_regrid_mask.fits')

# export the mask to fits file
exportfits(imagename=maskDir+'center_SFR_co32.mask',fitsimage=maskDir+'center_SFR_co32_mask.fits')

exportfits(imagename=maskDir+'hinge_SFR_co32.mask',fitsimage=maskDir+'hinge_SFR_co32_mask.fits')

exportfits(imagename=maskDir+'south_SFR_co32.mask',fitsimage=maskDir+'south_SFR_co32_mask.fits')

exportfits(imagename=maskDir+'whole_SFR_co32.mask',fitsimage=maskDir+'whole_SFR_co32_mask.fits')

exportfits(imagename=maskDir+'rest_SFR_co32.mask',fitsimage=maskDir+'rest_SFR_co32_mask.fits')
