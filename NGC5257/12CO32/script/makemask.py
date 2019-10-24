Dir='/1/home/heh15/workingspace/Arp240/NGC5257/12CO32/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'

inpimage=Dir+'NGC5257_12CO10_noise42.mom0'
region=regionDir+'hinge_init.crtf'
output=Dir+'NGC5257_12CO10_noise42.mom0:mask0'
makemask(mode='copy',inpimage=inpimage,inpmask=region,output=output)

# makemask(mode='delete',inpmask=image+':mask0')


imsmooth(imagename='NGC5257_12CO10_noise42.mom0',
         major='3.789arcsec',
         minor='2.989arcsec',
         pa='-17.357deg',
         outfile='NGC5257_12CO10_noise42_masked_smooth.mom0')

image=imageDir+'NGC5257co32_all.map40r.mom0'
ia.open(image)
ia.calcmask("'"+image+"'>0")
ia.close()
