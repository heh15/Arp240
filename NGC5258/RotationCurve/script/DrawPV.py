'''
Jan 17th, 2019

center position: 13h39m57.692s, 00d49m50.838s
'''

############################################################
# directory
Dir='/home/heh15/workingspace/Arp240/NGC5258/RotationCurve/'
imageDir=Dir+'Image/'
scriptDir=Dir+'script/'
testDir=scriptDir+'PV-diagram'

############################################################
# main program

imagename=imageDir+'NGC5258_12CO10_combine_contsub.image/'
outfile=imageDir+'NGC5258_12CO10_pv.image'
mask=imageDir+'NGC5258_12CO10_combine_contsub.mask'
mode='length'
center=["13h39m57.692s","00d49m50.838s"]
length = {"value": 31.5, "unit": "arcsec"}
pa={"value": 213.3, "unit": "deg"}

impv(imagename=imagename,outfile=outfile,mode=mode,center=center,length=length,pa=pa,mask="'"+mask+"'")
exportfits(imagename=outfile,fitsimage=imageDir+'NGC5258_12CO10_pv.fits',velocity=True)

# draw the 12CO21 pv image. 
imagename=imageDir+'NGC5258_12CO21_combine_noise45.image'
outfile=imageDir+'NGC5258_12CO21_pv.image'
mask=imageDir+'NGC5258_12CO21_combine_noise45.mask'
mode='length'
center=["13h39m57.692s","00d49m50.838s"]
length = {"value": 31.5, "unit": "arcsec"}
pa={"value": 213.3, "unit": "deg"}

impv(imagename=imagename,outfile=outfile,mode=mode,center=center,length=length,pa=pa,mask="'"+mask+"'")
exportfits(imagename=outfile,fitsimage=imageDir+'NGC5258_12CO21_pv.fits',velocity=True)
