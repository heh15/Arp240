'''
May 3rd
'''

Dir='/1/home/heh15/workingspace/Arp240/ratio/'
image_12CO10='NGC5258_12CO10_combine_uvrange_smooth_regrid21_masked'
image_12CO21='NGC5258_12CO21_combine_uvtaper_smooth_masked'
image_ratio='NGC5258_2110_ratio_uvtaper'
origDir=os.getcwd()
imageDir=Dir+'NGC5258/2110/uvtaper/'
measureDir=Dir+'NGC5258/scatter_plot/'
os.chdir(measureDir)

# beam size: 2.186 arcsec, 1.896 arcsec.


factor=[22,19]
imhead(imagename=imageDir+image_12CO10+'.image.mom0')
imrebin(imagename=imageDir+image_ratio+'.image',factor=factor,outfile=image_ratio+'_rebin.image')
imrebin(imagename=imageDir+image_12CO21+'.image.mom0',factor=factor,outfile=image_12CO21+'_rebin.image')
imrebin(imagename=imageDir+image_12CO10+'.image.mom0',factor=factor,outfile=image_12CO10+'_rebin.image')

os.chdir(origDir)

# convert it to the fits file. 
imagelist=glob.glob('*.image*')
for image in imagelist:
    fitsimage=re.sub('\.image.*$','.fits',image)
    exportfits(imagename=image,fitsimage=fitsimage)

