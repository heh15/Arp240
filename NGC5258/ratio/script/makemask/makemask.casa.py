inpimage='NGC5257_12CO10_combine_contsub_region.image'
mode='copy'
inpmask='anomaly.crtf'
output='anomaly.mask'

makemask(inpimage=inpimage,mode=mode,inpmask=inpmask,output=output)

immath(imagename=['NGC5257_12CO10_combine_contsub_region.mask',
                  'anomaly.mask'],
       expr='IM0-IM1',
       outfile='temp.mask')

inpimage='temp.mask'
inpmask='NGC5257_12CO10_combine_contsub_region_sub_tmp.mask > 0.5'
output='mask0'

ia.open('temp.mask')
ia.calcmask('"temp.mask">0.5',name='mask0')
ia.close()

output='NGC5257_12CO10_combine_contsub_region_sub.mask'
makemask(mode=mode, inpimage='temp.mask', inpmask='temp.mask:mask0', output=output)

rmtables('temp.mask')

