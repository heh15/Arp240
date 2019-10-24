'''
Oct. 18th, 2018
'''
import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re

############################################################
# basic setting

Dir='/home/heh15/workingspace/Arp240/continuum/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'
logDir=Dir+'log/'

ratio=1.5

origDir = os.getcwd()
os.chdir(imageDir)

index=['107GHz','95GHz','33GHz']
header=['flux','error']

continuum=np.empty((3,2))
center=np.empty((3,2))

############################################################
# function

def beam_get(imagename,region_init):
    beam=imfit(imagename=imagename,region=region_init)
    x_value=beam['results']['component0']['shape']\
             ['direction']['m0']['value']
    y_value=beam['results']['component0']['shape']\
             ['direction']['m1']['value']
    bmaj_value=beam['results']['component0']\
                ['shape']['majoraxis']['value']
    bmin_value=beam['results']['component0']['shape']\
                ['minoraxis']['value']
    pa_value=beam['results']['component0']['shape']\
              ['positionangle']['value']
    x=str(x_value)+'rad'
    y=str(y_value)+'rad'
    bmaj=str(bmaj_value*ratio)+'arcsec'
    bmin=str(bmin_value*ratio)+'arcsec'
    pa=str(pa_value)+'deg'
    region='ellipse[['+x+','+y+'],['+bmaj+','+bmin+'],'+pa+']'

    return region

def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)

array_round=np.vectorize(round_sig)

def table_output(results,header,index):
    results_tmp=array_round(results,sig=3)
    output=np.column_stack((index,results_tmp))
    header_tmp=copy.copy(header)
    header_tmp.insert(0,' ')
    output=np.vstack((header_tmp,output))

    return output

##################################################
# main program

# south source

# 107GHz flux and uncertainty
imagename='NGC5257_13CO10_cont_spw01.image.image/'
region_init=regionDir+'south_source.crtf'
region_rms=regionDir+'emission_free.crtf'
region=beam_get(imagename,region_init)
flux=imstat(imagename=imagename,region=region)\
                    ['flux']['0']
peak=imstat(imagename=imagename)['max']['0']
rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
error_rel=rms/peak
error=flux*error_rel
continuum[0][0]=flux
continuum[0][1]=error

# 95GHz flux and uncertainty
imagename='NGC5257_13CO10_cont_spw23.image.image/'
region_init=regionDir+'south_source.crtf'
region_rms=regionDir+'emission_free.crtf'
region=beam_get(imagename,region_init)
flux=imstat(imagename=imagename,region=region)\
                    ['flux']['0']
peak=imstat(imagename=imagename)['max']['0']
rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
error_rel=rms/peak
error=flux*error_rel
continuum[1][0]=flux
continuum[1][1]=error

# 33 GHz flux and uncertainty
imagename='NGC5257_33GHz_smooth_regrid.image/'
region_init=regionDir+'south_source.crtf'
region_rms=regionDir+'emission_free.crtf'
region=beam_get(imagename,region_init)
flux=imstat(imagename=imagename,region=region)\
                    ['flux']['0']
peak=imstat(imagename=imagename)['max']['0']
rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
error_rel=rms/peak
error=flux*error_rel
continuum[2][0]=flux
continuum[2][1]=error

# central region

# 107GHz flux and uncertainty
imagename='NGC5257_13CO10_cont_spw01.image.image/'
region_init=regionDir+'center.crtf'
region_rms=regionDir+'emission_free.crtf'
region=beam_get(imagename,region_init)
flux=imstat(imagename=imagename,region=region)\
                    ['flux']['0']
peak=imstat(imagename=imagename)['max']['0']
rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
error_rel=rms/peak
error=flux*error_rel
center[0][0]=flux
center[0][1]=error

# 95 GHz flux and uncertainty
imagename='NGC5257_13CO10_cont_spw23.image.image/'
region_init=regionDir+'center.crtf'
region_rms=regionDir+'emission_free.crtf'
region=beam_get(imagename,region_init)
flux=imstat(imagename=imagename,region=region)\
                    ['flux']['0']
peak=imstat(imagename=imagename)['max']['0']
rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
error_rel=rms/peak
error=flux*error_rel
center[1][0]=flux
center[1][1]=error

# 33 GHz flux and uncertainty
imagename='NGC5257_33GHz_smooth_regrid.image/'
region_init=regionDir+'center.crtf'
region_rms=regionDir+'emission_free.crtf'
region=beam_get(imagename,region_init)
flux=imstat(imagename=imagename,region=region)\
                    ['flux']['0']
peak=imstat(imagename=imagename)['max']['0']
rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
error_rel=rms/peak
error=flux*error_rel
center[2][0]=flux
center[2][1]=error


os.chdir(origDir)

# save the results.
output=table_output(continuum,header,index)
fmt='%5s %10s %10s'

filename=logDir+'south_flux.txt'
with open(filename,'w') as file:
    file.write('south region flux \n')
    np.savetxt(file,output,fmt)
    
output_center=table_output(center,header,index)
fmt='%5s %10s %10s'
filename=logDir+'south_flux.txt'
with open(filename,'a') as file:
    file.write('\n')
    file.write('center flux \n')
    np.savetxt(file,output_center,fmt)

