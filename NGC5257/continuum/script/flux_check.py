'''
Oct. 19th, 2018

check the flux in different spectral windows of 12CO 1-0. 

'''
import numpy as np
import os,glob
import scipy.ndimage as sni
import sys
import re
import copy

############################################################
# basic setting

Dir='/home/heh15/workingspace/Arp240/continuum/'
image12='/home/heh15/workingspace/Arp240/12CO10/continuum/'
image13='/home/heh15/workingspace/Arp240/13CO10/continuum/'
regionDir=Dir+'region/'
logDir=Dir+'log/'

ratio=1.0

origDir = os.getcwd()


index=['spw0','spw1','spw2','spw3']
header=['flux','error']

center=np.empty((4,2))
center13=np.empty((4,2))

# float_formatter = lambda x: "%.2f" % x
# np.set_printoptions(formatter={'float': lambda x: format(x, '3.2E')})

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

# 12CO 1-0
# central region
spws=np.array(range(4))

for spw in spws:
    imagename=image12+'NGC5257_12CO10_cont_12m_spw'+str(spw)+'.image'
    region_init=regionDir+'center.crtf'
    region_rms=regionDir+'emission_free.crtf'
    region=beam_get(imagename,region_init)
    flux=imstat(imagename=imagename,region=region)\
                        ['flux']['0']
    peak=imstat(imagename=imagename)['max']['0']
    rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
    error_rel=rms/peak
    error=flux*error_rel
    center[spw][0]=flux
    center[spw][1]=error


# 13CO 1-0
# central region
spws=np.array(range(4))

for spw in spws:
    imagename=image13+'NGC5257_cont_12m_13CO_spw'+str(spw)+'.image'
    region_init=regionDir+'center.crtf'
    region_rms=regionDir+'emission_free.crtf'
    region=beam_get(imagename,region_init)
    flux=imstat(imagename=imagename,region=region)\
                        ['flux']['0']
    peak=imstat(imagename=imagename)['max']['0']
    rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
    error_rel=rms/peak
    error=flux*error_rel
    center13[spw][0]=flux
    center13[spw][1]=error

filename=Dir+'log/frequency.txt'
frequency=np.loadtxt(filename,skiprows=1)


# save the results.
output=table_output(center,header,index)
fmt='%5s %10s %10s'

filename=logDir+'center_flux.txt'
with open(filename,'w') as file:
    file.write('12CO 1-0 flux \n')
    np.savetxt(file,output,fmt)


output13=table_output(center13,header,index)
fmt='%5s %10s %10s'
filename=logDir+'center_flux.txt'
with open(filename,'a') as file:
    file.write('\n')
    file.write('13CO 1-0 flux \n')
    np.savetxt(file,output13,fmt)
