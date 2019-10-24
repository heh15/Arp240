############################################################
# basic setting

Dir='/home/heh15/workingspace/Arp240/NGC5257/SFR/'
imageDir=Dir+'image/'
regionDir=Dir+'region/'
logDir=Dir+'log/'
scriptDir=Dir+'script'


############################################################
# function

def beam_get(beam,region_init,ratio=4.0/2.355):
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

############################################################
# main program
imagename=imageDir+'NGC5257_33GHz.image'
pbcorimage=imageDir+'NGC5257_33GHz_pbcor.image'

regions=['center','south']
values=['xvalue','yvalue','flux','error']
dicts=dict.fromkeys(regions)
for key in dicts.keys():
    dicts[key]=dict.fromkeys(values)

### south continuum source 

region_init=regionDir+'southern.crtf'
region_rms=regionDir+'emission_free.crtf'
beam=imfit(imagename=imagename,region=region_init)
region=beam_get(beam,region_init)
flux=imstat(imagename=pbcorimage,region=region)\
                    ['flux']['0']
peak=imstat(imagename=imagename)['max']['0']
rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
error_rel=rms/peak
error=flux*error_rel

x_value=beam['results']['component0']['shape']\
             ['direction']['m0']['value']
y_value=beam['results']['component0']['shape']\
             ['direction']['m1']['value']
dicts['south']['flux']=flux
dicts['south']['error']=error
dicts['south']['xvalue']=x_value
dicts['south']['yvalue']=y_value

# ### center 

# region_init=regionDir+'center.crtf'
# region_rms=regionDir+'emission_free.crtf'
# beam=imfit(imagename=imagename,region=region_init)
# region=beam_get(beam,region_init)
# flux=imstat(imagename=pbcorimage,region=region)\
#                     ['flux']['0']
# peak=imstat(imagename=imagename)['max']['0']
# rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
# error_rel=rms/peak
# error=flux*error_rel

# x_value=beam['results']['component0']['shape']\
#              ['direction']['m0']['value']
# y_value=beam['results']['component0']['shape']\
#              ['direction']['m1']['value']

# dicts['center']['flux']=flux
# dicts['center']['error']=error
# dicts['center']['xvalue']=x_value
# dicts['center']['yvalue']=y_value
