
############################################################
# basic setting
ratio=1

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


############################################################
# main program

imagename='spitzer_24um_regrid_minus.image/'
region_init='south_init.crtf'
region_rms='emission_free.crtf'
region=beam_get(imagename,region_init)
flux=imstat(imagename=imagename,region=region)\
                    ['sum']['0']
peak=imstat(imagename=imagename)['max']['0']
rms=imstat(imagename=imagename,region=region_rms)['rms']['0']
error_rel=rms/peak
error=flux*error_rel

