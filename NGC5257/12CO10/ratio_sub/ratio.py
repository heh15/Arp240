############################################################
# section 1

# task: (12m-7m) subtraction map
immath(imagename=['NGC_5257_CO10_7m.image.mom0',
                  'NGC_5257_CO10_12m.image.mom0'],
       expr='IM1-(IM0/45)',
       outfile='NGC_5257_127_sub.image.mom0')
# The subtraction map is similar to the 12m map.

# task: (12m/7m) ratio map
immath(imagename=['NGC_5257_CO10_12m.image.mom0',
                  'NGC_5257_CO10_7m.image.mom0'],
       expr='IM0*45/IM1',
       outfile='NGC_5257_127_ratio.image.flux')
# Does beam size matter? No. 
# The beam size given and measured is different? According to the given. 
imstat('NGC_5257_712_ratio.image',region='ellipse [[204.97001930deg, 0.84030271deg], [17.3516arcsec, 15.5251arcsec], 90.00000000deg] coord=ICRS, corr=[I], linewidth=1, linestyle=-, symsize=1, symthick=1, color=magenta, font="DejaVu Sans", fontsize=11, fontstyle=normal, usetex=false')
imsubimage('NGC_5257_712_ratio.image',mask='NGC_5257_712_ratio.image>0',outfile='NGC_5257_712_ratio.image.pb')
imsubimage('NGC_5257_127_ratio.image',mask='NGC_5257_127_ratio.image>0',outfile='NGC_5257_127_ratio.image.pb')
# The 712 ratio map has a region similar to the 7m source region,with median value to be 1.32.(indeed overlaped,99%)
# The 127 ratio map has no such region, median value of the same region is 0.76


# task:12m/(12m-7m) 
immath(imagename=['NGC_5257_127_sub.image.mom0',
                  'NGC_5257_CO10_12m.image.mom0'],
       expr='IM1/IM0',
       outfile='NGC_5257_12712_ratio.image.mom0')
# The ratio map would show the galaxy structure with 95%.
imsubimage('NGC_5257_12712_ratio.image.mom0',mask='NGC_5257_12712_ratio.image.mom0>0',outfile='NGC_5257_12712_ratio.image.mom0.cor')
# The masked map shows the structrue similar to the source in 12m map
imstat('NGC_5257_12712_ratio.image.mom0.cor',region='centralregion.crtf')
# The median value of centra region is 1.67
imstat('NGC_5257_12712_ratio.image.mom0.cor',region='region.crtf')
# The median value is 0.67(12/127)
# 12/(combine-7m) ?

#task:12m/combine and 7m/combine
immath(imagename=['NGC_5257_CO10_12m.image.mom0',
                  'NGC_5257_CO10_combine.image.mom0'],
       expr='IM0/IM1',
       outfile='NGC_5257_12c_ratio.image.mom0')
imsubimage('NGC_5257_12c_ratio.image.mom0',mask='NGC_5257_12c_ratio.image.mom0>0',outfile='NGC_5257_12c_ratio.image.mom0.cor')
imstat('NGC_5257_12c_ratio.image.mom0.cor',region='centralregion.crtf')['median'][0]
# median value of central region: 0.91
imstat('NGC_5257_12c_ratio.image.mom0.cor',region='region.crtf')['median'][0]
# median value of the whole region: 0.81
# there is no obvious structure, a blury silouette of the whole. 

immath(imagename=['NGC_5257_CO10_7m.image.mom0',
                  'NGC_5257_CO10_combine.image.mom0'],
       expr='IM0/(IM1*45)',
       outfile='NGC_5257_7c_ratio.image.mom0')
imsubimage('NGC_5257_7c_ratio.image.mom0',mask='NGC_5257_7c_ratio.image.mom0>0',outfile='NGC_5257_7c_ratio.image.mom0.cor')
imstat('NGC_5257_7c_ratio.image.mom0.cor',region='centralregion.crtf')['median'][0]
# median value of central region: 0.43
imstat('NGC_5257_7c_ratio.image.mom0.cor',region='region.crtf')['median'][0]
# median value of the whole region: 1.13
# a clearer structure of the whole galaxy.


############################################################
# section 2

#task: combine-7m map.
immath(imagename=['NGC_5257_CO10_combine.image.mom0',
                  'NGC_5257_CO10_7m.image.mom0'],
       expr='IM0-(IM1/45)',
       outfile='NGC_5257_c7_sub.image.mom0')
# similar to 12m map. 

immath(imagename=['NGC_5257_CO10_combine.image.mom0',
                  'NGC_5257_CO10_12m.image.mom0'],
       expr='IM0-IM1',
       outfile='NGC_5257_c12_sub.image.mom0')
# a elliptical region remsembles the 7m map
