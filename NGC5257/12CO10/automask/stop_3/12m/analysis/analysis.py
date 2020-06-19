immoments(imagename='NGC_5257_CO10_12m_auto.image/',moments=[0],outfile='NGC_5257_CO10_12m_auto.image.mom0')


immath(imagname=['NGC_5257_CO10_12m_auto.image.mom0/',
                  'NGC_5257_CO10_combine_Nat.image.mom0/'],
       expr='IM1-IM0',
       outfile='NGC5257_c12_sub.mom0')
