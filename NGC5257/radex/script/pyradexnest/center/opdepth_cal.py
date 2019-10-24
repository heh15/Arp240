import math
import numpy as np

############################################################
# basic settings

beammaj=2.186
beammin=1.896
freq=110.2

############################################################
# function

def colDen(opdepth,j):
    h=6.62e-34;k=1.38e-23;mu=0.122*3.16e-25
    gj=2*j+1
    colDen=opdepth*3*h*gj/(8*(math.pi)**3*mu**2*j*(math.exp(1)-1))

    return colDen

# method 2, from Cormier et al. 2018, equ.3)
def colden(tau13,I13):
    N13=3e14/(1-math.exp(-1))*tau_13/(1-math.exp(-tau13))*I13

    return N13

def Jy2K(jansky,beammaj,beammin,freq):
    kelvin=jansky/(0.0109*beammaj*beammin*(freq/115.271)**2)

    return kelvin

ratio=13.2322
abun=50
linewidth=51.1

tau_13=1.0/ratio
tau_12=tau_13*abun

jansky=0.44
kelvin=Jy2K(jansky,beammaj,beammin,freq)

I13=kelvin
N13=colden(tau_13,I13)

N12=N13*abun

N=np.log10(N12)

print(N)
