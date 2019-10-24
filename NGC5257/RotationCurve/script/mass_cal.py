import numpy as np 
from scipy.interpolate import interp1d

# function 

def mass_cal(radius_kpc,vel,error):
    vel_2=np.multiply(vel,vel)
    mass=2.33*10**5*np.multiply(radius_kpc,vel_2)
    mass_error=mass*2*error/vel
    return mass, mass_error

def mass_minus(mass_radius):
    mass_ring=np.empty(shape=(0,0))
    mass_ring=np.append(mass_ring,mass_radius[0])
    for i in range(mass_radius.shape[0]-1):
        mass=mass_radius[i+1]-mass_radius[i]
        mass_ring= np.append(mass_ring,mass)
    return mass_ring

def error_add(mass_error):
    error_ring=np.empty(shape=(0,0))
    error_ring=np.append(error_ring,mass_error[0])
    for i in range(mass_error.shape[0]-1):
        error=np.sqrt(mass_error[i+1]**2+mass_error[i]**2)
        error_ring= np.append(error_ring,error)
    return error_ring


# main program

####################
# calculate the mass

# import the velocity curve.
filename='log/vel_12CO10.txt' 
input10=np.transpose(np.loadtxt(filename))
radius10_arcsec=input10[0]
vel10=input10[1]
error10=input10[2]

filename='log/vel_12CO21.txt'
input21=np.transpose(np.loadtxt(filename))
radius21_arcsec=input21[0]
vel21=input21[1]
error21=input21[2]

# calculate the mass
radius10_kpc=radius10_arcsec*0.48
mass_12CO10= mass_cal(radius10_kpc,vel10,error10)

M_12CO10=np.log10(mass_12CO10[0])
M_12CO10=np.c_[M_12CO10,0.434*np.divide(mass_12CO10[1],mass_12CO10[0])]
M_12CO10=np.transpose(M_12CO10)

radius21_kpc=radius21_arcsec*0.48
mass_12CO21= mass_cal(radius21_kpc,vel21,error21)

M_12CO21=np.log10(mass_12CO21[0])
M_12CO21=np.c_[M_12CO21,0.434*np.divide(mass_12CO21[1],mass_12CO21[0])]
M_12CO21=np.transpose(M_12CO21)

# fig=plt.figure()
# plt.errorbar(radius_kpc,M_12CO10[0],M_12CO10[1],linestyle='none',label='12CO10')
# plt.errorbar(radius_kpc,M_12CO21[0],M_12CO21[1],linestyle='none',label='12CO21')
# plt.xlabel('radius (kpc)')
# plt.ylabel('log(mass) (solar mass)')

####################
# compare with other mass

#########
# compare with 12CO10 mass. 

# import the stellar mass data. 
filename='../spitzer/log/stellarmass.txt'
M_star=np.transpose(np.loadtxt(filename))[1]
# star=plt.scatter(radius_kpc,M_star,label='stellar mass')

# import the molecular mass data
filename='../spitzer/log/H2mass.txt'
M_gas=np.transpose(np.loadtxt(filename))[1]

# mol=plt.scatter(radius_kpc,M_gas,color='red',label='molecular mass')
# plt.legend(loc='lower right')
# plt.savefig('picture/NGC5257_mass_all.png')

##  stellar fraction within certain radius. 
m_star=np.power(10, M_star)
m_star=m_star[1:7]
error_star=0.3*np.sqrt(m_star*10**6)

mass_total=mass_12CO10[0][1:]
radius_ksub=radius10_kpc[1:]
error_total=mass_12CO10[1][1:]
fraction=np.divide(m_star,mass_total)
fraction_err=np.sqrt((error_total/mass_total)**2+(error_star/m_star)**2)*fraction

fig=plt.figure()
CO10=plt.errorbar(radius_ksub,fraction,fraction_err,color='red',label='CO10 mass')
plt.xlabel('radius (kpc)')
plt.ylabel('stellar mass fraction')


## surface brightness of the stellar mass and total mass. 
mass_total=mass_12CO10[0][1:7]
mass_error=mass_12CO10[1][1:7]

areas=np.empty(shape=(0,0))
area=math.pi*radius_ksub[0]**2*1000**2
areas=np.append(areas,area)

for i in range(mass_total.shape[0]-1):
    area=math.pi*(radius_ksub[i+1]**2-radius_ksub[i]**2)*1000**2
    areas=np.append(areas,area)

m_star=np.power(10, M_star)
m_star=m_star[1:7]
mass_ring=mass_minus(mass_total)
error_ring=error_add(mass_error)

mstar_ring=mass_minus(m_star)
m_sSB=np.divide(mstar_ring,areas)
mass_SB=np.divide(mass_ring,areas)
error_SB= error_ring/areas

# import the surface density and error. 
filename='../spitzer/log/SB_star.txt'
SB_star=np.transpose(np.loadtxt(filename))
SB_star=SB_star[:,1:7]

# fig=plt.figure()
# plt.errorbar(radius_ksub,mass_SB,error_SB, color='red',marker='o')
# plt.errorbar(radius_ksub,SB_star[0],SB_star[1],color='blue',marker='s' )
# plt.xlabel('radius (kpc)')
# plt.ylabel('surface brightness( solar mass/kpc^2)')
# plt.savefig('picture/SB_total.png')

#########
# compare with 12CO21 mass. 

# import the stellar mass data. 
filename='../spitzer/log/stellarmass.txt'
M_star=np.transpose(np.loadtxt(filename))[1]
# star=plt.scatter(radius_kpc,M_star,label='stellar mass')

# import the molecular mass data
filename='../spitzer/log/H2mass.txt'
M_gas=np.transpose(np.loadtxt(filename))[1]

# mol=plt.scatter(radius_kpc,M_gas,color='red',label='molecular mass')
# plt.legend(loc='lower right')
# plt.savefig('picture/NGC5257_mass_all.png')

## stellar fraction within certain radius. 
m_star=np.power(10, M_star)
m_star=m_star[0:7]
error_star=0.3*np.sqrt(m_star*10**6)

mass_total=mass_12CO21[0][1::3]
radius_ksub=radius21_kpc[1::3]
error_total=mass_12CO21[1][1::3]
fraction=np.divide(m_star,mass_total)
fraction_err=np.sqrt((error_total/mass_total)**2+(error_star/m_star)**2)*fraction

CO21=plt.errorbar(radius_ksub,fraction,fraction_err,color='blue',label='CO21 mass')
plt.xlabel('radius (kpc)')
plt.ylabel('stellar mass fraction')
plt.legend()
plt.savefig('picture/fraction_stellar.png')


## surface brightness of the stellar mass and total mass. 
mass_total=mass_12CO21[0][1::3]
mass_error=mass_12CO21[1][1::3]

areas=np.empty(shape=(0,0))
area=math.pi*radius_ksub[0]**2*1000**2
areas=np.append(areas,area)

for i in range(mass_total.shape[0]-1):
    area=math.pi*(radius_ksub[i+1]**2-radius_ksub[i]**2)*1000**2
    areas=np.append(areas,area)

m_star=np.power(10, M_star)
m_star=m_star[0:7]
mass_ring=mass_minus(mass_total)
error_ring=error_add(mass_error)

mstar_ring=mass_minus(m_star)
m_sSB=np.divide(mstar_ring,areas)
mass21_SB=np.divide(mass_ring,areas)
error21_SB= error_ring/areas

# import the surface density and error. 
filename='../spitzer/log/SB_star.txt'
SB_star=np.transpose(np.loadtxt(filename))
SB_star=SB_star[:,0:7]

fig=plt.figure()
plt.errorbar(radius_ksub,mass21_SB,error21_SB, color='green',marker='o')
plt.errorbar(radius_ksub,SB_star[0],SB_star[1],color='blue',marker='s' )
plt.xlabel('radius (kpc)')
plt.ylabel('surface brightness( solar mass/kpc^2)')
plt.savefig('picture/SB_all.png')
