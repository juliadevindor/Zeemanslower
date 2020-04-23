import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

mu_0=4*np.pi*1e-7 # magnetic field constant

def B_coil(I_coil, N_coil, L_coil, R_coil, z0_coil, pos): # magnetic field of a single coil
    return 1e4*mu_0*N_coil*I_coil/(2*L_coil)*((pos-z0_coil+L_coil/2)/np.sqrt(R_coil**2+(pos-z0_coil+L_coil/2)**2) - (pos-z0_coil-L_coil/2)/np.sqrt(R_coil**2+(pos-z0_coil-L_coil/2)**2)) # magnetic field of a single coil


I_1=10
N_1=150
L_1=0.025
R_1=0.043
z0_1=-0.09-L_1/2

with open("Measurements/MOT_new.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split()[0]) for line in lines])
    y = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot((x-33.5)/100, y*10, label="MOT w/o screws")
    plt.vlines(0.05,-350,750)
plt.vlines(0.05-0.117,-350,750)
plt.vlines(-0.1,-350,750)


with open("Measurements/Zeemanslower_old_standalone_Bz.txt", "r") as g:
    lines = g.readlines()
    x4 = np.asarray([float(line.split()[0]) for line in lines])
    y4 = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot((x4-32.5)/100, y4 ,label="Messung Zeemanslower von uns")
xourslower=np.arange(0.0,0.62,0.01)
funcourslower=interpolate.interp1d((x4-32.5)/100,y4, bounds_error=False, fill_value=0)
#plt.plot(xourslower-0.5, funcourslower(xourslower)*10)
funcourmot=interpolate.interp1d((x-33.5)/100,y*10, bounds_error=False, fill_value=0)
#plt.plot(xourslower-0.5, funcourmot(xourslower-0.5))
file= open("./TEST.txt","w+")
for i in range(len(xourslower)):
     file.write("{};{}\n".format(xourslower[i]-0.55, funcourslower(xourslower)[i]*10+funcourmot(xourslower-0.5)[i]+B_coil(I_1,N_1,L_1,R_1,z0_1,xourslower-0.5)[i]))
file.close()

with open("./TEST.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x,y,".",label="Messung gesamt von uns")
plt.plot(x,B_coil(I_1,N_1,L_1,R_1,z0_1,xourslower-0.5))


with open("Measurements/example_magnetic_field.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x, y, label="good field")


plt.grid()
plt.xlabel("distance from right edge of slower/ m", fontsize=15)
plt.ylabel("B/ mT", fontsize=15)
plt.legend(loc=1, prop={'size': 13})
plt.show()
plt.close()

