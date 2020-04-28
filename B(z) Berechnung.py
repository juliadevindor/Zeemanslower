import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import os
import glob
import scipy.constants as scc


hPlanck=scc.h #6.62607004e-34 Js
muB=scc.physical_constants['Bohr magneton'][0] # 9.274009994e-24 J/T

def slower_field(pos, L0,v0,wavelength,omega, omega0):
  # length of slower L0 in m
    B0 = hPlanck*v0/(wavelength*muB) #0.08 #0.19185 # y-Achsenabschnitt
    #print(B0)
    Bbias = -5e-2 #0#(hPlanck/(2*np.pi))*(omega-omega0)/(muB)#-1.478*1e-9  # shiftet nur den Plot nach oben/ unten
    B = np.empty(len(pos))

    for i in range(0,len(pos)):
        print(pos[i])
        #if pos[i]<0.0:
        #    B[i]=0
        if 1-pos[i]/L0 < 0: # to avoid "nan"
            B[i]=0
        elif pos[i]>L0:
            B[i]=0
        else:
            B[i] = B0 * np.sqrt(1-pos[i]/L0) + Bbias
    return B

L0=0.828#m
v0_0=1200 #m/s
v0_1=800 #m/s
v0_2=500 #m/s
wavelength=670.992e-9 #m (Gehm)
omega=446799923264221.4 #Hz (mein Wert) #446799900000000 #Hz (Andis Wert)
omega0=446.7896e12 #Hz (Gehm)

pos=np.linspace(-0.0,L0+0.2,num=1000)
B = np.empty(len(pos))
B=slower_field(pos,L0,v0_0,wavelength,omega,omega0)
file=open("B(z).txt","w+")
for i in range(len(pos)):
    file.write("{};{}\n".format(pos[i]-0.5,B[i]*1e4))
file.close()

fig, ax = plt.subplots()
with open("B(z)_1m.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.6,y,label="andi")
with open("B(z).txt","r") as f: # plot measured magnetic field
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.6,y, label="PLOT", color="black")
ax.plot(pos+0.1,B*1e4,label="v0=1000m/s")
plt.annotate("z in m", xy=(1.01, 0), ha='left', va='top', xycoords='axes fraction', fontsize=15)
plt.annotate("B in Gauss", xy=(-0.05, 1.05), ha='left', va='top', xycoords='axes fraction', fontsize=15)
#plt.hlines(0,pos[0]+0.1,pos[-1]+0.1)
ax.spines['left'].set_position('zero')
# set the y-spine
ax.spines['bottom'].set_position('zero')
#plt.xlabel("Position / m")
plt.rcParams.update({'font.size': 22})
#plt.ylabel("B / T")
plt.legend()
plt.grid()
plt.show()

