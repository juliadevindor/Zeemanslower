#######################################################################################################################
# Calculation of an ideal magnetic field for a Zeeman slower                                                          #
#######################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as scc

#constants
hPlanck=scc.h #6.62607004e-34 Js
muB=scc.physical_constants['Bohr magneton'][0] # 9.274009994e-24 J/T

def slower_field(pos, L0,v0,wavelength,omega, omega0): #ideal field of a Zeeman slower
  # length of slower L0 in m
    B0 = hPlanck*v0/(wavelength*muB) # y-axis interception
    Bbias = 0 #(hPlanck/(2*np.pi))*(omega-omega0)/(muB) # shifts the field up/downwards
    B = np.empty(len(pos))

    for i in range(0,len(pos)):
        if 1-pos[i]/L0 < 0: # to avoid "nan"
            B[i]=0
        elif pos[i]>L0:
            B[i]=0
        else:
            B[i] = B0 * np.sqrt(1-pos[i]/L0) + Bbias #formula for the ideal Zeeman slower field
    return B

L0=0.5#m length of the slower
v0=800 #m/s #maximum velocity, that gets decelerated
wavelength=670.992e-9 #m (Gehm)
omega=446799923264221.4 #Hz
omega0=446.7896e12 #Hz (Gehm)

pos=np.linspace(-0.0,L0+0.2,num=1000) #position on the z-axis
B = np.empty(len(pos))

#calculate B-field
B=slower_field(pos,L0,v0,wavelength,omega,omega0)

#Write to file
file=open("B(z).txt","w+")
for i in range(len(pos)):
    if B[i]!=0:
        file.write("{};{}\n".format(pos[i]-0.5,B[i]*1e4))
file.close()

#Plotting
fig, ax = plt.subplots()
ax.plot(pos+0.1,B*1e4,label="v0=1000m/s")
plt.annotate("z in m", xy=(1.01, 0), ha='left', va='top', xycoords='axes fraction', fontsize=22)
plt.annotate("B in Gauss", xy=(-0.05, 1.10), ha='left', va='top', xycoords='axes fraction', fontsize=22)
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
plt.rcParams.update({'font.size': 15})
xticks = ax.xaxis.get_major_ticks()
xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.legend()
plt.grid()
plt.show()


