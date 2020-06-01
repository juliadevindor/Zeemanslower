##########################################################################################################################################
# returns the magnetic field B(z) of an experimental setup which contains 8 coils followed by a MOT consisting of an anti Helmholtz coil #
# author: Julia Winter                                                                                                                   #
##########################################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import scipy.constants as scc
from scipy.optimize import curve_fit

num = 500

hPlanck=scc.h #6.62607004e-34 Js
muB=scc.physical_constants['Bohr magneton'][0] # 9.274009994e-24 J/T
# constants
mu_B=9.274*1e-24 #J/T
hbar=1.055*1e-34 #Js
k_var=2*np.pi/(671*1e-9) #m
Gamma=2*np.pi*5.87*1e6 #Hz
m_Li=6.015*1.66*1e-27 #kg
mu_0=4*np.pi*1e-7 # magnetic field constant
R= 0.043 # inner radius of Zeeman-coils in m (not approved)
d_wire=0.001# thickness of the wire in m
# distances
dist_oven_slower = 0.08
dist_coils_small = 0.002
dist_coils_large = 0.004

def B_coil(pos,N_coil): # magnetic field of a single coil
    coils = 8
    I_coil = 4.8
    L = 0.05
    R_coil = R

    z0 = np.empty([coils])  # center of the coils
    N_wires = np.empty([coils])  # number of wires per layer for each coil
    M = np.empty([coils])  # number of wire layers for each coil
    B=np.empty([num])

    for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
        if j == 0:
            z0[j] = L / 2 + dist_oven_slower  # z0 for the first Zeeman coil
        if j == 1 or j == coils - 1:
            z0[j] = z0[j - 1] + L/ 2 + L / 2 + dist_coils_large  # z0 for second and last Zeeman coil
        if j != 0 and j != 1 and j != coils - 1:
            z0[j] = z0[j - 1] + L / 2 + L/ 2 + dist_coils_small  # z0 for all the other Zeeman coils
    for p in range(0, coils):  # loop over all Zeeman slower coils to obtain the length of all coils
        N_wires[p] = L/ d_wire  # number of wires in each layer
        M[p] = np.abs(round(N_coil[p] / N_wires[p], 0))  # number of layers
    M = M.astype(int)

    for z in range(num):
        for j in range(0, coils):  # loop over Zeeman coils
            for mi in range(1, M[j] + 1):  # loop over all layers of each coil
                R_coil=R + mi * d_wire
                B[z] += mu_0*N_wires[p]*I_coil/(2*L)*((pos[z]-z0[j]+L/2)/np.sqrt(R_coil**2+(pos[z]-z0[j]+L/2)**2) - (pos[z]-z0[j]-L/2)/np.sqrt(R_coil**2+(pos[z]-z0[j]-L/2)**2))

    return B

fig, ax = plt.subplots()
print("plotting")

with open("magnetic_field_real.txt", 'r') as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x,y,label="Spin-flip field")
    #ax.plot(x+0.5,y,".",label="Real slower field")

with open("B(z)_0_5m.txt", "r") as g:  # plot measured magnetic field
    lines = g.readlines()
    xnew = np.asarray([float(line.split(";")[0]) for line in lines])
    ynew = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(xnew+0.5, ynew, label="Ideal slower field of length 0.5m")

N=[1100,900,800,700,600,500,300,100]
L_slower = 0.51  ##??
pos=np.linspace(0,L_slower,num=num)
#plt.plot(pos,1e4*B_coil(pos,N),".")
popt, pcov = curve_fit(B_coil, xnew, ynew,p0=N)
#plt.plot(xdata, func(xdata, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

plt.xlabel("Position in m", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.grid()
plt.rcParams.update({'font.size': 22})
xticks = ax.xaxis.get_major_ticks()
xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.legend(prop={'size': 15})
#plt.ylim(0,1400)
#plt.xlim(0,0.75)
plt.show()



