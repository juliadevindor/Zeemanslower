##########################################################################################################################################
# returns the magnetic field B(z) of an experimental setup which contains 8 coils followed by a MOT consisting of an anti Helmholtz coil #
# author: Julia Winter                                                                                                                   #
##########################################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import scipy.constants as scc
from scipy.optimize import curve_fit

num = 5000

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
d_wire=0.001# length of the wire in m
b_wire=0.001# thickness of the wire in m
# distances
dist_coils_small = 0.002
dist_coils_large = 0.004

def B_coil_single(z,I,coil):

    z0 = np.empty([coils])  # center of the coils
    N_wires = np.empty([coils])  # number of wires per layer for each coil
    M = np.empty([coils])  # number of wire layers for each coil
    B = 0
    B_0tot=0

    for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
        if j == 0:
            z0[j] = L / 2   # z0 for the first Zeeman coil
        if j == 1 or j == coils - 1:
            z0[j] = z0[j - 1] + L / 2 + L / 2 + dist_coils_large  # z0 for second and last Zeeman coil
        if j != 0 and j != 1 and j != coils - 1:
            z0[j] = z0[j - 1] + L / 2 + L / 2 + dist_coils_small  # z0 for all the other Zeeman coils
    for p in range(0, coils):  # loop over all Zeeman slower coils to obtain the length of all coils
        N_wires[p] = L / d_wire  # number of wires in each layer
        M[p] = np.abs(round(N_coil / N_wires[p], 0))  # number of layers
    M = M.astype(int)

    for j in range(0, coils):  # loop over Zeeman coils
        for mi in range(1, M[j] + 1):  # loop over all layers of each coil
            R_coil = R + mi * b_wire
            B += mu_0 * N_wires[j] * I / (2 * L) * ((z - z0[j] + L / 2) / np.sqrt(R_coil ** 2 + (z - z0[j] + L / 2) ** 2) - (z - z0[j] - L / 2) / np.sqrt(R_coil ** 2 + (z - z0[j] - L / 2) ** 2))
            if j==coil:
                B_0tot += mu_0 * N_wires[j] * I / (2 * L) * ((z - z0[j] + L / 2) / np.sqrt(R_coil ** 2 + (z - z0[j] + L / 2) ** 2) - (z - z0[j] - L / 2) / np.sqrt(R_coil ** 2 + (z - z0[j] - L / 2) ** 2))


    return 1e4 * B_0tot


def B_coil(z, I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11):#,I12,I13,I14,I15,I16,I17): # magnetic field of a single coil

    I_coil=np.array([I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11])#,I12,I13,I14,I15,I16,I17])

    z0 = np.empty([coils])  # center of the coils
    N_wires = np.empty([coils])  # number of wires per layer for each coil
    M = np.empty([coils])  # number of wire layers for each coil
    B=0

    for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
        if j == 0:
            z0[j] = L / 2   # z0 for the first Zeeman coil
        if j == 1 or j == coils - 1:
            z0[j] = z0[j - 1] + L/ 2 + L / 2 + dist_coils_large  # z0 for second and last Zeeman coil
        if j != 0 and j != 1 and j != coils - 1:
            z0[j] = z0[j - 1] + L / 2 + L/ 2 + dist_coils_small  # z0 for all the other Zeeman coils
        #print(j,z0[j])

    for p in range(0, coils):  # loop over all Zeeman slower coils to obtain the length of all coils
        N_wires[p] = L/ d_wire  # number of wires in each layer
        M[p] = np.abs(round(N_coil / N_wires[p], 0))  # number of layers
    M = M.astype(int)

    for j in range(0, coils):  # loop over Zeeman coils
        for mi in range(1, M[j] + 1):  # loop over all layers of each coil
            R_coil=R + mi * d_wire
            B += mu_0*N_wires[j]*I_coil[j]/(2*L) *((z-z0[j]+L/2)/np.sqrt(R_coil**2+(z-z0[j]+L/2)**2) - (z-z0[j]-L/2)/np.sqrt(R_coil**2+(z-z0[j]-L/2)**2))

    return 1e4*B

fig, ax = plt.subplots()
print("plotting")

#I=np.array([5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0])
coils = 11
L = 0.05
N_coil=50
L_field=0.6
#with open("sim_setup/real_magn_field_0_5m.txt", "r") as g:  # plot measured magnetic field
with open("B(z)_0_6m.txt", "r") as g:  # plot measured magnetic field
    lines = g.readlines()
    xnew = np.asarray([float(line.split(";")[0]) for line in lines])
    ynew = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(xnew+0.5, ynew, label="Ideal slower field of length {}m".format(L_field))
L_slower = xnew[-1]+L_field  ##??
pos=np.linspace(0,L_field,num=num)
#plt.plot(pos,B_coil(pos,1100,900,800,700,600,500,300,100),".",label="old real field")
popt, pcov = curve_fit(B_coil, xnew+0.5, ynew,method="trf",bounds=(0,100))
#popt, pcov = curve_fit(B_coil, xnew, ynew)
#print(pcov)
print(popt)
plt.plot(pos,B_coil(pos,*popt),label="Fit: real field of length {}m".format(L_field))

for i in range(coils):
    plt.plot(pos,B_coil_single(pos,popt[i],i),color="black")

file=open("B(z)_fit.txt","w+")
for zpos in pos:
    file.write("{};{}\n".format(zpos-0.5, B_coil(zpos,*popt)))
file.close()

print(0.05*coils+2*dist_coils_large+(coils-3)*dist_coils_small)

ss_res=0
ss_tot=0
for i in range(len(ynew)):
    # residual sum of squares
    ss_res += np.sum((ynew[i] - B_coil(xnew[i],*popt)) ** 2)
    # total sum of squares
    ss_tot += np.sum((ynew[i] - np.mean(ynew)) ** 2)

# r-squared
r2 = 1 - (ss_res / ss_tot)
print("goodness of fit", r2)

#for i in pos:
#    print(i,B_coil(i,1100,900,800,700,600,500,300,100),B_coil(i,*popt))

plt.xlabel("Position in m", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.grid()
plt.rcParams.update({'font.size': 22})
xticks = ax.xaxis.get_major_ticks()
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.legend(prop={'size': 15})
#plt.ylim(0,1400)
#plt.xlim(0,0.75)

ax.spines['left'].set_position('zero')
# set the y-spine
ax.spines['bottom'].set_position('zero')
ax.yaxis.get_major_ticks()[1].label1.set_visible(False)

plt.show()



