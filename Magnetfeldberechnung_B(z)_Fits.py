#################################################################################################################
# returns the magnetic field B(z) of an experimental setup which contains 8 coils followed by a MOT consisting  #
#################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def B_coil_single(z,I,coil): # magnetic field of a single slower coil
    z0 = np.empty([coils])  # center of the coils
    N_wires = np.empty([coils])  # number of wires per layer for each coil
    M = np.empty([coils])  # number of wire layers for each coil
    B = 0
    B_0tot=0

    for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
        if j == 0:
            z0[j] = L[j] / 2  # z0 for the first Zeeman coil
        if j == 1 or j == coils - 1:
            z0[j] = z0[j - 1] + L[j] / 2 + L[j] / 2 + dist_coils_large  # z0 for second and last Zeeman coil
        if j != 0 and j != 1 and j != coils - 1:
            z0[j] = z0[j - 1] + L[j] / 2 + L[j] / 2 + dist_coils_small  # z0 for all the other Zeeman coils
    for p in range(0, coils):  # loop over all Zeeman slower coils to obtain the length of all coils
        N_wires[p] = L[p] / d_wire  # number of wires in each layer
        M[p] = np.abs(round(N_coil[p] / N_wires[p], 0))  # number of layers
    M = M.astype(int)
    for j in range(0, coils):  # loop over Zeeman coils
        for mi in range(1, M[j] + 1):  # loop over all layers of each coil
            R_coil = R + mi * b_wire
            B += mu_0 * N_wires[j] * I / (2 * L[j]) * ((z - z0[j] + L[j] / 2) / np.sqrt(R_coil ** 2 + (z - z0[j] + L[j] / 2) ** 2) - (z - z0[j] - L[j] / 2) / np.sqrt(R_coil ** 2 + (z - z0[j] - L[j] / 2) ** 2))
            if j==coil:
                B_0tot += mu_0 * N_wires[j] * I / (2 * L[j]) * ((z - z0[j] + L[j] / 2) / np.sqrt(R_coil ** 2 + (z - z0[j] + L[j] / 2) ** 2) - (z - z0[j] - L[j] / 2) / np.sqrt(R_coil ** 2 + (z - z0[j] - L[j]/ 2) ** 2))
    return 1e4 * B_0tot


def B_coil(z, I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,I14,I15,I16,I17,I18,I19): # magnetic field all coils

    I_coil=np.array([I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,I14,I15,I16,I17,I18,I19])
    z0 = np.empty([coils])  # center of the coils
    N_wires = np.empty([coils])  # number of wires per layer for each coil
    M = np.empty([coils])  # number of wire layers for each coil
    B=0
    for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
        if j == 0:
            z0[j] = L[j] / 2   # z0 for the first Zeeman coil
        if j == 1 or j == coils - 1:
            z0[j] = z0[j - 1] + L[j]/ 2 + L[j] / 2 + dist_coils_large  # z0 for second and last Zeeman coil
        if j != 0 and j != 1 and j != coils - 1:
            z0[j] = z0[j - 1] + L[j] / 2 + L[j]/ 2 + dist_coils_small  # z0 for all the other Zeeman coils
    for p in range(0, coils):  # loop over all Zeeman slower coils to obtain the length of all coils
        N_wires[p] = L[p]/ d_wire  # number of wires in each layer
        M[p] = np.abs(round(N_coil[p] / N_wires[p], 0))  # number of layers
    M = M.astype(int)
    for j in range(0, coils):  # loop over Zeeman coils
        for mi in range(1, M[j] + 1):  # loop over all layers of each coil
            R_coil=R + mi * b_wire
            B += mu_0*N_wires[j]*I_coil[j]/(2*L[j]) \
                 *((z-z0[j]+L[j]/2)/np.sqrt(R_coil**2+(z-z0[j]+L[j]/2)**2) - (z-z0[j]-L[j]/2)/np.sqrt(R_coil**2+(z-z0[j]-L[j]/2)**2))
    return 1e4*B

fig, ax = plt.subplots()
print("plotting")

##### Enter all parameters here #####
coils = 19 #number of coils
L = np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]) #lengths of the coils
N_coil=np.array([320,290,260,260,270,260,250,250,240,230,220,190,190,160,160,140,130,130,100]) #windings of the coils
L_field=1.15 # length of the field (not of the slower!)
num = 5000 #number of points
R= 0.045/2 # inner radius of Zeeman-coils in m (not approved)
d_wire=0.005# length of the wire in m
b_wire=0.001# thickness of the wire in m
dist_coils_small = 0.002 #distance of coils (small)
dist_coils_large = 0.004 #distance of coils (large)

mu_0=4*np.pi*1e-7 # magnetic field constant

with open("fields/B(z)_1_0m_full_B.txt", "r") as g:  # plot ideal magnetic field
    lines = g.readlines()
    xnew = np.asarray([float(line.split(";")[0]) for line in lines])
    ynew = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(xnew+0.5, ynew, label="Ideal slower field of length {}m".format(round(L_field,2)))

L_slower = xnew[-1]+L_field
pos=np.linspace(0,L_field,num=num)

popt, pcov = curve_fit(B_coil, xnew+0.5, ynew,method="trf",bounds=(0,24.5)) #fit

for i in range(coils):
    print(round(popt[i],3),end=",")
print(" ")
plt.plot(pos,B_coil(pos,*popt),label="Fit: realistic field of length {}m".format(round(L_field,2)))

for i in range(coils): #plot fields of single coils
    plt.plot(pos,B_coil_single(pos,popt[i],i),color="black")

#write fit to file
file=open("B(z)_fit.txt","w+")
for zpos in pos:
    file.write("{};{}\n".format(zpos-0.5, B_coil(zpos,*popt)))
file.close()

#Plotting
plt.xlabel("Position in m", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.grid()
plt.rcParams.update({'font.size': 22})
xticks = ax.xaxis.get_major_ticks()
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
ax.yaxis.get_major_ticks()[1].label1.set_visible(False)

plt.legend(prop={'size': 15})
plt.show()



