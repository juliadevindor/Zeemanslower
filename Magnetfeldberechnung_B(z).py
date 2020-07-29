##########################################################################################################################################
# returns the magnetic field B(z) of an experimental setup which contains 8 coils followed by a MOT consisting of an anti Helmholtz coil #
##########################################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import scipy.constants as scc

hPlanck=scc.h #6.62607004e-34 Js
muB=scc.physical_constants['Bohr magneton'][0] # 9.274009994e-24 J/T

def B_coil(I_coil, N_coil, L_coil, R_coil, z0_coil, pos): # magnetic field of a single coil
    return mu_0*N_coil*I_coil/(2*L_coil)*((pos-z0_coil+L_coil/2)/np.sqrt(R_coil**2+(pos-z0_coil+L_coil/2)**2) - (pos-z0_coil-L_coil/2)/np.sqrt(R_coil**2+(pos-z0_coil-L_coil/2)**2)) # magnetic field of a single coil

def integrand_K(x,k): # integrand for calculating the magnetic field of the HH coils
    return 1/(np.sqrt(1-k*np.sin(x)**2))

def integrand_E(x,k):  # integrand for calculating the magnetic field of the HH coils
    return np.sqrt(1-k*np.sin(x)**2)

def B_HHcoils(pos, z0_HHcoils, R_HHcoils, I_HHcoils, N_HHcoils, a_HHcoils): # magnetic field of one HH coil
    k=4*R_HHcoils*(pos-z0_HHcoils)/((R_HHcoils+pos-z0_HHcoils)**2+(a_HHcoils)**2) #=k^2
    int1=quad(integrand_K, 0, 0.5*np.pi, args=(k))
    int2=quad(integrand_E, 0, 0.5*np.pi, args=(k))
    zpos=pos-z0_HHcoils
    return mu_0*I_HHcoils*N_HHcoils/(2*np.pi) * (-a_HHcoils)/(zpos*np.sqrt((zpos+R_HHcoils)**2+(a_HHcoils)**2))*(-int1[0]+(R_HHcoils**2+(zpos)**2+(a_HHcoils)**2)/((R_HHcoils-zpos)**2+(a_HHcoils)**2)*int2[0])

# constants
mu_B=9.274*1e-24 #J/T
hbar=1.055*1e-34 #Js
k_var=2*np.pi/(671*1e-9) #m
Gamma=2*np.pi*5.87*1e6 #Hz
m_Li=6.015*1.66*1e-27 #kg
mu_0=4*np.pi*1e-7 # magnetic field constant

num=10000 #number of points
R= 0.045/2 # inner radius of Zeeman-coils in m
d_wire=0.001 # thickness of the wire in m
b_wire=0.005 # width of the wire in m
s0 = 5  # saturation parameter

# HH coils
R_HH=0.078 # inner radius of the HH coils (should be 0.087m)
R_HH_outer=0.117 # outer radius of the HH coils
I_HH=40 # current in A
N_HH=99 # number of  turns of each HH coil
d_HH=0.0493 # distance from one coil to the MOT center
d_wire_HH=0.001 # thickness of the wire in m
M_HH=4 # 4 layers of wire for each HH coil
L_MOT=0.6 # relevant length of MOT field

# distances in m
dist_oven_slower=0
dist_slower_MOT=0.2 # not approved
dist_coils_small=0.002
dist_coils_large=0.004

coils = 16 #number of coils
print("Number of coils:", coils)

# initialize arrays
B_tot = np.empty([num])  # total magnetic field of all Zeeman coils
B_HHtot = np.empty([num])  # magnetic field of the MOT
z = np.empty([num])  # position on the beam axis for plotting (starting from coil 1 next to the oven)
z0 = np.empty([coils])  # center of the coils
N_wires = np.empty([coils])  # number of wires per layer for each coil
M = np.empty([coils])  # number of wire layers for each coil
L = np.empty([coils])  # length of coils
array_coils = np.array([coils])  # write coil number to array for file name

L_slower = (coils-2) * dist_coils_small + 2 * dist_coils_large + dist_oven_slower#+dist_slower_MOT#  total length
print("L_slower w/o L_coils",L_slower)

N = np.array([320,270,240, 260,240, 250,230,200, 220,180,140, 200,100, 10,70,320])
I = np.array([14.0,14.0,14.0, 12.5,12.5, 10.7,10.7,10.7, 8.6,8.6,8.6, 4.2,4.2, -9.0,-9.0,-9.0])
L = np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])

for p in range(0, coils): # loop over all Zeeman slower coils to obtain the length of all coils
    N_wires[p] = L[p]/b_wire # number of wires in each layer
    M[p] = np.abs(round(N[p]/N_wires[p],0)) # number of layers
M=M.astype(int)

for q in range(0, coils): # loop over coils to determine the total length and the length of the slower
    L_slower += L[q]

z0_HH = L_slower+dist_slower_MOT  # center of the HH coils (not approved as the distance MOT-slower is missing)
L_ges=0.93 #total length of the B-field in m

for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
    if j==0:
        z0[j] = L[j]/2 + dist_oven_slower # z0 for the first Zeeman coil
    if j==1 or j==coils-1:
        z0[j] = z0[j-1] + L[j]/2 + L[j-1]/2 + dist_coils_large # z0 for second and last Zeeman coil
    if j!=0 and j!=1 and j!=coils-1:
        z0[j] = z0[j-1] + L[j]/2 + L[j-1]/2 + dist_coils_small # z0 for all the other Zeeman coils

file = open("magnetic_field_real.txt","w+")  # open file

for o in range(0, num):  # loop over z (along the beam axis)
    B = 0  # startvalue for magnetic field
    if o == 0: # set initial value of z
        z[o] = 0.0
    else:
        z[o]=z[o-1]+L_ges/num

    for j in range(0, coils):  # loop over Zeeman coils
        for mi in range(0, M[j]): # loop over all layers of each coil
            R_layer=R+mi*d_wire
            B += B_coil(I[j], N_wires[j], L[j], R_layer, z0[j], z[o]) # calculate magnetic field and add to total magnetic field

    for mi_HH in range(0,M_HH): # loop over all layers of HH coils
        B_HHtot[o] += (B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, I_HH, N_HH/M_HH, d_HH) + B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, -I_HH, N_HH/M_HH, -d_HH))  # calculating HH field

    #B+= B_HHtot[o] # add HH coils to total magnetic field

    B_tot[o] = B  # store total magnetic field in B_tot

    file.write(str(round(z[o]-0.5, 5)))  # shift by 0.5m so that Andreas programm can use the txt file normally
    file.write(";")
    file.write(str(round(B_tot[o]*1e4, 5))) # Tesla -> Gauss
    file.write("\n")

file.close()

fig, ax = plt.subplots()

print("plotting")

with open("magnetic_field_real.txt","r") as g: # plot measured magnetic field
    lines = g.readlines()
    xnew = np.asarray([float(line.split(";")[0]) for line in lines])
    ynew = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(xnew+0.5, ynew, label="Realistic slower field of length 0.8m")

with open("fields/B(z)_fit_0_8m_SF.txt","r") as g: # plot measured magnetic field
    lines = g.readlines()
    xnew = np.asarray([float(line.split(";")[0]) for line in lines])
    ynew = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(xnew+0.5, ynew, label="Fit slower field of length 0.8m")

plt.xlabel("Position in m", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.grid()
plt.rcParams.update({'font.size': 22})
xticks = ax.xaxis.get_major_ticks()
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.legend(prop={'size': 15})
plt.show()



