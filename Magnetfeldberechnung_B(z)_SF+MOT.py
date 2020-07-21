##########################################################################################################################################
# returns the magnetic field B(z) of an experimental setup which contains 8 coils followed by a MOT consisting of an anti Helmholtz coil #
# author: Julia Winter                                                                                                                   #
##########################################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import scipy.constants as scc
from scipy.optimize import curve_fit

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
    return mu_0*I_HHcoils*N_HHcoils/(2*np.pi)*(-a_HHcoils)/(zpos*np.sqrt((zpos+R_HHcoils)**2+(a_HHcoils)**2))*(-int1[0]+(R_HHcoils**2+(zpos)**2+(a_HHcoils)**2)/((R_HHcoils-zpos)**2+(a_HHcoils)**2)*int2[0])

fig, ax = plt.subplots()
print("plotting")

coils = 16
shift=550
dist=0#0.04
L = np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])
N = np.array([320,270,250,240,250,250,230,200,200,170,170,150,150,150,100,150])
I = np.array([14.0,14.0,13.605,13.16,11.847,10.889,10.623,10.61,9.448,8.974,7.221,5.509,2.895,-0.862,-5.944,-17.341])
L_field=1.5 #0.93 + dist#+time for increase of  magnetic field
num = 5000
mu_0=4*np.pi*1e-7 # magnetic field constant
R= 0.045/2 # inner radius of Zeeman-coils in m (not approved)
d_wire=0.005# length of the wire in m
b_wire=0.001# thickness of the wire in m
dist_coils_small = 0.002
dist_coils_large = 0.002 #0.004

R_HH=0.078 # inner radius of the HH coils (should be 0.087m)
R_HH_outer=0.117 # outer radius of the HH coils
I_HH=40 # current in A
N_HH=99 # number of  turns of each HH coil
d_HH=0.0493 # distance from one coil to the MOT center
d_wire_HH=0.001 # thickness of the wire in m
M_HH=4 # 4 layers of wire for each HH coil
L_MOT=0.6 # relevant length of MOT field
z0_HH = 1.0  # center of the HH coils (not approved as the distance MOT-slower is missing)

pos=np.linspace(-dist,L_field-dist,num=num)
#plt.plot(pos,B_coil(pos,*I_coil),label="Real field of length {}m".format(round(L_field,2)))
B_tot = np.empty([num])  # total magnetic field of all Zeeman coils
z = np.empty([num])  # position on the beam axis for plotting (starting from coil 1 next to the oven)
z0 = np.empty([coils])  # center of the coils
B_HHtot = np.empty([num])  # magnetic field of the HH coils
N_wires = np.empty([coils])  # number of wires per layer for each coil
M = np.empty([coils])  # number of wire layers for each coil
for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
    if j == 0:
        z0[j] = L[j] / 2 - dist  # z0 for the first Zeeman coil
    if j == 1 or j == coils - 1:
        z0[j] = z0[j - 1] + L[j] / 2 + L[j] / 2 + dist_coils_large  # z0 for second and last Zeeman coil
    if j != 0 and j != 1 and j != coils - 1:
        z0[j] = z0[j - 1] + L[j] / 2 + L[j] / 2 + dist_coils_small  # z0 for all the other Zeeman coils
for p in range(0, coils):  # loop over all Zeeman slower coils to obtain the length of all coils
    N_wires[p] = L[p] / d_wire  # number of wires in each layer
    M[p] = np.abs(round(N[p] / N_wires[p], 0))  # number of layers
M = M.astype(int)
for o in range(0, num):  # loop over z (along the beam axis)
    B = 0  # startvalue for magnetic field
    if o == 0: # set initial value of z
        z[o] = 0.0
    else:
        z[o]=z[o-1]+L_field/num
    for j in range(0, coils):  # loop over Zeeman coils
        for mi in range(0, M[j]): # loop over all layers of each coil
            R_layer=R+mi*b_wire
            B += B_coil(I[j], N_wires[j], L[j], R_layer, z0[j], z[o]) # calculate magnetic field and add to total magnetic field
    for mi_HH in range(0,M_HH): # loop over all layers of HH coils
        B_HHtot[o] += (B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, -I_HH, N_HH/M_HH, d_HH) + B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, I_HH, N_HH/M_HH, -d_HH))  # calculating HH field

    B+= B_HHtot[o] # add HH coils to total magnetic field
    B_tot[o] = B

plt.plot(z,B_tot*1e4,label="Real spin-flip slower L=0.8m + MOT",color="blue")

with open("fields/B(z)_0_8m_SF_long.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x + 0.5, y, label="Ideal spin-flip slower L=0.8m", color="red")
with open("Measurements/Feld_MOT_01_12_2017_40A.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split()[0]) for line in lines])
    y = np.asarray([float(line.split()[1]) for line in lines])
    #ax.plot(x/100 +0.863, y*10, label="MOT measurement", color="green")
with open("sim_setup/example_magnetic_field.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y*1e4, label="FIELD", color="green")
file=open("B(z)_slower_MOT.txt","w+")
i=0
for zpos in pos:
    file.write("{};{}\n".format(zpos-0.5+dist, B_tot[i]*1e4))
    i+=1
file.close()
with open("B(z)_slower_MOT.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y*1e4, label="measurement", color="green")

plt.xlabel("Position in m", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.grid()
plt.rcParams.update({'font.size': 22})
xticks = ax.xaxis.get_major_ticks()
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#ax.spines['left'].set_position('zero')
#ax.spines['bottom'].set_position('zero')
#ax.yaxis.get_major_ticks()[1].label1.set_visible(False)


plt.legend(prop={'size': 15})
plt.show()



