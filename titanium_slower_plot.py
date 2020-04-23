import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import os
import glob

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


def slower_field(pos, L0):
  # length of slower in m
    B0 = 0.08 #0.19185 # y-Achsenabschnitt
    Bbias = -1.478*1e-9  # shiftet nur den Plot nach oben/ unten
    B = np.empty(len(pos))

    for i in range(0,len(pos)):
        if 1-pos[i]/L0 < 0: # to avoid "nan"
            B[i]=0
        else:
            B[i] = B0 * np.sqrt(1-pos[i]/L0) + Bbias

    return B

# constants
mu_B=9.274*1e-24 #J/T
hbar=1.055*1e-34 #Js
k_var=2*np.pi/(671*1e-9) #m
Gamma=2*np.pi*5.87*1e6 #Hz
m_Li=6.015*1.66*1e-27 #kg
mu_0=4*np.pi*1e-7 # magnetic field constant
# parameters of the experimental setup
duration="short" #"short"

if duration=="long":
    num=11001 # number of steps along the vacuum tube/ beam axis z (z=-0.6 to 0.2)
else:
    num=110
#coils=10 # number of coils
samples=1 # number of samples (magnetic fields) to be generated
sample_count=0
R= 0.043 # inner radius of Zeeman-coils in m (not approved)
d_wire=0.001# thickness of the wire in m
#b_wire=0.005
s0 = 5  # saturation parameter

# HH coils
R_HH=0.078 # inner radius of the HH coils (should be 0.087m)
R_HH_outer=0.117 # outer radius of the HH coils
I_HH=40 # current in A
N_HH=99 # number of  turns of each HH coil
d_HH=0.0493 # distance from one coil to the MOT center
d_wire_HH=0.001 # thickness of the wire in m
M_HH=4 # 4 layers of wire for each HH coil

# distances
dist_oven_slower=0.0811
dist_slower_MOT=0.0 # not approved
dist_coils_small=0.002
dist_coils_large=0.004

coils = 8 #randrange(8, 13)  # random number of coils
    # initialize arrays
B_tot = np.empty([num])  # total magnetic field of all Zeeman coils
z = np.empty([num])  # position on the beam axis for plotting (starting from coil 1 next to the oven)
z0 = np.empty([coils])  # center of the coils
B_HHtot = np.empty([num])  # magnetic field of the HH coils
B_1tot = np.empty([num])  # magnetic field of the first coil
B_2tot = np.empty([num])  # magnetic field of the second coil
B_3tot = np.empty([num])  # magnetic field of the third coil
B_4tot = np.empty([num])  # magnetic field of the fourth coil
B_5tot = np.empty([num])  # magnetic field of the fifth coil
B_6tot = np.empty([num])  # magnetic field of the sixth coil
B_7tot = np.empty([num])  # magnetic field of the seventh coil
B_8tot = np.empty([num])  # magnetic field of the eighth coil
B_9tot = np.empty([num])  # magnetic field of the first additional coil
B_10tot = np.empty([num])  # magnetic field of the second additional coil
B_11tot = np.empty([num])  # magnetic field of the third additional coil
# B_sum=np.empty([num]) # for testing
N_wires = np.empty([coils])  # number of wires per layer for each coil
M = np.empty([coils])  # number of wire layers for each coil
L = np.empty([coils])  # length of coils
array_coils = np.array([coils])  # write coil number to array for file name

L_slower = (coils-3) * dist_coils_small + 2 * dist_coils_large # length of slower
L_ges = R_HH_outer + L_slower + dist_oven_slower+dist_slower_MOT#  total length

N = np.array([647,415,386,362,336,311,268,266]) # field like the one that has been measured
I = np.array([4.1,4.8,4.8,4.8,4.8,4.8,4.8,4.2]) # field like the one that has been measured
L = np.array([0.034, 0.032, 0.035, 0.0344, 0.0342, 0.0353, 0.032, 0.033])  # real length values for Zeeman coils

for p in range(0, coils): # loop over all Zeeman slower coils to obtain the length of all coils
    N_wires[p] = L[p]/d_wire # number of wires in each layer
    M[p] = round(N[p]/N_wires[p],0) # number of layers
M=M.astype(int)

for q in range(0, coils): # loop over coils to determine the total length and the length of the slower
    L_ges += L[q] # add all coils
    L_slower += L[q]
L_ges=L_ges+5*dist_coils_small+2*dist_coils_large
z0_HH = L_ges  # center of the HH coils (not approved as the distance MOT-slower is missing)
dist = np.array([round(L_slower,4)])
print("Length of slower: ",L_slower, "Total length: ", L_ges)

for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
    if j==0:
        z0[j] = L[j]/2 + dist_oven_slower # z0 for the first Zeeman coil
    if j==1 or j==coils-1:
        z0[j] = z0[j-1] + L[j]/2 + L[j-1]/2 + dist_coils_large # z0 for second and last Zeeman coil
    if j!=0 and j!=1 and j!=coils-1:
        z0[j] = z0[j-1] + L[j]/2 + L[j-1]/2 + dist_coils_small # z0 for all the other Zeeman coils

B_HHtot=np.empty([num]) # array must be empty

for o in range(0, num):  # loop over z (along the beam axis)
    B = 0  # startvalue for magnetic field
    if o == 0:  # set initial value of z
        z[o] = 0.0
    else:
        if duration == "long":
            z[o] = z[o - 1] + 0.0001  # 0.0001374#  # step in z
        else:
            z[o] = z[o - 1] + 0.01
    for j in range(0, coils):  # loop over Zeeman coils
        for mi in range(1, M[j] + 1):  # loop over all layers of each coil
            B += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j],z[o])  # calculate magnetic field and add to total magnetic field

            if j == 0:
                B_1tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 1:
            #    B_2tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 2:
            #    B_3tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 3:
            #    B_4tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 4:
            #    B_5tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 5:
            #    B_6tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 6:
            #    B_7tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 7:
            #    B_8tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 8:
            #    B_9tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 9:
            #    B_10tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            #if j == 10:
            #    B_11tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])

    for mi_HH in range(0, M_HH):  # loop over all layers of HH coils
        B_HHtot[o] += (B_HHcoils(z[o], z0_HH, R_HH + 24 * d_wire_HH, I_HH, N_HH / M_HH, d_HH) + B_HHcoils(z[o], z0_HH, R_HH + 24 * d_wire_HH,-I_HH, N_HH / M_HH,-d_HH))  # calculating HH field

    #B+= B_HHtot[o] # add HH coils to total magnetic field

    B_tot[o] = B  # store total magnetic field in B_tot
    print(z[o]-0.5,";",B*1e4)

#plt.plot(z-0.5,B_HHtot*1e4,label="HH simulation")
plt.plot(z-0.5,B_tot*1e4,label="slower+HH simulation")
#plt.plot(z-0.5,B_1tot*1e4,label="coil1 simulation")

#with open("sim_setup/example_magnetic_field - Kopie.txt", "r") as f:
#    lines = f.readlines()
#    x = np.asarray([float(line.split(";")[0]) for line in lines])
#    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #plt.plot(x, y, label="Andis field")

with open("Titanium slower/titanium_slower_mot_on_new_shifted.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split()[0]) for line in lines])
    y = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot(x, y, label="MOT on NEW shifted")

with open("Titanium slower/titanium_slower_mot_off.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split()[0]) for line in lines])
    y = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot(-x/100+0.182, y, label="MOT off shifted")

with open("Measurements/Feld_MOT_Slower.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split()[0]) for line in lines])
    y = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot((x)/100-0.45, y*10, label="Messung Stefan ges")

plt.grid()
plt.legend()
plt.show()

