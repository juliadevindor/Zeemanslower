##########################################################################################################################################
# returns the magnetic field B(z) of an experimental setup which contains 8 coils followed by a MOT consisting of an anti Helmholtz coil #
# author: Julia Winter                                                                                                                   #
##########################################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import scipy.constants as scc
from scipy.optimize import curve_fit

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

num=1000
sample_count=0
R= 0.043 # inner radius of Zeeman-coils in m (not approved)
d_wire=0.001# thickness of the wire in m
b_wire=0.005
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

# distances
dist_oven_slower=0#0.08
dist_slower_MOT=0.1 # not approved
dist_coils_small=0.002
dist_coils_large=0.004

coils = 17# #randrange(8, 13)  # random number of coils
print("Number of coils:", coils)

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
B_12tot = np.empty([num])  # magnetic field of the third additional coil
B_13tot = np.empty([num])  # magnetic field of the third additional coil
N_wires = np.empty([coils])  # number of wires per layer for each coil
M = np.empty([coils])  # number of wire layers for each coil
L = np.empty([coils])  # length of coils
array_coils = np.array([coils])  # write coil number to array for file name

L_slower = (coils-2) * dist_coils_small + 2 * dist_coils_large + dist_oven_slower#+dist_slower_MOT#  total length
print("L_slower w/o L_coils",L_slower)
#N = np.array([320,320,270,270,250,220,200,170,160]) #np.array([800,700,650,600,550,450,350,300,250,200, 50,100,650]) # field like the one that has been measured
#I = np.array([15.625,11.9,11.668,11.541,10.997,11.352,10.995,10.396,10.122]) # field like the one that has been measured
#N = np.array([50,50,50,50,50,50,50,50,50]) #np.array([800,700,650,600,550,450,350,300,250,200, 50,100,650]) # field like the one that has been measured
#I = np.array([100.0,76.16286526,63.00735423,62.31878558,54.98645356,49.94757199,43.9799914,35.34623289,32.3912794]) # field like the one that has been measured
N = np.array([320,320,300,300,280,280,260,260,240,240,200,200,180,180,160,160,140]) #np.array([800,700,650,600,550,450,350,300,250,200, 50,100,650]) # field like the one that has been measured
I = np.array([25.0,25.0,22.567,22.013,22.554,21.628,22.247,21.185,21.738,20.478,22.937,21.253,21.413,19.358,18.239,16.134,11.441])
# field like the one that has been measured
L = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05,0.05,0.05,0.05,0.05, 0.05, 0.05, 0.05, 0.05,0.05,0.05,0.05])  # real length values for Zeeman coils

for p in range(0, coils): # loop over all Zeeman slower coils to obtain the length of all coils
    N_wires[p] = L[p]/d_wire # number of wires in each layer
    M[p] = np.abs(round(N[p]/N_wires[p],0)) # number of layers
M=M.astype(int)

for q in range(0, coils): # loop over coils to determine the total length and the length of the slower
    L_slower += L[q]

z0_HH = L_slower+dist_slower_MOT  # center of the HH coils (not approved as the distance MOT-slower is missing)
dist = np.array([round(L_slower,4)])
L_ges=L_slower#+L_MOT+dist_slower_MOT
print("Length of slower: ",L_slower, "total length of field", L_ges)


for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
    if j==0:
        z0[j] = L[j]/2 + dist_oven_slower # z0 for the first Zeeman coil
    if j==1 or j==coils-1:
        z0[j] = z0[j-1] + L[j]/2 + L[j-1]/2 + dist_coils_large # z0 for second and last Zeeman coil
    if j!=0 and j!=1 and j!=coils-1:
        z0[j] = z0[j-1] + L[j]/2 + L[j-1]/2 + dist_coils_small # z0 for all the other Zeeman coils
    print(z0[j])
file = open("magnetic_field_real.txt","w+")  # open file

for o in range(0, num):  # loop over z (along the beam axis)
    B = 0  # startvalue for magnetic field
    if o == 0: # set initial value of z
        z[o] = 0.0
    else:
        z[o]=z[o-1]+L_ges/num

    for j in range(0, coils):  # loop over Zeeman coils
        for mi in range(1, M[j]+1): # loop over all layers of each coil
            B += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j], z[o]) # calculate magnetic field and add to total magnetic field

            ##
            if j==0:
                B_1tot[o] += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j], z[o])
            if j==1:
                B_2tot[o] += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j],z[o])
            if j==2:
                B_3tot[o] += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j],z[o])
            if j==3:
                B_4tot[o] += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j],z[o])
            if j==4:
                B_5tot[o] += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j],z[o])
            if j==5:
                B_6tot[o] += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j],z[o])
            if j==6:
                B_7tot[o] += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j],z[o])
            if j==7:
                B_8tot[o] += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j], z[o])
            if j==8:
                B_9tot[o] += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j], z[o])
            if j == 9:
                B_10tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            if j == 10:
                B_11tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            if j == 11:
                B_12tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            if j == 12:
                B_13tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
            ##
    for mi_HH in range(0,M_HH): # loop over all layers of HH coils
        B_HHtot[o] += (B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, I_HH, N_HH/M_HH, d_HH) + B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, -I_HH, N_HH/M_HH, -d_HH))  # calculating HH field

    #B+= B_HHtot[o] # add HH coils to total magnetic field

    #if z[o]>0.6:
    #    B_tot[o]=0.0
    #else:
    B_tot[o] = B  # store total magnetic field in B_tot

    #if z[o]<=0.2+L_ges:
    file.write(str(round(z[o]-0.5, 5)))  # shift by 0.5m so that Andreas programm can use the txt file normally
    file.write(";")
    file.write(str(round(B_tot[o]*1e4, 5))) # Tesla -> Gauss
    file.write("\n")

file.close()

fig, ax = plt.subplots()

print("plotting")
#ax.plot(z,B_1tot*1e4, color="orange")
#ax.plot(z,B_2tot*1e4, color="orange")
#ax.plot(z,B_3tot*1e4, color="orange")
#ax.plot(z,B_4tot*1e4, color="orange")
#ax.plot(z,B_5tot*1e4, color="orange")
#ax.plot(z,B_6tot*1e4, color="orange")
#ax.plot(z,B_7tot*1e4, color="orange")
#ax.plot(z,B_8tot*1e4, color="orange")
#ax.plot(z,B_9tot*1e4, color="orange")
#ax.plot(z,B_10tot*1e4, color="orange")
#ax.plot(z,B_11tot*1e4, color="orange")
#ax.plot(z,B_12tot*1e4, color="orange")
#plt.plot(z-0.5,B_13tot*1e4, color="orange")

#plt.plot(z-0.5,B_HHtot*1e4, color="orange")

with open("magnetic_field_real.txt", 'r') as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x,y,label="Spin-flip field")
    ax.plot(x+0.5,y,".",label="Real slower field")

with open("B(z)_0_9m.txt","r") as g: # plot measured magnetic field
    lines = g.readlines()
    xnew = np.asarray([float(line.split(";")[0]) for line in lines])
    ynew = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(xnew+0.5, ynew, label="Ideal slower field of length 0.5m")

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



