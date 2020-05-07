##########################################################################################################################################
# returns the magnetic field B(z) of an experimental setup which contains 8 coils followed by a MOT consisting of an anti Helmholtz coil #
# author: Julia Winter                                                                                                                   #
##########################################################################################################################################

from random import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import scipy.constants as scc
import os
import glob

hPlanck=scc.h #6.62607004e-34 Js
muB=scc.physical_constants['Bohr magneton'][0] # 9.274009994e-24 J/T

def slower_field(pos, L0,v0,wavelength,omega, omega0):
  # length of slower L0 in m
    B0 = hPlanck*v0/(wavelength*muB) #0.08 #0.19185 # y-Achsenabschnitt
    #print(B0)
    Bbias = -3.5e-2 #0#(hPlanck/(2*np.pi))*(omega-omega0)/(muB)#-1.478*1e-9  # shiftet nur den Plot nach oben/ unten
    B = np.empty(len(pos))

    for i in range(0,len(pos)):
        #print(pos[i])
        #if pos[i]<0.0:
        #    B[i]=0
        if 1-pos[i]/L0 < 0: # to avoid "nan"
            B[i]=0
        elif pos[i]>L0:
            B[i]=0
        else:
            B[i] = B0 * np.sqrt(1-pos[i]/L0) + Bbias
    return B

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

num=10000
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
L_MOT=0.6 # relevant length of MOT field

# distances
dist_oven_slower=0.08
dist_slower_MOT=0.1 # not approved
dist_coils_small=0.002
dist_coils_large=0.004

coils = 13 #randrange(8, 13)  # random number of coils
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
N = np.array([750,700,600,550,500,400,300,250,200,150, 50,100,650]) #np.array([800,700,650,600,550,450,350,300,250,200, 50,100,650]) # field like the one that has been measured
I = np.array([4.8,4.8,4.8,4.8,4.8,4.8,4.8, 4.8,4.8,4.8,-4.8,-4.8,-4.8]) # field like the one that has been measured
L = np.array([0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,0.05, 0.05,0.05, 0.04])  # real length values for Zeeman coils

for p in range(0, coils): # loop over all Zeeman slower coils to obtain the length of all coils
    N_wires[p] = L[p]/d_wire # number of wires in each layer
    M[p] = np.abs(round(N[p]/N_wires[p],0)) # number of layers
M=M.astype(int)

for q in range(0, coils): # loop over coils to determine the total length and the length of the slower
    L_slower += L[q]

z0_HH = L_slower+dist_slower_MOT  # center of the HH coils (not approved as the distance MOT-slower is missing)
dist = np.array([round(L_slower,4)])
L_ges=L_slower+dist_slower_MOT+L_MOT
print("Length of slower: ",L_slower, "total length of field", L_ges)


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

   # if z[o]>L_slower:
   #     B_tot[o]=0.0
    #else:
    B_tot[o] = B  # store total magnetic field in B_tot

    #if z[o]<=0.2+L_ges:
    file.write(str(round(z[o]-0.5, 5)))  # shift by 0.5m so that Andreas programm can use the txt file normally
    file.write(";")
    file.write(str(round(B_tot[o]*1e4, 5))) # Tesla -> Gauss
    file.write("\n")

file.close()


   # Bmax = max(B_tot) # look for maximum value of the magnetic field
   # index = np.where(B_tot == Bmax) # at which z is Bmax
    #grad = np.gradient(B_tot) # gradient of the magnetic field
   # if Bmax * 1e4 > 550:  # maximum B should be greater than 550 Gauss
  #      if z[index] - L_ges > -(L_ges - 0.05) and z[index] - L_ges < -(L_ges - 0.25):  # maximum B should be in specific range (at the beginning of the slower)
  #          for i in range(0, len(B_tot) - 1): # loop over all point of magnetic field
  #              if z[i] - L_ges > z[index] - L_ges and z[i] - L_ges < 0.0:  # after maximum B and before 0
  #                  if B_tot[i] > B_tot[i + 1]: # B should be decreasing
  #                      hey=1

print("plotting")
#plt.plot(z-0.5,B_1tot*1e4, color="orange")
#plt.plot(z-0.5,B_2tot*1e4, color="orange")
#plt.plot(z-0.5,B_3tot*1e4, color="orange")
#plt.plot(z-0.5,B_4tot*1e4, color="orange")
#plt.plot(z-0.5,B_5tot*1e4, color="orange")
#plt.plot(z-0.5,B_6tot*1e4, color="orange")
#plt.plot(z-0.5,B_7tot*1e4, color="orange")
#plt.plot(z-0.5,B_8tot*1e4, color="orange")
#plt.plot(z-0.5,B_9tot*1e4, color="orange")
#plt.plot(z-0.5,B_10tot*1e4, color="orange")
#plt.plot(z-0.5,B_11tot*1e4, color="orange")
#plt.plot(z-0.5,B_12tot*1e4, color="orange")
#plt.plot(z-0.5,B_13tot*1e4, color="orange")

#plt.plot(z-0.5,B_HHtot*1e4, color="orange")
#plt.plot(z-L_ges,np.gradient(B_HHtot*1e4, z-L_ges), color="green", label="ha")
#plt.plot(z-L_ges,B_HH1*1e4, color="black")
#plt.plot(z-L_ges,B_HH2*1e4, color="grey")
#print("MAX", max(np.gradient(B_HHtot*1e4, z-L_ges)), min(np.gradient(B_HHtot*1e4, z-L_ges)))
#print("z", z[np.where(np.gradient(B_HHtot*1e4, z-L_ges)==max(np.gradient(B_HHtot*1e4, z-L_ges)))]-0.5, z[np.where(np.gradient(B_HHtot*1e4, z-L_ges)==min(np.gradient(B_HHtot*1e4, z-L_ges)))]-0.5)

with open("magnetic_field_real.txt", 'r') as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x,y,label="calculated field")

with open("sim_setup/example_magnetic_field.txt","r") as g: # plot measured magnetic field
    lines = g.readlines()
    xnew = np.asarray([float(line.split(";")[0]) for line in lines])
    ynew = np.asarray([float(line.split(";")[1]) for line in lines])
    #plt.plot(xnew, ynew, label="perfect field", color="black") #x+0.6093

L0=0.828#m
v0_0=1500 #m/s
v0_1=800 #m/s
v0_2=500 #m/s
wavelength=670.992e-9 #m (Gehm)
omega=446799923264221.4 #Hz (mein Wert) #446799900000000 #Hz (Andis Wert)
omega0=446.7896e12 #Hz (Gehm)

pos=np.linspace(0.0,L0+0.2,num=100)
B = np.empty(len(pos))
B=slower_field(pos,L0,v0_1,wavelength,omega,omega0)
#plt.plot(pos-0.4,B*1e4,".",label="v0=1000m/s")

plt.xlabel("Position in m", fontsize=15)
plt.ylabel("Magnetic field in Gauss", fontsize=15)
plt.grid()
plt.legend(prop={'size': 15})
plt.show()



