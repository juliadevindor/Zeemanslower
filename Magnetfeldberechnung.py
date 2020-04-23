##########################################################################################################################################
# returns the magnetic field B(z) of an experimental setup which contains 8 coils followed by a MOT consisting of an anti Helmholtz coil #
# author: Julia Winter                                                                                                                   #
##########################################################################################################################################

from random import *
import matplotlib.pyplot as plt
import numpy as np
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

# additional coils
#R_add1=0.0225
#R_add2=0.0158
#R_add3=0.0167
#N_add1=50
#N_add2=16.2
#N_add3=-16.7
#I_add=4.0
#L_add1=0.03
#L_add2=0.012
#L_add3=0.011
#z0_add1= -0.2+0.6093
#z0_add2= -0.085+0.5
#z0_add3= -0.071+0.5
#N_wires_add1 = L_add1/d_wire ###???
#M_add1 = round(N_add1/N_wires_add1,0)
#M_add1 = int(M_add1)

#N = np.array([740, 415, 386, 362, 336, 311, 290, 266]) # field like the one that has been measured
#I = np.array([4.1, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.2]) # field like the one that has been measured

#N = np.array([500, 450, 415, 386, 362, 336, 311, 290, 270, 266]) # field with 9 coils for testing
#I = np.array([4.1, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.2]) # field with 9 coils for testing

#L_real = np.array([0.034, 0.032, 0.035, 0.0344, 0.0342, 0.0353, 0.032, 0.033]) # real length values for Zeeman coils
#n=np.arange(200,750,1) # random sample of turns
#ii=np.arange(4.0,5.0,0.1) # random sample of currents

for i in range(0,100000): # loop for creating different magnetic fields
    if sample_count==samples: break # break the loop if number of samples is reached

    coils = 10 #randrange(8, 13)  # random number of coils
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
    # B_sum=np.empty([num]) # for testing
    N_wires = np.empty([coils])  # number of wires per layer for each coil
    M = np.empty([coils])  # number of wire layers for each coil
    L = np.empty([coils])  # length of coils
    array_coils = np.array([coils])  # write coil number to array for file name

    L_slower = (coils-3) * dist_coils_small + 2 * dist_coils_large # length of slower
    L_ges = R_HH_outer + L_slower + dist_oven_slower+dist_slower_MOT#  total length
   # L_ges =L_slower + dist_oven_slower+0.2#  total length

    #N_rand=np.random.choice(n,coils) # pick random value for the turns of each Zeeman coil
    #I_rand=np.random.choice(ii,coils) # pick random value for the current of each Zeeman coil
    #N_sorted=np.sort(N_rand) # sorting
    #I_sorted=np.sort(I_rand)
    #N=N_sorted[::-1] # sorting sothat the magnetic field becomes decreasing
    #I=I_sorted[::-1]

    N = np.array([680,560,520,520,480,480,440,320,280,200]) # field like the one that has been measured
    I = np.array([4.1,4.8,4.8,4.8,4.8,4.8,4.8,4.8,4.8,4.2]) # field like the one that has been measured
    L = np.array([0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04])  # real length values for Zeeman coils
    #N = [700, 450, 430, 420, 415, 386, 362, 336, 311, 270, 250] # so klappt es mit 11 spulen
    #I = [4.1, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.2]
    #L = np.array([0.034, 0.033, 0.033, 0.033, 0.032, 0.035, 0.0344, 0.0342, 0.0353, 0.032, 0.033])  # real length values for Zeeman coils

    for p in range(0, coils): # loop over all Zeeman slower coils to obtain the length of all coils
    #    if coils>8: # if there are more than 8 coils
    #        if p>0 and p<coils-7: # put additional coils between first and second coil
    #            L[p] = 0.034
    #        if p==0:
    #            L[p]=L_real[p]
    #        else:
    #            L[p]=L_real[p+(8-coils)]
    #    else:
    #        L[p] = L_real[p]
        N_wires[p] = L[p]/d_wire # number of wires in each layer
        M[p] = round(N[p]/N_wires[p],0) # number of layers
    M=M.astype(int)

    for q in range(0, coils): # loop over coils to determine the total length and the length of the slower
        L_ges += L[q] # add all coils
        L_slower += L[q]
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

    all_vals = np.concatenate([N, I, M, dist, array_coils], axis=0) # concatenate all important setup parameters
    filename0 = str(all_vals[0:coils]) # prepare name of the txt file
    filename1 = str(all_vals[coils:2*coils])
    filename2 = str(all_vals[2*coils:3*coils])
    filename3 = str(all_vals[3*coils])
    filename4 = str(all_vals[3*coils+1])
    filename0 = filename0.replace("]", "")
    filename1 = filename1.replace("]", "")
    filename2 = filename2.replace("]", "")
    filename0 = filename0.replace("[", "")
    filename1 = filename1.replace("[", "")
    filename2 = filename2.replace("[", "")
    filename3 = filename3.replace("]", "")
    filename3 = filename3.replace("[", "")
    filename0 = filename0.replace(".", "")
    filename1 = filename1.replace(".", "")
    filename2 = filename2.replace(".", "")
    filename3 = filename3.replace(".", "")
    filename4 = filename4.replace(".0", "")
    filename0=filename0.replace(" ", "_")
    filename1=filename1.replace(" ", "_")
    filename2=filename2.replace(" ", "_")
    filename0=filename0.replace("__", "_")
    filename1=filename1.replace("__", "_")
    filename2=filename2.replace("__", "_")

    file = open("sim_setup/magnetic_fields/Bfield_N_{}_I_{}_M_{}_length_{}_coils_{}.txt".format(filename0, filename1, filename2, filename3, filename4),"w+")  # open file
    print("i={}, N={}, I={}, M={}, length={}, coils={}".format(i, filename0, filename1, filename2, filename3, filename4))

    HHfile = open("HH_coils.txt","w+")  # open file for HH coils
    B_HHtot=np.empty([num]) # array must be empty

    for o in range(0, num):  # loop over z (along the beam axis)
        print(o)
        B = 0  # startvalue for magnetic field
        if o == 0: # set initial value of z
            z[o] = 0.0
        else:
            if duration=="long":
                z[o] = z[o-1] + 0.0001#0.0001374#  # step in z
            else:
                z[o] = z[o-1] + 0.01
        for j in range(0, coils):  # loop over Zeeman coils
            for mi in range(1, M[j]+1): # loop over all layers of each coil
                #if(z[o]<0.6):
                B += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j], z[o]) # calculate magnetic field and add to total magnetic field
                #else:
                #    B=0
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

        #for mi_add in range(0,M_add1):
        #    B_9tot[o]  += B_coil(I_add, N_add1, L_add1, R_add1, z0_add1,z[o])
        #B_10tot[o] += B_coil(I_add, N_add2, L_add2, R_add2, z0_add2,z[o])
        #B_11tot[o] += B_coil(I_add, N_add3, L_add3, R_add3, z0_add3,z[o])

        for mi_HH in range(0,M_HH): # loop over all layers of HH coils
            B_HHtot[o] += (B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, I_HH, N_HH/M_HH, d_HH) + B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, -I_HH, N_HH/M_HH, -d_HH))  # calculating HH field
        HHfile.write("{};{}\n".format(z[o]-L_ges,B_HHtot[o]*1e4))

        #B+= B_HHtot[o] # add HH coils to total magnetic field

        #B+= B_9tot[o]
        #B+= B_10tot[o]
        #B+= B_11tot[o]

        B_tot[o] = B  # store total magnetic field in B_tot

        if z[o]<=0.2+L_ges:
            file.write(str(round(z[o]-L_ges, 5)))  # shift by 0.5m so that Andreas programm can use the txt file normally
            file.write(";")
            file.write(str(round(B_tot[o]*1e4, 5))) # Tesla -> Gauss
            file.write("\n")

    file.close()
    HHfile.close()
    #plt.plot(z,B_tot*1e4,label="my Bfield")

    Bmax = max(B_tot) # look for maximum value of the magnetic field
    index = np.where(B_tot == Bmax) # at which z is Bmax
    #grad = np.gradient(B_tot) # gradient of the magnetic field
    if Bmax * 1e4 > 550:  # maximum B should be greater than 550 Gauss
        if z[index] - L_ges > -(L_ges - 0.05) and z[index] - L_ges < -(L_ges - 0.25):  # maximum B should be in specific range (at the beginning of the slower)
            for i in range(0, len(B_tot) - 1): # loop over all point of magnetic field
                if z[i] - L_ges > z[index] - L_ges and z[i] - L_ges < 0.0:  # after maximum B and before 0
                    if B_tot[i] > B_tot[i + 1]: # B should be decreasing
                        hey=1
                        # if func1(B_tot,i,z)+500<=func2(s0,B_tot[i]) : hey=1 #B_tot in T?
                        # else:
                        # print("func1>func2 for func1=", func1(B_tot,i,z),"func2=", func2(s0,B_tot[i]),"z=", z[i]-0.5)
                        #    if os.path.isfile("sim_setup/magnetic_fields/Bfield_NI_{}_M_{}_length_{}.txt".format(filename1, filename2, filename3)): os.remove("sim_setup/magnetic_fields/Bfield_NI_{}_M_{}_length_{}.txt".format(filename1, filename2, filename3))
                    # add else conditions so that txt file containing magnetic field is deleted if it does not fulfill all the conditions
      #              else:
     #                   if os.path.isfile(
     #                       "sim_setup/magnetic_fields/Bfield_N_{}_I_{}_M_{}_length_{}_coils_{}.txt".format(filename0, filename1, filename2, filename3, filename4)): os.remove(
     #                       "sim_setup/magnetic_fields/Bfield_N_{}_I_{}_M_{}_length_{}_coils_{}.txt".format(filename0, filename1, filename2, filename3, filename4))
    #    else:
    #        if os.path.isfile(
    #            "sim_setup/magnetic_fields/Bfield_N_{}_I_{}_M_{}_length_{}_coils_{}.txt".format(filename0, filename1, filename2, filename3, filename4)): os.remove(
    #            "sim_setup/magnetic_fields/Bfield_N_{}_I_{}_M_{}_length_{}_coils_{}.txt".format(filename0, filename1, filename2, filename3, filename4))

   # else:
   #     if os.path.isfile("sim_setup/magnetic_fields/Bfield_N_{}_I_{}_M_{}_length_{}_coils_{}.txt".format(filename0, filename1, filename2, filename3, filename4)): os.remove("sim_setup/magnetic_fields/Bfield_N_{}_I_{}_M_{}_length_{}_coils_{}.txt".format(filename0, filename1, filename2, filename3, filename4))

    if os.path.isfile(
        "sim_setup/magnetic_fields/Bfield_N_{}_I_{}_M_{}_length_{}_coils_{}.txt".format(filename0, filename1, filename2,
                                                                                        filename3, filename4)): sample_count+=1 # add 1 to sample_count if the magnetic field is saved
print("plotting")
#plt.plot(z,B_1tot*1e4, color="orange")
#plt.plot(z,B_2tot*1e4, color="orange")
#plt.plot(z,B_3tot*1e4, color="orange")
#plt.plot(z,B_4tot*1e4, color="orange")
#plt.plot(z,B_5tot*1e4, color="orange")
#plt.plot(z,B_6tot*1e4, color="orange")
#plt.plot(z,B_7tot*1e4, color="orange")
#plt.plot(z,B_8tot*1e4, color="orange")
#plt.plot(z,B_9tot*1e4, color="orange")
#plt.plot(z,B_10tot*1e4, color="orange")
#plt.plot(z,B_11tot*1e4, color="orange")

#plt.plot(z,B_HHtot*1e4, color="orange")
#plt.plot(z-L_ges,np.gradient(B_HHtot*1e4, z-L_ges), color="green", label="ha")
#plt.plot(z-L_ges,B_HH1*1e4, color="black")
#plt.plot(z-L_ges,B_HH2*1e4, color="grey")
#print("MAX", max(np.gradient(B_HHtot*1e4, z-L_ges)), min(np.gradient(B_HHtot*1e4, z-L_ges)))
#print("z", z[np.where(np.gradient(B_HHtot*1e4, z-L_ges)==max(np.gradient(B_HHtot*1e4, z-L_ges)))]-0.5, z[np.where(np.gradient(B_HHtot*1e4, z-L_ges)==min(np.gradient(B_HHtot*1e4, z-L_ges)))]-0.5)

path = "C:/Users/jwinte02/PycharmProjects/atomic_levels_sim-master/sim_setup/magnetic_fields/*.txt" # path for plotting
for filename in glob.glob(path): # plot all txt files in folder magnetic_fields
    with open(filename, 'r') as f:
        filename=filename.replace("__","_")
        name_magn_field = filename.split("_")
        del name_magn_field[0:6]
        numberofcoils = name_magn_field[-1]
        numberofcoils = numberofcoils.replace(".txt", "")
        numberofcoils = int(numberofcoils)
        N = name_magn_field[0:numberofcoils]
        I = name_magn_field[numberofcoils+1:2*numberofcoils+1]
        M = name_magn_field[2*numberofcoils+2:3*numberofcoils+2]
        dist = name_magn_field[3*numberofcoils+3]
        dist=dist.replace(".txt", "")
        I=[w.replace("4", "4.") for w in I]
        dist=dist.replace("0", "0.",1)
        N = str(N)
        I = str(I)
        M = str(M)
        numberofcoils = str(numberofcoils)
        N = N.replace("[", "")
        N = N.replace("]", "")
        I = I.replace("[", "")
        I = I.replace("]", "")
        M = M.replace("[", "")
        M = M.replace("]", "")
        M = M.replace("'", "")
        N = N.replace("'", "")
        I = I.replace("'", "")
        lines = f.readlines()
        x = np.asarray([float(line.split(";")[0]) for line in lines])
        y = np.asarray([float(line.split(";")[1]) for line in lines])
        plt.plot(x+L_ges,y, label="N={}, \n I={}, \n M={}, \n dist={}, coils={}".format(N, I, M, dist, numberofcoils))#, color="blue") #x+0.6093

with open("sim_setup/example_magnetic_field.txt","r") as g: # plot measured magnetic field
    lines = g.readlines()
    xnew = np.asarray([float(line.split(";")[0]) for line in lines])
    ynew = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(xnew+0.5, ynew, label="andreas", color="black") #x+0.6093

    B_zeeman=slower_field(z, L_ges)
    #for i in range(np.where(B_zeeman==0)[0][0], len(B_zeeman)-1):
    #    B_zeeman[i]+=B_HHtot[i]*1e4

    #plt.plot(z,B_zeeman*1e4, label="slower field", color="red")
    #print(z-L_ges, slower_field(z))

    file = open("zeeman_field.txt", "w+")
    for i in range(0, len(B_zeeman)-1):
        file.write(str(z[i]-L_ges))
        file.write(";")
        file.write(str(round(B_zeeman[i]*1e4, 4)))
        file.write("\n")

    file.close()

plt.grid()
plt.legend(prop={'size': 6})
plt.show()



