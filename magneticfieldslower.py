import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as scc
from scipy.integrate import quad
import os
import glob

k_var=2*np.pi/(671*1e-9) #m
Gamma=2*np.pi*5.87*1e6 #Hz
m_Li=6.015*1.66*1e-27 #kg
mu_0=4*np.pi*1e-7 # magnetic field constant

def slower_field(det,v_0,lambda_res,pos, L0):
  # length of slower in m
    B0 = h_bar*2*np.pi*v_0/(lambda_res*mb) # y-Achsenabschnitt
    Bbias = h_bar/mb *det  # shiftet nur den Plot nach oben/ unten
    B = np.empty(len(pos))
    for i in range(0,len(pos)):
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
    if zpos==0: return 0.0
    else: return mu_0*I_HHcoils*N_HHcoils/(2*np.pi) * (-a_HHcoils)/(zpos*np.sqrt((zpos+R_HHcoils)**2+(a_HHcoils)**2))*(-int1[0]+(R_HHcoils**2+(zpos)**2+(a_HHcoils)**2)/((R_HHcoils-zpos)**2+(a_HHcoils)**2)*int2[0])


if __name__ == '__main__':
    file=open("TEST.txt","w+")
    h_bar = scc.hbar  # Plank constant/2pi
    mb = scc.physical_constants['Bohr magneton'][0]  # Bohr magneton

    L0=0.61 #meters
    det=-130*1e6 #detuning omega-omega_0 in Hz
    v_0=1000#1470 #most probable velocity in beam of 6Li in m/s (value by Andreas)
    lambda_res= 670.977*1e-9#resonance wavelength in meters
    pos = np.arange(0, L0+0.4, 0.01)
    Bfield=np.empty([len(pos)])
    Bfield_ex=np.empty([5000])
    pos_ex=np.empty([5000])
    Bfield=slower_field(det,v_0,lambda_res,pos,L0)*1e4 #in Gauss
    for i in range(len(pos)):
        file.write("{};{}\n".format(pos[i]-0.95,Bfield[i]-50))

    with open("Measurements/MOT_new.txt", "r") as g:
        lines = g.readlines()
        x = np.asarray([float(line.split()[0]) for line in lines])
        y = np.asarray([float(line.split()[1]) for line in lines])
        for i in range(len(x)):
            file.write("{};{}\n".format((x[i]-38.5)/100, y[i]*10))
    file.close()

    with open("Test.txt", "r") as g:
        lines = g.readlines()
        xtest = np.asarray([float(line.split(";")[0]) for line in lines])
        ytest = np.asarray([float(line.split(";")[1]) for line in lines])
        #plt.plot(xtest,ytest)

    z=np.empty([len(pos)])
    B_HHtot=np.empty([len(pos)])
    B_tot=np.empty([len(pos)])
    B_1tot=np.empty([len(pos)])
    B_2tot=np.empty([len(pos)])
    B_3tot=np.empty([len(pos)])
    B_4tot=np.empty([len(pos)])
    B_5tot=np.empty([len(pos)])
    B_6tot=np.empty([len(pos)])
    B_7tot=np.empty([len(pos)])
    B_8tot=np.empty([len(pos)])
    B_9tot=np.empty([len(pos)])
    B_10tot=np.empty([len(pos)])


    coils=10
    L_ges=0.5
    N = np.array([980, 860, 820, 820, 760, 750, 600, 520, 500, 450])  # field like the one that has been measured
    I = np.array([4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8])  # field like the one that has been measured
    L = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])  # real length values for Zeeman coils
    M = np.empty([coils])
    R = 0.043  # inner radius of Zeeman-coils in m (not approved)
    d_wire = 0.001  # thickness of the wire in m
    s0 = 5  # saturation parameter
    N_wires = np.empty([coils])
    z0 = np.empty([coils])
    dist_oven_slower = 0.0811
    dist_coils_small = 0.002
    dist_coils_large = 0.004

    for j in range(0, coils):  # calculate z0(center of coils) for all Zeeman coils
        if j==0:
            z0[j] = L[j]/2 + dist_oven_slower # z0 for the first Zeeman coil
        if j==1 or j==coils-1:
            z0[j] = z0[j-1] + L[j]/2 + L[j-1]/2 + dist_coils_large # z0 for second and last Zeeman coil
        if j!=0 and j!=1 and j!=coils-1:
            z0[j] = z0[j-1] + L[j]/2 + L[j-1]/2 + dist_coils_small # z0 for all the other Zeeman coils
    for p in range(0, coils): # loop over all Zeeman slower coils to obtain the length of all coils
        N_wires[p] = L[p]/d_wire # number of wires in each layer
        M[p] = round(N[p]/N_wires[p],0) # number of layers
    M=M.astype(int)

    # HH coils
    R_HH = 0.078  # inner radius of the HH coils (should be 0.087m)
    R_HH_outer = 0.117  # outer radius of the HH coils
    I_HH = 40  # current in A
    N_HH = 99  # number of  turns of each HH coil
    d_HH = 0.0493  # distance from one coil to the MOT center
    d_wire_HH = 0.001  # thickness of the wire in m
    M_HH = 4  # 4 layers of wire for each HH coil
    z0_HH = 0.85

    #Bfield of coils
    for o in range(len(pos)):  # loop over z (along the beam axis)
        B = 0  # startvalue for magnetic field
        z[o]=pos[o]
        B_1tot[o]=0
        B_2tot[o]=0
        B_3tot[o]=0
        B_4tot[o]=0
        B_5tot[o] = 0
        B_6tot[o] = 0
        B_7tot[o] = 0
        B_8tot[o] = 0
        B_9tot[o] = 0
        B_10tot[o] = 0
        B_HHtot[o]=0

        for j in range(0, coils):  # loop over Zeeman coils
            for mi in range(1, M[j]+1): # loop over all layers of each coil
                B += B_coil(I[j], N_wires[j], L[j], R+mi*d_wire, z0[j], z[o]) # calculate magnetic field and add to total magnetic field
                if j == 0:
                    B_1tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
                if j == 1:
                    B_2tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
                if j == 2:
                    B_3tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
                if j == 3:
                    B_4tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
                if j == 4:
                    B_5tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
                if j == 5:
                    B_6tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
                if j == 6:
                    B_7tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
                if j == 7:
                    B_8tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
                if j == 8:
                    B_9tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])
                if j == 9:
                    B_10tot[o] += B_coil(I[j], N_wires[j], L[j], R + mi * d_wire, z0[j], z[o])

        for mi_HH in range(0,M_HH): # loop over all layers of HH coils
            B_HHtot[o] += (B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, I_HH, N_HH/M_HH, d_HH) + B_HHcoils(z[o], z0_HH, R_HH+24*d_wire_HH, -I_HH, N_HH/M_HH, -d_HH))  # calculating HH field
        B+= B_HHtot[o] # add HH coils to total magnetic field

        B_tot[o] = B  # store total magnetic field in B_tot

    plt.plot(z-0.85,B_tot*1e4)

    file=open("./SLOWER.txt","w+")
    for i in range(len(z)):
        file.write("{};{}\n".format(z[i]-0.85,B_tot[i]*1e4))
    file.close
    #plt.plot(z,B_1tot*1e4)
    #plt.plot(z,B_2tot*1e4)
    #plt.plot(z,B_3tot*1e4)
    #plt.plot(z,B_4tot*1e4)
    #plt.plot(z,B_5tot*1e4)
    #plt.plot(z,B_6tot*1e4)
    #plt.plot(z,B_7tot*1e4)
    #plt.plot(z,B_8tot*1e4)
    #plt.plot(z,B_9tot*1e4)
    #plt.plot(z,B_10tot*1e4)
    #plt.plot(z,B_HHtot*1e4,color="black")

    plt.ylabel("B in Gauss")
    plt.xlabel("position in m")
    plt.grid()
    plt.show()