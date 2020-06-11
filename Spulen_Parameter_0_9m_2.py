import math
import numpy as np

def I_parallel(I_coils,num):
    I_tot=0
    for i in range(num):
        I_tot+=I_coils[i]
    return I_tot

def R_parallel(R_coils,num):
    R_tot=0
    for i in range(num):
        R_tot+=1/R_coils[i]
    return 1/R_tot

def U_parallel(U_coils,num):
    return U_coils[0]

def I_reihe(I_coils,num):
    return I_coils[0]

def R_reihe(R_coils,num):
    R_tot=0
    for i in range(num):
        R_tot+=R_coils[i]
    return R_tot

def U_reihe(U_coils,num):
    U_tot=0
    for i in range(num):
        U_tot+=U_coils[i]
    return U_tot


L=0.9 #m
coils=17
L_coil=0.05
R_coil=0.045/2
d_wire=0.005
b_wire=0.001
N_layer=L_coil/d_wire #=10 windings per layer
I_fit=np.array([25.0,25.0,23.26,24.916,25.0,24.205,23.232,22.975,22.836,24.504,25.0,24.993,24.109,24.394,22.694,18.602,17.01]) #A
Temp=100
N_real=np.array([320,300,300,260,250,250,250,240,230,200,180,170,160,140,130,130,100]) #Vorschlag
M_real=np.empty([coils]) # number of layers
Resistance_real=np.empty([coils]) #Ohm
Power_real=np.empty([coils]) #Watt
#z0=np.array([0.025,0.07900000000000001,0.131,0.183,0.235,0.28700000000000003,0.3390000000000001,0.3910000000000001,0.4450000000000002]) #m

for i in range(coils):
    M_real[i]=np.abs(round(N_real[i]/N_layer,0))
    M_real=(M_real).astype(int)
    Resistance_real[i]=0
    for j in range(M_real[i]):
        Resistance_real[i]+=(0.017e-6*2*np.pi*N_layer*(R_coil+j*b_wire))/(d_wire*b_wire)
    Resistance_real[i]=Resistance_real[i]*(1+4.3e-3*(Temp-20))

    #X=0
    #z=z0[i] # pos in m
    #for k in range(int(M_fit)): # number of layers used in fit
    #    rad=R_coil+k*b_wire
    #    X+=((z-z0[i]+L/2)/np.sqrt(rad**2+(z-z0[i]+L/2)**2) - (z-z0[i]-L/2)/np.sqrt(rad**2+(z-z0[i]-L/2)**2))

    #Y=0
    #for k in range(M_real[i]):
    #    rad=R_coil+k*b_wire
    #    Y+=((z-z0[i]+L/    2) / np.sqrt(rad ** 2 + (z - z0[i] + L / 2) ** 2) - (z - z0[i] - L / 2) / np.sqrt(
    #    rad ** 2 + (z - z0[i] - L / 2) ** 2))

    #I_fit[i]=I_fit[i]*X/Y
    #I_fit[i]=N_fit*I_fit[i]/N_real[i]

    Power_real[i]=I_fit[i]**2*Resistance_real[i]

print("#",end=" ")
print("I_fit/A",end=" ")
print("N_real",end=" ")
print("Durchmesser_Spule/cm",end=" ")
print("layers_M",end=" ")
print("R_coil/Ohm",end=" ")
print("P_coil/Watt")

for i in range(coils):
    print(i,end=" ")
    print(round(I_fit[i],3),end=" ")
    print(N_real[i],end=" ")
    print(round(100*2*(R_coil+M_real[i]*b_wire),3),end=" ")
    print(M_real[i],end=" ")
    print(round(Resistance_real[i],3),end=" ")
    print(round(Power_real[i],3))

print(" ")

for i in range(coils):
    print(round(I_fit[i],3),end=",")

#Schaltkreis durchrechnen
print(" ")

#erste Spule einzeln
print("first coil with individual power supply")
print("I_tot/A", round(I_fit[0],3), "R_tot/Ohm", round(Resistance_real[0],3),"U_tot/V",round(I_fit[0]*Resistance_real[0] ,3),"P_tot/W", round(Power_real[0],3))
print(" ")

#Rest des Schaltkreises gepaart parallel und dann in Reihe

num_coils=3
coils_use=np.array([[2,17],[3,16],[4,15]])
print("{} following coils in parallel in pairs of two on one power supply".format(num_coils*2))
R=np.empty([num_coils])
I=np.empty([num_coils])
for i in range(num_coils):
    R[i]=R_parallel([Resistance_real[coils_use[i,0]-1],Resistance_real[coils_use[i,1]-1]],2)
    print("R_{}_{}=".format(coils_use[i,0],coils_use[i,1]),round(R[i],3))
    I[i]=I_parallel([I_fit[coils_use[i,0]-1],I_fit[coils_use[i,1]-1]],2)
    print("I_{}_{}=".format(coils_use[i,0],coils_use[i,1]),round(I[i],3))
R_ges=R_reihe(R,num_coils)
I_ges=I_reihe(I,num_coils)
print("R_ges/Ohm:",round(R_ges,3),"I_ges/A:",round(I_ges,3),"U_ges/V:",round(R_ges*I_ges,3),"P_ges/Watt:",round(R_ges*I_ges**2,3))
print(" ")

num_coils=3
coils_use=np.array([[5,14],[6,13],[7,12]])
print("{} following coils in parallel in pairs of two on one power supply".format(num_coils*2))
R=np.empty([num_coils])
I=np.empty([num_coils])
for i in range(num_coils):
    R[i]=R_parallel([Resistance_real[coils_use[i,0]-1],Resistance_real[coils_use[i,1]-1]],2)
    print("R_{}_{}=".format(coils_use[i,0],coils_use[i,1]),round(R[i],3))
    I[i]=I_parallel([I_fit[coils_use[i,0]-1],I_fit[coils_use[i,1]-1]],2)
    print("I_{}_{}=".format(coils_use[i,0],coils_use[i,1]),round(I[i],3))
R_ges=R_reihe(R,num_coils)
I_ges=I_reihe(I,num_coils)
print("R_ges/Ohm:",round(R_ges,3),"I_ges/A:",round(I_ges,3),"U_ges/V:",round(R_ges*I_ges,3),"P_ges/Watt:",round(R_ges*I_ges**2,3))
print(" ")

num_coils=2
coils_use=np.array([[8,11],[9,10]])
print("{} following coils in parallel in pairs of two on one power supply".format(num_coils*2))
R=np.empty([num_coils])
I=np.empty([num_coils])
for i in range(num_coils):
    R[i]=R_parallel([Resistance_real[coils_use[i,0]-1],Resistance_real[coils_use[i,1]-1]],2)
    print("R_{}_{}=".format(coils_use[i,0],coils_use[i,1]),round(R[i],3))
    I[i]=I_parallel([I_fit[coils_use[i,0]-1],I_fit[coils_use[i,1]-1]],2)
    print("I_{}_{}=".format(coils_use[i,0],coils_use[i,1]),round(I[i],3))
R_ges=R_reihe(R,num_coils)
I_ges=I_reihe(I,num_coils)
print("R_ges/Ohm:",round(R_ges,3),"I_ges/A:",round(I_ges,3),"U_ges/V:",round(R_ges*I_ges,3),"P_ges/Watt:",round(R_ges*I_ges**2,3))