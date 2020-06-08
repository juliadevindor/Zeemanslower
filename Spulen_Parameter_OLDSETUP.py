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


L=0.5 #m
coils=8
N_fit=np.array([647,415,386,362,336,311,268,266])
L_coil=0.036
R_coil=0.045/2
d_wire=0.001#0.005
b_wire=0.001
N_layer=L_coil/d_wire #=10 windings per layer
I_fit=np.array([4.1,4.8,4.8,4.8,4.8,4.8,4.8,4.2]) #A
Resistance_real=np.array([3.0,1.7,1.56,1.46,1.32,1.19,1.04,1.0]) #Ohm
Temp=100

#I_real=np.empty([coils]) #A
#M_real=np.empty([coils]) # number of layers
M_fit=np.empty([coils]) # number of layers
#Resistance_real=np.empty([coils]) #Ohm
Power_real=np.empty([coils]) #Watt

for i in range(coils):
    M_fit[i]=N_fit[i]/N_layer #=5 Layers
    #M_real[i]=np.abs(round(N_real[i]/N_layer,0))
    #M_real=(M_real).astype(int)
    #Resistance_real[i]=0
    #for j in range(0,int(M_fit[i])):
    #    Resistance_real[i]+=(0.017e-6*2*np.pi*N_layer*(R_coil+j*b_wire))/((d_wire/2)**2*np.pi)
    Resistance_real[i]=Resistance_real[i]*(1+4.3e-3*(Temp-20))
    Power_real[i]=I_fit[i]**2*Resistance_real[i]

print("#",end=" ")
print("N_fit",end=" ")
print("I_fit/A",end=" ")
print("Durchmesser_Spule/cm",end=" ")
print("layers_M",end=" ")
print("R_coil/Ohm",end=" ")
print("P_coil/Watt")

for i in range(coils):
    print(i,end=" ")
    print(N_fit[i],end=" ")
    print(round(I_fit[i],3),end=" ")
    print(round(100*2*(R_coil+M_fit[i]*b_wire),3),end=" ")
    print(M_fit[i],end=" ")
    print(round(Resistance_real[i],3),end=" ")
    print(round(Power_real[i],3))

#Schaltkreis durchrechnen
print(" ")

#erste Spule einzeln
print("first coil with individual power supply")
print("I_tot/A", round(I_fit[0],3), "R_tot/Ohm", round(Resistance_real[0],3),"U_tot/V",round(I_fit[0]*Resistance_real[0] ,3),"P_tot/W", round(Power_real[0],3))
print(" ")

#achte (letzte) Spule einzeln
print("eight (last) coil with individual power supply")
print("I_tot/A", round(I_fit[7],3), "R_tot/Ohm", round(Resistance_real[7],3),"U_tot/V",round(I_fit[7]*Resistance_real[7] ,3),"P_tot/W", round(Power_real[7],3))
print(" ")

#Rest des Schaltkreises gepaart parallel und dann in Reihe

num_coils=3
coils_use=np.array([[2,7],[3,6],[4,5]])
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