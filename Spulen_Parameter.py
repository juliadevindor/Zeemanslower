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
coils=9
N_fit=50
L_coil=0.05
R_coil=0.043
d_wire=0.005
b_wire=0.001
N_layer=L_coil/d_wire #=10 windings per layer
I_fit=np.array([100.0,76.16286526,63.00735423,62.31878558,54.98645356,49.94757199,43.9799914,35.34623289,32.3912794]) #A
M_fit=N_fit/N_layer #=5 Layers

N_real=np.array([320,320,270,270,250,220,200,170,160]) #Vorschlag

I_real=np.empty([coils]) #A
M_real=np.empty([coils]) # number of layers
Resistance_real=np.empty([coils]) #Ohm
Power_real=np.empty([coils]) #Watt
z0=np.array([0.025,0.07900000000000001,0.131,0.183,0.235,0.28700000000000003,0.3390000000000001,0.3910000000000001,0.4450000000000002]) #m

for i in range(coils):
    M_real[i]=np.abs(round(N_real[i]/N_layer,0))
    M_real=(M_real).astype(int)
    Resistance_real[i]=0
    for j in range(M_real[i]):
        Resistance_real[i]+=(0.017e-6*2*np.pi*N_layer*(R_coil+j*b_wire))/(d_wire*b_wire)

    X=0
    z=z0[i] # pos in m
    for k in range(int(M_fit)): # number of layers used in fit
        rad=R_coil+k*b_wire
        X+=((z-z0[i]+L/2)/np.sqrt(rad**2+(z-z0[i]+L/2)**2) - (z-z0[i]-L/2)/np.sqrt(rad**2+(z-z0[i]-L/2)**2))

    Y=0
    for k in range(M_real[i]):
        rad=R_coil+k*b_wire
        Y+=((z-z0[i]+L/    2) / np.sqrt(rad ** 2 + (z - z0[i] + L / 2) ** 2) - (z - z0[i] - L / 2) / np.sqrt(
        rad ** 2 + (z - z0[i] - L / 2) ** 2))

    I_real[i]=I_fit[i]*X/Y

    Power_real[i]=I_real[i]**2*Resistance_real[i]

print("#",end=" ")
print("N_fit",end=" ")
print("I_fit/A",end=" ")
print("(N*I)_fit",end=" ")
print("N_real",end=" ")
print("Durchmesser_Spule/cm",end=" ")
print("I_real/A",end=" ")
print("layers_M",end=" ")
print("R_coil/Ohm",end=" ")
print("P_coil/Watt")

for i in range(coils):
    print(i,end=" ")
    print(N_fit,end=" ")
    print(round(I_fit[i],3),end=" ")
    print(round(N_fit * I_fit[i],3),end=" ")
    print(N_real[i],end=" ")
    print(round(100*2*(R_coil+M_real[i]*b_wire),3),end=" ")
    print(round(I_real[i],3),end=" ")
    print(M_real[i],end=" ")
    print(round(Resistance_real[i],3),end=" ")
    print(round(Power_real[i],3))
print(I_real)
#Schaltkreis durchrechnen
print(" ")

#erste Spule einzeln
print("first coil with individual power supply")
print("I_tot/A", round(I_real[0],3), "R_tot/Ohm", round(Resistance_real[0],3),"U_tot/V",round(I_real[0]*Resistance_real[0] ,3),"P_tot/W", round(Power_real[0],3))
print(" ")

#Rest des Schaltkreises gepaart parallel und dann in Reihe
print("Remaining coils in parallel in pairs of two on one power supply")
R_29=R_parallel([Resistance_real[1],Resistance_real[8]],2)
R_38=R_parallel([Resistance_real[2],Resistance_real[7]],2)
R_47=R_parallel([Resistance_real[3],Resistance_real[6]],2)
R_56=R_parallel([Resistance_real[4],Resistance_real[5]],2)
print("R_29/Ohm:",round(R_29,3),"R_38/Ohm:",round(R_38,3),"R_47/Ohm:",round(R_47,3),"R_56/Ohm:",round(R_56,3))

I_29=I_parallel([I_real[1],I_real[8]],2)
I_38=I_parallel([I_real[2],I_real[7]],2)
I_47=I_parallel([I_real[3],I_real[6]],2)
I_56=I_parallel([I_real[4],I_real[5]],2)
print("I_29/A:",round(I_29,3),"I_38/A:",round(I_38,3),"I_47/A:",round(I_47,3),"I_56/A:",round(I_56,3))

print("Then put them in a row")
R_ges=R_reihe([R_29,R_38,R_47,R_56],4)
I_ges=I_reihe([I_29,I_38,I_47,I_56],4)
print("R_ges/Ohm:",round(R_ges,3),"I_ges/A:",round(I_ges,3),"U_ges/V:",round(R_ges*I_ges,3),"P_ges/Watt:",round(R_ges*I_ges**2,3))