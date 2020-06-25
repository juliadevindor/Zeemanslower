import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import os
import glob
import scipy.constants as scc


hPlanck=scc.h #6.62607004e-34 Js
muB=scc.physical_constants['Bohr magneton'][0] # 9.274009994e-24 J/T

def slower_field(pos, L0,v0,wavelength,omega, omega0):
  # length of slower L0 in m
    B0 = hPlanck*v0/(wavelength*muB) #0.08 #0.19185 # y-Achsenabschnitt
    #print(B0)
    Bbias = 0#-3e-2 #0#(hPlanck/(2*np.pi))*(omega-omega0)/(muB)#-1.478*1e-9  # shiftet nur den Plot nach oben/ unten
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

L0=1.0#m
v0_0=1200 #m/s
v0_1=800 #m/s
v0_2=500 #m/s
wavelength=670.992e-9 #m (Gehm)
omega=446799923264221.4 #Hz (mein Wert) #446799900000000 #Hz (Andis Wert)
omega0=446.7896e12 #Hz (Gehm)

pos=np.linspace(-0.0,L0+0.2,num=1000)
B = np.empty(len(pos))
B=slower_field(pos,L0,1580,wavelength,omega,omega0)
file=open("B(z).txt","w+")

for i in range(len(pos)):
    if B[i]!=0:
        file.write("{};{}\n".format(pos[i]-0.5,B[i]*1e4))
file.close()

fig, ax = plt.subplots()

with open("B(z)_0_5m_changed.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="blue",label="decr. field slower of length 0.5m")
    print("L=",0.5,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.5, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.5,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_0_6m_changed.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="orange",label="decr. field slower of length 0.6m")
    print("L=",0.6,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.6, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.6,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_0_7m_changed.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="green",label="decr. field slower of length 0.7m")
    print("L=",0.7,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.7, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.7,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_0_8m_changed.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="red",label="decr. field slower of length 0.8m")
    print("L=",0.8,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.8, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.8,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_0_9m_changed.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="purple",label="decr. field slower of length 0.9m")
    print("L=",0.9,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.9, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.9,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_1_0m_changed.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="brown",label="decr. field slower of length 1m")
    print("L=",1,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 1, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",1,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))

with open("B(z)_fit_0_5m.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="blue")
    print("FIT: L=",0.5,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.5, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.5,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_fit_0_6m.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="orange")
    print("L=",0.6,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.6, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.6,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_fit_0_7m.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="green")
    print("L=",0.7,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.7, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.7,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_fit_0_8m.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="red")
    print("L=",0.8,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.8, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.8,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_fit_0_9m.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="purple")
    print("L=",0.9,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 0.9, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",0.9,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))
with open("B(z)_fit_1_0m.txt","r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5,y,color="brown")
    print("L=",1,"m: ",round(max(y),3),round(x[np.argmax(y)]+0.5,3))
    print("L=", 1, "m: ", round(min(y), 3), round(x[np.argmin(y)] + 0.5, 3))
    print("L=",1,"m: ",round(x[np.argmin(y)] + 0.5, 3)-round(x[np.argmax(y)]+0.5,3))

with open("B(z).txt","r") as f: # plot measured magnetic field
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5,y, label="new", color="black")
mass_li=9.988142325599999e-27
gamma=36882297.753144175
lambda_=6.70977e-07
k_=2*np.pi/lambda_
wavelength_laser=6.709770033520598e-07
amax=hPlanck/(2*np.pi)*gamma*k_/(2*mass_li)*0.5
#print("vmax",np.sqrt(2*amax*0.5))
#print("vmax",np.sqrt(2*amax*0.6))
#print("vmax",np.sqrt(2*amax*0.7))
#print("vmax",np.sqrt(2*amax*0.8))
#print("vmax",np.sqrt(2*amax*0.9))
#print("vmax",np.sqrt(2*amax*1.0))

#ax.plot(pos+0.1,B*1e4,label="v0=1000m/s")
plt.annotate("z in m", xy=(1.01, 0), ha='left', va='top', xycoords='axes fraction', fontsize=22)
plt.annotate("B in Gauss", xy=(-0.05, 1.10), ha='left', va='top', xycoords='axes fraction', fontsize=22)
#plt.hlines(0,pos[0]+0.1,pos[-1]+0.1)
ax.spines['left'].set_position('zero')
# set the y-spine
ax.spines['bottom'].set_position('zero')
#plt.xlabel("Position / m")
plt.rcParams.update({'font.size': 15})
#plt.ylabel("B / T")
xticks = ax.xaxis.get_major_ticks()
xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.legend()
plt.grid()
plt.show()


