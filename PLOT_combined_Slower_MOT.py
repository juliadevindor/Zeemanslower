import matplotlib.pyplot as plt
import numpy as np

mu_0=4*np.pi*1e-7 # magnetic field constant

def B_coil(I_coil, N_coil, L_coil, R_coil, z0_coil, pos): # magnetic field of a single coil
    return 1e4*mu_0*N_coil*I_coil/(2*L_coil)*((pos-z0_coil+L_coil/2)/np.sqrt(R_coil**2+(pos-z0_coil+L_coil/2)**2) - (pos-z0_coil-L_coil/2)/np.sqrt(R_coil**2+(pos-z0_coil-L_coil/2)**2)) # magnetic field of a single coil


I_1=10
N_1=90
L_1=0.005
R_1=0.043
z0_1=0.383+0.02-L_1/2#+0.035/2
print("AT",z0_1)

I_2=-10
N_2=90
L_2=0.005
R_2=0.043
z0_2=0.383+0.015-L_2/2

shift=0.5

with open("Measurements/combined_data.txt", "r") as g:
    lines = g.readlines()
    xMOT = np.asarray([float(line.split(";")[0]) for line in lines])
    yMOT = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(xMOT+shift, yMOT, label="combined 1")
    #plt.plot(xMOT+shift, yMOT+B_coil(I_1,N_1,L_1,R_1,z0_1,xMOT+0.5), label="combined 1 with add. coils") #+B_coil(I_2,N_2,L_2,R_2,z0_2,xMOT+0.5)
    plt.vlines(0.383,-200,800)# MOT Rand
    plt.vlines(0.383+0.0355,-200,800) #Slower Ende?
    plt.vlines(0.383+0.02-L_1/2,-200,800,color="blue")# spule 1
    plt.vlines(0.383+0.015-L_2/2,-200,800,color="red")# spule 2

plt.plot(xMOT+shift, B_coil(I_1,N_1,L_1,R_1,z0_1,xMOT+0.5)+B_coil(I_2,N_2,L_2,R_2,z0_2,xMOT+0.5))
#plt.plot(xMOT+shift, B_coil(I_2,N_2,L_2,R_2,z0_2,xMOT+0.5))

file = open("field+coil1.txt", "w+")
for i in range(len(xMOT)):
    file.write("{};{}\n".format(xMOT[i],yMOT[i]+B_coil(I_1,N_1,L_1,R_1,z0_1,xMOT[i]+0.5)+B_coil(I_2,N_2,L_2,R_2,z0_2,xMOT[i]+0.5)))
file.close()

with open("./field+coil1.txt", "r") as g:
    lines = g.readlines()
    xcomb = np.asarray([float(line.split(";")[0]) for line in lines])
    ycomb = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(xcomb+0.5, ycomb, label="combined with coil")

with open("Measurements/combined_data_2.txt", "r") as g:
    lines = g.readlines()
    xMOT2 = np.asarray([float(line.split(";")[0]) for line in lines])
    yMOT2 = np.asarray([float(line.split(";")[1]) for line in lines])
    #plt.plot(xMOT2+0.5, yMOT2, label="combined 2")

with open("Measurements/example_magnetic_field.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x+shift, y, label="good field")

with open("Measurements/Zeemanslower_old_standalone_Bz.txt", "r") as g:
    lines = g.readlines()
    x4 = np.asarray([float(line.split()[0]) for line in lines])
    y4 = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot((x4-33)/100, y4*10, ".",label="Messung Zeemanslower von uns")


with open("Measurements/Zeemanslower_old_standalone_Bz_vactube_new.txt", "r") as g:
    lines = g.readlines()
    x4 = np.asarray([float(line.split()[0]) for line in lines])
    y4 = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot((x4-15)/100, y4*10, ".",label="Messung Zeemanslower von uns")

plt.grid()
plt.legend()
plt.show()