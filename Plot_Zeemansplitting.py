import matplotlib.pyplot as plt
import matplotlib.pyplot
import numpy as np
from cmath import sqrt
from numba import vectorize

import cmath

@vectorize(["complex128(complex128)"])

def SQRT(x):
    return cmath.sqrt(x)



matplotlib.rcParams.update({'font.size': 15})

fig, ((ax1), (ax2)) = plt.subplots(2, 1, sharex=False, sharey=False)
plt.subplots_adjust(hspace=0.5)
ax1.set(xlabel="Bfield/ G", ylabel="Energy shfit/ MHz")
ax2.set_title("Groundstate 2S1/2")
B=np.linspace(0,160,1000)
p1=[76.0684,1.40098,0.982695,-38.704,0.000318876,13482.,110.431,2.03522]

ax2.plot(B,p1[0] +p1[1]*B, label="F=3/2, mf=3/2")
ax2.plot(B,p1[2]*(p1[3]-p1[4]*B+np.sqrt(p1[5]+p1[6]*B+p1[7]*B**2)), label="F=3/2, mf=1/2")
ax2.plot(B,p1[2]*(p1[3]+p1[4]*B+np.sqrt(p1[5]-p1[6]*B+p1[7]*B**2)), label="F=3/2, mf=-1/2")
ax2.plot(B,p1[0] -p1[1]*B, label="F=3/2, mf=-3/2")
ax2.plot(B,p1[2]*(p1[3]+p1[4]*B-np.sqrt(p1[5]-p1[6]*B+p1[7]*B**2)), label="F=1/2, mf=-1/2")
ax2.plot(B,p1[2]*(p1[3]-p1[4]*B-np.sqrt(p1[5]+p1[6]*B+p1[7]*B**2)), label="F=1/2, mf=1/2")

ax2.set(xlabel="Bfield/ G",ylabel="Energy shfit/ MHz")
ax1.set_title("Excited state 2P3/2")
B=np.linspace(0,6,1000)

p2=[-1.65,2.80287,1.30906,-0.210075,1.42751,1.10329,\
    0.300034,0.509956,-0.81682,-0.897791,1.14407,\
    0.514564,0.891251j,-2.46846,1.02748,1.74637,\
    -3.23379,0.691848,4.70363,1.77636e-15,-49.707,\
    79.6037,188.909,104.134,90.3401,37.6035,\
    21.3044,1.02913,0.64831,8.88178e-16,0.324155,\
    0.561453j]

ax1.plot(B,( p2[8]*(p2[9]-p2[10]*B)-(p2[27]*(p2[13]-p2[14]*B-p2[15]*B**2))/(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)+p2[28]*(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="F=1/2, mf=1/2" )
ax1.plot(B,( p2[8]*(p2[9]+p2[10]*B)-(p2[27]*(p2[13]+p2[14]*B-p2[15]*B**2))/(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)+p2[28]*(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="F=1/2, mf=-1/2" )
ax1.plot(B,p2[2]*(p2[3]+p2[4]*B+np.sqrt(p2[5] +p2[6]*B+p2[7]*B**2)), label="F=3/2 mf=3/2")
ax1.plot(B,( p2[8]*(p2[9]-p2[10]*B)+((p2[11] -p2[12])*(p2[13]-p2[14] *B-p2[15] *B**2))/(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)-(p2[30] +p2[31])*(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3) ).real, label="F=3/2 mf=1/2")
ax1.plot(B,( p2[8]*(p2[9]+p2[10]*B)+((p2[11] -p2[12])*(p2[13]+p2[14] *B-p2[15] *B**2))/(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)-(p2[30] +p2[31])*(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="F=3/2 mf=-1/2" )
ax1.plot(B,p2[2]*(p2[3]-p2[4]*B+np.sqrt(p2[5] -p2[6]*B+p2[7]*B**2)), label="F=3/2 mf=-3/2")
ax1.plot(B,p2[0]+p2[1]*B, color="yellow",label="F=5/2 mf=5/2")
ax1.plot(B,p2[2]*(p2[3]+p2[4]*B-np.sqrt(p2[5] +p2[6]*B+p2[7]*B**2)), label="F=5/2 mf=3/2")
ax1.plot(B,( p2[8]*(p2[9]-p2[10]*B)+((p2[11] +p2[12])*(p2[13]-p2[14] *B-p2[15] *B**2))/(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)-(p2[30] -p2[31])*(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="F=5/2 mf=1/2" )
ax1.plot(B,( p2[8]*(p2[9]+p2[10]*B)+((p2[11] +p2[12])*(p2[13]+p2[14] *B-p2[15] *B**2))/(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)-(p2[30] -p2[31])*(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="F=5/2 mf=-1/2" )
ax1.plot(B,p2[2]*(p2[3]-p2[4]*B-np.sqrt(p2[5]-p2[6]*B+p2[7]*B**2)), label="F=5/2 mf=-3/2")
ax1.plot(B,p2[0]-p2[1]*B, color="black", label="F=5/2 mf=-5/2")

#ax1.grid()
#ax2.grid()
ax1.legend()
ax2.legend()
plt.show()

B=np.linspace(0,600,10000)

plt.plot(B,p1[0] +p1[1]*B, label="F=3/2, mf=3/2",color="red")
plt.plot(B,p1[2]*(p1[3]-p1[4]*B+np.sqrt(p1[5]+p1[6]*B+p1[7]*B**2)), label="F=3/2, mf=1/2",color="red")
plt.plot(B,p1[2]*(p1[3]+p1[4]*B+np.sqrt(p1[5]-p1[6]*B+p1[7]*B**2)), label="F=3/2, mf=-1/2",color="red")
plt.plot(B,p1[0] -p1[1]*B, label="F=3/2, mf=-3/2",color="red")
plt.plot(B,p1[2]*(p1[3]+p1[4]*B-np.sqrt(p1[5]-p1[6]*B+p1[7]*B**2)), label="F=1/2, mf=-1/2",color="red")
plt.plot(B,p1[2]*(p1[3]-p1[4]*B-np.sqrt(p1[5]+p1[6]*B+p1[7]*B**2)), label="F=1/2, mf=1/2",color="red")


plt.plot(B,( p2[8]*(p2[9]-p2[10]*B)-(p2[27]*(p2[13]-p2[14]*B-p2[15]*B**2))/(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)+p2[28]*(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="1",color="black")
plt.plot(B,( p2[8]*(p2[9]+p2[10]*B)-(p2[27]*(p2[13]+p2[14]*B-p2[15]*B**2))/(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)+p2[28]*(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="2",color="black")
plt.plot(B,p2[2]*(p2[3]+p2[4]*B+np.sqrt(p2[5] +p2[6]*B+p2[7]*B**2)), label="3",color="black")
plt.plot(B,( p2[8]*(p2[9]-p2[10]*B)+((p2[11] -p2[12])*(p2[13]-p2[14] *B-p2[15] *B**2))/(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)-(p2[30] +p2[31])*(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3) ).real, label="4",color="black")
plt.plot(B,( p2[8]*(p2[9]+p2[10]*B)+((p2[11] -p2[12])*(p2[13]+p2[14] *B-p2[15] *B**2))/(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)-(p2[30] +p2[31])*(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="5",color="black")
plt.plot(B,p2[2]*(p2[3]-p2[4]*B+np.sqrt(p2[5] -p2[6]*B+p2[7]*B**2)), label="6",color="black")
plt.plot(B,p2[0]+p2[1]*B,label="7",color="black")
plt.plot(B,p2[2]*(p2[3]+p2[4]*B-np.sqrt(p2[5] +p2[6]*B+p2[7]*B**2)), label="8",color="black")
plt.plot(B,( p2[8]*(p2[9]-p2[10]*B)+((p2[11] +p2[12])*(p2[13]-p2[14] *B-p2[15] *B**2))/(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)-(p2[30] -p2[31])*(p2[16]+p2[17] *B+p2[18] *B**2-p2[19] *B**3+SQRT(p2[20]-p2[21] *B-p2[22] *B**2-p2[23] *B**3-p2[24] *B**4-p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="9",color="black")
plt.plot(B,( p2[8]*(p2[9]+p2[10]*B)+((p2[11] +p2[12])*(p2[13]+p2[14] *B-p2[15] *B**2))/(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)-(p2[30] -p2[31])*(p2[16]-p2[17] *B+p2[18] *B**2+p2[29] *B**3+SQRT(p2[20]+p2[21] *B-p2[22] *B**2+p2[23] *B**3-p2[24] *B**4+p2[25] *B**5-p2[26] *B**6))**(1/3)).real, label="10",color="black")
plt.plot(B,p2[2]*(p2[3]-p2[4]*B-np.sqrt(p2[5]-p2[6]*B+p2[7]*B**2)), label="11",color="black")
plt.plot(B,p2[0]-p2[1]*B, color="black", label="12")
plt.xlabel("B in Gauss")
plt.grid()
plt.ylabel("Energyshift in MHz")
plt.show()