from light_atom_interaction import lorentzian_probability_2, lorentzian_probability
import numpy as np
import matplotlib.pyplot as plt

num=500
x=np.empty(num)
y=np.empty(num)
laserfreq=446800257000000
laserdet=18e9
for i in range(num):
    x[i]=laserfreq+laserdet
    y[i]=lorentzian_probability_2(446818663232638.56, x[i], 0, 36882297.753144175, 0.11, 0.022)
    laserdet+=1e9/(num)

plt.plot(x,y,".",label="andi")
#plt.show()

laserfreq=446799923264221.4
laserdet=-1e9
for i in range(num):
    x[i]=laserfreq+laserdet
    y[i]=lorentzian_probability(1.0,2,5,11,0.05370822448026744,2807329413780815.0, x[i], 0, 36882297.753144175, 10000, 4.9657465828447025, 1.8146475014265426, 587.7135659441951, 0, 0, -1, 6.709769684152644e-07)
    laserdet+=1e9/(num)

plt.grid()
plt.plot(x,y,".",label="meins")
plt.legend()
plt.show()
