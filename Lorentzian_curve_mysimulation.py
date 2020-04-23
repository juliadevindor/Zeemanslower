import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as scc
from Position import Position
from trans_strength import trans_strength
from light_atom_interaction import lorentzian_probability
from scipy.optimize import curve_fit

def lorentzian(x, x0, a, gam):
    return a * gam**2 / ( gam**2 + (x - x0)**2)

laser_det=0
intensity = 1 #mW/cm^2
pol=2
GS=5
ES=11
Bfield=1e-37#10e-4 #Tesla
natural_line_width= 2 * np.pi * 5.87E6
vx=0
vy=0
vz=0 #m/s
kx=0
ky=0
kz=-1

laser_det=-1000e6
laser_freq=446799923264221.4
wavelength=scc.c / laser_freq

num=10
rho_exc=np.empty(num)
frequency=np.empty(num)

for i in range(num):
    frequency[i]=laser_freq+laser_det
    rho_exc[i]=lorentzian_probability(pol, GS, ES, Bfield, Position(GS, ES, pol, Bfield),laser_freq, laser_det, natural_line_width,intensity, vx, vy, vz, kx, ky, kz,wavelength)
    laser_det+=2*1000e6/num


plt.plot(frequency,rho_exc,".")
plt.grid()
plt.show()