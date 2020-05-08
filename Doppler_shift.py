import matplotlib.pyplot as plt
import numpy as np
import math
from light_atom_interaction import lorentzian_probability_2, doppler_shift

with open("velocity_atom.txt", 'r') as f:
    lines = f.readlines()
    vx = np.asarray([float(line.split()[0]) for line in lines])
    vy = np.asarray([float(line.split()[1]) for line in lines])
    vz = np.asarray([float(line.split()[2]) for line in lines])

num=100
lorentz=np.empty([len(vx)])
lorentz_sum=np.empty(num)
laser_det=np.empty(num)

kx=0
ky=0
kz=-1
wavelength=6.709770033520598e-07
laser_frequency=2*math.pi*446799900000000
laser_detuning=-1e9#-920000000.0
natural_line_width=2 * math.pi * 5.87E6
sat=1 #I/I_sat
init_freq= 446799978232118.25-0#excitation_frequency -freq_shift_splitting[current_groundstate]-->0 oder -228E6

for det in range(num):
    print(det)
    laser_det[det]=laser_detuning
    for i in range(len(vx)):
        atom_freq=doppler_shift(vx[i], vy[i], vz[i], init_freq, kx, ky, kz, wavelength)
        lorentz[i]=lorentzian_probability_2(atom_freq, laser_frequency, laser_detuning, natural_line_width, sat,sat)
        lorentz_sum[det]+=lorentz[i]
    laser_detuning+=2e9/num

plt.plot(laser_det,lorentz_sum,".")
plt.show()