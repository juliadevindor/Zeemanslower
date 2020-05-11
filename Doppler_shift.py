import matplotlib.pyplot as plt
import numpy as np
import math
from light_atom_interaction import lorentzian_probability_2, doppler_shift

with open("velocity_atom.txt", 'r') as f:
    lines = f.readlines()
    vx = np.asarray([float(line.split()[0]) for line in lines])
    vy = np.asarray([float(line.split()[1]) for line in lines])
    vz = np.asarray([float(line.split()[2]) for line in lines])

num=10000
lorentz=np.empty([len(vx)])
lorentz_sum=np.empty(num)
laser_det=np.empty(num)

kx=0
ky=0
kz=-1
wavelength=6.709770033520598e-07
laser_frequency=2*math.pi*446799900000000
laser_detuning=-500e6#-920000000.0
natural_line_width=2 * math.pi * 5.87E6
sat=1 #I/I_sat
laser_intensity=sat
laser_sat_intensity=sat
init_freq= 2*math.pi*(446799978232118.25-0)#excitation_frequency -freq_shift_splitting[current_groundstate]-->0 oder -228E6

for det in range(num):
    print(det)
    laser_det[det]=laser_detuning
    lorentz_sum[det]=0
    for i in range(len(vx)):
        atom_freq= init_freq - (math.sqrt(vx[i]**2+vy[i]**2+vz[i]**2)*math.sqrt(kx**2+ky**2+kz**2)*math.cos(0)) * 2 * (math.pi / wavelength) #cos in rad
        laser_freq_modus = laser_frequency - 2 * math.pi * laser_detuning  # ohne 2 pi und Vorzeichen geändert
        excitation_probability = 0.5 * (laser_intensity / laser_sat_intensity) * (natural_line_width ** 2) / (
                    4 * (laser_freq_modus - atom_freq) ** 2 + (natural_line_width ** 2) * (
                        1 + laser_intensity / laser_sat_intensity))

        lorentz[i]=excitation_probability
        lorentz_sum[det]+=lorentz[i]
    laser_detuning+=3000e6/num

plt.plot(laser_det*1e-6,lorentz_sum)
plt.xlabel("Laser Detuning in MHz", fontsize=15)
plt.ylabel("Sum of Lorentzian Probabilities", fontsize=15)
plt.grid()

plt.show()