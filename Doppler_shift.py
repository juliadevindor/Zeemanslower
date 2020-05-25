import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import math
from light_atom_interaction import lorentzian_probability_2, doppler_shift
import scipy.constants as scc

vel_light=scc.c
def Dopplershift_simple(nu, alpha, vx, vy, vz) :
    lambda1=scc.c/nu
    scalarproduct=math.cos(alpha)*vx + 0*vy - math.sin(alpha)*vz
    return nu-scalarproduct*2*math.pi/lambda1

def Dopplershift(nu , alpha , vx , vy , vz ) :
   vel=math.sqrt(vx**2+vy**2+vz**2)
   beta=vel/vel_light
   gamma=1/(np.sqrt(1-beta*2))
   alpha += np.arcsin(vx/vel)
   return nu*gamma*(1+beta*np.sin(alpha))

def Probability(nu , nu_res , Gamma , I , I_sat):
    return 0.5 * (I / I_sat) * (Gamma**2) / (4 * (nu - nu_res)**2 + (Gamma ** 2) * (1 + I/I_sat))

with open("velocity_atom_start_SF.txt", 'r') as f:
    lines = f.readlines()
    vx = np.asarray([float(line.split()[0]) for line in lines])
    vy = np.asarray([float(line.split()[1]) for line in lines])
    vz = np.asarray([float(line.split()[2]) for line in lines])
    gs = np.asarray([int(line.split()[3]) for line in lines])

num=1000
lorentz=np.empty([len(vx)])
lorentz_sum=np.empty(num)
laser_det=np.empty(num)

kx=0
ky=0
kz=-1
wavelength=6.709770033520598e-07
#laser_frequency=2*math.pi*446799900000000
#laser_detuning=-1000e6#-920000000.0
natural_line_width= 5.87E6
sat=1 #I/I_sat
laser_intensity=sat
laser_sat_intensity=sat
init_freq= np.array([446799978232118.25 , 446799978232118.25-228e6])#excitation_frequency -freq_shift_splitting[current_groundstate]-->0 oder -228E6
alpha_0=6*math.pi/180 #degree to rad

nu0 = 446799978232118.25 #2*math.pi*446799900000000 #in Hz
Gamma = natural_line_width #in Hz
nu = np.linspace(nu0 - 350*Gamma , nu0 + 50*Gamma , 50000)
spectrum = np.zeros(nu.size)

for i in range(len(vx)):
    #nu_shifted = Dopplershift(nu, alpha_0, vx[i], vy[i], vz[i])
    nu_shifted=Dopplershift_simple(nu,alpha_0,vx[i],vy[i],vz[i])
    prob = Probability(nu_shifted, init_freq[gs[i]], Gamma , sat, sat)
    spectrum += prob

print(np.sum(spectrum))
fig=plt.figure()
ax = fig.add_subplot(111)
ax.plot(nu*1e-12,spectrum,".")
plt.xlabel("Laser Frequency in THz", fontsize=22)
plt.ylabel("Sum of Lorentzian Probabilities", fontsize=22)
plt.grid()
#plt.ylim(-10,2100)
ax = plt.gca()
ax.ticklabel_format(useOffset=False)

plt.show()
