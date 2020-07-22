#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.ticker
import math
from light_atom_interaction import lorentzian_probability_2, doppler_shift
import scipy.constants as scc
##########################################
## Create Li6 Spectrum Plot of hot beam ##
##########################################

vel_light=scc.c
def Dopplershift_simple(nu, alpha, vx, vy, vz) :
    lambda1=scc.c/nu
    scalarproduct=math.cos(alpha)*vx + 0*vy - math.sin(alpha)*vz
    return nu-scalarproduct*1/lambda1

def Dopplershift(nu , alpha , vx , vy , vz ) :
   vel=math.sqrt(vx**2+vy**2+vz**2)
   beta=vel/vel_light
   gamma=1/(np.sqrt(1-beta*2))
   alpha += np.arcsin(vx/vel)
   return nu*gamma*(1+beta*np.sin(alpha))

def Probability(nu , nu_res , Gamma , I , I_sat):
    return 0.5 * (I / I_sat) * (Gamma**2) / (4 * (nu - nu_res)**2 + (Gamma ** 2) * (1 + I/I_sat))

with open("velocity_atom_start_all.txt", 'r') as f:
    lines = f.readlines()
    vx = np.asarray([float(line.split()[0]) for line in lines])
    vy = np.asarray([float(line.split()[1]) for line in lines])
    vz = np.asarray([float(line.split()[2]) for line in lines])
    gs = np.asarray([int(line.split()[3]) for line in lines])

with open("velocity_atom_all.txt", 'r') as f:
    lines = f.readlines()
    vx2 = np.asarray([float(line.split()[0]) for line in lines])
    vy2 = np.asarray([float(line.split()[1]) for line in lines])
    vz2 = np.asarray([float(line.split()[2]) for line in lines])
    gs2 = np.asarray([int(line.split()[3]) for line in lines])


# Define function which reads voltage value measured by an Oscilloscope from a file
def load(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    V = []
    for line in lines[85005:181005]:
        words = line.split("\t")
        V.append(float(words[1]))
    V = np.array(V)
    V -= V.min()
    return V

# Define function describing changed spectrum due to Doppler shift
def dopplerSpectrum(x, norm, angle, x0, freqPerTime):
    deltaFreq = (x-x0) * freqPerTime #transform x into delta freq
    freq = -deltaFreq + 446799677 #add resonance freq in MHz
    mc2 = 5603.051e6 #m_Li*c^2 in eV(?)
    kb = 8.61733e-5 #kB in eV
    temperature = 273.15+500 #500°C into K
    returnValue = norm * (deltaFreq/freq)**3/freq*np.exp( -mc2/2/kb/temperature*(deltaFreq/freq/np.sin(angle*np.pi/180))**2 )
    #this function is explained in the Mathematica notebook
    returnValue[ returnValue < 0 ] = 0
    return returnValue

def totalSpectrum(x, norm1, norm2, angle, x0, freqPerTime):
    term1 = dopplerSpectrum(x, norm1, angle, x0, freqPerTime)
    term2 = dopplerSpectrum(x, norm2, angle, x0 + 228/freqPerTime, freqPerTime)
    return term1 + term2 #sum of two peaks

# Fit function to measured spectrum here
V = load("UnslowedBeam.txt")

xx = np.linspace(0,V.size,V.size)
yy = totalSpectrum(xx, 1/2, 1, 5, 5835, 228/(55000-27000))
popt, pcov = curve_fit( totalSpectrum, xx, V/V.max(),
                        p0=( 1/2/yy.max(), 1/yy.max(), 4.5, 0, 228/(55000-27000)))
with np.printoptions(precision=3, suppress=True):
    print(pcov)

print("")
print("norm1\t=\t{0}".format(popt[0]))
print("norm2\t=\t{0}".format(popt[1]))
print("angle\t=\t{0}".format(popt[2]))
print("x0\t=\t{0}".format(popt[3]))
print("fPerT\t=\t{0}".format(popt[4]))
print("")
#with x0 and fperT the x-axis can be transformed into delta nu
xx = (xx-popt[3]) * popt[4]

#########################################################
## Create Li6 Spectrum Plot of hot beam vs slowed beam ##
#########################################################

V_slowed = load("SlowedBeam.txt")
plt.figure(figsize=(5,4))

# Find minimum residual for V_slowed. Therefore, compare end tails (right end) of V and V_slowed
V_slowed /= V.max()
V /= V.max()

width = 20000
minAt = 0
minResidual = 0
for run in range(10):
    minAt = 0
    minResidual = ((V_slowed[-width:] - V[-width:])**2).sum()
    for i in range(1,15000):
        residual = ((V_slowed[-width:] - V[-width-i:-i])**2).sum()
        if residual < minResidual:
            minResidual = residual
            minAt = i
    V_slowed *= (V_slowed[-width:] * V[-width-minAt:-minAt]).sum()/(V_slowed[-width:]**2).sum()

ax1 = plt.subplot(1,1,1)
ax1.plot( xx, V, "k", label="No Zeeman laser")
ax1.plot( xx[:-minAt], V_slowed[minAt:] , color="darkorange", label="Laser on")
ax1.spines["top"].set_visible(False)
ax1.set_xlabel(r"$\Delta\nu$ [-MHz]", fontsize=15)
ax1.set_ylabel("Intensity [A.U.]", fontsize=15)
ax1.set_ylim(-0.1,1.7)

# Add two x-axes above
ax2 = ax1.twiny()
newlabels = np.array([0,1000,2000,3000,4000,5000,6000])
newpos = 446799677 * (1-1/(1+newlabels/2.998e8*np.sin(popt[2]*np.pi/180)))
ax2.plot( [newpos[0], newpos[0]], ax1.get_ylim(), "r--", alpha=0.3 )
ax2.set_xticks(newpos)
ax2.set_xticklabels(newlabels)
ax2.xaxis.set_ticks_position("top")
ax2.tick_params(axis="x", colors="red")
ax2.xaxis.set_label_position("top")
ax2.xaxis.label.set_color("red")
ax2.spines["top"].set_color("red")
ax2.set_xlabel("Velocity of F=1/2 [m/s]", fontsize=15)
ax2.set_xlim(ax1.get_xlim())

ax3 = ax1.twiny()
newlabels = np.array([0,1000,2000,3000,4000])
newpos = 446799677-(446799677-228)/(1+newlabels/2.998e8*np.sin(popt[2]*np.pi/180))
ax3.plot( [newpos[0], newpos[0]], ax1.get_ylim(), "b--", alpha=0.3 )
ax3.set_xticks(newpos)
ax3.set_xticklabels(newlabels)
ax3.xaxis.set_ticks_position("top")
ax3.tick_params(axis="x", colors="blue")
ax3.xaxis.set_label_position("top")
ax3.xaxis.label.set_color("blue")
ax3.spines["top"].set_position(("outward",40))
ax3.spines["top"].set_color("blue")
ax3.spines["top"].set_bounds(newpos[0], ax1.get_xlim()[-1])
ax3.set_xlabel("Velocity of F=3/2 [m/s]", fontsize=15)
ax3.set_xlim(ax1.get_xlim())




num=len(xx)
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
nu0 = 446799978232118.25 #2*math.pi*446799900000000 #in Hz
Gamma = natural_line_width #in Hz
nu = np.linspace(nu0 - 60*Gamma , nu0 + 100*Gamma , len(xx))
spectrum = np.zeros(nu.size)

#ax4 = ax1.twiny()

j=4.468405424956707 #alpha in degrees

spectrum_simple = np.zeros(nu.size)
spectrum_simple_s = np.zeros(nu.size)
alpha_0=-j*math.pi/180 #degree to rad
for i in range(len(vx)):
    nu_shifted_simple=Dopplershift_simple(nu,alpha_0,vx[i],vy[i],vz[i])
    prob_simple = Probability(nu_shifted_simple, init_freq[gs[i]], Gamma, sat, sat)
    spectrum_simple += prob_simple
    nu_shifted_simple_s = Dopplershift_simple(nu, alpha_0, vx2[i], vy2[i], vz2[i])
    prob_simple_s = Probability(nu_shifted_simple_s, init_freq[gs2[i]], Gamma, sat, sat)
    spectrum_simple_s += prob_simple_s

max1=np.where(spectrum_simple==np.amax(spectrum_simple[:int(len(spectrum_simple)/2)]))
max2=np.where(spectrum_simple==np.amax(spectrum_simple[int(len(spectrum_simple)/2)+1:]))
max1_old=np.where(V==np.amax(V[:int(len(V)/2)-10000]))
max2_old=np.where(V==np.amax(V))
max1_old_ind=max1_old[0][1]
max2_old_ind=max2_old[0][1]

print("max simu_old",xx[max1],xx[max2],xx[max1]-xx[max2])
print("max data_old",xx[max1_old_ind],xx[max2_old_ind],xx[max1_old_ind]-xx[max2_old_ind])
print("max simu_nu",nu[max1]*1e-12,nu[max2]*1e-12,(nu[max1]-nu[max2])*1e-6)

xxnew=[0.0]#xx*(228.00051256/504.18579242)

for i in range(len(xx)-1):
    if i<len(xx[int(max1[0]):]):
        xxnew.append(xx[int(max1[0])+i])
    else:
        xxnew.append(xxnew[i-1]+0.01)

xxnew = xxnew[::-1]
print("NEW")
#print("max simu_new",xxnew[max1],xxnew[max2],xxnew[max1]-xxnew[max2])
print("max simu_new_NU",nu[max1]*1e-12,nu[max2]*1e-12,(nu[max1]-nu[max2])*1e-6)

ax1.plot(xxnew,spectrum_simple/spectrum_simple.max(),".",label="unslowed beam for alpha={}°".format(j))
#ax1.plot(xxnew,spectrum_simple_s/spectrum_simple.max(),label="slowed beam for alpha={}°".format(j))

#ax4.plot(xx,totalSpectrum(xx, popt[0], popt[1], popt[2], 0, 1),label="FIT")

#plt.xlabel("Laser Frequency in THz", fontsize=22)
#plt.ylabel("Sum of Lorentzian Probabilities", fontsize=22)
#plt.grid()
#plt.ylim(-10,2100)
#ax = plt.gca()
#plt.legend()
#ax.ticklabel_format(useOffset=False)
#plt.show()


ax1.legend(framealpha=1,fontsize=10)
ax1.grid()
plt.show()
plt.savefig("SlowedBeam.pdf",bbox_inches='tight', pad_inches=0.02, dpi=150)
