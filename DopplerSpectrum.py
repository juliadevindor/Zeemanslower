#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

##########################################
## Create Li6 Spectrum Plot of hot beam ##
##########################################

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

# Drawing begins here
plt.figure(figsize=(5,4))

ax1 = plt.subplot(1,1,1)
ax1.plot( xx, V/V.max(), "k", label="Data")
ax1.plot( xx, totalSpectrum(xx, popt[0], popt[1], popt[2], 0, 1), color="orange", label="Fit")
ax1.plot( xx, dopplerSpectrum(xx, popt[0], popt[2], 0, 1), "r--",label="F=1/2")
ax1.plot( xx, dopplerSpectrum(xx, popt[1], popt[2], 228, 1), "b--",label="F=3/2")
ax1.spines["top"].set_visible(False)
ax1.set_xlabel(r"$\Delta\nu$ [-MHz]", fontsize=15)
ax1.set_ylabel("Intensity [A.U.]", fontsize=15)
ax1.set_ylim(-0.1,1.1)
ax1.legend(fontsize=10)

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

plt.savefig("UnslowedBeam.pdf",bbox_inches='tight', pad_inches=0.02, dpi=150)

#########################################################
## Create Li6 Spectrum Plot of hot beam vs slowed beam ##
#########################################################
plt.close()

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
ax1.legend(framealpha=1,fontsize=10)

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

plt.savefig("SlowedBeam.pdf",bbox_inches='tight', pad_inches=0.02, dpi=150)
