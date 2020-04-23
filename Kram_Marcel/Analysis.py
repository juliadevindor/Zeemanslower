#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit

###############################################################################
# Define loading of file.
# The oscilloscope measurements contain a preamble of 5 lines. These are skipped.
# Each line contains the time of the measurement and the measured voltage.
# Load both informations and return a numpy array for each.
###############################################################################

def load(fileName):
    file = open(fileName)
    lines = file.readlines()
    file.close()
    t, V = [], []
    for line in lines[5:]:
        words = line.split('\t')
        t.append( float( words[0] ) )
        V.append( float( words[1] ) )
    return np.array(t), np.array(V)

###############################################################################
# Define Doppler shifted spectrum.
# The following function was gained by convoluting a Dirac delta distribution
# and a modified (v**3) Doppler distribution using Mathematica.
# For further information contact the author.
###############################################################################

def spectrum(freqL, freq0, T, angle, norm):
    omegaL = freqL * const.pi * 2
    omega0 = freq0 * const.pi * 2
    c = const.speed_of_light
    kB = const.Boltzmann
    m = 9.988346e-27
    partA = 2*c*omegaL**2*np.cos(angle)
    partB = np.sqrt( 2*c**2*omega0**2*( 2*omega0**2-omegaL**2+omegaL**2*np.cos(2*angle) ) )
    partC = omega0**2+omegaL**2*np.cos(angle)**2
    condition = (partA + partB > 0)
    returnValue = condition * norm * np.exp( -m*(partA+partB)**2/(8*kB*T*partC**2) ) * (partA + partB)**2/(2*np.sqrt(2*const.pi)*(kB*T)**(3/2)*partC**2)
    return returnValue

###############################################################################
# Define Doppler shifted spectrum of two near transitions.
# The transitions are shifted by 228 MHz.
# Each spectrum can be renormalized independently.
###############################################################################

def doubleSpectrum(freqL, freq0, T, angle, norm_1, norm_2):
    returnValue  = spectrum(freqL, freq0      , T, angle, norm_1)
    returnValue += spectrum(freqL, freq0+228e6, T, angle, norm_2)
    return returnValue/np.max(returnValue)



if __name__ == '__main__':

    ###########################################################################
    # Load different measurements.
    # The numbers (00000 and above) represent the measurement.
    # t_* is always the respective time, V_* is always the voltage.
    # M1 arrays should always be the same. They contain the memorized
    # measurement of an unslowed atom beam.
    # C2 contains the ramp voltage of the measurement. A falling/positive
    # value stands for falling frequencies. A rising/negative value stands
    # for rising frequencies.
    # C3 contains a real measurement.
    ###########################################################################

    t_M1_00000, V_M1_00000 = load('21092018-Zeeman/M1--Slower210918--00000.txt')
    V_M1_00000 -= np.min(V_M1_00000)
    V_M1_00000 /= np.max(V_M1_00000)

    t_C2_00000, V_C2_00000 = load('21092018-Zeeman/C2--Slower210918--00000.txt')
    V_C2_00000 /= np.max(V_C2_00000)

    t_C3_00000, V_C3_00000 = load('21092018-Zeeman/C3--Slower210918--00000.txt')
    V_C3_00000 -= np.min(V_C3_00000)
    V_C3_00000 /= np.max(V_C3_00000)

    ###########################################################################
    # To get better spectra, only measurements with negative V_C2 are used.
    # This means, that we only look at a continously rising frequency.
    # The ends are cut of to guarantee this.
    ###########################################################################

    condition  = V_C2_00000 < 0
    t_M1_00000 = t_M1_00000[condition][30:-30]
    V_M1_00000 = V_M1_00000[condition][30:-30]
    t_C2_00000 = t_C2_00000[condition][30:-30]
    V_C2_00000 = V_C2_00000[condition][30:-30]
    t_C3_00000 = t_C3_00000[condition][30:-30]
    V_C3_00000 = V_C3_00000[condition][30:-30]

    ###########################################################################
    # Draw the spectrum.
    # Everything below this is TO DO.
    # The time of the measurement has to be adjusted so that measurment
    # and function are drawn on top of each other.
    ###########################################################################

    wavelength = 670.977e-9
    frequency = const.speed_of_light/wavelength
    f = np.linspace(frequency-1e9, frequency+0.25e9,1000)
    p = doubleSpectrum(f, frequency, 700, 94/180*const.pi, 3, 1)
    plt.plot(t_M1_00000, V_M1_00000)
    plt.plot(t_C2_00000, V_C2_00000)
    plt.plot(t_C3_00000, V_C3_00000)
    plt.show()
