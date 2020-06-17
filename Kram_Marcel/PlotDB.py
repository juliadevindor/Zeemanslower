#!/usr/bin/python3

# Load the necessary modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define function of a peak
def peak( nu , N , nu0 , Gamma ):
    return N/( ((nu-nu0)/Gamma)**2 + 1 )

# Define function describing the background
def background( nu , offset , N , nu0 , sigma ):
    return offset - N*np.exp(-(nu-nu0)**2/2/sigma**2)

# Define function describing six peaks plus background
def func( nu , offset ,
          N_1 , nu0_1 , Gamma_1 ,
          N_2 , nu0_2 , Gamma_2 ,
          N_3 , nu0_3 , Gamma_3 ,
          N_4 , nu0_4 , Gamma_4 ,
          N_5 , nu0_5 , Gamma_5 ,
          N_6 , nu0_6 , Gamma_6 ,
          N_bkg , nu0_bkg , Gamma_bkg ):
    returnValue  = 0
    returnValue += peak( nu , N_1 , nu0_1 , Gamma_1 )
    returnValue += peak( nu , N_2 , nu0_2 , Gamma_2 )
    returnValue += peak( nu , N_3 , nu0_3 , Gamma_3 )
    returnValue += peak( nu , N_4 , nu0_4 , Gamma_4 )
    returnValue += peak( nu , N_5 , nu0_5 , Gamma_5 )
    returnValue += peak( nu , N_6 , nu0_6 , Gamma_6 )
    returnValue += background( nu , offset , N_bkg , nu0_bkg , Gamma_bkg )
    return returnValue

# Load the Excel file and read the column titles
df = pd.read_excel( "Daten Dopplerbreite.xlsx" )
cols = df.columns

# Load specific columns of the Excel file
# This has to be done in a funny way by identifying
# a column by its name.
# Additionaly the first and last 10 percent are cut
# away to only measure the dip
freq  = df[cols[0]].values
Vspec = df[cols[1]].values
Vhot  = df[cols[2]].values
freq_short  = freq[ 10*freq.size//100 : 90*freq.size//100 ]
Vspec_short = Vspec[ 10*Vspec.size//100 : 90*Vspec.size//100 ]
Vhot_short  = Vhot[ 10*Vhot.size//100 : 90*Vhot.size//100 ]

# N_1: height of peak
# nu0: position of peak
# Gamma: width of peak
# Fit a function with multiple peaks to Vspec_short
offset    = 1
N_1       = 0.02
nu0_1     = 384.22906
Gamma_1   = 0.00001
N_2       = 0.09
nu0_2     = 384.22909
Gamma_2   = 0.00001
N_3       = 0.08
nu0_3     = 384.229115
Gamma_3   = 0.00001
N_4       = 0.25
nu0_4     = 384.229145
Gamma_4   = 0.000008
N_5       = 0.37
nu0_5     = 384.229165
Gamma_5   = 0.00001
N_6       = 0.17
nu0_6     = 384.229215
Gamma_6   = 0.00001
N_bkg     = 1.07
nu0_bkg   = 384.22919
Gamma_bkg = 0.0002
params = ( offset ,
           N_1 , nu0_1 , Gamma_1 ,
           N_2 , nu0_2 , Gamma_2 ,
           N_3 , nu0_3 , Gamma_3 ,
           N_4 , nu0_4 , Gamma_4 ,
           N_5 , nu0_5 , Gamma_5 ,
           N_6 , nu0_6 , Gamma_6 ,
           N_bkg , nu0_bkg , Gamma_bkg )
popt, pcov = curve_fit( func , freq_short , Vspec_short ,
                        p0 = params)
print( "Offset:\t" , popt[0] )
print( "N_1:\t" , popt[1] )
print( "nu0_1:\t" , popt[2] )
print( "Gamma_1:\t" , popt[3] )
print( "N_2:\t" , popt[4] )
print( "nu0_2:\t" , popt[5] )
print( "Gamma_2:\t" , popt[6] )
print( "N_3:\t" , popt[7] )
print( "nu0_3:\t" , popt[8] )
print( "Gamma_3:\t" , popt[9] )
print( "N_4:\t" , popt[10] )
print( "nu0_4:\t" , popt[11] )
print( "Gamma_4:\t" , popt[12] )
print( "N_5:\t" , popt[13] )
print( "nu0_5:\t" , popt[14] )
print( "Gamma_5:\t" , popt[15] )
print( "N_6:\t" , popt[16] )
print( "nu0_6:\t" , popt[17] )
print( "Gamma_6:\t" , popt[18] )
print( "N_bkg:\t" , popt[19] )
print( "nu0_bkg:\t" , popt[20] )
print( "Gamma_bkg:\t" , popt[21] )




# =============================================================================
# # Define function of a peak
# def maxbol(v, N, m, v_1, k, T):
#     return N*v**3*np.exp((-m*v_1**2)/(2*k*T))
# 
# # Define third function describing the angle shift
# def shift(v_0 ,l , theta):
#     return v_0 * l * np.sin(theta)
# 
# # Define function describing the background
# #def background(ofsset, B, su, su0_bkg, zeta):
# #    return ofsset - B*np.exp(-(su-su0_bkg)**2/2/zeta**2)
# 
# # Define function describing velocity distribution + angle shift + background
# def funk2(v, N, m, v_1, k, T, v_0, l, theta, ofsset, B, su, su0_bkg , zeta):
#     returnValue  = 0
#     returnValue += maxbol(v, N, m, v_1, k, T)
#     returnValue *= shift(v_0 ,l , theta)
# #    returnValue += background(ofsset, B, su, su0_bkg, zeta)
#     return returnValue
# 
# # v: velocity of atoms
# # m: mass of atoms
# # k: Boltzmann constant
# # T: Temperature beam
# # C: Normierungskonstante
# # l: wave vektor
# # start parameter of maxwell boltzmann with doppler broadening
# N       = 0
# m       = 1.4192261*10**(-22)
# v_1     = 1
# k       = 1.380649*10**(-23)
# T       = 600
# C       = 1
# v_0     = 1
# l       = 1
# theta   = 3
# su      = 1
# ofsset  = 0
# B       = 1.07
# su0_bkg     = 384.22919
# zeta    = 0.0002
# 
# param = ( N, m, k, T, C, l , theta,  su , 
#          su, ofsset, B, su0_bkg, zeta)
# pop, pcov = curve_fit( funk2 , freq_short , Vhot_short , p0 = param)
# print( "N:\t" , pop[0] )
# print( "m:\t" , pop[1] )
# print( "k:\t" , pop[2] )
# print( "T:\t" , pop[3] )
# print( "l:\t" , pop[4] )
# print( "theta:\t" , pop[5] )
# =============================================================================




# =============================================================================
# #Define function of a peak
# def doppler(N, f_R, f_L, m, k, T, omega):
#     k_1 = 2 * np.pi * 384.22919
#     v_0 = ((f_R - f_L)*2*np.pi )/(k_1*np.sin(omega))
#     print( v_0 )
#     print(k_1)
#     return N * v_0**3 * np.exp( -(m*v_0**2)/(2*k*T) )
# 
# # Define function describing velocity distribution + angle shift
# def funk2(f_L, N, f_R, T, omega):
#     m       = 1.4192261e-25
#     k       = 1.380649e-23
#     returnValue = doppler(N, f_R, f_L, m, k, T, omega)
#     print( returnValue )
#     return returnValue
# 
# # N = N*2*np.pi/(k*np.sin(omega)) : Normierungskonstante
# # f_R:      Resonanzfrequenz atom   
# # f_L:      Frequenz Laser      
# # m:        mass of atoms
# # k:        Boltzmann constant
# # T:        Temperature beam
# # omega:    angle beams
#     
# # start parameter of maxwell boltzmann with doppler broadening
# N       = 10e13
# f_R     = 384.2291
# T       = 600
# omega   = 0.1
# 
# param = (N, f_R, T, omega)
# pop, pcov = curve_fit( funk2 , freq_short , Vhot_short , p0 = param)
# #print( "N:\t" , pop[0] )
# #print( "f_R:\t" , pop[1] )
# #print( "T:\t" , pop[2] )
# #print( "omega:\t" , pop[3] )
# =============================================================================


# omega:    Frequenz Atom
# omegaL:   Frequenz Laser
# v:        Geschwindigkeit Atome
# alpha:    Winkel zwischen Atom- und Laserstrahl
# m:        mass atom
# T:        temperature beam
# kb:       Boltzmann constant
# c:        light speed
# norm:     Normierungskonstante
def boltzmann(omega, omegaL, alpha, c, m, T, kb, norm):
    v = (2*c*omegaL**2*np.sin(alpha) + np.sqrt(2) * np.sqrt(c**2*omega**2 * (2*omega**2 - omegaL**2 + omegaL**2 *np.sin(2*alpha))))/(2*(omega**2 + omegaL**2 *np.sin(alpha)**2))
    exponent = - v**2 * m / (2 * kb * T)
    shift = np.abs(omegaL * ((v)/(c**2*(1-(v**2)/(c**2))**(3/2))
                             -(v**2*np.sin(alpha))/(c**3*(1-(v**2)/(c**2))**(3/2)) 
                             -(np.sin(alpha))/(c*np.sqrt(1-(v**2/c**2)))))
    omegaL*(1/(np.sqrt(1-(v/c)**2))*np.sin(alpha)-1)
    print("v", v)
    print("exponent", exponent)
    print("shift", shift)
    plt.plot(omegaL, v)
    plt.show()
    return norm * v**3 * np.exp(exponent) / shift



# Define function describing velocity distribution + angle shift
def funk(omegaL, alpha, T, omega1, omega2, omega3, norm1, norm2, norm3):
    m       = 1.4192261e-25
    kb      = 1.380649e-23
    c       = 299792458
    returnValue  = 0
    returnValue += boltzmann(omega1, omegaL, alpha, c, m, T, kb, norm1)  
    returnValue += boltzmann(omega2, omegaL, alpha, c, m, T, kb, norm2)  
    returnValue += boltzmann(omega3, omegaL, alpha, c, m, T, kb, norm3)  
    print(returnValue)
    return returnValue

# start parameter of maxwell boltzmann with doppler broadening
omega1       = 384.2291
omega2       = 384.22925
omega3       = 384.22935
alpha        = 0.1
T            = 600
norm1        = 10e13
norm2        = 10e13
norm3        = 10e13

param = (alpha, T, omega1, omega2, omega3, norm1, norm2, norm3)
#pop, pcov = curve_fit( funk , freq_short , Vhot_short , p0 = param)
#print( "omega1:\t" , pop[0] )

funk( freq_short , *param )

# =============================================================================
# # Plot the spectra vs frequencies
# fig, ax = plt.subplots()
# ax.ticklabel_format(useOffset=False)
# #plt.plot( freq_short , Vspec_short , '.' , label="DFSS")
# #plt.plot( freq_short , func( freq_short , *params ) , '-' , label="start parameter")
# #plt.plot( freq_short , func( freq_short , *popt) , 'r-' , label="Fit by PC")
# plt.plot( freq_short , Vhot_short  , '.' , label="Hot atoms")
# plt.plot( freq_short , funk( freq_short , *param ) , 'r' , label="Hot Fit start parameter")
# #plt.plot( freq_short , funk( freq_short , *pop) , 'y' , label="Hot atom fit")
# plt.xlabel("frequency f [THz]")
# plt.ylabel("measured spectrum [A.U.]")
# plt.legend()
# #plt.show()
# # get current figure
# figure = plt.gcf()  
# figure.set_size_inches(8,5)
# plt.savefig("name.png")
# =============================================================================
