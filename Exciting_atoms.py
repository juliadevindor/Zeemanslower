import math
from numba import jit
import random
from Position import Position
from trans_strength import trans_strength
from light_atom_interaction import lorentzian_probability,laser_intensity_gauss
import scipy.constants as scc
from magnetic_field import max_step_length_fit_field_function

@jit(nopython=True)
def Excitation(Debug_flag, current_groundstate,counter,mass_lithium_6,GS_quantum_number,Bfield, z_pos, zeeman_distance, mean_free_path_bigger, max_step_fit_function, target_center_z, maximum_dist, vx,vy,vz,kx,ky,kz,wavelength,polarization,excitation_counter,mean_free_path_smaller,laser_frequency,x_y_pos_component_squared,laser_intensity,laser_detuning,natural_line_width,threshold):
    #print("nu_laser",laser_frequency-(kz*1128/wavelength+Position(5,11,2,1200e-4)/(2*math.pi)))
    Bfield+=1e-30
    h_bar = scc.hbar
    deex_rate = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    excitation_rate = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    rho_ex = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    rho_deex = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    excitation_rate_sum = 0.0
    GS = GS_quantum_number  
    if Debug_flag==1: print("=====================================")
    if Debug_flag==1: print("zpos:", z_pos, "GS:", GS)

    if z_pos > zeeman_distance:
        for a in range(0,12): #loop over all excited states
            for pol in range(0, 3):  # loop over all polarizations: 0: sigmin; 1:pi; 2:sigplus
                rho_ex[a]=polarization[pol]*lorentzian_probability(pol,GS,a,Bfield,Position(GS, a, pol, Bfield), laser_frequency, laser_detuning,natural_line_width,laser_intensity_gauss(x_y_pos_component_squared,z_pos, laser_intensity), vx, vy, vz, kx, ky, kz, wavelength)
                ### repumper laser ###
                #rho_ex[a]+=polarization[pol]*lorentzian_probability(pol,GS,a,Bfield,Position(GS, a, pol, Bfield), laser_frequency, laser_detuning+228e6,natural_line_width,laser_intensity_gauss(x_y_pos_component_squared,z_pos, laser_intensity), vx, vy, vz, kx, ky, kz, wavelength)
                ### repumper laser ###
                if rho_ex[a]>1e-20 and Debug_flag==1: print("ES", a, "sigma", pol, "rho ex", rho_ex[a])
                excitation_rate[a]+=rho_ex[a]*natural_line_width #Gamma_i
                if rho_ex[a]>1e-20 and Debug_flag==1: print("exc rate", excitation_rate[a])
            excitation_rate_sum += excitation_rate[a] #Sum over all Gamma_i
            if Debug_flag==1: print("exc rate sum", excitation_rate_sum)

    excitation_rate_sum+=1e-35 #to avoid diving by 0

    rnd_number = random.random()
    excitation_time_step = ((-1 / (excitation_rate_sum + 1e-35)) * math.log(rnd_number + 1e-35))
    if Debug_flag==1: print("exc time step", excitation_time_step)
    maximum_path_length = max_step_length_fit_field_function(max_step_fit_function, z_pos - target_center_z,maximum_dist)
    if Debug_flag==1: print("max path length", maximum_path_length)
    atom_path_length = excitation_time_step*vz
    if Debug_flag==1: print("atom path length",atom_path_length)

    if atom_path_length < maximum_path_length:  # Scatter a photon!
        if Debug_flag==1: print("scatter")
        random_exc_state = random.random()  # random number for excited state
        prob1 = 0  # start of probability-intervall
        prob2 = 0  # end of probability-intervall
        for i in range(0, 12):  # loop over all excited states
            if i > 0:
                prob1 += excitation_rate[i - 1] / excitation_rate_sum
            prob2 += excitation_rate[i] / excitation_rate_sum
            if prob1 < random_exc_state <= prob2:
                ES_quantum_number = i  # excited state after excitation from ground state

        deex_rate_sum = 0.0
        ES = ES_quantum_number
        if Debug_flag==1: print("exc state", ES)

        for GS in range(0, 6): #loop over all ground states
            for pol in range(0, 3):  # loop over all polarizations: 0: sigmin; 1:pi; 2:sigplus
                rho_deex[GS]=trans_strength(GS, ES, pol, Bfield)
                deex_rate[GS]+=rho_deex[GS] * natural_line_width
                if rho_deex[GS]>1e-30 and Debug_flag==1: print("GS", GS, "pol", pol, "rho_deex", rho_deex[GS], "deex rate",deex_rate[GS])
            deex_rate_sum += deex_rate[GS]
        if Debug_flag==1: print("deex_rate_sum",deex_rate_sum)

        random_exc_state = random.random()  # random number for new GS
        prob1=0
        prob2=0
        for i in range(0, 6):
            if i > 0:
                prob1 += deex_rate[i-1] / deex_rate_sum
            prob2 += deex_rate[i] / deex_rate_sum
            if random_exc_state > prob1 and random_exc_state <= prob2:
                GS_quantum_number = i
                if GS_quantum_number==0 or GS_quantum_number==2: #lower GS
                    current_groundstate=0
                else:
                    current_groundstate=1 #upper GS

        if Debug_flag==1: print("groundstate",GS_quantum_number)
        # de-excitation
        excitation_counter += 1
        mean_free_path_smaller += 1

        atom_time_step = excitation_time_step
        vz -= h_bar * ((2 * math.pi) / (wavelength * mass_lithium_6))
        # RIGHT: take two random numbers u and v between 0 and 1
        u = random.random()
        v = random.random()
        # for equally distributed points on a sphere use this two formulas. Source: http://mathworld.wolfram.com/SpherePointPicking.html
        random_number_theta = 2 * math.pi * u
        random_number_phi = math.acos(2 * v - 1)
        # Assigning velocity components as in the spherical coordinate system
        vx -= math.cos(random_number_theta) * math.sin(random_number_phi) * scc.h / (
                mass_lithium_6 * wavelength)
        vy -= math.sin(random_number_theta) * math.sin(random_number_phi) * scc.h / (
                mass_lithium_6 * wavelength)
        vz -= math.cos(random_number_phi) * scc.h / (mass_lithium_6 * wavelength)
        counter+=1

    else: # Do not scatter, just continue
        atom_path_length = maximum_path_length
        mean_free_path_bigger += 1

    return vz,vy,vx,excitation_time_step, GS_quantum_number, current_groundstate,counter,atom_path_length
