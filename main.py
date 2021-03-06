# module for calculating runtime of simulation
from datetime import datetime
# module for creating pathes, directories and folders
import pathlib
import numpy as np
import scipy.constants as scc
# math module provides mathematical operations like power, squareroot and so on
import math
# random needed for drawing random numbers which are e.g. uniformly or gaussian distributed
import random
# argparse is used for parsing command line arguments
import argparse
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# import of self written modules, first one handles the needed calculations around the Maxwell Boltzmann disribution
from distributions import number_sampling
# atomic material module holds informations (energy gap between ground and excited state, mass, einstein coefficents,...) about the used atom sorts
from atomic_material import Lithium_6
# util contains all functions which are needed in some special case but do not belong to any superordinate topic
from util import calculate_p_max
# this module contains all informations and all calculations regarding the magnetic field
from magnetic_field import spline_fit_field_function, max_step_length_fit_field_function
# import specific plotting methods used in this masterthesis
from plotting import line_plotting, hist_plotting, slice_plotting, eval_plotting, state_occupation_development_plot
# script for logging simulation information to a log file
from log_file import create_log_file
# module for loading all settings
from load_settings import load_files
import os
# module for excitation and de-excitation processes
from Exciting_atoms import Excitation

# import numba for JIT compiling
from numba import jit, int32, float64, void, cuda

# timestep function does the simulation steps for discrete time steps
@jit(nopython=True)
def timestep(pol,laser_frequency,laser_detuning, atom_count, p_max, v_min, v_max, x_min, x_max, y_min, y_max, wavelength,
             laser_beam_radius, zeeman_distance, target_center_z, spline_fit,target_radius, intensity, max_step_fit_function,
             slicing_position_array):
    if atom_count<2000:
        z_step=0.000000001
    else:
        z_step=0.01

    cutoff_number = 10001
    r_target = target_radius
    loop_counter = 0
    # initialize mot_counter for counting atoms entering mot, the indices of the atoms in
    # the MOT and three lists for saving velocity components of these atoms + one list for saving velocity
    atoms_in_mot = 0
    vel_x_atoms_in_mot = []
    vel_y_atoms_in_mot = []
    vel_z_atoms_in_mot = []
    vel_z_atoms_in_mot_upper_groundstate = []
    vel_z_atoms_in_mot_lower_groundstate = []
    dead_atoms_upper_groundstate = []
    dead_atoms_lower_groundstate = []

    plane_slice_pos = slicing_position_array
    plane_slice_flags = []
    plane_slice_upper_groundstate = [[0.0],[0.0],[0.0],[0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0]]
    plane_slice_lower_groundstate = [[0.0],[0.0],[0.0],[0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0],[0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0]]

    for i in range(len(plane_slice_pos)):
        plane_slice_flags.append(0)

    # counter for counting how many times the atoms were in total excited
    excitation_counter = 0

    capture_count_z_vel = 0
    wavevector_x = 0
    wavevector_y = 0
    wavevector_z = -1

    groundstate_upper_lower_start=[]
    groundstate_upper_lower=[]
    start_x_vel_atoms_in_mot = []
    start_x_vel = []
    start_y_vel_atoms_in_mot = []
    start_y_vel = []
    start_z_vel_atoms_in_mot = []
    start_z_vel = []
    start_z_vel_upper_state_atoms_in_mot = []
    start_z_vel_lower_state_atoms_in_mot = []
    dead_atoms_pos=[]
    dead_atoms_vx=[]
    dead_atoms_vy=[]
    dead_atoms_vz=[]

    z_histogram=[[0.0],[0.0],[0.0],[0.0],[0.0],[0.0]]
    vel_z_histo=[[[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0], [0.0], [0.0], [0.0]], [[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0], [0.0], [0.0], [0.0]], [[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0], [0.0], [0.0], [0.0]], [[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0], [0.0], [0.0], [0.0]], [[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0], [0.0], [0.0], [0.0]], [[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0],[0.0], [0.0], [0.0], [0.0]]]
    # initializing lists with 1 for observing several variables. Needed for compiling through Numba
    observing_z_vel = [0.0]
    observing_z_pos = [0.0]
    observing_magnetic_field = [0.0]
    excitation_freq_development = [0.0]
    excitation_probability_development = [0.0]
    zeeman_shift = [0.0]

    light_beam_radius_squared = laser_beam_radius ** 2
    threshold = 1E-35
    mean_free_path_bigger = 0
    mean_free_path_smaller = 0
    # for every atom do the simulation steps
    counter=0
    counter_dead=0

    for i in range(0, atom_count):
        z_pos_histo= [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0, 0.0, 0.0,0.0],
                       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0],
                       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0]]

        counter=0
        if i % 1000 == 0:
            print("Atom number", i)
        mot_flag = 0
        for f in range(len(plane_slice_flags)):
            plane_slice_flags[f] = 0
        x_pos = random.uniform(x_min, x_max + threshold)
        y_pos = random.uniform(y_min, y_max + threshold)
        z_pos = 0.0
        z_pos_before=z_pos

        # only z velocity component
        x_velocity = 0.0
        y_velocity = 0.0
        if atom_count>100:
            z_velocity = number_sampling(1, p_max, v_min, v_max, mass_lithium_6, temperature)
        else:
            z_velocity=2000-i*20

        start_x_velocity = x_velocity
        start_y_velocity = y_velocity
        start_z_velocity = z_velocity

        start_x_vel.append(start_x_velocity)
        start_y_vel.append(start_y_velocity)
        start_z_vel.append(start_z_velocity)

        atom_dead = 1

        random_ground_state = random.random()  # choose GS (all equally likely)
        if allgs==1:
            for state in range(0,7):
                if random_ground_state<(state+1)/6 and random_ground_state>=(state)/6:
                    if state==0 or state==2:
                        current_groundstate=0
                        initial_state=0
                    else:
                        current_groundstate=1
                        initial_state=1
                    GS_quantum_number=state

        elif allgs==0:
            initial_state=1
            current_groundstate=1
            GS_quantum_number=5

        groundstate_upper_lower_start.append(current_groundstate)

        pos_value = 0.0
        pos_index = 0

        z_pos_histo_binning=z_pos
        # if the atom is not dead do the steps
        while atom_dead != 0:
            for iii in range(0,23):
                if iii==pos_index:
                    pos_value=plane_slice_pos[iii]
                    break
            current_excitation_freq = [0.0, 0.0, 0.0, 0.0]
            loop_counter += 1
            if atom_count < cutoff_number:
                if z_pos-z_pos_before>=z_step:
                    observing_z_pos.append(z_pos)
                    observing_z_vel.append(z_velocity)
                    z_pos_before=z_pos
            # derive the squared x- and y-position for the following comparisons if the position in the xy-plane is beyond some (geometrical) limit
            x_y_pos_component_squared = x_pos ** 2 + y_pos ** 2
            # check if atoms move outside the laser beam or hit the wall of the MOT
            if x_y_pos_component_squared > light_beam_radius_squared or z_pos > target_center_z+0.12 or z_pos < 0.0 or z_velocity < 0.0:
                dead_atoms_pos.append(z_pos)
                dead_atoms_vx.append(x_velocity)
                dead_atoms_vy.append(y_velocity)
                dead_atoms_vz.append(z_velocity)
                # set the index to 0 (=dead)
                atom_dead = 0
                counter_dead+=1
                if current_groundstate == 0:
                    dead_atoms_lower_groundstate.append(z_velocity)
                if current_groundstate == 1:
                    dead_atoms_upper_groundstate.append(z_velocity + 2000)
                if atom_count < cutoff_number:
                    observing_z_pos.append(z_pos)
                    observing_z_vel.append(z_velocity)
                    observing_z_pos.append(-1)
                    observing_z_vel.append(-1)

                continue
            # is the atom is not hitting the rear wall of the MOT
            if z_pos < target_center_z + target_radius:
                # check if atom is inside the beam of the laser and if the frequency for exciting the atom is equal to the frequency of the laser
                if x_y_pos_component_squared < light_beam_radius_squared:
                    excitation_index=0

                    # Excitation and De-excitation
                    z_velocity_new,y_velocity_new,x_velocity_new,excitation_time_step, GS_quantum_number, current_groundstate,counter,atom_path_length=\
                        Excitation(repumper,current_groundstate,counter,mass_lithium_6, GS_quantum_number, spline_fit_field_function(spline_fit, z_pos - target_center_z,maximum_distance,target_center_z), z_pos, zeeman_distance,
                                   mean_free_path_bigger, max_step_fit_function, target_center_z, maximum_distance, x_velocity, y_velocity,
                                   z_velocity, wavevector_x, wavevector_y, wavevector_z, wavelength, pol, excitation_counter, mean_free_path_smaller,
                                   laser_frequency, x_y_pos_component_squared, intensity, laser_detuning, natural_line_width)

                    if z_pos==0.0 or z_pos>=z_pos_histo_binning+0.005:
                        z_histogram[GS_quantum_number].append(z_pos)
                        z_pos_histo_binning=z_pos
                    if z_pos>=pos_value and z_pos_histo[GS_quantum_number][pos_index]==0.0:
                        z_pos_histo[GS_quantum_number][pos_index]=pos_value
                        vel_z_histo[GS_quantum_number][pos_index].append(z_velocity_new)
                        if pos_index+1<len(plane_slice_pos): pos_index+=1
                    atom_time_step=atom_path_length/(math.sqrt(z_velocity**2+y_velocity**2+x_velocity**2+1e-30))

                    # position updating for all atoms before calculating new velocity depending of effects occuring
                    x_pos += x_velocity * atom_time_step
                    y_pos += y_velocity * atom_time_step
                    z_pos += z_velocity * atom_time_step

                    if z_pos>=target_center_z-0.01 and MOT_field==1: #1cm before MOT center: not transverse heating
                        x_velocity=0.0
                        y_velocity=0.0
                    else:
                        x_velocity=x_velocity_new
                        y_velocity=y_velocity_new

                    z_velocity=z_velocity_new

                    if atom_count < cutoff_number:
                        if z_pos - z_pos_before >= z_step:
                            excitation_freq_development.append(current_excitation_freq[excitation_index])
                            observing_z_pos.append(z_pos)
                            observing_z_vel.append(z_velocity)
                            observing_magnetic_field.append(spline_fit_field_function(spline_fit, z_pos - target_center_z,maximum_distance,target_center_z))
                            z_pos_before=z_pos
                    for ii in range(len(plane_slice_pos)):
                        if plane_slice_pos[ii] < z_pos and plane_slice_flags[ii] != 1:
                            if current_groundstate == 0:
                                plane_slice_lower_groundstate[ii].append(z_velocity)
                            if current_groundstate == 1:
                                plane_slice_upper_groundstate[ii].append(z_velocity + 2000)
                            plane_slice_flags[ii] = 1
                else:
                    maximum_path_length = max_step_length_fit_field_function(max_step_fit_function,z_pos - target_center_z, maximum_distance)
                    v_ges=math.sqrt(x_velocity**2+y_velocity**2+z_velocity**2+1e-35)
                    x_pos+=x_velocity/v_ges *maximum_path_length
                    y_pos+=y_velocity/v_ges *maximum_path_length
                    z_pos+=z_velocity/v_ges *maximum_path_length

            # atoms entering mot, therefor save number of entering atoms and their velocity components
            if target_center_z <= z_pos <= target_center_z + r_target and x_y_pos_component_squared < r_target ** 2 and mot_flag == 0:
                mot_flag = 1
                atoms_in_mot += 1
                vel_x_atoms_in_mot.append(x_velocity)
                start_x_vel_atoms_in_mot.append(start_x_velocity)
                vel_y_atoms_in_mot.append(y_velocity)
                start_y_vel_atoms_in_mot.append(start_y_velocity)
                vel_z_atoms_in_mot.append(z_velocity)
                start_z_vel_atoms_in_mot.append(start_z_velocity)
                groundstate_upper_lower.append(current_groundstate)

                if initial_state == 1:
                    start_z_vel_upper_state_atoms_in_mot.append((start_z_velocity + 2000))
                else:
                    start_z_vel_lower_state_atoms_in_mot.append(start_z_velocity)

                if current_groundstate == 1:
                    vel_z_atoms_in_mot_upper_groundstate.append(z_velocity + 2000)
                else:
                    vel_z_atoms_in_mot_lower_groundstate.append(z_velocity)

                if atom_count < cutoff_number:
                    observing_z_pos.append(-1)
                    observing_z_vel.append(-1)

                atom_dead = 0
                continue

    print("counter",counter)
    print("dead atoms in lower and upper gs:",counter_dead)

    return dead_atoms_pos,dead_atoms_vx,dead_atoms_vy,dead_atoms_vz,start_x_vel,start_y_vel,start_z_vel,groundstate_upper_lower_start,groundstate_upper_lower,z_histogram, vel_z_histo, atoms_in_mot, excitation_counter, capture_count_z_vel, observing_z_pos, observing_magnetic_field, excitation_freq_development, excitation_probability_development, observing_z_vel, vel_x_atoms_in_mot, vel_y_atoms_in_mot, vel_z_atoms_in_mot, vel_z_atoms_in_mot_upper_groundstate, vel_z_atoms_in_mot_lower_groundstate, start_z_vel_atoms_in_mot, start_z_vel_upper_state_atoms_in_mot, start_z_vel_lower_state_atoms_in_mot, zeeman_shift, loop_counter, start_z_vel, plane_slice_upper_groundstate, plane_slice_lower_groundstate, dead_atoms_upper_groundstate, dead_atoms_lower_groundstate



MOT_field=1 # 1: with MOT field; 0: without
repumper="off" # on / off

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Load necessary information from files')
    parser.add_argument('--sim_params_file',
                        help='Specify the JSON-file containing information about the simulation parameters')
    parser.add_argument('--exp_params_file',
                        help='Specify the JSON-file containing information about the experiment setup')
    parser.add_argument('--atom_params_file',
                        help='Specify the JSON-file containing information about the atom parameters')
    parser.add_argument('--raw_magnetic_field_data',
                        help='Specify the location of the TXT-file containing the raw magnetic field data')
    parser.add_argument('--max_step_size',
                        help='Specify the location of the TXT-file containing the maximum step size for a given magnetic field')
    args = parser.parse_args()
    sim_param_data, exp_param_data, atomic_data, spline_fit, file_name_magnetic_field, max_step_length_file, maximum_distance = load_files(args.sim_params_file,
                                                                                                                                           args.exp_params_file,
                                                                                                                                           args.atom_params_file,
                                                                                                                                           args.raw_magnetic_field_data,
                                                                                                                                           args.max_step_size)
    # mass of observed atom
    mass_lithium_6 = Lithium_6.mass_lithium_six
    # temperature at which atom species vaporises
    temperature = sim_param_data['temperature']
    # number of observed atoms
    n = sim_param_data['particle_number']
    allgs= 1 #0=gs5 only, 1=all gs
    # minimal considered velocity
    v_min = sim_param_data['velocity_min']
    # maximal considered velocity
    v_max = sim_param_data['velocity_max']
    # threshold for limits of used distributions
    threshold = 0.000001
    # minimal and maximal starting positions of atoms
    y_min = exp_param_data["center_atomic_source"] - 0.5 * exp_param_data["y_expansion_atomic_source"]
    y_max = exp_param_data["center_atomic_source"] + 0.5 * exp_param_data["y_expansion_atomic_source"]
    x_min = exp_param_data["center_atomic_source"] - 0.5 * exp_param_data["x_expansion_atomic_source"]
    x_max = exp_param_data["center_atomic_source"] + 0.5 * exp_param_data["x_expansion_atomic_source"]

    # zeeman slower distance
    zeeman_distance = exp_param_data["zeeman_slower_distance"]
    target_center_z = exp_param_data["mot_distance"] #equal to length of the slower
    target_radius = exp_param_data["mot_radius"]
    # educated guessing
    bin_count = 80

    # laser properties
    laser_det = (sim_param_data["slower_laser_detuning"]) # Frequency in 1/s (c/lambda)
    laser_freq = (sim_param_data["slower_laser_frequency"])  # Frequency in 1/s (c/lambda)
    laser_pol = (sim_param_data["laser_polarisation"])  # laser pol: sigminus, pi, sigplus
    wavelength = scc.c / laser_freq  # change wavelength, as its connected to f
    natural_line_width = 2 * math.pi * 5.87E6
    intensity = sim_param_data['slower_laser_intensity'] #in params of sat intensity
    laser_beam_radius = sim_param_data['slower_laser_diameter'] / 2
    slicing_positions = sim_param_data['positions_for_slicing']

    # test area, new JSON parameters read in
    print("slicing positions",slicing_positions)
    print("laser polarisation", sim_param_data['laser_polarisation'])
    print("Frequency offset", atomic_data['frequency_offset_ground_state'])

    p_max = calculate_p_max(n, v_min, v_max, mass_lithium_6, temperature)
    # function call of timestep
    startTime = datetime.now()

    # line breaks used for better readability
    dead_pos,dead_vx,dead_vy,dead_vz,start_vx,start_vy,start_vz,gs_upper_lower_start,gs_upper_lower,z_histo,v_z_histo,atoms_in_mot, excitation_counter, capture_count_z_velocity, observing_z_position, \
    observing_magnetic_field, excitation_freq_development, excitation_probability_development, observing_z_velocity, \
    vel_x_atoms_in_mot, vel_y_atoms_in_mot, vel_z_atoms_in_mot, vel_upper_groundstate, vel_lower_groundstate, \
    start_z_vel_atoms_in_mot, start_vel_upper_state, start_vel_lower_state, zeeman_shift, loop_Counter, start_vel_z, \
    vel_z_plane_slices_upper_gs, vel_z_plane_slices_lower_gs, vel_dead_atoms_upper, \
    vel_dead_atoms_lower = timestep(laser_pol,laser_freq,laser_det, n, p_max, v_min, v_max, x_min, x_max, y_min, y_max,
                                    wavelength, laser_beam_radius, zeeman_distance, target_center_z,
                                    spline_fit, target_radius, intensity,max_step_length_file, slicing_positions)

    runtime = datetime.now() - startTime
    print("runtime", runtime)

    print("Plotting...")
    # create output folder or if folder already exists, create log folder for simulation run
    create_log_file(file_name_magnetic_field, atoms_in_mot, excitation_counter, start_vel_z,
                    capture_count_z_velocity, exp_param_data, sim_param_data, runtime, startTime,
                    bin_count, v_max)
    pathlib.Path('simulation_results/' + str(startTime.strftime("%Y_%m_%d %H_%M_%S"))).mkdir(parents=True,
                                                                                             exist_ok=True)
    with open('simulation_results/' + str(startTime.strftime("%Y_%m_%d %H_%M_%S")) + '/mot_vel_distribution.txt',
              'w') as mot_file:
        for i in range(len(vel_x_atoms_in_mot)):
            mot_file.write(
                str(vel_x_atoms_in_mot[i]) + ";" + str(vel_y_atoms_in_mot[i]) + ";" + str(vel_z_atoms_in_mot[i]) + "\n")

    #plot velocity evolutions of single atoms
    line_plotting(observing_z_position, observing_z_velocity, 'z position', 'z velocity', 0.0,target_center_z+0.005,0.0, 2000.0, startTime, False)

    #plot dead atoms
    fig, ax = plt.subplots()
    ax.hist(dead_pos, bins=100)
    print("number of dead atoms", len(dead_pos))
    plt.xlim(target_center_z-0.5,target_center_z+0.001)
    plt.xlabel("Position in m", fontsize=22)
    plt.ylabel("Number of dead atoms", fontsize=22)
    plt.title("Total number of dead atoms: {} of {}".format(len(dead_pos),n),fontsize=22)
    plt.rcParams.update({'font.size': 22})
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.show()

    #plot velocity distributions for different ground states
    pathlib.Path('simulation_results/v_distr').mkdir(parents=True, exist_ok=True)
    positions=slicing_positions
    pos_i=0
    for pos in positions:
        fig, ax = plt.subplots()
        labels=["GS 0","GS 1","GS 2","GS 3","GS 4","GS 5"]
        colors=["red","cyan","orange","blue","green","purple"]
        bin_size = 10
        min_edge = 0
        max_edge = 5700
        N = int((max_edge-min_edge)/bin_size)
        bin_list = np.linspace(min_edge, max_edge, N+1)
        res, bins, patches=ax.hist([v_z_histo[0][pos_i], v_z_histo[1][pos_i], v_z_histo[2][pos_i],v_z_histo[3][pos_i], v_z_histo[4][pos_i], v_z_histo[5][pos_i]], bins=bin_list, stacked=True,color=colors, label=labels)
        plt.legend(loc="upper right",fontsize=22)
        plt.xlabel("v_z in m/s", fontsize=22)
        plt.ylabel("Atoms in GS", fontsize=22)
        plt.title("Atoms at z={}m".format(round(pos,3)),fontsize=22)
        plt.rcParams.update({'font.size': 22})
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        figure = plt.gcf()  # get current figure
        figure.set_size_inches(13.66, 6.71)
        #plt.ylim(0,350)
        v_z_histo[5][pos_i].sort()
        plt.savefig('simulation_results/' + "v_distr" + "/" + "vz" + "_Histo_pos" + str(round(pos,5)).replace('.', '_') + "_allGS" + ".png")
        plt.close()
        pos_i+=1

    # delete pickle-files
    pickle_path=os.listdir("./sim_setup/")
    for file in pickle_path:
        if file.endswith(".pickle"):
            filepath="./sim_setup/"+file
            os.remove(filepath)