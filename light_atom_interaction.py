import math
from numba import jit
from trans_strength import trans_strength


# function for calculating the Doppler shift
@jit(nopython=True)
def doppler_shift(velocity_x, velocity_y, velocity_z, init_freq, wavevector_x, wavevector_y, wavevector_z, wavelength):
    '''
    :param velocity_x: Float, x-component of the atom's velocity vector.
    :param velocity_y: Float, y-component of the atom's velocity vector.
    :param velocity_z: Float, z-component of the atom's velocity vector.
    :param init_freq: Float, the initial frequency of the atom. Depends on the current ground state of the atom.
    :param wavevector_x: Float, x-component of the wavevector of the laser light.
    :param wavevector_y: Float, y-component of the wavevector of the laser light.
    :param wavevector_z: Float, z-component of the wavevector of the laser light.
    :param wavelength: Float, the wavelength of the laser.
    :return: Float, the Doppler-shifted frequency.

    Calculates the Doppler shift for a moving atom.
    '''
    reduced_freq = init_freq - (wavevector_x * abs(velocity_x) + wavevector_y * abs(velocity_y) + wavevector_z * abs(
        velocity_z)) * 2 * (math.pi / wavelength)

    return reduced_freq

# lorentzian probability distribution for calculating the scatter rate
@jit(nopython=True)
def lorentzian_probability(pol, GS, a, Bfield, atom_freq, laser_frequency, laser_detuning,
                           natural_line_width, laser_intensity, velocity_x, velocity_y, velocity_z, wavevector_x,
                           wavevector_y, wavevector_z, wavelength):
    '''
    :param pol: Int, polariztion (0,1,2 for sigmin, pi, sigplus).
    :param GS: Int, the GS quantum number.
    :param a: Int, the excited state quantum number.
    :param Bfield: Float, the Bfield strength.
    :param atom_freq: Float, the freuqency of the atom.
    :param laser_frequency: Float, the frequency of the laser.
    :param laser_detuning: Float, the detuning of the laser.
    :param natural_line_width: Float, the natural line width of the excited state.
    :param laser_intensity: Float, the intensity of the laser at the current position of the atom.
    :param laser_sat_intensity: Float, the saturation intensity of the atom.
    :return: Float, excitation probability of the atom.

    Calculates the scatter rate using the lorentzian probability.
    '''
    laser_freq_modus = 2 * math.pi * laser_frequency + 2 * math.pi * laser_detuning # Include Laser detuning
    reduced_freq = laser_freq_modus - (
                wavevector_x * abs(velocity_x) + wavevector_y * abs(velocity_y) + wavevector_z * abs(
            velocity_z)) * 2 * (math.pi / wavelength)  # Doppler shift
    delta = reduced_freq - atom_freq
    ##### FACTOR OF 10 in SQRT is necessary but needs to be explained
    omega = 2 * math.pi * 1e6 * 11.925*trans_strength(GS, a, pol, Bfield) * 4.37 * math.sqrt(0.1*laser_intensity) # Rabi frequency mit 11.925Debye
    s_val = 0.5 * omega ** 2 / (delta ** 2 + (natural_line_width ** 2) / 4)
    excitation_rate = s_val / (2 * (1 + s_val))

    return excitation_rate

# lorentzian probability distribution for calculating the probability of excitment
@jit(nopython=True)
def lorentzian_probability_TEST(polarization, pol, GS, a, Bfield, atom_freq, laser_frequency, laser_detuning,
                                natural_line_width, laser_intensity, velocity_x, velocity_y, velocity_z, wavevector_x,
                                wavevector_y, wavevector_z, wavelength): #just for testing stuff
    '''
    :param atom_freq: Float, the Doppler shifted freuqency of the atom.
    :param laser_frequency: Float, the frequency of the laser.
    :param laser_detuning: Float, the detuning of the laser.
    :param natural_line_width: Float, the natural line width of the excited state.
    :param laser_intensity: Float, the intensity of the laser at the current position of the atom.
    :param laser_sat_intensity: Float, the saturation intensity of the atom.
    :return: Float, excitation probability of the atom.

    Calculates the scatter rate using the lorentzian probability.
    '''
    laser_freq_modus = 2 * math.pi * laser_frequency + 2 * math.pi * laser_detuning
    reduced_freq = laser_freq_modus - (
                wavevector_x * abs(velocity_x) + wavevector_y * abs(velocity_y) + wavevector_z * abs(
            velocity_z)) * 2 * (math.pi / wavelength)  # doppler shift
    delta = reduced_freq - atom_freq
    omega = 2 * math.pi * 1e6 * trans_strength(GS, a, pol, Bfield) * 4.37 * math.sqrt(laser_intensity)
    s_val = 0.5 * omega ** 2 / (delta ** 2 + (natural_line_width ** 2) / 4)
    excitation_probability = s_val / (2 * (1 + s_val))

    return reduced_freq, atom_freq, delta


@jit(nopython=True)
def lorentzian_probability_2(atom_freq, laser_frequency, laser_detuning, natural_line_width, laser_intensity,
                             laser_sat_intensity): #old main (Andis simulation)
    '''
    :param atom_freq: Float, the Doppler shifted freuqency of the atom.
    :param laser_frequency: Float, the frequency of the laser.
    :param laser_detuning: Float, the detuning of the laser.
    :param natural_line_width: Float, the natural line width of the excited state.
    :param laser_intensity: Float, the intensity of the laser at the current position of the atom.
    :param laser_sat_intensity: Float, the saturation intensity of the atom.
    :return: Float, excitation probability of the atom.

    Calculates the scatter rate using the lorentzian probability.
    '''

    laser_freq_modus = laser_frequency - 2 * math.pi * laser_detuning  # ohne 2 pi und Vorzeichen geändert

    excitation_probability = 0.5 * (laser_intensity / laser_sat_intensity) * (natural_line_width ** 2) / (4 * (laser_freq_modus - atom_freq) ** 2 + (natural_line_width ** 2) * (1 + laser_intensity / laser_sat_intensity))

    return excitation_probability


@jit(nopython=True)
def laser_intensity_gauss(radius_squared, distance, start_intensity):
    '''
    :param radius_squared: Float, the radial position of the atom in the beam.
    :param distance: Float, the distance of the atom to the laser origin.
    :param start_intensity: Float, the initial intensity of the laser.
    :return: Float, the intensity the atom receives at the current position.

    Calculates the intensity of the laser assuming a Gaussian beam profile.
    '''

    width = distance * (0.0098 / 0.59) + 0.002
    intensity = start_intensity * math.exp(-((radius_squared) / width ** 2))

    return intensity


def scatter_rate(natural_line_width, atom_freq, laser_frequency, rabi_freq, velocity_z, wavelength):
    '''
    :param natural_line_width: Float, the natural line width of the excited state.
    :param atom_freq: Float, the Doppler shifted freuqency of the atom.
    :param laser_frequency: Float, the frequency of the laser.
    :param rabi_freq: Float, the current Rabi frequency of the atom.
    :param velocity_z: Float, z-component of the atom's velocity vector.
    :param wavelength: Float, the wavelength of the laser.
    :return: Flaot, scatter rate.

    Calculates the scatter rate using formula 9.3 in "Atomic physics" by C. J. Foot
    '''
    return (natural_line_width / 2) * (rabi_freq ** 2 / 2) / (
                (laser_frequency - atom_freq + ((2 * math.pi * velocity_z) / wavelength)) ** 2 + (
                    rabi_freq ** 2 / 2) + (natural_line_width ** 2 / 4))
