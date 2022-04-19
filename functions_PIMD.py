# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import time

# Conversions
conv = 2625.499638   # [1 Hartree = 2625.499638 kJ/mol]    or the other way arround 3.8088E-4
conv_1 = 96.4853     # eV to kJ/mol
# Constants
part = 3                     # Quantum Particles
kb = 3.16683 * 10 ** (-6)    # Boltzman in Hartree/K
hbar, m = 1, 1               # hbar and mass
# For Exitons                                            from [ ... ]
k_exitons = 9.64931080e-07
w = math.sqrt(k_exitons / 0.84*m)
hw_exitons = 2.579052069  #  (26.73mev) -  0.02673eV * conv_1   [kJ/mol]

#                                                       until here [ ... ]
# For Barak's Fermions Article                           from [ ... ]
k_harmonic_fermions = 1.21647924 * 10 ** (-8)  # spring constant Hartree/Bohr^2
w = math.sqrt(k_harmonic_fermions / m)
hw = hbar*w  # 0.0001104536101718726   # [sqrt(Hartree /kg*K]
hw_2 = 0.2894557624503526   # 0.003 * conv_1  # 0.289458 kJ/mol    3meV as in the Article (to have same number as Barak)
#                                                       until here [ ... ]
# Cut the data until energy minimization
cut_data = 10000  # for data
cut_log = cut_data - 1       # for pimdb log the minus 1 is for the data to be with same index
cutoff = 0  # for block avg

#################################################################################################################

def hw_to_k_const(hw_ev):
    '''
    :param hw: hbar omega of the harmonic oscillator in meV
    conv_1 = 96.4853  -> eV to kJ/mol
    conv = 2625.499638   [1 Hartree = 2625.499638 kJ/mol]
    :return: k - spring constant in Hartree / Bohr^2 (for LAMMPS unit: electron)
    '''
    hw_ev = hw_ev * 0.001  # from meV to eV
    w0 = hw_ev * conv_1 / conv
    k = m * w0**2
    return k


def temp_for_lammps(hw_ev, C):
    '''
    :param C: hw_ev in meV ; number of bhw desired (for example bhw = 6 = C)
    :return: Returns the temperature [K] to put in LAMMPS and beta for python code [[1 / kJ/mol]
    '''
    # Conversion from Hartree to kJ/mol
    hw_ev = hw_ev * 0.001   # from meV to eV
    kb1 = kb * conv  # [kJ/mol/K]
    hw_kj_mol = hw_ev * conv_1
    b = C / hw_kj_mol # [1 / kJ/mol]
    T = 1/(kb1*b)
    return T, b


def statistical_error_estimation(fermionic_energy, wj):
    '''
    :param fermionic_energy: list of the calculated energy of same simulation with different seed
    :param wj: list of the calculated Wj of same simulation with different seed
    :return: The average energy and stdv
    '''
    wj_2 = wj**2
    sum_wj = sum(wj)
    avg_energy = sum(fermionic_energy * wj) / sum_wj
    n = sum_wj**2 / sum(wj_2)
    numerator = 0
    for j in range(len(fermionic_energy)):
        numerator += wj[j] * (fermionic_energy[j] - avg_energy)**2
    stdv = (n/(n-1)) * (numerator / sum_wj)
    energy_error = np.sqrt(stdv) / np.sqrt(n)
    return avg_energy, energy_error


def block_averaging(cutoff, block_size, data):
    '''

    :param cutoff: Cut the cutoff amount of data from the beggining of the simulation since it has not achieved equilibrium.
    :param block_size: obtained from a STDV vs Blocksize simulation to get the corresponding STDV
    :param data: data that is evaluated
    :return: the average and stdv of the data.
    '''
    data_cut = data[cutoff:]
    data_len = len(data_cut)
    number_of_blocks = int(data_len/block_size)
    average_array = np.zeros(number_of_blocks)
    for i in range(number_of_blocks):
        average_array[i] = (data_cut[i * block_size:(i + 1) * block_size]).mean()
    averages = average_array.mean()
    stdv = np.std(average_array, ddof=1) / np.sqrt(number_of_blocks)

    return number_of_blocks, averages, stdv


def make_block_average(cutoff, data):
    block_size_array = np.linspace(1, 10000, 200).astype(int)
    avg_array = np.zeros(len(block_size_array))
    stdv_array = np.zeros(len(block_size_array))
    number_of_blocks_array = np.zeros(len(block_size_array))

    for i, value in enumerate(block_size_array):
        number_of_blocks, avg, stdv = block_averaging(cutoff, block_size=value, data=data)
        avg_array[i] = avg
        stdv_array[i] = stdv
        number_of_blocks_array[i] = number_of_blocks


    figblock = plt.figure()
    plt.plot(block_size_array, stdv_array, label="stdv", color="red")
    plt.xlabel("Block Size")
    plt.ylabel("STDV")
    plt.legend()
    plt.show()
# To run block averaging
    # potenital_energy, time_step, trap, newvir, trapvir = etot_b_2d_harmonic(number_of_files, path)
    # make_block_average(cutoff, pot)


def read_file_pot(filename, cut, s_rows, s_footer):
    '''
    :param filename: input log.lammps.{} file
    :param cut: The first "cut" values will be eliminated.
    :param s_rows: skiprows in log.lammps file since there is no data
    :param s_footer: skipfooter in log.lammps file since there is no data
    :return: potenital_energy, time_step, trap, newvir, trapvir
    '''
    df = pd.read_csv(filename, sep='\s+', engine='python', skiprows=s_rows, skipfooter=s_footer) # delimiter='\s+'
    time_step = df['Time'][cut:]
    potenital_energy = df['PotEng'][cut:]
    trap = df['c_trap'][cut:]
    trapvir = df['c_trapvir'][cut:]
    newvir = df['v_newvir'][cut:]
    return potenital_energy, time_step, trap, newvir, trapvir


def avg_sign(filename, beta, cut):
    '''
    :param filename: pimdb.log file
    :param beta: beta in 1 / kJ/mol
    :param cut: he first "cut" values will be eliminated.
    :return: the average sign value
    '''
    # Here the pimdb.log file is in units of kJ/mol hence beta has to also be in same units
    df = pd.read_csv(filename, delimiter='\s+')
    df = df.iloc[cut:-1] # I added this and removed "cut from all [] example  e_1_1 = df.iloc[cut:, 0]
    v_n_column = df.iloc[:, -1]

    w_b = np.exp(- beta * v_n_column)         #  W(3)_Bosons
    print(len(w_b))
    e_1_1 = df.iloc[:, 0]
    e_2_2 = df.iloc[:, 1]
    e_2_1 = df.iloc[:, 2]
    e_3_3 = df.iloc[:, 3]
    e_3_2 = df.iloc[:, 4]
    e_3_1 = df.iloc[:, 5]
                                             # W(3) Fermions
    w_f_0 = 1
    w_f_1 = np.exp(- beta * e_1_1)
    w_f_2 = 0.5 * ((np.exp(- beta * e_2_1) * w_f_1) - (np.exp(- beta * e_2_2) * w_f_0))
    w_f_3 = (1/3) * ((np.exp(- beta * e_3_1) * w_f_2) - (np.exp(- beta * e_3_2) * w_f_1) + (np.exp(- beta * e_3_3) * w_f_0))

                                             #  W(3)_Bosons second way
    # w_b_0 = 1
    # w_b_1 = np.exp(- beta * e_1_1)
    # w_b_2 = 0.5 * ((np.exp(- beta * e_2_1) * w_b_1) + (np.exp(- beta * e_2_2) * w_b_0))
    # w_b_3 = (1/3) * ((np.exp(- beta * e_3_1) * w_b_2) + (np.exp(- beta * e_3_2) * w_b_1) +
    #                  (np.exp(- beta * e_3_3) * w_b_0))

    s = w_f_3 / w_b
    # s1 = w_f_3 / w_b_3
    return s


def etot_b_2d_harmonic(number_of_files, path):  # Bosons
    '''
    In the Boson Ensemble
    :param number_of_files: number of beads
    :param path: the file path
    :return: time_step, pot, avg_pot_exp, stdv_pot, trap, trapvir, newvir
    '''
    file = [path+'log.lammps.{}'.format(k) for k in range(0, number_of_files)]  # Makes a list of log.lammps.{}
    pot, trap, newvir, trapvir = 0, 0, 0, 0
    time_step = 0
    for i in range(0, number_of_files):
        try:
            potential, time_step, trap1, newvir1, trapvir1 = read_file_pot(file[i], cut_data, 153, 38)  # Harmonic (153,38) Auxiliary(167, 38)
        except:
            potential, time_step, trap1, newvir1, trapvir1 = read_file_pot(file[i], cut_data, 167, 38)
        pot += potential
        trap += trap1
        newvir += newvir1
        trapvir += trapvir1
    number_of_blocks, avg_pot_exp, stdv_pot = block_averaging(cutoff, block_size=900, data=pot)

    return time_step, pot, avg_pot_exp, stdv_pot, trap*conv, trapvir*conv, newvir*conv


def etot_f_2d_harmonic(number_of_files, path, beta):   # Fermions
    '''
    :param number_of_files: number of beads
    :param path: location of path
    :param beta: beta in 1 / kJ/mol
    :return: Fermionic Energy for Harmonic Oscilator potentail
    '''
    # Fermions Sign Re-Weighting # since pimdb.log is bigger file than log.lammps then I reset indexes
    file_pimdb = path + 'pimdb.log'
    time_step, pot, avg_pot_exp, stdv_pot, trap, trapvir, newvir = etot_b_2d_harmonic(number_of_files, path)
    sign_array = avg_sign(file_pimdb, beta, cut_log).reset_index(drop=True)   # Sign average
    wj = sum(sign_array)    # Weights for sign average of each simulation
    e_s_num_h = (((trap.reset_index(drop=True) / number_of_files)) + newvir.reset_index(drop=True)) * sign_array    # (orange) for harmonic potential
    e_s_num_h = np.mean(e_s_num_h)
    s_denom = sign_array
    s_denom = np.mean(s_denom)

    avg_etot_f_h = e_s_num_h / s_denom

    return time_step,  avg_etot_f_h, sign_array, wj


def etot_f_2d_harmonic_aux(number_of_files, path, beta):   # Fermions
    '''
    :param number_of_files: number of beads
    :param path: location of path
    :param beta: beta in 1 / kJ/mol
    :return: Fermionic Energies for the Auxiliary system and the Bogoliubov ineq. term
    '''
    # Fermions Sign Re-Weighting # since pimdb.log is bigger file than log.lammps then I reset indexes
    file_pimdb = path + 'pimdb.log'
    time_step, pot, avg_pot_exp, stdv_pot, trap, trapvir, newvir = etot_b_2d_harmonic(number_of_files, path)
    sign_array = avg_sign(file_pimdb, beta, cut_log).reset_index(drop=True)               # Sign average
    wj = sum(sign_array)   # Weights for sign average of each simulation

                           # Below - a - Auxiliary: (grey) for sum of potentials that are position dependent
    e_s_num_a = ((pot.reset_index(drop=True) / number_of_files) + newvir.reset_index(drop=True)) * sign_array
    e_s_num_a = np.mean(e_s_num_a)
                           # Below - b - bogoliubov: (purpule) <V_g>_H'   (needs to be subtracted! )
    e_s_num_b = (pot.reset_index(drop=True) - trap.reset_index(drop=True)) * sign_array
    e_s_num_b = np.mean(e_s_num_b)

    s_denom = sign_array
    s_denom = np.mean(s_denom)

                           # Below a is for <E>_H' and b for <V_g>_H'
    avg_etot_f_a = (e_s_num_a ) / s_denom   # Auxiliary
    avg_etot_f_b = (e_s_num_b ) / s_denom   # Bogoliubov



    return time_step, sign_array, wj, avg_etot_f_a, avg_etot_f_b


def ana_3_fermions(beta, hbarw):
    '''
    :return: Analytical result for three Fermions in Canonical Ensemble <E> in harmonic potential
    '''
    exp = np.exp(beta * hbarw)
    exp2 = np.exp(2 * beta * hbarw)
    exp3 = np.exp(3 * beta * hbarw)
    exp4 = np.exp(4 * beta * hbarw)
    exp5 = np.exp(5 * beta * hbarw)
    exp6 = np.exp(6 * beta * hbarw)
    exp7 = np.exp(7 * beta * hbarw)
    exp8 = np.exp(8 * beta * hbarw)
    num = hbarw*(5*exp6 + 31*exp5 + 47*exp4 + 50*exp3 + 47*exp2 + 31*exp + 5)
    denom = (exp-1)*(exp+1)*(exp2+exp+1)*(exp2+4*exp+1)
    etot_3_f = num/denom
    return etot_3_f

# def etot_f_2d_harmonic_aux(number_of_files, path, beta):   # Fermions
#     '''
#     :param number_of_files: number of beads
#     :param path: location of path
#     :param beta: beta in 1 / kJ/mol
#     :return: Fermionic Energies for the Auxiliary system and the Bogoliubov ineq. term
#     '''
#     # Fermions Sign Re-Weighting # since pimdb.log is bigger file than log.lammps then I reset indexes
#     file_pimdb = path + 'pimdb.log'
#     time_step, pot, avg_pot_exp, stdv_pot, trap, trapvir, newvir = etot_b_2d_harmonic(number_of_files, path)
#     sign_array = avg_sign(file_pimdb, beta, cut_log).reset_index(drop=True)               # Sign average
#     pot_array = pot.reset_index(drop=True)
#     wj = sum(sign_array)   # Weights for sign average of each simulation
#     # Below - a - Auxiliary: (grey) for sum of potentials that are position dependent
#     e_s_num_a = ((pot.reset_index(drop=True) / number_of_files) + newvir.reset_index(drop=True)) * sign_array
#     number_of_blocks, avg_pot_f_num_a, stdv_pot_f_num = block_averaging(cutoff, block_size=900, data=e_s_num_a)
#     # Below - b - bogoliubov: (purpule) <V_g>_H'   (needs to be subtracted! )
#     e_s_num_b = (pot.reset_index(drop=True) - trap.reset_index(drop=True)) * sign_array
#     number_of_blocks, avg_pot_f_num_b, stdv_pot_f_num = block_averaging(cutoff, block_size=900, data=e_s_num_b)
#
#     s_denom = sign_array
#     number_of_blocks, avg_sign_f_denom, stdv_pot_f_denom = block_averaging(cutoff, block_size=900, data=s_denom)
#     # I am multipling by the conversion from Hartree to kJ/mol since the potential comes from the lammps.log files
#     # and these are in units of Hartree.
#     # Below a is for <E>_H' and b for <V_g>_H'
#     avg_etot_f_a = (avg_pot_f_num_a ) / avg_sign_f_denom   # Auxiliary
#     avg_etot_f_b = (avg_pot_f_num_b ) / avg_sign_f_denom   # Bogoliubov
#
#     number_of_blocks, avg_trap_f_num, stdv_trap_f_num = block_averaging(cutoff, block_size=900, data=trap)
#     number_of_blocks, avg_trapvir_f_num, stdv_trapvir_f_num = block_averaging(cutoff, block_size=900, data=trapvir)
#     number_of_blocks, avg_newvir_f_num, stdv_newvir_f_num = block_averaging(cutoff, block_size=900, data=newvir)
#
#     return time_step, sign_array, wj, avg_etot_f_a, avg_etot_f_b, avg_trap_f_num, \
#            avg_trapvir_f_num, avg_newvir_f_num

