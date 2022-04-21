# Attacking the Fermion Sign Problem
from functions_PIMD import *
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import fmin
start = time.time()



bhw = [1.5]
temp = [23.2091]  # [K]
beta = 5.184154
# Path of Folder and how many Beads (number of files) (For Harmonic num_of_files is 12 for Auxiliary is 72)
number_of_files = 32
path = '/home/netanelb2/Desktop/Netanel/Research/sgnprob/'     # Auxiliary
file_pimdb = path + 'pimdb.log'


#######################################################################################################################

#                                 ###      2 Bosons /  Fermions under Harmonic Trap hw=3meV         ###

#                                 ###  Bosons   ###

time_step, pot_b, avg_pot_exp, stdv_pot, trap_b, trapvir_b, newvir_b = etot_b_2d_harmonic(number_of_files, path)
print("Total Energy Fermions / P: ", (2 * avg_pot_exp / (number_of_files * hw)))  # [Kj/mol]


#                                 ###   Fermions   ###
time_step,  avg_etot_f_h, sign_array, wj = etot_f_2d_harmonic(number_of_files, path, beta)
print("<S>_B: ", np.mean(sign_array))  # [Kj/mol]
print("W_j: ", wj)  # [Kj/mol]
 #division here by hw_2 is only for graph purposes [Kj/mol]
print(" Harmonic - Re-weighted Energy Fermions: ", avg_etot_f_h / hw_2)
















# plt.figure(1)
# plt.rcParams.update({'font.size': 13})
# n, bins, patches = plt.hist(pos, 100, density=True, facecolor='b', alpha=0.75)
# plt.xlabel('Position')
# plt.ylabel('Counts')
# plt.title('Position of Bead')
# plt.xlim([-10, 10])
# plt.grid(True)
# plt.show()























stop = time.time()
duration = stop - start
print("Time of execution:", duration)