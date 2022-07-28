# # EXITON

from functions_PIMD import *
from scipy import optimize
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import math
from scipy.optimize import fmin
start = time.time()
#################################################### CONSTANTS ################################################
nmtom = 10**-9
mevtoJ = 1.60217656535 * 10 ** -22

moire_period = 19  # [nm] From Tran article
V_e = 18  # [meV] From Tran article
moire_period = moire_period * nmtom  # [m]
V_e = V_e * mevtoJ     # [J]
conv_1 = 96.4853  # eV to kJ/mol

hbar = 1.05457182 * 10 ** -34   # J * s
mass_electron_kg = 9.10938e-31     # kg
omega = 6.088 * 10 ** -13 / hbar   # omega = 5.772959114344625e+21 1 /s (3.8MeV for exiting electon)

lamb = 1.7978199e-12   #  Hartree / Bohr^4
lamb1 = 1.5e14   #  J / m^4
lamb2 = 0.5e14   #  J / m^4
mass_exciton = mass_electron_kg * 0.84     # kg  mass of exciton

pi = math.pi
psi = pi
hbaromega = 0
mass_e_amu = 0.00054858 * 0.84  # mass of single exiton
kB = 0.0083144621 #Boltzmann const in kJ/mol/K

# bhw_val = [1, 3, 4, 6, 10, 30, 60]
# beta_array = np.array([0.3877,  1.1632,  1.5510,  2.3264,  3.8774, 11.6322, 23.2644])
bhw_val = [10]
# beads = [1, 8, 16, 32, 64, 80, 96]
# beads = [4, 8, 12, 16]
beads = [4, 8, 12, 16, 64]
# beads = [1, 2, 4, 8, 12, 16]
# beads = [64]
hw_ev = 26.73  # eV of hbaromega
k_spring = 0.84*m*(((hw_ev/1000)*conv_1/conv)/hbar)**2  # Hartree/Bohr^2  8.10542107E-7
# print(k_spring)

Benergy= np.zeros(len(beads))
Bstdv = np.zeros(len(beads))
for bhw in bhw_val:
    for i, p in enumerate(beads):
        temperature, beta = temp_for_lammps(hw_ev, bhw)
        print ("temp : ", temperature)
        # path = "/home/netanelb2/Desktop/Netanel/Research/exiton/anharmonic/boson3/bhw10/{}beads/".format(p)    # Anharmonic  # REMEMBER TO CHANGE N = 10 / 3 in functions_PIMD
        path = "/home/netanelb2/Desktop/Netanel/Research/exiton/harmonic/boson3/bhw10/{}beads/".format(p)   # Harmonic
        print("Path: ", path)

        time_step, pot, avg_pot_exp, stdv_pot, trap, stdv_trap,trapvir, newvir, l, p_l = etot_b_2d_exitons(p, path, beta)  # [ kJ/mol ]
        EB = np.mean ((trap / p) + newvir)  # it is divided by hw in hartree + conversion
        EB_Err = stdv_pot / (p)
        # Benergy[i] = EB / hw_exitons
        Benergy[i] = EB * 1000 / conv_1  # [meV]
        Bstdv[i] = EB_Err * 1000 / conv_1  # [meV]
        print("Energy Bosons / P: ", EB * 1000 / conv_1, "+-", EB_Err * 1000 / conv_1)  # [meV]
    print("Boson Energy:", list(Benergy))
    print("Boson Energy STDV :", list(Bstdv))


# bhw_energy_beta = np.array([1, 3, 4, 6, 10, 30, 60])
# beta_array = np.array([0.3877,  1.1632,  1.5510,  2.3264,  3.8774, 11.6322, 23.2644])
# tit = "bhw{} - Energy vs Beads".format(bhw_val[0])
# # plt.axhline(y=3, color='r', linestyle='-', label="Convergence") # 3.36
# # plt.axvline(x = 12, color='r', linestyle='-')
# plt.plot(beads, Benergy, 'o', color="blue")
# plt.errorbar(beads, Benergy, yerr=Bstdv, fmt="o", ecolor="black")
# plt.title(tit)
# # plt.ylim([2,4])
# plt.xlabel("beads")
# plt.ylabel("Total Energy [meV] ")
# plt.legend(loc='lower right')
# plt.show()
#
# Plot for Energy vs Beta
# fig = plt.figure()
# # q = np.linspace(0.9, 61, 1000)
# # p = ana_3_bosons(q, 1)
# # plt.axhline(y=10, color='r', linestyle='-', label="Convergence")
# # plt.plot(q, p, color="black", label="Analytical Result")
# # plt.plot(bhw_energy_beta, energy_beta, 'o', color="blue")
# # plt.errorbar(bhw_energy_beta, energy_beta, yerr=stdv_energy_beta, fmt="o", ecolor="black")
# plt.plot(beta_array, energy_beta, 'o', color="blue")
# plt.errorbar(beta_array, energy_beta, yerr=stdv_energy_beta, fmt="o", ecolor="black")
# plt.title("Energy vs bhw - 10 Bosons")
# plt.ylabel('<E>_B / hw')
# plt.xlabel('bhw')
# plt.legend()
# plt.show()

###########################################  CONDENSATE  PLOT > #############################################
# fig = plt.figure(figsize=(10, 7))
# plt.bar(l, p_l / (1 / N))
# plt.xlabel('Permutation Length')
# plt.ylabel('Normalized Permutation Probability')
# plt.show()
###########################################  CONDENSATE  PLOT ^ #############################################

#################################### Plot of Periodic Potential > ################################################
x_moire_lims = 100
units = 'hartreebohr'  # or 'meVmeter'
# moire_plot(units, x_moire_lims, hbaromega, moire_period, V_e, nmtom, mevtoJ, psi, pi)
# ########################################## Plot 1D Slice of 4 Potentials########################################
x_harm_lims = 1.5 * 10 ** 13 / 6    # harmonic
x_an_lims = 1.5 * 10 ** 13   / 6  # Anharm
x_m_lims = 1.5 * 10 ** 13   # Moire

bohrconv = 5.2918e-11
V_e = 0.000661487
moire_period = 359.048
x_moire_lims = 1.5 * 10 ** 13
lamb=9.5*10**-11
x_h = np.linspace(-x_harm_lims * bohrconv, x_harm_lims * bohrconv, 5024)
x_an = np.linspace(-x_an_lims * bohrconv, x_an_lims * bohrconv, 5024)
x_m = np.linspace(-x_m_lims * bohrconv, x_m_lims * bohrconv, 5024)
y_m = 0
# four_potential_slice(x_h, x_an, x_m, y_m, V_e, moire_period, pi, psi, lamb, hbaromega)
########################################3D Contour of Harmonic and Anharmonic potential##############################

x_lim = 200
lamb1 = 9.5*10**-11    # Hartree / bohr^4    ( first LAMBDA =1.7978 *10**-12 H/bohr^4)
k = 8.10279e-07      # Hartree / Bohr^2
V = 0.0006614875     # Hartree
pi = math.pi
# three_d_contour_potential(x_lim, lamb1, k, V, pi)
######################################################## ^ ###########################################################


stop = time.time()
duration = stop - start
print("Time of execution:", duration)
