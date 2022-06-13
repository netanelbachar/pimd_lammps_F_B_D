# # EXITON
from functions_PIMD import *
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import fmin
start = time.time()
#################################################### CONSTANTS ################################################
nmtom = 10**-9
mevtoJ = 1.60217656535 * 10 ** -22

moire_period = 19  # [nm] From Tran article
V_e = 18  # [meV] From Tran article
moire_period = moire_period * nmtom
V_e = V_e * mevtoJ

hbar = 1.05457182 * 10 ** -34   # J * s
mass_electron_kg = 9.10938e-31     # kg
omega = 6.088 * 10 ** -13 / hbar   # omega = 5.772959114344625e+21 1 /s (3.8MeV for exiting electon)

lamb = 2.893e14   #  J / m^4
lamb1 = 1.5e14   #  J / m^4
lamb2 = 0.5e14   #  J / m^4
mass_exciton = mass_electron_kg * 0.84     # kg  mass of exciton

pi = math.pi
psi = pi
hbaromega = 0
mass_e_amu = 0.00054858 * 0.84  # mass of single exiton
kB = 0.0083144621 #Boltzmann const in kJ/mol/K
N = 3
# bhw_val = [3, 4, 6, 10, 30, 60]
bhw_val = [3]
# beads = [1, 8, 16, 32, 80, 96]
beads = [4, 8, 12, 16, 32, 64]
# beads = [16]
hw_ev = 26.73  # eV of hbaromega
k_spring = 0.84*m*(((hw_ev/1000)*conv_1/conv)/hbar)**2  # Hartree/Bohr^2  8.10542107E-7
print(k_spring)
## bhw =1 (310.186, 0.3877393605270402) (temp, beta)

# number_of_files = 72
Benergy= np.zeros(len(beads))
Bstdv = np.zeros(len(beads))
for bhw in bhw_val:
    for i, p in enumerate(beads):
        temperature, beta = temp_for_lammps(hw_ev, bhw)
        print ("temp : ", temperature)
        path = "/home/netanelb2/Desktop/Netanel/Research/exiton/harmonic/boson3/bhw3/{}beads/".format(p)
        print("Path: ", path)

        time_step, pot, avg_pot_exp, stdv_pot, trap, stdv_trap,trapvir, newvir, l, p_l = etot_b_2d_exitons(p, path, beta)
        # print("1) Energy Fermions / P: ", (2 * avg_pot_exp / (p*hw_exitons)), "+-", stdv_pot / hw_exitons)  # [Kj/mol]
        EB = np.mean(  ((trap / p) + newvir)  /  hw_exitons)
        EB_Err = stdv_pot / (p*hw_exitons)
        Benergy[i] = EB
        Bstdv[i] = EB_Err
        print("Energy Bosons / P: ", EB, "+-", EB_Err)  # [Kj/mol]
    print("Boson Energy:", list(Benergy))
    print("Boson Energy STDV :", list(Bstdv))

###### Boson3
# energy_beta = np.array([5.688106496436944, 3.178857302590968, 3.063968997580119, 2.993909363080621, 2.970251897644674, 2.9643799268188094, 2.8668079996664977])
# stdv_energy_beta = np.array([0.0029970987701447133, 0.0010822574671667015, 0.004603890811087005, 0.003535004388812486, 0.002729430135437658, 0.0015357912207428484, 0.001220797219627618])
# bhw_energy_beta = np.array([1, 3, 4, 6, 10, 30, 60])
###### Boson10
# energy_beta = np.array([0, 0, 0, 0, 0, 0])
# stdv_energy_beta = np.array([0, 0, 0, 0, 0, 0])
# bhw_energy_beta = np.array([0, 0, 0, 0, 0, 0])

#
tit = "bhw{} - Energy vs Beads".format(bhw_val[0])
plt.axhline(y=3.2, color='r', linestyle='-', label="Convergence") # 3.36
# plt.axvline(x = 12, color='r', linestyle='-')
plt.plot(beads, Benergy, 'o', color="blue")
plt.errorbar(beads, Benergy, yerr=Bstdv, fmt="o", ecolor="black")
plt.title(tit)
# plt.ylim([2,4])
plt.xlabel("beads")
plt.ylabel("Total Energy")
plt.legend(loc='lower right')
plt.show()

# Plot for Energy vs Beta
# fig = plt.figure()
# q = np.linspace(0.9, 61, 1000)
# p = ana_3_bosons(q, 1)
# plt.plot(q, p, color="black", label="Analytical Result")
# plt.plot(bhw_energy_beta, energy_beta, 'o', color="blue")
# plt.errorbar(bhw_energy_beta, energy_beta, yerr=stdv_energy_beta, fmt="o", ecolor="black")
# plt.title("Energy vs bhw - 3 Bosons")
# plt.ylabel('<E>_B / hw')
# plt.xlabel('bhw')
# plt.legend()
# plt.show()

###########################################  CONDENSATE  > #############################################
# fig = plt.figure(figsize=(10, 7))
# plt.bar(l, p_l / (1 / N))
# plt.xlabel('Permutation Length')
# plt.ylabel('Normalized Permutation Probability')
# plt.show()
###########################################  CONDENSATE  ^ #############################################

stop = time.time()
duration = stop - start
print("Time of execution:", duration)


#################################### Plot of Periodic Potential > ####################################
# x_moire_lims = 35
#
# X, Y = np.meshgrid(np.linspace(-x_moire_lims*nmtom, x_moire_lims*nmtom, 1024), np.linspace(-x_moire_lims*nmtom, x_moire_lims*nmtom, 1024)) # * 10 ** -9  # nm
# Z = hbaromega + 2 * V_e * (np.cos((-2*pi/(math.sqrt(3)*moire_period))*X - (2*pi/moire_period)*Y - psi) +
#                            np.cos((-2*pi/(math.sqrt(3)*moire_period))*X + (2*pi/moire_period)*Y - psi) +
#                            np.cos(((4*pi)/(math.sqrt(3)*moire_period))*X - psi))
# Z = Z/mevtoJ
# # # Plot Moire Potential 2D
# levels = np.linspace(Z.min(), Z.max(), 50)
# fig, ax = plt.subplots()
# plt.set_cmap('coolwarm')
# graph = ax.contourf(X, Y, Z, levels=levels)
# ax.set_title('Moire Potential (meV)')
# plt.xlabel('x [m]')
# plt.ylabel('y [m]')
# plt.colorbar(graph)
# plt.show()
#
#
# ########################################## Plot 1D Slice of 4 Potentials########################################
# x_h_lims = 5      # Anharm with same difference in G.S energy as Harmonic
# x_lims = 7      # An harm similar to harmonic
# x_m_lims = 20   # Moire
# fig = plt.figure()
# x = np.linspace(-x_m_lims*nmtom, x_m_lims*nmtom, 100000)
# x1 = np.linspace(-x_lims*nmtom, x_lims*nmtom, 100000)
# x2 = np.linspace(-x_h_lims*nmtom, x_h_lims*nmtom, 100000)
# # y = np.linspace(-30,30, 100000)
# y = 0
#
# ## In Joul ###########################
# z = hbaromega + 2 * V_e * (np.cos((-2 * pi / (math.sqrt(3) * moire_period)) * x - (2 * pi / moire_period) * y - psi) +
#                          np.cos((-2 * pi / (math.sqrt(3) * moire_period)) * x + (2 * pi / moire_period) * y - psi) +
#                          np.cos(((4 * pi) / (math.sqrt(3) * moire_period)) * x - psi))
# z_harm = ((16*np.pi**2*V_e) / (2*moire_period**2)) * x1**2 - 108*mevtoJ  # 108 is the minimum of z
# z_anharm = 0.25 * lamb * x2**4 - 108*mevtoJ
# z_anharm1 = 0.25 * lamb1 * x1**4 - 108*mevtoJ
# z_anharm2 = 0.25 * lamb2 * x1**4 - 108*mevtoJ
#
# ######################################
# ## In meV ############################
# z=z/mevtoJ
# z_harm = z_harm/mevtoJ
# z_anharm = z_anharm/mevtoJ
# z_anharm1 = z_anharm1/mevtoJ
# z_anharm2 = z_anharm2/mevtoJ
# #####################################
#
# plt.title('')
# plt.xlabel('position at x axis [m]')
# plt.ylabel('energy [meV]')
# plt.plot(x, z, color="black", label="Moire Potential at y = 0")
# plt.plot(x1, z_harm, color="blue", label="Harmonic Potential ")
# plt.plot(x1, z_anharm, color="purple", label="Anharmonic Potential".format(lamb))
# # plt.plot(x1, z_anharm, color="purple", label="Anharmonic $\lambda$ = {:.2e} J/nm^4 ".format(lamb))
# # plt.plot(x1, z_anharm1, color="red", label="Anharmonic1 $\lambda$ = {:.2e} J/nm^4 ".format(lamb1))
# # plt.plot(x1, z_anharm2, color="orange", label="Anharmonic2 $\lambda$ = {:.2e} J/nm^4".format(lamb2))
# # plt.title("Slice of Moire Potential at y = 0")
# plt.legend(loc="lower right")
# plt.show()

######################################################## ^ ###########################################################