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

hbar = 1.05457182 * 10 ** -34   # J * s
mass_electron_kg = 9.10938e-31     # kg
omega = 6.088 * 10 ** -13 / hbar   # omega = 5.772959114344625e+21 1 /s (3.8MeV for exiting electon)

lamb = 1e13   #  J / m^4
lamb1 = 1.5e14   #  J / m^4
lamb2 = 0.5e14   #  J / m^4
mass_exciton = mass_electron_kg * 0.84     # kg  mass of exciton

pi = math.pi
psi = pi
hbaromega = 0
mass_e_amu = 0.00054858 * 0.84  # mass of single exiton
kB = 0.0083144621 #Boltzmann const in kJ/mol/K

# bhw_val = [3, 4, 6, 10, 30, 60]
bhw_val = [60]
beads = [1, 8, 16, 64, 80, 96]
# beads = [4, 8, 12, 16]
# beads = [4, 8, 12, 16, 64]
# beads = [1, 2, 4, 8, 12, 16]
# beads = [32]
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
        path = "/home/netanelb2/Desktop/Netanel/Research/exiton/anharmonic/boson3/bhw60/{}beads/".format(p)    # Anharmonic  # REMEMBER TO CHANGE N = 10 / 3
        # path = "/home/netanelb2/Desktop/Netanel/Research/exiton/harmonic/boson10/bhw30/{}beads/".format(p)   # Harmonic
        print("Path: ", path)

    #     time_step, pot, avg_pot_exp, stdv_pot, trap, stdv_trap,trapvir, newvir, l, p_l = etot_b_2d_exitons(p, path, beta)
    #     # print("1) Energy Fermions / P: ", (2 * avg_pot_exp / (p*hw_exitons)), "+-", stdv_pot / hw_exitons)  # [Kj/mol]
    #     EB = np.mean ((trap / p) + newvir)  # it is divided by hw in hartree + conversion
    #     EB_Err = stdv_pot / (p)
    #     Benergy[i] = EB / hw_exitons
    #     Bstdv[i] = EB_Err / hw_exitons
    #     print("Energy Bosons / P: ", EB / hw_exitons, "+-", EB_Err / hw_exitons)  # [Kj/mol]
    # print("Boson Energy:", list(Benergy))
    # print("Boson Energy STDV :", list(Bstdv))

##### Boson3 Harmonic
# energy_beta = np.array([5.6560275685720605, 3.147383429177774, 3.0648197524823098, 3.0149694389179937, 3.005403343416991, 2.965474462494367 , 2.8645714276422196])
# stdv_energy_beta = np.array([0.0023915946511268516, 0.0007634896205101604, 0.0018825207740906726, 0.0015105669121401432 ,0.001638165972034178, 0.0009084307395090711, 0.0006890692807012612])
# bhw_energy_beta = np.array([1, 3, 4, 6, 10, 30, 60])
###### Boson10 Harmonic
energy_beta = np.array([13.646026202596689, 10.2375407447492, 10.073433336077665, 10.062382145802628, 10.010431952912294, 9.89156989769864, 9.545800537485356])
stdv_energy_beta = np.array([0.010058492324181866, 0.003955722773832558, 0.0033504818644061484, 0.003901105319106849, 0.00295066843408104, 0.0017469206883813514, 0.0012489266779385198])
bhw_energy_beta = np.array([1, 3, 4, 6, 10, 30, 60])


tit = "bhw{} - Energy vs Beads".format(bhw_val[0])
# plt.axhline(y=3, color='r', linestyle='-', label="Convergence") # 3.36
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
fig = plt.figure()
# q = np.linspace(0.9, 61, 1000)
# p = ana_3_bosons(q, 1)
# plt.axhline(y=10, color='r', linestyle='-', label="Convergence")
# plt.plot(q, p, color="black", label="Analytical Result")
plt.plot(bhw_energy_beta, energy_beta, 'o', color="blue")
plt.errorbar(bhw_energy_beta, energy_beta, yerr=stdv_energy_beta, fmt="o", ecolor="black")
plt.title("Energy vs bhw - 10 Bosons")
plt.ylabel('<E>_B / hw')
plt.xlabel('bhw')
plt.legend()
plt.show()

###########################################  CONDENSATE  PLOT > #############################################
# fig = plt.figure(figsize=(10, 7))
# plt.bar(l, p_l / (1 / N))
# plt.xlabel('Permutation Length')
# plt.ylabel('Normalized Permutation Probability')
# plt.show()
###########################################  CONDENSATE  PLOT ^ #############################################

stop = time.time()
duration = stop - start
print("Time of execution:", duration)


#################################### Plot of Periodic Potential > ####################################
x_moire_lims = 35

X, Y = np.meshgrid(np.linspace(-x_moire_lims*nmtom, x_moire_lims*nmtom, 1024), np.linspace(-x_moire_lims*nmtom, x_moire_lims*nmtom, 1024)) # * 10 ** -9  # nm
Z = hbaromega + 2 * V_e * (np.cos((-2*pi/(math.sqrt(3)*moire_period))*X - (2*pi/moire_period)*Y - psi) +
                           np.cos((-2*pi/(math.sqrt(3)*moire_period))*X + (2*pi/moire_period)*Y - psi) +
                           np.cos(((4*pi)/(math.sqrt(3)*moire_period))*X - psi))
Z = Z/mevtoJ
# # Plot Moire Potential 2D
levels = np.linspace(Z.min(), Z.max(), 50)
fig, ax = plt.subplots()
plt.set_cmap('coolwarm')
graph = ax.contourf(X, Y, Z, levels=levels)
ax.set_title('Moire Potential (meV)')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.colorbar(graph)
plt.show()
#
#
# ########################################## Plot 1D Slice of 4 Potentials########################################
x_lims = 6      # harmonic
x_h_lims = 10      # Anharm
x_m_lims = 20   # Moire

fig = plt.figure()
x = np.linspace(-x_m_lims*nmtom, x_m_lims*nmtom, 100000)
x1 = np.linspace(-x_lims*nmtom, x_lims*nmtom, 100000)
y1 = np.linspace(-x_lims*nmtom, x_lims*nmtom, 100000)
x2 = np.linspace(-x_h_lims*nmtom, x_h_lims*nmtom, 100000)
y2 = np.linspace(-x_h_lims*nmtom, x_h_lims*nmtom, 100000)
# y = np.linspace(-30,30, 100000)
y = 0

## In Joul ###########################
# z = hbaromega + 2 * V_e * (np.cos((-2 * pi / (math.sqrt(3) * moire_period)) * x - (2 * pi / moire_period) * y - psi) +
#                          np.cos((-2 * pi / (math.sqrt(3) * moire_period)) * x + (2 * pi / moire_period) * y - psi) +
#                          np.cos(((4 * pi) / (math.sqrt(3) * moire_period)) * x - psi))
#
# z_harm = - 6 * V_e + ((16*np.pi**2*V_e) / (2*moire_period**2)) * x1**2   # 108 is the minimum of z
# z_anharm = - 6 * V_e + 0.25 * lamb * x2**4
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

#################################### Plot of Anharmonic potential > ####################################

# x_anharm = 5
# x_m_lims = 20   # Moire
#
# X, Y = np.meshgrid(np.linspace(-x_anharm*nmtom, x_anharm*nmtom, 2000), np.linspace(-x_anharm*nmtom, x_anharm*nmtom, 2000)) # [nm]
# # Z = - 6*V_e + 0.25 * lamb * (X**4 + Y**4)
# Z = - 6*V_e + 0.25 * lamb * ((X**2 + Y**2)*(X**2 + Y**2))
# Z = Z/mevtoJ  # [meV]
# print (Z.min())
#
# # # Plot Moire Potential 2D
# levels = np.linspace(Z.min(), Z.max(), 50)
# fig, ax = plt.subplots()
# plt.set_cmap('coolwarm')
# graph = ax.contourf(X, Y, Z, levels=levels)
# # ax.set_title('Anharmonic Potential (x^4+y^4 (meV))')
# ax.set_title('Anharmonic Potential ((x^2+y^2)^2 (meV))')
# plt.xlabel('x [m]')
# plt.ylabel('y [m]')
# plt.colorbar(graph)
# plt.show()
######################################################## ^ ###########################################################
########################################3D Contour #####################################################
# x_anharm = 19
# x_m_lims = 20   # Moire
# lamb1 = 1e12
# X, Y = np.meshgrid(np.linspace(-x_anharm*nmtom, x_anharm*nmtom, 500), np.linspace(-x_anharm*nmtom, x_anharm*nmtom, 500)) # [nm]
# # Z1 = - 6*V_e + 0.25 * lamb * (X**4+Y**4)
# # Z1 = Z1/mevtoJ  # [meV]
#
# # Z = - 6*V_e + 0.25 * lamb1 * ((X**2 + Y**2)*(X**2 + Y**2))
# Z = -6*V_e +  lamb1 * (X**4 + Y**4) + 2*lamb1 * (X**2*Y**2)
# Z_harm = -6*V_e + ((16*np.pi**2*V_e) / (2*moire_period**2)) * (X**2 + Y**2)
# # Z_moire = hbaromega + 2 * V_e * (np.cos((-2 * pi / (math.sqrt(3) * moire_period)) * X - (2 * pi / moire_period) * Y - psi) +
# #                          np.cos((-2 * pi / (math.sqrt(3) * moire_period)) * X + (2 * pi / moire_period) * Y - psi) +
# #                          np.cos(((4 * pi) / (math.sqrt(3) * moire_period)) * X - psi))
#
# Z = Z/mevtoJ  # [meV]
# Z_harm = Z_harm/mevtoJ  # [meV]
# # Z_moire = Z_moire/mevtoJ  # [meV]
#
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.contour3D(X, Y, Z, 200, cmap='winter')  # Anharmonic (Blue)
# ax.contour3D(X, Y, Z_harm, 500, cmap='autumn')  # Harmonic (Red)
# # ax.contour3D(X, Y, Z_moire, 50, cmap='binary')  # Moire
# ax.set_zlim([-109, 45])
# # ax.set_title('Anharmonic Potential (x^4+y^4 (meV))')
# ax.set_title('Blue - Anharmonic / Red - Harmonic')
# ax.set_xlabel('x [m]')
# ax.set_ylabel('y [m]')
# plt.legend(["h","k"])
# plt.show()
######################################################## ^ ###########################################################