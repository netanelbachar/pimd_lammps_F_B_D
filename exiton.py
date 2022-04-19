# # EXITON
from functions_PIMD import *
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import fmin
start = time.time()


# Moire Potential =  2V sum Cos (bj * r - psi )
# Constants
v = 18  # meV
a = 19  # nm
pi = math.pi
psi = pi
hbaromega = 0
mass = 0.00054858 * 0.84  # mass of single exiton
kB = 0.0083144621 #Boltzmann const in kJ/mol/K
bhw_val = [3, 4, 6, 10, 30, 60]
number_of_files = 72
hw_ev = 26.73  # eV of hbaromega

# Harmonic Siumlation Constants

temp = []
betas =[]
for bhw in bhw_val:
    temperature, beta = temp_for_lammps(hw_ev, bhw)
    path = "/home/netanelb2/Desktop/Netanel/Research/PIMD/runs/fermions/five2/bhw1p25/"

    time_step, avg_etot_f_h, sign_array, wj = etot_f_2d_harmonic(number_of_files, path, beta)

    print("<S>_B: ", np.mean(sign_array))  # [Kj/mol]
    print("W_j: ", wj)  # [Kj/mol]
    # division here by hw_2 is only for graph purposes [Kj/mol]
    print(" Harmonic - Re-weighted Energy Fermions: ", avg_etot_f_h / hw_2)





















stop = time.time()
duration = stop - start
print("Time of execution:", duration)

                                    ###       Plot    of Periodic Potential     ###

# X, Y = np.meshgrid(np.linspace(-25, 25, 512), np.linspace(-25, 25, 512))
# Z = hbaromega + 2* v * (np.cos((-2*pi/(math.sqrt(3)*a))*X - (2*pi/a)*Y - psi) +
#                            np.cos((-2*pi/(math.sqrt(3)*a))*X + (2*pi/a)*Y - psi) +
#                            np.cos(((4*pi)/(math.sqrt(3)*a))*X - psi))
# # Plot Moire Potential 2D
# levels = np.linspace(Z.min(), Z.max(), 50)
# fig, ax = plt.subplots()
# plt.set_cmap('coolwarm')
# graph = ax.contourf(X, Y, Z, levels=levels)
# ax.set_title('Moire Potential (meV)')
# plt.colorbar(graph)
# plt.show()
#
# # Plot Moire Potential 1D Slicr
# fig = plt.figure()
# x = np.linspace(-30, 30, 512)
# y = 10
# z = hbaromega + 2 * v * (   np.cos((-2*pi/(math.sqrt(3)*a))*x - (2*pi/a)*y - psi) +
#                             np.cos((-2*pi/(math.sqrt(3)*a))*x + (2*pi/a)*y - psi) +
#                             np.cos(((4*pi)/(math.sqrt(3)*a))*x - psi)   )
# plt.title('')
# plt.xlabel('X axis')
# plt.ylabel('Z axis')
# plt.plot(x, z, color="black")
# plt.title("Slice of Moire Potential at y = 10")
# plt.legend(loc="upper left")
# plt.show()









