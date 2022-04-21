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
# bhw_val = [3, 4, 6, 10, 30, 60]
bhw_val = [60]
# beads = [4, 8, 12, 16, 32, 64]
# beads = [4, 8, 12, 16, 32]
beads = [1, 8, 16, 32]
hw_ev = 26.73  # eV of hbaromega
k_spring = 0.84*m*(((hw_ev/1000)*conv_1/conv)/hbar)**2  # Hartree/Bohr^2  8.10542107E-7
print(k_spring)
# Harmonic Siumlation Constants



# number_of_files = 72
for bhw in bhw_val:
    for p in beads:
        temperature, beta = temp_for_lammps(hw_ev, bhw)
        path = "/home/netanelb2/Desktop/Netanel/Research/exiton/harmonic/boson3/bhw60/{}beads/".format(p)
        print ("Path: ", path)

        # time_step, pot, avg_pot_exp, stdv_pot, trap, trapvir, newvir = etot_b_2d_exitons(p, path)
        #
        # print("1) Energy Fermions / P: ", (2 * avg_pot_exp / (p*hw_exitons)), "+-", stdv_pot / hw_exitons)  # [Kj/mol]
        # print("2) Energy Fermions / P: ", np.mean(  ((trap / p) + newvir)  /  hw_exitons), "+-", stdv_pot / hw_exitons)  # [Kj/mol]


# ### bhw60 ### [10 M Iterations printfreq = 1000
# energy_array = np.array([0.10006391075648322, 0.77703294702903, 1.410437119600736,  2.193401746855274])
# stdv_energy_array = np.array([0.0003840310171667427, 0.009679779037786999, 0.019406191722290946, 0.04068149961357886])
# ### bhw60 ### ]

# ### bhw30 ### [10 M Iterations printfreq = 1000
energy_array = np.array([0.20009499411574647, 1.4125716306809173, 2.199205615167362,  2.7169993965385744])
stdv_energy_array = np.array([0.0006667631743926251, 0.011722411672378636, 0.036450812219989026,  0.051901870666212414])
# ### bhw30 ### ]

# ### bhw10 ### [10 M Iterations printfreq = 1000
# energy_array = np.array([1.8696828627901545, 2.553335327677357, 0,  2.8608938468767238, 2.970251897644674])
# stdv_energy_array = np.array([0.010971150570047529, 0.026133841124951283, 0, 0.045527554966455726, 0.08734176433400506])
# ### bhw10 ### ]

# ### bhw6 ### [10 M Iterations printfreq = 1000
# energy_array = np.array([2.41059346982529, 2.8164069323021854 , 2.927500169221537,  2.954762032805305, 2.993909363080621])
# stdv_energy_array = np.array([0.014428673054994368, 0.030519542018269803, 0.04222784361948804,0.05599016449365126,  0.11312014044199956])
# ### bhw6 ### ]

# ### bhw4 ### [10 M Iterations printfreq = 1000
# energy_array = np.array([2.7237398546906246, 2.9599930191533717 , 2.999988987051284,  3.0276072355708123, 3.063968997580119])
# stdv_energy_array = np.array([0.018623912765726017, 0.0357161550133291, 0.05617213434647265,0.072999610932241,  0.14732450595478416])
# ### bhw4 ### ]


## bhw3 ### [10M Iterations printfreq = 1000
# energy_array = np.array([3.1358126429287956, 3.314385368638948, 3.361266695984942,   3.366668146863846, 3.4146613864528774])
# stdv_energy_array = np.array([0.021532260953199937, 0.036715536215920365, 0.0609103765057502, 0.09042567249403877,  0.22697133791789711])
## bhw3 ### ]

plt.axhline(y=3.0, color='r', linestyle='-', label="Convergence (16,2.11)")
# plt.axvline(x = 12, color='r', linestyle='-')
plt.plot(beads, energy_array, 'o', color="blue")
plt.errorbar(beads, energy_array, yerr=stdv_energy_array, fmt="o", ecolor="black")
plt.title("bhw - Tot Energy vs Beads")
# plt.ylim([2,4])
plt.xlabel("beads")
plt.ylabel("Total Energy")
plt.legend(loc='lower right')
plt.show()

# Plot for Energy vs Beta
fig = plt.figure()
q = np.linspace(0.5, 11, 1000)
p = ana_3_bosons(q, 1)
plt.plot(q, p, color="black", label="Analytical Result")
plt.ylabel('<E>_B / hw')
plt.xlabel('bhw')
plt.legend()
plt.show()




















stop = time.time()
duration = stop - start
print("Time of execution:", duration)

#                                     ###       Plot    of Periodic Potential     ###
#
# X, Y = np.meshgrid(np.linspace(-33, 33, 512), np.linspace(-28.5, 28.5, 512))
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









