import matplotlib.pyplot as plt
import numpy as np
import math
#
#This code solves 1D Schrodinger Equation for any 1D Potential
# This work quite fine for harmonic potential with a = -6 b = 6 and N = 1001
# Constant - Only to see results
m, w, hbar_1 = 1, 1, 1
de, re, alpha = 1, 1, 1
########################
hbar = 1.05457182 * 10 ** -34   # J * s
mass_electron = 9.10938e-31     # kg
omega = 6.088 * 10 ** -13 / hbar   # omega = 5.772959114344625e+21 1 /s (3.8MeV for exiting electon)
moire_period = 19 # [nm] From Tran article
V_e = 18  # [meV] From Tran article

# lamb = 2.893e14  # J / nm^4
lamb1 = 1.5e14   #  J / m^4
lamb2 = 0.5e14   #  J / m^4

a = float(input('enter lower limit of the domain: '))
b = float(input('enter upper limit of the domain: '))
N = int(input('enter number of grid points: '))
potential = str(input('Which potential: '))  # harm, anharm, harm_exc
x = np.linspace(a, b, N) * 10 ** -9
h = x[1]-x[0]


V_e = V_e * 1.60217656535 * 10 ** -22   # convert from meV to J
moire_period = moire_period * 10 ** -9  # convert from nm to m
mass_exciton = mass_electron * 0.84     # kg  mass of exciton


def harmonic_potential(m, w, x):
    return 0.5 * m * w**2 * x**2


def harmonic_potential_exciton(V_e, a_m , x):
    energy = -6 * V_e + (16 * math.pi**2 * V_e / (2 * a_m**2)) * x ** 2  # J
    return energy


def anharmonic_potential(lamb, x):
    # -97.65480657 this is the number for "hbarw/2" since without -6V I obtain 10.34519343 for the first state
    energy = -6 * V_e + 0.25 * lamb * x**4
    return energy


def morse_potential(x):
    return de*(1-np.exp(-alpha*(x-re)))**2

def schrodinger_1D(num_eigen=4, potential=None, mass=None, omega=None, lamb=None, V_e=None, moireperiod=None):
    T = np.zeros((N-2)**2).reshape(N-2, N-2)
    for i in range(N-2):
        for j in range(N-2):
            if i==j:
                T[i, j] = -2
            elif np.abs(i-j) == 1:
                T[i, j] = 1
            else:
                T[i, j] = 0


    V = np.zeros((N-2)**2).reshape(N-2, N-2)
    for i in range(N-2):
        for j in range(N-2):
            if i==j:
                if potential == "harm":
                    V[i, j] = harmonic_potential(mass, omega, x[i+1])
                elif potential == "anharm":
                    V[i, j] = anharmonic_potential(lamb, x[i + 1])
                elif potential == "harm_exc":
                    V[i, j] = harmonic_potential_exciton(V_e, moireperiod, x[i+1])
                else:
                    print("Choose a Potential")
            else:
                V[i, j] = 0

    H = -(hbar**2/(2*mass_exciton*h**2))*T + V

    val, vec = np.linalg.eig(H)
    z = np.argsort(val)
    z = z[0:num_eigen]
    # energies = (val[z]/val[z][0]) # Energies normalized to the first state
    energies = (val[z])
    energies_difference = np.zeros(len(energies)-1)
    for i in range(0, len(energies)-1):
        energies_difference[i] = energies[i+1] - energies[i]

    print("E of state", energies / (1.60217656535*10**-22), "meV")
    # print(energies)
    print("dE", energies_difference / (1.60217656535*10**-22), "meV")
    # print(energies_difference )
    return z, vec


z, vec = schrodinger_1D(num_eigen=5, potential=potential, mass=m, omega=w, lamb=lamb, V_e=V_e, moireperiod=moire_period)

plt.figure(figsize=(12, 10))
for i in range(len(z)):
    y = []
    y = np.append(y, vec[:, z[i]])
    y = np.append(y, 0)
    y = np.insert(y, 0, 0)
    plt.plot(x, y, lw=3, label="{} ".format(i))
    plt.xlabel('x', size=14)
    plt.ylabel('$\psi$(x)',size=14)
plt.legend()
plt.title('normalized wavefunctions for a harmonic oscillator using finite difference method', size=14)
plt.show()


#### Energy hw(nx+ny+1) for harmonic exciton from article
def hbar_omega_nx_ny(nx, ny, V_e, hbar, mass, a_m):
    hw_article = np.sqrt((16*math.pi**2*V_e*hbar**2)/(mass*a_m**2))*(nx + ny + 1) # in J
    return hw_article


hw_article = hbar_omega_nx_ny(0, 0, V_e, hbar, mass_exciton, moire_period)
print("hw_2", hw_article / (1.60217656535*10**-22), "meV")

## Perturbation Theory for Anharmonic potential
# num = 5
# pert_energy = np.zeros(num)
# for n in range(0, num):
#     E = hbar * omega * (lamb / 4) * (6*n** 2 + 6*n + 3)
#     pert_energy[n] = E
# print("perturbation energy: ", pert_energy / (1.60217656535 * 10 ** -22), "meV")


