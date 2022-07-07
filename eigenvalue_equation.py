import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import eigs
from scipy import sparse
import os

#
path = "/home/netanelb2/Desktop/Netanel/Research/exiton/harmonic_periodic/try1/log.lammps"

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

lamb = 2.893e14  # J / nm^4
lamb1 = 1.5e14   #  J / nm^4
lamb2 = 0.5e14   #  J / nm^4

a = -25 #float(input('enter lower limit of the domain: '))
b = 25  #float(input('enter upper limit of the domain: '))
N = 1001 #int(input('enter number of grid points: '))
potential = "harm_exc" #str(input('Which potential: '))  # harm, anharm, harm_exc
x = np.linspace(a, b, N) * 10 ** -9  # nano-meters
y = x
h = x[1]-x[0]


V_e = V_e * 1.60217656535 * 10 ** -22   # convert from meV to [J]
moire_period = moire_period * 10 ** -9  # convert from nm to [m]
mass_exciton = mass_electron * 0.84     # [kg] mass of exciton


def harmonic_potential(m, w, x):
    return 0.5 * m * w**2 * (x**2)


def harmonic_potential_2D(m, w, x, y):
    return 0.5 * m * w**2 * (x**2 + y**2)


def harmonic_potential_exciton(V_e, a_m , x):
    energy = -6 * V_e + (16 * math.pi**2 * V_e / (2 * a_m**2)) * (x ** 2)  # J
    return energy


def harmonic_potential_exciton_2d(V_e, a_m , x, y):
    return (16 * math.pi**2 * V_e / (2 * a_m**2)) * (x ** 2 + y ** 2)

def harmonic_dimless_2d(x, y):
    return  0.5 * (x ** 2 + y ** 2)



#### Energy hw(nx+ny+1) for harmonic exciton from article
def hbar_omega_nx_ny(nx, ny, V_e, hbar, mass, a_m):
    hw_article = np.sqrt((16*math.pi**2*V_e*hbar**2)/(mass*a_m**2))*(nx + ny + 1) # in J
    return hw_article


def anharmonic_potential(lamb, x):
    # -97.65480657 this is the number for "hbarw/2" since without -6V I obtain 10.34519343 for the first state
    energy = -6 * V_e + 0.25 * lamb * x**4
    return energy


def morse_potential(x):
    return de*(1-np.exp(-alpha*(x-re)))**2


def laplacian_1D(n):
    D = -2.0 * np.eye(n)
    for i in range(n-1):
        D[i, i+1] = D[i+1, i] = 1
    return D


def laplacian_2D(n):
    D1 = laplacian_1D(n)
    return np.kron(D1, np.eye(n)) + np.kron(np.eye(n), D1)


def read_file_pot_classic(filename, cut, s_rows, s_footer):
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
    kinetic_energy = df['c_kin'][cut:]
    tot_energy = potenital_energy + kinetic_energy

    return time_step, potenital_energy, kinetic_energy, tot_energy


# Perturbation Theory for Anharmonic potential
def perturbation_anharmonic():
    num = 5
    pert_energy = np.zeros(num)
    for n in range(0, num):
        E = hbar * omega * (lamb / 4) * (6*n** 2 + 6*n + 3)
        pert_energy[n] = E
    print("perturbation energy: ", pert_energy / (1.60217656535 * 10 ** -22), "meV")


def schrodinger_1D(num_eigen=4, potential=None, m=None, w=None, lamb=None, V_e=None, moireperiod=None):
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
                    V[i, j] = harmonic_potential(m, w, x[i+1])
                    H = -(hbar_1**2/(2*m*h**2))*T + V     # SI units    m, hbar = 1, 1
                elif potential == "anharm":
                    V[i, j] = anharmonic_potential(lamb, x[i + 1])
                elif potential == "harm_exc":
                    V[i, j] = harmonic_potential_exciton(V_e, moireperiod, x[i+1])
                else:
                    print("Choose a Potential")
            else:
                V[i, j] = 0

    H = -(hbar**2/(2*mass_exciton*h**2))*T + V     # SI units

    val, vec = np.linalg.eig(H)
    z = np.argsort(val)
    z = z[0:num_eigen]
    # energies = (val[z]/val[z][0]) # Energies normalized to the first state
    energies = (val[z])
    energies_difference = np.zeros(len(energies)-1)
    for i in range(0, len(energies)-1):
        energies_difference[i] = energies[i+1] - energies[i]

    # print("E of state", energies , "meV")
    # # print(energies)
    # print("dE", energies_difference , "meV")
    # # print(energies_difference )

    print("E of state", energies / (1.60217656535*10**-22), "meV")
    # print(energies)
    print("dE", energies_difference / (1.60217656535*10**-22), "meV")
    # print(energies_difference )
    return z, vec


z, vec = schrodinger_1D(num_eigen=4, potential=potential, m=m, w=w, lamb=lamb, V_e=V_e, moireperiod=moire_period)

# hw_article = hbar_omega_nx_ny(0, 0, V_e, hbar, mass_exciton, moire_period)
# print("hw_2", hw_article / (1.60217656535*10**-22), "meV")

# plt.figure(figsize=(12, 10))
# for i in range(len(z)):
#     y = []
#     y = np.append(y, vec[:, z[i]])
#     y = np.append(y, 0)
#     y = np.insert(y, 0, 0)
#     plt.plot(x, y, lw=3, label="{} ".format(i))
#     plt.xlabel('x', size=14)
#     plt.ylabel('$\psi$(x)',size=14)
# plt.legend()
# plt.title('normalized wavefunctions for a harmonic oscillator using finite difference method', size=14)
# plt.show()


##Create meshgrid for x and y
# based on the article: Solving 2D Time Independent Schrodinger Equation Using Numerical Method
def schrodinger_2d(N=None, L=None, potential=None, V_e=None, moireperiod=None):
    X, Y = np.meshgrid(np.linspace(-L/2, L/2, N, dtype=float), np.linspace(-L/2, L/2, N, dtype=float))
    #Create matrix
    diag = np.ones([N])
    diags = np.array([diag, -2*diag, diag])
    D = sparse.spdiags(diags, np.array([-1, 0, 1]), N, N)

    if potential == "harm2":
        V = harmonic_potential_2D(m, w, X, Y)
        T = - (hbar_1**2 / (2 * m)) * sparse.kronsum(D, D)
    elif potential == "harm_exc_2":
        V = harmonic_potential_exciton_2d(V_e, moireperiod, X, Y)
        T = - (hbar**2 / (2 * mass_exciton)) * sparse.kronsum(D, D)
    elif potential == "dimless":
        # X = X * (mass_exciton*omega)/hbar
        # Y = Y * (mass_exciton*omega)/hbar
        V = harmonic_dimless_2d(X, Y)
        T = -( 1 / 2 ) * sparse.kronsum(D, D)

    else:
        print("Choose a Potential")

    U = sparse.diags(V.reshape(N**2), (0))
    H = T + U
    #solve eigen vector and value
    eigenvalues, eigenvectors = eigsh(H, k=10, which='SM')
    e_v = np.zeros(len(eigenvalues))
    for i, v in enumerate(eigenvalues):
        e_v[i] = round(v/eigenvalues[0])

    return eigenvalues, eigenvectors, e_v

# llist = [1]  # harmonic exciton - harm_exc_2 - 0.5*10**-9
# for i in llist:
#     eigenvalues, eigenvectors, e_v = schrodinger_2d(N=100, L=i, potential="dimless", V_e=V_e, moireperiod=moire_period)
#     print("THIS is L", i)
#     print("EigenValues - rounded & normalized: ", e_v)
#     print("EigenValues normalized: ", eigenvalues / eigenvalues[0])
#     print("EigenValues in Joul: ", eigenvalues )
#     print("EigenValues in meV: ", eigenvalues / (1.60217656535*10**-22), "meV")



############## conservation of energy of a classical Simulation from log.lammps input ########

t, pot, kin, tot_e = read_file_pot_classic(path, 0, 66, 27)

plt.figure(figsize=(12, 10))
plt.plot(t, pot, color='b')
plt.plot(t, kin, color='r')
plt.plot(t, tot_e, color='k')
plt.xlabel('Time [fs]', size=14)
plt.ylabel('Total energy [Hartree]', size=14)
plt.legend()
plt.title('Conservation of Energy', size=14)
plt.show()

plt.figure(figsize=(12, 10))
plt.plot(t, (((tot_e-tot_e[0])*100)/tot_e[0]), color='k')
plt.xlabel('Time [fs]', size=14)
plt.ylim(-0.1, 0.1)
plt.ylabel('E(t)-E(0)/E(0)  [%]', size=14)
plt.legend()
plt.title('Energy Conservation', size=14)
plt.show()