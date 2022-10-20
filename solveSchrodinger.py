import numpy as np
import matplotlib.pyplot as plt
# This codes solves the finite schrodinger equation for any potential given. Here the Harmonic and An harmonic
# potentials are written.
hbar = 1             # Hartree
m = 0.84             # amu
L = 120              # bohr
n = 51
pot = 'harm'

def laplacian_1D(n):
    D = -2.0 * np.eye(n)
    for i in range(n-1):
        D[i, i+1] = D[i+1, i] = 1
    return D


def laplacian_2D(n):
    D1 = laplacian_1D(n)
    return np.kron(D1, np.eye(n)) + np.kron(np.eye(n), D1)


def two_d_schrodinger(L, n, s, pot):
    '''
    :param L: length of X, Y coord
    :param n: dimension of matrix
    :param s: number of Energy state
    :return: eigenvalues , eigenvector
    '''
    # Coordinated
    x = np.linspace(-L, L, n)
    dx = x[1] - x[0]
    X = np.kron(x, np.ones(n))
    Y = np.kron(np.ones(n), x)
    # Energies
    T = -hbar ** 2 / (2 * m) / (dx ** 2) * laplacian_2D(n)  # Kinetic Energy
    if pot == 'harm':
        omega = 0.00098215 # 1/s- This number is obtained from sqrt(k/m) where k = 8.1045*10^-7 and m = 0.84 (26.73 meV)
        V = np.diag( 0.5 * m * omega**2 * (X**2 + Y**2))        # Potential
    elif pot == 'anharm':
        l = 9.5*10**-11        # Hartree / bohr**4    - 1.7978 *10**-12   - 9.5*10**-11
        V = np.diag(l * (X**4 + Y**4) + 2*l*((X**2)*(Y**2)))
    H = T + V                                               # Hamiltonian
    # Eigenvalues and Vectors
    E, U = np.linalg.eigh(H)                                # Eigenvalue/vector
    psi = U[:, s].reshape(n, n)

    return E, psi


E, psi = two_d_schrodinger(L, n, 0, pot)
E1 = E[0:7]
energies_difference = np.zeros(6)
for i in range(0, 6):
    energies_difference[i] = E1[i+1] - E1[i]

print(E[0:15])
print(energies_difference)

plt.xlabel(r'$x_1$')
plt.ylabel(r'$x_2$')
shw = plt.imshow(psi[:, ::-1], interpolation='bilinear', extent=[-L, L, -L, L])
bar = plt.colorbar(shw)
plt.show()

