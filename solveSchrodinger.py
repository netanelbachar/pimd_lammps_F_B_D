import numpy as np
import matplotlib.pyplot as plt

def laplacian_1D(n):
    D = -2.0 * np.eye(n)
    for i in range(n-1):
        D[i, i+1] = D[i+1, i] = 1
    return D

def laplacian_2D(n):
    D1 = laplacian_1D(n)
    return np.kron(D1, np.eye(n)) + np.kron(np.eye(n), D1)

hbar = 1
m = 0.84
L = 50
omega = 0.00095548246

n = 101

x = np.linspace(-L, L, n)
dx = x[1] - x[0]
X1 = np.kron(x, np.ones(n))
X2 = np.kron(np.ones(n), x)

T = -hbar ** 2 / (2 * m) / (dx ** 2) * laplacian_2D(n)
# T = -hbar ** 2 / 2 * laplacian_2D(n)
V = np.diag( 0.5 * m * omega**2 * ( X1**2 + X2**2 ) )

H = T + V

E, U = np.linalg.eigh(H)

psi0 = U[:,0].reshape(n, n)
plt.xlabel(r'$x_1$')
plt.ylabel(r'$x_2$')
shw = plt.imshow(psi0[:,::-1], interpolation='bilinear', extent=[-L, L, -L, L])
bar = plt.colorbar(shw)

plt.show();

print(E[0:10])