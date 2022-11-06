# This code SHOULD provide a Winding_number analysis for B.E.C
# the code takes system_01.xyz and system_P.xyz files (trajectories of beads 0 and P-1 of each particle)
# for each timestep, for each particle - the code finds the minimum distance between beads (r(P)i-r(0)i or r(P)i-r(0)j)
# then, the code closes the "rings" according to the minimums
# later the code computes the winding number FOR the specific particle in EACH coordinate (x,y)
# finally, the code finds the average winding number out of all timesteps, and returns it

'''
# This code SHOULD provide a Winding_number analysis for B.E.C

# the code takes system_01.xyz and system_P.xyz files (trajectories of beads 0 and P-1 of each particle)
# for each timestep, for each particle - the code finds the minimum distance between beads (r(P)i-r(0)i or r(P)i-r(0)j)
# then, the code closes the "rings" according to the minimums
# later the code computes the winding number FOR the specific particle in EACH coordinate (x,y)
# finally, the code finds the average winding number out of all timesteps, and returns it
'''
# Imports
import networkx as nx
import numpy as np
import os
from tqdm import tnrange, tqdm, tqdm_notebook
from IPython.display import Image
from numpy import linalg as LA
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Functions:

def block_average(datastream, isplt=True, maxBlock=0):
    npoints = len(datastream)
    minBlock = 1;
    if (maxBlock == 0):
        maxBlock = int(npoints / 4)

    nblocks = maxBlock - minBlock

    mu_block = np.zeros(nblocks)
    std_block = np.zeros(nblocks)

    counter = 0
    for iblock in range(minBlock, maxBlock):
        ngrid_block = int(npoints / iblock)
        container = np.zeros(ngrid_block)
        for i in range(1, ngrid_block + 1):
            istart = (i - 1) * iblock
            iend = istart + iblock
            container[i - 1] = np.mean(datastream[istart:iend])

        mu_block[counter] = np.mean(container)
        std_block[counter] = np.std(container) / np.sqrt(ngrid_block - 1)
        counter += 1

    block = np.arange(minBlock, maxBlock)

    if isplt:
        plt.subplot(2, 1, 1)
        plt.plot(block, std_block, 'ro-', lw=1)
        plt.xlabel('Block')
        plt.ylabel('std')

        plt.subplot(2, 1, 2)
        plt.errorbar(block, mu_block, std_block)
        plt.xlabel('Block')
        plt.ylabel('mean')

        plt.show()
        print(len(mu_block))
        print('"<x> =" {} "+/-" {}'.format(mu_block, std_block))

    mean_mu = np.mean(mu_block)
    print(mean_mu)

    return block, mu_block, std_block, mean_mu

def to_graph(pair):
    g_connected =[]
    G = nx.Graph()
    G.add_edges_from(pair)
    g_connected=list(nx.connected_components(G))
    ring_len_list=[]

    for i in range(len(g_connected)):
        ring_len_list.append(len(g_connected[i]))

    # G.nodes
    # nx.draw(G, with_labels=True, font_weight='bold')

    return g_connected, ring_len_list


# a function that helps pandas remove rows with index 5i or 5i+1
def remove_nan(index, natoms):
    if natoms == 3:  # for N = 3
        if index % 5 == 0 or index % 5 == 1:  # 5 for 3 particles 12 for 10 particles
            return True
    else:          # for N = 10
        if index % 12 == 0 or index % 12 == 1:    # 5 for 3 particles 12 for 10 particles
            return True
    return False

#Constants
bhw = [1, 3, 4, 6, 10, 30, 60]

conv_m_to_bohr = 5.29177249*10**-11
m = 0.84*9.1093837e-31   # mass of exciton
hbar = 1.05457182e-34  #MKS
kb = 1.380649e-23   # MKS
kB = 0.0083144621 #Boltzmann const in kJ/mol/K
omega = 0.02673  # [eV]                                # only relevant for moire
omega = omega/27.2114  # [Hartree]                     # only relevant for moire
omega_kJmol = omega/3.8088E-4  # [kJmol]               # only relevant for moire


nstart = 270000
nstart = 1948714
# for bhw1 N=10  itr = 23624568

#only read first and last bead
nbeads=2

first_bead = '01'
last_bead = 32
positions_skipped = []
unitcells=[]
#nskip=0
#nstepend=16000
natoms = 10 # Number of Bosons in simulation

bhw = 1
T = 1 / (kB*(bhw / omega_kJmol))
print("Temperature:", T)
p = last_bead
path_ = "/home/netanelb2/Desktop/Netanel/Research/exiton/moire_one_{}p/".format(natoms)
# path_ = "/home/netanelb2/Desktop/Netanel/Research/exiton/kfir_results/harmonic/"
path = path_ + 'boson{}/bhw{}/{}beads/'.format(natoms, bhw, p)

# The only relevant files are 'system_01.xyz' and 'system_P.xyz' hence only 2 types of beads (1 and P)


files = [path+'system_{}.xyz'.format(k) for k in [first_bead, last_bead]]

#######################################READ FILE ###############################################
for inx, file in enumerate(files):
    length = 0
    print('File: ', file)
    if inx == 0:
        pos_1 = pd.read_table(file, skiprows=lambda x: remove_nan(x, natoms), delim_whitespace=True, names=['x', 'y', 'z']).to_numpy()
    elif inx == 1:
        pos_P = pd.read_table(file, skiprows=lambda x: remove_nan(x, natoms), delim_whitespace=True, names=['x', 'y', 'z']).to_numpy()
        length += len(pos_1)

positions = np.append(pos_1, pos_P)
positions = np.reshape(positions, (nbeads*length, 3))  # since there are 2 list that are appended

ndata=int(positions.shape[0]/nbeads/natoms)
natom_frame=int(ndata*natoms)
beads=np.zeros((nbeads, natom_frame, 3))

for i in range(nbeads):
    beads[i] = positions[i*natom_frame:(i+1)*natom_frame, :]
    # print(beads)
# print('nbeads:', nbeads, 'ndata:', ndata, 'natoms:', natoms)
beads = np.reshape(beads, (nbeads, ndata, natoms, 3))

# print(beads.shape)

############################################################################################

Lfac = [718.096]  # Lammps unit cell
cond_L = np.zeros(len(Lfac))
for q, L_size in enumerate(Lfac):
    print("L size: ", L_size)
    lx, ly, lz = L_size, L_size/2.30940107659, 1500
    # lx, ly, lz = 718.019, 310 .945, 1500  # 19nm in bohr
    L = np.array([lx, ly, lz])
    W2, dW = 0.0, 0.0
    ring_dist = []
    w2list = []
    wlistx = []
    wlisty = []
    for t in range(nstart, ndata):  # t goes frin nstart to all time steps
        pair = []
        Wtmp = np.zeros(3)
        if t % 10000 == 0:
            print("t", t)   # This is only to see progess of run.
        for i in range(natoms):  # natoms
            delW = beads[-1, t, :, :] - beads[0, t, i, :]    # -1 is for bead "P" and 0 is for bead "1"
            # need to consider the PBC!
            # x
            delx = delW[:, 0]
            delx = np.where(delx < L[0] / 2, delx, delx - L[0])
            delx = np.where(delx >= -L[0] / 2, delx, delx + L[0])
            delW[:, 0] = delx
            # y
            dely = delW[:, 1]
            dely = np.where(dely < L[1] / 2, dely, dely - L[1])
            dely = np.where(dely >= -L[1] / 2, dely, dely + L[1])
            delW[:, 1] = dely
            # z
            delz = delW[:, 2]
            delz = np.where(delz < L[2] / 2, delz, delz - L[2])
            delz = np.where(delz >= -L[2] / 2, delz, delz + L[2])
            delW[:, 2] = delz

            j = LA.norm(delW, axis=-1).argmin()
            # min1 = LA.norm(delW, axis=-1)[j]

            if (i != j):
                # print(i,j)

                pair += [(i, j)]
                dW = (beads[0, t, j, :] - beads[0, t, i, :])   # 0 means first bead and time t of atom j and i
                # dW = (beads[-1, t, i, :] - beads[0, t, i, :])
                #
                # PBC
                # x
                delx = dW[0]
                delx = np.where(delx < L[0] / 2, delx, delx - L[0])
                delx = np.where(delx >= -L[0] / 2, delx, delx + L[0])
                dW[0] = delx
                # y
                dely = dW[1]
                dely = np.where(dely < L[1] / 2, dely, dely - L[1])
                dely = np.where(dely >= - L[1] / 2, dely, dely + L[1])
                dW[1] = dely
                # z
                delz = dW[2]
                delz = np.where(delz < L[2] / 2, delz, delz - L[2])
                delz = np.where(delz >= - L[2] / 2, delz, delz + L[2])
                dW[2] = delz

                Wtmp += dW / L
            else:
                ring_dist += [1]


        # W2 += np.sum(Wtmp * Wtmp)
        W2 += Wtmp * Wtmp
        gconnected, ring_len_dist = to_graph(pair)
        ring_dist += ring_len_dist
        w2list.append(np.sum(Wtmp * Wtmp))
        wlistx.append(Wtmp[0])
        wlisty.append(Wtmp[1])
        pass
        # pdb.set_trace()  # for DEBUGGING

    print("W2", W2)
    print(W2 / (ndata - nstart))
    ww = W2 / (ndata - nstart)


    lambda_D = np.sqrt(2*np.pi*hbar**2/(m*kb*T)) / conv_m_to_bohr


    # cond = 2*np.pi*(lx/lambda_D)**2*(ww)/(3*natoms)   # 3D
    cond = 2*np.pi*np.sum((L/lambda_D)**2*(ww)/(2*natoms))     # 2D
    print("condensation frac ", cond)
    cond_L[q] = cond

    pass
print("Condesation as function of L: ", cond_L, "L: ", Lfac)



temp = [310, 103, 77, 51, 31, 10, 5]
#
# fig = plt.figure(figsize=(10, 7))
# tit = "Moire One-well - 3 Excitons"
# # plt.axhline(y=3, color='r', linestyle='-', label="Convergence") # 3.36
# # plt.axvline(x = 12, color='r', linestyle='-')
# plt.plot(temp, superfluid_frac, 'o', color="blue")
# # plt.plot(temp, ww1, 'o', color="red")
# plt.title(tit)
# plt.xlabel("temp [K]")
# plt.ylabel("Superfluid Fraction ")
# plt.legend(loc='lower right')
# plt.show()

fig = plt.figure(figsize=(10, 7))
tit = "Moire One-well - 3 Excitons - Temp {}".format(round(T,3))
# plt.axhline(y=3, color='r', linestyle='-', label="Convergence") # 3.36
# plt.axvline(x = 12, color='r', linestyle='-')
plt.plot(Lfac, cond_L, 'o', color="blue")
# plt.plot(temp, ww1, 'o', color="red")
plt.title(tit)
plt.xlabel("L [bohr]")
plt.ylabel("Superfluid Fraction ")
plt.legend(loc='lower right')
plt.show()

fig = plt.figure(figsize=(10, 7))
tit = "Moire One-well - 3 Excitons - Temp {}".format(round(T,3))
# plt.axhline(y=3, color='r', linestyle='-', label="Convergence") # 3.36
# plt.axvline(x = 12, color='r', linestyle='-')
# plt.plot(Lfac, cond_L100, 'o', color="blue")
# plt.plot(Lfac, cond_L77, 'o', color="red")
# plt.plot(Lfac, cond_L51, 'o', color="yellow")
# plt.plot(Lfac, cond_L61, 'o', color="orange")
# plt.plot(Lfac, cond_L10, 'o', color="green")
# plt.plot(temp, ww1, 'o', color="red")
plt.title(tit)
plt.xlabel("L [bohr]")
plt.ylabel("Superfluid Fraction ")
plt.legend(loc='lower right')
plt.show()