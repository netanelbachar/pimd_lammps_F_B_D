import networkx as nx
import numpy as np
import os
from tqdm import tnrange, tqdm, tqdm_notebook
from IPython.display import Image
from numpy import linalg as LA
import pdb
import matplotlib.pyplot as plt
#Constants




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



T_list =[ 1.8]

conv_m_to_bohr = 5.29177249*10**-11
m = 6.6464731e-27
hbar = 1.05457182e-34
kb = 1.380649e-23
# T = 1.6

W2 = 0.0
ring_dist = []
w2list = []
wlist = []
probw = []
W2_prob = []
nstart = 10000

lx , ly, lz = 14.318839687622253, 14.318839687622253, 14.318839687622253
L = np.array([lx, ly, lz])


#only read first and last bead
nbeads=2
path ="/home/netanelb2/Desktop/Netanel/Research/exiton/superfluid/"
filepath = [path+'t{}/dump/'.format(k) for i, k in enumerate(T_list)]
positions_skipped=[]
unitcells=[]
#nskip=0
#nstepend=16000
natoms=64

#######################################READ FILE ###############################################
cond_array = np.zeros(len(filepath))
for q, file in enumerate(filepath):
    cond = 0
    for count, filename in tqdm(enumerate(sorted(os.listdir(file)))):
        src=filepath[q]+filename
        with open(src,"r") as f:
            #store each bead info.
            #print(src)
            positions=[]
            for lines in f:
                line=lines.split()
                if(len(line)==3):
                    unitcells.append([float(u) for u in line])
                if(len(unitcells)%3==0):
                    unitcell=[]
                    for u in unitcells[-3:]:
                        unitcell+=[u[1]-u[0]]
                    #print(unitcell)
                if(len(line)==5):
                    position=[float(pos) for pos in line]
                    position=position[2:]
                    #print(position)
                    positions.append([a*b for a,b in zip(position,unitcell)])
            positions_skipped+=positions[:]

    atoms=np.asarray(positions_skipped)
    atoms.shape
    #atoms=np.reshape(atoms, (-1,3))

    ndata=int(atoms.shape[0]/nbeads/natoms)
    natom_frame=int(ndata*natoms)
    beads=np.zeros((nbeads,natom_frame,3))
    for i in range(nbeads):
        beads[i]=atoms[i*natom_frame:(i+1)*natom_frame,:]
    print('nbeads:', nbeads, 'ndata:', ndata, 'natoms:', natoms)
    beads=np.reshape(beads,(nbeads,ndata,natoms,3))  # ARRANGE FILE

    for t in range(nstart, ndata):  # ndata
        pair = []
        Wtmp = np.zeros(3)

        # print("t", t)
        for i in range(natoms):  # natoms
            delW = beads[0, t, :, :] - beads[-1, t, i, :]
            # need to consider the PBC!
            # x
            delx = delW[:, 0]
            delx = np.where(delx < lx / 2, delx, delx - lx)
            delx = np.where(delx > -lx / 2, delx, delx + lx)
            delW[:, 0] = delx
            # y
            dely = delW[:, 1]
            dely = np.where(dely < ly / 2, dely, dely - ly)
            dely = np.where(dely > -ly / 2, dely, dely + ly)
            delW[:, 1] = dely
            # z
            delz = delW[:, 2]
            delz = np.where(delz < lz / 2, delz, delz - lz)
            delz = np.where(delz > -lz / 2, delz, delz + lz)
            delW[:, 2] = delz

            j = LA.norm(delW, axis=-1).argmin()
            # min1 = LA.norm(delW, axis=-1)[j]
            # tmptmp = np.delete(LA.norm(delW, axis=-1), j)
            # jj = tmptmp.argmin()
            # min2 = LA.norm(delW, axis=-1)[jj]

            # print(i,j,jj)
            # print(delW)
            if (i != j):
                # print(i,j)
                # j=jj
                pair += [(i, j)]
                dW = (beads[0, t, j, :] - beads[0, t, i, :])   # 0 means first bead and time t of atom j and i
                # dW = (beads[0, t, i, :] - beads[-1, t, i, :])

                # PBC
                # x
                delx = dW[0]
                delx = np.where(delx < lx / 2, delx, delx - lx)
                delx = np.where(delx > -lx / 2, delx, delx + lx)
                dW[0] = delx
                # y
                dely = dW[1]
                dely = np.where(dely < ly / 2, dely, dely - ly)
                dely = np.where(dely > -ly / 2, dely, dely + ly)
                dW[1] = dely
                # z
                delz = dW[2]
                delz = np.where(delz < lz / 2, delz, delz - lz)
                delz = np.where(delz > -lz / 2, delz, delz + lz)
                dW[2] = delz

                # print(dW)
                Wtmp += dW / L
            else:
                ring_dist += [1]



        W2 += np.sum(Wtmp * Wtmp)
        gconnected, ring_len_dist = to_graph(pair)
        ring_dist += ring_len_dist
        pass
        # w2list.append(np.sum(Wtmp * Wtmp))
        # wlist.append(Wtmp[0])
        # pdb.set_trace()  # for DEBUGGING

    print("W2", W2)
    print(W2 / (ndata - nstart))
    ww = W2 / (ndata - nstart)

    temp=np.array([1.0,1.6,1.8,2.0,2.4, 3.0])


    lambda_D_kfir = np.sqrt(2*np.pi*hbar**2/(m*kb*T_list[q])) / conv_m_to_bohr
    lambda_D=15.413053134916042/np.sqrt(T_list[q])

    cond = 2*np.pi*(lx/lambda_D)**2*(ww)/(3*natoms)

    print (cond)
    cond_array[q] = cond

print(cond_array)















###### Recreated result of Chang woo's He 4 using his code



# Rho is the recreated result by NBS
# rho=np.array([1.0681773371920247,0.7814648339266727, 0.37933121092481636, 0.13475290148913563, 0.0])
#
# temp        =np.array([1.00176179982669,1.597057080646505,1.9988541977687397,2.4120079258140006,2.994164552830728])
# temp_e      =np.array([0.0016712484331259305,0.005740745617858237,0.0083819039724898,0.0079194179181487,0.00799851670935463])
# winding     =np.array([1.07569455,0.8086872,0.3806735,0.13496873,0.0])
# winding_e   =np.array([1.12649116-1.07569455,0.84857185-0.8086872,0.38904519-0.38,0.13521991-0.13496873,0.0])
#
# temp_pimc        =np.array([1.18,1.6, 1.82,2.0, 2.35,2.5,2.86])
# winding_pimc     =np.array([0.98,0.88,0.66,0.44,0.16,0.08,0.03])
#
# fig, ax1 = plt.subplots()
#
# ax1.set_title('Superfluid fraction of He-4', fontsize=18)
#
# ax1.tick_params(axis="y", labelsize=14, direction='in')
# ax1.tick_params(axis="x", labelsize=14, direction='in')
#
# ax1.set_xlabel('Temperature (K)', fontsize=16)
# ax1.set_ylabel(r'$\rho_s\:/\:\rho$', fontsize=16)
#
# ax1.set_xlim(0.9,3.1)
# ax1.set_ylim(-0.1,1.2)
# ax1.errorbar(temp, winding, yerr=winding_e, fmt='ro--', ecolor='black', capsize=5,markersize=8,
#             markeredgewidth=1.0, markeredgecolor='k',label="Boson PIMD")
# # ax1.plot(temp_pimc,winding_pimc,'co--',lw=1,markersize=8,
# #          markeredgewidth=1.0, markeredgecolor='k',label="Pollock & Ceperley (1987) - PIMC")
# ax1.plot(temp,rho,'co--',lw=1,markersize=8, markeredgewidth=1.0, markeredgecolor='k',label="Boson PIMD - Recreate")
#
# ax1.legend()
























