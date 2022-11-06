#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:12:40 2019
This code will allow analysing bosonic simulations.
@author: hirshb and updated by NBS
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from auxfunctions import CalcPhiEstimator, CalcPhiEstimator_from_PLUMED, CalcVirEstimator_from_PLUMED, \
    CalcUintEstimator_from_PLUMED, permutation_prob_3, calc_Wn,  permutation_prob_10
correct = False


#HARMONIC FORCE CONSTANT!
# omega = 0.003 #in eV
# omega = omega/27.2114 #eV to Ha
# omega_kJmol = omega/3.8088E-4 #Ha to kJmol

#Constants
kB = 0.0083144621 #Boltzmann const in kJ/mol/K
d = 2

#Changing parameters
# seeds = [98743501, 269451, 666472, 782943, 1239451]
seeds = [1239451]
Natoms = 10

# bhw_energy_beta = np.array([1, 3, 4, 6, 10, 30, 60])
# beta_array = np.array([0.3877,  1.1632,  1.5510,  2.3264,  3.8774, 11.6322, 23.2644])
# bhw_val = [10, 30, 60]
bhw_val = [30]

step_start = 10000  # line to start from colvar.
step_end = 300000  # 300000  # 272745 #

Nbeads = [4, 8, 12, 16, 32, 64]
Nbeads = [4, 8, 12, 16, 32]
Nbeads_ = [1, 8, 16, 32, 64, 80, 96]
# Nbeads = [32]
Nbeads_ = [96]


gs = [0]  # Interaction between particles


skiprows_in_Phi = 197  # 134 moire_onewell_3    197 moire_onewell_10   harmonic_10 = 209  #anharmonic 107  #153


#MOIRE FORCE CONSTANT!
omega = 0.02673  # [eV]
omega = omega/27.2114  # [Hartree]
omega_kJmol = omega/3.8088E-4  # [kJmol]

#For BSE, Min block size, max block size, step and whether to plot convergence vs. block size or not.
Nblocks = 5
maxsize = int(np.floor((step_end-step_start)/Nblocks))
minsize = maxsize
step = 1 #dummy
plot = False


count = 0
for g in gs:
    for ind, seed in enumerate(seeds):
        print("g: ", g, "seed: ", seed)
        for bhw in bhw_val:
            T = 1 / (kB*(bhw / omega_kJmol)) # [K]
            beta = 1/(kB*T)   # kJ/mol/K
            if T > 45:     # for higher temperatures dt is smaller hence more iterations  45 Kelvin
                step_end = 10*step_end
                print("number of steps*10 because high temp (more itr lower dt)")
            if T < 15:
                Nbeads = Nbeads_

            Benergy = np.zeros(len(Nbeads))  # To zero the energy for the next temp loop
            save_data = pd.DataFrame(columns=['Nbeads', 'g', 'seed', 'bhw', 'EB', 'l', 'p_l'])
            for i, p in enumerate(Nbeads):
                # if p == 16:
                #     step_end = 4213583
                # else:
                #     step_end = 1958664
                path_ = '/home/netanelb2/Desktop/Netanel/Research/exiton/moire_one_{}p/'.format(Natoms)
                # path_ = "/home/netanelb2/Desktop/Netanel/Research/exiton/kfir_results/harmonic/"
                # path_ = '/home/netanelb2/Desktop/Netanel/Research/exiton/dt05/harmonic/'
                path = path_+'boson{}/bhw{}/{}beads/'.format(Natoms, bhw, p)
                # path = "/home/netanelb2/Desktop/Netanel/Research/exiton/dt05/harmonic/boson10/bhw1/16beads/"

                print("path: ", path, "Temp [K]:", round(T, 3), "beta [kJ/mol/K]:", round(beta), "omega_kjmol: [kJ/mol]", round(omega_kJmol, 3))

                #Get Pot energy
                Trap = CalcPhiEstimator_from_PLUMED(path, 'log.lammps', p, step_start, step_end, skip=skiprows_in_Phi, potlabel='c_trap')
                # NBS: I obtain Trap when I add all the TRAP from all logfiles.
                Phi = CalcPhiEstimator_from_PLUMED(path, 'log.lammps', p, step_start, step_end, skip=skiprows_in_Phi) # Pot Estimator
                # Get kin energy
                Vir = CalcVirEstimator_from_PLUMED(path, 'log.lammps', p, step_start, step_end, skip=skiprows_in_Phi) # Kinetic Estimator

                Trap = Trap / 3.8088E-4  # Ha to kJmol
                Phi = Phi / 3.8088E-4  # Ha to kJmol
                Vir = Vir / 3.8088E-4  # Ha to kJmol

                #Get RP energies
                fname = "pimdb.log"

                try:
                    data = pd.read_csv(path + fname, sep='\s+', dtype='float64', header=None)
                except:
                    print("IOError: Problem reading file " + path+fname)
                    raise
                if p == max(Nbeads):  # Permutation Probability (Condensation)
                    if Natoms == 3:
                        l, p_l = permutation_prob_3(path + fname, beta, step_start, Natoms)
                    elif Natoms == 10:
                        l, p_l = permutation_prob_10(path+fname, beta, step_start, Natoms)
                else:
                    l, p_l = 0, 0

                # This below create title for pimdb.log file ("data")
                all_labels_e = []
                save_label_n = []
                for n in range(1, Natoms+1):
                    save_label_e = []
                    for k in range(n, 0, -1):
                        lis = [ '{:02d}'.format(int(x)) for x in range(n-k+1,n+1)]
                        label=''.join(map(str, lis))
                        #        print(label)
                        if k==1:
                            index = str(lis[0])
                        else:
                            index = str(lis[0]) + '-' + str(lis[-1])

                        save_label_e.append(label)
                        all_labels_e.append(label)
                    save_label_n.append(str(n))
                data.columns = ['E' + x for x in all_labels_e] + ['VB' + str(x) for x in range(Natoms+1)] # gives name E01    E0102      E02   E010203  ...  VB0       VB1      VB2      VB3

                data = data.iloc[step_start:step_end] # Reads the pimdb.log file from step_start untill step_end

                vir_EB = (Phi + Vir)
                EB3 = np.mean(vir_EB)  # [kJ/mol]
                print('Boson Energy: ' + str(EB3), "kJ/mol")
                save_data.loc[count, 'Nbeads'] = p

                save_data.loc[count,'seed'] = seed
                save_data.loc[count,'g'] = g
                save_data.loc[count,'bhw'] = bhw
                save_data.loc[count, 'l'] = l
                save_data.loc[count, 'p_l'] = p_l
                count += 1

                # convert from kJ/mol to meV
                conv_1 = 96.4853  # kJ/mol to eV
                Benergy = np.array(list(save_data['EB'])) * 1000 / conv_1  # meV
                save_data.loc[count, 'EB'] = Benergy  # meV

#################################### PLOTS ###############################################
            r_beta = round(beta)  # beta rounded
            np.savetxt(path_+'Figures/output_bosons{}bhw{}'.format(Natoms, bhw), save_data, fmt="%s", header='Nbeads g seed bhw EB l p_l', comments="")
            path_save_plot_betaEB = path_+'Figures/betaEB_bosons{}bhw{}Nbeads{}.png'.format(Natoms, bhw, Nbeads[-1])
            path_save_plot_cond = path_+'Figures/cond_bosons{}bhw{}Nbeads{}.png'.format(Natoms, bhw, Nbeads[-1])

            fig = plt.figure(figsize=(10, 7))
            tit = "Beta{} - Energy vs Beads".format(r_beta)
            # plt.axhline(y=3, color='r', linestyle='-', label="Convergence") # 3.36
            # plt.axvline(x = 12, color='r', linestyle='-')
            plt.plot(Nbeads, Benergy, 'o', color="blue")
            plt.title(tit)
            plt.xlabel("Number of Beads")
            plt.ylabel("Boson Energy [meV] ")
            plt.legend(loc='lower right')
            plt.savefig(path_save_plot_betaEB)
            # plt.show()

            # Plot for Condensation

            fig1 = plt.figure(figsize=(10, 7))
            plt.bar(l, p_l / (1 / Natoms))
            plt.xlabel('Permutation Length', fontsize=20)
            plt.ylabel('Normalized Permutation Probability', fontsize=20)
            # plt.ylim(0.99, 1.01)
            plt.savefig(path_save_plot_cond)
            # plt.show()

            print("All data:", save_data)  # ['g', 'seed', 'bhw', 'EB', 'l', 'p_l']
            print("Beads: ", Nbeads, "Boson Energy: [meV]", Benergy)

        # meanB = np.mean(save_data['EB'])
        # BSEB = np.std(save_data['EB'])/np.sqrt(len(save_data['EB']))

# Results

# pl1 = np.array([0.14931201, 0.12929172, 0.12086855, 0.11247752, 0.10071464,
#        0.09195319, 0.084679  , 0.07788191, 0.07016435, 0.0626571 ])
# pl3 = np.array([0.10060294, 0.1004373 , 0.10030757, 0.10023092, 0.10012034,
#        0.10000366, 0.09985377, 0.09966185, 0.09949627, 0.09928539])
# pl4 = np.array([0.10035573, 0.10030799, 0.10021925, 0.10016867, 0.10007645,
#        0.09998898, 0.0999106 , 0.09977306, 0.09962452, 0.09957475])
# pl6 = np.array([0.10019791, 0.10016776, 0.10013255, 0.10009298, 0.10006604,
#        0.09999076, 0.09993939, 0.09986757, 0.09977976, 0.09976527])
# pl10 = np.array([0.10002745, 0.10002405, 0.10002039, 0.10001356, 0.10000751,
#        0.10000231, 0.09999131, 0.09997929, 0.099967  , 0.09996712])
# pl60 = np.array([0.10000565, 0.10000537, 0.10000364, 0.10000253, 0.10000319,
#        0.10000024, 0.09999879, 0.09999473, 0.09999359, 0.09999228])
#
# fig1 = plt.figure(figsize=(10, 7))
# # tit = "Permutation Probability - Temp {} K".format(77)
# plt.bar(l, pl6 / (1 / 10))
# plt.title(tit,  fontsize=20)
# plt.xlabel('Permutation Length', fontsize=20)
# plt.ylabel('Normalized Permutation Probability', fontsize=20)
# plt.ylim(0.99, 1.01)
# plt.show()