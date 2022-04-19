from functions_PIMD import *
start = time.time()

# Harmonic Siumlation Constants
# kB = 0.0083144621 #Boltzmann const in kJ/mol/K     # Barak
# T = 34.8/bhw  # (bhw = 1.25)                       # Barak
# beta = 1/(kB*T)                                    # Barak
bhw = [1.25, 1.5, 1.75, 2.0, 2.5]
beta_list_5_v2 = [4.320127964728478, 5.184153557674173, 6.048179150619869, 6.912204743565564, 8.640255929456956]  # [1 / kJ/mol ]
temp = [27.8509, 23.2091, 19.8935, 17.4068, 13.9254]  # [K]
# Auxiliary and Bogoliubov Simulation Constants
bhw_g = [3, 4, 5, 6]
temp_g = [11.6, 8.7, 6.96, 5.8]  # [K]
beta_list_g = [10.368307115348346, 13.824409487131128, 17.28051185891391, 20.73661423069669]


# Path of Folder and how many Beads (number of files) (For Harmonic num_of_files is 12 for Auxiliary is 72)
number_of_files = 72
# path = '/home/netanelb2/Desktop/Netanel/Research/PIMD/runs/fermions/five2/bhw1p25/sim5/'  # Harmonic
path = '/home/netanelb2/Desktop/Netanel/Research/PIMD/runs/fermions/gauss/bogo/bhw3/sim1/'     # Auxiliary
file_pimdb = path + 'pimdb.log'
beta = 10.368307115348346

#######################################################################################################################

#                                 ###      3 Fermions under Harmonic potential run         ###

# time_step,  avg_etot_f_h, sign_array, wj = etot_f_2d_harmonic(number_of_files, path, beta)
#
# print("<S>_B: ", np.mean(sign_array))  # [Kj/mol]
# print("W_j: ", wj)  # [Kj/mol]
#  #division here by hw_2 is only for graph purposes [Kj/mol]
# print(" Harmonic - Re-weighted Energy Fermions: ", avg_etot_f_h / hw_2)

# #                               ###     3 Fermions under Auxiliary potential run          ###

# Re-weighted results: <(pot_est+kin_est)*s> / <s>  = <(pot/P + newvir)*s>_B / <s>_B
# Bogoliubov Results: E_H' - <pot - harmonic(trap)>_H' = E_H' - <V_g>_H'>  >=  E_H

time_step, sign_array, wj, avg_etot_f_a, avg_etot_f_b = etot_f_2d_harmonic_aux(number_of_files, path, beta)
print("<S>_B: ", np.mean(sign_array))  # [Kj/mol]
print("W_j: ", wj)  # [Kj/mol]
print("Auxiliary - Re-weighted Energy Fermions: ", avg_etot_f_a / hw_2)
print("Auxiliary - Bogoliubov: ",  avg_etot_f_b / (number_of_files * hw_2))
print("Re-weighted - Vg", (avg_etot_f_a - (avg_etot_f_b / number_of_files)) / hw_2)  # E_H' - <V_g>_H'







stop = time.time()
duration = stop - start
print("Time of execution:", duration)


















#                                     ### Results for bhw = {3,4,5,6} ### Auxiliary System and Bogoliubov
#bhw 3
wj_3 = np.array([210004.8593726606, 186185.55305947014, 213982.27978131443, 202573.24525250704, 183733.2605096042])
etot_f_3_a_array = np.array([5.676597176932612, 5.739254266235507, 5.678254283657318, 5.679397044460362, 5.742916118912327])
etot_f_3_b_array = np.array([5.383706690956485, 5.48518149484957, 5.384307947084825, 5.391548967415173, 5.490810443012252])
sign_avg_3 = np.array([0.2837903505035954, 0.2516020987290137, 0.2891652429477222, 0.2737476287196041, 0.24828818987784354])
#bhw 4
wj_4 = np.array([106580.16795465646, 117527.41639917549, 78411.4400749824, 88679.25590672411, 118002.68890790107])
etot_f_4_a_array = np.array([5.539805274644795, 5.524756050228152, 5.5854505238310885, 5.6412809860425135, 5.455452642488769])
etot_f_4_b_array = np.array([5.278729613216011, 5.239719413476682, 5.39723880811313, 5.457426175882458, 5.1456864815264])
sign_avg_4 = np.array([0.144027253992779, 0.15882083297185878, 0.10596140550673297, 0.11983683230638394, 0.15946309311878523])
#bhw 5
wj_5 = np.array([44076.777059272594, 75709.22267712229, 68376.33525409011, 33845.90619840927, 57058.63030484316])
etot_f_5_a_array = np.array([ 5.634939395614105, 5.40020354085211, 5.428324918794425, 5.715037549546594, 5.561358368004263])
etot_f_5_b_array = np.array([5.495353844147165, 5.086822357472333, 5.14568998659368, 5.697222685235016, 5.36024423474])
sign_avg_5 = np.array([0.059563212242260265, 0.10230976037448958, 0.09240045304606771, 0.045737711078931445,  0.07710625716870698])
#bhw6
wj_6 = np.array([57984.19628909118, 36641.15233128122, 50682.70690898466, 24550.539287341282, 54689.07166212774])
etot_f_6_a_array = np.array([5.227673814242167, 5.302847477041122, 5.308631963219681, 5.674839387092791, 5.2778622640190616])
etot_f_6_b_array = np.array([4.842119677836193, 5.007890288326537, 4.94528093720046, 5.546773077706502, 4.917727469932478])
sign_avg_6 = np.array([0.07835702201228538, 0.0495150707179476, 0.0684901444716009, 0.033176404442353084, 0.07390415089476722])

# Final Results gaussian and Bogo
sign_avg_a_b = np.array([0.26931870215555576, 0.137621883579308, 0.07542347878209119, 0.060688558507790834])
etot_f_array_a = np.array([5.701457314550066, 5.5414850955346395, 5.515302665885774, 5.3193271535079045])
stdv_f_array_a = np.array([0.015272990105292414, 0.03133796446517137, 0.05964174674504432, 0.06728098590472299])
etot_f_array_b = np.array([5.424137999200333, 5.288264040424372, 5.2957061054480885, 4.987910812140755])
stdv_f_array_b = np.array([0.024639556703329606,  0.05579659762332275, 0.10987127588158654, 0.10653525451071216])

# #                            try2         ### Results for bhw = {1.25,1.5,1.75,2,2.5} ### Only Harmonic
# #bhw1.25
# wj_1p25 = np.array([112864.16604383902, 117606.2297878852, 118604.76613394207, 115077.77326353025, 125637.66042745464])
# etot_f_1p25_array = np.array([6.823037821416585,  6.850215389410304, 6.795925346156208, 6.873700071067444, 6.700930320876907])
# sign_avg_1p25 = np.array([0.15251914330248517, 0.1589273375511962, 0.1589273375511962, 0.15551050441017603, 0.16978062219926301])
# #bhw1.5
# wj_1p5 = np.array([62734.91903474599, 63377.748534289414, 64940.63148946411, 67218.10356636555, 64170.46041154753])
# etot_f_1p5_array = np.array([6.313795800561833, 6.331024235908322, 6.223979873880679, 6.173050415572744, 6.229694759976889])
# sign_avg_1p5 = np.array([0.08477691761452161, 0.08564560612741813, 0.08775761012089744, 0.09083527508968317, 0.08671683839398316])
# #bhw1.75
# wj_1p75 = np.array([40004.92207657017, 39004.27577724805, 33770.499426682836, 39013.97948954744, 34774.6139840806])
# etot_f_1p75_array = np.array([5.828133724793112,  5.659618629414346, 6.012260612605722, 5.760196121018797, 5.936577617774265])
# sign_avg_1p75 = np.array([0.054060705508878606,  0.05270848078006493, 0.045635810036057885, 0.05272159390479384, 0.046992721600108915])
# #bhw2
# wj_2 = np.array([16967.33652565199, 25162.97481631888, 18542.542643326957, 26069.307393955667, 20747.484344292334])
# etot_f_2_array = np.array([5.9653755294430475, 5.012753871764524,  5.663783907526129, 5.098676598797878, 5.372140426371725])
# sign_avg_2 = np.array([0.02292883314277296, 0.03400402002205254, 0.02505749005854994, 0.03522879377561577, 0.028037141005800452])
# #bhw2.5
# wj_2p5 = np.array([5657.982104567379, 5787.995613428307, 7703.345914436861,  2699.0748836396215, 7482.558285767174])
# etot_f_2p5_array = np.array([4.9533534368382295, 4.638495539123845, 4.875233215145961, 7.755674085126569, 4.9699248959868605])
# sign_avg_2p5 = np.array([0.00764592176292889, 0.007821615693822035, 0.010409926911401163, 0.003647398491404894, 0.010111565251036722])

### TOTAL harmonic potential
sign_avg_5 = [0.159402852, 0.087146427678, 0.0504238, 0.0290510885, 0.010799586]
etot_f_array_5 = [6.806878494086993, 6.2530150950758525,  5.832238527064522, 5.365639185421901, 5.132805031301638]
stdv_f_array_5 = [0.030484528845417365, 0.029700493444146128,  0.06184914896734548, 0.17514134943286927, 0.4463006923704222]

sign_avg_5_v2 = np.array([0.15955686287327467, 0.08724413765519949, 0.05048910486412432, 0.029095192602032777, 0.007946735077753436])
etot_f_array_5_v2 = np.array([6.901289992143137, 6.383315530057243, 6.001863829786661, 5.711517982801404, 5.746343105765936])
stdv_f_array_5_v2 = np.array([ 0.02798230127157498, 0.029509423362973766, 0.06765394958776051, 0.20068033218394465, 0.5243678741776379])

#                            try1         ### Results for bhw = {1.25,1.5,1.75,2,2.5} ### Only Harmonic
# #bhw1.25
# wj_1p25 = np.array([112864.166, 117606.2297, 118604.766133, 115077.773263, 125637.66042])
# etot_f_1p25_array = np.array([6.8230378, 6.850215, 6.7959253, 6.8737000, 6.700930])
# sign_avg_1p25 = [0.15251914, 0.1589273, 0.1602767, 0.1555105,  0.16978062]
# #bhw1.5
# wj_1p5 = np.array([62734.91903, 63377.748534, 64940.631489, 67218.103566, 64170.46041])
# etot_f_1p5_array = np.array([6.3137958, 6.331024, 6.2239798, 6.173050, 6.2296947])
# sign_avg_1p5 = [0.0847769,  0.08564560, 0.0877576,  0.0908352, 0.08671683839]
# #bhw1.75
# wj_1p75 = np.array([40004.922076, 39004.275777, 33770.499426, 39013.97948, 34774.61398])
# etot_f_1p75_array = np.array([5.8281337,  5.65961862, 6.0122606, 5.7601961, 5.936577617])
# sign_avg_1p75 = [0.0540607, 0.0527084, 0.0456358, 0.0527215939, 0.0469927]
# #bhw2
# wj_2 = np.array([16967.3365, 25162.974816, 18542.54264, 26069.3073939, 20747.4843442])
# etot_f_2_array = np.array([5.9653755, 5.012753, 5.6637839, 5.09867, 5.3721404])
# sign_avg_2 = [0.022928, 0.03400402, 0.02505749, 0.035228793, 0.02803714]
# #bhw2.5
# # wj_2p5 = np.array([8465.32979, 7525.056499, 9785.34025, 7888.28519, 6295.1803964])
# # etot_f_2p5_array = np.array([5.0299618, 4.789912, 4.5190857, 4.64861188, 5.05696040])
# # sign_avg_2p5 = [0.01143963, 0.0101689, 0.0132234, 0.010659, 0.008507]
# wj_2p5 = np.array([5657.98210, 5787.99561, 7703.345914, 2699.074883, 7482.55828])
# etot_f_2p5_array = np.array([4.95335, 4.63849, 4.8752332, 7.755674, 4.96992])
# sign_avg_2p5 = [0.0076459, 0.00782, 0.0104099, 0.003647398, 0.0101115]

#                                      Results for first 5 points (Only harmonic potential)
sign_avg_5_try1 = [0.159402852, 0.087146427678, 0.0504238, 0.0290510885, 0.010799586]
etot_f_array_5_try1 = [6.8068783, 6.2530149,  5.8322385, 5.365637, 5.132802]
stdv_f_array_5_try1 = [0.0304845, 0.0297005, 0.06184914, 0.17514208, 0.44630122]


#                                        ############ Barak Data ##############

# #bhw1.25
# wj_1p25_barak = np.array([103264.19756, 101075.215348, 101163.337928, 102849.908, 102390.9363])
# etot_f_1p25_array_barak = np.array([6.91262412, 6.940551, 6.9571200, 6.91561865, 6.96245745])
# sign_avg_1p25_barak = [0.158873, 0.15550595, 0.1556294, 0.15823888, 0.157520]
# #bhw1.5
# wj_1p5_barak = np.array([57812.3495, 58005.62921257, 57173.91370, 58093.27251, 57899.4352763])
# etot_f_1p5_array_barak = np.array([6.42738909, 6.4094888, 6.44295, 6.40783615, 6.387069032])
# sign_avg_1p5_barak = [0.0889360, 0.089240, 0.087956, 0.089373, 0.0890702]
# #bhw1.75
# wj_1p75_barak = np.array([29442.8146,29441.3204, 31422.4637, 31655.970031, 33391.6106127])
# etot_f_1p75_array_barak = np.array([6.22525393, 6.1482325, 6.1021057, 6.01349610, 5.909392])
# sign_avg_1p75_barak = [0.04529363,0.045295148, 0.048347,0.0486994, 0.051378489]
# #bhw2
# wj_2_barak = np.array([15657.3754801, 18275.584461, 15986.0143289, 17739.5295075, 18435.9103375])
# etot_f_2_array_barak = np.array([6.12199098, 5.8776883, 6.09188909, 5.74251589, 5.72096018])
# sign_avg_2_barak = [0.02408663, 0.0281190291, 0.024596491, 0.02728800, 0.02836]
# #bhw2.5
# wj_2p5_barak = np.array([6002.020729958988, 5643.323090707206,  4485.547580785226, 2915.760388031396, 7095.686301233578])
# etot_f_2p5_array_barak = np.array([5.385968571985375, 5.5266261470974625, 6.069622557267915, 7.493737880005943, 4.873965362086076])
# sign_avg_2p5_barak = [0.00922770815, 0.0086894861, 0.006901265, 0.004491869, 0.0109132339]
# #bhw3
# wj_3_barak = np.array([184977.3512652572, 184831.61767780365, 157468.9732848415, 185986.329866, 151628.11917911714])
# etot_f_3_array_barak = np.array([5.705093350017439, 5.707361641287298, 5.743720552496348, 5.660455943308668, 5.82292526590855])
# sign_avg_3_barak = [0.2845750381907103, 0.284363136987,0.242257207287930, 0.286142708, 0.233271209729]
# #bhw4
# wj_4_barak = np.array([85245.914543502, 89021.0652735, 81207.909253, 106117.4592])
# etot_f_4_array_barak = np.array([5.50345780765321, 5.542641647, 5.59843929, 5.34207670])
# sign_avg_4_barak = [0.14208271154, 0.1483646660, 0.135355980, 0.17686890]
# #bhw5
# wj_5_barak = np.array([75218.23284305, 36657.1608508, 46238.832556, 87446.944005])
# etot_f_5_array_barak = np.array([5.324355655, 5.677672714, 5.55642997736, 5.32250848731])
# sign_avg_5_barak = [0.1157166338038, 0.056403758, 0.07114490724, 0.1345301797]
# #bhw6
# wj_6_barak = np.array([34339.488792, 34339.355887, 30027.815102, 22171.530207, 32965.38598])
# etot_f_6_array_barak = np.array([5.32455984926, 5.3287294439, 5.2844063732, 5.3733643311, 5.3366618738])
# sign_avg_6_barak = [0.0686894870, 0.0472559, 0.060066226, 0.044333521, 0.040470853]

# Barak results for 3 fermions - Harmonic potential and Auxiliary.
# 5 is for first "5" points and "g" for last 4 points in fig6
etot_f_array_5_barak = [6.93755744, 6.4148657, 6.074709, 5.9004776433, 5.6297522]
stdv_f_array_5_barak = [0.01030199, 0.0094601, 0.05510046, 0.083725350, 0.40240035]
sign_f_array_5_barak = [0.15715344, 0.08891504, 0.047802733, 0.026490030, 0.00804471243]
etot_f_array_g_barak = [5.723669675071802, 5.48707494, 5.4201400039, 5.32728000]
stdv_f_array_g_barak = [0.02653834, 0.0575122, 0.086961, 0.013275099]
sign_f_array_g_barak = [0.26612186003, 0.1506680643, 0.09444886968, 0.0521631974]
# Unitl here barak resuls



#      Plot of the Fermionic Energy as well as the Average Sign to replicate Fig6
fig, (ax1, ax2) = plt.subplots(2)
q = np.linspace(0.8, 11, 1000)
p = ana_3_fermions(q, 1)
ax1.plot(q, p, 'g',  label="Analytical Result")
ax1.errorbar(bhw, etot_f_array_5_v2, yerr=stdv_f_array_5_v2, fmt='^', color='orange', label='H_original')
ax1.errorbar(bhw_g, etot_f_array_a, yerr=stdv_f_array_a, fmt='H', color='grey', label='H_aux')
ax1.errorbar(bhw, etot_f_array_5_barak, yerr=stdv_f_array_5_barak, fmt='x', color='green', label='Baraks results')
ax1.errorbar(bhw_g, etot_f_array_b, yerr=stdv_f_array_b, fmt='X', color='purple', label='Bogoliubov')
ax1.errorbar(bhw_g, etot_f_array_g_barak, yerr=stdv_f_array_g_barak, fmt='x', color='green')
ax1.set_ylabel('<E>_F')
ax1.set_ylim([1, 8])
ax1.set_xlim([1, 6.2])
ax2.plot(bhw, sign_avg_5_v2, 'D', color='red')
ax2.plot(bhw, sign_f_array_5_barak, 'D', color='green')
ax2.plot(bhw_g, sign_avg_a_b, 'D',  color='red', label='<s>_B')
ax2.plot(bhw_g, sign_f_array_g_barak, 'D',  color='green', label='Baraks - <s>_B')
ax2.set_ylabel('<S>_B')
ax2.set_xlabel('bhw')
ax2.set_ylim([-0.01, 0.28])
ax2.set_xlim([1, 6.2])
ax1.legend()
ax2.legend()
plt.show()
















# To calcualte E avgerage and stdv according to fermion supporting info
# avg_energy, energy_error = statistical_error_estimation(etot_f_1p25_array, wj_1p25)
# print("<E_F>: ", avg_energy, "+-", energy_error)

#                                                         ### Bosons run ###

# time_step, pot, avg_pot_exp, stdv_pot = etot_b_2d_harmonic(number_of_files, path)
# print("Total Energy Fermions / P: ", (2 * avg_pot_exp / (number_of_files * hw)))  # [Kj/mol]


                                                         # Results #
# Results from Barak's artilce:
# Natoms = 3
# bhws = [1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 6]
# ##bhws = [3, 4, 5, 6]
# analytic = [6.89333779, 6.36220084, 6.01473798, 5.77315087, 5.46533286, 5.28532944, 5.10755903, 5.04008246, 5.01482427]
# g,#seed#,EF,Err_EF,EB,Err_EB,sign,Err_sign,neff_reW,bhw,Wj,EF_corr,Err_EF_corr
# 0.0,#98743501.0#,5.385968571985375,2.1132931542097952,3.216831456368679,0.0008102108362961339,0.009227708151542006,0.0034864803536090398,94.60365877312594,2.5,6002.020729958988,5.385968571985375,2.1132931542097952
# 0.0,#269451.0#,5.5266261470974625,1.4976033641370217,3.2127540136006094,0.001769823082885048,0.00868948612488037,0.0023017135189826454,83.85152951973014,2.5,5643.323090707206,5.5266261470974625,1.4976033641370217
# 0.0,#666472.0#,6.069622557267915,1.5428141445133405,3.2132450833359405,0.001427607025719007,0.006901265636481244,0.0014369515996920433,52.915969065364244,2.5,4485.547580785226,6.069622557267915,1.5428141445133405
# 0.0,#782943.0#,7.493737880005943,3.5688294211287,3.2141862512354002,0.001736622325001616,0.004491869374804708,0.0018347522174338323,22.364885294766964,2.5,2915.760388031396,7.493737880005943,3.5688294211287
# 0.0,#1239451.0#,4.873965362086076,1.6017255669153663,3.2146077683841368,0.0016959293541149814,0.010913233936965382,0.0033463418722573953,132.43399266019878,2.5,7095.686301233578,4.873965362086076,1.6017255669153663

# Bosons 2 Quantum Particles

# Total Energy vs number of Beads - FOR BHW = 3 find converagance number of beads (All data below is for 2M itr)

# beads_array = np.array([2, 4, 8, 16, 32, 64, 72])
# energy_array = np.array([1.7395061436968315, 1.9979284032744322, 2.0681892824853216, 2.1298897342298404, 2.0976698115269805, 2.124903807277683, 2.0974433474243304])
# stdv_energy_array = np.array([0.014384479212790001, 0.012366058535679911, 0.010588211817062711, 0.01120683054549245, 0.010863020424765689, 0.012543236748566056, 0.014649020291975198])
# plt.axhline(y=2.11, color='r', linestyle='-', label="Convergence (16,2.11)")
# plt.axvline(x = 16, color='r', linestyle='-')
# plt.plot(beads_array, energy_array, 'o', color="blue")
# plt.errorbar(beads_array, energy_array, yerr=stdv_energy_array, fmt="o", ecolor="black")
# plt.title("bhw=3 - Tot Energy vs Beads")
# plt.xlabel("beads")
# plt.ylabel("Total Energy")
# plt.legend(loc='lower right')
# plt.show()


# Fig7 Replication BOSONS partilces

# w = 1
# bhw = [2, 3, 4, 5, 6]
# energy_array_hw = [2.4096314412021793, 2.0974433474243304, 2.0392088383062266, 2.008478942194873, 1.9993265561758102]
# stdv_energy_array = [0.024837541563550877, 0.014649020291975198, 0.0034618631722831725, 0.014645367890929661, 0.009282708812675358]
# figtotenergy = plt.figure()
# plt.rcParams.update({'font.size': 13})
# q = np.linspace(0.2, 11, 1000)
# p = part * hbar * w * (np.exp(hbar*w*q*4)+np.exp(hbar*w*q*3)+4*np.exp(hbar*w*q*2)+np.exp(hbar*w*q)+1)/(np.exp(hbar*w*q*4)-1)
# plt.plot(q, p, 'g',  label="Analytical Result")
# plt.plot(bhw, energy_array_hw, '.', label="Simulation", color="blue")
# plt.errorbar(bhw, energy_array_hw, yerr=stdv_energy_array, ecolor="black", fmt='o', markersize=3)
# plt.title("Total Energy vs bhw")
# plt.xlabel("bhw")
# plt.ylabel("<E>/hw")
# plt.legend(loc='upper right')
# plt.show()

# Fig7 Replication Distinguishable partilces

# w=1
# bhw = [2, 3, 4, 5, 6]
# energy_array_hw = [2.616543754337255, 2.197529044650856, 2.072392746224948, 2.005743796543706, 1.9996729203280086]
# stdv_energy_array = [0.014703183110991508, 0.016797341507523517, 0.003752243273667098, 0.009754321389733773, 0.008868449634375707]
# figtotenergy = plt.figure()
# plt.rcParams.update({'font.size': 13})
# q = np.linspace(0.2, 11, 1000)
# p = part * hbar * w * (1 + ((2 * np.exp(- hbar * w * q))/(1 - np.exp(- hbar * w * q))))
# plt.plot(q, p, 'g',  label="Analytical Result")
# plt.plot(bhw, energy_array_hw, '.', label="Mean Total Energy", color="blue")
# plt.errorbar(bhw, energy_array_hw, yerr=stdv_energy_array, ecolor="black", fmt='o', markersize=3)
# plt.xlabel("bhw")
# plt.ylabel("<E>/hw")
# plt.legend(loc='upper right')
# plt.show()

# Figure 7 - BOTH DISTINGUISHABLE and BOSONS results - b for bosons   d for distinguishable

# bhw = [2, 3, 4, 5, 6]
# w=1
# figtotenergy = plt.figure()
# plt.rcParams.update({'font.size': 13})
# energy_array_hw_b = [2.4096314412021793, 2.0974433474243304, 2.0392088383062266, 2.008478942194873, 1.9993265561758102]
# stdv_energy_array_b = [0.024837541563550877, 0.014649020291975198, 0.0034618631722831725, 0.014645367890929661, 0.009282708812675358]
# energy_array_hw_d = [2.616543754337255, 2.197529044650856, 2.072392746224948, 2.005743796543706, 1.9996729203280086]
# stdv_energy_array_d = [0.014703183110991508, 0.016797341507523517, 0.003752243273667098, 0.009754321389733773, 0.008868449634375707]
# q = np.linspace(0.2, 11, 1000)
# p = part * hbar * w * (np.exp(hbar*w*q*4)+np.exp(hbar*w*q*3)+4*np.exp(hbar*w*q*2)+np.exp(hbar*w*q)+1)/(np.exp(hbar*w*q*4)-1)
# plt.plot(q, p, 'k',  label="Bosons - Analytical Result")
# plt.plot(bhw, energy_array_hw_b, '.', label="Boson - Simulation", color="blue")
# plt.errorbar(bhw, energy_array_hw_b, yerr=stdv_energy_array_b, ecolor="blue", fmt='o', markersize=3)
#
# x1 = np.linspace(0.2, 11, 1000)
# y1 = part * hbar * w * (1 + ((2 * np.exp(- hbar * w * x1))/(1 - np.exp(- hbar * w * x1))))
# plt.plot(x1, y1, 'r',  label="Distinguishable - Analytical Result", dashes=[6, 2])
# plt.plot(bhw, energy_array_hw_d, '.', label="Distinguishable - Simulation", color="green")
# plt.errorbar(bhw, energy_array_hw_d, yerr=stdv_energy_array_d, ecolor="red", fmt='o', markersize=3, color="green")
# plt.xlabel("bhw")
# plt.ylabel("<E>/hw")
# plt.legend(loc='upper right')
# plt.show()

# Block Averaging - STDV of Data vs Block Size

# block_size_array = np.linspace(1, 1000, 200).astype(int)
# avg_array = np.zeros(len(block_size_array))
# stdv_array = np.zeros(len(block_size_array))
# number_of_blocks_array = np.zeros(len(block_size_array))
#
# for i, value in enumerate(block_size_array):
#     number_of_blocks, avg, stdv = block_averaging(5000, block_size=value, data=pot)
#     avg_array[i] = avg
#     stdv_array[i] = stdv
#     number_of_blocks_array[i] = number_of_blocks
#
# figblock = plt.figure()
# plt.plot(block_size_array, stdv_array, label="stdv", color="red")
# plt.xlabel("Block Size")
# plt.ylabel("STDV")
# plt.legend()
# plt.show()
#
# number_of_blocks1, avg, stdv = block_averaging(cutoff, block_size=165, data=pot)
# print ("number of blocks", number_of_blocks1)

