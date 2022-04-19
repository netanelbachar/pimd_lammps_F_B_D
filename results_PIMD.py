import numpy as np
from functions_PIMD import statistical_error_estimation


#bhw1.25
wj_1p25 = np.array([112972.63583981787,  117720.93886400807,118717.85161879819, 115192.1142230634,125756.85208542869])
etot_f_1p25_array = np.array([6.917397164559292, 6.935955205243089,6.892726323388188, 6.965293944217862, 6.8038274590554195])
sign_avg_1p25 = np.array([0.152665724107862, 0.15908234981622713, 0.16042952921459214, 0.15566501922035594, 0.16994169200733608])
#bhw1.5
wj_1p5 = np.array([62805.32736347416, 63451.087888912254,65012.87271373798,  67291.74204346964, 64242.279314644125])
etot_f_1p5_array = np.array([6.4566064408197095, 6.452530450469957, 6.334465776710926, 6.31658790472761, 6.362632255608484])
sign_avg_1p5 = np.array([0.08487206400469481,  0.08574471336339494,  0.08785523339694322, 0.09093478654522924, 0.0868138909657353])
#bhw1.75
wj_1p75 = np.array([ 40054.281755032374, 39052.8282878969, 33817.84457254546, 39064.10674372506,  34820.62663806022])
etot_f_1p75_array = np.array([5.973790662082808, 5.826473291662583, 6.195400586017178, 5.948898113344316, 6.13109881947882])
sign_avg_1p75 = np.array([0.05412740777707078, 0.05277409228094175, 0.045699789962899265, 0.0527893334374663,  0.04705490086224354])
#bhw2
wj_2 = np.array([16997.20435468659,25199.035401843823, 18571.19735856291, 26105.963743996865,  20778.811768431075])
etot_f_2_array = np.array([6.264032185310658, 5.188877473066963,  5.930085107873297, 5.284286838386023,5.62762867611058])
sign_avg_2 = np.array([ 0.022969195073900796, 0.034052750543032195, 0.025096212646706636,  0.03527832938377955, 0.028079475362744695])
#bhw2.5
wj_2p5 = np.array([5671.401659373476,  5802.093433770617, 7718.613276236949, 2713.267116107566, 7497.544302199108])
etot_f_2p5_array = np.array([5.52155457455794, 5.450467678932403, 5.336345227914472, 8.849354668036758, 5.444494772734039])
sign_avg_2p5 = np.array([0.007664056296450643, 0.007840666802392725, 0.010430558481401283, 0.003666577183929143, 0.01013181662459339])

sign_avg_5_v2 = np.array([0.15955686287327467, 0.08724413765519949, 0.05048910486412432, 0.029095192602032777, 0.007946735077753436])
etot_f_array_5_v2 = np.array([6.901289992143137, 6.383315530057243, 6.001863829786661, 5.711517982801404, 5.746343105765936])
stdv_f_array_5_v2 = np.array([ 0.02798230127157498, 0.029509423362973766, 0.06765394958776051, 0.20068033218394465, 0.5243678741776379])


e1p25, s_e1p25 = statistical_error_estimation(etot_f_1p25_array, wj_1p25)
e1p5, s_e1p5 = statistical_error_estimation(etot_f_1p5_array, wj_1p5)
e1p75, s_e1p75 = statistical_error_estimation(etot_f_1p75_array, wj_2)
e2, s_e2 = statistical_error_estimation(etot_f_2_array, wj_2p5)
e2p5, s_e2p5 = statistical_error_estimation(etot_f_2p5_array, wj_2p5)

print (" EF 1.25 : ", e1p25, "+-", s_e1p25, "Average Sign: ", np.mean(sign_avg_1p25))
print (" EF 1.5 : ", e1p5, "+-", s_e1p5, "Average Sign: ", np.mean(sign_avg_1p5))
print (" EF 1.75 : ", e1p75, "+-", s_e1p75, "Average Sign: ", np.mean(sign_avg_1p75))
print (" EF 2 : ", e2, "+-", s_e2, "Average Sign: ", np.mean(sign_avg_2))
print (" EF 2.5 : ", e2p5, "+-", s_e2p5, "Average Sign: ", np.mean(sign_avg_2p5))















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
# etot_f_array_5_barak = [6.93755744, 6.4148657, 6.074709, 5.9004776433, 5.6297522]
# stdv_f_array_5_barak = [0.01030199, 0.0094601, 0.05510046, 0.083725350, 0.40240035]
# sign_f_array_5_barak = [0.15715344, 0.08891504, 0.047802733, 0.026490030, 0.00804471243]
# etot_f_array_g_barak = [5.723669675071802, 5.48707494, 5.4201400039, 5.32728000]
# stdv_f_array_g_barak = [0.02653834, 0.0575122, 0.086961, 0.013275099]
# sign_f_array_g_barak = [0.26612186003, 0.1506680643, 0.09444886968, 0.0521631974]
# Unitl here barak resuls