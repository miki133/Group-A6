import numpy as np
import scipy as sp

CD0_0deg = 0.004713
CD0_10deg = 0.046930
CM0_0deg = -0.486076
CM0_10deg = -1.38817
max_load_density = 1.225
max_load_velocity = 193
Cl_des = 0.6645720488089126
Cd_des = (0.647562 * CD0_10deg + CD0_0deg) / (1 + 0.647562)
Cm_des = (0.647562 * CM0_10deg + CM0_0deg) / (1 + 0.647562)
CL0_0deg = 0.355876
CL0_10deg = 1.141277

y_list_10deg = np.array([1.0079,3.0168,5.005,6.9591,8.8656,10.7116,12.4844,14.172,15.7627,17.2457,18.611,19.8491,20.9516,21.911,22.7207,23.3752,23.8701,24.2019,24.3683])
lift_coefficient_list_0deg = np.array([0.282096,0.304281, 0.321909, 0.336951, 0.350082, 0.36165,  0.371837, 0.380699,
 0.388162, 0.394007, 0.397804, 0.398821, 0.395875, 0.387164, 0.370173, 0.341672,
 0.297375, 0.231283, 0.145097])
lift_coefficient_list_10deg = np.array([0.883134,0.961389,1.019523,1.069302,1.112873,1.151257,1.184958,1.214071,1.238245,1.256577,1.267348,1.267608,1.252565,1.214896,1.1447,1.031093,0.865262, 0.644307, 0.393162])
induced_drag_coefficient_list_10deg = np.array( [0.025671,0.038217,0.042483,0.044196,0.044749,0.044783,0.044705,0.044849,0.045546,0.047199,0.050294,0.055398,0.062996,0.073024,0.084219,0.093937,0.098426 ,0.091332, 0.078368])
induced_drag_coefficient_list_0deg = np.array( [0.002736, 0.003955, 0.004272, 0.004359, 0.004356, 0.004316, 0.004274, 0.004259,
 0.004299, 0.004433, 0.004707, 0.005187, 0.005951, 0.007086, 0.008682, 0.010871,
 0.013667, 0.015753, 0.01499 ])
pitching_moment_quarterchord_list_10deg = np.array( [-0.121113,-0.105958,-0.100657, -0.098559, -0.09773,  -0.097395, -0.097202,-0.09697,  -0.096552, -0.095776, -0.094341, -0.091759, -0.087201, -0.079401,-0.066986, -0.049743, -0.03042,  -0.017746,  0.022791])
pitching_moment_quarterchord_list_0deg = np.array( [-0.107224, -0.105085, -0.104375, -0.104116, -0.104059, -0.104083, -0.104125,
 -0.10415,  -0.10412,  -0.103995, -0.103694, -0.103085, -0.101931, -0.099836,
 -0.096236, -0.090524, -0.082149, -0.071452, -0.046055])
chord_list_10deg = np.array( [8.2522, 7.7842, 7.3209, 6.8657, 6.4215, 5.9914, 5.5784, 5.1852, 4.8146, 4.4691,4.151,  3.8626, 3.6057, 3.3822, 3.1936, 3.0411, 2.9258, 2.8485, 2.8097])
chord_list_0deg = np.array( [8.2522, 7.7842, 7.3209, 6.8657, 6.4215, 5.9914, 5.5784, 5.1852, 4.8146, 4.4691,4.151,  3.8626, 3.6057, 3.3822, 3.1936, 3.0411, 2.9258, 2.8485, 2.8097])



chord_function_10deg = sp.interpolate.interp1d(y_list_10deg, chord_list_10deg, kind='cubic', fill_value="extrapolate")
CL_function_10deg = sp.interpolate.interp1d(y_list_10deg, lift_coefficient_list_10deg, kind='cubic', fill_value="extrapolate")
CD_function_10deg = sp.interpolate.interp1d(y_list_10deg, induced_drag_coefficient_list_10deg, kind='cubic', fill_value="extrapolate")
CM_function_10deg = sp.interpolate.interp1d(y_list_10deg, pitching_moment_quarterchord_list_10deg, kind='cubic', fill_value="extrapolate")
CL_function_0deg = sp.interpolate.interp1d(y_list_10deg, lift_coefficient_list_0deg, kind='cubic', fill_value="extrapolate")
CD_function_0deg = sp.interpolate.interp1d(y_list_10deg, induced_drag_coefficient_list_10deg, kind='cubic', fill_value="extrapolate")
CM_function_0deg = sp.interpolate.interp1d(y_list_10deg, pitching_moment_quarterchord_list_10deg, kind='cubic', fill_value="extrapolate")


def lift_force_function(z):
    return 0.5 * max_load_velocity ** 2 * max_load_density * chord_function_10deg(z) * (CL_function_0deg(z) + (Cl_des - CL0_0deg) / (CL0_10deg - Cl_des) * (CL_function_10deg(z) - CL_function_0deg(z)))


def drag_force_function(z):
    return 0.5 * max_load_velocity ** 2 * max_load_density * chord_function_10deg(z) * (CD_function_0deg(z) + (
                Cd_des - CD0_0deg) / (CD0_10deg - Cd_des) * (CD_function_10deg(z) - CD_function_0deg(z)))


def moment_function(z):
    return 0.5 * max_load_velocity ** 2 * max_load_density * chord_function_10deg(z) * (CM_function_0deg(z) + (
                Cm_des - CM0_0deg) / (CM0_10deg - Cm_des) * (CM_function_10deg(z) - CM_function_0deg(z)))