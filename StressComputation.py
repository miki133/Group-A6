from InternalForcesAE2111 import z_values, b
from InternalForcesAE2111 import shearyz_values, momentxy_values, norm_values
from DeflectionComputation import geometry_lst, area_lst, q_design_op
from DeflectionComputation import Point3Stress, Point4Stress, Point2Stress, Point1Stress
from DeflectionComputation import nr_of_stringers, nr_of_ribs, t_ratio, spanwisesplit
import numpy as np
import matplotlib.pyplot as plt

E = 68.9*10**9
v = 0.33

#Iterables
N_str = nr_of_stringers
N_rib = nr_of_ribs
thickness = t_ratio
web_end = spanwisesplit

a = 0.5*b/(N_rib+1)



def k_c_k_s_eq(A):
    if A <1.6:
        k_c = 8.5975*A**4 - 41.785*A**3 + 76.944*A**2 - 62.87*A + 23.189
    elif 1.5<A<2.0:
        k_c = -1.73*np.log(A)+5.21
    else:
        k_c = 4

    if A <3:
        k_s = -1.2667*A**3 + 9.0857 * A**2 - 21.826*A + 23.68
    elif A>8:
        k_s = 4
    else:
        k_s = -0.25*A + 6.5333
    return k_c, k_s


def normal_stress(z):
    ind = int(z * 10)
    return norm_values[ind]/area_lst[ind][4]


def tau_ave(z):
    ind = int(z*10)
    r, h_f, h_r, h_c, k, d, e, t, theta, beta = geometry_lst[ind][0:10]
    t_f = t
    t_c = t
    t_r = t

    ind = int(z * 10)
    tau_ave = shearyz_values[ind]/(h_f*t_f + h_r*t_r + h_c*t_c)
    return tau_ave*1.5


def tau_cr(z):
    ind = int(z*10)
    r, b_1, b_2, b_3, k, d, e, t, theta, beta = geometry_lst[ind][0:10]
    t = t

    A1 = a/b_1
    A2 = a/b_2
    if b_3 != 0:
        A3 = a/b_3
    else:
        A3 = 0

    k_s1 = k_c_k_s_eq(A1)[1]
    k_s2 = k_c_k_s_eq(A2)[1]
    k_s3 = k_c_k_s_eq(A3)[1]

    tau_cr_1 = np.pi**2 *k_s1 * E/(12*(1-v**2)) *(t/b_1)**2
    tau_cr_2 = np.pi**2 *k_s2 * E/(12*(1-v**2)) *(t/b_2)**2
    if b_3 == 0:
        tau_cr_3 = 0
    else:
        tau_cr_3 = (np.pi**2 *k_s3 * E/(12*(1-v**2)) *(t/b_3)**2)
    return tau_cr_1, tau_cr_2, tau_cr_3


def sigma_cr(z):
    ind = int(z*10)
    r, b_1, b_2, b_3, k, d, e, t, theta, beta, h_1 = geometry_lst[ind][0:11]

    if q_design_op == 2:
        h_1 = h_1/(N_str/2 + 1)
    else:
        if b_3 != 0:
            h_1 = 0.5 * h_1/(N_str/4 + 1)
        else:
            h_1 = h_1/(N_str / 4 + 1)

    A = a / h_1
    k_c1 = k_c_k_s_eq(A)[0]

    sigma_cr_1 = np.pi**2 * k_c1 * E/(12*(1-v**2)) * (t/h_1)**2
    return sigma_cr_1


def torque_flow_stress(z):
    ind = int(z*10)
    t = geometry_lst[ind][7]
    q = momentxy_values[ind]/(2*area_lst[ind][3])
    return q/t


def column_stress():
    A_s = 2 * 10 ** -3
    L_s = np.sqrt(5 * A_s)
    t_s = L_s / 10  # for thin walled approximation to be true
    K = 1 / 4
    L = a
    I_xx_stiffner = 19 / 48 * L_s ** 3 * t_s
    sigma_cr_bk = (K * np.pi ** 2 * E * I_xx_stiffner) / (L ** 2 * A_s)
    return sigma_cr_bk


# Preallocate arrays
tau_ave1_values = np.empty(len(z_values))
tau_ave2_values = np.empty(len(z_values))
tau_cr1_values = np.empty(len(z_values))
tau_cr2_values = np.empty(len(z_values))
tau_cr3_values = np.empty(len(z_values))
sigma_cr_values = np.empty(len(z_values))
sigma_cr_br = np.empty(len(z_values))


# Calculate Stresses
for i in range(len(z_values)):
    tau_ave1_values[i] = tau_ave(z_values[i]) - torque_flow_stress(z_values[i])
    tau_ave2_values[i] = tau_ave(z_values[i]) + torque_flow_stress(z_values[i])

    tau_cr1, tau_cr2, tau_cr3 = tau_cr(z_values[i])
    tau_cr1_values[i] = tau_cr1
    tau_cr2_values[i] = tau_cr2
    tau_cr3_values[i] = tau_cr3

    sigma_cr_values[i] = sigma_cr(z_values[i])
    sigma_cr_br[i] = column_stress()
    Point1Stress[i] = Point1Stress[i] - normal_stress(z_values[i])
    Point2Stress[i] = Point2Stress[i] - normal_stress(z_values[i])
    Point3Stress[i] = Point3Stress[i] - normal_stress(z_values[i])
    Point4Stress[i] = Point4Stress[i] - normal_stress(z_values[i])


fig, axs = plt.subplots(2, 3, figsize=(13, 5), layout='constrained')

axs.flat[0].plot(z_values, tau_ave1_values, color='red')
axs.flat[0].plot(z_values, tau_cr1_values, color='blue')

axs.flat[1].plot(z_values, tau_ave2_values, color='red')
axs.flat[1].plot(z_values, tau_cr2_values, color='blue')

axs.flat[2].plot(z_values, tau_cr3_values, color='blue')
axs.flat[2].plot(z_values, tau_ave2_values, color='red')

axs.flat[3].plot(z_values, Point1Stress, color='red')
axs.flat[3].plot(z_values, Point2Stress, color='red')
axs.flat[3].plot(z_values, Point3Stress, color='red')
axs.flat[3].plot(z_values, Point4Stress, color='red')
axs.flat[3].plot(z_values, sigma_cr_values, color='blue')

axs.flat[4].plot(z_values, sigma_cr_br, color='blue')
axs.flat[4].plot(z_values, Point3Stress, color='red')
axs.flat[4].plot(z_values, Point4Stress, color='red')

plt.show()

####################################################

# exceeds_limit2 = False
#
# for i in range(len(z_values)):
#     tau_ave2 = tau_ave2_values[i]
#     tau_cr2 = tau_cr2_values[i]
#     z_location = z_values[i]
#
#     if tau_ave2 > tau_cr2:
#         print(f"At z = {z_location}, tau_ave ({tau_ave}) exceeds tau_cr ({tau_cr})")
#         exceeds_limit2 = True
#
# if exceeds_limit2:
#     print("At least one pair of tau_ave2 exceeds tau_cr2")
# else:
#     print("No pair of tau_ave2 exceeds tau_cr2")
#
# exceeds_limit1 = False
#
# for i in range(len(z_values)):
#     tau_ave1 = tau_ave1_values[i]
#     tau_cr1 = tau_cr1_values[i]
#     z_location = z_values[i]
#
#     if tau_ave1 > tau_cr1:
#         print(f"At z = {z_location}, tau_ave1 ({tau_ave1}) exceeds tau_cr1 ({tau_cr1})")
#         exceeds_limit1 = True
#         break
#
# if exceeds_limit1:
#     print("Limit exceeded for tau_ave1 and tau_cr1")
# else:
#     print("No pair of tau_ave1 exceeds tau_cr1")
#
# exceeds_limit3 = False
# half_length = len(z_values) // 2  # Assuming len(z_values) is even
#
# for i in range(half_length):
#     tau_ave3 = tau_ave2_values[i]
#     tau_cr3 = tau_cr3_values[i]
#     z_location = z_values[i]
#
#     if tau_ave3 > tau_cr3:
#         print(f"At z = {z_location}, tau_ave3 ({tau_ave3}) exceeds tau_cr3 ({tau_cr3})")
#         exceeds_limit3 = True
#         break
#
# if exceeds_limit3:
#     print("Limit exceeded for tau_ave3 and tau_cr3 in the first half of the span")
# else:
#     print("No pair of tau_ave3 exceeds tau_cr3 in the first half of the span")


# Perform comparisons using vectorized operations

def failure_test():
    tau_ave1_fail = tau_ave1_values > tau_cr1_values
    tau_ave2_fail = tau_ave2_values > tau_cr2_values
    tau_ave3_fail = np.logical_and(tau_ave2_values > tau_cr3_values, tau_cr3_values != 0)
    point4_sigma_fail = Point4Stress > sigma_cr_values
    point3_sigma_fail = Point3Stress > sigma_cr_values

    Failure = tau_ave1_fail.any() or tau_ave2_fail.any() or tau_ave3_fail.any() or \
              point3_sigma_fail.any() or point4_sigma_fail.any()
    return Failure

# print(Failure)
