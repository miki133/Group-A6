from InternalForcesAE2111 import b
from InternalForcesAE2111 import shearyz_values, momentxy_values, norm_values
from DeflectionComputation import q_design_op
import numpy as np

E = 68.9*10**9
v = 0.33


def k_c_k_s_eq(A):
    if q_design_op !=2:
        if A <1.6:
            k_c = 8.5975*A**4 - 41.785*A**3 + 76.944*A**2 - 62.87*A + 23.189
        elif 1.5<A<2.0:
            k_c = -1.73*np.log(A)+5.21
        else:
            k_c = 4
    else:
        if A != 0:
            k_c = 2.0896*A**(-1.155)
        else:
            k_c = 0

    if A <3:
        k_s = -1.2667*A**3 + 9.0857 * A**2 - 21.826*A + 23.68
    elif A>8:
        k_s = 4
    else:
        k_s = -0.25*A + 6.5333
    return k_c, k_s


def normal_stress(z, area_lst):
    ind = int(z * 10)
    return norm_values[ind]/area_lst[ind][4]


def tau_ave(z, geometry_lst):
    ind = int(z*10)
    r, h_f, h_r, h_c, k, d, e, t, theta, beta = geometry_lst[ind][0:10]
    t_f = t
    t_c = t
    t_r = t

    ind = int(z * 10)
    tau_ave = shearyz_values[ind]/(h_f*t_f + h_r*t_r + h_c*t_c)
    return tau_ave*1.5


def tau_cr(z, geometry_lst, N_rib):
    ind = int(z*10)
    a = 0.5 * b / (N_rib + 1)

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


def sigma_cr(z, geometry_lst, N_rib, N_str):
    ind = int(z*10)
    a = 0.5 * b / (N_rib + 1)
    r, b_1, b_2, b_3, k, d, e, t, theta, beta, h_1 = geometry_lst[ind][0:11]

    # if q_design_op == 2:
    #     h_1 = h_1/(N_str/2 + 1)
    # else:
    if b_3 != 0:
        h_1 = 0.5 * h_1/(N_str/4 + 1)
    else:
        h_1 = h_1/(N_str / 4 + 1)

    A = a / h_1
    k_c1 = k_c_k_s_eq(A)[0]

    sigma_cr_1 = np.pi**2 * k_c1 * E/(12*(1-v**2)) * (t/h_1)**2
    return sigma_cr_1


def torque_flow_stress(z, geometry_lst, area_lst):
    ind = int(z*10)
    t = geometry_lst[ind][7]
    q = momentxy_values[ind]/(2*area_lst[ind][3])
    return q/t


def column_stress(N_rib):
    a = 0.5 * b / (N_rib + 1)
    A_s = 2 * 10 ** -3
    L_s = np.sqrt(5 * A_s)
    t_s = L_s / 10  # for thin walled approximation to be true
    K = 1 / 4
    L = a
    I_xx_stiffner = 19 / 48 * L_s ** 3 * t_s
    sigma_cr_bk = (K * np.pi ** 2 * E * I_xx_stiffner) / (L ** 2 * A_s)
    return sigma_cr_bk


def failure_test(tau_ave1_values, tau_ave2_values, tau_cr1_values, tau_cr2_values, tau_cr3_values,
                 Point4Stress, Point3Stress, sigma_cr_values, sigma_cr_br):
    # Define conditions directly in the `any()` function to avoid intermediate variables
    return any([
        np.any(tau_ave1_values > tau_cr1_values),
        np.any(tau_ave2_values > tau_cr2_values),
        np.any(np.logical_and(tau_ave2_values > tau_cr3_values, tau_cr3_values != 0)),
        np.any(Point3Stress > sigma_cr_values),
        np.any(Point4Stress > sigma_cr_values),
        np.any(Point3Stress > sigma_cr_br),
        np.any(Point4Stress > sigma_cr_br)
    ])
