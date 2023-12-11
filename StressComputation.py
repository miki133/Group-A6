from InternalForcesAE2111 import z_values, shearfuncyz, torquefunc, b
from DeflectionComputation import a_ratio, p_ratio, d_ratio, t_ratio, k_ratio, q_design_op, Jcalculation
import numpy as np
import matplotlib.pyplot as plt

E = 68.9*10**9
v = 0.33

def r(z):
    r = 0.5 * (8.487 - 0.233 * z)
    return r

def geometry(z):
    c = (8.487 - 0.233 * z)
    r = 0.5 * c
    #Front/Rearweb
    b_1 = c * p_ratio
    b_2 = c * a_ratio
    d = c * d_ratio
    t = t_ratio*(0.003125 * r - 0.00235)
    k = k_ratio * r
    e = b_1 - d - b_2
    theta = np.arctan(d / r)
    beta = np.arctan(e / r)

    # b_3 = b_2 + (r - k) * np.tan(beta) + (r - k) * np.tan(theta)
    if q_design_op == 2 or z > b / 4:
        b_3 = 0
    else:
        b_3 = b_2 + (r - k) * np.tan(beta) + (r - k) * np.tan(theta)
    return r, b_1, b_2, b_3, k, d, e, t, theta, beta


def enclosed_area(z):
    r, a, d, e = geometry(z)[0], geometry(z)[2], geometry(z)[5], geometry(z)[6]
    A_1 = a*r
    A_2 = 1/2*r*d
    A_3 = 1/2*r*e
    return A_3 + A_2 + A_1


def tau_ave(z):
    r, h_f, h_r, h_c, k, d, e, t, theta, beta = geometry(z)
    t_f = t
    t_c = t
    t_r = t

    tau_ave = shearfuncyz(z)/(h_f*t_f + h_r*t_r + h_c*t_c)
    return tau_ave*1.5


def tau_cr(z):
    #beam 1 is LE, beam 2 is TE, beam 3 is middle
    # r, b_1, b_2, b_3, k, d, e, t, theta, beta = geometry(z)
    # print("taucr", r, b_1, b_2, b_3, k, d, e, t, theta, beta)

    t = t_ratio * (0.003125 * r(z) - 0.00235)

    b_1 = p_ratio*r(z)          #height LE beam at point z
    b_2 = a_ratio*r(z)          #height TE beam at point z

    d = r(z) * d_ratio          #"gap" above TE beam at point z
    e = b_1 - d - b_2           #"gap" below TE beam at point z
    k = k_ratio*r(z)            #Distance in x direction of centre beam

    theta = np.arctan(d / r(z))
    beta = np.arctan(e / r(z))

    b_3 = b_2 + (r(z) - k) * np.tan(beta) + (r(z) - k) * np.tan(theta)

    k_s = 9


    if q_design_op == 2 or z >b/4:
        on_off_factor = 0
    else:
        on_off_factor = 1

    tau_cr_1 = np.pi**2 *k_s * E/(12*(1-v**2)) *(t/b_1)**2
    tau_cr_2 = np.pi**2 *k_s * E/(12*(1-v**2)) *(t/b_2)**2
    tau_cr_3 = (np.pi**2 *k_s * E/(12*(1-v**2)) *(t/b_3)**2)*on_off_factor
    return tau_cr_1, tau_cr_2, tau_cr_3


def torque_flow_stress(z):
    t = geometry(z)[7]
    q = torquefunc(z)/(2*enclosed_area(z))
    return q/t

tau_ave1_values = np.empty((0, 1))
tau_ave2_values = np.empty((0, 1))
tau_cr1_values = np.empty((0, 1))
tau_cr2_values = np.empty((0, 1))
tau_cr3_values = np.empty((0, 1))


for i in range(len(z_values)):
    tau_ave1_values = np.append(tau_ave1_values, tau_ave(z_values[i]) - torque_flow_stress(z_values[i]))
    tau_ave2_values = np.append(tau_ave2_values, tau_ave(z_values[i]) + torque_flow_stress(z_values[i]))
    tau_cr1, tau_cr2, tau_cr3 = tau_cr(z_values[i])
    tau_cr1_values = np.append(tau_cr1_values, tau_cr1)
    tau_cr2_values = np.append(tau_cr2_values, tau_cr2)
    tau_cr3_values = np.append(tau_cr3_values, tau_cr3)






fig, axs = plt.subplots(2, 3, figsize=(13, 5), layout='constrained')
# # axs[0][2].set_visible(False)
#
#
#
axs.flat[2].plot(z_values, tau_cr3_values, color='blue')
axs.flat[2].plot(z_values, tau_ave2_values, color='red')

axs.flat[1].plot(z_values, tau_ave2_values, color='red')
axs.flat[1].plot(z_values, tau_cr2_values, color='blue')

axs.flat[0].plot(z_values, tau_ave1_values, color='red')
axs.flat[0].plot(z_values, tau_cr1_values, color='blue')


# axs.flat[2].set_title("Internal Normal Force")
# axs.flat[2].set(ylabel='Internal Normal Load [N]')
# axs.flat[2].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# axs.flat[2].axhline(0, color='black')
#
#
# axs.flat[0].plot(z_values, shearyz_values, color='blue')
# axs.flat[0].set_title('Internal Shear yz-plane')
# # axs.flat[0].set(ylabel='Internal Shear [N]')
# # axs.flat[0].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# # # axs.flat[0].axhline(0, color='black')
# #
# # axs.flat[1].plot(z_values, shearxz_values, color='blue')
# # axs.flat[1].set_title('Internal Shear xz-plane')
# # # axs.flat[1].axhline(0, color='black')
# # axs.flat[1].set(ylabel='Internal Shear [N]')
# # # axs.flat[1].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
#
#
# axs.flat[2].plot(z_values, tau_cr1_values, color='orange', label='Leading Edge')
# axs.flat[3].set_title('Internal Moment yz-plane')
# axs.flat[3].set(ylabel='Internal Moment [Nm]')
# axs.flat[3].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# # axs.flat[3].axhline(0, color='black')
#
#
# axs.flat[2].plot(z_values, tau_cr2_values, color='blue', label='Trailing Edge')
# axs.flat[4].set_title('Internal Moment xz-plane')
# axs.flat[4].set(ylabel='Internal Moment [Nm]')
# axs.flat[4].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# # axs.flat[4].axhline(0, color='black')
#
#
# axs.flat[2].plot(z_values, tau_cr3_values, color='green')
# axs.flat[5].set_title('Internal Torque')
# axs.flat[5].set(ylabel='Internal Moment [Nm]')
# # axs.flat[5].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# # axs.flat[5].axhline(0, color='black')

plt.show()
