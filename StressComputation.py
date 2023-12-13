from InternalForcesAE2111 import z_values, shearfuncyz, torquefunc, normalfunc, b
from DeflectionComputation import a_ratio, p_ratio, d_ratio, t_ratio, k_ratio, q_design_op
# from MaxiMagic import Point1Stress, Point2Stress, Point3Stress, Point4Stress
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

    h_1 = np.sqrt(d ** 2 + r ** 2)
    h_2 = np.sqrt(e ** 2 + r ** 2)
    return r, b_1, b_2, b_3, k, d, e, t, theta, beta, h_1, h_2


def enclosed_area(z):
    r, a, d, e = geometry(z)[0], geometry(z)[2], geometry(z)[5], geometry(z)[6]
    A_1 = a*r
    A_2 = 1/2*r*d
    A_3 = 1/2*r*e
    return A_3 + A_2 + A_1

def cross_sec_area(z):
    r, b_1, b_2, b_3, k, d, e, t, theta, beta, h_1, h_2 = geometry(z)
    A_str = 2*10**-3
    if q_design_op == 1:
        N = 8
    elif q_design_op == 2:
        N = 20
    else:
        N = 4
    A = (b_1 + b_2 + b_3 + h_1 + h_2)*t + N*A_str
    return A

def normal_stress(z):
    return normalfunc(z)/cross_sec_area(z)


def tau_ave(z):
    r, h_f, h_r, h_c, k, d, e, t, theta, beta = geometry(z)[0:10]
    t_f = t
    t_c = t
    t_r = t

    tau_ave = shearfuncyz(z)/(h_f*t_f + h_r*t_r + h_c*t_c)
    return tau_ave*1.5


def tau_cr(z):
    #beam 1 is LE, beam 2 is TE, beam 3 is middle
    r, b_1, b_2, b_3, k, d, e, t, theta, beta = geometry(z)[0:10]
    t = 2*t

    k_s = 9

    tau_cr_1 = np.pi**2 *k_s * E/(12*(1-v**2)) *(t/b_1)**2
    tau_cr_2 = np.pi**2 *k_s * E/(12*(1-v**2)) *(t/b_2)**2
    if b_3 == 0:
        tau_cr_3 = 0
    else:
        tau_cr_3 = (np.pi**2 *k_s * E/(12*(1-v**2)) *(t/b_3)**2)#*on_off_factor
    return tau_cr_1, tau_cr_2, tau_cr_3

def sigma_cr(z):
    #beam 1 is LE, beam 2 is TE, beam 3 is middle
    r, b_1, b_2, b_3, k, d, e, t, theta, beta = geometry(z)[0:10]

    h_1 = np.sqrt(d ** 2 + r ** 2)
    h_2 = np.sqrt(e ** 2 + r ** 2)
    k_c = 4

    sigma_cr_1 = np.pi**2 *k_c * E/(12*(1-v**2)) *(t/h_1)**2
    sigma_cr_2 = np.pi**2 *k_c * E/(12*(1-v**2)) *(t/h_2)**2
    return sigma_cr_1, sigma_cr_2


def torque_flow_stress(z):
    t = geometry(z)[7]
    q = torquefunc(z)/(2*enclosed_area(z))
    return q/t

tau_ave1_values = np.empty((0, 1))
tau_ave2_values = np.empty((0, 1))
tau_cr1_values = np.empty((0, 1))
tau_cr2_values = np.empty((0, 1))
tau_cr3_values = np.empty((0, 1))
sigma_cr1_values = np.empty((0, 1))
sigma_cr2_values = np.empty((0, 1))


sigma_values = np.empty((0, 1))


for i in range(len(z_values)):
    tau_ave1_values = np.append(tau_ave1_values, tau_ave(z_values[i]) - torque_flow_stress(z_values[i]))
    tau_ave2_values = np.append(tau_ave2_values, tau_ave(z_values[i]) + torque_flow_stress(z_values[i]))

    tau_cr1, tau_cr2, tau_cr3 = tau_cr(z_values[i])
    tau_cr1_values = np.append(tau_cr1_values, tau_cr1)
    tau_cr2_values = np.append(tau_cr2_values, tau_cr2)
    tau_cr3_values = np.append(tau_cr3_values, tau_cr3)

    sigma_cr1, sigma_cr2 = sigma_cr(z_values[i])
    sigma_cr1_values = np.append(sigma_cr1_values, sigma_cr1)
    sigma_cr2_values = np.append(sigma_cr2_values, sigma_cr2)

    sigma_values = np.append(sigma_values, normal_stress(z_values[i]))






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

# axs.flat[3].plot(z_values, Point1Stress, color='red')
# axs.flat[3].plot(z_values, Point2Stress, color='red')
# axs.flat[3].plot(z_values, -Point3Stress, color='red')
# axs.flat[3].plot(z_values, -Point4Stress, color='red')
axs.flat[3].plot(z_values, sigma_cr1_values, color='blue')
axs.flat[3].plot(z_values, sigma_cr2_values, color='blue')


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


#stiffner calculations for colum buckling
A_s = 3*10**-3
L_s = np.sqrt(5*A_s)
t_s = L_s/10 #for thin walled approximation to be true
K = 1/4
L = 24.41
I_xx_stiffner = 19/48*L_s**3*t_s
sigma_cr_bk = (K*np.pi**2*E*I_xx_stiffner)/(L**2*A_s)

exceeds_limit2 = False

for i in range(len(z_values)):
    tau_ave2 = tau_ave2_values[i]
    tau_cr2 = tau_cr2_values[i]
    z_location = z_values[i]

    if tau_ave2 > tau_cr2:
        print(f"At z = {z_location}, tau_ave ({tau_ave}) exceeds tau_cr ({tau_cr})")
        exceeds_limit2 = True

if exceeds_limit2:
    print("At least one pair of tau_ave2 exceeds tau_cr2")
else:
    print("No pair of tau_ave2 exceeds tau_cr2")
    
    
exceeds_limit1 = False

for i in range(len(z_values)):
    tau_ave1 = tau_ave1_values[i]
    tau_cr1 = tau_cr1_values[i]
    z_location = z_values[i]

    if tau_ave1 > tau_cr1:
        print(f"At z = {z_location}, tau_ave1 ({tau_ave1}) exceeds tau_cr1 ({tau_cr1})")
        exceeds_limit1 = True
        break

if exceeds_limit1:
    print("Limit exceeded for tau_ave1 and tau_cr1")
else:
    print("No pair of tau_ave1 exceeds tau_cr1")

exceeds_limit3 = False
half_length = len(z_values) // 2  # Assuming len(z_values) is even

for i in range(half_length):
    tau_ave3 = tau_ave2_values[i]
    tau_cr3 = tau_cr3_values[i]
    z_location = z_values[i]

    if tau_ave3 > tau_cr3:
        print(f"At z = {z_location}, tau_ave3 ({tau_ave3}) exceeds tau_cr3 ({tau_cr3})")
        exceeds_limit3 = True
        break

if exceeds_limit3:
    print("Limit exceeded for tau_ave3 and tau_cr3 in the first half of the span")
else:
    print("No pair of tau_ave3 exceeds tau_cr3 in the first half of the span")
