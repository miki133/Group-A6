import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from Wingloading import Lift_func, Drag_func, Moment_func, chord_function_0deg, n, alpha

# n = defined in wingloading
b = 48.82
P_eng = n*135787
T = 431000
M_eng = 0
# AoA = np.radians(12.044974098211046)
AoA = alpha
c4angle = np.radians(20.05)
eng_loc = 8.91
eng_loc_vert = 5


def weightfunc(z):
    weight = n*(6688.42 - 183.56*z)
    return weight


def vertforce(z):
    F_v = (Lift_func(z)*np.cos(AoA) - Drag_func(z)*np.sin(AoA) - weightfunc(z))*1.5
    return F_v


def horforce(z):
    F_h = (Lift_func(z)*np.sin(AoA) + Drag_func(z)*np.cos(AoA))*1.5
    return F_h

def normalfunc(z):
    N = (-T*np.sin(c4angle)*(1 - np.heaviside(z - eng_loc, 1)))*1.5
    return N


def shearfuncyz(z):
    ind = np.where(z_values == z)[0][0]
    new_z = [j for j in z_values if j >= z]
    new_vert = vert_values[ind:]
    shear = -P_eng * (1 - np.heaviside(z - eng_loc, 1)) + sp.integrate.simps(new_vert, new_z)
    return shear


def shearfuncxz(z):
    ind = np.where(z_values == z)[0][0]
    new_z = [j for j in z_values if j >= z]
    new_hor = hor_values[ind:]
    shear = -T * np.cos(c4angle) * (1 - np.heaviside(z - eng_loc, 1)) + sp.integrate.simps(new_hor, new_z)
    return shear


def momentfuncyz(z):
    ind = np.where(z_values == z)[0][0]
    new_z = [j for j in z_values if j >= z]
    new_shear = shearyz_values[ind:]
    moment = - sp.integrate.simps(new_shear, new_z) - M_eng*(1 - np.heaviside(z - eng_loc, 1))
    return moment


def momentfuncxz(z):
    ind = np.where(z_values == z)[0][0]
    new_z = [j for j in z_values if j >= z]
    new_shear = shearxz_values[ind:]
    moment = - sp.integrate.simps(new_shear, new_z) - M_eng*(1 - np.heaviside(z - eng_loc, 1))
    return moment   


def torquefunc(z): #Fix moment of the thrust
    ind = np.where(z_values == z)[0][0]
    new_z = [j for j in z_values if j >= z]
    new_vert_mom_arm = vert_mom_arm[ind::]
    # new_vert = vert_values[ind::]*np.tan(c4angle)*new_z
    new_vert = vert_values[ind::]*new_vert_mom_arm*np.cos(AoA)
    # new_hor = hor_values[ind::]*np.tan(np.radians(5))/np.cos(c4angle)*new_z
    # print(new_vert
    # print(new_vert_mom_arm)
    new_hor = hor_values[ind::]*0
    new_mom = mom_values[ind::]
    moment = sp.integrate.simps(new_vert, new_z) \
    - T * np.cos(c4angle) * (0.412 + 1.5) * (1 - np.heaviside(z - eng_loc, 1)) #/ (
        # np.cos(c4angle)) \
    - sp.integrate.simps(new_hor, new_z) \
    + P_eng * (8.25 * 0.25 + 1.416) * (1 - np.heaviside(z - eng_loc, 1)) \
    + sp.integrate.simps(new_mom, new_z)

    # - T* np.cos(c4angle)*(eng_loc*np.tan(np.radians(5))-1.912)*(1 - np.heaviside(z - eng_loc, 1))/(np.cos(c4angle)) \
    # + P_eng*(-8.25*0.25 - 1.65 - 0.5*5.812 + eng_loc*np.tan(c4angle))*(1 - np.heaviside(z - eng_loc, 1))  \
    return moment


z_values = np.arange(0, b/2 + 0.1, 0.1)

vert_values = np.empty((0, 1))
hor_values = np.empty((0, 1))
mom_values = np.empty((0, 1))
norm_values = np.empty((0, 1))

drag_values = np.empty((0, 1))
lift_values = np.empty((0, 1))

vert_mom_arm = np.empty((0, 1))
hor_mom_arm = np.empty((0, 1))

shearyz_values = np.empty((0, 1))
shearxz_values = np.empty((0, 1))
momentyz_values = np.empty((0, 1))
momentxz_values = np.empty((0, 1))
momentxy_values = np.empty((0, 1))

#Calculation Distributed Forces at each point
for i in range(len(z_values)):
    vert_values = np.append(vert_values, vertforce(z_values[i]))
    hor_values = np.append(hor_values, horforce(z_values[i]))
    mom_values = np.append(mom_values, Moment_func(z_values[i]))
    lift_values = np.append(lift_values, Lift_func(z_values[i]))
    drag_values = np.append(drag_values, Drag_func(z_values[i]))
    norm_values = np.append(norm_values, normalfunc(z_values[i]))
    vert_mom_arm = np.append(vert_mom_arm, -0.25*chord_function_0deg(z_values[i]) + chord_function_0deg(z_values[i])*0.433)

#Calculation Shear Forces
for i in range(len(z_values)):
    shearyz_values = np.append(shearyz_values, shearfuncyz(z_values[i]))
    shearxz_values = np.append(shearxz_values, shearfuncxz(z_values[i]))

#Calculation Moments
for i in range(len(z_values)):
    momentyz_values = np.append(momentyz_values, momentfuncyz(z_values[i]))
    momentxz_values = np.append(momentxz_values, momentfuncxz(z_values[i]))
    momentxy_values = np.append(momentxy_values, torquefunc(z_values[i]))


fig, axs = plt.subplots(2, 3, figsize=(13, 5), layout='constrained')
# axs[0][2].set_visible(False)



axs.flat[2].plot(z_values, norm_values, color='red')
axs.flat[2].set_title("Internal Normal Force")
axs.flat[2].set(ylabel='Internal Normal Load [N]')
axs.flat[2].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# axs.flat[2].axhline(0, color='black')


axs.flat[0].plot(z_values, shearyz_values, color='blue')
axs.flat[0].set_title('Internal Shear yz-plane')
axs.flat[0].set(ylabel='Internal Shear [N]')
axs.flat[0].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# axs.flat[0].axhline(0, color='black')


axs.flat[1].plot(z_values, shearxz_values, color='blue')
axs.flat[1].set_title('Internal Shear xz-plane')
# axs.flat[1].axhline(0, color='black')
axs.flat[1].set(ylabel='Internal Shear [N]')
# axs.flat[1].ticklabel_format(style='sci', axis='y', scilimits=(6,6))


axs.flat[3].plot(z_values, momentyz_values, color='orange')
axs.flat[3].set_title('Internal Moment yz-plane')
axs.flat[3].set(ylabel='Internal Moment [Nm]')
axs.flat[3].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# axs.flat[3].axhline(0, color='black')


axs.flat[4].plot(z_values, momentxz_values, color='orange')
axs.flat[4].set_title('Internal Moment xz-plane')
axs.flat[4].set(ylabel='Internal Moment [Nm]')
axs.flat[4].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# axs.flat[4].axhline(0, color='black')


axs.flat[5].plot(z_values, momentxy_values, color='green')
axs.flat[5].set_title('Internal Torque')
axs.flat[5].set(ylabel='Internal Moment [Nm]')
# axs.flat[5].ticklabel_format(style='sci', axis='y', scilimits=(6,6))
# axs.flat[5].axhline(0, color='black')


for ax in axs.flat:
    ax.set(xlabel='Spanwise Location [m]')
    # ax.label_outer()

plt.show()
