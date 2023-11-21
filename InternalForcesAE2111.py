#Code written by Robin van Veldhuizen, 5727227
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from Mikimagic import lift_force_function, drag_force_function, moment_function

b = 48.82
P_eng = 135787
T = 431000
M_eng = 0
AoA = np.radians(6.456449)
c4angle = np.radians(20.05)
eng_loc = 8.91
eng_loc_vert = 5


def weightfunc(z):
    weight = 6859 - 188.2*z
    return weight


def vertforce(z):
    F_v = lift_force_function(z)*np.cos(AoA) - drag_force_function(z)*np.sin(AoA) - weightfunc(z)
    return F_v


def horforce(z):
    F_h = lift_force_function(z)*np.sin(AoA) + drag_force_function(z)*np.cos(AoA)
    return F_h


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
    shear = -T * (1 - np.heaviside(z - eng_loc, 1)) + sp.integrate.simps(new_hor, new_z)
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


def torquefunc(z):
    ind = np.where(z_values == z)[0][0]
    new_z = [j for j in z_values if j >= z]
    new_vert = vert_values[ind::]*np.tan(c4angle)*new_z
    new_hor = hor_values[ind::]*np.tan(np.radians(5))/np.cos(c4angle)*z
    new_mom = mom_values[ind::]
    moment = sp.integrate.simps(new_vert, new_z) \
    + T*(eng_loc*np.tan(np.radians(5))-1.912)*(1 - np.heaviside(z - eng_loc, 1))/(np.cos(c4angle)) \
    - sp.integrate.simps(new_hor, new_z) \
    - P_eng*(-8.25*0.25 - 1.65 - 0.5*5.812 + eng_loc*np.tan(c4angle))*(1 - np.heaviside(z - eng_loc, 1))  \
    + sp.integrate.simps(new_mom, new_z)
    return moment


z_values = np.arange(0, b/2 + 0.1, 0.1)

vert_values = np.empty((0, 1))
hor_values = np.empty((0, 1))
mom_values = np.empty((0, 1))
shearyz_values = np.empty((0, 1))
shearxz_values = np.empty((0, 1))
momentyz_values = np.empty((0, 1))
momentxz_values = np.empty((0, 1))
momentxy_values = np.empty((0, 1))

#Calculation Distributed Forces at each point
for i in range(len(z_values)):
    vert_values = np.append(vert_values, vertforce(z_values[i]))
    hor_values = np.append(hor_values, horforce(z_values[i]))
    mom_values = np.append(mom_values, moment_function(z_values[i]))

#Calculation Shear Forces
for i in range(len(z_values)):
    shearyz_values = np.append(shearyz_values, shearfuncyz(z_values[i]))
    shearxz_values = np.append(shearxz_values, shearfuncxz(z_values[i]))

#Calculation Moments
for i in range(len(z_values)):
    momentyz_values = np.append(momentyz_values, momentfuncyz(z_values[i]))
    momentxz_values = np.append(momentxz_values, momentfuncxz(z_values[i]))
    momentxy_values = np.append(momentxy_values, torquefunc(z_values[i]))

print(sp.integrate.quad(horforce, 0, b/2))
print(sp.integrate.quad(vertforce, 0, b/2))
print(sp.integrate.quad(lift_force_function, 0, b/2))


fig, axs = plt.subplots(2, 3, figsize=(15, 5), layout='constrained')

axs.flat[0].plot(z_values, shearyz_values)
axs.flat[0].set_title('Internal Shear yz-plane')

axs.flat[1].plot(z_values, shearxz_values)
axs.flat[1].set_title('Internal Shear xz-plane')

axs.flat[3].plot(z_values, momentyz_values)
axs.flat[3].set_title('Internal Moment yz-plane')

axs.flat[4].plot(z_values, momentxz_values)
axs.flat[4].set_title('Internal Moment xz-plane')

axs.flat[2].plot(z_values, momentxy_values)
axs.flat[2].set_title('Internal Moment xy-plane')

plt.show()