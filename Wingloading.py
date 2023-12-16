import scipy as sp
import numpy as np
import os
import matplotlib.pyplot as plt

y_list_0deg = np.empty((0, 1))
chord_list_0deg = np.empty((0, 1))
induced_angle_list_0deg = np.empty((0, 1))
lift_coefficient_list_0deg = np.empty((0, 1))
induced_drag_coefficient_list_0deg = np.empty((0, 1))
pitching_moment_quarterchord_list_0deg = np.empty((0, 1))
y_list_10deg = np.empty((0, 1))
chord_list_10deg = np.empty((0, 1))
induced_angle_list_10deg = np.empty((0, 1))
lift_coefficient_list_10deg = np.empty((0, 1))
induced_drag_coefficient_list_10deg = np.empty((0, 1))
pitching_moment_quarterchord_list_10deg = np.empty((0, 1))

# ------------ XFLR5 Data -----------
# Open The Text File
script_directory = os.path.dirname(os.path.abspath(__file__))
file_path_0deg = os.path.join(script_directory, '0 deg wing A6.txt')
with open(file_path_0deg, 'r') as file:
    lines_0deg = file.readlines()
file_path_10deg = os.path.join(script_directory, '10 deg wing A6.txt')
with open(file_path_10deg, 'r') as file:
    lines_10deg = file.readlines()
# Append the values to the Lists
start_line = 21
for i in range(19):
    line = lines_0deg[start_line + i + 19].split('   ')

    y_list_0deg = np.append(y_list_0deg, float(line[1]))
    chord_list_0deg = np.append(chord_list_0deg, float(line[2]))
    induced_angle_list_0deg = np.append(induced_angle_list_0deg, float(line[3]))
    lift_coefficient_list_0deg = np.append(lift_coefficient_list_0deg, float(line[4]))
    induced_drag_coefficient_list_0deg = np.append(induced_drag_coefficient_list_0deg, float(line[6]))
    pitching_moment_quarterchord_list_0deg = np.append(pitching_moment_quarterchord_list_0deg, float(line[8]))

    line = lines_10deg[start_line + i + 19].split('   ')

    y_list_10deg = np.append(y_list_10deg, float(line[1]))
    chord_list_10deg = np.append(chord_list_10deg, float(line[2]))
    induced_angle_list_10deg = np.append(induced_angle_list_10deg, float(line[3]))
    lift_coefficient_list_10deg = np.append(lift_coefficient_list_10deg, float(line[4]))
    induced_drag_coefficient_list_10deg = np.append(induced_drag_coefficient_list_10deg, float(line[6]))
    pitching_moment_quarterchord_list_10deg = np.append(pitching_moment_quarterchord_list_10deg, float(line[8]))

chord_function_0deg = sp.interpolate.interp1d(y_list_0deg, chord_list_0deg, kind='cubic', fill_value="extrapolate")
chord_function_10deg = sp.interpolate.interp1d(y_list_10deg, chord_list_10deg, kind='cubic', fill_value="extrapolate")
CL_function_0deg = sp.interpolate.interp1d(y_list_0deg, lift_coefficient_list_0deg, kind='cubic', fill_value="extrapolate")
CL_function_10deg = sp.interpolate.interp1d(y_list_10deg, lift_coefficient_list_10deg, kind='cubic', fill_value="extrapolate")
CD_function_0deg = sp.interpolate.interp1d(y_list_0deg, induced_drag_coefficient_list_0deg, kind='cubic', fill_value="extrapolate")
CD_function_10deg = sp.interpolate.interp1d(y_list_10deg, induced_drag_coefficient_list_10deg, kind='cubic', fill_value="extrapolate")
CM_function_0deg = sp.interpolate.interp1d(y_list_0deg, pitching_moment_quarterchord_list_0deg, kind='cubic', fill_value="extrapolate")
CM_function_10deg = sp.interpolate.interp1d(y_list_10deg, pitching_moment_quarterchord_list_10deg, kind='cubic', fill_value="extrapolate")
induced_AOA_function_0deg = sp.interpolate.interp1d(y_list_0deg, induced_angle_list_0deg, kind='cubic', fill_value="extrapolate")
induced_AOA_function_10deg = sp.interpolate.interp1d(y_list_10deg, induced_angle_list_10deg, kind='cubic', fill_value="extrapolate")

b = 48.82
A = 8.65
e = 0.8
g = 9.80556
n_max = 2.5
n_min = -1
S = 275.5

CL0_0deg = 0.355876
CL0_10deg = 1.141277
CD0_0deg = 0.004713
CD0_10deg = 0.046930
CM0_0deg = -0.486076
CM0_10deg = -1.38817
CD0 = 0.019827

WOE = 90261.9
ZFW = 109215
MTOW = 170382
Weights = np.array([WOE, ZFW, MTOW])*g

rho_s = 1.225
rho_cr = 0.441
Densities = np.array([rho_s, rho_cr])

def CL_des_func(z): # Cl per unit length
    return CL_function_0deg(z) + ((CL_des - CL0_0deg)/(CL0_10deg - CL0_0deg))*(CL_function_10deg(z)-CL_function_0deg(z))


# def CD_des_func(z): # Cd per unit length
#     return CD_function_0deg(z) + ((CL_des - CL0_0deg)/(CL0_10deg - CL0_0deg))**2 * (CD_function_10deg(z)-CD_function_0deg(z))

def CD_des_func(z): # Cd per unit length using weird math I did
    return CD_function_0deg(z) + (2*CL_des)/(np.pi*A*e) * ((CL_des - CL0_0deg)/(CD0_10deg - CD0_0deg)) * (CD_function_10deg(z)-CD_function_0deg(z))


# def CD_des_func(z): # Cd per unit length using drag polar
#     CD_des = CD0 + CL_des_func(z)**2/(np.pi * A * e)
#     return CD_function_0deg(z) + ((CD_des - CD0_0deg)/(CD0_10deg - CD0_0deg)) * (CD_function_10deg(z)-CD_function_0deg(z))


#CHECK!!! From rsearchgate "Lift coefficient change vs pitching moment coefficient"
# dCL = -4dCM
def CM_des_func(z):
    return CM_function_0deg(z) - 0.25 * ((CL_des - CL0_0deg)/(CM0_10deg - CM0_0deg)) * (CM_function_10deg(z) - CM_function_0deg(z))


def Lift_func(z): # Lift per unit lenght (L')
    return 0.5 * rho * V**2 * chord_function_0deg(z) * CL_des_func(z)


def Drag_func(z): # Drag per unit lenght (D')
    return 0.5 * rho * V**2 * chord_function_0deg(z) * CD_des_func(z)

def Moment_func(z):
    return 0.5 * rho * V**2 * chord_function_0deg(z)**2 * CM_des_func(z)


def alpha_func(CL_des):
    sin_alpha = ((CL_des - CL0_0deg)/(CL0_10deg - CL0_0deg))*np.sin(np.radians(10))
    return np.arcsin(sin_alpha)


# fig, axs = plt.subplots(2, 2, figsize=(11, 5), layout='constrained')
# Velocities1 = np.linspace(138, 193, 100)
# Velocities2 = np.linspace(87, 155, 100)
#
# i = 0
# j = 2
#
# for rho in Densities:
#         Drag1_array = np.empty((0, 1))
#         Drag2_array = np.empty((0, 1))
#         Lift_array = np.empty((0, 1))
#         CLdes1_array = np.empty((0, 1))
#         CLdes2_array = np.empty((0, 1))
#
#         print(i)
#         for V in Velocities1:
#             CL_des = n_max*Weights[2]/(0.5 * rho * V**2 * S)
#             Lift = sp.integrate.quad(Lift_func, 0, b / 2)[0] * 2
#             Drag1 = sp.integrate.quad(Drag_func, 0, b / 2)[0] * 2
#             Drag1_array = np.append(Drag1_array, Drag1)
#             Lift_array = np.append(Lift_array, Lift)
#             CLdes1_array = np.append(CLdes1_array, CL_des)
#
#         for V in Velocities2:
#             CL_des = n_min*Weights[2]/(0.5 * rho * V**2 * S)
#             # Lift = sp.integrate.quad(Lift_func, 0, b / 2)[0] * 2
#             # print(np.degrees(alpha_func(CL_des)))
#             Drag2 = sp.integrate.quad(Drag_func, 0, b / 2)[0] * 2
#             Drag2_array = np.append(Drag2_array, Drag2)
#             CLdes2_array = np.append(CLdes2_array, CL_des)
#
#         try:
#             ValidClind1 = np.where(CLdes1_array<1.8)[0][0]
#         except IndexError:
#             ValidClind1 = -1
#         try:
#             ValidClind2 = np.where(CLdes2_array>-0.5)[0][0]
#         except IndexError:
#             ValidClind2 = -1
#
#
#         axs.flat[i].plot(Velocities1, Drag1_array) #label = f'Density = {rho}, n = {n_max}')
#         if ValidClind1 != -1:
#             axs.flat[i].plot(Velocities1[ValidClind1], Drag1_array[ValidClind1], 'ro', label = f'Drag = {round(Drag1_array[ValidClind1])} [N] CL = {round(CLdes1_array[ValidClind1],1)}')
#         # axs.flat[i].set_title(f'Weight = {round(Weights[2])},  CL_des = {round(CLdes1_array[0], 3)}')
#         if rho == 1.225:
#             axs.flat[i].set_title(f'Drag to Velocity at sea level and n = {n_max}')
#         else:
#             axs.flat[i].set_title(f'Drag to Velocity at cruise altitude and n = {n_max}')
#         axs.flat[i].legend()
#
#
#         axs.flat[j].plot(Velocities2, Drag2_array) #, label = f'Density = {rho}, n = {n_min}')
#         if ValidClind2 != -1:
#             axs.flat[j].plot(Velocities2[ValidClind2], Drag2_array[ValidClind2], 'ro', label = f'Drag = {round(Drag2_array[ValidClind2])} [N] CL = {round(CLdes2_array[ValidClind2],1)}')
#             print("Velocity", Velocities2[ValidClind2])
#         # axs.flat[j].set_title(f'Weight = {round(Weights[2])}, CL_des = {round(CLdes2_array[0], 3)}')
#         if rho == 1.225:
#             axs.flat[j].set_title(f'Drag to Velocity at sea level and n = {n_min}')
#         else:
#             axs.flat[j].set_title(f'Drag to Velocity at cruise altitude and n = {n_min}')
#         axs.flat[j].legend()
#
#         for ax in axs.flat:
#             ax.set(xlabel='Velocity [m/s]')
#             ax.set(ylabel='Drag [N]')
#             ax.ticklabel_format(style='sci', axis='y', scilimits=(6, 6))
#
#         i = i+1
#         j = j+1
#         print('CLdes', CLdes1_array[0])
#         print('alpha', np.degrees(alpha_func(CLdes1_array[0])))
#         print('Weight', Weights[2], 'Density', rho)
#         print('Lift', Lift_array[0])
#         print('L/D', Lift_array[0]/max(Drag1_array))
#         print(Lift_array[0]/(Weights[2]*n_max))
#         print('--------------------------')
#
# plt.show()

#Critical case values (to lazy to automate this)

# Which n are we using, change here!!!!!!!!!!!!
n = 2.5

W = Weights[2]
rho = Densities[0]

if n == 2.5:
    V = 138
if n == -1:
    V = 141.262626262

CL_des = n*W/(0.5 * rho * V**2 * S)

# # V = Velocities1[0]
# V = 141.26262626
# W = Weights[2]
# rho = Densities[0]
# # CL_des = n_max*W/(0.5 * rho * V**2 * S)

alpha = np.degrees(alpha_func((CL_des)))

# print(alpha)
#
# print(CL_des)
# print(np.degrees(alpha_func(CL_des)))



