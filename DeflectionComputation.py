import numpy as np
import scipy
from sympy import symbols, integrate
from InternalForcesAE2111 import z_values, momentyz_values, momentxz_values, momentxy_values as torque_values
import matplotlib.pyplot as plt

np.seterr(divide='ignore', invalid='ignore')
b = 48.82

nr_of_stringers = 10
nr_of_ribs = 10
spanwisesplit = 0.5
t_ratio = 1

q_design_op = 1
if q_design_op == 1:
    a_ratio = 0.0841
    p_ratio = 0.127
    d_ratio = 0.01
    # t_ratio = 1.85
    k_ratio = 0.5

    # spanwisesplit = 0.85#Ratio of span
    # nr_of_stringers = 12
    # nr_of_ribs = 22
    A_str = 2*10**(-3)

    #Change the numbers here
elif q_design_op == 2:
    a_ratio = 0.0841
    p_ratio = 0.0841
    d_ratio = 0
    # t_ratio = 1.2
    k_ratio = 0
    # spanwisesplit = 0#Ratio of span
    # nr_of_stringers = 20
    # nr_of_ribs = 22
    A_str = 2*10**(-3)

    #Change the numbers here
elif q_design_op == 3:
    a_ratio = 0.0841
    p_ratio = 0.127
    d_ratio = 0.01
    # t_ratio = 1.4
    k_ratio = 0.5
    spanwisesplit = 0.5#Ratio of span
    # nr_of_stringers = 4
    # nr_of_ribs = 22
    A_str = 2*10**(-3)

    #Change the numbers here


def geometry(z):
    c = (8.487 - 0.233 * z)
    r = 0.5 * c
    b_1 = c * p_ratio
    b_2 = c * a_ratio
    d = c * d_ratio
    t = t_ratio*(0.003125 * r - 0.00235)
    k = k_ratio * r
    e = b_1 - d - b_2
    theta = np.arctan(d / r)
    beta = np.arctan(e / r)

    if q_design_op == 2 or z > (b / 2)*spanwisesplit:
        b_3 = 0
    else:
        b_3 = b_2 + (r - k) * np.tan(beta) + (r - k) * np.tan(theta)

    h_1 = np.sqrt(d ** 2 + r ** 2)
    h_2 = np.sqrt(e ** 2 + r ** 2)
    return r, b_1, b_2, b_3, k, d, e, t, theta, beta, h_1, h_2


def enclosed_area(z):
    ind = int(z*10)
    r, a, d, e = geometry_lst[ind][0], geometry_lst[ind][2], geometry_lst[ind][5], geometry_lst[ind][6]
    A_1 = a*r
    A_2 = 1/2*r*d
    A_3 = 1/2*r*e
    A = A_3 + A_2 + A_1
    return A_1, A_2, A_3, A


def cross_sec_area(z):
    ind = int(z * 10)
    r, b_1, b_2, b_3, k, d, e, t, theta, beta, h_1, h_2 = geometry_lst[ind]
    A_cross = (b_1 + b_2 + b_3 + h_1 + h_2)*t + nr_of_stringers*A_str
    return A_cross


def Icalculation(z):
    ind = int(z*10)
    [r, p, a, l, k, d, e, t, theta, beta, h_1, h_2] = geometry_lst[ind]

    G = 26e9
    E = 68.9e9

    # A_1 = a * r
    # A_2 = 0.5 * r * d
    # A_3 = 0.5 * r * e
    A_stringer = 2 * 10**-3
    A_1, A_2, A_3, A_tot = area_lst[ind][0:4]

    # Centroids
    C_1_x = 0.5 * r
    C_1_z = 0.5 * a + d

    C_2_x = r / 3
    C_2_z = 2 / 3 * d

    C_3_x = r / 3
    C_3_z = a + d + e / 3

    C_t_x = (A_1 * C_1_x + A_2 * C_2_x + A_3 * C_3_x) / (A_1 + A_2 + A_3)
    C_t_z = (A_1 * C_1_z + A_2 * C_2_z + A_3 * C_3_z) / (A_1 + A_2 + A_3)


    Stress_Points = [[-C_t_x, -C_t_z],
        [r - C_t_x, e - C_t_z],
        [-C_t_x, p - C_t_z],
        [r - C_t_x, e + a - C_t_z]]

    if z > spanwisesplit * b / 2:
        A = 0.5 * (a + p) * r
        Divider = nr_of_stringers / 2 + 1  # half the actual nr of stringers

        # Stringer calculations using vectorized operations
        dist_between = r / Divider
        dist_from_left_wall = np.arange(1, int(Divider) + 1) * dist_between
        h_stringer_theta = C_t_z - dist_from_left_wall * np.tan(theta)
        h_stringer_beta = p - C_t_z - dist_from_left_wall * np.tan(beta)

        I_xx_Stringer = A_stringer * (np.sum(h_stringer_theta ** 2) + np.sum(h_stringer_beta ** 2))
        I_xy_Stringer = A_stringer * (np.sum(-h_stringer_theta * (-C_t_x + dist_from_left_wall)) +
                                      np.sum(h_stringer_beta * (-C_t_x + dist_from_left_wall)))
        I_yy_Stringer = A_stringer * (np.sum(2 * (C_t_x - dist_from_left_wall) ** 2))

        # Wall calculations
        I_xx_walls = [
            (t * h_1 ** 3 * np.sin(-beta) ** 2) / 12 + (t * h_1) * (p - d / 2 - C_t_z) ** 2,
            (t * a ** 3) / 12 + (t * a) * (e + a / 2 - C_t_z) ** 2,
            (t * h_2 ** 3 * np.sin(theta) ** 2) / 12 + (t * h_2) * (e / 2 - C_t_z) ** 2,
            (t * p ** 3) / 12 + (t * p) * (p / 2 - C_t_z) ** 2
        ]

        I_yy_walls = [
            (t * h_1 ** 3 * np.cos(-beta) ** 2) / 12 + (t * h_1) * (r / 2 - C_t_x) ** 2,
            (t * a) * (r - C_t_x) ** 2,
            (t * h_2 ** 3 * np.cos(theta) ** 2) / 12 + (t * h_2) * (r / 2 - C_t_x) ** 2,
            (t * p) * C_t_x ** 2
        ]

        I_xy_walls = [
            (t * h_1 ** 3 * np.sin(-beta) * np.cos(-beta)) / 12 + (t * h_1) * (p - d / 2 - C_t_z) * (r / 2 - C_t_x),
            t * a * (r - C_t_x) * (e + a / 2 - C_t_z),
            (t * h_2 ** 3 * np.sin(theta) * np.cos(theta)) / 12 + (t * h_2) * (e / 2 - C_t_z) * (r / 2 - C_t_x),
            (t * p) * (-C_t_x) * (p / 2 - C_t_z)
        ]

        # Total moments of inertia
        I_xx = np.sum(I_xx_walls) + I_xx_Stringer
        I_yy = np.sum(I_yy_walls) + I_yy_Stringer
        I_xy = np.sum(I_xy_walls) + I_xy_Stringer
    else:
        DividerA1 = np.ceil(nr_of_stringers / 2) + 1
        DividerA2 = np.floor(nr_of_stringers / 2) + 1

        # Initialize stringer variables
        I_xx_Stringer = 0
        I_yy_Stringer = 0
        I_xy_Stringer = 0

        # Loop for stringer calculations for DividerA1
        for i in range(int(DividerA1)):
            dist_between = k / DividerA1
            dist_from_left_wall = (i + 1) * dist_between

            # Calculate stringer heights
            hstringer_theta = C_t_z - dist_from_left_wall * np.tan(theta)
            hstringer_beta = p - C_t_z - dist_from_left_wall * np.tan(beta)

            # Update stringer moments of inertia
            I_xx_Stringer += A_stringer * (hstringer_theta ** 2 + hstringer_beta ** 2)
            I_yy_Stringer += A_stringer * 2 * (C_t_x - dist_from_left_wall) ** 2
            I_xy_Stringer += A_stringer * (-hstringer_theta * (-C_t_x + dist_from_left_wall) + hstringer_beta * (
                        -C_t_x + dist_from_left_wall))

        # Loop for stringer calculations for DividerA2
        for i in range(int(DividerA2)):
            dist_between = (r - k) / DividerA2
            dist_from_left_wall = (i + 1) * dist_between

            # Calculate stringer heights
            hstringer_theta = C_t_z - dist_from_left_wall * np.tan(theta)
            hstringer_beta = p - C_t_z - dist_from_left_wall * np.tan(beta)

            # Update stringer moments of inertia
            I_xx_Stringer += A_stringer * (hstringer_theta ** 2 + hstringer_beta ** 2)
            I_yy_Stringer += A_stringer * 2 * (C_t_x - dist_from_left_wall) ** 2
            I_xy_Stringer += A_stringer * (-hstringer_theta * (-C_t_x + dist_from_left_wall) + hstringer_beta * (
                        -C_t_x + dist_from_left_wall))


        # Calculate walls moments of inertia
        I_xx_walls = [
            (t * h_1 ** 3 * np.sin(-beta) ** 2) / 12 + (t * h_1) * (p - d / 2 - C_t_z) ** 2,
            (t * a ** 3) / 12 + (t * a) * (e + a / 2 - C_t_z) ** 2,
            (t * h_2 ** 3 * np.sin(theta) ** 2) / 12 + (t * h_2) * (e / 2 - C_t_z) ** 2,
            (t * p ** 3) / 12 + (t * p) * (p / 2 - C_t_z) ** 2,
            (t * l ** 3) / 12 + (t * l) * (l / 2 + e / 2 - C_t_z) ** 2
        ]

        I_yy_walls = [
            (t * h_1 ** 3 * np.cos(-beta) ** 2) / 12 + (t * h_1) * (r / 2 - C_t_x) ** 2,
            (t * a) * (r - C_t_x) ** 2,
            (t * h_2 ** 3 * np.cos(theta) ** 2) / 12 + (t * h_2) * (r / 2 - C_t_x) ** 2,
            (t * p) * C_t_x ** 2,
            (t * l) * (k - C_t_x) ** 2
        ]

        I_xy_walls = [
            (t * h_1 ** 3 * np.sin(-beta) * np.cos(-beta)) / 12 + (t * h_1) * (p - d / 2 - C_t_z) * (r / 2 - C_t_x),
            t * a * (r - C_t_x) * (e + a / 2 - C_t_z),
            (t * h_2 ** 3 * np.sin(theta) * np.cos(theta)) / 12 + (t * h_2) * (e / 2 - C_t_z) * (r / 2 - C_t_x),
            (t * p) * (-C_t_x) * (p / 2 - C_t_z),
            t * l * (e / 2 + l / 2 - C_t_z) * (k - C_t_x)
        ]

        # Calculate total moments of inertia
        I_xx = np.sum(I_xx_walls) + I_xx_Stringer
        I_yy = np.sum(I_yy_walls) + I_yy_Stringer
        I_xy = np.sum(I_xy_walls) + I_xy_Stringer

    return -I_xx * E, I_xx, I_yy, I_xy, Stress_Points



def Jcalculation(z):
    ind = int(z*10)
    [r, p, a, l, k, d, e, t, theta, beta, h_1, h_2] = geometry_lst[ind]
    G = 26e9
    E = 68.9e9
    #Multicell = True
    T = 1
    #if not Multicell:
        #print("now single cell")

    if z > spanwisesplit*b/2:
        A = 1/2 * (a + p) * r

        J = 4 * A ** 2 * t / (a + p + h_1 + h_2)
        #print(J)
    else:
        # length of middle bar

        #Area of Sections with multicell
        Area_1 = 0.5 * (p + l) * k
        Area_2 = 0.5 * (l + a) * (r - k)

        #print(Area_1, Area_2)
        # System of equations to find q1 q2 and d0/dy
        #Perimeter_1 = b + l + h_1 - (k / np.cos(theta))+ h_2 - (k / np.cos(beta))
        #Perimeter_2 = l + a + (k / np.cos(theta)) + (k / np.cos(beta))
        Perimeter_1 = p + l + h_1*(k/r) + h_2*(k/r)
        Perimeter_2 = l + a + h_1*(1 - k/r) + h_2*(1 - k/r)
        A = np.array([[2*Area_1, 2*Area_2, 0],[Perimeter_1/(t * G), -l/(t * G), -2*Area_1],[-l/(t * G), Perimeter_2/ (t * G), -2*Area_2]])
        B = np.array([T, 0, 0])

        solution = np.linalg.solve(A, B)

        #print(f"q1 is {solution[0]} q2 is {solution[1]} d0/dy is {solution[2]}")

        J = T/(solution[2] * G)
        #print(J)
    return J * G, J

# Moment and Torque functions
#M = input("Moment Function")
#T = input("Torque Function")
# v prime prime
#v = M / (E * I_xx)
# phi prime
#phi = T / (G * J)

def angle_function(z):
    ind = int(z * 10)
    new_z = [j for j in z_values if j <= z]
    new_moment = momentyz_values[:ind+1]
    new_ivalues = ivalues[:ind+1]
    EI_values = np.divide(new_moment, new_ivalues)
    return scipy.integrate.simps(EI_values, new_z)


def deflect_function(z):
    ind = int(z * 10)
    new_z = [j for j in z_values if j <= z]
    new_anglevalues = angle_values[:ind + 1]
    return scipy.integrate.simps(new_anglevalues, new_z)


def twist_function(z):
    ind = int(z * 10)
    new_z = [j for j in z_values if j <= z]
    new_torque = torque_values[:ind + 1]
    new_jvalues = jvalues[:ind+1]
    GJ_values = np.divide(new_torque, new_jvalues)
    return scipy.integrate.simps(GJ_values, new_z)

# Define the variable
z = symbols('z')

#Momentfuncyz
# Compute the indefinite integral
# indefinite_integral = integrate(divfunction, z)
# deflection = integrate(indefinite_integral, z)
# print(dblintegral)

geometry_lst = []
area_lst = []
for z in z_values:
    r, b_1, b_2, b_3, k, d, e, t, theta, beta, h_1, h_2 = geometry(z)
    geometry_lst.append([r, b_1, b_2, b_3, k, d, e, t, theta, beta, h_1, h_2])
    A1, A2, A3, A_tot = enclosed_area(z)
    A_cross = cross_sec_area(z)
    area_lst.append([A1, A2, A3, A_tot, A_cross])


ivalues, Ixx_lst, Iyy_lst, Ixy_lst = [], [], [], []
stress_point_lst = []

jvalues = []
torsionalstiffness = []
# angle_values = np.empty((0, 1))
# deflect_values = np.empty((0, 1))
# twist_values = np.empty((0, 1))


for z in z_values:
    ivalue, Ixx, Iyy, Ixy, Stress_points = Icalculation(z)

    # Append values to respective lists
    ivalues.append(ivalue)
    Ixx_lst.append(Ixx)
    Iyy_lst.append(Iyy)
    Ixy_lst.append(Ixy)
    stress_point_lst.append(Stress_points)

    jvalue, torsionalstiffness_value = Jcalculation(z)

    # Append values to respective lists
    jvalues.append(jvalue)
    torsionalstiffness.append(torsionalstiffness_value)

# Convert lists to NumPy arrays if needed
ivalues = np.array(ivalues)
Ixx_lst = np.array(Ixx_lst)
Iyy_lst = np.array(Iyy_lst)
Ixy_lst = np.array(Ixy_lst)
jvalues = np.array(jvalues)
torsionalstiffness = np.array(torsionalstiffness)



# for i in range(len(z_values)):
#     angle_values = np.append(angle_values, angle_function(z_values[i]))
#     deflect_values = np.append(deflect_values, deflect_function(z_values[i]))
#     twist_values = np.append(twist_values, twist_function(z_values[i]))

"""
fig, axs = plt.subplots(1, 2, figsize=(10, 5), layout='constrained')


#axs.flat[0].plot(z_values, angle_values)
#axs.flat[0].set_title('Angle')

axs.flat[0].plot(z_values, deflect_values)
axs.flat[0].set_title('Deflection')
axs.flat[0].set(ylabel=r'Deflection [m]', xlabel='Spanwise Location [m]')

axs.flat[1].plot(z_values, twist_values)
axs.flat[1].set_title('Twist')
axs.flat[1].set(ylabel=r'Wing Twist [rad]', xlabel='Spanwise Location [m]')

axs.flat[0].plot(z_values, momentofinertia)
axs.flat[0].set_title('Moment of Inertia')
axs.flat[0].set(ylabel=r'Moment of Inertia [$m^{4}$]', xlabel='Spanwise Location [m]')

axs.flat[1].plot(z_values, torsionalstiffness)
axs.flat[1].set_title('Torsional Stiffness')
axs.flat[1].set(ylabel=r'Torsional Stiffness [$m^{4}$]', xlabel='Spanwise Location [m]')

plt.show()
"""


def Stress_Analysis(z, point):
    ind = int(z * 10)  # Calculate index directly without np.where
    stress_point = stress_point_lst[ind][int(point)]
    moment_yz = momentyz_values[ind]
    moment_xz = momentxz_values[ind]
    ixx = Ixx_lst[ind]
    iyy = Iyy_lst[ind]
    ixy = Ixy_lst[ind]

    stress = ((moment_yz * iyy - moment_xz * ixy) * stress_point[1]
              + (moment_xz * ixx - moment_yz * ixy) * stress_point[0]) \
             / (ixx * iyy - ixy ** 2)

    return stress

num_points = 4  # Number of stress points
stress_arrays = [[] for _ in range(num_points)]

for i in range(len(z_values)):
    stresses = [Stress_Analysis(z_values[i], point) for point in range(num_points)]
    for idx, stress in enumerate(stresses):
        stress_arrays[idx].append(stress)

Point1Stress = np.array(stress_arrays[0])
Point2Stress = np.array(stress_arrays[1])
Point3Stress = np.array(stress_arrays[2])
Point4Stress = np.array(stress_arrays[3])

########################################################################

# def principal_axis(z):
#     try:
#         result = np.degrees(-np.arctan((momentfuncxz(z) * Icalculation(z)[1] - momentfuncyz(z) * Icalculation(z)[3]) / (
#                     momentfuncyz(z) * Icalculation(z)[2] - momentfuncxz(z) * Icalculation(z)[3])))
#     except RuntimeWarning:
#         result = 0
#     return result
#
# principal_axis_values = np.empty((0, 1))
#
# for i in range(len(z_values)):
#     principal_axis_values = np.append(principal_axis_values, principal_axis(z_values[i]))
#
# def Moment_Angle(z):
#     try:
#         result = np.degrees(np.arctan(momentfuncxz(z)/momentfuncyz(z)))
#     except RuntimeWarning:
#         result = 0
#     return result
#
# Moment_Angle_values = np.empty((0, 1))
#
# for i in range(len(z_values)):
#     Moment_Angle_values = np.append(Moment_Angle_values, Moment_Angle(z_values[i]))

#############################################################################################

# Moment_valuesX = np.empty((0, 1))
# Moment_valuesY = np.empty((0, 1))
#
# for i in range(len(z_values)):
#     Moment_valuesX = np.append(Moment_valuesX, momentfuncyz(z_values[i]))
#     Moment_valuesY = np.append(Moment_valuesY, momentfuncxz(z_values[i]))



#Axis
# fig, axs = plt.subplots(2)
# axs[0].plot(z_values, principal_axis_values)
# axs[1].plot(z_values, Moment_Angle_values)
# axs[2].plot(z_values, Moment_valuesX)
# axs[2].plot(z_values, Moment_valuesY)
# print(max(max(Point1Stress), max(Point2Stress)), min(min(Point3Stress), min(Point4Stress)))
#print(max(Point1Stress), max(Point2Stress), min(Point3Stress), min(Point4Stress))


#########################################################################################
# positive_load = [251231315.16698426, 54296088.62229689, -221555531.49262363, -258785594.87434104]
# negative_load = [-84791413.52927507, -36782737.48649467, 96954942.68921526, 85982636.67852694]
# #Stress at points
# """
# fig, axs = plt.subplots(2)
# axs[1].plot(z_values, Point1Stress, label = "Left")
# axs[1].plot(z_values, Point2Stress, label = "Right")
# axs[1].set_title('Bottom Points')
# plt.legend()
# axs[0].plot(z_values, Point3Stress, label = "Left")
# axs[0].plot(z_values, Point4Stress, label = "Right")
# axs[0].set_title('Top Points')
# plt.legend()
# plt.show()
# """


# # Cyclic loading calculations
# stress_max_tension = 1
# A = 4.3378e-7
# m = 2.6183
# K1c = 29e6
# c_min = 1.27e-3
# positive_load = np.array(positive_load)
# negative_load = np.array(negative_load)
# delta_sigma = max(abs(positive_load - negative_load))
# c_crit = K1c ** 2 / np.pi / delta_sigma ** 2
# delta_sigma = delta_sigma/ 10e6
# c = c_min
# """
# cycles = 0
# plotx = []
# ploty = []
# while c <= c_crit:
#     cycles += 1
#     delta_K = 1.1 * delta_sigma * np.sqrt(np.pi * c)
#     dcdN = A * delta_K ** m
#     c += dcdN / 10e3
#     plotx.append(cycles)
#     ploty.append(c)
#     if cycles % 100000 == 0:
#         print(cycles, c, c_crit)
# print(cycles)
# plt.plot(plotx, ploty)
# plt.show()
# """
##############################################################################



Point1Stress = abs(Point1Stress)
Point2Stress = abs(Point2Stress)
Point3Stress = abs(Point3Stress)
Point4Stress = abs(Point4Stress)

margin_of_safety = []
for z in range(len(z_values)):
    margin = 276e6 / max(Point1Stress[z], Point2Stress[z], Point3Stress[z], Point4Stress[z])
    if margin < 10:
        margin_of_safety.append(margin)
    else:
        margin_of_safety.append(10)



# print(margin_of_safety[0])
plt.plot(z_values, margin_of_safety)
plt.show()

