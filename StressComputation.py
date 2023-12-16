import numpy as np
import scipy
from sympy import symbols, integrate
from InternalForcesAE2111 import momentfuncyz, torquefunc, momentfuncxz
import matplotlib.pyplot as plt

# q_design_op = int(input("Select design 1, 2 or 3"))
q_design_op = 1
if q_design_op == 1 :
    a_ratio = 0.0841
    p_ratio = 0.127
    d_ratio = 0.01
    t_ratio = 1
    k_ratio = 0.5
    spanwisesplit = 0.9#Ratio of span
    nr_of_stringers = 30
    #Change the numbers here
elif q_design_op == 2:
    a_ratio = 0.0841
    p_ratio = 0.0841
    d_ratio = 0
    t_ratio = 1.2
    k_ratio = 0
    spanwisesplit = 0#Ratio of span
    nr_of_stringers = 20
    #Change the numbers here
elif q_design_op == 3:
    a_ratio = 0.0841
    p_ratio = 0.127
    d_ratio = 0.01
    t_ratio = 1.4
    k_ratio = 0.5
    spanwisesplit = 0.5#Ratio of span
    nr_of_stringers = 4
    #Change the numbers here


b = 48.82
span = 48.82



def Icalculation(z):
    r = 0.5 * (8.487 - 0.233 * z)
    c = (8.487 - 0.233 * z)
    a = a_ratio * c
    p = p_ratio * c
    d = d_ratio * c
    t = t_ratio * (0.003125 * r - 0.00235)
    #print(t)
    k = k_ratio * r
    G = 26e9
    E = 68.9e9
    #T = 1
    e = p-d-a
    A_1 = a*r
    A_2 = 1/2*r*d
    A_3 = 1/2*r*e
    A_stringer = 2*10**-3
    #print(t)
    #centroid A_1
    C_1_x = 1/2*r
    C_1_z = 1/2*a+d

    #centroid A_2
    C_2_x = r/3
    C_2_z = 2/3*d

    #centroid A_3
    C_3_x = r/3
    C_3_z = a+d+e/3

    #centroid x and z
    C_t_x = (A_1*C_1_x+A_2*C_2_x+A_3*C_3_x)/(A_1+A_2+A_3)
    C_t_z = (A_1*C_1_z+A_2*C_2_z+A_3*C_3_z)/(A_1+A_2+A_3)
    #print(C_t_x, C_t_z)
    #angles triangles
    theta = np.arctan(d/r)
    beta = np.arctan(e/r)

    #ipotenusa
    h_1 = np.sqrt(d**2+r**2)
    h_2 = np.sqrt(e**2+r**2)

    # Stress_Points = np.empty((0,2))
    # Point 1, Bottom Left
    x_1 = - C_t_x
    y_1 = - C_t_z
    # Stress_Points = np.append(Stress_Points, (np.array([x_1, y_1])).reshape(2,1))
    # Point 2, Bottom right
    x_2 = r - C_t_x
    y_2 = e - C_t_z
    # Stress_Points = np.append(Stress_Points, (np.array([x_2, y_2])).reshape(2,1))
    # Point 3, Top left
    x_3 = - C_t_x
    y_3 = p - C_t_z
    # Stress_Points = np.append(Stress_Points, (np.array([x_3, y_3])).reshape(2,1))
    # Point 4, Top right
    x_4 = r - C_t_x
    y_4 = e + a - C_t_z
    # Stress_Points = np.append(Stress_Points, (np.array([x_4, y_4])).reshape(2,1))
    Stress_Points = [[x_1, y_1], [x_2, y_2], [x_3, y_3], [x_4, y_4]]




    #if Multicell == False:
    if z > spanwisesplit*span/2:
        A = 1/2 * (a + p) * r
        Divider = nr_of_stringers / 2 + 1 # half the actual nr of stringers
        # stringer calculations
        I_xx_Stringer = 0
        I_yy_Stringer = 0
        I_xy_Stringer = 0
        for i in range(int(Divider)):
            distbetween = r / Divider
            distfromleftwall = (i + 1) * distbetween
            hstringer = C_t_z - distfromleftwall * np.tan(theta)
            I_xx_Stringer = I_xx_Stringer + A_stringer * hstringer ** 2
            I_xy_Stringer = I_xy_Stringer + A_stringer * -hstringer * (-C_t_x + distfromleftwall)
            hstringer = p - C_t_z - distfromleftwall * np.tan(beta)
            I_xx_Stringer = I_xx_Stringer + A_stringer * hstringer ** 2
            I_yy_Stringer = I_yy_Stringer + A_stringer * 2 * (C_t_x - distfromleftwall) ** 2
            I_xy_Stringer = I_xy_Stringer + A_stringer * hstringer * (-C_t_x + distfromleftwall)


        # centroids of walls
        # upper wall
        I_xx_1 = (t * h_1 ** 3 * np.sin(-beta) ** 2) / 12 + (t * h_1) * (p - d / 2 - C_t_z) ** 2
        # print("I_xx_1:", I_xx_1)
        I_yy_1 = (t * h_1 ** 3 * np.cos(-beta) ** 2) / 12 + (t * h_1) * (r / 2 - C_t_x) ** 2
        # print("I_yy_1:", I_yy_1)
        I_xy_1 = (t * h_1 ** 3 * np.sin(-beta) * np.cos(-beta)) / 12 + (t * h_1) * (p - d / 2 - C_t_z) * (r / 2 - C_t_x)
        # print("I_xy_1:", I_xy_1)
        # right wall
        I_xx_2 = (t * a ** 3) / 12 + (t * a) * (e + a / 2 - C_t_z) ** 2
        # print("I_xx_2:", I_xx_2)
        I_yy_2 = (t * a) * (r - C_t_x) ** 2
        # print("I_yy_2:", I_yy_2)
        I_xy_2 = t * a * (r - C_t_x) * (e + a / 2 - C_t_z)
        # print("I_xy_2:", I_xy_2)
        # bottom wall
        I_xx_3 = (t * h_2 ** 3 * np.sin(theta) ** 2) / 12 + (t * h_2) * (e / 2 - C_t_z) ** 2
        # print("I_xx_3:", I_xx_3)
        I_yy_3 = (t * h_2 ** 3 * np.cos(theta) ** 2) / 12 + (t * h_2) * (r / 2 - C_t_x) ** 2
        # print("I_yy_3:", I_yy_3)
        I_xy_3 = (t * h_2 ** 3 * np.sin(theta) * np.cos(theta)) / 12 + (t * h_2) * (e / 2 - C_t_z) * (r / 2 - C_t_x)
        # print("I_xy_3:", I_xy_3)
        # left wall
        I_xx_4 = (t * p ** 3) / 12 + (t * p) * (p / 2 - C_t_z) ** 2
        # print("I_xx_4:", I_xx_4)
        I_yy_4 = (t * p) * C_t_x ** 2
        # print("I_yy_4:", I_yy_4)
        I_xy_4 = (t * p) * (- C_t_x) * (p / 2 - C_t_z)
        # print("I_xy_4:", I_xy_4)
        # Total
        I_xx = I_xx_1 + I_xx_2 + I_xx_3 + I_xx_4 + I_xx_Stringer
        I_yy = I_yy_1 + I_yy_2 + I_yy_3 + I_yy_4 + I_yy_Stringer
        I_xy = I_xy_1 + I_xy_2 + I_xy_3 + I_xy_4 + I_xy_Stringer

    #if Multicell == True:
    else:
        DividerA1 = np.ceil(nr_of_stringers / 2) + 1
        DividerA2 = np.floor(nr_of_stringers / 2) + 1
        # stringer calculations
        I_xx_Stringer = 0
        I_yy_Stringer = 0
        I_xy_Stringer = 0
        for i in range(int(DividerA1)):
            distbetween = k / DividerA1
            distfromleftwall = (i + 1) * distbetween
            hstringer = C_t_z - distfromleftwall * np.tan(theta)
            I_xx_Stringer = I_xx_Stringer + A_stringer * hstringer ** 2
            I_xy_Stringer = I_xy_Stringer + A_stringer * -hstringer * (-C_t_x + distfromleftwall)
            hstringer = p - C_t_z - distfromleftwall * np.tan(beta)
            I_xx_Stringer = I_xx_Stringer + A_stringer * hstringer ** 2
            I_yy_Stringer = I_yy_Stringer + A_stringer * 2 * (C_t_x - distfromleftwall) ** 2
            I_xy_Stringer = I_xy_Stringer + A_stringer * hstringer * (-C_t_x + distfromleftwall)
        for i in range(int(DividerA2)):
            distbetween = (r - k) / DividerA2
            distfromleftwall = (i + 1) * distbetween
            hstringer = C_t_z - distfromleftwall * np.tan(theta)
            I_xx_Stringer = I_xx_Stringer + A_stringer * hstringer ** 2
            I_xy_Stringer = I_xy_Stringer + A_stringer * -hstringer * (-C_t_x + distfromleftwall)
            hstringer = p- C_t_z - distfromleftwall * np.tan(beta)
            I_xx_Stringer = I_xx_Stringer + A_stringer * hstringer ** 2
            I_yy_Stringer = I_yy_Stringer + A_stringer * 2 * (C_t_x - distfromleftwall) ** 2
            I_xy_Stringer = I_xy_Stringer + A_stringer * hstringer * (-C_t_x + distfromleftwall)

        # length of middle bar
        l = a + (r - k) * np.tan(beta) + (r - k) * np.tan(theta)
        # centroids of walls
        # upper wall
        I_xx_1 = (t * h_1 ** 3 * np.sin(-beta) ** 2) / 12 + (t * h_1) * (p - d / 2 - C_t_z) ** 2
        # print("I_xx_1:", I_xx_1)
        I_yy_1 = (t * h_1 ** 3 * np.cos(-beta) ** 2) / 12 + (t * h_1) * (r / 2 - C_t_x) ** 2
        # print("I_yy_1:", I_yy_1)
        I_xy_1 = (t * h_1 ** 3 * np.sin(-beta) * np.cos(-beta)) / 12 + (t * h_1) * (p - d / 2 - C_t_z) * (
                    r / 2 - C_t_x)
        # print("I_xy_1:", I_xy_1)
        # right wall
        I_xx_2 = (t * a ** 3) / 12 + (t * a) * (e + a / 2 - C_t_z) ** 2
        # print("I_xx_2:", I_xx_2)
        I_yy_2 = (t * a) * (r - C_t_x) ** 2
        # print("I_yy_2:", I_yy_2)
        I_xy_2 = t * a * (r - C_t_x) * (e + a / 2 - C_t_z)
        # print("I_xy_2:", I_xy_2)
        # bottom wall
        I_xx_3 = (t * h_2 ** 3 * np.sin(theta) ** 2) / 12 + (t * h_2) * (e / 2 - C_t_z) ** 2
        # print("I_xx_3:", I_xx_3)
        I_yy_3 = (t * h_2 ** 3 * np.cos(theta) ** 2) / 12 + (t * h_2) * (r / 2 - C_t_x) ** 2
        # print("I_yy_3:", I_yy_3)
        I_xy_3 = (t * h_2 ** 3 * np.sin(theta) * np.cos(theta)) / 12 + (t * h_2) * (e / 2 - C_t_z) * (
                    r / 2 - C_t_x)
        # print("I_xy_3:", I_xy_3)
        # left wall
        I_xx_4 = (t * p ** 3) / 12 + (t * p) * (p / 2 - C_t_z) ** 2
        # print("I_xx_4:", I_xx_4)
        I_yy_4 = (t * p) * C_t_x ** 2
        # print("I_yy_4:", I_yy_4)
        I_xy_4 = (t * p) * (- C_t_x) * (p / 2 - C_t_z)
        # print("I_xy_4:", I_xy_4)
        # mid wall
        I_xx_5 = (t * l ** 3) / 12 + (t * l) * (l / 2 + e / 2 - C_t_z) ** 2
        # print("I_xx_5:", I_xx_5)
        I_yy_5 = (t * l) * (k - C_t_x) ** 2
        # print("I_yy_5:", I_yy_5)
        I_xy_5 = t * l * (e / 2 + l / 2 - C_t_z) * (k - C_t_x)
        # print("I_xy_5:", I_xy_5)

        # Total
        I_xx = I_xx_1 + I_xx_2 + I_xx_3 + I_xx_4 + I_xx_5 + I_xx_Stringer
        I_yy = I_yy_1 + I_yy_2 + I_yy_3 + I_yy_4 + I_yy_5 + I_yy_Stringer
        I_xy = I_xy_1 + I_xy_2 + I_xy_3 + I_xy_4 + I_xy_5 + I_xy_Stringer

    return -I_xx*E, I_xx, I_yy, I_xy, Stress_Points


def Jcalculation(z):
    r = 0.5 * (8.487 - 0.233 * z)
    a = a_ratio * r
    p = p_ratio * r
    d = d_ratio * r
    t = t_ratio*(0.003125 * r - 0.00235)
    k = k_ratio * r
    G = 26e9
    E = 68.9e9
    #Multicell = True
    T = 1
    e = p - d - a
    #angles triangles
    theta = np.arctan(d/r)
    beta = np.arctan(e/r)

    #ipotenusa
    h_1 = np.sqrt(d**2+r**2)
    h_2 = np.sqrt(e**2+r**2)
    #if not Multicell:
        #print("now single cell")

    if z > spanwisesplit*span/2:
        A = 1/2 * (a + p) * r

        J = 4 * A ** 2 * t / (a + p + h_1 + h_2)
        #print(J)
    else:
        # length of middle bar

        l = a + (r - k) * np.tan(beta) + (r - k) * np.tan(theta)
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
    ind = np.where(z_values == z)[0][0]
    new_z = [j for j in z_values if j <= z]
    new_moment = momentyzvalues[:ind+1]
    new_ivalues = ivalues[:ind+1]
    EI_values = np.divide(new_moment, new_ivalues)
    return scipy.integrate.simps(EI_values, new_z)

def deflect_function(z):
    ind = np.where(z_values == z)[0][0]
    new_z = [j for j in z_values if j <= z]
    new_anglevalues = angle_values[:ind + 1]
    return scipy.integrate.simps(new_anglevalues, new_z)

def twist_function(z):
    ind = np.where(z_values == z)[0][0]
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

z_values = np.arange(0, b/2 + 0.1, 0.1)
momentofinertia = np.empty((0,1))
torsionalstiffness = np.empty((0,1))
ivalues = np.empty((0, 1))
angle_values = np.empty((0, 1))
momentyzvalues = np.empty((0, 1))
deflect_values = np.empty((0, 1))
torque_values = np.empty((0, 1))
jvalues = np.empty((0, 1))
twist_values = np.empty((0, 1))
for i in range(len(z_values)):
    ivalues = np.append(ivalues, Icalculation(z_values[i])[0])
    momentyzvalues = np.append(momentyzvalues, momentfuncyz(z_values[i]))
    momentofinertia = np.append(momentofinertia, Icalculation(z_values[i])[1])

for i in range(len(z_values)):
    angle_values = np.append(angle_values, angle_function(z_values[i]))

for i in range(len(z_values)):
    deflect_values = np.append(deflect_values, deflect_function(z_values[i]))
# print(deflect_values)

for i in range(len(z_values)):
    jvalues = np.append(jvalues, Jcalculation(z_values[i])[0])
    torque_values = np.append(torque_values, torquefunc(z_values[i]))
    torsionalstiffness = np.append(torsionalstiffness, Jcalculation(z_values[i])[1])

for i in range(len(z_values)):
    twist_values = np.append(twist_values, twist_function(z_values[i]))
# print(twist_values)
# fig, axs = plt.subplots(1, 2, figsize=(10, 5), layout='constrained')

"""
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


# Below this are computations for tensile stress, some can also be used for compression as well
# Analysis of points bottom left: 1, bottom right: 2, top left: 3, top right: 4
#print(Icalculation((10))[4])
def Stress_Analysis(z, point):
    return (((momentfuncyz(z) * Icalculation(z)[2] - momentfuncxz(z) * Icalculation(z)[3]) * Icalculation(z)[4][int(point)][1] + (momentfuncxz(z) * Icalculation(z)[1] - momentfuncyz(z) * Icalculation(z)[3]) * Icalculation(z)[4][int(point)][0]) / (Icalculation(z)[1] * Icalculation(z)[2] - Icalculation(z)[3] ** 2))


Point1Stress = np.empty((0, 1))
Point2Stress = np.empty((0, 1))
Point3Stress = np.empty((0, 1))
Point4Stress = np.empty((0, 1))


for i in range(len(z_values)):
    Point1Stress = np.append(Point1Stress, Stress_Analysis(z_values[i], 0))
    Point2Stress = np.append(Point2Stress, Stress_Analysis(z_values[i], 1))
    Point3Stress = np.append(Point3Stress, Stress_Analysis(z_values[i], 2))
    Point4Stress = np.append(Point4Stress, Stress_Analysis(z_values[i], 3))


def principal_axis(z):
    return np.degrees(-np.arctan((momentfuncxz(z) * Icalculation(z)[1] - momentfuncyz(z) * Icalculation(z)[3] ) / (momentfuncyz(z) * Icalculation(z)[2] - momentfuncxz(z) * Icalculation(z)[3])))

principal_axis_values = np.empty((0, 1))

for i in range(len(z_values)):
    principal_axis_values = np.append(principal_axis_values, principal_axis(z_values[i]))

def Moment_Angle(z):
    return np.degrees(np.arctan(momentfuncxz(z)/momentfuncyz(z)))

Moment_Angle_values = np.empty((0, 1))

for i in range(len(z_values)):
    Moment_Angle_values = np.append(Moment_Angle_values, Moment_Angle(z_values[i]))

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

positive_load = [251231315.16698426, 54296088.62229689, -221555531.49262363, -258785594.87434104]
negative_load = [-84791413.52927507, -36782737.48649467, 96954942.68921526, 85982636.67852694]
#Stress at points
"""
fig, axs = plt.subplots(2)
axs[1].plot(z_values, Point1Stress, label = "Left")
axs[1].plot(z_values, Point2Stress, label = "Right")
axs[1].set_title('Bottom Points')
plt.legend()
axs[0].plot(z_values, Point3Stress, label = "Left")
axs[0].plot(z_values, Point4Stress, label = "Right")
axs[0].set_title('Top Points')
plt.legend()
plt.show()
"""
# Cyclic loading calculations
stress_max_tension = 1
A = 4.3378e-7
m = 2.6183
K1c = 29e6
c_min = 1.27e-3
positive_load = np.array(positive_load)
negative_load = np.array(negative_load)
delta_sigma = max(abs(positive_load - negative_load))
c_crit = K1c ** 2 / np.pi / delta_sigma ** 2
delta_sigma = delta_sigma/ 10e6
c = c_min
"""
cycles = 0
plotx = []
ploty = []
while c <= c_crit:
    cycles += 1
    delta_K = 1.1 * delta_sigma * np.sqrt(np.pi * c)
    dcdN = A * delta_K ** m
    c += dcdN / 10e3
    plotx.append(cycles)
    ploty.append(c)
    if cycles % 100000 == 0:
        print(cycles, c, c_crit)
print(cycles)
plt.plot(plotx, ploty)
plt.show()
"""
Point1Stress = abs(Point1Stress)
Point2Stress = abs(Point2Stress)
Point3Stress = abs(Point3Stress)
Point4Stress = abs(Point4Stress)

margin_of_safety = []
for z in range(len(z_values)):
    margin = 310e6 / max(Point1Stress[z], Point2Stress[z], Point3Stress[z], Point4Stress[z])
    if margin < 10:
        margin_of_safety.append(margin)
    else:
        margin_of_safety.append(10)
# print(margin_of_safety[0])
#plt.plot(z_values, margin_of_safety)
#plt.show()


# plt.plot(z_values, Point3Stress)
# plt.show()

# print(Point1Stress)
# print(Point4Stress)
