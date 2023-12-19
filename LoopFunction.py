import numpy as np
import matplotlib.pyplot as plt

from InternalForcesAE2111 import z_values, b
from DeflectionComputation import Icalculation, geometry, enclosed_area, Stress_Analysis
from StressComputation import tau_ave, tau_cr, torque_flow_stress, sigma_cr, normal_stress, column_stress, failure_test

rho = 2700
#Iterables
nr_of_stringers = 8
nr_of_ribs = 20

def Weight_func(geometry_lst):
    total_vol = 0
    for i in range(len(z_values)):
        total_vol += area_lst[i][4]*geometry_lst[i][7]
    total_vol += 2*10**(-3)*b/2*nr_of_stringers
    rib_spacing = b/2/(nr_of_ribs-1)
    for rib in range(nr_of_ribs):
        z_loc = rib_spacing*rib
        t_root = geometry_lst[int(z_loc*10)][7]
        area_rib = area_lst[int(z_loc*10)][3]
        total_vol += area_rib*t_root
    return total_vol*rho

W_min = 10**12

spanwisesplit = 0.9
t_ratio = 2.1

# spanwisesplit = 0.75
# while spanwisesplit <= 1:
#     print(spanwisesplit, nr_of_ribs, nr_of_stringers)
#     t_ratio = 1
#     Failure = True
#     while Failure and t_ratio<25:
        # Compute geometry values for all z_values at once
geometry_lst = np.array([geometry(z, t_ratio, spanwisesplit) for z in z_values])
area_lst = np.array([enclosed_area(z, geometry_lst, nr_of_stringers) for z in z_values])

ivalues, Ixx_lst, Iyy_lst, Ixy_lst = [], [], [], []
stress_point_lst = []
for z in z_values:
    ivalue, Ixx, Iyy, Ixy, Stress_points = Icalculation(z, geometry_lst, area_lst, nr_of_stringers, spanwisesplit)

    # Append values to respective lists
    ivalues.append(ivalue)
    Ixx_lst.append(Ixx)
    Iyy_lst.append(Iyy)
    Ixy_lst.append(Ixy)
    stress_point_lst.append(Stress_points)
# ivalues, Ixx_lst, Iyy_lst, Ixy_lst, stress_point_lst = Icalculation(z_values, geometry_lst, area_lst, nr_of_stringers, spanwisesplit)


# Convert lists to NumPy arrays if needed
ivalues = np.array(ivalues)
Ixx_lst = np.array(Ixx_lst)
Iyy_lst = np.array(Iyy_lst)
Ixy_lst = np.array(Ixy_lst)


num_points = 4  # Number of stress points
stress_arrays = [[] for _ in range(num_points)]

for i in range(len(z_values)):
    stresses = [Stress_Analysis(z_values[i], point, stress_point_lst, Ixx_lst, Iyy_lst, Ixy_lst) for point in range(num_points)]
    for idx, stress in enumerate(stresses):
        stress_arrays[idx].append(stress)

# for i in range(len(z_values)):
#     stresses = [Stress_Analysis(z_values[i], point, stress_point_lst[i], Ixx_lst[i], Iyy_lst[i], Ixy_lst[i]) for
#                 point in range(num_points)]
#     for idx, stress in enumerate(stresses):
#         stress_arrays[idx].append(stress)

Point1Stress = abs(np.array(stress_arrays[0]))
Point2Stress = abs(np.array(stress_arrays[1]))
Point3Stress = abs(np.array(stress_arrays[2]))
Point4Stress = abs(np.array(stress_arrays[3]))

tau_ave1_values = np.empty(len(z_values))
tau_ave2_values = np.empty(len(z_values))
tau_cr1_values = np.empty(len(z_values))
tau_cr2_values = np.empty(len(z_values))
tau_cr3_values = np.empty(len(z_values))
sigma_cr_values = np.empty(len(z_values))
sigma_cr_br = np.empty(len(z_values))
yield_stress = np.array([276e6]*len(z_values))

for i in range(len(z_values)):
    tau_ave1_values[i] = tau_ave(z_values[i], geometry_lst) - torque_flow_stress(z_values[i], geometry_lst, area_lst)
    tau_ave2_values[i] = tau_ave(z_values[i], geometry_lst) + torque_flow_stress(z_values[i], geometry_lst, area_lst)

    tau_cr1, tau_cr2, tau_cr3 = tau_cr(z_values[i], geometry_lst, nr_of_ribs)
    tau_cr1_values[i] = tau_cr1
    tau_cr2_values[i] = tau_cr2
    tau_cr3_values[i] = tau_cr3

    sigma_cr_values[i] = sigma_cr(z_values[i], geometry_lst, nr_of_ribs, nr_of_stringers)
    sigma_cr_br[i] = column_stress(nr_of_ribs)
    Point1Stress[i] = Point1Stress[i] - normal_stress(z_values[i], area_lst)
    Point2Stress[i] = Point2Stress[i] - normal_stress(z_values[i], area_lst)
    Point3Stress[i] = Point3Stress[i] - normal_stress(z_values[i], area_lst)
    Point4Stress[i] = Point4Stress[i] - normal_stress(z_values[i], area_lst)


margin_web_1 = tau_cr1_values/tau_ave1_values
margin_web_2 = tau_cr2_values/tau_ave2_values
margin_web_3 = tau_cr3_values/tau_ave2_values

print(tau_cr1_values)
print(tau_ave1_values)

margin_top_plate = sigma_cr_values/Point4Stress
margin_column = sigma_cr_br/Point4Stress
margin_yield = yield_stress/Point1Stress



# Create subplots
fig, axs = plt.subplots(2, 3, figsize=(13, 5), constrained_layout=True)

# Data for plotting grouped together
# plot_data = [
#     (axs[0, 0], z_values, tau_ave1_values, 'red'),
#     (axs[0, 0], z_values, tau_cr1_values, 'blue'),
#     (axs[0, 1], z_values, tau_ave2_values, 'red'),
#     (axs[0, 1], z_values, tau_cr2_values, 'blue'),
#     (axs[0, 2], z_values, tau_cr3_values, 'blue'),
#     (axs[0, 2], z_values, tau_ave2_values, 'red'),
#     # (axs[1, 0], z_values, Point1Stress, 'red'),
#     # (axs[1, 0], z_values, Point2Stress, 'red'),
#     # (axs[1, 0], z_values, Point3Stress, 'red'),
#     (axs[1, 0], z_values, Point4Stress, 'red'),
#     (axs[1, 0], z_values, sigma_cr_values, 'blue'),
#     (axs[1, 1], z_values, sigma_cr_br, 'blue'),
#     (axs[1, 1], z_values, Point4Stress, 'red'),
#     (axs[1, 2], z_values, Point1Stress, 'red'),
#     # (axs[1, 2], z_values, Point2Stress, 'red'),
#     # (axs[1, 2], z_values, Point3Stress, 'red'),
#     # (axs[1, 2], z_values, Point4Stress, 'red'),
#     (axs[1, 2], z_values, yield_stress, 'blue')]
plot_data = [
    (axs[0, 0], z_values, margin_web_1, 'blue'),
    (axs[0, 1], z_values, margin_web_2, 'blue'),
    (axs[0, 2], z_values, margin_web_3, 'blue'),
    (axs[1, 0], z_values, margin_top_plate, 'blue'),
    (axs[1, 1], z_values, margin_column, 'blue'),
    (axs[1, 2], z_values, margin_yield, 'blue')]

# Plot all data at once
for ax, x_data, y_data, color in plot_data:
    ax.plot(x_data, y_data, color=color)
    ax.set_ylim(0, 10)

plt.show()



    #     Failure = failure_test(tau_ave1_values, tau_ave2_values, tau_cr1_values, tau_cr2_values, tau_cr3_values,
    #                      Point4Stress, Point3Stress, sigma_cr_values, sigma_cr_br)
    #     # print(Failure, t_ratio, spanwisesplit)
    #     if Failure:
    #         t_ratio += 0.01
    #     else:
    #         W = Weight_func(geometry_lst)
    #         print("W", W, t_ratio, spanwisesplit, nr_of_ribs, nr_of_stringers)
    #         if W < W_min:
    #             W_min = W
    #             best_lst = [t_ratio, spanwisesplit, nr_of_ribs, nr_of_stringers, W]
    #             print('best', best_lst)
    # spanwisesplit += 0.05
    # #     nr_of_ribs += 1
    # # nr_of_stringers += 4

# print(best_lst)

print(Ixx_lst[0])
