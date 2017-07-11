import numpy as np
from sympy import symbols, pi, simplify, sqrt, atan2, acos

from IK import RADS_AT_REST_JOINT3, RADS_AT_REST_JOINT2, get_roll, get_pitch, get_yaw
from IK import generate_rrpy_matrix, get_wrist_coordinates, get_theta_3, get_alpha, get_beta, get_theta_2, \
     T0_3, rotate_z, rotate_y, T_total
from sample_data import given_values, correct_values

import os.path, json


# create symbols for variables
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

# DH Parameters

s = {alpha0: 0, a0: 0, d1: 0.75,
     alpha1: -pi / 2, a1: 0.35, d2: 0, q2: q2 - pi / 2,
     alpha2: 0, a2: 1.25, d3: 0,
     alpha3: -pi / 2, a3: -0.054, d4: 1.50,
     alpha4: pi / 2, a4: 0, d5: 0,
     alpha5: -pi / 2, a5: 0, d6: 0,
     alpha6: 0, a6: 0, d7: 0.303, q7: 0
     }

# open the saved error list, if file does not exist, create it
# load the json from the error list
# open error list from disk
error_list_path = "./error_list.json"

if os.path.isfile(error_list_path):
    with open(error_list_path, "r") as error_list_file:
        json_error_list = json.loads(error_list_file.read(), encoding='utf-8')
        px_error_list = json_error_list["px"]
        py_error_list = json_error_list["py"]
        pz_error_list = json_error_list["pz"]
        roll_error_list = json_error_list["roll"]
        pitch_error_list = json_error_list["pitch"]
        yaw_error_list = json_error_list["yaw"]
else:
    # if file does not exist, create an empty file and initialize lists
    with open(error_list_path, "w") as error_list_file:
        json_error_list = {
            "px": [],
            "py": [],
            "pz": [],
            "roll": [],
            "pitch": [],
            "yaw": []
        }

        px_error_list = json_error_list["px"]
        py_error_list = json_error_list["py"]
        pz_error_list = json_error_list["pz"]
        roll_error_list = json_error_list["roll"]
        pitch_error_list = json_error_list["pitch"]
        yaw_error_list = json_error_list["yaw"]

        json.dump(json_error_list, error_list_file)

for shelf in given_values:
    print(shelf)
    roll = given_values[shelf]["roll"]
    pitch = given_values[shelf]["pitch"]
    yaw = given_values[shelf]["yaw"]

    Rrpy = generate_rrpy_matrix(roll, pitch, yaw)

    px = given_values[shelf]["px"]
    py = given_values[shelf]["py"]
    pz = given_values[shelf]["pz"]

    # should be s[d7] + s[d6] or .303, however .193 works better
    dist_wrist_to_effector = 0.303

    wx, wy, wz = get_wrist_coordinates(Rrpy, px, py, pz, dist_wrist_to_effector)

    expected_wx = correct_values[shelf]["wx"]
    expected_wy = correct_values[shelf]["wy"]
    expected_wz = correct_values[shelf]["wz"]
    # assert (round(wx, 1) == round(expected_wx, 1)), (str(shelf) + " wx mismatch " + wx + " " + expected_wx)
    # assert (round(wy, 1) == round(expected_wy, 1)), (str(shelf) + " wy mismatch " + wy + " " + expected_wy)
    # assert (round(wz, 1) == round(expected_wz, 1)), (str(shelf) + " wz mismatch " + wz + " " + expected_wz)

    if round(wx, 1) != round(expected_wx, 1):
        print(str(shelf) + " wx mismatch ", wx, expected_wx)
    if round(wy, 1) != round(expected_wy, 1):
        print(str(shelf) + " wy mismatch ", wy, expected_wy)
    if round(wz, 1) != round(expected_wz, 1):
        print(str(shelf) + ' wz mismatch ', wz, expected_wz)

    # once we're done checking if we have sane wrist values, we recompute to use .193 for theta computation
    # dist_wrist_to_effector = 0.193
    wx, wy, wz = get_wrist_coordinates(Rrpy, px, py, pz, dist_wrist_to_effector)

    # compute theta for joint 1
    theta1 = atan2(wy, wx)
    theta1 = theta1.evalf()
    expected_theta1 = correct_values[shelf]["joint_1"]
    # print("expected theta1", expected_theta1)
    # print("theta1         ", theta1)
    # assert round(theta1, 2) == round(expected_theta1, 2), "theta1 mismatch"

    if round(theta1, 2) != round(expected_theta1, 2):
        print("theta1 mismatch", theta1, expected_theta1)

    # compute theta3
    # compute the side adjacents of theta 3
    distance_joint_2to3 = s[a2]
    distance_joint_3to5 = sqrt((s[a3] ** 2) + (s[d4] ** 2))

    # get the offsets to the origin of joint 2
    x_offset_to_joint2 = s[a1]
    z_offset_to_joint2 = s[d1]

    # compute for the 2 possible values of theta 3
    unadjusted_theta_3 = get_theta_3(distance_joint_2to3, distance_joint_3to5, wx, wz, theta1, x_offset_to_joint2,
                                     z_offset_to_joint2)

    # choose the first of the two theta_3 results, which is the elbow up result
    theta3 = unadjusted_theta_3 - RADS_AT_REST_JOINT3
    expected_theta3 = correct_values[shelf]["joint_3"]
    # assert round(theta3, 2) == round(expected_theta3, 2), "theta3 mismatch"

    if round(theta3, 2) != round(expected_theta3, 2):
        print("theta3 mismatch", theta3, expected_theta3)

    # compute the parts used for theta 2
    alpha = get_alpha(distance_joint_2to3, distance_joint_3to5, wx, wz, theta1, joint_2_x_offset=s[a1],
                      joint_2_z_offset=s[d1])
    beta = get_beta(wx, wz, theta1, s[a1], s[d1])

    unadjusted_theta_2 = get_theta_2(alpha, beta).evalf()

    # compute theta 2 value
    theta2 = RADS_AT_REST_JOINT2 - unadjusted_theta_2
    theta2 = theta2.evalf()
    expected_theta2 = correct_values[shelf]["joint_2"]
    # assert round(theta2, 2) == round(expected_theta2, 2), "theta3 mismatch"

    if round(theta2, 2) != round(expected_theta2, 2):
        print("theta2 mismatch", theta2, expected_theta2)

    theta1 = correct_values[shelf]["joint_1"]
    theta2 = correct_values[shelf]["joint_2"]
    theta3 = correct_values[shelf]["joint_3"]

    roll = given_values[shelf]["roll"]
    pitch = given_values[shelf]["pitch"]
    yaw = given_values[shelf]["yaw"]

    R0_3 = simplify(T0_3[0:3, 0:3])
    R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
    # convert sympy matrix to numpy matrix to avoid errors
    # how to convert lifted from: https://stackoverflow.com/a/37491889
    R0_3 = np.array(R0_3).astype(np.float64)

    Rrpy = (generate_rrpy_matrix(roll, pitch, yaw) * rotate_y(pi / 2)[0:3, 0:3] * rotate_z(pi)[0:3, 0:3])
    R3_6 = (np.linalg.inv(R0_3)) * Rrpy

    theta4 = atan2(R3_6[2, 2], -R3_6[0, 2])
    theta5 = acos(R3_6[1, 2])
    # correct solution for theta6
    theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])
    # solution is near, but not as close
    # theta6b = atan2(-R3_6[2,2], R3_6[1,0])

    # if the value of theta4 is greater than or less than pi, then reverse the rotation of theta4
    # by subtracting it from pi, and reversing the sign. Also reverse the rotations of theta5 and 6
    # if theta4 < -pi / 2:
    #     theta4 = -(-pi - theta4)
    #     theta5 = -theta5
    #     theta6 = -(pi - theta6)
    # elif theta4 > pi / 2:
    #     theta4 = -(pi - theta4)
    #     theta5 = -theta5
    #     theta6 = -(-pi - theta6)

    expected_theta4 = correct_values[shelf]["joint_4"]
    expected_theta5 = correct_values[shelf]["joint_5"]
    expected_theta6 = correct_values[shelf]["joint_6"]

    if round(theta4, 2) != round(expected_theta4, 2):
        print("theta4 mismatch", theta4, expected_theta4)

    if round(theta5, 2) != round(expected_theta5, 2):
        print("theta5 mismatch", theta5, expected_theta5)

    if round(theta6, 2) != round(expected_theta6, 2):
        print("theta6 mismatch", theta6, expected_theta6)

    FK_matrix = T_total.evalf(
        subs={q1: theta1, q2: theta2, q3: theta3, q4: theta4, q5: theta5, q6: theta6})

    FK_px = FK_matrix[0, 3]
    FK_py = FK_matrix[1, 3]
    FK_pz = FK_matrix[2, 3]
    # print("expected px  ", given_values[shelf]["px"])
    # print("FK_px        ", FK_px)
    # print("expected py  ", given_values[shelf]["py"])
    # print("FK_py        ", FK_py)
    # print("expected pz  ", given_values[shelf]["pz"])
    # print("FK_pz        ", FK_pz)

    FK_roll = get_roll(FK_matrix)
    FK_pitch = get_pitch(FK_matrix)
    FK_yaw = get_yaw(FK_matrix)
    # print("expected roll    ", roll)
    # print("FK_roll          ", FK_roll)
    # print("expected pitch   ", pitch)
    # print("FK_pitch         ", FK_pitch)
    # print("expected yaw     ", yaw)
    # print("FK_yaw           ", FK_yaw)

    # compute the error values for the px, py, pz and rpy
    px_error = (px - FK_px) # ** 2
    py_error = (py - FK_py) # ** 2
    pz_error = (pz - FK_pz) #** 2

    roll_error = (roll - FK_roll)  # ** 2
    pitch_error = (pitch - FK_pitch)  #** 2
    yaw_error = (yaw - FK_yaw) # ** 2

    # append each of the error values to their respective lists
    px_error_list.append(float(px_error))
    py_error_list.append(float(py_error))
    pz_error_list.append(float(pz_error))
    roll_error_list.append(float(roll_error))
    pitch_error_list.append(float(pitch_error))
    yaw_error_list.append(float(yaw_error))

# once all the operations have been done, append the new set of errors to disk
with open(error_list_path, "w") as error_list_file:
    json.dump(json_error_list, error_list_file, encoding='utf-8')
