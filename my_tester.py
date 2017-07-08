from sample_data import given_values, correct_values
from IK import generate_rrpy_matrix, get_wrist_coordinates, get_theta_3
from IK import RADS_AT_REST_JOINT3

from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2, acos

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

for shelf in given_values:
    print(shelf)
    roll = given_values[shelf]["roll"]
    pitch = given_values[shelf]["pitch"]
    yaw = given_values[shelf]["yaw"]

    Rrpy = generate_rrpy_matrix(roll, pitch, yaw)

    px = given_values[shelf]["px"]
    py = given_values[shelf]["py"]
    pz = given_values[shelf]["pz"]

    # should be s[d7] + s[d6] or .303, however .197 works better
    dist_wrist_to_effector = 0.303

    wx, wy, wz = get_wrist_coordinates(Rrpy, px, py, pz, dist_wrist_to_effector)

    expected_wx = correct_values[shelf]["wx"]
    expected_wy = correct_values[shelf]["wy"]
    expected_wz = correct_values[shelf]["wz"]

    # print("expected wx ", expected_wx)
    # print("wx          ", wx)
    # print("expected wy ", expected_wy)
    # print("wy          ", wy)
    # print("expected wz ", expected_wz)
    # print("wz          ", wz)
    # assert (round(wx, 1) == round(expected_wx, 1)), (str(shelf) + " wx mismatch " + wx + " " + expected_wx)
    # assert (round(wy, 1) == round(expected_wy, 1)), (str(shelf) + " wy mismatch " + wy + " " + expected_wy)
    # assert (round(wz, 1) == round(expected_wz, 1)), (str(shelf) + " wz mismatch " + wz + " " + expected_wz)

    if round(wx, 1) != round(expected_wx, 1):
        print(str(shelf) + " wx mismatch ",  wx, expected_wx)
    if round(wy, 1) != round(expected_wy, 1):
        print(str(shelf) + " wy mismatch ",  wy, expected_wy)
    if round(wz, 1) != round(expected_wz, 1):
        print(str(shelf) + ' wz mismatch ', wz, expected_wz)

    # compute theta for joint 1
    theta1 = atan2(wy, wx)
    theta1 = theta1.evalf()
    expected_theta1 = correct_values[shelf]["joint_1"]
    # print("expected theta1", expected_theta1)
    # print("theta1         ", theta1)
    # assert round(theta1, 2) == round(expected_theta1, 2), "theta1 mismatch"

    if round(theta1, 2) == round(expected_theta1, 2):
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
    theta3 = unadjusted_theta_3[0].evalf() - RADS_AT_REST_JOINT3
    expected_theta3 = correct_values[shelf]["joint_3"]
    print("expected theta3", expected_theta3)
    print("theta3         ", theta3)
