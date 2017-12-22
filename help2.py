from sympy import atan2
from sympy.matrices import Matrix
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2, acos, asin
import numpy as np
from numpy.linalg import inv
# from numpy import cos, sin

import math

def rot_alpha(rotation_matrix):
    r11 = rotation_matrix[0, 0]
    r12 = rotation_matrix[0, 1]
    r13 = rotation_matrix[0, 2]

    r21 = rotation_matrix[1, 0]
    r22 = rotation_matrix[1, 1]
    r23 = rotation_matrix[1, 2]

    r31 = rotation_matrix[2, 0]
    r32 = rotation_matrix[2, 1]
    r33 = rotation_matrix[2, 2]

    alpha_rotation = atan2(r21, r11)

    return alpha_rotation


def rot_beta(rotation_matrix):
    r11 = rotation_matrix[0, 0]
    r12 = rotation_matrix[0, 1]
    r13 = rotation_matrix[0, 2]

    r21 = rotation_matrix[1, 0]
    r22 = rotation_matrix[1, 1]
    r23 = rotation_matrix[1, 2]

    r31 = rotation_matrix[2, 0]
    r32 = rotation_matrix[2, 1]
    r33 = rotation_matrix[2, 2]

    beta_rotation = atan2(-r31, sqrt((r11 * r11) + (r21 * r21)))

    return beta_rotation


def rot_gamma(rotation_matrix):
    r11 = rotation_matrix[0, 0]
    r12 = rotation_matrix[0, 1]
    r13 = rotation_matrix[0, 2]

    r21 = rotation_matrix[1, 0]
    r22 = rotation_matrix[1, 1]
    r23 = rotation_matrix[1, 2]

    r31 = rotation_matrix[2, 0]
    r32 = rotation_matrix[2, 1]
    r33 = rotation_matrix[2, 2]

    gamma_rotation = atan2(r32, r33)

    return gamma_rotation

def rotate_x(rads):
    rotated = Matrix([
        [1, 0, 0, 0],
        [0, cos(rads), -sin(rads), 0],
        [0, sin(rads), cos(rads), 0],
        [0, 0, 0, 1]
    ])

    # rotated = np.asarray([
    #     [1, 0, 0, 0],
    #     [0, cos(rads), -sin(rads), 0],
    #     [0, sin(rads), cos(rads), 0],
    #     [0, 0, 0, 1]
    # ])

    return rotated

def rotate_z(rads):

    rotated = Matrix([
        [cos(rads), -sin(rads), 0, 0],
        [sin(rads), cos(rads), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    #
    # rotated = np.asarray([
    #     [cos(rads), -sin(rads), 0, 0],
    #     [sin(rads), cos(rads), 0, 0],
    #     [0, 0, 1, 0],
    #     [0, 0, 0, 1]
    # ])

    return rotated

def rotate_y(rads):
    rotated = Matrix([
        [cos(rads), 0, sin(rads), 0],
        [0, 1, 0, 0],
        [-sin(rads), 0, cos(rads), 0],
        [0, 0, 0, 1]
    ])

    # rotated = np.asarray([
    #     [cos(rads), 0, sin(rads), 0],
    #     [0, 1, 0, 0],
    #     [-sin(rads), 0, cos(rads), 0],
    #     [0, 0, 0, 1]
    # ])

    return rotated


def generate_rrpy_matrix(roll, pitch, yaw):
    Rot_X = rotate_x(roll)[0:3, 0:3]
    Rot_Y = rotate_y(pitch)[0:3, 0:3]
    Rot_Z = rotate_z(yaw)[0:3, 0:3]

    # Rrpy = Rot_Z * Rot_Y * Rot_X
    Rrpy = Rot_Z * Rot_Y * Rot_X

    return Rrpy


def quarternion2rotation(q):
    # Normalize quarternion in case it is not nomralized
    q = q / np.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2 + q[3] ** 2)
    q0 = q[3]
    qx = q[0]
    qy = q[1]
    qz = q[2]

    R = np.array([[1 - 2 * (qy ** 2 + qz ** 2), 2 * (qx * qy - q0 * qz), 2 * (q0 * qy + qx * qy)],
                  [2 * (qx * qy + q0 * qz), 1 - 2 * (qx ** 2 + qy ** 2), 2 * (qy * qz - q0 * qx)],
                  [2 * (qx * qz - q0 * qy), 2 * (q0 * qx + qy * qz), 1 - 2 * (qx ** 2 + qy ** 2)]])
    return R

# correct values for the different shelves
correct_values = {
    "shelf_3": {
        'joint_1': -0.44293154231873455,
        'joint_2': 0.5386775324322572,
        'joint_3': 0.16542712045273245,
        'joint_4': -0.6314924490169878,
        'joint_5': -0.8115364743375029,
        'joint_6': 0.46657506948980654,
        "wx": 1.897,
        "wy": -0.89986,
        "wz": 0.8108
    },

    # correct values for shelf 4
    "shelf_4": {
        'joint_1': 0.4430142917290958,
        'joint_2': 0.2410005481971762,
        'joint_3': -0.018890976127179826,
        'joint_4': 1.1347509560756794,
        'joint_5': -0.4926981628581206,
        'joint_6': -1.0838895504725423,
        "wx": 1.897,
        "wy": 0.89994,
        "wz": 1.5808
    },

    # correct values for shelf 7
    "shelf_7": {
        'joint_1': 0.4429112339320298,
        'joint_2': 0.2405030554081513,
        'joint_3': -0.5330847563438992,
        'joint_4': -1.0234459000219864,
        'joint_5': 0.5247176752176532,
        'joint_6': 0.9579212364217895,
        "wx": 1.897,
        "wy": 0.89982,
        "wz": 2.345
    },

    # correct values for shelf 9
    "shelf_9": {
        'joint_1': -0.44294603521585785,
        'joint_2': 0.2405436964408425,
        'joint_3': -0.5332512472834132,
        'joint_4': 1.0228278416487209,
        'joint_5': 0.5252679315532607,
        'joint_6': -0.956034634397418,
        "wx": 1.8969,
        "wy": -0.89988,
        "wz": 2.3451
    }

}

given_values = {
    # shelf 3 given values
    "shelf_3": {
        'roll': -0.00014556016151799498,
        'pitch': 0.000606189206810885,
        'yaw': -0.000422844112202551,

        'px': 2.2,
        'py': -0.89999,
        'pz': 0.81098
    },

    # shelf 4 given values
    "shelf_4": {
        'roll': 6.226668809756833e-06,
        'pitch': -0.0004581857291535871,
        'yaw': -0.0006602558977705618,

        'px': 2.0900268663882255,
        'py': 0.9000621976697936,
        'pz': 1.5810337949811457
    },

    # shelf 7 given values
    "shelf_7": {
        'roll': 0.000660472914656566,
        'pitch': -1.613033347887964e-05,
        'yaw': 0.0008755346172453829,

        "px": 2.2,
        "py": 0.90008,
        "pz": 2.345
    },

    # shelf 9 given values
    "shelf_9": {
        'roll': 0.0006378107154360105,
        'pitch': -0.0004715761181352515,
        'yaw': -0.000639290366401144,

        "px": 2.1999,
        "py": -0.90008,
        "pz": 2.345
    }
}
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

# Homogeneous transforms
T0_1 = Matrix([[cos(q1), -sin(q1), 0, a0],
               [sin(q1) * cos(alpha0), cos(q1) * cos(alpha0), -sin(alpha0), -sin(alpha0) * d1],
               [sin(q1) * sin(alpha0), cos(q1) * sin(alpha0), cos(alpha0), cos(alpha0) * d1],
               [0, 0, 0, 1]])

T0_1 = T0_1.subs(s)

T1_2 = Matrix([[cos(q2), -sin(q2), 0, a1],
               [sin(q2) * cos(alpha1), cos(q2) * cos(alpha1), -sin(alpha1), -sin(alpha1) * d2],
               [sin(q2) * sin(alpha1), cos(q2) * sin(alpha1), cos(alpha1), cos(alpha1) * d2],
               [0, 0, 0, 1]])

T1_2 = T1_2.subs(s)

T2_3 = Matrix([[cos(q3), -sin(q3), 0, a2],
               [sin(q3) * cos(alpha2), cos(q3) * cos(alpha2), -sin(alpha2), -sin(alpha2) * d3],
               [sin(q3) * sin(alpha2), cos(q3) * sin(alpha2), cos(alpha2), cos(alpha2) * d3],
               [0, 0, 0, 1]])

T2_3 = T2_3.subs(s)

T0_2 = T0_1 * T1_2
T0_3 = T0_2 * T2_3

for shelf in given_values:
    print shelf
    theta1 = correct_values[shelf]["joint_1"]
    theta2 = correct_values[shelf]["joint_2"]
    theta3 = correct_values[shelf]["joint_3"]
    
    roll = given_values[shelf]["roll"]
    pitch = given_values[shelf]["pitch"]
    yaw = given_values[shelf]["yaw"]
    
    R0_3 = simplify(T0_3[0:3, 0:3])
    
    # print(R0_3)
    
    R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
    # convert sympy matrix to numpy matrix to avoid errors
    # how to convert lifted from: https://stackoverflow.com/a/37491889
    R0_3 = np.array(R0_3).astype(np.float64)

    Rrpy = (generate_rrpy_matrix(roll, pitch, yaw) * rotate_y(pi/2)[0:3, 0:3] * rotate_z(pi)[0:3, 0:3])
    R3_6 = (inv(R0_3)) * Rrpy
    
    # theta4 = atan2(R3_6[2, 1], R3_6[2, 2])
    theta4 = atan2(R3_6[2,2], -R3_6[0,2])
    # theta4 = atan2(R3_6[2,2], 1-R3_6[0,2])
    theta5 = acos(R3_6[1, 2])
    # theta5 = atan2(-R3_6[2, 0], sqrt((R3_6[0, 0] * R3_6[0, 0]) + (R3_6[1, 0] * R3_6[1, 0])))
    # theta6c = atan2(R3_6[1, 0], R3_6[0, 0])
    # theta6c = atan2(R3_6[1,2], -R3_6[0,2])

    # correct solution for theta6
    theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])
    # solution is near, but not as close
    # theta6b = atan2(-R3_6[2,2], R3_6[1,0])

    # if the value of theta4 is greater than or less than pi, then reverse the rotation of theta4
    # by subtracting it from pi, and reversing the sign. Also reverse the rotations of theta5 and 6
    if theta4 < -pi/2:
        theta4 = -(-pi - theta4)
        theta5 = -theta5
        theta6 = -(pi - theta6)
    elif theta4 > pi/2:
        theta4 = -(pi - theta4)
        theta5 = -theta5
        theta6 = -(-pi - theta6)

    expected_theta4 = correct_values[shelf]["joint_4"]
    expected_theta5 = correct_values[shelf]["joint_5"]
    expected_theta6 = correct_values[shelf]["joint_6"]
    
    print("expected theta4 ", expected_theta4)
    print("computed theta4 ", theta4.evalf())
    print("expected theta5 ", expected_theta5)
    print("computed theta5 ", theta5.evalf())
    print("expected theta6 ", expected_theta6)
    print("computed theta6 ", theta6.evalf())

    # function I used to permute through
    # for i in range(3):
    #     for j in range(3):
    #         for k in range(3):
    #             for l in range(3):
    #                 mytheta = atan2(-R3_6[i,j], R3_6[k,l])
    #                 # print mytheta
    #                 if abs(round(mytheta, 1)) == abs(round(expected_theta6, 1)):
    #                     print mytheta
    #                     print(i, j, k, l, "\n")


