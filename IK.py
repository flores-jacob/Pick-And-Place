import numpy as np
from numpy import array, ndarray
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2, acos
from sympy.matrices import Matrix

from sample_data import given_values, correct_values

# These are the rads of joint 3 in resting position
RADS_AT_REST_JOINT3 = 1.60678078687695

# These are the rads of joint 2 in resting position
RADS_AT_REST_JOINT2 = pi/2 # 1.57079632679490


SHELF = "shelf_3"


# Law of cosines
def get_theta_3(len_link2_3, len_link3_5, wx, wz, q1, joint_2_x_offset=0, joint_2_z_offset=0):
    # TODO enable arm to bend backwards
    # This solution was adapted from: https://www.youtube.com/watch?v=llUBbpWVPQE time at 4:49
    x_dist_from_joint2 = (wx / cos(q1)) - joint_2_x_offset
    z_dist_from_joint2 = wz - joint_2_z_offset

    numerator = (x_dist_from_joint2 ** 2) + (z_dist_from_joint2 ** 2) - (len_link2_3 ** 2) - (len_link3_5 ** 2)
    denominator = (2 * len_link2_3 * len_link3_5)

    cos_theta = numerator / denominator

    if cos_theta > 1:
        mod_cos_theta = cos_theta - 1
        theta_3_elbow_up = atan2(+sqrt(1 - (mod_cos_theta ** 2)), mod_cos_theta) - (pi / 2)
        theta_3_elbow_down = atan2(-sqrt(1 - (mod_cos_theta ** 2)), mod_cos_theta) - (pi / 2)
    elif cos_theta < 0:
        theta_3_elbow_up = atan2(+sqrt(1 - (cos_theta ** 2)), cos_theta)
        theta_3_elbow_down = atan2(-sqrt(1 - (cos_theta ** 2)), cos_theta)
    else:
        theta_3_elbow_up = atan2(+sqrt(1 - (cos_theta ** 2)), cos_theta)
        theta_3_elbow_down = atan2(-sqrt(1 - (cos_theta ** 2)), cos_theta)

    # maybe we can also try to return the value that is positive, instead of assuming that the elbow_up values is the
    # positive one
    theta_3_elbow_up = theta_3_elbow_up.evalf()
    theta_3_elbow_down = theta_3_elbow_down.evalf()

    return theta_3_elbow_up


def get_beta(wx, wz, q1, a1=0, d1=0):
    x_dist_from_joint2 = (wx / cos(q1)) - a1
    z_dist_from_joint2 = wz - d1

    beta = atan2(z_dist_from_joint2, x_dist_from_joint2)

    return beta


def get_alpha(len_link2_3, len_link3_5, wx, wz, q1, joint_2_x_offset=0, joint_2_z_offset=0):
    x_dist_from_joint2 = (wx / cos(q1)) - joint_2_x_offset
    z_dist_from_joint2 = wz - joint_2_z_offset

    side_adj1 = len_link2_3
    side_adj2 = sqrt((x_dist_from_joint2 ** 2) + (z_dist_from_joint2 ** 2))
    side_opp = len_link3_5

    numerator = (side_opp ** 2) - (side_adj1 ** 2) - (side_adj2 ** 2)
    denominator = -2 * side_adj1 * side_adj2
    cos_alpha = numerator / denominator

    if cos_alpha > 1:
        mod_cos_alpha = cos_alpha - 1
        alpha_up = atan2(+sqrt(1 - (mod_cos_alpha ** 2)), mod_cos_alpha) - (pi / 2)
        alpha_down = atan2(-sqrt(1 - (mod_cos_alpha ** 2)), mod_cos_alpha) - (pi / 2)
    elif cos_alpha < 0:
        alpha_up = atan2(+sqrt(1 - (cos_alpha ** 2)), cos_alpha)
        alpha_down = atan2(-sqrt(1 - (cos_alpha ** 2)), cos_alpha)
    else:
        alpha_up = atan2(+sqrt(1 - (cos_alpha ** 2)), cos_alpha)
        alpha_down = atan2(-sqrt(1 - (cos_alpha ** 2)), cos_alpha)

    return alpha_up


def get_theta_2(alpha, beta):
    theta_2 = beta + alpha

    return theta_2


def rotate_x(rads):
    rotated = Matrix([
        [1, 0, 0, 0],
        [0, cos(rads), -sin(rads), 0],
        [0, sin(rads), cos(rads), 0],
        [0, 0, 0, 1]
    ])

    return rotated


def rotate_z(rads):
    rotated = Matrix([
        [cos(rads), -sin(rads), 0, 0],
        [sin(rads), cos(rads), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

    return rotated


def rotate_y(rads):
    rotated = Matrix([
        [cos(rads), 0, sin(rads), 0],
        [0, 1, 0, 0],
        [-sin(rads), 0, cos(rads), 0],
        [0, 0, 0, 1]
    ])

    return rotated


def get_wrist_coordinates(Rrpy, px, py, pz, dist_wrist_to_effector):
    lx = Rrpy[0, 0]
    ly = Rrpy[1, 0]
    lz = Rrpy[2, 0]

    # d6 is from the lesson. we use lx, ly, lz for an transform on the x axis instead of nx, ny, nz which is a transform
    # across the z axis
    wx = px - (dist_wrist_to_effector * lx)
    wy = py - (dist_wrist_to_effector * ly)
    wz = pz - (dist_wrist_to_effector * lz)

    return wx, wy, wz


def generate_rrpy_matrix(roll, pitch, yaw):
    Rot_X = rotate_x(roll)[0:3, 0:3]
    Rot_Y = rotate_y(pitch)[0:3, 0:3]
    Rot_Z = rotate_z(yaw)[0:3, 0:3]

    Rrpy = Rot_Z * Rot_Y * Rot_X

    return Rrpy

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

print("computing homogeneous transforms")

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

T3_4 = Matrix([[cos(q4), -sin(q4), 0, a3],
               [sin(q4) * cos(alpha3), cos(q4) * cos(alpha3), -sin(alpha3), -sin(alpha3) * d4],
               [sin(q4) * sin(alpha3), cos(q4) * sin(alpha3), cos(alpha3), cos(alpha3) * d4],
               [0, 0, 0, 1]])

T3_4 = T3_4.subs(s)

print("halfway through")

T4_5 = Matrix([[cos(q5), -sin(q5), 0, a4],
               [sin(q5) * cos(alpha4), cos(q5) * cos(alpha4), -sin(alpha4), -sin(alpha4) * d5],
               [sin(q5) * sin(alpha4), cos(q5) * sin(alpha4), cos(alpha4), cos(alpha4) * d5],
               [0, 0, 0, 1]])

T4_5 = T4_5.subs(s)

T5_6 = Matrix([[cos(q6), -sin(q6), 0, a5],
               [sin(q6) * cos(alpha5), cos(q6) * cos(alpha5), -sin(alpha5), -sin(alpha5) * d6],
               [sin(q6) * sin(alpha5), cos(q6) * sin(alpha5), cos(alpha5), cos(alpha5) * d6],
               [0, 0, 0, 1]])

T5_6 = T5_6.subs(s)

T6_G = Matrix([[cos(q7), -sin(q7), 0, a6],
               [sin(q7) * cos(alpha6), cos(q7) * cos(alpha6), -sin(alpha6), -sin(alpha6) * d7],
               [sin(q7) * sin(alpha6), cos(q7) * sin(alpha6), cos(alpha6), cos(alpha6) * d7],
               [0, 0, 0, 1]])

T6_G = T6_G.subs(s)

print("Done with individual matrices, proceeding to multiply")

T0_2 = T0_1 * T1_2
T0_3 = T0_2 * T2_3
# compose rotation matrix
# procedure taken from http://nghiaho.com/?page_id=846


# insert test cases and given data here the first objective is to produce Rrpy that is consistent with
# the one in Gazebo we will first obtain the givens in gazebo, roll, pitch, yaw and position coords we will
# then proceed to construct the Rrpy
roll = given_values[SHELF]["roll"]
pitch = given_values[SHELF]["pitch"]
yaw = given_values[SHELF]["yaw"]

px = given_values[SHELF]["px"]
py = given_values[SHELF]["py"]
pz = given_values[SHELF]["pz"]

Rrpy = generate_rrpy_matrix(roll, pitch, yaw)

dist_wrist_to_effector = s[d7]
d6_val = s[d6]

wx, wy, wz = get_wrist_coordinates(Rrpy, px, py, pz, dist_wrist_to_effector)

# compute theta for joint 1
theta1 = atan2(wy, wx)
theta1 = theta1.evalf()

# compute the side adjacents of theta 3
distance_joint_2to3 = s[a2]
distance_joint_3to5 = sqrt((s[a3] ** 2) + (s[d4] ** 2))

# get the offsets to the origin of joint 2
x_offset_to_joint2 = s[a1]
z_offset_to_joint2 = s[d1]

# compute for the elbow up value of theta 3
unadjusted_theta_3 = get_theta_3(distance_joint_2to3, distance_joint_3to5, wx, wz, theta1, x_offset_to_joint2,
                                 z_offset_to_joint2)

# choose the first of the two theta_3 results, which is the elbow up result
theta3 = unadjusted_theta_3.evalf() - RADS_AT_REST_JOINT3

# compute the parts used for theta 2
alpha = get_alpha(distance_joint_2to3, distance_joint_3to5, wx, wz, theta1, joint_2_x_offset=s[a1],
                  joint_2_z_offset=s[d1])
beta = get_beta(wx, wz, theta1, s[a1], s[d1])

# compute for theta2 without the offset
unadjusted_theta_2 = get_theta_2(alpha, beta).evalf()

# compute theta 2 value
theta2 = RADS_AT_REST_JOINT2 - unadjusted_theta_2
theta2 = theta2.evalf()

expected_theta1 = correct_values[SHELF]["joint_1"]
expected_theta2 = correct_values[SHELF]["joint_2"]
expected_theta3 = correct_values[SHELF]["joint_3"]

print("expected theta1 ", expected_theta1)
print("theta1          ", theta1)
print("expected theta2 ", expected_theta2)
print("theta2          ", theta2)
print("expected theta3 ", expected_theta3)
print("theta3          ", theta3)

R0_3 = T0_3[0:3, 0:3]
R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
# convert sympy matrix to numpy matrix to avoid errors
# how to convert lifted from: https://stackoverflow.com/a/37491889
R0_3 = np.array(R0_3).astype(np.float64)

Rrpy = (generate_rrpy_matrix(roll, pitch, yaw) * rotate_y(pi / 2)[0:3, 0:3] * rotate_z(pi)[0:3, 0:3])
R3_6 = (np.linalg.inv(R0_3)) * Rrpy

theta4 = atan2(R3_6[2, 2], -R3_6[0, 2])
theta4 = theta4.evalf()
theta5 = acos(R3_6[1, 2])
theta5 = theta5.evalf()
theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])
theta6 = theta6.evalf()

# if the value of theta4 is greater than or less than pi, then reverse the rotation of theta4
# by subtracting it from pi, and reversing the sign. Also reverse the rotations of theta5 and 6
# may be unecessary however
# if theta4 < -pi / 2:
#     theta4 = -(-pi - theta4)
#     theta5 = -theta5
#     theta6 = -(pi - theta6)
# elif theta4 > pi / 2:
#     theta4 = -(pi - theta4)
#     theta5 = -theta5
#     theta6 = -(-pi - theta6)

expected_theta4 = correct_values[SHELF]["joint_4"]
expected_theta5 = correct_values[SHELF]["joint_5"]
expected_theta6 = correct_values[SHELF]["joint_6"]

print("expected theta4 ", expected_theta4)
print("theta4          ", theta4)
print("expected theta5 ", expected_theta5)
print("theta5          ", theta5)
print("expected theta6 ", expected_theta6)
print("theta6          ", theta6)