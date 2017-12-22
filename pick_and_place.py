import random

import numpy as np
from numpy import array
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2, acos
from sympy.matrices import Matrix, Transpose

from sample_data import correct_values

# These are the rads of joint 3 in resting position
# RADS_AT_REST_JOINT3 = 1.60509488746813
# RADS_AT_REST_JOINT3 = 1.60926888235184
RADS_AT_REST_JOINT3 = 1.60678078687695

# These are the rads of joint 2 in resting position
# RADS_AT_REST_JOINT2 = 0.673308313432651
# RADS_AT_REST_JOINT2 = 0.674740942223551
# RADS_AT_REST_JOINT2 = 0.674740942223554
# RADS_AT_REST_JOINT2 = 0.671755266989842
RADS_AT_REST_JOINT2 = pi/2 # 1.57079632679490

RADS_AT_REST_JOINT4 = pi/2 # 1.57079632679490

RADS_AT_REST_JOINT6 = pi/2 # 1.57079632679490


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

    return theta_3_elbow_up, theta_3_elbow_down


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

    # alpha = acos(cos_alpha)


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


def get_rpy(rotation_matrix):
    # solution adapted from this page http://nghiaho.com/?page_id=846
    r11 = rotation_matrix[0, 0]
    r12 = rotation_matrix[0, 1]
    r13 = rotation_matrix[0, 2]

    r21 = rotation_matrix[1, 0]
    r22 = rotation_matrix[1, 1]
    r23 = rotation_matrix[1, 2]

    r31 = rotation_matrix[2, 0]
    r32 = rotation_matrix[2, 1]
    r33 = rotation_matrix[2, 2]

    theta_x = atan2(r32, r33)
    theta_y = atan2(-r31, sqrt((r32 ** 2) + (r33 ** 2)))
    theta_z = atan2(r21, r11)


    # theta5 = atan2(-R3_6[2, 0], sqrt((R3_6[0, 0] * R3_6[0, 0]) + (R3_6[1, 0] * R3_6[1, 0])))
    # theta4 = atan2(R3_6[2, 1], R3_6[2, 2])
    # theta6 = atan2(R3_6[1, 0], R3_6[0, 0])

    return (theta_x, theta_y, theta_z)


def get_roll(rotation_matrix):
    r32 = rotation_matrix[2, 1]
    r33 = rotation_matrix[2, 2]

    roll = atan2(r32, r33)

    return roll


def get_pitch(rotation_matrix):
    r31 = rotation_matrix[2, 0]
    r32 = rotation_matrix[2, 1]
    r33 = rotation_matrix[2, 2]

    pitch = atan2(-r31, sqrt((r32 ** 2) + (r33 ** 2)))

    return pitch


def get_yaw(rotation_matrix):
    # solution adapted from this page http://nghiaho.com/?page_id=846
    r11 = rotation_matrix[0, 0]
    r21 = rotation_matrix[1, 0]

    yaw = atan2(r21, r11)

    return yaw

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

# create symbols for variables
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

# DH Parameters

s = {alpha0:        0,  a0:       0,  d1: 0.75,
     alpha1:    -pi/2,  a1:    0.35,  d2:    0,    q2: q2 - pi/2,
     alpha2:        0,  a2:    1.25,  d3:    0,
     alpha3:    -pi/2,  a3:  -0.054,  d4: 1.50,
     alpha4:     pi/2,  a4:       0,  d5:    0,
     alpha5:    -pi/2,  a5:       0,  d6:    0,
     alpha6:        0,  a6:       0,  d7:0.303,     q7: 0
     }

print("computing homogeneous transforms")

# Homogeneous transforms
T0_1 = Matrix([ [               cos(q1),              -sin(q1),            0,                a0],
                [ sin(q1) * cos(alpha0), cos(q1) * cos(alpha0), -sin(alpha0), -sin(alpha0) * d1],
                [ sin(q1) * sin(alpha0), cos(q1) * sin(alpha0),  cos(alpha0),  cos(alpha0) * d1],
                [                     0,                     0,            0,                 1]])

T0_1 = T0_1.subs(s)

T1_2 = Matrix([ [               cos(q2),              -sin(q2),            0,                a1],
                [ sin(q2) * cos(alpha1), cos(q2) * cos(alpha1), -sin(alpha1), -sin(alpha1) * d2],
                [ sin(q2) * sin(alpha1), cos(q2) * sin(alpha1),  cos(alpha1),  cos(alpha1) * d2],
                [                     0,                     0,            0,                 1]])

T1_2 = T1_2.subs(s)

T2_3 = Matrix([ [               cos(q3),              -sin(q3),            0,                a2],
                [ sin(q3) * cos(alpha2), cos(q3) * cos(alpha2), -sin(alpha2), -sin(alpha2) * d3],
                [ sin(q3) * sin(alpha2), cos(q3) * sin(alpha2),  cos(alpha2),  cos(alpha2) * d3],
                [                     0,                     0,            0,                 1]])

T2_3 = T2_3.subs(s)

T3_4 = Matrix([ [               cos(q4),              -sin(q4),            0,                a3],
                [ sin(q4) * cos(alpha3), cos(q4) * cos(alpha3), -sin(alpha3), -sin(alpha3) * d4],
                [ sin(q4) * sin(alpha3), cos(q4) * sin(alpha3),  cos(alpha3),  cos(alpha3) * d4],
                [                     0,                     0,            0,                 1]])

T3_4 = T3_4.subs(s)

T4_5 = Matrix([ [               cos(q5),              -sin(q5),            0,                a4],
                [ sin(q5) * cos(alpha4), cos(q5) * cos(alpha4), -sin(alpha4), -sin(alpha4) * d5],
                [ sin(q5) * sin(alpha4), cos(q5) * sin(alpha4),  cos(alpha4),  cos(alpha4) * d5],
                [                     0,                     0,            0,                 1]])

T4_5 = T4_5.subs(s)

T5_6 = Matrix([ [               cos(q6),              -sin(q6),            0,                a5],
                [ sin(q6) * cos(alpha5), cos(q6) * cos(alpha5), -sin(alpha5), -sin(alpha5) * d6],
                [ sin(q6) * sin(alpha5), cos(q6) * sin(alpha5),  cos(alpha5),  cos(alpha5) * d6],
                [                     0,                     0,            0,                 1]])

T5_6 = T5_6.subs(s)

T6_G = Matrix([ [               cos(q7),              -sin(q7),            0,                a6],
                [ sin(q7) * cos(alpha6), cos(q7) * cos(alpha6), -sin(alpha6), -sin(alpha6) * d7],
                [ sin(q7) * sin(alpha6), cos(q7) * sin(alpha6),  cos(alpha6),  cos(alpha6) * d7],
                [                     0,                     0,            0,                 1]])

T6_G = T6_G.subs(s)

print("Done with individual matrices, proceeding to multiply")

T0_2 = simplify(T0_1 * T1_2)
T0_3 = simplify(T0_2 * T2_3)
T0_4 = simplify(T0_3 * T3_4)
T0_5 = simplify(T0_4 * T4_5)
T0_6 = simplify(T0_5 * T5_6)
T0_G = simplify(T0_6 * T6_G)

# T0_2 = T0_1 * T1_2
# T0_3 = T0_2 * T2_3
# T0_4 = T0_3 * T3_4
# T0_5 = T0_4 * T4_5
# T0_6 = T0_5 * T5_6
# T0_G = T0_6 * T6_G

# Gripper link  orientation correction so that URDF values are in accordance with DH Convention

print("computing orientations")

# Rotate z by 180 degrees
R_z = Matrix([
    [   cos(pi),     -sin(pi),        0,      0],
    [   sin(pi),      cos(pi),        0,      0],
    [         0,            0,        1,      0],
    [         0,            0,        0,      1]
])

# Rotate along the y axis by -90 degrees
R_y = Matrix([
    [   cos(-pi/2),      0,      sin(-pi/2),      0],
    [            0,      1,               0,      0],
    [  -sin(-pi/2),      0,      cos(-pi/2),      0],
    [            0,      0,               0,      1]
])

# R_corr_4_6_and_gripper = simplify(R_z * R_y)
R_corr_4_6_and_gripper = R_z * R_y

# Rotate along the z axis by 90 degrees
R_z90 = Matrix([
    [   cos(pi/2),     -sin(pi/2),        0,      0],
    [   sin(pi/2),      cos(pi/2),        0,      0],
    [         0,            0,        1,      0],
    [         0,            0,        0,      1]
])

R_x = Matrix([
    [   1,              0,              0,      0],
    [   0,      cos(pi/2),     -sin(pi/2),      0],
    [   0,      sin(pi/2),      cos(pi/2),      0],
    [   0,              0,              0,      1]
])

# R_corr_3_and_5 = simplify(R_z90 * R_x)
R_corr_3_and_5 = R_z90 * R_x

# q1_val = 0.92
# q2_val = -0.51
# q3_val = 1.01
# q4_val = -1.42
# q1_val = -.4
# q2_val = .5
# q3_val = -.12
# q1_val = .8
# q2_val = -.1
# q3_val = .25
q1_val = 0
q2_val = 0
q3_val = 0
# q4_val = 0.0
q4_val = 0.0
q5_val = 0
q6_val = 0.0
# q1_val = .4
# q2_val = -.8
# q3_val = -.11
# q4_val = -1.56
# q5_val = 1.9
# q6_val = 0
# q1_val = 0.4431106733822272
# q2_val = 0.24109151255939043
# q3_val = -0.019039282352903975
# q4_val = 1.1356308927170824
# q5_val = -0.493364849487004
# q6_val = -1.1356308927170824
q1_val = -0.44293154231873455
q2_val =  0.5386775324322572
q3_val =  0.16542712045273245
q4_val =  -0.6314924490169878
# q5_val =  -0.8115364743375029
q6_val =  0.46657506948980654

# q1_val = correct_values["shelf_3"]["joint_1"]
# q2_val =  correct_values["shelf_3"]["joint_2"]
# q3_val = correct_values["shelf_3"]["joint_3"]
# q4_val =  correct_values["shelf_3"]["joint_4"]
# q5_val =  correct_values["shelf_3"]["joint_5"]
# q6_val =  correct_values["shelf_3"]["joint_6"]

print("T0_1 = ", T0_1.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_2 = ", T0_2.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_3 = ", T0_3.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_4 = ", T0_4.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_5 = ", T0_5.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_6 = ", T0_6.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_G = ", T0_G.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

# T_total = simplify(T0_G * R_corr_4_6_and_gripper)
T_total = T0_G * R_corr_4_6_and_gripper

print("T_total = ", T_total.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

Rrpy = T_total.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val})

# print(Rrpy)

# px = T_total[0, 3]
# py = T_total[1, 3]
# pz = T_total[2, 3]
# 
# nx = T_total[0, 2]
# ny = T_total[1, 2]
# nz = T_total[2, 2]
# 
# lx = T_total[0, 0]
# ly = T_total[1, 0]
# lz = T_total[2, 0]

px = Rrpy[0, 3]
py = Rrpy[1, 3]
pz = Rrpy[2, 3]

nx = Rrpy[0, 2]
ny = Rrpy[1, 2]
nz = Rrpy[2, 2]

lx = Rrpy[0, 0]
ly = Rrpy[1, 0]
lz = Rrpy[2, 0]

dist_wrist_to_effector = s[d7]
d6_val = s[d6]

# d6 is from the lesson. we use lx, ly, lz for an transform on the x axis instead of nx, ny, nz which is a transform
# across the z axis
wx = px - ((d6_val + dist_wrist_to_effector) * lx)
wy = py - ((d6_val + dist_wrist_to_effector) * ly)
wz = pz - ((d6_val + dist_wrist_to_effector) * lz)


# wx = px - (.193 * lx)
# wy = py - (.193 * ly)
# wz = pz - (.193 * lz)

# get wirst coordinate values
# wx = wx.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val})
# wy = wy.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val})
# wz = wz.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val})

# wx = wx.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val})
# wy = wy.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val})
# wz = wz.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val})


# wx = 1.1868
#
# wy = 1.2207
#
# wz = 1.7115

print("wx ", wx)
print("wy ", wy)
print("wz ", wz)

##############


# compute theta for joint 1
theta_1 = atan2(wy, wx)
theta1 = theta_1.evalf()

# compute the side adjacents of theta 3
distance_joint_2to3 = s[a2]
distance_joint_3to5 = sqrt((s[a3] ** 2) + (s[d4] ** 2))

# get the offsets to the origin of joint 2
x_offset_to_joint2 = s[a1]
z_offset_to_joint2 = s[d1]

# compute for the 2 possible values of theta 3
unadjusted_theta_3 = get_theta_3(distance_joint_2to3, distance_joint_3to5, wx, wz, theta_1, x_offset_to_joint2,
                                 z_offset_to_joint2)

# choose the first of the two theta_3 results, which is the elbow up result
theta_3 = unadjusted_theta_3[0].evalf() - RADS_AT_REST_JOINT3
theta3 = theta_3

# compute the parts used for theta 2
alpha = get_alpha(distance_joint_2to3, distance_joint_3to5, wx, wz, theta_1, joint_2_x_offset=s[a1],
                  joint_2_z_offset=s[d1])
beta = get_beta(wx, wz, theta_1, s[a1], s[d1])

unadjusted_theta_2 = get_theta_2(alpha, beta).evalf()

# compute theta 2 value
theta_2 = RADS_AT_REST_JOINT2 - unadjusted_theta_2
theta2 = theta_2.evalf()

########


# print("theta_1 ", theta_1)
# print("theta_2 ", theta_2)
# print("theta_3 ", adjusted_theta_3)

# print("theta_1 ", theta_1.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta_2 ", theta_2.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta_3 ", theta_3.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
#
# print("unadjusted theta_2 ",
#       unadjusted_theta_2.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("unadjusted theta_3 ",
#       theta_3.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

# R0_1 = T0_1[0:3, 0:3].subs(s)
# R1_2 = T1_2[0:3, 0:3].subs(s)
# R2_3 = T2_3[0:3, 0:3].subs(s)
# R3_4 = T3_4[0:3, 0:3].subs(s)
# R4_5 = T4_5[0:3, 0:3].subs(s)
# R5_6 = T5_6[0:3, 0:3].subs(s)
# R6_G = T6_G[0:3, 0:3].subs(s)


R0_1 = T0_1[0:3, 0:3]
R1_2 = T1_2[0:3, 0:3]
R2_3 = T2_3[0:3, 0:3]
R3_4 = T3_4[0:3, 0:3]
R4_5 = T4_5[0:3, 0:3]
R5_6 = T5_6[0:3, 0:3]
R6_G = T6_G[0:3, 0:3]

# # Rotation from base to the 3rd joint
# R0_3 = R0_1 * R1_2 * R2_3 * R3_4 * R_corr[0:3, 0:3]
#
# # Rotation from the base to the gripper
# R0_G = R0_1 * R1_2 * R2_3 * R3_4 * R4_5 * R5_6 * R6_G
# R_Total = R0_G * R_corr[0:3, 0:3]

# Rotation from base to the 3rd joint
# R0_3 = R0_1 * R1_2 * R2_3 * R_corr[0:3, 0:3]
# # R0_3 = simplify(T0_3[0:3, 0:3] * R_corr_3_and_5[0:3, 0:3])
# R0_3 = T0_3[0:3, 0:3]
# # R0_3_corr = R0_3 * R_corr_3_and_5[0:3, 0:3]
# # R0_3_corr = (R0_3 * rotate_x(-pi/2)[0:3, 0:3]) * rotate_z(-pi/2)[0:3, 0:3]
# R0_3_corr = R0_3 * R_corr_3_and_5[0:3, 0:3]
# R0_3_corr = R0_3 * rotate_z(pi/2)[0:3, 0:3]
# R0_3_corr = R0_3_corr * rotate_x(pi/2)[0:3, 0:3]
# R_check = R0_3 * R_corr_3_and_5[0:3, 0:3]
# R_check = R_check.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

# R0_3_corr = R0_3_corr.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

# check = R0_3_corr == R_check
#
# print("check ", check)
#
# print R0_3_corr
# print R_check


# get the matrix for base link to joint 3, using computed theta1, 2, and 3 values
R0_3 = T0_3[0:3, 0:3]
R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
# convert sympy matrix to numpy matrix to avoid errors
# how to convert lifted from: https://stackoverflow.com/a/37491889
# R0_3 = np.array(R0_3).astype(np.float64)

# correct the orientation of the Rrpy matrix, make its axes align with the DH axes of the gripper
Rrpy = Rrpy[0:3, 0:3]  * rotate_y(pi / 2)[0:3, 0:3] * rotate_z(pi)[0:3, 0:3]
# get the matrix values for the wrist
# R3_6 = (np.linalg.inv(R0_3)) * Rrpy
R3_6 = Transpose(R0_3) * Rrpy


# compute values for thetas 4, 5, and 6
theta4 = atan2(R3_6[2, 2], -R3_6[0, 2])
theta4 = theta4.evalf()
theta5 = acos(R3_6[1, 2])
theta5 = theta5.evalf()
theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])
theta6 = theta6.evalf()

if theta4 < -pi / 2:
    theta4 = -(-pi - theta4)
    theta5 = -theta5
    theta6 = -(pi - theta6)
elif theta4 > pi / 2:
    theta4 = -(pi - theta4)
    theta5 = -theta5
    theta6 = -(-pi - theta6)

# Rotation from base to the 4th joint
# R0_4 = R0_1 * R1_2 * R2_3 * R3_4 * R_corr[0:3, 0:3]
# R0_4 = simplify(T0_4[0:3, 0:3] * R_corr_4_6_and_gripper[0:3, 0:3])
# R0_4 = T0_4[0:3, 0:3] * R_corr_4_6_and_gripper[0:3, 0:3]

# Rotation from base to the 5th joint
# R0_5 = R0_1 * R1_2 * R2_3 * R3_4 * R4_5 * R_corr[0:3, 0:3]
# R0_5 = simplify(T0_5[0:3, 0:3] * R_corr_3_and_5[0:3, 0:3])
# R0_5 = T0_5[0:3, 0:3] * R_corr_3_and_5[0:3, 0:3]

# Rotation from base to the 6th joint
# R0_6 = R0_1 * R1_2 * R2_3 * R3_4 * R4_5 * R5_6 * R_corr[0:3, 0:3]
# R0_6 = simplify(T0_6[0:3, 0:3] * R_corr_4_6_and_gripper[0:3, 0:3])
# R0_6 = T0_6[0:3, 0:3]
# R0_6_corr = R0_6 * R_corr_4_6_and_gripper[0:3, 0:3]
# R0_6_corr = R0_6_corr.evalf(subs={q1:theta1, q2:theta2, q3:theta3})
# Rotation from base to the gripper joint
# R0_G = R0_1 * R1_2 * R2_3 * R3_4 * R4_5 * R5_6 * R6_G * R_corr[0:3, 0:3]
# R0_G = simplify(T0_G[0:3, 0:3] * R_corr[0:3, 0:3])

# # Rotation from the base to the gripper
# R0_G = R0_1 * R1_2 * R2_3 * R3_4 * R4_5 * R5_6 * R6_G
# R_Total = R0_G * R_corr[0:3, 0:3]

# R_Total = R0_G
# R_Total = R0_6


# Rotation from R3 to the End Effector
# R_product = R0_3.inv() * R_Total
# R_product = simplify(R0_3.inv() * R_Total)

# Rotation from the third joint to the 4th joint
# R_3to4 = R0_3.inv() * R0_4
# R_3to4 = simplify(R0_3.inv() * R0_4)

# Rotation from the 4th joint to the 5th joint
# R_4to5 = R0_4.inv() * R0_5
# R_4to5 = simplify(R0_4.inv() * R0_5)

# Rotation from the 5th joint to the 6th joint
# R_5to6 = R0_5.inv() * R0_6
# R_5to6 = simplify(R0_5.inv() * R0_6)

# R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
#
# R3_6 = R0_3.inv() * T_total[0:3, 0:3]

# R0_3 = R0_3_corr.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
#
# Rrpy_corr = Rrpy[0:3, 0:3] * R_corr_4_6_and_gripper[0:3, 0:3].evalf()
# # Rrpy_corr = T_total[0:3, 0:3]
# R3_6 = R0_3.inv() * Rrpy_corr
# # R3_6 = (R3_6 * rotate_x(-pi/2)[0:3, 0:3]) * rotate_z(-pi/2)[0:3, 0:3]
#
# # TODO create a for loop that will test each and every possible permutation
# # R3_6 = (R3_6 * rotate_x(-pi/2)[0:3, 0:3]) * rotate_z(-pi/2)[0:3, 0:3]
#
#
# # # TODO the value for theta 4 could simply be the negative value of roll
# # # this is to counteract the roll of the arm
# # unadjusted_theta4 = get_roll(R3_6).evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4:q4_val, q5:q5_val, q6:q6_val})
# # # theta4 = RADS_AT_REST_JOINT4 - unadjusted_theta4
# # theta4 = unadjusted_theta4
# # # TODO the value for theta 5 could simply be the value for pitch
# # theta5 = get_pitch(R3_6).evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4:q4_val, q5:q5_val, q6:q6_val})
# # unadjusted_theta6 = get_yaw(R3_6).evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4:q4_val, q5:q5_val, q6:q6_val})
# # # theta6 = RADS_AT_REST_JOINT6 - unadjusted_theta6
# # # TODO theta6 could simply be the opposite of theta4
# # theta6 = unadjusted_theta6
#
#
# # Experimental code that should work
# R3_6_corr = ((R3_6 * rotate_y(pi / 2)[0:3, 0:3]) * rotate_z(-pi / 2)[0:3, 0:3]) #* rotate_z(pi)[0:3, 0:3]
#
# # R3_6_corr = (R3_6 * rotate_x(-pi/2)[0:3, 0:3]) * rotate_z(-pi/2)[0:3, 0:3]
#
# # Why -pitch for theta4 and roll for theta5? When theta4 is roll, and theta5 is pitch? This is because these are the
# # are the formulas that return the desired values base on the different permutations we have tried
# theta4 = -get_pitch(R3_6_corr).evalf(subs={q1: theta1, q2: theta2, q3: theta3})
# theta5 = get_roll(R3_6_corr).evalf(subs={q1: theta1, q2: theta2, q3: theta3})
# # theta6 is simply the revers of theta4
# theta6 = get_yaw(R3_6_corr).evalf(subs={q1: theta1, q2: theta2, q3: theta3})


# alpha = get_roll(R3_6_corr).evalf(subs={q1: theta1, q2: theta2, q3: theta3})
# beta = get_pitch(R3_6_corr).evalf(subs={q1: theta1, q2: theta2, q3: theta3})
# gamma = get_yaw(R3_6_corr).evalf(subs={q1: theta1, q2: theta2, q3: theta3})
#
# alpha = atan2(R3_6[1, 0], R3_6[0, 0])  # rotation about Z-axis
# beta = atan2(-R3_6[2, 0],
#              sqrt((R3_6[0, 0] * R3_6[0, 0]) + (R3_6[1, 0] * R3_6[1, 0])))  # rotation about Y-axis
# gamma = atan2(R3_6[2, 1], R3_6[2, 2])  # rotation about X-axis

# theta4 = gamma.evalf(subs=({q1:theta1, q2:theta2, q3:theta3}))
# theta5 = beta.evalf(subs=({q1:theta1, q2:theta2, q3:theta3}))
# theta6 = alpha.evalf(subs=({q1:theta1, q2:theta2, q3:theta3}))

# print("q4_val ", q4_val)
# print("q5_val ", q5_val)
# print("q6_val ", q6_val)
print("expected theta1 ", q1_val)
print("theta1          ", theta1)
print("expected theta2 ", q2_val)
print("theta2          ", theta2)
print("expected theta3 ", q3_val)
print("theta3          ", theta3)
print("expected theta4 ", q4_val)
print("theta4          ", theta4.evalf())
print("expected theta5 ", q5_val)
print("theta5          ", theta5.evalf())
print("expected theta6 ", q6_val)
print("theta6          ", theta6.evalf())


# print("alpha ", alpha)
# print("beta ", beta)
# print("gamma ", gamma)


# theta4_error_count = 0
# theta5_error_count = 0
#
# for i in range(1000):
#     q4_val = round(random.uniform(-3.14159265359/2, 3.14159265359/2), 2)
#     q5_val = round(random.uniform(-3.14159265359/2, 3.14159265359/2), 2)
#
#     # The correction value here was taken based on testing using different available permutations of the axes
#     R3_6_corr = (R3_6 * rotate_y(pi / 2)[0:3, 0:3]) * rotate_z(-pi / 2)[0:3, 0:3]
#     R3_6 = R3_6_corr
#
#
#     theta4 = -get_pitch(R3_6).evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4:q4_val, q5:q5_val, q6:q6_val})
#     theta5 = get_roll(R3_6).evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4: q4_val, q5: q5_val, q6: q6_val})
#
#     if round(theta4, 5) != q4_val:
#         print("q4_val not matching. expecting", q4_val, "got ", theta4, "instead ")
#         theta4_error_count +=1
#     if round(theta5, 5) != q5_val:
#         print("q5_val not matching. expecting", q5_val, "got ", theta5, "instead ")
#         theta5_error_count +=1
#
#
# print ("theta4 errors ", theta4_error_count)
# print ("theta5 errors ", theta5_error_count)




# theta_x3EE, theta_y3EE, theta_z3EE = get_rpy(R_product)
# # theta_x34, theta_y34, theta_z34 = get_rpy(R_3to4)
# theta_4 = get_roll(R_3to4)
# print("halfway done")
# # theta_x45, theta_y45, theta_z45 = get_rpy(R_4to5)
# theta_5 = get_pitch(R_4to5)
# print("computing 5 to 6")
# # theta_x56, theta_y56, theta_z56 = get_rpy(R_5to6)
# theta_6 = get_roll(R_5to6)

# print("computing 6 to G")
# theta_x6G, theta_y6G, theta_z6G = get_rpy(R_6toG)


# print("theta x 3EE", theta_x3EE.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta y 3EE", theta_y3EE.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta z 3EE", theta_z3EE.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

# print("joint 4 ", theta_4.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("joint 5 ", theta_5.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("joint 6 ", theta_6.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

# print("joint 4 ", theta_4.evalf(subs={q1: theta_1, q2: theta_2, q3: adjusted_theta_3, q4: theta_4, q5: theta_5, q6: theta_6}))
# print("joint 5 ", theta_5.evalf(subs={q1: theta_1, q2: theta_2, q3: adjusted_theta_3, q4: theta_4, q5: theta_5, q6: theta_6}))
# print("joint 6 ", theta_6.evalf(subs={q1: theta_1, q2: theta_2, q3: adjusted_theta_3, q4: theta_4, q5: theta_5, q6: theta_6}))

# print("joint 4 ", theta_4.evalf(subs={q1: theta_1, q2: theta_2, q3: adjusted_theta_3}))
# print("joint 5 ", theta_5.evalf(subs={q1: theta_1, q2: theta_2, q3: adjusted_theta_3}))
# print("joint 6 ", theta_6.evalf(subs={q1: theta_1, q2: theta_2, q3: adjusted_theta_3}))

# print("theta x 34", theta_x34.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta y 34", theta_y34.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta z 34", theta_z34.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
#
# print("theta x 45", theta_x45.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta y 45", theta_y45.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta z 45", theta_z45.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
#
# print("theta x 56", theta_x56.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta y 56", theta_y56.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta z 56", theta_z56.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

# print("theta x 6G", theta_x6G.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta y 6G", theta_y6G.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta z 6G", theta_z6G.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

