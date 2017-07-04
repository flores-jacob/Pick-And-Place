import numpy as np
from numpy import array
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2, acos
from sympy.matrices import Matrix

# These are the rads of joint 3 in resting position
# RADS_AT_REST_JOINT3 = 1.60509488746813
# RADS_AT_REST_JOINT3 = 1.60926888235184
RADS_AT_REST_JOINT3 = 1.60678078687695

# These are the rads of joint 2 in resting position
# RADS_AT_REST_JOINT2 = 0.673308313432651
# RADS_AT_REST_JOINT2 = 0.674740942223551
# RADS_AT_REST_JOINT2 = 0.674740942223554
# RADS_AT_REST_JOINT2 = 0.671755266989842
RADS_AT_REST_JOINT2 = 1.57079632679490
# 1.57324022779343


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
    cos_alpha = numerator/denominator

    alpha = acos(cos_alpha)

    return alpha


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

    return (theta_x, theta_y, theta_z)

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

print("halfway through")

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

R_corr_4_6_and_gripper = simplify(R_z * R_y)

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

R_corr_3_and_5 = simplify(R_z90 * R_x)


# q1_val = 0.92
# q2_val = -0.51
# q3_val = 1.01
# q1_val = .4
# q2_val = .4
# q3_val = -.11
# q4_val = -1.42
# q1_val = -.4
# q2_val = .5
# q3_val = -.12
q1_val = .8
q2_val = -.1
q3_val = .25
# q1_val = 0.0
# q2_val = 0.0
# q3_val = 0.0
q4_val = 0.0
q5_val = 0.0
q6_val = 0.0
# q4_val = 1
# q5_val = -1
# q6_val = .3


print("T0_1 = ", T0_1.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_2 = ", T0_2.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_3 = ", T0_3.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_4 = ", T0_4.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_5 = ", T0_5.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_6 = ", T0_6.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_G = ", T0_G.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

T_total = simplify(T0_G * R_corr_4_6_and_gripper)
# T_total = T0_G * R_corr

print("T_total = ", T_total.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

px = T_total[0, 3]
py = T_total[1, 3]
pz = T_total[2, 3]

nx = T_total[0, 2]
ny = T_total[1, 2]
nz = T_total[2, 2]

lx = T_total[0, 0]
ly = T_total[1, 0]
lz = T_total[2, 0]

dist_wrist_to_effector = s[d7]
d6_val = s[d6]

# d6 is from the lesson. we use lx, ly, lz for an transform on the x axis instead of nx, ny, nz which is a transform
# across the z axis
wx = px - ((d6_val + dist_wrist_to_effector) * lx)
wy = py - ((d6_val + dist_wrist_to_effector) * ly)
wz = pz - ((d6_val + dist_wrist_to_effector) * lz)

# get wirst coordinate values
wx = wx.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val})
wy = wy.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val})
wz = wz.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val})

# wx = 1.1868
#
# wy = 1.2207
#
# wz = 1.7115

print("wx ", wx)
print("wy ", wy)
print("wz ", wz)

# compute theta for joint 1
theta_1 = atan2(wy, wx)

# compute the side adjacents of theta 3
distance_joint_2to3 = s[a2]
distance_joint_3to5 = sqrt((s[a3] ** 2) + (s[d4] ** 2))

# get the offsets to the origin of joint 2
x_offset_to_joint2 = s[a1]
z_offset_to_joint2 = s[d1]

# compute for the 2 possible values of theta 3
theta_3 = get_theta_3(distance_joint_2to3, distance_joint_3to5, wx, wz, theta_1, x_offset_to_joint2, z_offset_to_joint2)

# choose the first of the two theta_3 results, which is the elbow up result
adjusted_theta_3 = theta_3[0].evalf() - RADS_AT_REST_JOINT3

# compute the parts used for theta 2
alpha = get_alpha(distance_joint_2to3, distance_joint_3to5, wx, wz, theta_1, joint_2_x_offset=s[a1], joint_2_z_offset=s[d1])
beta = get_beta(wx, wz, theta_1, s[a1], s[d1])

unadjusted_theta_2 = get_theta_2(alpha, beta).evalf()

# compute theta 2 value
theta_2 = RADS_AT_REST_JOINT2 - unadjusted_theta_2

# print("theta_1 ", theta_1)
# print("theta_2 ", theta_2)
# print("theta_3 ", adjusted_theta_3)

print("theta_1 ", theta_1.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta_2 ", theta_2.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta_3 ", adjusted_theta_3.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

print("unadjusted theta_2 ", unadjusted_theta_2.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("unadjusted theta_3 ", theta_3[0].evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))


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
R0_3 = simplify(T0_3[0:3, 0:3] * R_corr_3_and_5[0:3, 0:3])
# R0_3 = T0_3[0:3, 0:3] * R_corr[0:3, 0:3]


# Rotation from base to the 4th joint
# R0_4 = R0_1 * R1_2 * R2_3 * R3_4 * R_corr[0:3, 0:3]
R0_4 = simplify(T0_4[0:3, 0:3] * R_corr_4_6_and_gripper[0:3, 0:3])
# R0_4 = T0_4[0:3, 0:3] * R_corr[0:3, 0:3]


# Rotation from base to the 5th joint
# R0_5 = R0_1 * R1_2 * R2_3 * R3_4 * R4_5 * R_corr[0:3, 0:3]
R0_5 = simplify(T0_5[0:3, 0:3] * R_corr_3_and_5[0:3, 0:3])
# R0_5 = T0_5[0:3, 0:3] * R_corr[0:3, 0:3]


# Rotation from base to the 6th joint
# R0_6 = R0_1 * R1_2 * R2_3 * R3_4 * R4_5 * R5_6 * R_corr[0:3, 0:3]
R0_6 = simplify(T0_6[0:3, 0:3] * R_corr_4_6_and_gripper[0:3, 0:3])
# R0_6 = T0_6[0:3, 0:3] * R_corr[0:3, 0:3]


# Rotation from base to the gripper joint
# R0_G = R0_1 * R1_2 * R2_3 * R3_4 * R4_5 * R5_6 * R6_G * R_corr[0:3, 0:3]
# R0_G = simplify(T0_G[0:3, 0:3] * R_corr[0:3, 0:3])

# # Rotation from the base to the gripper
# R0_G = R0_1 * R1_2 * R2_3 * R3_4 * R4_5 * R5_6 * R6_G
# R_Total = R0_G * R_corr[0:3, 0:3]

# R_Total = R0_G
R_Total = R0_6

print("applying rotations")

# Rotation from R3 to the End Effector
# R_product = R0_3.inv() * R_Total
R_product = simplify(R0_3.inv() * R_Total)

# Rotation from the third joint to the 4th joint
# R_3to4 = R0_3.inv() * R0_4
R_3to4 = simplify(R0_3.inv() * R0_4)

# Rotation from the 4th joint to the 5th joint
# R_4to5 = R0_4.inv() * R0_5
R_4to5 = simplify(R0_4.inv() * R0_5)


print("almost done")
# Rotation from the 5th joint to the 6th joint
# R_5to6 = R0_5.inv() * R0_6
# R_5to6 = simplify(R0_5.inv() * R0_6)


# print("computing 6 to G")
# Rotation from the 6th joint to the gripper joint
# R_6toG = R0_6.inv() * R0_G

print("done with multiplications")

# print("R total 1", R_product.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("R total 2", R_product.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: 0.44, q5: q5_val, q6: q6_val}))

# R3_EE = simplify(T0_3.inv() * T_total)
# R3_EE = T0_3.inv() * T_total
#
# print(R3_EE[0:3, 0:3].evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

print("computing rpy values")

theta_x3EE, theta_y3EE, theta_z3EE = get_rpy(R_product)
theta_x34, theta_y34, theta_z34 = get_rpy(R_3to4)
print("halfway done")
theta_x45, theta_y45, theta_z45 = get_rpy(R_4to5)
print("computing 5 to 6")
# theta_x56, theta_y56, theta_z56 = get_rpy(R_5to6)
# print("computing 6 to G")
# theta_x6G, theta_y6G, theta_z6G = get_rpy(R_6toG)


print("theta x 3EE", theta_x3EE.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta y 3EE", theta_y3EE.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta z 3EE", theta_z3EE.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))


print("theta x 34", theta_x34.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta y 34", theta_y34.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta z 34", theta_z34.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

print("theta x 45", theta_x45.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta y 45", theta_y45.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta z 45", theta_z45.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

# print("theta x 56", theta_x56.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta y 56", theta_y56.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta z 56", theta_z56.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

# print("theta x 6G", theta_x6G.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta y 6G", theta_y6G.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
# print("theta z 6G", theta_z6G.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))


