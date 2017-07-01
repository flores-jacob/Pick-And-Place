import numpy as np
from numpy import array
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2
from sympy.matrices import Matrix

# These are the rads of joint 3 in resting position
RADS_AT_REST_JOINT3 = 1.60509488746813
# These are the rads of joint 2 in resting position
RADS_AT_REST_JOINT2 = 0.673308313432651


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


def get_alpha(len_link2_3, len_link3_5, theta_3):
    # TODO fix ugly hack of manually assigning alpha with a negative value
    alpha_side_adj = cos(theta_3) * len_link3_5
    alpha_side_opp = sin(theta_3) * len_link3_5

    alpha = -atan2(alpha_side_opp, len_link2_3 + alpha_side_adj)
    return alpha


def get_theta_2(alpha, beta):
    theta_2 = beta - alpha

    return theta_2


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
    [            0,               0,        1,      0],
    [            0,               0,        0,      1]
])

# Rotate along the y axis by -90 degrees
R_y = Matrix([
    [   cos(-pi/2),      0,      sin(-pi/2),      0],
    [               0,      1,                  0,      0],
    [  -sin(-pi/2),      0,      cos(-pi/2),      0],
    [               0,      0,                  0,      1]
])

R_corr = simplify(R_z * R_y)
# R_corr = R_z * R_y

# q1_val = 0.92
# q2_val = -0.51
# q3_val = 1.01
# q1_val = .4
# q2_val = .4
# q3_val = -2.95
# q4_val = -1.42
q1_val = -.4
q2_val = .5
q3_val = -.12
q4_val = 0
q5_val = 0.0
q6_val = 0.0

print("T0_1 = ", T0_1.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_2 = ", T0_2.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_3 = ", T0_3.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_4 = ", T0_4.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_5 = ", T0_5.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_6 = ", T0_6.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("T0_G = ", T0_G.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

T_total = simplify(T0_G * R_corr)
# T_total = T0_G * R_corr

# print("T_total = ", T_total.evalf(subs={q1: pi/4, q2: pi/4, q3:pi/4, q4:pi/4, q5:pi/4, q6:pi/4}))
print("T_total = ", T_total.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))

# T_total_evaluated = T_total.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val})

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
alpha = get_alpha(distance_joint_2to3, distance_joint_3to5, adjusted_theta_3)
beta = get_beta(wx, wz, theta_1, s[a1], s[d1])

# compute theta 2 value
theta_2 = RADS_AT_REST_JOINT2 - get_theta_2(alpha, beta).evalf()

print("theta_1 ", theta_1)
print("theta_2 ", theta_2)
print("theta_3 ", adjusted_theta_3)

print("theta_1 ", theta_1.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta_2 ", theta_2.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
print("theta_3 ", adjusted_theta_3.evalf(subs={q1: q1_val, q2: q2_val, q3: q3_val, q4: q4_val, q5: q5_val, q6: q6_val}))
