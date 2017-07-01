import numpy as np
from numpy import array
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2
from sympy.matrices import Matrix


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


# # Correct wrist orientation
# W_z = Matrix([
#     [   cos(pi/2),     -sin(pi/2),        0,      0],
#     [   sin(pi/2),      cos(pi/2),        0,      0],
#     [            0,             0,        1,      0],
#     [            0,             0,        0,      1]
# ])
#
# W_x = Matrix([
#     [   1,              0,              0,      0],
#     [   0,      cos(pi/2),     -sin(pi/2),      0],
#     [   0,      sin(pi/2),      cos(pi/2),      0],
#     [   0,              0,              0,      1]
# ])
#
# # Rotate joint 5 reference frame via the z axis(90 deg), then via the x axis (90 deg)
#
# W_corr = simplify(W_z * W_x)

# q1_val = 0.92
# q2_val = -0.51
# q3_val = 1.01
# q1_val = .4
# q2_val = .4
# q3_val = -2.95
# q4_val = -1.42
q1_val = 0
q2_val = .5
q3_val = 0
q4_val = 0
q5_val = 0.0
q6_val = 0.0

print("T0_1 = ", T0_1.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("T0_2 = ", T0_2.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("T0_3 = ", T0_3.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("T0_4 = ", T0_4.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("T0_5 = ", T0_5.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("T0_6 = ", T0_6.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("T0_G = ", T0_G.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))

T_total = simplify(T0_G * R_corr)
# T_total = T0_G * R_corr

# print("T_total = ", T_total.evalf(subs={q1: pi/4, q2: pi/4, q3:pi/4, q4:pi/4, q5:pi/4, q6:pi/4}))
print("T_total = ", T_total.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))

# T_total_evaluated = T_total.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val})

px = T_total[0, 3]
py = T_total[1, 3]
pz = T_total[2, 3]

nx = T_total[0, 2]
ny = T_total[1, 2]
nz = T_total[2, 2]

# px = T0_G[0, 3]
# py = T0_G[1, 3]
# pz = T0_G[2, 3]
#
# nx = T0_G[0, 2]
# ny = T0_G[1, 2]
# nz = T0_G[2, 2]

# l = symbols("l")
# l = s[d7]

# T0_5_corr = T0_5 * W_corr

T0_5_corr = T0_5 # * W_corr

l_dist_x = px - T0_5_corr[0, 3]
l_dist_y = py - T0_5_corr[1, 3]
l_dist_z = pz - T0_5_corr[2, 3]


# wx = px - ((s[d6]) * nx)
# wy = py - ((s[d6] + l_dist_y) * ny)
# wz = pz - ((s[d6] + l_dist_z) * nz)

# wx = px - ((s[d7]) * nx)
# wy = py - ((s[d7]) * ny)
# wz = pz - ((s[d7]) * nz)

wx = px - (.303 * nx)
wy = py - (.303 * ny)
wz = pz - (.303 * nz)

# wx = px * nx
# wy = py * ny
# wz = pz * nz

# print("wx ", wx.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))

theta_1 = atan2(wy, wx)

x_adj_theta_1 = cos(theta_1) * s[a1]
y_adj_theta_1 = sin(theta_1) * s[a1]
print("px ", px.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("py ", py.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("pz ", pz.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))

print("wx ", wx.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("wy ", wy.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("wz ", wz.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
# print("lx ", lx.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
# print("ly ", ly.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
# print("lz ", lz.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
# print("lx_j5 ", lx_j5.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
# print("ly_j5 ", ly_j5.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
# print("lz_j5 ", lz_j5.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("nx ", nx.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("ny ", ny.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("nz ", nz.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))


r_coord = sqrt(((wx - x_adj_theta_1) ** 2) + ((wy - y_adj_theta_1) ** 2))
# r_coord = sqrt((wx**2) + (wy**2))
s_coord = wz - s[d1]

# TODO get theta 2 to return q2
theta_2 = atan2(s_coord, r_coord)

x_adj_theta_2 = cos(theta_1 * (s[a1] + (cos(theta_2) * s[a2])))
y_adj_theta_2 = sin(theta_1 * (s[a1] + (sin(theta_2) * s[a2])))
z_adj_theta_2 = sin(theta_2) * s[a2]

r_2_coord = sqrt(((wx - x_adj_theta_2) ** 2) + ((wy - y_adj_theta_2) ** 2))
s_2_coord = wz - s[d1] - z_adj_theta_2

theta_3 = atan2(s_2_coord, r_2_coord)


# print("simplifying theta values")
# theta_1 = simplify(theta_1)
# theta_2 = simplify(theta_2)
# theta_3 = simplify(theta_3)


print("theta_1 ", theta_1)
print("theta_2 ", theta_2)
print("theta_3 ", theta_3)

print("theta_1 ", theta_1.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("theta_2 ", theta_2.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))
print("theta_3 ", theta_3.evalf(subs={q1: q1_val, q2: q2_val, q3:q3_val, q4:q4_val, q5:q5_val, q6:q6_val}))