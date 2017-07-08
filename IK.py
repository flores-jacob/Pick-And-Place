import numpy as np
from numpy import array, ndarray
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2, acos
from sympy.matrices import Matrix

from sample_data import given_values, correct_values

# These are the rads of joint 3 in resting position
# RADS_AT_REST_JOINT3 = 1.60509488746813
# RADS_AT_REST_JOINT3 = 1.60926888235184
RADS_AT_REST_JOINT3 = 1.60678078687695

# These are the rads of joint 2 in resting position
RADS_AT_REST_JOINT2 = 1.57079632679490

RADS_AT_REST_JOINT4 = 1.57079632679490

RADS_AT_REST_JOINT6 = 1.57079632679490
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

    theta_x = atan2(r32, r33)
    theta_y = atan2(-r31, sqrt((r32 ** 2) + (r33 ** 2)))
    theta_z = atan2(r21, r11)

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

    # Rot_X = rot_gamma(roll)[0:3, 0:3]
    # Rot_Y = rot_beta(pitch)[0:3, 0:3]
    # Rot_Z = rot_alpha(yaw)[0:3, 0:3]

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

# Gripper link  orientation correction so that URDF values are in accordance with DH Convention

print("computing orientations")

# Rotate z by 180 degrees
R_z = Matrix([
    [cos(pi), -sin(pi), 0, 0],
    [sin(pi), cos(pi), 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

# Rotate along the y axis by -90 degrees
R_y = Matrix([
    [cos(-pi / 2), 0, sin(-pi / 2), 0],
    [0, 1, 0, 0],
    [-sin(-pi / 2), 0, cos(-pi / 2), 0],
    [0, 0, 0, 1]
])

# R_corr_4_6_and_gripper = simplify(R_z * R_y)
R_corr_4_6_and_gripper = R_z * R_y


# Rotate along the z axis by 90 degrees
R_z90 = Matrix([
    [cos(pi / 2), -sin(pi / 2), 0, 0],
    [sin(pi / 2), cos(pi / 2), 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

R_x = Matrix([
    [1, 0, 0, 0],
    [0, cos(pi / 2), -sin(pi / 2), 0],
    [0, sin(pi / 2), cos(pi / 2), 0],
    [0, 0, 0, 1]
])

# R_corr_3_and_5 = simplify(R_z90 * R_x)
R_corr_3_and_5 = R_z90 * R_x

# compose rotation matrix
# procedure taken from http://nghiaho.com/?page_id=846


# insert test cases and given data here the first objective is to produce Rrpy that is consistent with
# the one in Gazebo we will first obtain the givens in gazebo, roll, pitch, yaw and position coords we will
# then proceed to construct the Rrpy

# generate_rot_matrix = yaw * pitch * roll


# Rrpy for shelf 4
Rrpy = np.asarray([[  9.99999455e-01,  -9.98744949e-04,   3.02840662e-04,
          0.00000000e+00],
       [  9.98888264e-04,   9.99999389e-01,  -4.73457027e-04,
          0.00000000e+00],
       [ -3.02367614e-04,   4.73759274e-04,   9.99999842e-01,
          0.00000000e+00],
       [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          1.00000000e+00]])

px = 2.0899794436124215
py = 0.9001324988518173
pz = 1.5810668259285738




# Rrpy for shelf 4 ryzy version
# computed using ryzy
Rrpy = np.asarray([[  9.99999681e-01,   4.58185713e-04,  -6.54029113e-04,
          0.00000000e+00],
       [ -4.58185613e-04,   9.99999895e-01,   3.02519797e-07,
          0.00000000e+00],
       [  6.54029183e-04,  -2.85297069e-09,   9.99999786e-01,
          0.00000000e+00],
       [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          1.00000000e+00]])

px = 2.0900268663882255
py = 0.9000621976697936
pz = 1.5810337949811457


Rrpy = np.asarray([
[    0.999999677064026, 0.00066025299681594, -0.000458189724437495, 0],
[-0.000660255780493678,   0.999999782013581,  -5.92414765511909e-6, 0],
[ 0.000458185713122114, 6.22666815612136e-6,     0.999999895013535, 0],
[                    0,                   0,                     0, 1]])




# Rrpy = array([[  9.99999947e-01,   7.64667591e-05,  -3.15630087e-04,
#           0.00000000e+00],
#        [ -7.67083401e-05,   9.99999704e-01,  -7.65451711e-04,
#           0.00000000e+00],
#        [  3.15571462e-04,   7.65475882e-04,   9.99999657e-01,
#           0.00000000e+00],
#        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#           1.00000000e+00]])
# px = 2.0900516509914198
# py = 0.9000068604315312
# pz = 1.5809190609211612

# Rrpy = array([[ 0.99323148, -0.07350677, -0.08993326,  0.        ],
#        [ 0.07069361,  0.99691572, -0.03408009,  0.        ],
#        [ 0.092161  ,  0.02749172,  0.99536453,  0.        ],
#        [ 0.        ,  0.        ,  0.        ,  1.        ]])
# px = 2.1296616660546173
# py = 0.6032108013179218
# pz = 1.7270503654194076
# ('wx ', 1.828712528984358)
# ('wy ', 0.58179063779495377)
# ('wz ', 1.6991255826751228)
# ('theta1 ', 1.26277951687629)
# ('theta2 ', 2.97607953081526 - 0.701546950012246*I)
# ('theta3 ', -3.17757711367185 + 2.60862850712765*I)
# ('theta4 ', 0)
# ('theta5 ', 0)
# ('theta6 ', 0)

# data for shelf 9
# Rrpy = np.asarray([[  9.99999447e-01,   5.96149636e-04,   8.66090459e-04,
#           0.00000000e+00],
#        [ -5.96493507e-04,   9.99999743e-01,   3.96833721e-04,
#           0.00000000e+00],
#        [ -8.65853665e-04,  -3.97350119e-04,   9.99999546e-01,
#           0.00000000e+00],
#        [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#           1.00000000e+00]])
# px = 2.0899788602474105
# py = 0.9000462046181804
# pz = 1.5811458368578455


roll = given_values["shelf_4"]["roll"]
pitch = given_values["shelf_4"]["pitch"]
yaw = given_values["shelf_4"]["yaw"]

Rot_X = rotate_x(roll)
Rot_Y = rotate_y(pitch)
Rot_Z = rotate_z(yaw)

Rrpy = Rot_Z * Rot_Y * Rot_X

print("generated matrix", Rrpy)



lx = Rrpy[0, 0]
ly = Rrpy[1, 0]
lz = Rrpy[2, 0]

dist_wrist_to_effector = s[d7]
d6_val = s[d6]

# d6 is from the lesson. we use lx, ly, lz for an transform on the x axis instead of nx, ny, nz which is a transform
# across the z axis
wx = px - (.193 * lx)
wy = py - (.193 * ly)
wz = pz - (.193 * lz)

# compute theta for joint 1
theta_1 = pi/2 - atan2(wx, wy)
theta_1 = theta_1.evalf()
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

# print("theta1 ", theta1)
# print("expected theta1 ", correct_values["shelf_4"]["joint_1"])
# print("theta2 ", theta2)
# print("theta3 ", theta3)

R0_3 = T0_3[0:3, 0:3]

R0_3_corr = R0_3 * R_corr_3_and_5[0:3, 0:3]
# R0_3_corr = (R0_3 * rotate_x(-pi/2)[0:3, 0:3]) * rotate_z(-pi/2)[0:3, 0:3]
R0_3_corr = R0_3_corr.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

R0_3 = R0_3_corr.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

Rrpy_corr = Rrpy[0:3, 0:3] * R_corr_4_6_and_gripper[0:3, 0:3].evalf()

R3_6 = R0_3.inv() * Rrpy_corr

R3_6 = (R3_6 * rotate_x(-pi/2)[0:3, 0:3]) * rotate_z(-pi/2)[0:3, 0:3]

R3_6_corr = (R3_6 * rotate_y(pi / 2)[0:3, 0:3]) * rotate_z(-pi / 2)[0:3, 0:3]

# Why -pitch for theta4 and roll for theta5? When theta4 is roll, and theta5 is pitch? This is because these are the
# are the formulas that return the desired values base on the different permutations we have tried
theta4 = get_pitch(R3_6_corr).evalf(subs={q1: theta1, q2: theta2, q3: theta3})
theta5 = -get_roll(R3_6_corr).evalf(subs={q1: theta1, q2: theta2, q3: theta3})
# theta6 is simply the revers of theta4
theta6 = -theta4

# print("theta4 ", theta4)
# print("theta5 ", theta5)
# print("theta6 ", theta6)