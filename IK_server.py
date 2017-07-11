#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *

import numpy as np

import json, os.path


# These are the rads of joint 3 in resting position
RADS_AT_REST_JOINT3 = 1.60678078687695

# These are the rads of joint 2 in resting position
RADS_AT_REST_JOINT2 = pi/2 # 1.57079632679490


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
    # compose rotation matrix
    # procedure taken from http://nghiaho.com/?page_id=846
    rot_x = rotate_x(roll)[0:3, 0:3]
    rot_y = rotate_y(pitch)[0:3, 0:3]
    rot_z = rotate_z(yaw)[0:3, 0:3]

    Rrpy = rot_z * rot_y * rot_x

    return Rrpy


# roll function for forward kinematics
def get_roll(rotation_matrix):
    r32 = rotation_matrix[2, 1]
    r33 = rotation_matrix[2, 2]

    roll_val = atan2(r32, r33)

    return roll_val


# pitch function for forward kinematics
def get_pitch(rotation_matrix):
    r31 = rotation_matrix[2, 0]
    r32 = rotation_matrix[2, 1]
    r33 = rotation_matrix[2, 2]

    pitch_val = atan2(-r31, sqrt((r32 ** 2) + (r33 ** 2)))

    return pitch_val


# yaw function for forward kinematics
def get_yaw(rotation_matrix):
    # solution adapted from this page http://nghiaho.com/?page_id=846
    r11 = rotation_matrix[0, 0]
    r21 = rotation_matrix[1, 0]

    yaw_val = atan2(r21, r11)

    return yaw_val


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

T0_2 = T0_1 * T1_2
T0_3 = T0_2 * T2_3
T0_4 = T0_3 * T3_4
T0_5 = T0_4 * T4_5
T0_6 = T0_5 * T5_6
T0_G = T0_6 * T6_G

T_total = T0_G * rotate_z(pi) * rotate_y(-pi/2)


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []

        ######### Start code to open json error file ###########

        # open the saved error list. if file does not exist, create it
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

        ######### End code to open json error file ###########


        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Define DH param symbols



            # Joint angle symbols



            # Modified DH params



            # Define Modified DH Transformation matrix



            # Create individual transformation matrices



            # Extract end-effector position and orientation from request
	        # px,py,pz = end-effector position
	        # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            # create Rrpy matrix
            Rrpy = generate_rrpy_matrix(roll, pitch, yaw)

            dist_wrist_to_effector = s[d7]
            d6_val = s[d6]

            # compute wrist values
            wx, wy, wz = get_wrist_coordinates(Rrpy, px, py, pz, d6_val + dist_wrist_to_effector)

            # Calculate joint angles using Geometric IK method
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
            unadjusted_theta_3 = get_theta_3(distance_joint_2to3, distance_joint_3to5, wx, wz, theta1,
                                             x_offset_to_joint2,
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

            # get the matrix for base link to joint 3, using computed theta1, 2, and 3 values
            R0_3 = T0_3[0:3, 0:3]
            R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
            # convert sympy matrix to numpy matrix to avoid errors
            # how to convert lifted from: https://stackoverflow.com/a/37491889
            R0_3 = np.array(R0_3).astype(np.float64)

            # correct the orientation of the Rrpy matrix
            Rrpy = Rrpy * rotate_y(pi / 2)[0:3, 0:3] * rotate_z(pi)[0:3, 0:3]
            # get the matrix values for the wrist
            R3_6 = (np.linalg.inv(R0_3)) * Rrpy

            # compute values for thetas 4, 5, and 6
            theta4 = atan2(R3_6[2, 2], -R3_6[0, 2])
            theta4 = theta4.evalf()
            theta5 = acos(R3_6[1, 2])
            theta5 = theta5.evalf()
            theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])
            theta6 = theta6.evalf()

            # compute forward kinematics for error checking
            FK_matrix = T_total.evalf(
                subs={q1: theta1, q2: theta2, q3: theta3, q4: theta4, q5: theta5, q6: theta6})

            FK_px = FK_matrix[0, 3]
            FK_py = FK_matrix[1, 3]
            FK_pz = FK_matrix[2, 3]
            FK_roll = get_roll(FK_matrix)
            FK_pitch = get_pitch(FK_matrix)
            FK_yaw = get_yaw(FK_matrix)

            # compute the error values for the px, py, pz and rpy
            px_error = (px - FK_px)  # ** 2
            py_error = (py - FK_py)  # ** 2
            pz_error = (pz - FK_pz)  # ** 2

            roll_error = (roll - FK_roll)  # ** 2
            pitch_error = (pitch - FK_pitch)  # ** 2
            yaw_error = (yaw - FK_yaw)  # ** 2

            # append each of the error values to their respective lists
            px_error_list.append(float(px_error))
            py_error_list.append(float(py_error))
            pz_error_list.append(float(pz_error))
            roll_error_list.append(float(roll_error))
            pitch_error_list.append(float(pitch_error))
            yaw_error_list.append(float(yaw_error))



            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        # once once done with IK and FK operations, save the updated set of errors to disk
        with open(error_list_path, "w") as error_list_file:
            json.dump(json_error_list, error_list_file, encoding='utf-8')

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
