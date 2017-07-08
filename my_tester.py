from sample_data import given_values, correct_values
from IK import generate_rrpy_matrix, get_wrist_coordinates
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2, acos

for shelf in given_values:
    roll = given_values[shelf]["roll"]
    pitch = given_values[shelf]["pitch"]
    yaw = given_values[shelf]["yaw"]

    Rrpy = generate_rrpy_matrix(roll, pitch, yaw)

    px = given_values[shelf]["px"]
    py = given_values[shelf]["py"]
    pz = given_values[shelf]["pz"]

    # should be s[d7] + s[d6] or .303, however .197 works better
    dist_wrist_to_effector = 0.197

    wx, wy, wz = get_wrist_coordinates(Rrpy, dist_wrist_to_effector)

    # compute theta for joint 1
    theta1 = atan2(wy, wx)
    expected_theta1 = correct_values[shelf]["joint_1"]

    print(shelf)
    print("expected theta1", expected_theta1)
    print("theta1         ", theta1)
    assert(round(theta1, 2) == round(expected_theta1, 2), "theta1 mismatch")