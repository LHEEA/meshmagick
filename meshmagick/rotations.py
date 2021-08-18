#!/usr/bin/env python
#  -*- coding: utf-8 -*-

# Written by Francois Rongere (D-ICE Engineering)
#
# Octobre 2018
#

import numpy as np
from math import pi, radians, degrees
from math import sin, cos, acos, asin, atan2, sqrt

norm = np.linalg.norm
_pi_2 = 0.5*pi

"""
This module aims at converting 3D rotations between the following representations:

    * Unit quaternion
    * Rotation matrix
    * Axis angle
    * Cardan angles

Functions are written for any pair of representation, bidirectionally.

Every angles are expected in radians
Quaternions are numpy arrays of shape (4, 1) but functions that accept a quaternion as input may be any iterable

Note that no check is done yet on the rotation nature of rotation matrices

You may launch the test suite by running the test function in a python console or by lanching the script in a terminal,
eventually with the number of random rotation as a command line argument. Default is 1000.


Please report any bug you may found to this file to:

Francois Rongere <Francois.Rongere@dice-engineering.com>

This file is under the GPL v3 license.
Please report to : https://www.gnu.org/licenses/gpl-3.0.fr.html

Please contact Francois Rongere if any need of integration of this file into other product.
"""


def is_unit_vector(vector):
    vector = np.asarray(vector, dtype=float)
    return np.isclose(norm(vector), 1.)

def quat_to_axis_angle(quat):
    assert is_unit_vector(quat)
    assert quat.size == 4
    quat = np.asarray(quat, dtype=float)
    q0 = quat[0]
    d = 1. - q0*q0

    if np.isclose(d, 0.):
        angle = 0.
        axis = np.array((0, 0, 1), dtype=float)

    else:
        angle = 2. * acos(q0)
        axis = quat[1:] / sqrt(d)

    return axis, angle

def axis_angle_to_quat(axis, angle):
    assert is_unit_vector(axis)

    half_angle = 0.5 * angle
    c = cos(half_angle)
    s = sin(half_angle)

    return np.array((c, s*axis[0], s*axis[1], s*axis[2]), dtype=float)

def quat_to_rotmat(quat):
    assert is_unit_vector(quat)
    assert quat.size == 4
    quat = np.asarray(quat, dtype=float)

    q0q0 = quat[0] * quat[0]
    q1q1 = quat[1] * quat[1]
    q2q2 = quat[2] * quat[2]
    q3q3 = quat[3] * quat[3]
    q0q1 = quat[0] * quat[1]
    q0q2 = quat[0] * quat[2]
    q0q3 = quat[0] * quat[3]
    q1q2 = quat[1] * quat[2]
    q1q3 = quat[1] * quat[3]
    q2q3 = quat[2] * quat[3]

    rotmat = np.zeros((3, 3), dtype=float)
    rotmat[0, 0] = (q0q0 + q1q1) * 2 - 1
    rotmat[0, 1] = (q1q2 - q0q3) * 2
    rotmat[0, 2] = (q1q3 + q0q2) * 2
    rotmat[1, 0] = (q1q2 + q0q3) * 2
    rotmat[1, 1] = (q0q0 + q2q2) * 2 - 1
    rotmat[1, 2] = (q2q3 - q0q1) * 2
    rotmat[2, 0] = (q1q3 - q0q2) * 2
    rotmat[2, 1] = (q2q3 + q0q1) * 2
    rotmat[2, 2] = (q0q0 + q3q3) * 2 - 1

    return rotmat

def rotmat_to_quat(rotmat):

    # No test here for the properties of the matrix (orthogonal, symmetric...)

    # Real s, tr;
    # Real half = (Real)0.5;

    r00 = rotmat[0, 0]
    r01 = rotmat[0, 1]
    r02 = rotmat[0, 2]
    r10 = rotmat[1, 0]
    r11 = rotmat[1, 1]
    r12 = rotmat[1, 2]
    r20 = rotmat[2, 0]
    r21 = rotmat[2, 1]
    r22 = rotmat[2, 2]

    tr = r00 + r11 + r22

    if tr >= 0.:
        s = sqrt(tr + 1)
        q0 = 0.5 * s
        s = 0.5 / s
        q1 = (r21 - r12) * s
        q2 = (r02 - r20) * s
        q3 = (r10 - r01) * s

    else:
        i = 0

        if r11 > r00:
            i = 1
            if r22 > r11:
                i = 2
        else:
            if r22 > r00:
                i = 2

        if i == 0:
            s = sqrt(r00 - r11 - r22 + 1)
            q1 = 0.5 * s
            s = 0.5 / s
            q2 = (r01 + r10) * s
            q3 = (r20 + r02) * s
            q0 = (r21 - r12) * s

        elif i == 1:
            s = sqrt(r11 - r22 - r00 + 1)
            q2 = 0.5 * s
            s = 0.5 / s
            q3 = (r12 + r21) * s
            q1 = (r01 + r10) * s
            q0 = (r02 - r20) * s

        elif i == 2:
            s = sqrt(r22 - r00 - r11 + 1)
            q3 = 0.5 * s
            s = 0.5 / s
            q1 = (r20 + r02) * s
            q2 = (r12 + r21) * s
            q0 = (r10 - r01) * s

    return np.array((q0, q1, q2, q3), dtype=float)

def cardan_to_quat(phi, theta, psi):

    phi_2 = 0.5 * phi
    theta_2 = 0.5 * theta
    psi_2 = 0.5 * psi

    cphi_2 = cos(phi_2)
    sphi_2 = sin(phi_2)
    ctheta_2 = cos(theta_2)
    stheta_2 = sin(theta_2)
    cpsi_2 = cos(psi_2)
    spsi_2 = sin(psi_2)

    q0 = cphi_2 * ctheta_2 * cpsi_2 + sphi_2 * stheta_2 * spsi_2
    q1 = -cphi_2 * stheta_2 * spsi_2 + ctheta_2 * cpsi_2 * sphi_2
    q2 = cphi_2 * cpsi_2 * stheta_2 + sphi_2 * ctheta_2 * spsi_2
    q3 = cphi_2 * ctheta_2 * spsi_2 - sphi_2 * cpsi_2 * stheta_2

    return np.array((q0, q1, q2, q3), dtype=float)

def rotmat_to_cardan(rotmat):

    r00 = rotmat[0, 0]
    # r01 = rotmat[0, 1]
    # r02 = rotmat[0, 2]
    r10 = rotmat[1, 0]
    r11 = rotmat[1, 1]
    r12 = rotmat[1, 2]
    r20 = rotmat[2, 0]
    r21 = rotmat[2, 1]
    r22 = rotmat[2, 2]

    if r20 < 1.:
        if r20 > -1.:
            theta = asin(-r20)
            psi = atan2(r10, r00)
            phi = atan2(r21, r22)

        else:  # r20 = -1
            # Not a unique solution: phi -psi = atan2(-r12, r11)
            theta = _pi_2
            psi = -atan2(-r12, r11)
            phi = 0.

    else:  # r20 = +1
        # Not a unique solution: phi+psi = atan2(-r12, r11)
        theta = -_pi_2
        psi = atan2(-r12, r11)
        phi = 0.

    return phi, theta, psi

def cardan_to_rotmat(phi, theta, psi):

    cphi = cos(phi)
    sphi = sin(phi)
    ctheta = cos(theta)
    stheta = sin(theta)
    cpsi = cos(psi)
    spsi = sin(psi)

    rotmat = np.zeros((3, 3), dtype=float)

    rotmat[0, 0] = ctheta * cpsi
    rotmat[0, 1] = sphi * stheta * cpsi - cphi * spsi
    rotmat[0, 2] = cphi * stheta * cpsi + sphi * spsi
    rotmat[1, 0] = ctheta * spsi
    rotmat[1, 1] = sphi * stheta * spsi + cphi * cpsi
    rotmat[1, 2] = cphi * stheta * spsi - sphi * cpsi
    rotmat[2, 0] = -stheta
    rotmat[2, 1] = ctheta * sphi
    rotmat[2, 2] = ctheta * cphi

    return rotmat

def quat_to_cardan(quat):
    return rotmat_to_cardan(quat_to_rotmat(quat))

def rotmat_to_axis_angle(rotmat):
    return quat_to_axis_angle(rotmat_to_quat(rotmat))

def axis_angle_to_rotmat(axis, angle):
    return quat_to_rotmat(axis_angle_to_quat(axis, angle))

def cardan_to_axis_angle(phi, theta, psi):
    return quat_to_axis_angle(cardan_to_quat(phi, theta, psi))

def axis_angle_to_cardan(axis, angle):
    return quat_to_cardan(axis_angle_to_quat(axis, angle))


def test():
    # Random quaternion to initialize tests
    quat = np.random.rand(4)
    quat /= norm(quat)
    print("\n0 - Initial unit quaternion")
    print(quat)

    axis, angle = quat_to_axis_angle(quat)
    print("\n1 - quat to axis angle")
    print(axis, angle)

    quat1 = axis_angle_to_quat(axis, angle)
    print("\n2 - axis angle to quat")
    print(quat1)
    assert np.all(np.isclose(quat1, quat))

    rotmat = quat_to_rotmat(quat1)
    print("\n3 - quat to rotmat")
    print(rotmat)

    quat2 = rotmat_to_quat(rotmat)
    print("\n4 - rotmat to quat")
    print(quat2)
    assert np.all(np.isclose(quat2, quat))

    phi, theta, psi = rotmat_to_cardan(rotmat)
    print("\n5 - rotmat to cardan")
    print(phi, theta, psi)

    quat3 = cardan_to_quat(phi, theta, psi)
    print("\n6 - cardan to quat")
    print(quat3)
    assert np.all(np.isclose(quat3, quat))

    rotmat2 = cardan_to_rotmat(phi, theta, psi)
    print("\n7 - cardan to rotmat")
    print(rotmat2)
    assert np.all(np.isclose(rotmat2, rotmat))

    phi2, theta2, psi2 = quat_to_cardan(quat3)
    print("\n8 - quat to cardan")
    print(phi2, theta2, psi2)
    assert np.all(np.isclose((phi, theta, psi), (phi2, theta2, psi2)))

    axis2, angle2 = rotmat_to_axis_angle(rotmat2)
    print("\n9 - rotmat to axis angle")
    print(axis2, angle2)
    assert np.all(np.isclose(axis, axis2))
    assert np.isclose(angle2, angle)

    rotmat3 = axis_angle_to_rotmat(axis2, angle2)
    print("\n10 - axis angle to rotmat")
    print(rotmat3)
    assert np.all(np.isclose(rotmat3, rotmat))

    axis3, angle3 = cardan_to_axis_angle(phi2, theta2, psi2)
    print("\n11 - cardan to axis angle")  # FIXME ---> PB
    print(axis3, angle3)
    assert np.all(np.isclose(axis, axis3))
    assert np.isclose(angle, angle3)

    phi3, theta3, psi3 = axis_angle_to_cardan(axis3, angle3)
    print("\n12 - axis angle to cardan")
    print(phi3, theta3, psi3)
    assert np.all(np.isclose((phi, theta, psi), (phi3, theta3, psi3)))


if __name__ == '__main__':
    import sys
    try:
        n = int(sys.argv[1])
    except IndexError:
        n = 1000
    print("===============================================================")
    print("TESTING ROTATION CONVERSIONS ON %u RANDOM ROTATIONS" % n)
    print("===============================================================")

    iter = 1
    while iter <= n:
        print("\n\nITERATION %u :" % iter)
        test()
        iter += 1

    print("====================")
    print("-- TESTING PASSED --")
    print("====================")
