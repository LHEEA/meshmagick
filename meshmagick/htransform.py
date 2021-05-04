#!/usr/bin/env python
# -*- coding: utf-8 -*-


__author__ = "Francois Rongere"
__copyright__ = "Copyright 2021, D-ICE Engineering"
__credits__ = "Francois Rongere"
__licence__ = "GPLv3"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@dice-engineering.com"

import numpy as np
from math import pi, radians, degrees
from math import sin, cos, acos, asin, atan2, sqrt


# import mesh


def _rotmat_to_cardan(rotmat):
    r00 = rotmat[0, 0]
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
            theta = 0.5 * pi
            psi = -atan2(-r12, r11)
            phi = 0.

    else:  # r20 = +1
        # Not a unique solution: phi+psi = atan2(-r12, r11)
        theta = -0.5 * pi
        psi = atan2(-r12, r11)
        phi = 0.

    return phi, theta, psi


def _cardan_to_rotmat(phi, theta, psi):
    cphi = cos(phi)
    sphi = sin(phi)
    ctheta = cos(theta)
    stheta = sin(theta)
    cpsi = cos(psi)
    spsi = sin(psi)

    rotmat = np.zeros((3, 3), dtype=np.float)

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


class HTransform:

    def __init__(self, rotmat=np.eye(3), trans=np.zeros(3)):
        # TODO: si rotmat 3x3, matrice, sinon, 3x1, angles en radians...
        self.rotmat = np.asarray(rotmat)
        self.transvec = np.asarray(trans)
        assert self.rotmat.shape == (3, 3)
        assert self.transvec.shape == (3,)


    @property
    def trans(self):
        return self.transvec

    @trans.setter
    def trans(self, t):
        self.transvec = np.asarray(t)
        assert self.transvec.shape == (3,)

    @property
    def rot(self):
        return self.rotmat

    @rot.setter
    def rot(self, R):
        self.rotmat = np.asarray(R)
        assert self.rotmat.shape == (3, 3)

    @property
    def cardan_rad(self):
        return _rotmat_to_cardan(self.rotmat)

    @cardan_rad.setter
    def cardan_rad(self, angles_rad):
        phi_rad, theta_rad, psi_rad = angles_rad
        self.rotmat = _cardan_to_rotmat(phi_rad, theta_rad, psi_rad)

    @property
    def cardan_deg(self):
        return map(degrees, self.cardan_rad)

    @cardan_deg.setter
    def cardan_deg(self, angles_deg):
        self.cardan_rad = map(radians, angles_deg)

    # def set_cardan(self, phi_rad, theta_rad, psi_rad, trans):
    #     self.rotmat = _cardan_to_rotmat(phi_rad, theta_rad, psi_rad)
    #     self.transvec = trans

    # def set_rotmat(self, rotmat, trans):
    #     self.rotmat = np.asarray(rotmat)
    #     self.trans = np.asarray(trans)
    #     assert self.rotmat.shape == (3, 3)
    #     assert self.trans.shape == (3,)

    # def get_cardan(self):
    #     return _rotmat_to_cardan(self.rotmat)

    def inverse(self):
        t = HTransform()
        t.set_rotmat(self.rotmat.transpose(), -np.dot(self.rotmat.transpose(), self.transvec))
        return t

    def to_matrix44(self):
        mat = np.zeros((4, 4))
        mat[:3, :3] = self.rotmat
        mat[:3, 3] = self.transvec
        mat[3, :] = [0, 0, 0, 1]
        return mat

    def __mul__(self, other):
        if isinstance(other, HTransform):
            return HTransform(np.dot(self.rotmat, other.rotmat),
                              np.dot(self.rotmat, other.transvec) + self.transvec)

        elif isinstance(other, (np.ndarray, list,)):
            other = np.asarray(other)

            if other.ndim == 1:
                assert other.shape == (3,)
                return np.dot(self.rotmat, other) + self.transvec
            elif other.ndim == 2:
                assert other.shape[1] == 3
                return np.einsum('ij, kj -> kj', self.rotmat, other) + self.transvec  # FIXME: verifier !!
        else:
            assert False

    def __repr__(self):
        return self.to_matrix44().__repr__()

    def __str__(self):
        return self.to_matrix44().__str__()


if __name__ == '__main__':
    print("Testing homogeneous transform")

    t1 = HTransform()
    t2 = HTransform()
    print(t1)

    angles = [pi / 3, pi / 6, 0]
    t1.set_rotmat(_cardan_to_rotmat(0, 0, pi / 2),
                  [0, 0, 0])
    print(t1)

    print(t1.get_cardan())

    t3 = t1 * t2
    print(t3)

    v = t1 * [1, 1, 0]
    print(v)
    print(t1.inverse() * v)

    # print(t1.inverse())
