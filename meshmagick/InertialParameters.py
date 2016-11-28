#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np
from math import pi

# TODO: preciser egalement une frame d'expression des coords de cog et de point

# TODO: donner la possibilite de jouer avec une densite --> creer une sous-classe ???


class InertialParameters(np.ndarray):
    def __new__(cls, mass, cog, xx, yy, zz, yz, xz, xy, point=None):
        
        mat = np.array([[xx, -xy, -xz],
                        [-xy, yy, -yz],
                        [-xz, -yz, zz]], dtype=np.float)
        obj = mat.view(cls)
        
        obj._mass = float(mass)
        
        assert len(cog) == 3
        obj._cog = np.asarray(cog, dtype=np.float)

        if point is None:
            obj._point = obj._cog
        else:
            assert len(point) == 3
            obj._point = np.asarray(point, dtype=np.float)

        return obj
        
    def _get_huygens_correction(self):
        p_g = self._cog - self._point
        return self._mass * (np.dot(p_g, p_g) * np.eye(3) - np.outer(p_g, p_g))
    
    def is_at_cog(self):
        return np.all(self._point == self._cog)
    
    @property
    def reduction_point(self):
        return self._point

    @reduction_point.setter
    def reduction_point(self, point):
        at_cog = self.at_cog
        assert len(point) == 3
        self._point = np.asarray(point, dtype=np.float)
        self[:] = at_cog + self._get_huygens_correction()
        return
    
    @property
    def gravity_center(self):
        return self._cog
    
    @property
    def at_cog(self):
        inertia = self - self._get_huygens_correction()
        inertia._point = inertia._cog
        return inertia
    
    @property
    def mass(self):
        return self._mass
    
    @property
    def xx(self):
        return self[0, 0]
    
    @property
    def yy(self):
        return self[1, 1]
    
    @property
    def zz(self):
        return self[2, 2]
    
    @property
    def yz(self):
        return -self[1, 2]
    
    @property
    def xz(self):
        return -self[0, 2]
    
    @property
    def xy(self):
        return -self[0, 1]
    
    
    def shift_at_cog(self):
        self[:] = self.at_cog
        self._point = self._cog
        return
    
    # The two following methods are taken from the numpy documentation
    # https://docs.scipy.org/doc/numpy/user/basics.subclassing.html#array-wrap
    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._mass = getattr(obj, '_mass', None)
        self._cog = getattr(obj, '_cog', None)
        self._point = getattr(obj, '_point', None)
    
    def __array_wrap__(self, out_arr, context=None):
        return super(InertialParameters, self).__array_wrap__(out_arr, context)
    
    def __str__(self):
        dtype = 'E'
        width = 12
        precision = 3
        str_repr = 'PRINCIPAL INERTIAL PARAMETERS:\n'
        str_repr += "\tMass: {:< {width}.{precision}{dtype}} kg\n".format(self._mass, width=width,
                                                                          precision=precision, dtype=dtype)
        str_repr += "\tGravity center:  {:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}\n".format(
            self._cog[0], self._cog[1],
                                                                                self._cog[2], dtype=dtype,
                                                                                 width=width, precision=precision)
        str_repr += "\tReduction point: {:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}\n".format(self._point[0],self._point[1],
                                                                               self._point[2], dtype=dtype,
                                                                                 width=width, precision=precision)
        str_repr += "\tInertia matrix:\n"
        str_repr += "\t\t{:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}\n".format(self[0, 0], self[0, 1], self[0, 2], dtype=dtype,
                                                                                 width=width, precision=precision)
        str_repr += "\t\t{:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}\n".format(self[1, 0], self[1, 1], self[1, 2], dtype=dtype,
                                                                                 width=width, precision=precision)
        str_repr += "\t\t{:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}{:< {width}.{precision}{dtype}}\n".format(self[2, 0], self[2, 1], self[2, 2], dtype=dtype,
                                                                                 width=width, precision=precision)
        return str_repr


# Here we define the analytical results for most of the current basic geometries
def cylinder_inertia(radius, height, density=1.):
    vol = pi*radius**2*height
    mass = density * vol
    
    Ixx = Iyy = mass * (radius**2/2. + height**2/12.)
    Izz = mass * radius**2/2.
    
    return InertialParameters(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


if __name__ == '__main__':

    inertia = InertialParameters(1, [0, 0, 0], 1, 2, 3, 4, 5, 6, point=[1, 2, 3])
    print inertia
    # print inertia.reduction_point
    print inertia.at_cog
    
    inertia.reduction_point = [2, 2, 2]
    print inertia
    print inertia.at_cog

    inertia.shift_at_cog()
    print inertia
    print inertia.at_cog
    
    inertia.reduction_point = [2, 5, 9]
    print inertia
    print inertia.at_cog

    cylinder = cylinder_inertia(5, 20, 1025)
    
    print cylinder
    cylinder.reduction_point = [1, 0, -10]
    print cylinder
    print cylinder.at_cog




