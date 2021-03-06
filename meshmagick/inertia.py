#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""This module defines a RigidBodyInertia class to handle 3D rotational inertia of rigid bodies"""

import numpy as np
from math import pi, sqrt
from copy import deepcopy

from . import densities

# FIXME: attention, changer les signes pour les produits d'inertie !
# TODO: ajouter la production d'inerties de solides connus --> utile pour comparaison !!
# TODO: indiquer une frame d'expression ...


class RigidBodyInertia(object):
    """
    Parameters
    ----------
    mass : float
        The mass of the body in kg
    cog : array_like
        The 3D coordinates of the center of gravity
    xx : float
        The principal inertia moment around x axis
    yy : float
        The principal inertia moment around y axis
    zz : float
        The principal inertia moment around z axis
    yz : float
        Inertia product
    xz : float
        Inertia product
    xy : float
        Inertia product
    point : array_like, optional
        The reduction point. If None, it is set to be cog
    """
    def __init__(self, mass, cog, xx, yy, zz, yz, xz, xy, point=None):
        
        self._mass = float(mass)
        
        assert len(cog) == 3
        self._cog = np.asarray(cog, dtype=float)
        
        if point is None:
            self._point = self._cog
        else:
            assert len(point) == 3
            self._point = np.asarray(point, dtype=float)
        
        self._3d_rotational_inertia = np.array([[float(xx), -float(xy), -float(xz)],
                                                [-float(xy), float(yy), -float(yz)],
                                                [-float(xz), -float(yz), float(zz)]], dtype=float)
    
    @property
    def mass(self):
        """The mass of the body"""
        return self._mass
    
    @property
    def gravity_center(self):
        """The position of the center of gravity"""
        return self._cog
    
    @property
    def inertia_matrix(self):
        """The 3D rotational inertia matrix"""
        return self._3d_rotational_inertia
    
    @property
    def reduction_point(self):
        """The reduction point of the inertia matrix
        
        Returns
        -------
        ndarray
        """
        return self._point
    
    @reduction_point.setter
    def reduction_point(self, point):
        """Set the reduction point"""
        mat_at_cog = self.at_cog.inertia_matrix
        assert len(point) == 3
        self._point = np.asarray(point, dtype=float)
        self._3d_rotational_inertia = mat_at_cog + self._huygens_transport()
    
    @property
    def at_cog(self):
        """Returns a new inertia object that is expressed at cog.
        
        It makes a copy of itself.
        
        Returns
        -------
        ndarray
        """
        inertia = deepcopy(self)
        inertia.shift_at_cog()
        return inertia
    
    def shift_at_cog(self):
        """Shift the inertia matrix internally at cog.
        
        The reduction point is then cog.
        """
        self._3d_rotational_inertia -= self._huygens_transport()
        self._point = self._cog
    
    def is_at_cog(self):
        """Returns whether the object is expressed at cog
        
        Returns
        -------
        bool
        """
        return np.all(self._point == self._cog)
    
    def _huygens_transport(self):
        p_g = self._cog - self._point
        return self._mass * (np.dot(p_g, p_g) * np.eye(3) - np.outer(p_g, p_g))
    
    @property
    def xx(self):
        """Get the principal inertia moment around x
        
        Returns
        -------
        float
        """
        return self._3d_rotational_inertia[0, 0]
    
    @property
    def yy(self):
        """Get the principal inertia moment around y

        Returns
        -------
        float
                """
        return self._3d_rotational_inertia[1, 1]
    
    @property
    def zz(self):
        """Get the principal inertia moment around z

        Returns
        -------
        float
        """
        return self._3d_rotational_inertia[2, 2]
    
    @property
    def yz(self):
        """Get the yz inertia product

        Returns
        -------
        float
        """
        return -self._3d_rotational_inertia[1, 2]
    
    @property
    def xz(self):
        """Get the xz inertia product

        Returns
        -------
        float
        """
        return -self._3d_rotational_inertia[0, 2]
    
    @property
    def xy(self):
        """Get the xy inertia product

        Returns
        -------
        float
        """
        return -self._3d_rotational_inertia[0, 1]
    
    def __str__(self):
        width = 15
        precision = 6
        dtype = 'E'
        fformat = '< {width:}.{precision:}{dtype:}'.format(width=width, precision=precision, dtype=dtype)
        
        str_repr = '\nPincipal inertial parameters:\n'
        str_repr += '\tMass: {mass:{fformat}} kg\n'.format(fformat=fformat, mass=self._mass)
        
        str_repr += '\tCOG:  {cog[0]:{fformat}}{cog[1]:{fformat}}{cog[2]:{fformat}}\n'.format(
            fformat=fformat, cog=self._cog)
        
        str_repr += '\tInertia matrix expressed at: {p[0]:{fformat}}{p[1]:{fformat}}{p[2]:{fformat}}\n'.format(
            fformat=fformat, p=self._point
        )
        str_repr += '\t\t{i[0]:{fformat}}{i[1]:{fformat}}{i[2]:{fformat}}\n'.format(
            fformat=fformat, i=self._3d_rotational_inertia[0]
        )
        str_repr += '\t\t{i[0]:{fformat}}{i[1]:{fformat}}{i[2]:{fformat}}\n'.format(
            fformat=fformat, i=self._3d_rotational_inertia[1]
        )
        str_repr += '\t\t{i[0]:{fformat}}{i[1]:{fformat}}{i[2]:{fformat}}\n'.format(
            fformat=fformat, i=self._3d_rotational_inertia[2]
        )
        
        return str_repr


# Principal geometrical shapes
# From "Handbook of equations for mass and area of various geometrical shapes, J.A. Myers, 1962"
def right_circular_cylinder(radius, length, density=1.):
    """Get the inertia of a right circular cylinder
    
    Returns
    -------
    RigidBodyInertia
    """
    
    vol = pi * radius ** 2 * length
    mass = density * vol
    
    Ixx = Iyy = mass * (3 * radius ** 2 + length ** 2) / 12.
    Izz = mass * radius ** 2 / 2.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def hollow_right_circular_cylinder(int_radius, ext_radius, length, density=1.):
    """Get the inertia of a hollow right circular cylinder

    Returns
    -------
    RigidBodyInertia
    """
    
    vol = pi * length * (ext_radius ** 2 - int_radius ** 2)
    mass = density * vol
    
    R2r2 = ext_radius ** 2 + int_radius ** 2
    Ixx = Iyy = mass * (3 * R2r2 + length ** 2) / 12.
    Izz = mass * R2r2 / 2.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def right_circular_cone(radius, length, density=1.):
    """Get the inertia of a right circular cone

    Returns
    -------
    RigidBodyInertia
    
    Note
    ----
    The center of gravity is at an altitude of z = H/4 over the circular basis center
    """
    
    vol = pi * radius ** 2 * length / 3.
    mass = density * vol
    
    Ixx = Iyy = 3 * mass * (radius ** 2 + length ** 2/4.) / 20.
    Izz = 3. * mass * radius ** 2 / 10.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def sphere(radius, density=1.):
    """Get the inertia of a sphere

    Returns
    -------
    RigidBodyInertia
    """
    
    vol = 4. * pi * radius ** 3 / 3.
    mass = density * vol
    
    Ixx = Iyy = Izz = 2 * mass * radius ** 2 / 5.
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0., 0., 0.)


def hollow_sphere(int_radius, ext_radius, density=1.):
    """Get the inertia of a hollow sphere

    Returns
    -------
    RigidBodyInertia
    """
    
    vol = 4. * pi * (ext_radius ** 3 - int_radius ** 3) / 3.
    mass = density * vol
    
    Ixx = Iyy = Izz = 2 * mass * (ext_radius ** 5 - int_radius ** 5) / (ext_radius ** 3 - int_radius ** 3) / 5
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0., 0., 0.)


def hemisphere(radius, density=1.):
    """Get the inertia of a hemisphere

    Returns
    -------
    RigidBodyInertia
    
    Note
    ----
    The center of gravity is situated at the altitude of z = 3R/8 over the circular basis center
    """
    
    vol = 2 * pi * radius ** 3 / 3.
    mass = density * vol
    
    Ixx = Iyy = Izz = 0.26 * mass * radius ** 2
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def elliptical_cylinder(a, b, length, density=1.):
    """Get the inertia of an elliptical cylinder

    Returns
    -------
    RigidBodyInertia
    
    Note
    ----
    * The center of gravity is located at an altitude of z=H/2 over the elliptical basis center
    * a is along x axis (ellipse semi axis)
    * b is along y axis (ellipse semi axis)
    * length is along z axis
    """
    
    vol = pi * a * b * length
    mass = density * vol
    
    Ixx = mass * (3 * b ** 2 + length ** 2) / 12.
    Iyy = mass * (3 * a ** 2 + length ** 2) / 12.
    Izz = mass * (a ** 2 + b ** 2) / 4.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def ellipsoid(a, b, c, density=1.):
    """Get the inertia of an ellipsoid

    Returns
    -------
    RigidBodyInertia
    
    Note
    ----
    * a is along z axis (ellipse semi axis)
    * b is along x axis (ellipse semi axis)
    * c is along y axis (ellipse semi axis)
    """
    
    vol = 4 * pi * a * b * c / 3.
    mass = density * vol
    
    Ixx = mass * (a ** 2 + c ** 2) / 5.
    Iyy = mass * (a ** 2 + b ** 2) / 5.
    Izz = mass * (b ** 2 + c ** 2) / 5.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def torus(chord_radius, tube_radius, density=1.):
    """Get the inertia of a torus

    Returns
    -------
    RigidBodyInertia
    """
    
    vol = 2*pi**2 * tube_radius**2 * chord_radius
    mass = density * vol
    
    Ixx = Izz = mass * (4*chord_radius**2 + 5*tube_radius**2) / 8.
    Iyy = mass * (4*chord_radius**2 + 3*tube_radius**2) / 4.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def right_angle_wedge(base, height, length, density=1.):
    """Get the inertia of a right angle wedge

    Returns
    -------
    RigidBodyInertia
    """
    
    vol = base*height*length / 2.
    mass = density * vol
    
    Ixx = mass * (2*height**2 + 3*length**2) / 36.
    Iyy = mass * (base**2 + height**2) / 18.
    Izz = mass * (2*base**2 + 3*length**2) / 36.

    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def isoceles_wedge(base, height, length, density=1.):
    """Get the inertia of an isocele wedge

    Returns
    -------
    RigidBodyInertia
    """
    
    vol = base * height * length / 2.
    mass = density * vol

    Ixx = mass * (2 * height**2 + 3 * length**2) / 36.
    Iyy = mass * (4*height**2 + 3*base**2) / 72.
    Izz = mass * (2 * length**2 + base**2) / 24.

    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def right_rectangular_pyramid(a, b, height, density=1.):
    """Get the inertia of a right rectangular pyramid

    Returns
    -------
    RigidBodyInertia
    
    Note
    ----
    The center of gravity is located at the altitude z=H/4 over the rectangular basis center
    """
    
    vol = a*b*height / 3.
    mass = density * vol
    
    Ixx = mass * (b**2 + 3*height**2/4.) / 20.
    Iyy = mass * (a**2 + 3*height**2/4.) / 20.
    Izz = mass * (a**2 + b**2) / 20.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def cube(a, density=1.):
    """Get the inertia of a cube

    Returns
    -------
    RigidBodyInertia
    """
    
    vol = a**3
    mass = density * vol
    
    Ixx = Iyy = Izz = mass * a**2 / 6.
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def rectangular_prism(a, b, h, density=1.):
    """Get the inertia of a rectangular prism

    Returns
    -------
    RigidBodyInertia
    
    Note
    ----
    * a is along x
    * b is along y
    * h is along z
    """
    
    vol = a*b*h
    mass = density * vol
    
    Ixx = mass * (b**2 + h**2) / 12.
    Iyy = mass * (a**2 + h**2) / 12.
    Izz = mass * (a**2 + b**2) / 12.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def circular_cone_shell(R, height, density=densities.get_density('STEEL'), thickness=0.02):
    """Get the inertia of a circular cone shell

    Returns
    -------
    RigidBodyInertia
    
    Note
    ----
    The center of gravity is located at an altitude of z=H/3 over the circular basis center
    """
    
    surface = pi * R * sqrt(R**2 + height**2)
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = mass * (R**2 + 2*height**2/9.) / 4.
    Izz = mass * R**2 / 2.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def frustrum_of_circular_cone_shell(r, R, height, density=densities.get_density('STEEL'), thickness=0.02):
    """Get the inertia of a frustrum of circular cone shell

    Returns
    -------
    RigidBodyInertia
    
    Note
    ----
    The center of gravity is located at an altitude of z=(H/3)*(2*r+R)/(r+R)
    """
    
    surface = pi * (R + r) * sqrt(height**2 + (R - r)**2)
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = mass * (R**2+r**2)/4. + mass * height**2 * (1+2*R*r/(R+r)**2) / 18.
    Izz = mass * (R**2+r**2) / 2.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def lateral_cylindrical_shell(R, H, density=densities.get_density('STEEL'), thickness=0.02):
    """Get the inertia of a lateral cylindrical shell

    Returns
    -------
    RigidBodyInertia
    """
    
    surface = 2*pi*R*H
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = mass * (R**2 + H**2/6.) / 2.
    Izz = mass * R**2
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def total_cylindrical_shell(R, H, density=densities.get_density('STEEL'), thickness=0.02):
    """Get the inertia of a total cylindrical shell

    Returns
    -------
    RigidBodyInertia
    """
    
    surface = 2 * pi * R * (R + H)
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = mass * (2*R**2*(R+2*H) + H**2*(3*R+H)) / 12. / (R+H)
    Izz = mass * R**2 * ((R+2*H) / (R+H)) / 2.
    
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def spherical_shell(R, density=densities.get_density('STEEL') ,thickness=0.02):
    """Get the inertia of a spherical shell

    Returns
    -------
    RigidBodyInertia
    """
    
    surface = 4*pi*R**2
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = Izz = 2 * mass*R**2 / 3.
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def hemispherical_shell(R, density=densities.get_density('STEEL'), thickness=0.02):
    """Get the inertia of a hemispherical shell

    Returns
    -------
    RigidBodyInertia
    
    Note
    ----
    The center of gravity is located at an altitude of z=R/2 over the circular basis center
    """
    
    surface = 2 * pi * R ** 2
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = 5 * mass * R ** 2 / 12.
    Izz = 2 * mass * R ** 2 / 3.
    return RigidBodyInertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


# TODO: placer cette classe tout en haut du module
# TODO: faire une methode pour templater a partir d'un array preexistant
# TODO: voir comment on gere les signes des produits d'inertie


class RotationalInertia3D(np.ndarray):
    
    __array_priority__ = 15
    
    def __new__(cls, xx, xy, yy, xz, yz, zz, point):
        
        lower_triangular = np.array([xx, xy, yy, xz, yz, zz], dtype=float)
        obj = lower_triangular.view(cls)
        
        return obj
    
    @property
    def array(self):
        print("generating full array")
        array = np.asarray(self, dtype=float)[[0, 1, 3, 1, 2, 4, 3, 4, 5]]
        array[[1, 2, 3, 5, 6, 7]] *= -1
        return array.reshape((3, 3))
    
    # def __mul__(self, other):
    #     print '\nOperation %s * %s' % (repr(self), repr(other))
    #
    #     if np.isscalar(other):
    #         return np.multiply(self, other)
    #
    #     elif isinstance(other, AngularVelocityVector):
    #         return np.dot(self.array, other)
    #
    #     else:
    #         raise RuntimeError('Operation not allowed: %s * %s') % (type(self), type(other))
    #
    #
    #     # return np.multiply(self, other)
    #
    # def __rmul__(self, other):
    #     print '\nRMUL %s * %s' % (repr(other), repr(self))
    # #
    # #     if isinstance(other, (float, int)):
    # #         return np.multiply(other, self)
    # #     elif isinstance(other, RotationalInertia3D):
    # #         print "I*I"
    # #         return
    # #     else:
    # #         return NotImplemented

    def __numpy_ufunc__(self, ufunc, method, i, inputs, **kwargs):
        print("In __numpy_ufunc__")
        return NotImplemented
    

    # def __array_finalize(self, obj):
    #     print "In __array_finalize__ :"
    #     print '\tself is %s' % repr(self)
    #     print '\tobj is %s' % repr(obj)
    #     if obj is None: return
    
    # def __array_prepare__(self, array, context=None):
    #     print "In __array_prepare__:"
    #     (func, args, _) = context
    #     print func
    #     print args


    # def __array_wrap__(self, obj, context=None):
    #     print 'In __array_wrap__ :'
    #     print '\tself is %s' % repr(self)
    #     print '\tobj: %s' % repr(obj)
    #     print '\tcontext: ', context
    #
    #     (func, args, _) = context
    #     print func, args
    #
    #
    #     return np.ndarray.__array_wrap__(self, obj, context)
    
    # def __array__(self):
    #     print "In __array__"
    
    # def __getitem__(self, key):
    #     print 'In __getitem__ with key :'
    #     print key
    #     return super(RotationalInertia3D, self).__getitem__(key)


class AngularVelocityVector(np.ndarray):
    __array_priority__ = 15
    def __new__(cls, array):
        assert len(array) == 3
        return np.asarray(array, dtype=float).view(cls)


if __name__ == '__main__':
    
    # inertia = RigidBodyInertia(1, [0, 0, 0], 1, 2, 3, 4, 5, 6)
    
    inertia = RotationalInertia3D(1, 2, 3, 4, 5, 6, [0, 0, 0])
    
    w = AngularVelocityVector([1, 1, 1])
    
    # print inertia.array
    # print inertia * 2
    # print inertia.__array_priority__
    # print w.__array_priority__
    print((inertia * 2))
    # print w*inertia
    # print inertia*2

    # print type(inertia.array)
    # print inertia.array
    # print inertia.array is inertia
    #
    # print inertia[:3]

    # inertia = RigidBodyInertia(1, [0, 0, 0], 1, 2, 3, 4, 5, 6, point=[1, 2, 3])
    
    # print inertia.inertia_matrix
    # print inertia.reduction_point
    # print inertia.at_cog.inertia_matrix
    #
    # inertia.reduction_point = [2, 2, 2]
    # print inertia.inertia_matrix
    # print inertia.at_cog.inertia_matrix
    #
    # inertia.shift_at_cog()
    # print inertia.inertia_matrix
    #
    # print inertia
    
    # R = 5
    # r = 4.5
    # R1 = 4.75
    # h = 20
    # rho = 8000
    #
    # e = 0.0001
    # r_point = [4, 3, -h/2]
    #
    # rcyl = hollow_right_circular_cylinder(r, R, h, density=rho)
    # rcyl.reduction_point = r_point
    #
    # lcyl = lateral_cylindrical_shell(R1, h, rho*(R**2-r**2) / (2*R1*e), e)
    # lcyl.reduction_point = r_point
    #
    # print rcyl
    # print lcyl
