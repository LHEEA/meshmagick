#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np
from math import pi, sqrt
from copy import deepcopy


# FIXME: attention, changer les signes pour les produits d'inertie !

# TODO: ajouter la production d'inerties de solides connus --> utile pour comparaison !!


densities = {'CONCRETE': 2300.,
             'REINFORCED_CONCRETE': 2400.,
             'FRESH_WATER': 1000.,
             'SALT_WATER': 1025.,
             'SAND': 1600.,
             'STEEL': 7850.,
             'ALUMINUM': 2700.,
             'LEAD': 11350.,
             'TITANIUM': 4500.,
             'POLYSTYRENE': 1050.,
             'GASOLINE': 750.,
             'DIESEL_FUEL': 850.,
             'ETHANOL': 789.,
             'AIR_20DEGC': 1.204,
             'BUTANE': 2.7,
             'PROPANE': 2.01,
             'HYDROGEN_-252DEGC': 70.,
             'NITROGEN_-195DEGC': 810.
             }


# TODO: indiquer une frame d'expression ...

class Inertia(object):
    def __init__(self, mass, cog, xx, yy, zz, yz, xz, xy, point=None):
        
        self._mass = float(mass)
        
        assert len(cog) == 3
        self._cog = np.asarray(cog, dtype=np.float)
        
        if point is None:
            self._point = self._cog
        else:
            assert len(point) == 3
            self._point = np.asarray(point, dtype=np.float)
        
        self._3d_rotational_inertia = np.array([[float(xx), -float(xy), -float(xz)],
                                                [-float(xy), float(yy), -float(yz)],
                                                [-float(xz), -float(yz), float(zz)]], dtype=np.float)
    
    @property
    def mass(self):
        return self._mass
    
    @property
    def gravity_center(self):
        return self._cog
    
    @property
    def inertia_matrix(self):
        return self._3d_rotational_inertia
    
    @property
    def reduction_point(self):
        return self._point
    
    @reduction_point.setter
    def reduction_point(self, point):
        """Change the reduction point"""
        mat_at_cog = self.at_cog.inertia_matrix
        assert len(point) == 3
        self._point = np.asarray(point, dtype=np.float)
        self._3d_rotational_inertia = mat_at_cog + self._huygens_transport()
        return
    
    @property
    def at_cog(self):
        inertia = deepcopy(self)
        inertia.shift_at_cog()
        return inertia
    
    def shift_at_cog(self):
        self._3d_rotational_inertia -= self._huygens_transport()
        self._point = self._cog
        return
    
    def is_at_cog(self):
        return np.all(self._point == self._cog)
    
    def _huygens_transport(self):
        p_g = self._cog - self._point
        return self._mass * (np.dot(p_g, p_g) * np.eye(3) - np.outer(p_g, p_g))
    
    @property
    def xx(self):
        return self._3d_rotational_inertia[0, 0]
    
    @property
    def yy(self):
        return self._3d_rotational_inertia[1, 1]
    
    @property
    def zz(self):
        return self._3d_rotational_inertia[2, 2]
    
    @property
    def yz(self):
        return -self._3d_rotational_inertia[1, 2]
    
    @property
    def xz(self):
        return -self._3d_rotational_inertia[0, 2]
    
    @property
    def xy(self):
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
    vol = pi * radius ** 2 * length
    mass = density * vol
    
    Ixx = Iyy = mass * (3 * radius ** 2 + length ** 2) / 12.
    Izz = mass * radius ** 2 / 2.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def hollow_right_circular_cylinder(int_radius, ext_radius, length, density=1.):
    vol = pi * length * (ext_radius ** 2 - int_radius ** 2)
    mass = density * vol
    
    R2r2 = ext_radius ** 2 + int_radius ** 2
    Ixx = Iyy = mass * (3 * R2r2 + length ** 2) / 12.
    Izz = mass * R2r2 / 2.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def right_circular_cone(radius, length, density=1.):
    """Note:
    The center of gravity is at an altitude of z = H/4 over the circular basis center
    """
    vol = pi * radius ** 2 * length / 3.
    mass = density * vol
    
    Ixx = Iyy = 3 * mass * (radius ** 2 + length ** 2/4.) / 20.
    Izz = 3. * mass * radius ** 2 / 10.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def sphere(radius, density=1.):
    vol = 4. * pi * radius ** 3 / 3.
    mass = density * vol
    
    Ixx = Iyy = Izz = 2 * mass * radius ** 2 / 5.
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0., 0., 0.)


def hollow_sphere(int_radius, ext_radius, density=1.):
    vol = 4. * pi * (ext_radius ** 3 - int_radius ** 3) / 3.
    mass = density * vol
    
    Ixx = Iyy = Izz = 2 * mass * (ext_radius ** 5 - int_radius ** 5) / (ext_radius ** 3 - int_radius ** 3) / 5
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0., 0., 0.)


def hemisphere(radius, density=1.):
    """Note
    The center of gravity is situated at the altitude of z = 3R/8 over the circular basis center
    """
    vol = 2 * pi * radius ** 3 / 3.
    mass = density * vol
    
    Ixx = Iyy = Izz = 0.26 * mass * radius ** 2
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def elliptical_cylinder(a, b, length, density=1.):
    """Note :
    The center of gravity is located at an altitude of z=H/2 over the elliptical basis center
    a is along x axis (ellipse semi axis)
    b is along y axis (ellipse semi axis)
    length is along z axis
    """
    
    vol = pi * a * b * length
    mass = density * vol
    
    Ixx = mass * (3 * b ** 2 + length ** 2) / 12.
    Iyy = mass * (3 * a ** 2 + length ** 2) / 12.
    Izz = mass * (a ** 2 + b ** 2) / 4.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def ellipsoid(a, b, c, density=1.):
    """Note
    a is along z axis (ellipse semi axis)
    b is along x axis (ellipse semi axis)
    c is along y axis (ellipse semi axis)
    """
    vol = 4 * pi * a * b * c / 3.
    mass = density * vol
    
    Ixx = mass * (a ** 2 + c ** 2) / 5.
    Iyy = mass * (a ** 2 + b ** 2) / 5.
    Izz = mass * (b ** 2 + c ** 2) / 5.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def torus(chord_radius, tube_radius, density=1.):
    vol = 2*pi**2 * tube_radius**2 * chord_radius
    mass = density * vol
    
    Ixx = Izz = mass * (4*chord_radius**2 + 5*tube_radius**2) / 8.
    Iyy = mass * (4*chord_radius**2 + 3*tube_radius**2) / 4.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)


def right_angle_wedge(base, height, length, density=1.):
    vol = base*height*length / 2.
    mass = density * vol
    
    Ixx = mass * (2*height**2 + 3*length**2) / 36.
    Iyy = mass * (base**2 + height**2) / 18.
    Izz = mass * (2*base**2 + 3*length**2) / 36.

    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def isoceles_wedge(base, height, length, density=1.):
    vol = base * height * length / 2.
    mass = density * vol

    Ixx = mass * (2 * height**2 + 3 * length**2) / 36.
    Iyy = mass * (4*height**2 + 3*base**2) / 72.
    Izz = mass * (2 * length**2 + base**2) / 24.

    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def right_rectangular_pyramid(a, b, height, density=1.):
    """Note
    The center of gravity is located at the altitude z=H/4 over the rectangular basis center
    """
    vol = a*b*height / 3.
    mass = density * vol
    
    Ixx = mass * (b**2 + 3*height**2/4.) / 20.
    Iyy = mass * (a**2 + 3*height**2/4.) / 20.
    Izz = mass * (a**2 + b**2) / 20.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def cube(a, density=1.):
    vol = a**3
    mass = density * vol
    
    Ixx = Iyy = Izz = mass * a**2 / 6.
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def rectangular_prism(a, b, h, density=1.):
    """Note:
    a is along x
    b is along y
    h is along z"""
    vol = a*b*h
    mass = density * vol
    
    Ixx = mass * (b**2 + h**2) / 12.
    Iyy = mass * (a**2 + h**2) / 12.
    Izz = mass * (a**2 + b**2) / 12.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def circular_cone_shell(R, height, density=densities['STEEL'], thickness=0.02):
    """Note
    The center of gravity is located at an altitude of z=H/3 over the circular basis center
    """
    surface = pi * R * sqrt(R**2 + height**2)
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = mass * (R**2 + 2*height**2/9.) / 4.
    Izz = mass * R**2 / 2.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

    

def frustrum_of_circular_cone_shell(r, R, height, density=densities['STEEL'], thickness=0.02):
    """Note:
    The center of gravity is located at an altitude of z=(H/3)*(2*r+R)/(r+R)
    """
    surface = pi * (R + r) * sqrt(height**2 + (R - r)**2)
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = mass * (R**2+r**2)/4. + mass * height**2 * (1+2*R*r/(R+r)**2) / 18.
    Izz = mass * (R**2+r**2) / 2.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def lateral_cylindrical_shell(R, H, density=densities['STEEL'], thickness=0.02):
    surface = 2*pi*R*H
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = mass * (R**2 + H**2/6.) / 2.
    Izz = mass * R**2
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def total_cylindrical_shell(R, H, density=densities['STEEL'], thickness=0.02):
    surface = 2 * pi * R * (R + H)
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = mass * (2*R**2*(R+2*H) + H**2*(3*R+H)) / 12. / (R+H)
    Izz = mass * R**2 * ((R+2*H) / (R+H)) / 2.
    
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def spherical_shell(R, density=densities['STEEL'] ,thickness=0.02):
    surface = 4*pi*R**2
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = Izz = 2 * mass*R**2 / 3.
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)

def hemispherical_shell(R, density=densities['STEEL'], thickness=0.02):
    """Note
    The center of gravity is located at an altitude of z=R/2 over the circular basis center
    """
    surface = 2 * pi * R ** 2
    sigma = density * thickness
    mass = sigma * surface
    
    Ixx = Iyy = 5 * mass * R ** 2 / 12.
    Izz = 2 * mass * R ** 2 / 3.
    return Inertia(mass, [0, 0, 0], Ixx, Iyy, Izz, 0, 0, 0)




if __name__ == '__main__':
    inertia = Inertia(1, [0, 0, 0], 1, 2, 3, 4, 5, 6, point=[1, 2, 3])
    
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
    
    R = 5
    r = 4.5
    R1 = 4.75
    h = 20
    rho = 8000
    
    e = 0.0001
    r_point = [4, 3, -h/2]
    
    rcyl = hollow_right_circular_cylinder(r, R, h, density=rho)
    rcyl.reduction_point = r_point
    
    lcyl = lateral_cylindrical_shell(R1, h, rho*(R**2-r**2) / (2*R1*e), e)
    lcyl.reduction_point = r_point
    
    print rcyl
    print lcyl
