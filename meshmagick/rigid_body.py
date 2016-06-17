#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from mesh import *
from nodes import *
from inertia import *
from itertools import count

# TODO: Make this class abstract
class RigidBody(object):
    _ids = count(0)
    def __init__(self, name=None, mass=None, cog=None): # TODO: doit avoir name, cog, frame, mass, inertia
        """
        Rigid Body constructor

        Parameters
        ----------
        mass : float
            mass of the body (should be specified in tons but is stored in kg)
        cog : array_like
            position of the center of gravity in the body reference frame
        name : str
            The body's name
        """

        if name:
            self._name = str(name)
        else:  # Default
            self._name = 'RigidBody_%u' % self._id

        self.mass = mass
        self.cog = cog

        self._id = self._ids.next()

        self._body_mesh = None

        # TODO: Passer a une classe
        self._spatial_inertia_tensor = None

        return

    def __str__(self):
        # TODO : regulariser la maniere de construire la representation
        if self._mass is not None and self._cog is not None:

            str_repr = \
            """
            --------------------------------------------
            \tRigid Body NAME : %s
            --------------------------------------------

            Mass : %.3f (tons)
            Center of gravity : %.3f, %.3f, %.3f (m)

            """ % (self._name,
                   self.mass,
                   self._cog[0],
                   self._cog[1],
                   self._cog[2]
                   )

        else:
            str_repr = \
            """
            Rigid body : %s
            """ % self._name

        return str_repr

    @property
    def body_mesh(self):
        return self._body_mesh

    @body_mesh.setter
    def body_mesh(self, value):
        if value is not None:
            if not isinstance(value, Mesh):
                raise TypeError, "RigidBody Mesh must be a Mesh object"
        self._body_mesh = value
        return

    @property
    def mass(self):
        return self._mass*1e-3

    @mass.setter
    def mass(self, value):
        if value:
            self._mass = float(value) * 1e3
        else:
            self._mass = value
        return

    @property
    def cog(self): # TODO: faire un guess
        return self._cog

    @cog.setter
    def cog(self, value):
        if value:
            self._cog = Node(value, name='_'.join((self._name, 'COG')))
        else:
            self._cog = None
        return

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = str(value)

    # @property
    # def spatial_inertia_tensor(self):
    #     return self._spatial_inertia_tensor

    # @spatial_inertia_tensor.setter
    # def spatial_inertia_tensor(self, value):
    #     value = np.asarray(value, dtype=np.float)
    #     if value.shape != (6, 6):
    #         raise ValueError, 'Rigid Body spatial inertia tensor must have shape (6, 6)'
    #
    #     # self._spatial_inertia_tensor = InertiaTensor(value) # TODO : Passer a une classe
    #     self._spatial_inertia_tensor = value
    #     return







class PlainRigidBody(RigidBody):

    def __init__(self, mass=None, cog=None, name=None, rho_material=1026.):

        RigidBody.__init__(self, mass=mass, cog=cog, name=name)
        self._rho_material = float(rho_material)

    def __str__(self):
        # TODO : Regulariser la maniere de representer la chose... (comme pour RigidBody)
        str_repr = RigidBody.__str__(self)
        str_rho_material = \
        """
            rho material : %.3f (kg/m3)
        """ % self._rho_material

        str_repr = ''.join((str_repr, str_rho_material))
        return str_repr

    @property
    def rho_material(self):
        return self._rho_material

    @rho_material.setter
    def rho_material(self, value):
        self._rho_material = float(value)
        return

    def _set_mass_from_properties(self):
        self._mass = self._rho_material * self._body_mesh.volume
        return

    def _set_cog_from_properties(self):

        volume = self._body_mesh.volume
        mesh_normals = self._body_mesh.faces_normals

        sigma_6_8 = self._body_mesh.get_surface_integrals()[6:9]
        cog = (mesh_normals.T * sigma_6_8).sum(axis=1) / (2*volume)
        cog[np.fabs(cog) < 1e-6] = 0.
        self._cog = cog
        return

    def _set_inertia_from_properties(self):
        # Following formula from technical report on inertia computations

        mesh_normals = self._body_mesh.faces_normals.T
        integrals = self._body_mesh.get_surface_integrals()[9:15]

        sigma9, sigma10, sigma11 = (mesh_normals * integrals[:3]).sum(axis=1)
        sigma12, sigma13, sigma14 = (mesh_normals * integrals[3:]).sum(axis=1)

        # vol = self._body_mesh.volume

        Ixx = (sigma10+sigma11) / 3. # TODO: voir pourquoi besoin de diviser par le volume !!!
        Iyy = (sigma9+sigma11) / 3.
        Izz = (sigma9+sigma10) / 3.
        Ixy = -sigma12 / 2.
        Ixz = -sigma14 / 2.
        Iyz = -sigma13 / 2.

        rotational_inertia = self._rho_material * \
            np.asarray(
                [[Ixx, Ixy, Ixz],
                 [Ixy, Iyy, Iyz],
                 [Ixz, Iyz, Izz]],
        dtype=np.float)

        # rotational_inertia[np.fabs(rotational_inertia) < 1e-6] = 0.

        cog = self._cog
        mass = self._mass

        # TODO:

        # Expressing the inertia matrix at cog
        rotational_inertia -= mass * (np.inner(cog, cog)*np.eye(3) - np.outer(cog, cog))

        print rotational_inertia
        # TODO: A finir, il faut maintenant stocker la matrice dans une classe InertiaTensor



        return

    def guess_mass_properties(self):
        self._set_cog_from_properties()
        self._set_mass_from_properties()
        self._set_inertia_from_properties()
        return

    # TODO: redefinir setters et getter pour mass et cog...

    # @property
    # def inertia_matrix(self):
    #     return self._inertia_matrix



class ShellRigidBody(RigidBody):

    def __init__(self, mass=None, cog=None, body_mesh=None, thickness=4., rho_material=7500.):
        # Attention thickness doit etre donne en cm !!
        RigidBody.__init__(self, mass=mass, cog=cog, body_mesh=body_mesh)
        self._thickness = thickness
        pass



if __name__ == '__main__':
    import mmio
    vertices, faces = mmio.load_VTP('Cylinder.vtp')

    mymesh = Mesh(vertices, faces)

    # Pour le cylindre
    plane = Plane()
    mymesh.symmetrize(plane)
    plane.normal = [0, 1, 0]
    mymesh.symmetrize(plane)

    mymesh.translate([20., 100, 3.2])

    mymesh.triangulate_quadrangles()
    # print mymesh
    # mymesh.show()
    # Fin pour cylindre

    body = PlainRigidBody(name='mybody', mass=0., cog=[0, 0, 0])

    body.cog = [0, 0, 0]


    body.body_mesh = mymesh
    body.name = "Mybody"

    # print body.body_mesh.volume

    body.guess_mass_properties()

    # print body.mass, body.body_mesh.volume*1026*1e-3
    # print body.cog


    # print type(body.spatial_inertia_tensor)

    # print body

