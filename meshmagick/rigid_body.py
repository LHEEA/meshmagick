#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import mesh
import nodes
import inertia
import frame
from itertools import count
from warnings import warn

# TODO: Make this class abstract
# TODO: make a system of frame attachment (observer pattern)
class RigidBody(object):
    _ids = count(0)
    def __init__(self, name=None, mass=None, cog=None, rb_inertia=None):

        self._id = self._ids.next()

        if name:
            self._name = str(name)
        else:  # Default
            self._name = 'RigidBody_%u' % self._id

        self.mass = mass
        self.cog = cog

        if rb_inertia:
            if not isinstance(rb_inertia, RigidBodyInertia): # TODO:
                raise ValueError, "RigidBody rb_inertia parameter must be a RigidBodyInertia object"

        # Setting the reference frame
        self._reference_frame = frame.Frame('_'.join((self._name, 'reference_frame')))



        self._body_mesh = None



        # # TODO: Passer a une classe
        # self._spatial_inertia_tensor = None

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
            if not isinstance(value, mesh.Mesh):
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
            self._cog = nodes.Node('_'.join((self._name, 'COG')), value)
        else:
            self._cog = None
        return

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = str(value)

    def show(self):
        if not self.body_mesh:
            warn("Unable to visualize %s as no mesh has been specified" % self._name)

        # TODO : FINIR
        # TODO: devrait montrer le cog, la masse, les axes principaux d'inertie...
        return


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

    def guess_mass_properties_from_mesh(self):
        # Following formula from technical report on inertia computations

        if not self.body_mesh:
            warn("Unable to compute %s's mass properties from mesh as no mesh has been specified" % self._name)
            return

        # Computing mass
        # --------------
        mesh_volume = self._body_mesh.volume
        mass = self._rho_material * mesh_volume
        self._mass = mass
        # TODO: voir si on ne supprime pas la ppte en tant que tel --> recuperer depuis l'objet Spatial Inertia

        # Computing COG
        # -------------
        mesh_normals = self._body_mesh.faces_normals.T
        integrals = self._body_mesh.get_surface_integrals()[6:15]
        sigma_6_8 = integrals[:3]

        cog = (mesh_normals * sigma_6_8).sum(axis=1) / (2*mesh_volume)
        # cog[np.fabs(cog) < 1e-6] = 0. # TODO: voir si utile, ou mettre en option d'affichage --> le mieux ?
        self._cog = cog

        # Computing the rotational inertia matrix
        # ---------------------------------------
        sigma9, sigma10, sigma11 = (mesh_normals * integrals[3:6]).sum(axis=1)
        sigma12, sigma13, sigma14 = (mesh_normals * integrals[6:10]).sum(axis=1)

        rho = self._rho_material

        Ixx = rho * (sigma10+sigma11) / 3.
        Iyy = rho * (sigma9+sigma11) / 3.
        Izz = rho * (sigma9+sigma10) / 3.
        Ixy = rho * sigma12 / 2.
        Ixz = rho * sigma14 / 2.
        Iyz = rho * sigma13 / 2.

        inertia_name = '_'.join((self._name, 'rotational_inertia'))
        self._inertia = inertia.RigidBodyInertia(inertia_name, Ixx=Ixx, Iyy=Iyy, Izz=Izz, Ixy=Ixy, Ixz=Ixz, Iyz=Iyz,
                                         parent=self)
        # print self._inertia

        return



class ShellRigidBody(RigidBody):

    def __init__(self, mass=None, cog=None, body_mesh=None, thickness=4., rho_material=7500.):
        # Attention thickness doit etre donne en cm !!
        RigidBody.__init__(self, mass=mass, cog=cog, body_mesh=body_mesh)
        self._thickness = thickness
        pass



if __name__ == '__main__':
    import mmio
    vertices, faces = mmio.load_VTP('Cylinder.vtp')

    mymesh = mesh.Mesh(vertices, faces)
    # mymesh.show()

    # Pour le cylindre
    plane = mesh.Plane()
    mymesh.symmetrize(plane)
    plane.normal = [0, 1, 0]
    mymesh.symmetrize(plane)

    # mymesh.translate([2000., 100, 3.2])

    mymesh.triangulate_quadrangles()
    # print mymesh
    # mymesh.show()
    # Fin pour cylindre

    body = PlainRigidBody(name='mybody', mass=0., cog=[0, 0, 0])

    body.body_mesh = mymesh
    body.name = "Mybody"

    body.guess_mass_properties_from_mesh()




