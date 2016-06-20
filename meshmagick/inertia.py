#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import mesh
import nodes
import frame
import rigid_body as rb

import numpy as np

# FIXME: attention, changer les signes pour les produits d'inertie !!!!!!!!!

# TODO: ajouter la production d'inerties de solides connus --> utile pour comparaison !!

# import sys
class RigidBodyInertia(np.ndarray):
    def __new__(cls, name, Ixx=0., Iyy=0., Izz=0., Ixy=0., Ixz=0., Iyz=0., node=[0., 0., 0.], parent=None):
        mat = np.array([[ Ixx, -Ixy, -Ixz],
                        [-Ixy,  Iyy, -Iyz],
                        [-Ixz, -Iyz,  Izz]], dtype=np.float)
        obj = mat.view(cls)
        obj._name = str(name)

        node_name = '_'.join((obj._name, 'reference_node'))
        obj._node = nodes.Node(node_name, node)

        if parent:
            if not isinstance(parent, rb.RigidBody):
                raise TypeError, "parent of a RigidBodyInertia must be a RigidBody object. Got %s instead" % type(
                    parent)
        obj._parent = parent
        return obj

    @property
    def Ixx(self):
        return self[0, 0]

    @property
    def Iyy(self):
        return self[1, 1]

    @property
    def Izz(self):
        return self[2, 2]

    @property
    def Ixy(self):
        return -self[0, 1]

    @property
    def Ixz(self):
        return -self[0, 2]

    @property
    def Iyz(self):
        return -self[1, 2]

    @Ixx.setter
    def Ixx(self, value):
        self[0, 0] = float(value)
        return

    @Iyy.setter
    def Iyy(self, value):
        self[1, 1] = float(value)
        return

    @Izz.setter
    def Izz(self, value):
        self[2, 2] = float(value)
        return

    @Ixy.setter
    def Ixy(self, value):
        self[0, 1] = self[1, 0] = -float(value)
        return

    @Ixz.setter
    def Ixz(self, value):
        self[0, 2] = self[2, 0] = -float(value)
        return

    @Iyz.setter
    def Iyz(self, value):
        self[1, 2] = self[2, 1] = -float(value)
        return

    # Some alternative constructors
    @classmethod
    def from_diag(cls, name, diag, node=[0., 0., 0.]):
        return cls.__new__(cls, name, Ixx=diag[0], Iyy=diag[1], Izz=diag[2], node=node)

    @classmethod
    def from_point_mass(cls, name, position, mass, node=[0., 0., 0.]):
        position = np.asarray(position, dtype=np.float) - np.asarray(node, dtype=np.float)
        x2, y2, z2 = position * position
        return cls.__new__(cls, name, Ixx=mass*(y2+z2), Iyy=mass*(x2+z2), Izz=mass*(x2+y2))

    # Class methods
    def shift_to_cog(self, cog, mass):
        """
        Assume that the current inertia is about the F frame's origin node, and expressed in F.

        Given the vector from the reference node to the body center of mass cog, and the mass of the body,
        we can shift the inertia to the center of mass. This produces a new inertia whose (implicit) frame F' is
        aligned with F but has origin at the cog and becmos a central inertial.

        Parameters
        ----------
        cog
        mass

        Returns
        -------

        """
        self -= mass * (np.inner(cog, cog)*np.eye(3) - np.outer(cog, cog))
        return

    def get_principal_axes(self):

        eig_values, eig_vectors = np.linalg.eigh(np.asarray(self))
        # u0, u1, u2 = map(np.asarray, eig_vectors.T)
        e0, e1, e2 = eig_vectors.T
        # print eig_values
        # print e0, e1, e2

        return eig_vectors.T



    def __str__(self):
        if self._parent:
            str_repr = "Rotational Inertia matrix of rigid body : '%s'" % self._parent.name
        else:
            str_repr = "Rotational Inertia matrix of a RigidBody"
        str_repr = '\n'.join((str_repr, "\t%10.3E, %10.3E, %10.3E" % (self.Ixx, self.Ixy, self.Ixz)))
        str_repr = '\n'.join((str_repr, "\t%10.3E, %10.3E, %10.3E" % (self.Ixy, self.Iyy, self.Iyz)))
        str_repr = '\n'.join((str_repr, "\t%10.3E, %10.3E, %10.3E\n" % (self.Ixz, self.Iyz, self.Izz)))
        str_repr = '\n'.join((str_repr, "Node '%s':" % (self._node.name))) # TODO: voir si il ne faut pas plutot definir
        # une frame ici !!!
        str_repr = '\n'.join((str_repr, "\t%10.3f, %10.3f, %10.3f" % (self._node[0], self._node[1], self._node[2])))
        return str_repr





    # @property
    # def rotational_inertia_matrix(self):
    #     return self[3:, 3:]

    # @property
    # def Ixx(self):
    #     return self.__internals__['inertia_matrix']







if __name__ == '__main__':

    # import mmio
    #
    # vertices, faces = mmio.load_VTP('SEAREV.vtp')
    # mymesh = Mesh(vertices, faces)

    # inertia = Inertia(mymesh)
    inertia = RigidBodyInertia.from_diag('essai', [1, 2, 3], [4, 5, 6])
    print inertia


    point = [1000, 20, 4]
    mass = 50

    inertia_point = RigidBodyInertia.from_point_mass('P0', point, mass)
    print inertia_point
