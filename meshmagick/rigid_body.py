#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from mesh import *

# TODO: Make this class abstract
class RigidBody(object):

    def __init__(self, mass=None, cog=None, body_mesh=None):
        self._mass = mass
        self._cog = cog
        self._body_mesh = body_mesh

class PlainRigidBody(RigidBody):
    def __init__(self, mass=None, cog=None, body_mesh=None):
        RigidBody.__init__(self, mass=mass, cog=cog, body_mesh=body_mesh)


    def get_cog_from_properties(self):

        volume = self._body_mesh.volume
        mesh_normals = self._body_mesh.faces_normals

        sigma_6_8 = self._body_mesh.get_surface_integrals()[6:9]

        return (mesh_normals.T * sigma_6_8).sum(axis=1) / (2*volume)



class ShellRigidBody(RigidBody):

    def __init__(self, mass=None, cog=None, body_mesh=None, thickness=4., rho_material=7500.):
        # Attention thickness doit etre donne en cm !!
        RigidBody.__init__(self, mass=mass, cog=cog, body_mesh=body_mesh)
        self._thickness = thickness
        pass

if __name__ == '__main__':
    import mmio
    vertices, faces = mmio.load_VTP('SEAREV.vtp')

    mymesh = Mesh(vertices, faces)

    body = PlainRigidBody(body_mesh=mymesh)

    print body.get_cog_from_properties()
