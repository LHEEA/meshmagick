#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from mesh import *

# TODO: Make this class abstract
class RigidBody(object):

    def __init__(self, mass=None, cog=None, body_mesh=None):
        self._mass = mass*1e3 # On entre la masse en tonnes !!!
        self._cog = cog
        self._body_mesh = body_mesh

class PlainRigidBody(RigidBody):
    def __init__(self, mass=None, cog=None, body_mesh=None, rho_material=1026.):
        RigidBody.__init__(self, mass=mass, cog=cog, body_mesh=body_mesh)
        self._rho_material = rho_material

    @property
    def mass(self):
        return self._mass*1e-3 # Attention: travailler en tonnes !!!

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

    @property
    def cog(self):
        return self._cog

    def _set_cog_from_properties(self):

        volume = self._body_mesh.volume
        mesh_normals = self._body_mesh.faces_normals

        sigma_6_8 = self._body_mesh.get_surface_integrals()[6:9]
        self._cog = (mesh_normals.T * sigma_6_8).sum(axis=1) / (2*volume)
        return

    def guess_mass_properties(self):
        self._set_cog_from_properties()
        self._set_mass_from_properties()
        return

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

    body.guess_mass_properties()
    print body.cog, body.mass
