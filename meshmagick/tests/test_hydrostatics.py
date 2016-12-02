#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import meshmagick.mmio as mmio
import meshmagick.hydrostatics as hs
from meshmagick.mesh import Mesh
from math import pi, fabs


# Importing cylinder mesh
vertices, faces = mmio.load_VTP('meshmagick/tests/data/Cylinder.vtp')
cylinder = Mesh(vertices, faces)
# cylinder.show()

hs_cylinder = hs.Hydrostatics(cylinder)

vertices, faces = mmio.load_VTP('meshmagick/tests/data/SEAREV.vtp')
searev = Mesh(vertices, faces)

def test_verbose():
    hs_cylinder.verbose
    hs_cylinder.verbose_off()
    hs_cylinder.verbose_off()
    return

def test_accessors():
    hs_cylinder.gravity
    hs_cylinder.gravity = 9.81
    hs_cylinder.rho_water = 1025.
    mass = hs_cylinder.mass
    hs_cylinder.mass = mass
    cg = hs_cylinder.gravity_center
    hs_cylinder.gravity_center = cg
    zg = hs_cylinder.zg
    hs_cylinder.zg = zg
    hs_cylinder.wet_surface_area
    hs_cylinder.displacement_volume
    hs_cylinder.displacement
    hs_cylinder.buoyancy_center
    hs_cylinder.flotation_surface_area
    hs_cylinder.flotation_center
    hs_cylinder.transversal_metacentric_height
    hs_cylinder.longitudinal_metacentric_height
    hs_cylinder.transversal_metacentric_radius
    hs_cylinder.longitudinal_metacentric_radius
    hs_cylinder.hydrostatic_stiffness_matrix
    hs_cylinder.S33
    hs_cylinder.S34
    hs_cylinder.S35
    hs_cylinder.S44
    hs_cylinder.S45
    hs_cylinder.S55
    reltol = hs_cylinder.reltol
    hs_cylinder.reltol = reltol
    tr = hs_cylinder.theta_relax
    hs_cylinder.theta_relax = tr
    zr = hs_cylinder.z_relax
    hs_cylinder.z_relax = zr
    max_iter = hs_cylinder.max_iterations
    hs_cylinder.max_iterations = max_iter
    max_restart = hs_cylinder.max_restart
    hs_cylinder.max_restart = max_restart
    allow_unstable = hs_cylinder.allow_unstable
    hs_cylinder.allow_unstable = allow_unstable
    hs_cylinder.allow_unstable_on()
    hs_cylinder.allow_unstable_off()
    hs_cylinder.reset()
    return

def test_set_displacement():
    disp = hs_cylinder.displacement
    hs_cylinder.set_displacement(1.1*disp)
    hs_cylinder.reset()
    # hs_cylinder.show()
    return

def test_equilibrate():
    hs_cylinder.gravity_center = [0, 0, 3]
    hs_cylinder.equilibrate()
    hs_cylinder.show()
    hs_cylinder.reset()
    # hs_cylinder.show()
    return

def test_hydrostatic_report():
    hs_cylinder.get_hydrostatic_report()
    return
    

def test_cylinder_base_hydrostatics():
    
    hydrostatics = hs.Hydrostatics(cylinder)
    hs_data = hydrostatics.hs_data
    
    # Analytical results
    d = 10.
    h = 10.
    waterplane_area = pi*d**2 / 4.
    immersed_volume = waterplane_area*h
    wet_surface_area = waterplane_area + pi*d*h

    assert fabs(hs_data['waterplane_area'] - waterplane_area) < 1
    assert fabs(hs_data['wet_surface_area'] - wet_surface_area) < 2
    assert fabs(hs_data['disp_volume'] - immersed_volume) < 10


def test_searev_base_hydrostatics():
    hydrostatics = hs.Hydrostatics(searev)
    hs_data = hydrostatics.hs_data

    assert fabs(hs_data['waterplane_area'] - 300) < 1
    assert fabs(hs_data['wet_surface_area'] - 550) < 1
    assert fabs(hs_data['disp_volume'] - 1177) < 1




# def test_displacement_equilibrium():
#
#     disp = 2000
#     CG = [-1, 0, -3]
#
#     cylinder_eq, CGc = hs.compute_equilibrium(searev, disp, CG, rho_water=1023, anim=True)
    
