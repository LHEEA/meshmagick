#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from mmio import load_VTP
from mesh import Mesh
from inertia_parameters import densities, right_circular_cylinder

# Comparaison pour un calcul au cog

vertices, faces = load_VTP('meshmagick/tests/data/Cylinder_fine_mesh.vtp')
cyl_mesh = Mesh(vertices, faces)
xmin, xmax, ymin, ymax, zmin, zmax = cyl_mesh.axis_aligned_bbox


# Analytique
R = (xmax-xmin) / 2.
H = zmax-zmin
cyl_anal_inertia = right_circular_cylinder(R, H, density=densities['SALT_WATER'])

# print cyl_mesh_inertia
# print cyl_anal_inertia.at_cog


# On decale
cyl_mesh.translate([1, 3, -10])
# cyl_mesh.show()

print "Numerique:"
print cyl_mesh.eval_plain_mesh_inertias(densities['SALT_WATER'])

cyl_anal_inertia.reduction_point = [-1, -3, 10]
print "analytique:"
print cyl_anal_inertia

