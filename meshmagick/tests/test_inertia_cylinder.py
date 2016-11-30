#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from meshmagick.mmio import load_VTP
from meshmagick.mesh import Mesh
from meshmagick.inertia import right_circular_cylinder
from meshmagick.densities import get_density

# Comparaison pour un calcul au cog

vertices, faces = load_VTP('meshmagick/tests/data/Cylinder.vtp')
cyl_mesh = Mesh(vertices, faces)
xmin, xmax, ymin, ymax, zmin, zmax = cyl_mesh.axis_aligned_bbox


# Analytique
R = (xmax-xmin) / 2.
H = zmax-zmin
cyl_anal_inertia = right_circular_cylinder(R, H, density=get_density('SALT_WATER'))

# print cyl_mesh_inertia
# print cyl_anal_inertia.at_cog


# On decale
cyl_mesh.translate([1, 3, -10])
# cyl_mesh.show()

print "Numerique:"
print cyl_mesh.eval_plain_mesh_inertias(get_density('SALT_WATER'))

cyl_anal_inertia.reduction_point = [-1, -3, 10]
print "analytique:"
print cyl_anal_inertia

