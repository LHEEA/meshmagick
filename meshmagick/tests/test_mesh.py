#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np

from meshmagick.mmio import load_VTP
from meshmagick.mesh import Mesh, Plane
from math import pi

from os.path import curdir

vertices, faces = load_VTP('meshmagick/tests/data/Cylinder.vtp')
cylinder = Mesh(vertices, faces)
print(cylinder)

def test_plane():
    plane = Plane()
    plane.set_normal_from_angles(0., 0.)
    plane.get_normal_orientation_wrt_z()
    
    plane.flip_normal()
    plane.coord_in_plane(cylinder.vertices)
    try:
        plane.get_edge_intersection([-1, -1, -1], [-2, -2, -2])
    except RuntimeError:
        pass
    plane.set_plane_parameters(1, 1, 1)
    plane.get_normal_orientation_wrt_z()
    plane.orthogonal_projection_on_plane(cylinder.vertices)
    plane.get_origin()
    return

def test_verbose():
    verbose = cylinder.verbose
    cylinder.verbose = verbose
    cylinder.verbose_off()
    cylinder.verbose_on()
    return

def test_heal_mesh():
    cylinder.heal_mesh()
    return

def test_faces():
    faces_tmp = cylinder.faces
    # assert np.all(faces == faces_tmp)
    cylinder.faces = faces_tmp
    cylinder.id
    return
    
def test_vertices():
    vertices_tmp = cylinder.vertices
    # assert np.all(vertices_tmp == vertices)
    cylinder.vertices = vertices_tmp
    return

def test_remove_internals():
    cylinder.triangles_ids
    cylinder._remove_faces_properties()
    cylinder._remove_triangles_quadrangles()
    cylinder._remove_connectivity()
    cylinder._remove_surface_integrals()
    return

def test_quality():
    cylinder.print_quality()
    return

# def test_show():
#     cylinder.show()

def test_mesh_inertia():
    cylinder.eval_plain_mesh_inertias()
    cylinder.eval_shell_mesh_inertias()
    cylinder._remove_surface_integrals()
    return

def test_bbox():
    assert (-5, 5, -5, 5, -10, 10) == cylinder.axis_aligned_bbox
    assert (-10, 10, -10, 10, -10, 10) == cylinder.squared_axis_aligned_bbox
    return

def test_is_mesh_conformal():
    cylinder.is_mesh_conformal()
    return
    
def test_is_mesh_closed():
    cylinder.is_mesh_closed()
    return
    
def test_rotate():
    cylinder.rotate_x(pi)
    cylinder.rotate_y(pi)
    cylinder.rotate_z(pi)
    return
    
def test_translate():
    cylinder.translate_x(0)
    cylinder.translate_y(0)
    cylinder.translate_z(0)
    cylinder.translate([0, 0, 0])
    return
    
def test_scale():
    cylinder.scalex(1)
    cylinder.scaley(1)
    cylinder.scalez(1)
    cylinder.scale(1)
    return

def test_merge_duplicate():
    cylinder.verbose_on()
    cylinder.merge_duplicates(atol=1e-5, return_index=True)
    return
    
def test_triangulate_quadrangles():
    cylinder.triangulate_quadrangles()
    return
    
def test_symmetrize():
    cylinder.symmetrize(Plane())
    return

def test_mirror():
    cylinder.mirror(Plane())
    
def test_quick_save():
    import os
    cylinder.quick_save()
    os.remove('quick_save.vtp')
    
