#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""This module allows to perform hydrostatics computations on meshes"""

import numpy as np
import math

from .mesh_clipper import MeshClipper
from .rotations import cardan_to_rotmat, rotmat_to_cardan
from .htransform import HTransform

__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-2015, Ecole Centrale de Nantes"
__credits__ = "Francois Rongere"
__licence__ = "CeCILL"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@ec-nantes.fr"
__status__ = "Development"


def compute_hydrostatics(mesh, rho_water, grav, rotmat_corr=np.eye(3), z_corr=0., zcog=0.):

    wmesh = mesh.copy()
    wmesh.rotate_matrix(rotmat_corr)
    wmesh.translate_z(z_corr)

    clipper = MeshClipper(wmesh, assert_closed_boundaries=True, verbose=False)
    cmesh = clipper.clipped_mesh

    wet_surface_area = cmesh.faces_areas

    inertia = cmesh.eval_plain_mesh_inertias(rho_water)

    xb, yb, zb = inertia.gravity_center

    disp_volume = inertia.mass / rho_water

    # Computing quantities from intersection polygons
    sigma0 = 0.  # \iint_{waterplane_area} dS = waterplane_area
    sigma1 = 0.  # \iint_{waterplane_area} x dS
    sigma2 = 0.  # \iint_{waterplane_area} y dS
    sigma3 = 0.  # \iint_{waterplane_area} xy dS
    sigma4 = 0.  # \iint_{waterplane_area} x^2 dS
    sigma5 = 0.  # \iint_{waterplane_area} y^2 dS

    xmin = []
    xmax = []
    ymin = []
    ymax = []

    polygons = clipper.closed_polygons
    for polygon in polygons:
        polyverts = clipper.clipped_crown_mesh.vertices[polygon]

        # TODO: voir si on conserve ce test...
        if np.any(np.fabs(polyverts[:, 2]) > 1e-3):
            print('The intersection polygon is not on the plane z=0')

        xi, yi = polyverts[0, :2]
        for (xii, yii) in polyverts[1:, :2]:
            dx = xii - xi
            dy = yii - yi
            px = xi + xii
            py = yi + yii

            sigma0 += dy * px
            sigma1 += dy * (px * px - xi * xii)
            sigma2 += dx * (py * py - yi * yii)
            sigma3 += dy * (py * px * px + yi * xi * xi + yii * xii * xii)
            sigma4 += dy * (xi * xi + xii * xii) * px
            sigma5 += dx * (yi * yi + yii * yii) * py

            xi, yi = xii, yii

        xmin.append(polyverts[:, 0].min())
        xmax.append(polyverts[:, 0].max())
        ymin.append(polyverts[:, 1].min())
        ymax.append(polyverts[:, 1].max())

    minx, maxx = [min(xmin), max(xmax)]
    miny, maxy = [min(ymin), max(ymax)]

    sigma0 /= 2
    sigma1 /= 6
    sigma2 /= -6
    sigma3 /= 24
    sigma4 /= 12
    sigma5 /= -12

    # Flotation surface
    waterplane_area = sigma0

    # Stiffness matrix coefficients that do not depend on the position of the gravity center
    rhog = rho_water * grav
    s33 = rhog * waterplane_area
    s34 = rhog * sigma2
    s35 = -rhog * sigma1
    s45 = -rhog * sigma3

    # Metacentric radius (Bouguer formulae)
    transversal_metacentric_radius = sigma5 / disp_volume  # Around Ox
    longitudinal_metacentric_radius = sigma4 / disp_volume  # Around Oy

    # Metacentric height
    a = zcog - zb  # BG
    gm_x = transversal_metacentric_radius - a
    gm_y = longitudinal_metacentric_radius - a

    # Stiffness matrix coefficients that depend on the position of the gravity center
    s44 = rhog * disp_volume * gm_x
    s55 = rhog * disp_volume * gm_y

    # Assembling stiffness matrix
    stiffness_matrix = np.array([[s33, s34, s35],
                                 [s34, s44, s45],
                                 [s35, s45, s55]], dtype=np.float)

    # Zeroing tiny coefficients
    stiffness_matrix[np.fabs(stiffness_matrix) < 1e-4] = 0.

    # Flotation center F:
    x_f = -s35 / s33
    y_f = s34 / s33
    # TODO: ajouter xf et yf dans le rapport hydro !!

    xmin, xmax, ymin, ymax, zmin, zmax = cmesh.axis_aligned_bbox

    # Storing data
    hs_data = dict()
    hs_data['wet_surface_area'] = wet_surface_area
    hs_data['disp_volume'] = disp_volume
    hs_data['disp_mass'] = rho_water * disp_volume
    hs_data['buoy_center'] = np.array([xb, yb, zb], dtype=np.float)
    hs_data['flotation_center'] = np.array([x_f, y_f, 0.], dtype=np.float)
    hs_data['waterplane_area'] = waterplane_area
    hs_data['transversal_metacentric_radius'] = transversal_metacentric_radius
    hs_data['longitudinal_metacentric_radius'] = longitudinal_metacentric_radius
    hs_data['gm_x'] = gm_x
    hs_data['gm_y'] = gm_y
    hs_data['stiffness_matrix'] = stiffness_matrix
    hs_data['lwl'] = maxx - minx
    hs_data['los'] = xmax - xmin
    hs_data['bos'] = ymax - ymin
    hs_data['draught'] = math.fabs(zmin)
    hs_data['fp'] = maxx
    hs_data['breadth'] = maxy - miny

    # TODO: we should better store the inertia object !
    inertia.shift_at_cog()
    hs_data['Ixx'] = inertia.xx
    hs_data['Iyy'] = inertia.yy
    hs_data['Izz'] = inertia.zz
    hs_data['Ixy'] = inertia.xy
    hs_data['Ixz'] = inertia.xz
    hs_data['Iyz'] = inertia.yz

    return hs_data


def disp_equilibrium(mesh, disp_tons, rho_water, grav, reltol=1e-6, verbose=False):
    if verbose:
        print("========================")
        print("Displacement equilibrium")
        print("========================")
        print("Target displacement:    {:.3f} tons".format(disp_tons))

    hs_data = dict()

    disp_kg = disp_tons * 1000

    itermax = 100
    iter = 1

    z_corr = 0.

    # Relaxation distance adapted to the size of the mesh
    _, _, _, _, zmin, zmax = mesh.axis_aligned_bbox
    z_relax = (zmax - zmin) / 5

    while True:

        if iter == itermax:
            print("No convergence after %s" % itermax)
            break

        hs_data = compute_hydrostatics(mesh, rho_water, grav, z_corr=z_corr)

        disp_volume = hs_data["disp_volume"]
        waterplane_area = hs_data["waterplane_area"]

        mass_residual = rho_water * disp_volume - disp_kg

        dz = mass_residual / (rho_water * waterplane_area)

        # Convergence criterion
        if math.fabs(mass_residual / disp_kg) < reltol:
            break

        # Relaxation of the correction: z correction is bounded by z_relax to avoid divergence
        if math.fabs(dz) > z_relax:
            dz = math.copysign(z_relax, dz)

        z_corr += dz

        iter += 1

    if verbose:
        print("Achieved displacement: {:.3f} tons ({:d} iterations)".format(hs_data["disp_mass"] / 1000., iter))

    return z_corr


def full_equilibrium(mesh, cog, disp_tons, rho_water, grav, reltol=1e-4, verbose=False):
    if verbose:
        print("========================")
        print("3DOF equilibrium")
        print("========================")
        print("Target diplacement:    {:.3f} tons".format(disp_tons))
        print("COG pos: ", cog)

    disp_kg = disp_tons * 1000

    rhog = rho_water * grav
    mg = disp_kg * grav

    wmesh = mesh.copy()

    # Initial equilibrium in displacement
    z_corr = disp_equilibrium(wmesh, disp_tons, rho_water, grav, reltol=reltol, verbose=False)
    rotmat_corr = np.eye(3)

    itermax = 100
    iter = 1

    while True:

        if iter == itermax:
            print("No convergence after %s" % itermax)
            break

        wcog = cog.copy()
        wcog = np.dot(rotmat_corr, wcog)
        wcog[2] += z_corr

        # Computing hydrostatics properties
        hs_data = compute_hydrostatics(mesh, rho_water, grav, rotmat_corr=rotmat_corr, z_corr=z_corr, zcog=wcog[2])

        # Computing the residue
        xg, yg, zg = wcog
        xb, yb, zb = hs_data["buoy_center"]
        disp_volume = hs_data['disp_volume']

        rhogv = rhog * disp_volume
        residue = np.array([rhogv - mg,
                            rhogv * yb - mg * yg,
                            -rhogv * xb + mg * xg])

        # Computing scale for nondimensionalization
        breadth = hs_data['breadth']
        lwl = hs_data['lwl']
        scale = np.array([mg, mg * breadth, mg * lwl])

        # Convergence criterion
        if np.all(np.fabs(residue / scale) < reltol):
            print("convergence in %u iterations" % iter)
            break

        # Solving for the next incremental correction
        tz, rx, ry = np.linalg.solve(hs_data["stiffness_matrix"], residue)

        # Updating corrections
        rotmat_corr = np.dot(cardan_to_rotmat(rx, ry, 0.), rotmat_corr)
        z_corr += tz

        iter += 1




    # rx, ry, _ = rotmat_to_cardan(rotmat_corr)
    # print("%f\t%f\t%f" % (z_corr, math.degrees(rx), math.degrees(ry)))


# TODO: reactiver plus tard !!! pour la relaxation
# rx_relax = 2 * math.pi / 180
# ry_relax = 2 * math.pi / 180
#
# # FIXME: manque le z_relax
#
# if math.fabs(rx) > rx_relax:
#     rx = math.copysign(rx_relax, rx)
#
# if math.fabs(ry) > ry_relax:
#     ry = math.copysign(ry_relax, ry)
# TODO: FIN reactiver plus tard !!!
