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


def compute_hydrostatics(mesh, rho_water, grav, zcog=0.):
    clipper = MeshClipper(mesh, assert_closed_boundaries=True, verbose=False)
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

        wmesh = mesh.copy()
        wmesh.translate_z(z_corr)

        hs_data = compute_hydrostatics(wmesh, rho_water, grav)

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


def full_equilibrium(mesh, cog, disp_tons, rho_water, grav, reltol=1e-6, verbose=False):
    if verbose:
        print("========================")
        print("3DOF equilibrium")
        print("========================")
        print("Target diplacement:    {:.3f} tons".format(disp_tons))
        print("COG pos: ", cog)

    hs_data_NW = dict()

    disp_kg = disp_tons * 1000

    itermax = 10000
    iter = 1

    wmesh_NW = mesh.copy()
    wcog_NW = cog.copy()

    z_corr_working = disp_equilibrium(wmesh_NW, disp_tons, rho_water, grav, reltol=reltol, verbose=True)
    z_corr_NW = z_corr_working

    rotmat_working = np.eye(3)

    transform_NW = HTransform(trans=[0, 0, z_corr_NW])

    while True:

        if iter == itermax:
            print("No convergence after %s" % itermax)
            break

        # WORKING >>>>>>>>>>>>>>>>>>
        wmesh_working = mesh.copy()
        wcog_working = cog.copy()

        wcog_working = np.dot(rotmat_working, wcog_working)
        wcog_working[2] += z_corr_working
        wmesh_working.rotate_matrix(rotmat_working)
        wmesh_working.translate_z(z_corr_working)

        hs_data_working = compute_hydrostatics(wmesh_working, rho_water, grav)

        xg_working, yg_working, zg_working = wcog_working
        xb_working, yb_working, zb_working = hs_data_working["buoy_center"]
        disp_volume_working = hs_data_working['disp_volume']

        rhogv_working = rho_water * grav * disp_volume_working
        mg = disp_kg * grav

        residue_working = np.array([rhogv_working - mg,
                                    rhogv_working * yb_working - mg * yg_working,
                                    -rhogv_working * xb_working + mg * xg_working])
        # <<<<<<<<<<<<<<<< END WORKING

        # NOT WORKING >>>>>>>>>>>>>>>>>>>
        wmesh_NW = mesh.copy()
        wcog_NW = cog.copy()

        wmesh_NW.transform(transform_NW)
        wcog_NW = transform_NW * wcog_NW

        hs_data_NW = compute_hydrostatics(wmesh_NW, rho_water, grav)

        xg_NW, yg_NW, _ = wcog_NW
        xb_NW, yb_NW, _ = hs_data_NW["buoy_center"]
        disp_volume_NW = hs_data_NW['disp_volume']

        rhogv_NW = rho_water * grav * disp_volume_NW
        mg = disp_kg * grav

        residue_NW = np.array([rhogv_NW - mg,
                               rhogv_NW * yb_NW - mg * yg_NW,
                               -rhogv_NW * xb_NW + mg * xg_NW])

        # >>>>>>>>>>>>> NOT WORKING

        # Computing scale [ COMMON ]
        breadth_COMM = hs_data_NW['breadth']
        lwl_COMM = hs_data_NW['lwl']
        scale_COMM = np.array([mg, mg * breadth_COMM, mg * lwl_COMM])

        # >>>>>>>>>>>>>>>WORKING
        stiffness_matrix_working = hs_data_working["stiffness_matrix"]
        tz_working, rx_working, ry_working = np.linalg.solve(stiffness_matrix_working, residue_working)

        if np.all(np.fabs(residue_working / scale_COMM) < reltol):
            print("convergence in %u iterations" % iter)
            break

        dR_WORKING = cardan_to_rotmat(rx_working, ry_working, 0.)
        rotmat_working = np.dot(dR_WORKING, rotmat_working)
        z_corr_working += tz_working

        # TODO: remoonter au dessus de la resolution de system
        # WORKING
        # if np.all(np.fabs(residue_working / scale_COMM) < reltol):
        #     print("convergence in %u iterations" % iter)
        #     break
        # >>>>>>>>>>>> WORKING



        # >>>>>>>>>>>>>> NOT WORKING
        if np.all(np.fabs(residue_NW / scale_COMM) < reltol):
            print("convergence in %u iterations" % iter)
            break

        stiffness_matrix_NW = hs_data_NW["stiffness_matrix"]
        tz_NW, rx_NW, ry_NW = np.linalg.solve(stiffness_matrix_NW, residue_NW)

        dt_NW = HTransform()
        # dt.set_cardan(rx, ry, 0., [0., 0., tz])
        dt_NW.set_rotmat(dR_WORKING, np.dot(np.transpose(dR_WORKING), np.array([0., 0., tz_NW])))
        # print(dt)
        # print(transform_NW)
        transform_NW = dt_NW * transform_NW
        # print(transform_NW)

        iter += 1

    print("WORKING:")
    rx_working, ry_working, _ = rotmat_to_cardan(rotmat_working)
    print("%f\t%f\t%f" % (z_corr_working, math.degrees(rx_working), math.degrees(ry_working)))

    print("NOT WORKING:")
    rx_NW, ry_NW, _ = transform_NW.get_cardan()
    # print(transform_NW.trans)
    print("%f\t%f\t%f" % (transform_NW.trans[2], math.degrees(rx_NW), math.degrees(ry_NW)))

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
