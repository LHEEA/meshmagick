#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""This module allows to perform hydrostatics computations on meshes"""

import numpy as np
import math

from .mesh_clipper import MeshClipper
from .rotations import cardan_to_rotmat, rotmat_to_cardan
from math import degrees, radians

from meshmagick import __version__ as version

__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-2015, Ecole Centrale de Nantes / D-ICE ENGINEERING"
__credits__ = "Francois Rongere"
__licence__ = "GPLv3"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@dice-engineering.com"
__status__ = "Development"


def compute_hydrostatics(mesh, cog, rho_water, grav, rotmat_corr=np.eye(3), z_corr=0., at_cog=False):
    wmesh = mesh.copy()
    wmesh.rotate_matrix(rotmat_corr)
    wmesh.translate_z(z_corr)

    wcog = cog.copy()
    wcog = np.dot(rotmat_corr, wcog)
    wcog[2] += z_corr

    clipper = MeshClipper(wmesh, assert_closed_boundaries=True, verbose=False)
    cmesh = clipper.clipped_mesh

    wet_surface_area = cmesh.faces_areas.sum()

    inertia = cmesh.eval_plain_mesh_inertias(rho_water)

    xb, yb, zb = inertia.gravity_center
    buoyancy_center = np.array([xb, yb, zb])

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

        if at_cog:
            # To get the stiffness expressed at cog
            polyverts[:, 0] -= xb
            polyverts[:, 1] -= yb

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
    zcog = wcog[2]
    a = zcog - zb  # BG
    transversal_metacentric_height = transversal_metacentric_radius - a
    longitudinal_metacentric_height = longitudinal_metacentric_radius - a

    # Stiffness matrix coefficients that depend on the position of the gravity center
    s44 = rhog * disp_volume * transversal_metacentric_height
    s55 = rhog * disp_volume * longitudinal_metacentric_height

    # Assembling stiffness matrix
    stiffness_matrix = np.array([[s33, s34, s35],
                                 [s34, s44, s45],
                                 [s35, s45, s55]], dtype=np.float)

    # Zeroing tiny coefficients
    stiffness_matrix[np.fabs(stiffness_matrix) < 1e-4] = 0.

    # Flotation center F:
    x_f = -s35 / s33
    y_f = s34 / s33

    if at_cog:
        # To get the flotation center expressed in the initial mesh frame
        x_f += xb
        y_f += yb

    waterplane_center = np.array([x_f, y_f, -z_corr])
    waterplane_center = np.dot(np.transpose(rotmat_corr), waterplane_center)

    # Buoyancy center in the initial frame
    buoyancy_center[2] -= z_corr
    buoyancy_center = np.dot(np.transpose(rotmat_corr), buoyancy_center)

    # bounding box
    xmin, xmax, ymin, ymax, zmin, zmax = cmesh.axis_aligned_bbox

    # Storing data
    hs_data = dict()
    hs_data['grav'] = grav
    hs_data['rho_water'] = rho_water
    hs_data["mesh"] = mesh.copy()
    hs_data["rotmat_eq"] = rotmat_corr
    hs_data["z_eq"] = z_corr
    hs_data["cog"] = cog.copy()

    hs_data['wet_surface_area'] = wet_surface_area
    hs_data['disp_volume'] = disp_volume
    hs_data['disp_mass'] = rho_water * disp_volume
    hs_data['buoyancy_center'] = buoyancy_center
    hs_data['waterplane_center'] = waterplane_center
    hs_data['waterplane_area'] = waterplane_area
    hs_data['transversal_metacentric_radius'] = transversal_metacentric_radius
    hs_data['longitudinal_metacentric_radius'] = longitudinal_metacentric_radius
    hs_data['transversal_metacentric_height'] = transversal_metacentric_height
    hs_data['longitudinal_metacentric_height'] = longitudinal_metacentric_height
    hs_data['stiffness_matrix'] = stiffness_matrix
    hs_data['lwl'] = maxx - minx
    hs_data['los'] = xmax - xmin
    hs_data['bos'] = ymax - ymin
    hs_data['draught'] = math.fabs(zmin)
    hs_data['fp'] = maxx
    hs_data['breadth'] = maxy - miny

    inertia.set_cog(wcog)
    inertia.shift_at_cog()
    hs_data['Ixx'] = inertia.xx
    hs_data['Iyy'] = inertia.yy
    hs_data['Izz'] = inertia.zz
    hs_data['Ixy'] = inertia.xy
    hs_data['Ixz'] = inertia.xz
    hs_data['Iyz'] = inertia.yz

    return hs_data


def disp_equilibrium(mesh, disp_tons, rho_water, grav, cog=np.zeros(3), reltol=1e-6, verbose=False):
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

        hs_data = compute_hydrostatics(mesh, cog, rho_water, grav, z_corr=z_corr)

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

    disp_kg = disp_tons * 1000

    rhog = rho_water * grav
    mg = disp_kg * grav

    wmesh = mesh.copy()

    # Initial equilibrium in displacement
    z_corr = disp_equilibrium(wmesh, disp_tons, rho_water, grav, reltol=reltol, verbose=False)
    rotmat_corr = np.eye(3)

    hs_data = dict()

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
        hs_data = compute_hydrostatics(mesh, cog, rho_water, grav, rotmat_corr=rotmat_corr, z_corr=z_corr)

        # Computing the residue
        xg, yg, zg = wcog
        buoyancy_center = hs_data["buoyancy_center"]
        buoyancy_center = np.dot(rotmat_corr, buoyancy_center)
        buoyancy_center[2] += z_corr
        xb, yb, zb = buoyancy_center
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
            print("\t>>> Convergence reached in %u iterations" % iter)
            break

        # Solving for the next incremental correction
        tz, rx, ry = np.linalg.solve(hs_data["stiffness_matrix"], residue)

        # Relaxation
        _, _, _, _, zmin, zmax = wmesh.axis_aligned_bbox
        z_relax = (zmax - zmin) / 5
        if math.fabs(tz) > z_relax:
            tz = math.copysign(z_relax, tz)

        rx_relax = ry_relax = radians(2)
        if math.fabs(rx) > rx_relax:
            rx = math.copysign(rx_relax, rx)
        if math.fabs(ry) > ry_relax:
            ry = math.copysign(ry_relax, ry)

        # Updating corrections
        rotmat_corr = np.dot(cardan_to_rotmat(rx, ry, 0.), rotmat_corr)
        z_corr += tz

        iter += 1

    return z_corr, rotmat_corr


def get_hydrostatic_report(hs_data):
    """Returns a hydrostatic report for the current configuration

    Returns
    -------
    str
    """

    def hspace():
        return '\n'

    def build_line(text, data, precision=3, dtype='f'):
        # TODO: ajouter unit
        textwidth = 40
        try:
            line = '\t{:-<{textwidth}}>  {:< .{precision}{dtype}}\n'.format(str(text).upper(), data,
                                                                            precision=precision,
                                                                            textwidth=textwidth,
                                                                            dtype=dtype
                                                                            )
        # ValueError for Python 2
        # TypeError for Python 3
        except (ValueError, TypeError):
            if isinstance(data, (np.ndarray, list, tuple)):
                data = np.asarray(data)
                if data.ndim == 1:
                    data_str = ''.join(['{:< 10.{precision}{dtype}}'.format(val, precision=precision, dtype=dtype)
                                        for val in data])
                    line = '\t{:-<{textwidth}}>  {}\n'.format(str(text).upper(), data_str, textwidth=textwidth)

                    # else:
                    #     print data

        return line

    print('\n\n')

    msg = "HYDROSTATIC REPORT (Meshmagick version %s)\n\n" % version

    msg += '\tCorrections made on initial mesh:\n'
    msg += build_line('Z correction (M)', hs_data['z_eq'], precision=3)
    heel, trim, _ = rotmat_to_cardan(hs_data['rotmat_eq'])
    msg += build_line('Heel correction (deg)', degrees(heel), precision=3)
    msg += build_line('Trim correction (deg)', degrees(trim), precision=3)

    msg += hspace()
    msg += build_line('Gravity acceleration (M/S**2)', hs_data['grav'], precision=2)
    msg += build_line('Density of water (kg/M**3)', hs_data['rho_water'], precision=1)

    msg += hspace()
    msg += build_line('Waterplane area (M**2)', hs_data['waterplane_area'], precision=1)
    msg += build_line('Waterplane center (M)', hs_data['waterplane_center'], precision=3)
    msg += build_line('Wetted Surface Area (M**2)', hs_data['wet_surface_area'], precision=1)
    msg += build_line('Displacement volume (M**3)', hs_data['disp_volume'], precision=3)
    msg += build_line('Displacement mass (tons)', hs_data['disp_mass'] / 1000., precision=3)
    msg += build_line('Buoyancy center (M)', hs_data["buoyancy_center"], precision=3)
    msg += build_line('Center of gravity (M)', hs_data['cog'], precision=3)

    msg += hspace()
    msg += build_line('Draft (M)', hs_data['draught'], precision=3)  # TODO
    msg += build_line('Length overall submerged (M)', hs_data['los'], precision=2)
    msg += build_line('Breadth overall submerged (M)', hs_data['bos'], precision=2)
    msg += build_line('Length at Waterline LWL (M)', hs_data['lwl'], precision=2)
    msg += build_line('Forward perpendicular FP (M)', hs_data['fp'], precision=2)

    msg += hspace()
    msg += build_line('Transversal metacentric radius (M)', hs_data['transversal_metacentric_radius'], precision=3)
    msg += build_line('Longitudinal metacentric radius (M)', hs_data['longitudinal_metacentric_radius'], precision=3)
    msg += build_line('Transversal metacentric height GMt (M)', hs_data['transversal_metacentric_height'], precision=3)
    msg += build_line('Longitudinal metacentric height GMl (M)', hs_data['longitudinal_metacentric_height'],
                      precision=3)
    if hs_data['transversal_metacentric_height'] < 0.:
        msg += '\t<<<  TRANSVERSALLY UNSTABLE  >>>\n'
    if hs_data['longitudinal_metacentric_height'] < 0.:
        msg += '\t<<<  LONGITUDINALLY UNSTABLE  >>>\n'

    msg += hspace()
    msg += '\tHYDROSTATIC STIFFNESS COEFFICIENTS (at COG horiz. Pos.):\n'
    msg += build_line('K33 (N/M)', hs_data['stiffness_matrix'][0, 0], precision=4, dtype='E')
    msg += build_line('K34 (N)', hs_data['stiffness_matrix'][0, 1], precision=4, dtype='E')
    msg += build_line('K35 (N)', hs_data['stiffness_matrix'][0, 2], precision=4, dtype='E')
    msg += build_line('K44 (N.M)', hs_data['stiffness_matrix'][1, 1], precision=4, dtype='E')
    msg += build_line('K45 (N.M)', hs_data['stiffness_matrix'][1, 2], precision=4, dtype='E')
    msg += build_line('K55 (N.M)', hs_data['stiffness_matrix'][2, 2], precision=4, dtype='E')

    # Il faut faire une correction avec le plan de la flottaison de certains coeffs
    msg += hspace()
    msg += '\tINERTIAS:\n'
    msg += build_line('Ixx', hs_data['Ixx'], precision=3, dtype='E')
    msg += build_line('Ixy', hs_data['Ixy'], precision=3, dtype='E')
    msg += build_line('Ixz', hs_data['Ixz'], precision=3, dtype='E')
    msg += build_line('Iyy', hs_data['Iyy'], precision=3, dtype='E')
    msg += build_line('Iyz', hs_data['Iyz'], precision=3, dtype='E')
    msg += build_line('Izz', hs_data['Izz'], precision=3, dtype='E')

    return msg
