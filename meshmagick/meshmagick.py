#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PYTHON_ARGCOMPLTETE_OK

# Python module to manipulate 2D meshes for hydrodynamics purposes

"""
This module contains utility function to manipulate, load, save and
convert surface mesh files used by the hydrodynamics community.
Two numpy arrays are manipulated in this module : _vertices and _faces.
_vertices is the array of nodes coordinates. It is an array of shape (nv, 3) where
nv is the number of nodes in the mesh.
_faces is the array of cell connectivities. It is an array of shape (nf, 4) where
nf is the number of cells in the mesh. Not that it has 4 columns as we consider
flat polygonal cells up to 4 edges (quads). Triangles are obtained by repeating
the first node at the end of the cell node ID list.

IMPORTANT NOTE:
IDs of _vertices are internally idexed from 0 in meshmagick. However, several mesh
file format use indexing starting at 1. This different convention might be transparent
to user and 1-indexing may not be present outside the I/O functions
"""

# TODO: move meshmagick.py at the root level of the project ?

from mesh import *
import mmio

import os, sys
import numpy as np
import numpy as np
import math
from datetime import datetime
from numpy.core.multiarray import dtype

__year__ = datetime.now().year

__author__     = "Francois Rongere"
__copyright__  = "Copyright 2014-%u, Ecole Centrale de Nantes" % __year__
__credits__    = "Francois Rongere"
__licence__    = "CeCILL"
__version__    = "1.0"
__maintainer__ = "Francois Rongere"
__email__      = "Francois.Rongere@ec-nantes.fr"
__status__     = "Development"


# _mult_surf = np.array([1/6., 1/6., 1/6., 1/12., 1/12., 1/12., 1/12., 1/12., 1/12., 1/20., 1/20., 1/20., 1/60., 1/60., 1/60.], dtype=float) # Defines the array coefficient to compute surface integrals efficiently
# def _get_surface_integrals(V, F, sum=True):
#     """_get_surface_integrals(_vertices, _faces, sum=True)
#
#     Internal function
#     Computes all the _faces' integrals that may be used in several computations such
#     as inertial properties of meshes. This function is partly based on the work of
#     David Eberly:
#     ...
#
#     Parameters:
#         V: ndarray
#             numpy array of the mesh's nodes coordinates
#         F: ndarray
#             numpy array of the mesh's _faces nodes connectivity
#         sum[optional]: bool
#             if let to True, the results will be summed over all _faces.
#             Otherwise, no sum will be performed and individual integral
#             on _faces will be returned
#
#     Return:
#         sint: ndarray
#             numpy array of shape (nf, 15) that contains different integrals
#             over the mesh _faces. If sum is False, nf is equal to the number
#             of _faces in the mesh. Otherwise, nf=1.
#
#             The different integrals that sint contains are:
#
#              sint[0]  = \int x dS
#              sint[1]  = \int y dS
#              sint[2]  = \int z dS
#              sint[3]  = \int yz dS
#              sint[4]  = \int xz dS
#              sint[5]  = \int xy dS
#              sint[6]  = \int x^2 dS
#              sint[7]  = \int y^2 dS
#              sint[8]  = \int z^2 dS
#              sint[9]  = \int x^3 dS
#              sint[10] = \int y^3 dS
#              sint[11] = \int z^3 dS
#              sint[12] = \int x^2y dS
#              sint[13] = \int y^2z dS
#              sint[14] = \int z^2x dS
#     """
#
#     # TODO : put reference to the Eberly's work in the docstring
#
#     nf = F.shape[0]
#
#     if sum:
#         sint = np.zeros(15, dtype=float)
#     else:
#         sint = np.zeros((nf, 15), dtype=float)
#
#     tri1 = [0, 1, 2]
#     tri2 = [0, 2, 3]
#
#     cross = np.cross
#
#     sint_tmp = np.zeros(15)
#     for (iface, face) in enumerate(F):
#         sint_tmp *= 0.
#         # sint_tmp = np.zeros(15) # FIXME : Essai, la version precedente serait mieux !
#
#         if face[0] == face[-1]:
#             nb = 1
#         else:
#             nb = 2
#
#         _vertices = V[face]
#
#         # Loop on triangles of the face
#         for itri in xrange(nb):
#             if itri == 0:
#                 triangle = tri1
#             else:
#                 triangle = tri2
#
#             V0, V1, V2 = _vertices[triangle]
#             x0, y0, z0 = V0
#             x1, y1, z1 = V1
#             x2, y2, z2 = V2
#
#             d0, d1, d2 = cross(V1-V0, V2-V0)
#             e1_c_e2 = math.sqrt(d0**2 + d1**2 + d2**2)
#
#             temp0 = V0 + V1
#             f1 = temp0 + V2
#             temp1 = V0*V0
#             temp2 = temp1 + V1*temp0
#             f2 = temp2 + V2*f1
#             f3 = V0*temp1 + V1*temp2 + V2*f2
#             g0 = f2 + V0*(f1+V0)
#             g1 = f2 + V1*(f1+V1)
#             g2 = f2 + V2*(f1+V2)
#
#             yz = z0*(4*y0 - y1 - y2)  - y0*(z1+z2) + 3*(y1*z1 + y2*z2)
#             xz = x0*(4*z0 - z1 - z2)  - z0*(x1+x2) + 3*(x1*z1 + x2*z2)
#             xy = y0*(4*x0 - x1 - x2)  - x0*(y1+y2) + 3*(x1*y1 + x2*y2)
#
#             # Update integrals
#             sint_tmp[0] += d0 * f1[0] # order 1 in vol, x in surf
#             sint_tmp[1] += d1 * f1[1] # order 1 in vol, y in surf
#             sint_tmp[2] += d2 * f1[2] # order 1 in vol, z in surf
#
#             sint_tmp[3] +=  d0 * yz # order yz in surf
#             sint_tmp[4] +=  d1 * xz # order xz in surf
#             sint_tmp[5] +=  d2 * xy # order xy in surf
#
#             sint_tmp[6] += d0 * f2[0] # order x in vol, x**2 in surf
#             sint_tmp[7] += d1 * f2[1] # order y in vol, y**2 in surf
#             sint_tmp[8] += d2 * f2[2] # order z in vol, z**2 in surf
#             sint_tmp[9] += d0 * f3[0] # order x**2 in vol, x**3 in surf
#             sint_tmp[10] += d1 * f3[1] # order y**2 in vol, y**3 in surf
#             sint_tmp[11] += d2 * f3[2] # order z**2 in vol, z**3 in surf
#             sint_tmp[12] += d0 * (y0*g0[0] + y1*g1[0] + y2*g2[0]) # order xy in vol, x**2*y in surf
#             sint_tmp[13] += d1 * (z0*g0[1] + z1*g1[1] + z2*g2[1]) # order yz in vol, y**2*z in surf
#             sint_tmp[14] += d2 * (x0*g0[2] + x1*g1[2] + x2*g2[2]) # order zx in vol, z**2*x in surf
#
#             if sum:
#                 sint += sint_tmp
#             else:
#                 sint[iface] = sint_tmp
#
#     if sum:
#         sint *= _mult_surf
#     else:
#         sint = np.array([sint[j]*_mult_surf for j in xrange(nf)], dtype=float)
#
#     return sint
#
#
# def get_mass_cog(V, F, rho=1.):
#     """get_mass_cog(_vertices, _faces, rho=1.)
#
#     Returns the mass and the center of gravity of a mesh
#
#     Parameters:
#         V: ndarray
#             numpy array of the coordinates of the mesh's nodes
#         F: ndarray
#             numpy array of the _faces' nodes connectivities
#         rho[optional]: float
#             specifies the density of the material enclosed by the mesh
#
#     Returns:
#
#     """
#     # TODO: allow to specify as mush as options as in get_inertial_properties... or remove this function !
#
#     return get_inertial_properties(V, F, rho=rho)[:2]
#
#
# _mult_vol = np.array([1., 1., 1., 1., 1., 1., 1/2., 1/2., 1/2., 1/3., 1/3., 1/3., 1/2., 1/2., 1/2.]) # Defines the array coefficient to compute volume integrals on meshes
# def get_inertial_properties(V, F, rho=7500., mass=None, thickness=None, shell=False, verbose=False):
#     """get_inertial_properties(_vertices, _faces, rho=7500., mass=None, thickness=None, shell=False, verbose=False)
#
#     Returns the inertial properties of a mesh. The mesh may be considred as being
#     filled with homogeneous material or as being a shell.
#
#     Parameters:
#         V: ndarray
#             numpy array of the coordinates of the mesh's nodes
#         F: ndarray
#             numpy array of the _faces' nodes connectivities
#         rho[optional]: float
#             the density of the material. By default, it is the steel density
#             (7500 kg/m**3)
#         mass[optional]: float
#             the mass of the mesh. If it is specified, it will overwrite
#             the density (even if explicitely given)
#         thickness[optional]: float
#             specifies the thickness of the hull. Used if shell is set to True.
#         shell[optional]: bool
#             if set to True, the mesh will be considered as a shell
#         verbose[optional]: bool
#             if set to True, the function will display a report on inertial
#             properties of the mesh
#
#     Returns:
#         mass: float
#             The mass of the mesh (computed or given...)
#         cog: ndarray
#             The coordinates of the center of gravity of the mesh
#         inertia_matrix: ndarray
#             The inertia matrix of the mesh expressed at the center
#              of gravity
#
#     """
#     # TODO : allow to specify a reduction point...
#     # The default density rho is that of steel
#
#     tol = 1e-8
#
#     # FIXME : la gestion des options n'est pas claire ici !!! Remettre les choses a plat
#
#     if shell: # The geometry is a shell with homogeneous thickness and density.
#         areas, normals = get_all_faces_properties(V, F)[:2]
#         St = areas.sum() # Total surface
#         # The geometry is considered as being a shell with thickness given
#         if mass is None:
#
#             # FIXME : may not work if thickness is not given !!
#
#             # if thickness is None:
#             #     # Assuming a standard thickness of 1cm
#             #     thickness = 1e-2
#
#             sigma = thickness * rho
#             mass = St*sigma
#
#         else:# A mass has been specified
#             sigma = mass/St # Surfacic density
#
#             if thickness is None:
#                 # Computing the equivalent thickness
#                 thickness = sigma / rho
#             else:
#                 # thickness has been given, overwriting the density of the medium accordingly
#                 rho_tmp = sigma / thickness
#                 rho = rho_tmp
#
#         # Getting surface integrals
#         sint = _get_surface_integrals(V, F, sum=False)
#
#         normals[normals==0.] = 1. # To avoid division by zero
#         # Correcting these integrals by normals
#         sint[:, :3] /= normals
#         sint[:, 3:6] /= normals
#         sint[:, 6:9] /= normals
#         sint[:, 9:12] /= normals
#         sint[:, 12:15] /= normals
#
#         nu = sint.sum(axis=0)
#
#         cog = np.array([nu[0], nu[1], nu[2]], dtype=np.float) * sigma / mass
#
#         inertia_matrix = np.array([
#             [nu[7]+nu[8] ,   -nu[5]   ,   -nu[4]],
#             [  -nu[5]    , nu[6]+nu[8],   -nu[3]],
#             [  -nu[4]    ,   -nu[3]   , nu[6]+nu[7]]
#         ], dtype=np.float) * sigma
#
#         if verbose:
#             print '\nPrincipal inertia parameters report:'
#             print '------------------------------------\n'
#             print 'Total surface         : %f m**2' % St
#             print 'Thickness             : %f m' % thickness
#             print 'Density               : %f kg/m**3' % rho
#             print 'Surface density       : %f kg/m**2' % sigma
#
#     else:
#         # The geometry is full
#         sint = _get_surface_integrals(V, F)
#
#         # Appliying multipliers
#         sint *= _mult_vol
#
#         # We take the mean of 3 possible computations from surface integrals
#         vol = (sint[0] + sint[1] + sint[2]) / 3.
#         cog = np.array([sint[6]/vol, sint[7]/vol, sint[8]/vol], dtype=float)
#
#         # Inertia matrix is expressed for the moment in O
#         # xx = sint[10] + sint[11] - vol*(cog[1]**2 + cog[2]**2)
#         # yy = sint[9] + sint[11] - vol*(cog[2]**2 + cog[0]**2)
#         # zz = sint[9] + sint[10] - vol*(cog[0]**2 + cog[1]**2)
#         # xy = -(sint[12] - vol*cog[0]*cog[1])
#         # yz = -(sint[13] - vol*cog[1]*cog[2])
#         # xz = -(sint[14] - vol*cog[2]*cog[0])
#
#         # Inertia matrix expressed in cog
#         xx = sint[10] + sint[11]
#         yy = sint[9] + sint[11]
#         zz = sint[9] + sint[10]
#         xy = -sint[12]
#         yz = -sint[13]
#         xz = -sint[14]
#
#         mass = rho * vol
#         # The inertia matrix is expressed in
#         inertia_matrix = rho * np.array(
#             [
#                 [xx, xy, xz],
#                 [xy, yy, yz],
#                 [xz, yz, zz]
#             ], dtype=np.float)
#
#     # Cleaning
#     cog[np.fabs(cog) < tol] = 0.
#     inertia_matrix[np.fabs(inertia_matrix) < tol] = 0.
#
#     if verbose:
#         print 'Mass                  : %f kg' % mass
#         print 'COG                   : (%f, %f, %f) m' % tuple(cog)
#         print 'Inertia matrix in COG : '
#         print '\t%E, %E, %E\n\t%E, %E, %E\n\t%E, %E, %E\n' % tuple(inertia_matrix.flatten())
#
#     return mass, cog, inertia_matrix
#
#
# def transport_inertia_matrix(mass, cog, Ig, point, rot=np.eye(3, dtype=np.float)):
#     """transport_inertia_matrix(mass, cog, Ig, point, rot)
#
#     Performs the transport of the inertia matrix of a mesh at an other
#     reduction point
#
#     Parameters:
#         mass: float
#             The mass of the mesh
#         cog: ndarray
#             The coordinates of the center of gravity
#         Ig: ndarray
#             The 3x3 inertia matrix of the mesh expressed at cog
#         point: ndarray
#             The coordinates of the reduction point
#         rot: ndarray
#             The rotation matrix defining the orientation of the
#             new axis system with respect to the current
#
#     Returns:
#         Ipoint:
#             The 3x3 inertia matrix of the mesh expressed at the
#             new reduction point and in a frame rotated by rot with
#             respect to the initial coordinate system
#     """
#
#     point_cog = cog - point
#     Ipoint = rot * Ig * rot.T + \
#              mass * (np.eye(3, dtype=float) * np.dot(point_cog, point_cog)
#                      - np.outer(point_cog, point_cog))
#     return Ipoint
#
#
# def get_volume(V, F):
#     """get_volume(_vertices, _faces)
#
#     Returns the volume of the mesh
#
#     Parameters:
#         V: ndarray
#             numpy array of the coordinates of the mesh's nodes
#         F: ndarray
#             numpy array of the _faces' nodes connectivities
#
#     Returns:
#         vol: float
#             The volume of the mesh
#     """
#     return _get_surface_integrals(V, F)[0]
#
#
# def get_COM(V, F):
#     """get_COM(_vertices, _faces)
#
#     Returns the center of mass (center of gravity) of the mesh
#
#     Parameters:
#         V: ndarray
#             numpy array of the coordinates of the mesh's nodes
#         F: ndarray
#             numpy array of the _faces' nodes connectivities
#
#     Returns:
#         com: ndarray
#             Coordinates of the center of gravity of the mesh
#     """
#     # FIXME: la sortie est tres etonnante. Ne doit pas faire ce que ca dit !!!
#
#     return _get_surface_integrals(V, F)


def merge_duplicates(vertices, tol=1e-8, return_index=False):

    nv = vertices.shape[0]

    index = np.arange(nv)
    levels = [[0, nv]]

    for idim in xrange(3):
        coord = vertices[index, idim].copy()

        for level in levels:
            isort = coord.argsort()

            output = np.unique(coord, return_index=True, return_counts=True)





    if return_index:
        return vertices, index
    else:
        return vertices



def merge_duplicates(V, F=None, verbose=False, tol=1e-8, return_index=False):
    """merge_duplicates(_vertices, _faces, verbose=False, tol=1e-8)

    Returns a new node array where close nodes have been merged into one node (following tol). It also returns
    the connectivity array _faces with the new node IDs.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the _faces' nodes connectivities
        verbose[optional]: bool
            if set to True, displays information on the merge procedure
        tol[optional]: float
            the tolerance used to define nodes that are coincident and
            that have to be merged

    Returns:
        _vertices: ndarray
            numpy array of the coordinates of the mesh's nodes where
            every node is different
        _faces: ndarray
            numpy array of the _faces' nodes connectivities, accordingly
            to the new node list that has been merged
    """
    # TODO: Refaire la documentation --> les entrees sorties ont change !!

    # TODO : Set a tolerance option in command line arguments

    # This function is a bottleneck in the clipping routines
    # TODO: use np.unique to cluster groups --> acceleration !!

    if verbose:
        print "* Removing duplicate _vertices:"
    nv, nbdim = V.shape

    levels = [0, nv]
    Vtmp = []
    iperm = np.arange(nv)

    for dim in range(nbdim):
        # Sorting the first dimension
        values = V[:, dim].copy()
        if dim > 0:
            values = values[iperm]
        levels_tmp = []
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            if istop-istart > 1:
                level_values = values[istart:istop]
                iperm_view = iperm[istart:istop]

                iperm_tmp = level_values.argsort()

                level_values[:] = level_values[iperm_tmp]
                iperm_view[:] = iperm_view[iperm_tmp]

                levels_tmp.append(istart)
                vref = values[istart]

                for idx in xrange(istart, istop):
                    cur_val = values[idx]
                    if np.abs(cur_val - vref) > tol:
                        levels_tmp.append(idx)
                        vref = cur_val

            else:
                levels_tmp.append(levels[ilevel])
        if len(levels_tmp) == nv:
            # No duplicate _vertices
            if verbose:
                print "\t -> No duplicate _vertices detected :)"
            break

        levels_tmp.append(nv)
        levels = levels_tmp

    else:
        # Building the new merged node list
        Vtmp = []
        newID = np.arange(nv)
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            Vtmp.append(V[iperm[istart]])
            newID[iperm[range(istart, istop)]] = ilevel
        V = np.array(Vtmp, dtype=float)
        # Applying renumbering to cells
        if F is not None:
            for cell in F:
                cell[:] = newID[cell]

        if verbose:
            nv_new = V.shape[0]
            print "\t -> Initial number of nodes : {:d}".format(nv)
            print "\t -> New number of nodes     : {:d}".format(nv_new)
            print "\t -> {:d} nodes have been merged".format(nv-nv_new)

    if F is not None:
        if return_index:
            return V, F, newID
        else:
            return V, F
    else:
        if return_index:
            return V, newID
        else:
            return V




#=======================================================================
#                         MESH MANIPULATION HELPERS
#=======================================================================


def _is_point_inside_polygon(point, poly):
    """_is_point_inside_polygon(point, poly)

    Internal function to Determine if a point is inside a given polygon.
    This algorithm is a ray casting method.

    Parameters:
        point: ndarray
            2D coordinates of the point to be tested
        poly: ndarray
            numpy array of shape (nv, 2) of 2D coordinates
            of the polygon
    """

    # FIXME : Do we have to repeat the first point at the last position
    # of polygon ???

    x = point[0]
    y = point[1]

    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in xrange(n):
        p2x, p2y = poly[i]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def generate_lid(V, F, max_area=None, verbose=False):
    """generate_lid(_vertices, _faces, max_area=None, verbose=False)

    Meshes the lid of a mesh with triangular _faces to be used in irregular frequency
    removal in BEM softwares. It clips the mesh againt the plane Oxy, extract the intersection
    polygon and relies on meshpy (that is a wrapper around the TRIANGLE meshing library).
    It is able to deal with moonpools.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the _faces' nodes connectivities
        max_area[optional]: float
            The maximum area of triangles to be generated
        verbose[optional]: bool
            If set to True, generates output along the processing

    """

    # TODO: Put the reference of TRIANGLE and meshpy (authors...) in the docstring

    # TODO: remove verbose mode and place it into the main of meshmagick !!!

    try:
        import meshpy.triangle as triangle
    except:
        raise ImportError, 'Meshpy has to be available to use the generate_lid() function'

    if verbose:
        print '\n--------------'
        print 'Lid generation'
        print '--------------\n'

    # Clipping the mesh with Oxy plane
    V, F, clip_infos = clip_by_plane(V, F, Plane(), infos=True)

    nv = V.shape[0]
    nf = F.shape[0]

    if max_area is None:
        max_area = get_all_faces_properties(V, F)[0].mean()

    # Analysing polygons to find holes
    polygons = clip_infos['PolygonsNewID']
    nb_pol = len(polygons)

    holes = []
    boundaries = []
    for ipoly, polygon in enumerate(polygons):
        points = V[polygon][:, :2]
        n = points.shape[0]
        # Testing the orientation of each polygon by computing the signed area of it
        sarea = np.array([points[j][0]*points[j+1][1] - points[j+1][0]*points[j][1] for j in xrange(n-1)],
                         dtype=np.float).sum()
        if sarea < 0.:
            holes.append(polygon)
        else:
            boundaries.append(polygon)

    nb_hole = len(holes)
    nb_bound = len(boundaries)

    hole_dict = dict([(j, []) for j in xrange(nb_bound)])
    if nb_hole > 0:
        if verbose:
            if nb_hole == 1:
                word = 'moonpool has'
            else:
                word = 'moonpools have'
            print '\t-> %u %s been detected' % (nb_hole, word)

        # TODO : getting a point inside the hole polygon

        def pick_point_inside_hole(hole):

            # First testing with the geometric center of the hole
            point = np.array(hole).sum(axis=0)/len(hole)
            if not _is_point_inside_polygon(point, hole):
                # Testing something else
                raise RuntimeError, 'The algorithm should be refined to more complex polygon topologies... up to you ?'

            return point

        # Assigning holes to boundaries
        if nb_bound == 1 and nb_hole == 1:
            # Obvious case
            hole_dict[0].append( ( 0, pick_point_inside_hole(V[holes[0]][:, :2]) ) )
        else:
            # We may do a more elaborate search
            for ihole, hole in enumerate(holes):
                P0 = V[hole[0]][:2]
                # Testing against all boundary polygons
                for ibound, bound in enumerate(boundaries):
                    if _is_point_inside_polygon(P0, V[bound][:, :2]):
                        hole_dict[ibound].append( ( ihole, pick_point_inside_hole(V[hole][:, :2]) ) )
                        break

    def round_trip_connect(start, end):
        return [(j, j+1) for j in xrange(start, end)] + [(end, start)]

    # Meshing every boundaries, taking into account holes
    for ibound, bound in enumerate(boundaries):

        nvp = len(bound)-1

        # Building the loop
        points = map(tuple, list(V[bound][:-1, :2]))

        edges = round_trip_connect(0, nvp-1)

        info = triangle.MeshInfo()

        if len(hole_dict[ibound]) > 0:
            for ihole, point in hole_dict[ibound]:
                hole = holes[ihole]
                points.extend(map(tuple, list(V[hole][:-1, :2])))
                edges.extend(round_trip_connect(nvp, len(points)-1))

                # Marking the point as a hole
                info.set_holes([tuple(point)])

        info.set_points(points)
        info.set_facets(edges)

        # Generating the lid
        mesh = triangle.build(info, max_volume=max_area, allow_boundary_steiner=False)

        mesh_points = np.array(mesh.points)
        nmp = len(mesh_points)
        mesh_tri = np.array(mesh.elements, dtype=np.int32)

        # Resizing
        nmt = mesh_tri.shape[0]
        mesh_quad = np.zeros((nmt, 4), dtype=np.int32)
        mesh_quad[:, :-1] = mesh_tri + nv
        mesh_quad[:, -1] = mesh_quad[:, 0]

        mesh_points_3D = np.zeros((nmp, 3))
        mesh_points_3D[:, :-1] = mesh_points

        # show(_vertices, _faces)
        # return

        # Adding the lid to the initial mesh
        V = np.append(V, mesh_points_3D, axis=0)
        nv += nmp
        F = np.append(F, mesh_quad, axis=0)
        nf += nmt

    # Merging duplicates
    V, F = merge_duplicates(V, F)

    if verbose:
        if nb_bound == 1:
            verb = 'lid has'
        else:
            verb = 'lids have'
        print "\n\t-> %u %s been added successfully\n" % (nb_bound, verb)

    return V, F

def fill_holes(V, F, verbose=False):

    import vtk

    if verbose:
        print "Filling holes"

    polydata = _build_vtkPolyData(V, F)

    fillHolesFilter = vtk.vtkFillHolesFilter()

    if vtk.VTK_MAJOR_VERSION <= 5:
        fillHolesFilter.SetInputConnection(polydata.GetProducerPort())
    else:
        fillHolesFilter.SetInputData(polydata)

    fillHolesFilter.Update()

    polydata_filled = fillHolesFilter.GetOutput()

    V, F = _dump_vtk(polydata_filled)

    if verbose:
        print "\t--> Done!"

    return V, F

def detect_features(V, F, verbose=True):

    mesh = Mesh(V, F, verbose=verbose)
    mesh.detect_features(verbose=verbose)

    return

def _build_polyline(curve):
    import vtk

    npoints = len(curve)

    points = vtk.vtkPoints()
    for point in curve:
        points.InsertNextPoint(point)

    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(npoints)

    for id in xrange(npoints):
        polyline.GetPointIds().SetId(id, id)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(cells)

    return polydata


# =======================================================================
#                         COMMAND LINE USAGE
# =======================================================================


def main():

    import argparse
    import sys

    try:
        import argcomplete

        acok = True
    except:
        acok = False

    parser = argparse.ArgumentParser(
        description="""  --  MESHMAGICK --
                    A python module and a command line utility to manipulate meshes from
                    different format used in hydrodynamics as well as for visualization.

                    The formats currently supported by meshmagick are :

                    *-----------*------------*-----------------*----------------------*
                    | File      | R: Reading | Software        | Keywords             |
                    | extension | W: writing |                 |                      |
                    *-----------*------------*-----------------*----------------------*
                    |   .mar    |    R/W     | NEMOH (1)       | nemoh, mar           |
                    |   .gdf    |    R/W     | WAMIT (2)       | wamit, gdf           |
                    |   .inp    |    R       | DIODORE (3)     | diodore-inp, inp     |
                    |   .DAT    |    W       | DIODORE (3)     | diodore-dat          |
                    |   .hst    |    R/W     | HYDROSTAR (4)   | hydrostar, hst       |
                    |   .nat    |    R/W     |    -            | natural, nat         |
                    |   .msh    |    R       | GMSH (5)        | gmsh, msh            |
                    |   .rad    |    R       | RADIOSS         | rad, radioss         |
                    |   .stl    |    R/W     |    -            | stl                  |
                    |   .vtu    |    R/W     | PARAVIEW (6)    | vtu                  |
                    |   .vtp    |    R/W     | PARAVIEW (6)    | vtp                  |
                    |   .vtk    |    R/W     | PARAVIEW (6)    | paraview-legacy, vtk |
                    |   .tec    |    R/W     | TECPLOT (7)     | tecplot, tec         |
                    *---------- *-----------------------------------------------------*

                    By default, Meshmagick uses the filename extensions to choose the
                    appropriate reader/writer. This behaviour might be bypassed using the
                    -ifmt and -ofmt optional arguments. When using these options, keywords
                    defined in the table above must be used as format identifiers.

                    (1) NEMOH is an open source BEM Software for seakeeping developped at
                        Ecole Centrale de Nantes (LHHEA)
                    (2) WAMIT is a BEM Software for seakeeping developped by WAMIT, Inc.
                    (3) DIODORE is a BEM Software for seakeeping developped by PRINCIPIA
                    (4) HYDROSTAR is a BEM Software for seakeeping developped by
                        BUREAU VERITAS
                    (5) GMSH is an open source meshing software developped by C. Geuzaine
                    and J.-_faces. Remacle
                    (6) PARAVIEW is an open source visualization software developped by
                        Kitware
                    (7) TECPLOT is a visualization software developped by Tecplot


                    """,
        epilog='--  Copyright 2014-2015  -  Francois Rongere  /  Ecole Centrale de Nantes  --',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('infilename', # TODO : voir pour un typ=file pour tester l'existence
                        help='path of the input mesh file in any format')

    parser.add_argument('-o', '--outfilename', type=str,
                        help="""path of the output mesh file. The format of
                         this file is determined from the extension given
                         """)

    parser.add_argument('-ifmt', '--input-format',
                        help="""Input format. Meshmagick will read the input file considering the
                         INPUT_FORMAT rather than using the extension
                         """)

    parser.add_argument('-ofmt', '--output-format',
                        help="""Output format. Meshmagick will write the output file considering
                        the OUTPUT_FORMAT rather than using the extension
                        """)

    parser.add_argument('-q', '--quiet',
                        help="""switch of verbosity of meshmagick""",
                        action='store_true')

    parser.add_argument('-i', '--info',
                        help="""extract informations on the mesh on the standard output""",
                        action='store_true')

    parser.add_argument('--quality',
                        help="""prints mesh quality""",
                        action='store_true')

    parser.add_argument('-t', '--translate',
                        nargs=3, type=float,
                        help="""translates the mesh in 3D
                        Usage -translate tx ty tz""")

    parser.add_argument('-tx', '--translatex', type=float,
                        help="""translates the mesh following the x direction""")

    parser.add_argument('-ty', '--translatey', type=float,
                        help="""translates the mesh following the y direction""")

    parser.add_argument('-tz', '--translatez', type=float,
                        help="""translates the mesh following the z direction""")

    parser.add_argument('-r', '--rotate',
                        nargs=3, type=float,
                        help="""rotates the mesh in 3D""")

    parser.add_argument('-rx', '--rotatex', type=float,
                        help="""rotates the mesh around the x direction""")

    parser.add_argument('-ry', '--rotatey', type=float,
                        help="""rotates the mesh around the y direction""")

    parser.add_argument('-rz', '--rotatez', type=float,
                        help="""rotates the mesh around the z direction""")

    parser.add_argument('-s', '--scale', type=float,
                        help="""scales the mesh. CAUTION : if used along
                         with a translation option, the scaling is done before
                        the translations. The translation magnitude should be set
                        accordingly to the newly scaled mesh.
                        """)

    parser.add_argument('-hn', '--heal_normals', action='store_true',
                        help="""Checks and heals the normals consistency and
                        verify if they are outward.
                        """)

    parser.add_argument('-fn', '--flip-normals', action='store_true',
                        help="""flips the normals of the mesh""")

    parser.add_argument('-p', '--plane', nargs='+', action='append',
                        help="""Defines a plane used by the --clip_by_plane and --symmetrize options.
                        It can be defined by the floats nx ny nz c where [nx, ny, nz]
                        is a normal vector to the plane and c defines its position
                        following the equation <N|X> = c with X a point belonging
                        to the plane.
                        It can also be defined by a string among [Oxy, Oxz, Oyz, \Oxy, \Oxz, \Oyz]
                        for quick definition. Several planes may be defined on the same command
                        line. Planes with a prepended '\' have normals inverted i.e. if Oxy has its
                        normal pointing upwards, \Oxy has its normal pointing downwards.
                        In that case, the planes are indexed by an integer starting by
                        0 following the order given in the command line.
                        """)

    parser.add_argument('-c', '--clip_by_plane', nargs='*', action='append',
                        help="""cuts the mesh with a plane. Is no arguments are given, the Oxy plane
                        is used. If an integer is given, it should correspond to a plane defined with
                        the --plane option. If a key string is given, it should be a valid key (see
                        help of --plane option for valid plane keys). A normal and a scalar could
                        also be given for the plane definition just as for the --plane option. Several
                        clipping planes may be defined on the same command line.""")

    parser.add_argument('-md', '--merge-duplicates', nargs='?', const='1e-8', default=None,
                        help="""merges the duplicate nodes in the mesh with the absolute tolerance
                        given as argument (default 1e-8)""")

    parser.add_argument('-tq', '--triangulate_quadrangles', action='store_true',
                        help="""Triangulate all quadrangle _faces by a simple splitting procedure.
                        Twho triangles are generated and from both solution, the one with the best
                        aspect ratios is kept. This option may be used in conjunction with a
                        mesh export in a format that only deal with triangular cells like STL format.""")

    parser.add_argument('-sym', '--symmetrize', nargs='*', action='append',
                        help="""Symmetrize the mesh by a plane defined wether by 4 scalars, i.e.
                        the plane normal vector coordinates and a scalar c such as N.X=c is the
                        plane equation (with X a point of the plane) or a string among Oxz, Oxy
                        and Oyz which are shortcuts for planes passing by the origin and whose
                        normals are the reference axes. Default is Oxz if only -y is specified.
                        Be careful that symmetry is applied before any rotation so as the plane
                        equation is defined in the initial frame of reference.""")

    # Arguments concerning the hydrostatics
    # TODO : permettre de specifier un fichier de sortie
    parser.add_argument('-hs', '--hydrostatics', action='store_true',
                        help="""Compute hydrostatics. When used with the --verbose option, hydrostatic
                        data will be shown on the command line.""") # TODO : completer l'aide avec les options

    # parser.add_argument('--mass', default=None, type=float,
    #                     help="""Specifies the mass of the device for hydrostatics calculations.
    #                     When used, an hydrostatic equilibrium is resolved. Depending on the join
    #                     use of --mass, --cog and --zcog options, several behavior may be obtained.
    #
    #                     1) --mass is given a value:
    #
    #                         This options trigs an equilibrium resolution.
    #
    #                         a) No options --cog or --zcog are given:
    #                             Equilibrium is searched in z only, no hydrostatic stiffness
    #                             matrix is computed but the stiffness in heave.
    #                         b) --zcog is given alone:
    #                             Equilibrium is searched in z only, the hydrostatics stiffness
    #                             matrix is given.
    #                         c) --cog is given alone:
    #                             Equilibrium is searched in 6 dof and the hydrostatics stiffness
    #                             matrix is given
    #
    #                         Note that if both options --cog and --zcog are given, --zcog will
    #                         be ignored but taken from --cog and the result will be the same as c).
    #
    #                     2) --mass is used:
    #
    #                         No equilibrium resolution is performed if no mass is specified.
    #                         Mesh position is considered as being the equilibrium position and
    #                         hydrostatics data are generated as-is. Mass of the device is then
    #                         an output and is considered to be the computed displacement * rho_water.
    #
    #                         a) No options --cog or --zcog are given:
    #                             Only the stiffness in heave can be calculated
    #                         b) --zcog is given alone:
    #                             The full stiffness matrix is calculated and the estimated cog
    #                             position is given.
    #                         c) --cog is given alone:
    #                             Idem b) but the residual may be not identically 0. as the mesh may
    #                             not be at equilibrium.
    #
    #                         Note that if both options --cog and --zcog are given, --zcog will
    #                         be ignored but taken from --cog and the result will be the same as c).
    #                     """)

    # parser.add_argument('--cog', nargs=3, default=None, type=float,
    #                     help="""Specifies the center of gravity to be used along with the
    #                     --hydrostatics option. See --mass for more information. See also the
    #                     --inertias option for side effects.
    #                     """)

    parser.add_argument('--zcog', default=None, type=float,
                        help="""Specifies the z position of the center of gravity to be used along
                        the --hydrostatics option.
                        """)

    # parser.add_argument('--zcog', default=None, type=float,
    #                     help="""Specifies the z position of the center of gravity to be used along
    #                     the --hydrostatics option. See --mass for more information. See also the
    #                     --inertias option for side effects.
    #                     """)

    # parser.add_argument('--anim', action='store_true',
    #                     help="""Generates animation files for paraview visualization of the equilibrium
    #                     resolution.
    #                     """)

    parser.add_argument('--rho-water', default=1023., type=float,
                        help="""Specifies the density of salt water. Default is 1023 kg/m**3.
                        """)

    parser.add_argument('-g', '--grav', default=9.81, type=float,
                        help="""Specifies the acceleration of gravity on the earth surface.
                        Default is 9.81 m/s**2.
                        """)

    # parser.add_argument('--rho-medium', default=7500., type=float,
    #                     help="""Specified the density of the medium used for the device. Default
    #                     is steel and is 7500 kg/m**3.
    #                     """)

    # parser.add_argument('--thickness', default=0.01, type=float,
    #                     help="""Specifies the thickness of the hull. This option is only used if
    #                     both the --inertias and --hull are used. Default is 0.01 m.
    #                     """)

    # parser.add_argument('-gz', '--gz-curves', nargs='?', const=5., default=None, type=float,
    #                     help=""" [EXPERIMENTAL] Computes the GZ curves with angle spacing given as argument.
    #                     Default is 5 degrees (if no argument given)
    #                     """)

    # TODO : permettre de rajouter des ballasts
    # parser.add_argument('--inertias', action='store_true', # TODO : specifier un point de calcul
    #                     help="""Compute the principal inertia properties of the mesh. By default,
    #                     the device is considered to be a hull. Then the --thickness and --rho-medium
    #                     options may be used to tune the properties.
    #                     If the --no-hull option is used, then the device will be considered to be
    #                     filled with the medium of density rho-medium. Be carefull that the default
    #                     medium is steel so that you may get a really heavy device. Please consider
    #                     to specify an other density with the --rho-medium option.
    #
    #                     Note that the inertia matrix is expressed at center of gravity
    #                     location.
    #
    #                     A side effect of this option is that for hydrostatics computations, the values
    #                     computed by this option will be used instead of other options --mass and --cog
    #                     that will be overriden. Be carefull that the device may need some ballast.
    #                     """)

    # parser.add_argument('--no-hull', action='store_true',
    #                     help="""Specifies that the device should be considered as being filled with
    #                     the material of density rho-medium. It is only used by the --inertias option.
    #                     """)

    parser.add_argument('--lid', nargs='?', const=1., default=None, type=float,
                        help="""Generate a triangle mesh lid on the mesh clipped by the Oxy plane.
                        """)

    parser.add_argument('--fill-holes', '-fh', action='store_true',
                        help="""Fill little holes by triangulation if any.
                        """)

    parser.add_argument('--show', action='store_true',
                        help="""Shows the input mesh in an interactive window""")

    parser.add_argument('--version', action='version',
                        version='meshmagick - version %s\n%s'%(__version__, __copyright__),
                        help="""Shows the version number and exit""")


    if acok:
        argcomplete.autocomplete(parser)

    # TODO : Utiliser des sous-commandes pour l'utilisation de meshmagick

    args, unknown = parser.parse_known_args()

    if args.quiet:
        verbose = False
    else:
        verbose = True

    if verbose:
        print '\n============================================='
        print 'meshmagick - version %s\n%s' % (__version__, __copyright__)
        print '============================================='

    # LOADING DATA FROM FILE
    if args.input_format is not None:
        format = args.input_format
    else:
        # Format based on extension
        _ ,ext = os.path.splitext(args.infilename)
        format = ext[1:].lower()
        if format == '':
            raise IOError, 'Unable to determine the input file format from its extension. Please specify an input format.'

    # Loading mesh elements from file
    if os.path.isfile(args.infilename):
        V, F = mmio.load_mesh(args.infilename, format)
        mesh = Mesh(V, F)
        # Ensuring triangles are following the right convention (last id = first id)
        mesh.heal_triangles()
        if verbose:
            mesh.verbose_on()
            print '%s successfully loaded' % args.infilename
    else:
        raise IOError, 'file %s not found'%args.infilename

    # Merge duplicate _vertices
    if args.merge_duplicates is not None:
        tol = float(args.merge_duplicates)
        mesh.merge_duplicates(tol=tol)

    # TODO : put that dict at the begining of the main function
    plane_str_list = {'Oxy':[0.,0.,1.],
                      'Oxz':[0.,1.,0.],
                      'Oyz':[1.,0.,0.],
                      '\Oxy':[0.,0.,-1.],
                      '\Oxz':[0.,-1.,0.],
                      '\Oyz':[-1.,0.,0.]}

    # Defining planes
    if args.plane is not None:
        nb_planes = len(args.plane)

        if verbose:
            if nb_planes == 1:
                verb = 'plane has'
            else:
                verb = 'planes have'
            print '\n%u %s been defined' % (nb_planes, verb)

        planes = [Plane() for i in xrange(nb_planes)]
        for (iplane, plane) in enumerate(args.plane):
            if len(plane) == 4:
                # plane is defined by normal and scalar
                try:
                    planes[iplane] = Plane(normal=map(float, plane[:3]), scalar=plane[3])
                    # planes[iplane].normal = np.array(map(float, plane[:3]), dtype=np.float)
                    # planes[iplane].c = float(plane[3])
                except:
                    raise AssertionError, 'Defining a plane by normal and scalar requires four scalars'

            elif len(plane) == 1:
                if plane_str_list.has_key(plane[0]):
                    planes[iplane].normal = np.array(plane_str_list[plane[0]], dtype=np.float)
                    planes[iplane].c = 0.
                else:
                    raise AssertionError, '%s key for plane is not known. Choices are [%s].' % (plane[0], ', '.join(plane_str_list.keys()) )
            else:
                raise AssertionError, 'Planes should be defined by a normal and a scalar or by a key to choose among [%s]' % (', '.join(plane_str_list.keys()))



    # Symmetrizing the mesh
    if args.symmetrize is not None:
        nb_sym = len(args.symmetrize)

        if verbose:
            if nb_sym == 1:
                verb = 'plane'
            else:
                verb = 'planes'
            print '\nMesh is being symmetrized by %u %s' % (nb_sym, verb)

        sym_plane = Plane()
        for plane in args.symmetrize:
            if len(plane) == 0:
                # Default symmetry by plane Oxz
                sym_plane.normal = np.array([0., 1., 0.], dtype=np.float)
                sym_plane.c = 0.
            elif len(plane) == 1:
                try:
                    # Plane ID
                    plane_id = int(plane[0])
                    if plane_id < nb_planes:
                        sym_plane = planes[plane_id]
                    else:
                        raise AssertionError, 'Plane with ID %u has not been defined with option --plane' % plane_id
                except:
                    # A key string
                    if plane_str_list.has_key(plane[0]):
                        sym_plane.normal = np.asarray(plane_str_list[plane[0]], dtype=np.float)
                        sym_plane.c = 0.
                    else:
                        raise AssertionError, 'Planes should be defined by a normal and a scalar or by a key to choose among [%s]' % (', '.join(plane_str_list.keys()))
            elif len(plane) == 4:
                # Plane defined by a normal and a scalar
                try:
                    sym_plane.normal = np.array(map(float, plane[:3]), dtype=np.float)
                    sym_plane.c = float(plane[3])
                except:
                    raise AssertionError, 'Defining a plane by normal and scalar requires four scalars'
            else:
                raise AssertionError, 'Unknown mean to define a plane for symmetry'
            mesh.symmetrize(sym_plane)
            if verbose:
                print '\t-> Done.'

    # Heal normals
    if args.heal_normals:
        if verbose:
            print '\nOPERATION: heal normals'
        mesh.heal_normals()
        if verbose:
            print '\t-> Done.'

    # Mesh translations
    if args.translate is not None:
        if verbose:
            print '\nOPERATION: Translation by [%f, %f, %f]' % tuple(args.translate)
        mesh.translate(args.translate)
        if verbose:
            print '\t-> Done.'

    if args.translatex is not None:
        if verbose:
            print '\nOPERATION: Translation by %f along X' % args.translatex
        mesh.translate_x(args.translatex)
        if verbose:
            print '\t-> Done.'

    if args.translatey is not None:
        if verbose:
            print '\nOPERATION: Translation by %f along Y' % args.translatey
        mesh.translate_y(args.translatey)
        if verbose:
            print '\t-> Done.'

    if args.translatez is not None:
        if verbose:
            print '\nOPERATION: Translation by %f along Z' % args.translatez
        mesh.translate_z(args.translatez)
        if verbose:
            print '\t-> Done.'

    # Mesh rotations
    if args.rotate is not None:
        if verbose:
            print '\nOPERATION: Rotation by [%f, %f, %f]' % tuple(args.rotate)
        mesh.rotate(args.rotate*math.pi/180.)
        if verbose:
            print '\t-> Done.'

    if args.rotatex is not None:
        if verbose:
            print '\nOPERATION: Rotation by %f around X (Roll)' % args.rotatex
        mesh.rotate_x(args.rotatex*math.pi/180.)
        if verbose:
            print '\t-> Done.'

    if args.rotatey is not None:
        if verbose:
            print '\nOPERATION: Rotation by %f around Y (Pitch)' % args.rotatey
        mesh.rotate_y(args.rotatey * math.pi / 180.)
        if verbose:
            print '\t-> Done.'

    if args.rotatez is not None:
        if verbose:
            print '\nOPERATION: Rotation by %f around Z (Yaw)' % args.rotatez
        mesh.rotate_z(args.rotatez * math.pi / 180.)
        if verbose:
            print '\t-> Done.'

    if args.scale is not None:
        if verbose:
            print '\nOPERATION: Scaling by %f' % args.scale
        mesh.scale(args.scale)
        if verbose:
            print '\t-> Done.'

    if args.flip_normals:
        if verbose:
            print '\nOPERATION: Flipping normals'
        mesh.flip_normals()
        if verbose:
            print '\t-> Done.'

    if args.triangulate_quadrangles:
        mesh.triangulate_quadrangles()

    # Clipping the mesh
    if args.clip_by_plane is not None:
        clipping_plane = Plane()
        nb_clip = len(args.clip_by_plane)

        if verbose:
            if nb_clip == 1:
                verb = 'plane'
            else:
                verb = 'planes'
            print '\nMesh is being clipped by %u %s' % (nb_clip, verb)

        for plane in args.clip_by_plane:
            if len(plane) == 0:
                # Default clipping plane Oxy
                clipping_plane.normal = np.array([0., 0., 1.], dtype=np.float)
                clipping_plane.c = 0.
            elif len(plane) == 1:
                try:
                    # Plane ID
                    plane_id = int(plane[0])
                    if plane_id < nb_planes:
                        clipping_plane = planes[plane_id]
                    else:
                        raise AssertionError, 'Plane with ID %u has not been defined with option --plane' % plane_id
                except:
                    # A key string
                    if plane_str_list.has_key(plane[0]):
                        clipping_plane.normal = np.asarray(plane_str_list[plane[0]], dtype=np.float)
                        clipping_plane.c = 0.
                    else:
                        raise AssertionError, 'Planes should be defined by a normal and a scalar or by a key to choose among [%s]' % (', '.join(plane_str_list.keys()))
            elif len(plane) == 4:
                # Plane defined by a normal and a scalar
                try:
                    clipping_plane.normal = np.array(map(float, plane[:3]), dtype=np.float)
                    clipping_plane.c = float(plane[3])
                except:
                    raise AssertionError, 'Defining a plane by normal and scalar requires four scalars'
            else:
                raise AssertionError, 'Unknown mean to define a plane for clipping'
            mesh = mesh.clip(clipping_plane)
            if verbose:
                print '\t-> Done.'

    # Compute principal inertia parameters
    # if args.inertias:
    #     # TODO : completer l'aide avec la logique de cette fonction !!
    #     if verbose:
    #         print '\n------------------'
    #         print 'Computing inertias'
    #         print '------------------'
    #     if args.no_hull:
    #         hull = False
    #     else:
    #         hull = True
    #
    #     mass, cog, inertia_matrix = get_inertial_properties(_vertices, _faces,
    #                                     rho=args.rho_medium,
    #                                     mass=args.mass,
    #                                     thickness=args.thickness,
    #                                     shell=hull,
    #                                     verbose=verbose)
    #     # Replacing values in command line arguments in the eventuality of hydrostatics computations
    #     args.mass = mass
    #     args.cog = cog
    #     if verbose:
    #         print '\t-> Done.'

    # if args.gz_curves is not None:
    #     if verbose:
    #         print '\n-------------------'
    #         print 'Computing GZ curves'
    #         print '-------------------'
    #
    #     spacing = args.gz_curves
    #     try:
    #         import hydrostatics as hs
    #     except:
    #         raise ImportError, '--gz-curves option relies on the hydrostatics module that can not be found'
    #
    #     # if args.hydrostatics:
    #     #     raise RuntimeError, """GZ computations can not be performed at the same time as a hydrostatics equilibrium
    #     #                            resolution as it needs a full mesh to perform clipping at different angles"""
    #
    #     if args.zcog is None:
    #         raise RuntimeError, 'For the GZ computations, the --zcog option is mandatory'
    #
    #     hsMesh = hs.HydrostaticsMesh(_vertices, _faces, rho_water=args.rho_water, g=args.grav)
    #     hs.get_GZ_curves(hsMesh, args.zcog,
    #                      spacing=spacing,
    #                      rho_water=args.rho_water,
    #                      g=args.grav,
    #                      verbose=verbose)
    #     if verbose:
    #         print '\t-> Done.'


    # Compute hydrostatics
    # ESSAI
    if args.hydrostatics: # TODO: A remettre en place
        try:
            import hydrostatics as hs
        except:
            raise ImportError, '--hydrostatics option relies on the hydrostatics module that can not be found'

        if args.zcog is None:
            raise ValueError, 'The hydrostatics option shall be used along with the --zcog option for the hydrostatic stiffness matrix to be computed'

        output_hs = hs.compute_hydrostatics(mesh, args.zcog, verbose=verbose)
        mesh = output_hs['mesh_hs']
        # _vertices = outputHS['Vc']
        # _faces = outputHS['Fc']

    # if args.hydrostatics:
    #     try:
    #         import hydrostatics as hs
    #     except:
    #         raise ImportError, '--hydrostatics option relies on the hydrostatics module that can not be found'
    #
    #     if args.anim:
    #         anim=True
    #     else:
    #         anim=False
    #
    #     if args.cog is None:
    #         cog = None
    #     else:
    #         cog = np.asarray(args.cog, dtype=np.float)
    #
    #     # TODO : Revoir la structure afin de ne jouer que sur l'objet !!
    #     hsMesh = hs.HydrostaticsMesh(_vertices, _faces, rho_water=args.rho_water, g=args.grav)
    #     _vertices, _faces = hs.get_hydrostatics(hsMesh,
    #                                mass=args.mass,
    #                                zcog=args.zcog,
    #                                cog=cog,
    #                                rho_water=args.rho_water,
    #                                g=args.grav,
    #                                anim=anim,
    #                                verbose=verbose)[:2]


    # Lid generation on a clipped mesh
    if args.lid is not None:# TODO: A remettre en place
        V, F = generate_lid(V, F, max_area=args.lid, verbose=verbose)

    if args.fill_holes:# TODO: A remettre en place
        V, F = fill_holes(V, F, verbose=verbose)


    # WARNING : No more mesh modification should be released from this point until the end of the main

    if args.info:
        print mesh

    if args.quality:
        mesh.print_quality()

    if args.show:
        mesh.show()

    if args.outfilename is None:
        base, ext = os.path.splitext(args.infilename)
        write_file = False
        # if write_file:
        #     args.outfilename = '%s_modified%s' % (base, ext)
        # Case where only the output format is given
        if args.output_format is not None:
            write_file = True
            args.outfilename = '%s.%s' % (base, args.output_format)
    else:
        write_file = True


    # Writing an output file
    if write_file:

        if args.output_format is not None:
            format = args.output_format
        else:
            if args.outfilename is None:
                # We base the output format on the input format used
                if args.input_format is not None:
                    format = args.input_format
                else:
                    format = os.path.splitext(args.infilename)[1][1:].lower()
                    if not mmio.know_extension(format):
                        raise IOError, 'Could not determine a format from input file extension, please specify an input format or an extension'
            else:
                format = os.path.splitext(args.outfilename)[1][1:].lower()

        if verbose:
            print 'Writing %s' % args.outfilename
        mmio.write_mesh(args.outfilename, mesh._vertices, mesh._faces, format)
        if verbose:
            print '\t-> Done.'

    if verbose:
        print '\n============================================================='
        print 'meshmagick - version %s\n%s' % (__version__, __copyright__)
        print 'Maintainer : %s <%s>' % (__maintainer__, __email__)
        print 'Good Bye!'
        print '============================================================='

if __name__ == '__main__':
    main()
