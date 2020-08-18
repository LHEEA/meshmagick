#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PYTHON_ARGCOMPLTETE_OK

# Python module to manipulate 2D meshes for hydrodynamics purposes
# TODO: Change the following docstring as it is no more up to date
"""
This module contains utility function to manipulate, load, save and
convert surface mesh files used by the hydrodynamics community.
Two numpy arrays are manipulated in this module : vertices and faces.
vertices is the array of nodes coordinates. It is an array of shape (nv, 3) where
nv is the number of nodes in the mesh.
faces is the array of cell connectivities. It is an array of shape (nf, 4) where
nf is the number of cells in the mesh. Not that it has 4 columns as we consider
flat polygonal cells up to 4 edges (quads). Triangles are obtained by repeating
the first node at the end of the cell node ID list.

IMPORTANT NOTE:
IDs of _vertices are internally idexed from 0 in meshmagick. However, several mesh
file format use indexing starting at 1. This different convention might be transparent
to user and 1-indexing may not be present outside the I/O functions
"""

# TODO: move meshmagick.py at the root level of the project ?

import os, sys
import numpy as np
import math
from datetime import datetime
from warnings import warn
from time import strftime

from .mesh import *
from . import mmio
from .mesh_clipper import MeshClipper
from . import hydrostatics as hs
import argparse
from . import densities

__year__ = datetime.now().year

__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-%u, Ecole Centrale de Nantes" % __year__
__credits__ = "Francois Rongere"
__licence__ = "GPLv3"
__version__ = "2.0.1"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@ec-nantes.fr"
__status__ = "Development"

__all__ = ['main']

# =======================================================================
#                         MESH MANIPULATION HELPERS
# =======================================================================
# TODO: those functions should disappear from this module

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
    # TODO: place this code into a utils module
    # FIXME : Do we have to repeat the first point at the last position
    # of polygon ???

    x = point[0]
    y = point[1]

    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n):
        p2x, p2y = poly[i]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

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

    # TODO: rely on the wrapper done for the triangle lib that has been done in cython and no more on meshpy.

    # TODO: Put the reference of TRIANGLE and meshpy (authors...) in the docstring

    # TODO: remove verbose mode and place it into the main of meshmagick !!!

    # TODO: Faire de cette fonction une methode dans mesh ??? --> non on a un autre module qui wrappe triangle avec
    # cython...

    try:
        import meshpy.triangle as triangle
    except:
        raise ImportError('Meshpy has to be available to use the generate_lid() function')

    if verbose:
        print('\n--------------')
        print('Lid generation')
        print('--------------\n')

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
        signed_area = np.array(
            [points[j][0] * points[j + 1][1] - points[j + 1][0] * points[j][1] for j in range(n - 1)],
            dtype=np.float).sum()
        if signed_area < 0.:
            holes.append(polygon)
        else:
            boundaries.append(polygon)

    nb_hole = len(holes)
    nb_bound = len(boundaries)

    hole_dict = dict([(j, []) for j in range(nb_bound)])
    if nb_hole > 0:
        if verbose:
            if nb_hole == 1:
                word = 'moonpool has'
            else:
                word = 'moonpools have'
            print(('\t-> %u %s been detected' % (nb_hole, word)))

        # TODO : getting a point inside the hole polygon

        def pick_point_inside_hole(hole):

            # First testing with the geometric center of the hole
            point = np.array(hole).sum(axis=0) / len(hole)
            if not _is_point_inside_polygon(point, hole):
                # Testing something else
                raise RuntimeError('The algorithm should be refined to more complex polygon topologies... up to you ?')

            return point

        # Assigning holes to boundaries
        if nb_bound == 1 and nb_hole == 1:
            # Obvious case
            hole_dict[0].append((0, pick_point_inside_hole(V[holes[0]][:, :2])))
        else:
            # We may do a more elaborate search
            for ihole, hole in enumerate(holes):
                P0 = V[hole[0]][:2]
                # Testing against all boundary polygons
                for ibound, bound in enumerate(boundaries):
                    if _is_point_inside_polygon(P0, V[bound][:, :2]):
                        hole_dict[ibound].append((ihole, pick_point_inside_hole(V[hole][:, :2])))
                        break

    def round_trip_connect(start, end):
        return [(j, j + 1) for j in range(start, end)] + [(end, start)]

    # Meshing every boundaries, taking into account holes
    for ibound, bound in enumerate(boundaries):

        nvp = len(bound) - 1

        # Building the loop
        points = list(map(tuple, list(V[bound][:-1, :2])))

        edges = round_trip_connect(0, nvp - 1)

        info = triangle.MeshInfo()

        if len(hole_dict[ibound]) > 0:
            for ihole, point in hole_dict[ibound]:
                hole = holes[ihole]
                points.extend(list(map(tuple, list(V[hole][:-1, :2]))))
                edges.extend(round_trip_connect(nvp, len(points) - 1))

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
        print(("\n\t-> %u %s been added successfully\n" % (nb_bound, verb)))

    return V, F


def fill_holes(V, F, verbose=False):
    import vtk

    if verbose:
        print("Filling holes")

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
        print("\t--> Done!")

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

    for id in range(npoints):
        polyline.GetPointIds().SetId(id, id)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(cells)

    return polydata


def list_medium():
    return ', '.join(densities.list_medium())

# =======================================================================
#                         COMMAND LINE USAGE
# =======================================================================


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

                +-----------+------------+-----------------+----------------------+
                | File      | R: Reading | Software        | Keywords             |
                | extension | W: writing |                 |                      |
                +===========+============+=================+======================+
                |   .mar    |    R/W     | NEMOH [#f1]_    | nemoh, mar           |
                +-----------+------------+-----------------+----------------------+
                |   .nem    |    R       | NEMOH [#f1]_    | nemoh_mesh, nem      |
                +-----------+------------+-----------------+----------------------+
                |   .gdf    |    R/W     | WAMIT [#f2]_    | wamit, gdf           |
                +-----------+------------+-----------------+----------------------+
                |   .inp    |    R       | DIODORE [#f3]_  | diodore-inp, inp     |
                +-----------+------------+-----------------+----------------------+
                |   .DAT    |    W       | DIODORE [#f3]_  | diodore-dat          |
                +-----------+------------+-----------------+----------------------+
                |   .hst    |    R/W     | HYDROSTAR [#f4]_| hydrostar, hst       |
                +-----------+------------+-----------------+----------------------+
                |   .nat    |    R/W     |    -            | natural, nat         |
                +-----------+------------+-----------------+----------------------+
                |   .msh    |    R       | GMSH [#f5]_     | gmsh, msh            |
                +-----------+------------+-----------------+----------------------+
                |   .rad    |    R       | RADIOSS         | rad, radioss         |
                +-----------+------------+-----------------+----------------------+
                |   .stl    |    R/W     |    -            | stl                  |
                +-----------+------------+-----------------+----------------------+
                |   .vtu    |    R/W     | PARAVIEW [#f6]_ | vtu                  |
                +-----------+------------+-----------------+----------------------+
                |   .vtp    |    R/W     | PARAVIEW [#f6]_ | vtp                  |
                +-----------+------------+-----------------+----------------------+
                |   .vtk    |    R/W     | PARAVIEW [#f6]_ | paraview-legacy, vtk |
                +-----------+------------+-----------------+----------------------+
                |   .tec    |    R/W     | TECPLOT [#f7]_  | tecplot, tec         |
                +-----------+------------+-----------------+----------------------+
                |   .med    |    R       | SALOME [#f8]_   | med, salome          |
                +-----------+------------+-----------------+----------------------+
                |   .obj    |    R       | WAVEFRONT       | obj                  |
                +-----------+------------+-----------------+----------------------+

                By default, Meshmagick uses the filename extensions to choose the
                appropriate reader/writer. This behaviour might be bypassed using the
                -ifmt and -ofmt optional arguments. When using these options, keywords
                defined in the table above must be used as format identifiers.
                
                .. rubric:: Footnotes
                
                .. [#f1] NEMOH is an open source BEM Software for seakeeping developped at
                         Ecole Centrale de Nantes (LHHEA)
                .. [#f2] WAMIT is a BEM Software for seakeeping developped by WAMIT, Inc.
                .. [#f3] DIODORE is a BEM Software for seakeeping developped by PRINCIPIA
                .. [#f4] HYDROSTAR is a BEM Software for seakeeping developped by
                         BUREAU VERITAS
                .. [#f5] GMSH is an open source meshing software developped by C. Geuzaine
                         and J.-_faces. Remacle
                .. [#f6] PARAVIEW is an open source visualization software developped by
                         Kitware
                .. [#f7] TECPLOT is a visualization software developped by Tecplot
                .. [#f8] SALOME-MECA is an open source software for computational mechanics
                         developped by EDF-R&D


                """,
    epilog='--  Copyright 2014-%u  -  Francois Rongere  /  Ecole Centrale de Nantes  --' % __year__,
    formatter_class=argparse.RawDescriptionHelpFormatter)

# TODO: ajouter option pour voir l'ensemble des formats de fichier geres par meshmagick avec une explication du logiciel utilise

parser.add_argument('infilename',  # TODO : voir pour un typ=file pour tester l'existence
                    help='path of the input mesh file in any supported format')

parser.add_argument('-o', '--outfilename', type=str,
                    help="""path of the output mesh file. The format of
                     this file is determined from the given extension.
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

parser.add_argument('-t', '--translate', metavar=('Tx', 'Ty', 'Tz'),
                    nargs=3, type=float,
                    help="""translates the mesh in 3D
                    Usage -translate tx ty tz""")

parser.add_argument('-tx', '--translatex', type=float, metavar='Tx',
                    help="""translates the mesh following the x direction""")

parser.add_argument('-ty', '--translatey', type=float, metavar='Ty',
                    help="""translates the mesh following the y direction""")

parser.add_argument('-tz', '--translatez', type=float, metavar='Tz',
                    help="""translates the mesh following the z direction""")

parser.add_argument('-r', '--rotate', metavar=('Rx', 'Ry', 'Rz'),
                    nargs=3, type=float,
                    help="""rotates the mesh in 3D following a rotation
                    coordinate vector. It is done around fixed axis. Angles
                    must be given in degrees.""")

parser.add_argument('-rx', '--rotatex', type=float, metavar='Rx',
                    help="""rotates the mesh around the x direction.
                    Angles must be given in degrees.""")

parser.add_argument('-ry', '--rotatey', type=float, metavar='Ry',
                    help="""rotates the mesh around the y direction.
                    Angles must be given in degrees.""")

parser.add_argument('-rz', '--rotatez', type=float, metavar='Rz',
                    help="""rotates the mesh around the z direction.
                    Angles must be given in degrees.""")

parser.add_argument('-s', '--scale', type=float, metavar='S',
                    help="""scales the mesh. CAUTION : if used along
                     with a translation option, the scaling is done after
                    the translations. The translation magnitude should be set
                    accordingly to the original mesh.
                    """)

parser.add_argument('-sx', '--scalex', type=float, metavar='Sx',
                    help="""scales the mesh along x axis. CAUTION : if used along
                     with a translation option, the scaling is done after
                    the translations. The translation magnitude should be set
                    accordingly to the original mesh.
                    """)

parser.add_argument('-sy', '--scaley', type=float, metavar='Sy',
                    help="""scales the mesh along y axis. CAUTION : if used along
                     with a translation option, the scaling is done after
                    the translations. The translation magnitude should be set
                    accordingly to the original mesh.
                    """)

parser.add_argument('-sz', '--scalez', type=float, metavar='Sz',
                    help="""scales the mesh along z axis. CAUTION : if used along
                     with a translation option, the scaling is done after
                    the translations. The translation magnitude should be set
                    accordingly to the original mesh.
                    """)

parser.add_argument('-hn', '--heal-normals', action='store_true',
                    help="""Checks and heals the normals consistency and
                    verify if they are outward.
                    """)

parser.add_argument('-fn', '--flip-normals', action='store_true',
                    help="""flips the normals of the mesh""")

parser.add_argument('-hm', '--heal-mesh', action='store_true',
                    help="""Applies the following sanity transformation on the
                    mesh: Removes unused vertices, Removes degenerated faces,
                    Merge duplicate vertices, Heal triangles description,
                    Heal normal orientations.
                    """)

parser.add_argument('-p', '--plane', nargs='+', action='append', metavar='Arg',
                    help="""Defines a plane used by the --clip_by_plane and --symmetrize options.
                    It can be defined by the floats nx ny nz c where [nx, ny, nz]
                    is a normal vector to the plane and c defines its position
                    following the equation <N|X> = c with X a point belonging
                    to the plane.
                    It can also be defined by a string among [Oxy, Oxz, Oyz, /Oxy, /Oxz, /Oyz]
                    for quick definition. Several planes may be defined on the same command
                    line. Planes with a prepended '/' have normals inverted i.e. if Oxy has its
                    normal pointing upwards, /Oxy has its normal pointing downwards.
                    In that case, the planes are indexed by an integer starting by
                    0 following the order given in the command line.
                    """)

parser.add_argument('-c', '--clip-by-plane', nargs='*', action='append', metavar='Arg',
                    help="""cuts the mesh with a plane. Is no arguments are given, the Oxy plane
                    is used. If an integer is given, it should correspond to a plane defined with
                    the --plane option. If a key string is given, it should be a valid key (see
                    help of --plane option for valid plane keys). A normal and a scalar could
                    also be given for the plane definition just as for the --plane option. Several
                    clipping planes may be defined on the same command line.""")

parser.add_argument('-md', '--merge-duplicates', nargs='?', const='1e-8',
                    default=None, metavar='Tol',
                    help="""merges the duplicate nodes in the mesh with the absolute tolerance
                    given as argument (default 1e-8). Tolerance must be lower than 1""")

parser.add_argument('-tq', '--triangulate-quadrangles', action='store_true',
                    help="""Triangulate all quadrangle _faces by a simple splitting procedure.
                    Twho triangles are generated and from both solution, the one with the best
                    aspect ratios is kept. This option may be used in conjunction with a
                    mesh export in a format that only deal with triangular cells like STL format.""")

parser.add_argument('-sym', '--symmetrize', nargs='*', action='append', metavar='Arg',
                    help="""Symmetrize the mesh by a plane defined wether by 4 scalars, i.e.
                    the plane normal vector coordinates and a scalar c such as N.X=c is the
                    plane equation (with X a point of the plane) or a string among ['Oxy',
                    'Oxz', 'Oyz', '/Oxy', '/Oxz', '/Oyz'] which are shortcuts for planes
                    passing by the origin and whose normals are the reference axes. Default
                    is Oxz if no argument is given to --sym option.
                    Be careful that symmetry is applied before any rotation so as the plane
                    equation is defined in the initial frame of reference.""")

parser.add_argument('--mirror', nargs='+', metavar='Arg',
                    help="""Mirror the mesh through the specified plane. Plane may be specified
                    with reference planes keys (see --plane option), or by 4 scalars, or by the
                    id of a plane defined with the --plane option. By default, the Oxy plane
                    is used when the option has no argument.""")

# FIXME: on devrait pouvoir laisser les valeurs par defaut --> creer une option --rho-medium
parser.add_argument('-pi', '--plain-inertia', action='store_true',
                    help="""Evaluates the inertia properties of the mesh condidering it as
                    uniformly plain of a medium of density rho_medium in kg/m**3. Default
                    is 1023 kg/m**3.""")

# TODO: creer une option --thickness
parser.add_argument('-si', '--shell-inertia', action='store_true',
                    help="""Evaluates the inertia properties of the mesh condidering it as
                    uniformly plain of a medium of density rho_medium in kg/m**3. Default
                    is 1023 kg/m**3.""")

parser.add_argument('--rho-medium', type=float,
                    help="""The density (in kg/m**3) of the medium used for evaluation of
                    inertia parameters of the mesh. For the hypothesis of plain homogeneous
                    mesh, the default is that of salt water (1023 kg/m**3) . For the
                    hypothesis of a shell, default is that of steel (7850 kg/m**3).

                    It is possible to specify medium by a name. Available medium are
                    currently: %s
                    """ % list_medium())

parser.add_argument('--list-medium', action='store_true',
                    help="""Lists the available medium keywords along with their density.
                    """
                    )

parser.add_argument('--thickness', type=float,
                    help="""The thickness of the shell used for the evaluation of inertia
                    parameters of the mesh. The default value is 0.02m.""")


# Arguments for hydrostatics computations
# TODO: l'option -hs devrait etre une sous-commande au sens de argparse
# TODO: completer l'aide de -hs
parser.add_argument('-hs', '--hydrostatics', action='store_true',
                    help="""Compute hydrostatics data and throws a report on the
                    command line. When used along with options --disp, --cog or
                    --zcog, the behaviour is different.""")

# TODO: replace disp by mass as it is more correct...
parser.add_argument('--disp', default=None, type=float, metavar='Disp',
                    help="""Specifies the mass of the mesh for hydrostatic computations.
                        It must be given in tons.
                        """)

parser.add_argument('--cog', nargs=3, metavar=('Xg', 'Yg', 'Zg'),
                    help="""Specifies the 3D position of the center of gravity.
                    The third coordinate given has priority over the value given
                    with the --zcog option.""")

parser.add_argument('--zcog', default=None, type=float, metavar='Zcog',
                    help="""Specifies the z position of the center of gravity. This
                        is the minimal data needed for hydrostatic stiffness matrix
                        computation. It is however overwriten by the third component
                        of cog when --cog option is used. If none of these two option
                        is given, zcog is set to 0.
                        """)

parser.add_argument('--rho-water', default=1023., type=float, metavar='Rho',
                    help="""Specifies the density of salt water. Default is 1023 kg/m**3.
                    """)

parser.add_argument('-g', '--grav', default=9.81, type=float, metavar='G',
                    help="""Specifies the acceleration of gravity on the earth surface.
                    Default is 9.81 m/s**2.
                    """)

# parser.add_argument('--hs_solver_params', nargs='+')

parser.add_argument('-af', '--absolute-force', nargs=6, action='append',
                    metavar=('X', 'Y', 'Z', 'Fx', 'Fy', 'Fz'),
                    help="""Add an additional absolute force applied to the mesh in
                    hydrostatics equilibrium computations. It is absolute as force
                    orientation does not change during hydrostatic equilibrium
                    computations, but application point follows the mesh motion.
                    The option is called with 6 arguments such that:
                    
                        -af x y z fx fy fz
                     
                    where [x, y, z] are the coordinates of the application point
                    (in meters) and [fx, fy, fz] are the components of the force.
                    Those are expressed in the initial mesh frame.
                     """)

parser.add_argument('-rf', '--relative-force', nargs=6, action='append', type=float,
                    metavar=('X', 'Y', 'Z', 'Fx', 'Fy', 'Fz'),
                    help="""Add an additional relative force applied to the mesh in
                    hydrostatics equilibrium computations. It is relative as force
                    orientation follows the orientation of the mesh during equilibrium
                    computations. The option is called with 6 arguments such that:
                     
                        -af x y z fx fy fz
                     
                    where [x, y, z] are the coordinates of the application point
                    of the force and [fx, fy, fz] are the components of the force.
                    Those are expressed in the initial mesh frame.
                     """)

parser.add_argument('--hs-report', type=str, metavar='filename',
                    help="""Write the hydrostatic report into the file given as an argument""")

# ARGUMENTS RELATED TO THE COMPUTATION OF INERTIA PARAMETERS
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

# parser.add_argument('--lid', nargs='?', const=1., default=None, type=float,
#                     help="""Generate a triangle mesh lid on the mesh clipped by the Oxy plane.
#                     """)

# parser.add_argument('--fill-holes', '-fh', action='store_true',
#                     help="""Fill little holes by triangulation if any.
#                     """)

parser.add_argument('--show', action='store_true',
                    help="""Shows the input mesh in an interactive window""")

parser.add_argument('--version', action='version',
                    version='meshmagick - version %s\n%s' % (__version__, __copyright__),
                    help="""Shows the version number and exit""")


def main():
    if acok:
        argcomplete.autocomplete(parser)

    # TODO : Utiliser des sous-commandes pour l'utilisation de meshmagick

    args, unknown = parser.parse_known_args()

    if args.quiet:
        verbose = False
    else:
        verbose = True

    if verbose:
        print('\n=============================================')
        print(('meshmagick - version %s\n%s' % (__version__, __copyright__)))
        print('=============================================')

    # LOADING DATA FROM FILE
    if args.input_format is not None:
        format = args.input_format
    else:
        # Format based on extension
        _, ext = os.path.splitext(args.infilename)
        format = ext[1:].lower()
        if format == '':
            raise IOError('Unable to determine the input file format from its extension. Please specify an input format.')

    # Loading mesh elements from file
    if os.path.isfile(args.infilename):
        V, F = mmio.load_mesh(args.infilename, format)

        # Give the name of the mesh the filename
        basename = os.path.basename(args.infilename)
        mesh_name, _ = os.path.splitext(basename)

        mesh = Mesh(V, F, name=mesh_name)
        # Ensuring triangles are following the right convention (last id = first id)
        mesh.heal_triangles()
        if verbose:
            mesh.verbose_on()
            print(('%s successfully loaded' % args.infilename))
    else:
        raise IOError('file %s not found' % args.infilename)

    # Merge duplicate _vertices
    if args.merge_duplicates is not None:
        tol = float(args.merge_duplicates)
        mesh.merge_duplicates(atol=tol)

    # TODO : put that dict at the begining of the main function or in the module
    plane_str_list = {'Oxy': [0., 0., 1.],
                      'Oxz': [0., 1., 0.],
                      'Oyz': [1., 0., 0.],
                      '/Oxy': [0., 0., -1.],
                      '/Oxz': [0., -1., 0.],
                      '/Oyz': [-1., 0., 0.]}

    if args.quality:
        mesh.print_quality()

    # Defining planes
    planes = []
    if args.plane is not None:
        nb_planes = len(args.plane)

        if verbose:
            if nb_planes == 1:
                verb = 'plane has'
            else:
                verb = 'planes have'
            print(('\n%u %s been defined:' % (nb_planes, verb)))
            # TODO: ajouter un recapitulatif des plans definis

        planes = [Plane() for i in range(nb_planes)]
        for (iplane, plane) in enumerate(args.plane):
            if len(plane) == 4:
                # plane is defined by normal and scalar
                try:
                    planes[iplane] = Plane(normal=list(map(float, plane[:3])), scalar=plane[3])
                except:
                    raise AssertionError('Defining a plane by normal and scalar requires four scalars')

            elif len(plane) == 1:
                if plane[0] in plane_str_list:
                    planes[iplane].normal = np.array(plane_str_list[plane[0]], dtype=np.float)
                    planes[iplane].c = 0.
                else:
                    raise AssertionError('%s key for plane is not known. Choices are [%s].'
                                         % (plane[0], ', '.join(list(plane_str_list.keys()))))
            else:
                raise AssertionError('Planes should be defined by a normal and a scalar '
                                     'or by a key to choose among [%s]' % (', '.join(list(plane_str_list.keys()))))
        if verbose:
            for plane_id, plane in enumerate(planes):
                print(("\t%u: %s" % (plane_id, plane)))

    # Mirroring the mesh
    if args.mirror is not None:
        sym_plane = Plane()
        print((args.mirror))

        if len(args.mirror) == 1:
            # May be a standard plane or a plane id
            if len(planes) > 0:
                try:
                    plane_id = int(args.mirror[0])
                    if plane_id >= 0 and plane_id < len(planes):
                        sym_plane = planes[plane_id]
                    else:
                        raise AssertionError('Only planes IDs from 0 to %u have been defined. %u is outside the range.'
                                             % (len(planes) - 1, plane_id))
                except ValueError:
                    # Cannot be converted to an int, it should be the key of a plane
                    try:
                        sym_plane.normal = plane_str_list[args.mirror[0]]
                    except KeyError as err:
                        raise KeyError('%s is not a standard plane identifier' % err)
            else:
                try:
                    sym_plane.normal = plane_str_list[args.mirror[0]]
                except KeyError as err:
                    raise KeyError('%s is not a standard plane identifier' % err)

        elif len(args.mirror) == 4:
            sym_plane.normal = args.mirror[:3]
            sym_plane.c = args.mirror[3]
        else:
            raise ValueError('Bad plane definition.')

        if verbose:
            print(('Mirroring the mesh by :\n\t%s' % sym_plane))

        mesh.mirror(sym_plane)
        if verbose:
            print('\t-> Done.')

    # Symmetrizing the mesh
    if args.symmetrize is not None:
        nb_sym = len(args.symmetrize)
        for iplane, plane in enumerate(args.symmetrize):
            if len(plane) == 0:
                args.symmetrize[iplane] = ['Oxz']

        if verbose:
            if nb_sym == 1:
                verb = 'plane'
            else:
                verb = 'planes'
            print(('\nMesh is being symmetrized by %u %s:' % (nb_sym, verb)))

        for plane in args.symmetrize:
            sym_plane = Plane()
            if len(plane) == 1:
                if len(planes) > 0:
                    try:
                        plane_id = int(plane[0])
                        if plane_id >= 0 and plane_id < len(planes):
                            sym_plane = planes[plane_id]
                        else:
                            raise AssertionError(
                                'Only planes IDs from 0 to %u have been defined. %u is outside the range.' % (
                                    len(planes) - 1, plane_id))
                    except ValueError:
                        try:
                            sym_plane.normal = plane_str_list[plane[0]]
                        except KeyError as err:
                            raise KeyError('%s is not a standard plane identifier' % err)
                else:
                    try:
                        sym_plane.normal = plane_str_list[plane[0]]
                    except KeyError as err:
                        raise KeyError('%s is not a standard plane identifier' % err)

            elif len(plane) == 4:
                sym_plane.normal = plane[:3]
                sym_plane.c = plane[3]
            else:
                raise ValueError('Bad plane definition.')

            if verbose:
                print(('\t%s' % sym_plane))
            mesh.symmetrize(sym_plane)

        if verbose:
            print('\t-> Done.')

    # Globally heal the mesh
    if args.heal_mesh:
        if verbose:
            print('\nOPERATION: heal the mesh')
        mesh.heal_mesh()
        if verbose:
            print('\tDone.')

    # Heal normals
    if args.heal_normals and not args.heal_mesh:
        if verbose:
            print('\nOPERATION: heal normals')
        mesh.heal_normals()
        if verbose:
            print('\t-> Done.')

    # Mesh translations
    if args.translate is not None:
        if verbose:
            print(('\nOPERATION: Translation by [%f, %f, %f]' % tuple(args.translate)))
        mesh.translate(args.translate)
        if verbose:
            print('\t-> Done.')

    if args.translatex is not None:
        if verbose:
            print(('\nOPERATION: Translation by %f along X' % args.translatex))
        mesh.translate_x(args.translatex)
        if verbose:
            print('\t-> Done.')

    if args.translatey is not None:
        if verbose:
            print(('\nOPERATION: Translation by %f along Y' % args.translatey))
        mesh.translate_y(args.translatey)
        if verbose:
            print('\t-> Done.')

    if args.translatez is not None:
        if verbose:
            print(('\nOPERATION: Translation by %f along Z' % args.translatez))
        mesh.translate_z(args.translatez)
        if verbose:
            print('\t-> Done.')

    # Mesh rotations
    if args.rotate is not None:
        if verbose:
            print(('\nOPERATION: Rotation by [%f, %f, %f] (degrees)' % tuple(args.rotate)))
        mesh.rotate(list(map(math.radians, args.rotate)))
        if verbose:
            print('\t-> Done.')

    if args.rotatex is not None:
        if verbose:
            print(('\nOPERATION: Rotation by %f around X (Roll)' % args.rotatex))
        mesh.rotate_x(math.radians(args.rotatex))
        if verbose:
            print('\t-> Done.')

    if args.rotatey is not None:
        if verbose:
            print(('\nOPERATION: Rotation by %f around Y (Pitch)' % args.rotatey))
        mesh.rotate_y(math.radians(args.rotatey))
        if verbose:
            print('\t-> Done.')

    if args.rotatez is not None:
        if verbose:
            print(('\nOPERATION: Rotation by %f around Z (Yaw)' % args.rotatez))
        mesh.rotate_z(math.radians(args.rotatez))
        if verbose:
            print('\t-> Done.')

    if args.scale is not None:
        if verbose:
            print(('\nOPERATION: Scaling by %f' % args.scale))
        mesh.scale(args.scale)
        if verbose:
            print('\t-> Done.')

    if args.scalex is not None:
        if verbose:
            print(('\nOPERATION: Scaling by %f along the x axis' % args.scalex))
        mesh.scalex(args.scalex)
        if verbose:
            print('\t-> Done.')

    if args.scaley is not None:
        if verbose:
            print(('\nOPERATION: Scaling by %f along the y axis' % args.scaley))
        mesh.scaley(args.scaley)
        if verbose:
            print('\t-> Done.')

    if args.scalez is not None:
        if verbose:
            print(('\nOPERATION: Scaling by %f along the z axis' % args.scalez))
        mesh.scalez(args.scalez)
        if verbose:
            print('\t-> Done.')

    if args.flip_normals:
        if verbose:
            print('\nOPERATION: Flipping normals')
        mesh.flip_normals()
        if verbose:
            print('\t-> Done.')

    if args.triangulate_quadrangles:
        mesh.triangulate_quadrangles()

    # Clipping the mesh
    if args.clip_by_plane is not None:
        nb_clip = len(args.clip_by_plane)
        for iplane, plane in enumerate(args.clip_by_plane):
            if len(plane) == 0:
                args.clip_by_plane[iplane] = ['Oxy']

        if verbose:
            if nb_clip == 1:
                verb = 'plane'
            else:
                verb = 'planes'
            print(('\nMesh is being clipped by %u %s' % (nb_clip, verb)))

        for plane in args.clip_by_plane:
            clipping_plane = Plane()
            if len(plane) == 1:
                if len(planes) > 0:
                    try:
                        plane_id = int(plane[0])
                        if plane_id >= 0 and plane_id < len(planes):
                            clipping_plane = planes[plane_id]
                        else:
                            raise AssertionError(
                                'Only planes IDs from 0 to %u have been defined. %u is outside the range.' % (
                                    len(planes) - 1, plane_id))
                    except ValueError:
                        try:
                            clipping_plane.normal = plane_str_list[plane[0]]
                        except KeyError as err:
                            raise KeyError('%s is not a standard plane identifier' % err)
                else:
                    try:
                        clipping_plane.normal = plane_str_list[plane[0]]
                    except KeyError as err:
                        raise KeyError('%s is not a standard plane identifier' % err)

            elif len(plane) == 4:
                clipping_plane.normal = plane[:3]
                clipping_plane.c = plane[3]
            else:
                raise ValueError('Bad plane definition.')

            if verbose:
                print(('\t%s' % clipping_plane))

            clipper = MeshClipper(mesh, plane=clipping_plane)
            mesh = clipper.clipped_mesh

        if verbose:
            print('\t-> Done.')

    # Listing available medium
    if args.list_medium:
        col_width = 22
        hline = '+{0:s}+{0:s}+\n'.format('-' * col_width)
        table = '\n' + hline
        table += '|{:<{n}s}|{:>{n}s}|\n'.format('NAME', 'DENSITY (KG/M**3)', n=col_width)
        table += hline
        for medium in densities.list_medium():
            table += '|{:<{n}s}|{:>{n}.3f}|\n'.format(medium, densities.get_density(medium), n=col_width)
            table += hline
        print(table)

    # Calculate the plain inertia
    if args.plain_inertia:
        if args.rho_medium is None:
            rho_medium = 1023.
        else:
            rho_medium = args.rho_medium

        inertia = mesh.eval_plain_mesh_inertias(rho_medium=rho_medium)

        if verbose:
            print(("\nInertial parameters for a uniform distribution of a medium of density %.1f kg/m**3 in the mesh:\n" \
                  % rho_medium))
            print(("\tMass = %.3f tons" % (inertia.mass / 1000.)))
            cog = inertia.gravity_center
            print(("\tCOG (m):\n\t\txg = %.3f\n\t\tyg = %.3f\n\t\tzg = %.3f" % (cog[0], cog[1], cog[2])))
            mat = inertia.inertia_matrix
            print("\tInertia matrix (SI):")
            print(("\t\t%.3E\t%.3E\t%.3E" % (mat[0, 0], mat[0, 1], mat[0, 2])))
            print(("\t\t%.3E\t%.3E\t%.3E" % (mat[1, 0], mat[1, 1], mat[1, 2])))
            print(("\t\t%.3E\t%.3E\t%.3E" % (mat[2, 0], mat[2, 1], mat[2, 2])))
            point = inertia.reduction_point
            print(("\tExpressed at point : \t\t%.3E\t%.3E\t%.3E" % (point[0], point[1], point[2])))

    if args.shell_inertia:
        if args.rho_medium is None:
            rho_medium = 7850.
        else:
            rho_medium = args.rho_medium
        if args.thickness is None:
            thickness = 0.02
        else:
            thickness = args.thickness

        inertia = mesh.eval_shell_mesh_inertias(rho_medium=rho_medium, thickness=thickness)

        if verbose:
            print(("\nInertial parameters for a shell distribution of a medium of density %.1f kg/m**3 and a thickness " \
                  "of %.3f m over the mesh:\n" % (rho_medium, thickness)))
            print(("\tMass = %.3f tons" % (inertia.mass / 1000.)))
            cog = inertia.gravity_center
            print(("\tCOG (m):\n\t\txg = %.3f\n\t\tyg = %.3f\n\t\tzg = %.3f" % (cog[0], cog[1], cog[2])))
            mat = inertia.inertia_matrix
            print("\tInertia matrix (SI):")
            print(("\t\t%.3E\t%.3E\t%.3E" % (mat[0, 0], mat[0, 1], mat[0, 2])))
            print(("\t\t%.3E\t%.3E\t%.3E" % (mat[1, 0], mat[1, 1], mat[1, 2])))
            print(("\t\t%.3E\t%.3E\t%.3E" % (mat[2, 0], mat[2, 1], mat[2, 2])))
            point = inertia.reduction_point
            print(("\tExpressed at point : \t\t%.3E\t%.3E\t%.3E" % (point[0], point[1], point[2])))

    additional_forces = []
    if args.relative_force is not None:
        for item in args.relative_force:
            force = hs.Force(point=list(map(float, item[:3])), value=list(map(float, item[3:])), mode='relative')
            additional_forces.append(force)

    if args.absolute_force is not None:
        for item in args.absolute_force:
            force = hs.Force(point=list(map(float, item[:3])), value=list(map(float, item[3:])), mode='absolute')
            additional_forces.append(force)

    if args.hydrostatics:
        grav = args.grav
        rho_water = args.rho_water

        hs_solver = hs.Hydrostatics(mesh, rho_water=rho_water, grav=grav, verbose=verbose)

        for force in additional_forces:
            hs_solver.add_force(force)

        if args.disp is not None:
            disp = args.disp
            has_disp = True
        else:
            disp = hs_solver.displacement
            has_disp = False

        if args.cog is not None:
            cog = list(map(float, args.cog))
            has_cog = True
        else:
            cog = hs_solver.gravity_center
            has_cog = False

        if args.zcog is not None:
            zcog = args.zcog
            has_zcog = True
        else:
            zcog = hs_solver.zg
            has_zcog = False

        # Overwrite zcog by cog[2] in case both have been given
        if has_cog and has_zcog:
            zcog = cog[2]

        case = (has_disp, has_cog, has_zcog)

        if case == (True, False, False) or case == (True, False, True):
            if verbose:
                print(("\nSetting mesh displacement to %.3f tons" % disp))
            hs_solver.zg = zcog
            hs_solver.set_displacement(disp)

        msg = """\nWARNING !! : Keep in mind that the mesh may have a new orientation in the Oxy plane so the hydrostatic
        stiffness matrix may be impractical as not expressed in a convenient axis system. Use the result
        with caution and have a look at the mesh by using the --show option.
        
        This should be fixed in a future release."""

        if case == (False, True, False) or case == (False, True, True):
            if verbose:
                print(("""
                \nSetting mesh position so that we are at the current displacement of %.3f tons and in a stable
                configuration with a gravity center located at [%.3f, %.3f, %.3f] meters.""" \
                      % (disp, cog[0], cog[1], cog[2])))

            hs_solver.gravity_center = cog
            hs_solver.mass = disp
            hs_solver.equilibrate(init_disp=False)
            warn(msg)

        if case == (False, False, True) or case == (False, False, False):
            if verbose:
                print(("\nGenerating hydrostatic data for a zcog of %.3f meters." % zcog))
            hs_solver.zg = zcog

        if case == (True, True, False) or case == (True, True, True):
            if verbose:
                print(("""
                \nSetting mesh position so that the displacement is equal to %.3f tons and in a stable
                configuration with a gravity center located at [%.3f, %.3f, %.3f] meters.""" \
                      % (disp, cog[0], cog[1], cog[2])))
            hs_solver.gravity_center = cog
            hs_solver.mass = disp
            hs_solver.equilibrate(init_disp=True)
            warn(msg)

        # TODO: voir pour une option pour sortir plutot le maillage coupe pour Nemoh
        mesh = hs_solver

        if verbose:
            print((hs_solver.get_hydrostatic_report()))

        if args.hs_report is not None:
            with open(args.hs_report, 'w') as f:
                f.write('==============================================\n')
                f.write('Hydrostatic report generated by Meshmagick\n')
                f.write('Meshfile: %s\n' % os.path.abspath(args.infilename))
                f.write('%s\n' % strftime('%c'))
                f.write('meshmagick - version %s\n%s\n' % (__version__, __copyright__))
                f.write('==============================================\n')
                f.write(hs_solver.get_hydrostatic_report())



    # WARNING : No more mesh modification should be released from this point until the end of the main

    if args.info:
        print(mesh)

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
                        raise IOError('Could not determine a format from input file extension, please specify an input format or an extension')
            else:
                format = os.path.splitext(args.outfilename)[1][1:].lower()

        if verbose:
            print(('Writing %s' % args.outfilename))
        mmio.write_mesh(args.outfilename, mesh.vertices, mesh.faces, format)
        if verbose:
            print('\t-> Done.')

    if verbose:
        print('\n=============================================================')
        print(('meshmagick - version %s\n%s' % (__version__, __copyright__)))
        print(('Maintainer : %s <%s>' % (__maintainer__, __email__)))
        print('Good Bye!')
        print('=============================================================')


if __name__ == '__main__':
    main()
