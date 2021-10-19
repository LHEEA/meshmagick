#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PYTHON_ARGCOMPLETE_OK

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

""" Pour activer l'aucompletion meshmagick:
installer argcomplete (conda install argcomplete)
Activer la completion globale:
    activate-global-python-argcomplete --dest=/home/<username>/
Ca installe un fichier /home/<username>/python-argcomplete
qu'il faut sourcer dans le .bashrc soit ajotuer la ligne suivante dans le bashrc:
    source /home/<username>/python-argcomplete
L'autocompletion devrait etre maintenant active pour meshmagick
"""

# TODO: move meshmagick.py at the root level of the project ?

import os, sys
# import numpy as np
# import math
from datetime import datetime
# from warnings import warn
from time import strftime
import argparse

from meshmagick.mesh import *
from meshmagick import mmio
from meshmagick.mesh_clipper import MeshClipper
from meshmagick import hydrostatics as hs

from meshmagick import densities
from meshmagick import __version__

__year__ = datetime.now().year

__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-%u, Ecole Centrale de Nantes / D-ICE Engineering" % __year__
__credits__ = "Francois Rongere"
__licence__ = "GPLv3"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@dice-engineering.com"

__all__ = ['main']


def list_medium():
    return ', '.join(densities.list_medium())


# =======================================================================
#                         COMMAND LINE USAGE
# =======================================================================


try:
    import argcomplete
    has_argcomplete = True
except:
    has_argcomplete = False

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

parser.add_argument('-cc', '--concatenate-file', type=str,
                    help="""Concatenate a mesh from the specified path. The file format has to be
                    the same as the input file.""")

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

parser.add_argument('-pn', '--project-name', default="NO_NAME", type=str, metavar='Project Name',
                    help="""The project name for hydrostatics ourput
                    """)

parser.add_argument('-hs', '--hydrostatics', action='store_true',
                    help="""Compute hydrostatics data and throws a report on the
                    command line. When used along with options -mdisp, --cog or
                    --zcog, the behaviour is different.""")

# TODO: replace disp by mass as it is more correct...
parser.add_argument('-mdisp', '--mass-displacement', default=None, type=float, metavar='Disp',
                    help="""Specifies the mass of the mesh for hydrostatic computations.
                        It MUST be given in tons.
                        """)

parser.add_argument('-cog', '--cog', nargs=3, metavar=('Xg', 'Yg', 'Zg'),
                    help="""Specifies the 3D position of the center of gravity.
                    The third coordinate given has priority over the value given
                    with the --zcog option.""")

parser.add_argument('-zg', '--zcog', default=None, type=float, metavar='Zcog',
                    help="""Specifies the z position of the center of gravity. This
                        is the minimal data needed for hydrostatic stiffness matrix
                        computation. It is however overwriten by the third component
                        of cog when --cog option is used. If none of these two option
                        is given, zcog is set to 0.
                        """)

parser.add_argument('-lpp', '--lpp', default=None, type=float, metavar='Lpp',
                    help="""Specifies the Lpp value as it cannot be calculated with
                        only the mesh as it depends on the AP position that is a 
                        rudder position dependent information that the mesh does not
                        enclose. It helps do better inertia (Iyy & Izz) approximations
                        using standard formulas. See also the alternative -AP formula.  
                        """)

parser.add_argument('-ap', '--orig-at-ap', action='store_true',
                    help="""Tell the solver that the origin is ar After perpendicular 
                    sot that lpp can be calculated from this information.""")

parser.add_argument('-wd', '--water-density', default=1025., type=float, metavar='Rho',
                    help="""Specifies the density of salt water. Default is 1025 kg/m**3.
                    """)

parser.add_argument('-g', '--grav', default=9.81, type=float, metavar='G',
                    help="""Specifies the acceleration of gravity on the earth surface.
                    Default is 9.81 m/s**2.
                    """)

# parser.add_argument('--hs_solver_params', nargs='+')


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

parser.add_argument('-sh', '--show', action='store_true',
                    help="""Shows the input mesh in an interactive window""")

parser.add_argument('-v', '--version', action='version',
                    version='meshmagick - version %s\n%s' % (__version__, __copyright__),
                    help="""Shows the version number and exit""")


def main():
    if has_argcomplete:
        argcomplete.autocomplete(parser)

    # Parse command line arguments
    args = parser.parse_args()

    verbose = True
    if args.quiet:
        verbose = False

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
            raise IOError(
                'Unable to determine the input file format from its extension. Please specify an input format.')

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

    if args.concatenate_file is not None:
        print('Concatenate %s with %s' % (args.infilename, args.concatenate_file))
        # Loading the file
        if os.path.isfile(args.concatenate_file):
            Vc, Fc = mmio.load_mesh(args.concatenate_file, format)

            # Give the name of the mesh the filename
            basename = os.path.basename(args.concatenate_file)
            mesh_name, _ = os.path.splitext(basename)

            mesh_c = Mesh(Vc, Fc, name=mesh_name)
            # Ensuring triangles are following the right convention (last id = first id)
            mesh_c.heal_triangles()
            if verbose:
                mesh_c.verbose_on()
                print(('%s successfully loaded' % args.concatenate_file))
        else:
            raise IOError('file %s not found' % args.concatenate_file)
        mesh += mesh_c

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


    has_disp = has_cog = has_zcog = False

    if args.mass_displacement is not None:
        disp = args.mass_displacement
        has_disp = True

    cog = np.zeros(3)
    if args.cog is not None:
        cog = list(map(float, args.cog))
        has_cog = True

    if args.zcog is not None:
        zcog = args.zcog
        has_zcog = True

    orig_at_ap = False
    if args.orig_at_ap:
        orig_at_ap = True

    # parameters
    water_density = args.water_density
    grav = args.grav
    lpp = args.lpp

    project_name = args.project_name


    # FIXME: il faut prendre en compte dans les z_corr et rotmat_corr les informations de transformation donnees
    # dans les operations de translation et de rotation precedemment !!!

    z_corr = 0.
    rotmat_corr = np.eye(3, 3)

    if args.hydrostatics:

        reltol = 1e-6
        z_corr = 0.
        rotmat_corr = np.eye(3, 3)
        hs_data = dict()

        if not has_disp and not has_cog:
            print(">>>> Performing hydrostatic computation on the current hull configuration considered at equilibrium")
            if not has_zcog:
                raise RuntimeError("zcog should at least be given for correct stiffness values computations")

            hs_data = hs.compute_hydrostatics(mesh, np.zeros(3), water_density, grav)
            xb, yb, _ = hs_data["buoyancy_center"]
            cog = np.array([xb, yb, zcog])
            hs_data = hs.compute_hydrostatics(mesh, cog, water_density, grav,
                                              at_cog=True, lpp=lpp, orig_at_ap=orig_at_ap)

        elif has_disp and not has_cog:
            print(">>>> Computing equilibrium of the hull for the given displacement of %f tons" % disp)
            if not has_zcog:
                raise RuntimeError("zcog should at least be given for correct stiffness values computations")

            z_corr = hs.displacement_equilibrium(mesh, disp, water_density, grav,
                                                 reltol=reltol, verbose=True)
            hs_data = hs.compute_hydrostatics(mesh, np.zeros(3), water_density, grav,
                                              z_corr=z_corr, at_cog=False)
            xb, yb, _ = hs_data["buoyancy_center"]
            cog = np.array([xb, yb, zcog])
            hs_data = hs.compute_hydrostatics(mesh, cog, water_density, grav,
                                              z_corr=z_corr, at_cog=True, lpp=lpp, orig_at_ap=orig_at_ap)

        elif has_disp and has_cog:
            print(">>>> Computing equilibrium in 3DOF for the given displacement and COG")
            if has_zcog:
                warn("zcog is redundant with cog, taking cog and ignoring zcog")
            z_corr, rotmat_corr = hs.full_equilibrium(mesh, cog, disp, water_density, grav, reltol=reltol, verbose=True)
            hs_data = hs.compute_hydrostatics(mesh, cog, water_density, grav,
                                              z_corr=z_corr, rotmat_corr=rotmat_corr, at_cog=True, lpp=lpp,
                                              orig_at_ap=orig_at_ap)

        elif not has_disp and has_cog:
            print(
                ">>>> Computing equilibrium in 3DOF for the given COG, considering the current configuration presents the "
                "target displacement")
            if has_zcog:
                warn("zcog is redundant with cog, taking cog and ignoring zcog")

            hs_data = hs.compute_hydrostatics(mesh, np.zeros(3), water_density, grav,
                                              at_cog=False)
            disp = hs_data['disp_mass'] / 1000
            z_corr, rotmat_corr = hs.full_equilibrium(mesh, cog, disp, water_density, grav,
                                                      reltol=reltol, verbose=True)
            hs_data = hs.compute_hydrostatics(mesh, cog, water_density, grav,
                                              z_corr=z_corr, rotmat_corr=rotmat_corr, at_cog=True, lpp=lpp,
                                              orig_at_ap=orig_at_ap)

        mesh.rotate_matrix(rotmat_corr)
        mesh.translate_z(z_corr)
        # cog = np.dot(rotmat_corr, cog)
        # cog[2] += z_corr
        has_cog = True
        has_disp = True
        has_zcog = False

        hs_report = hs.get_hydrostatic_report(hs_data)
        print(hs_report)

        if args.hs_report is not None:
            with open(args.hs_report, 'w') as f:
                f.write('==============================================\n')
                f.write('Hydrostatic report generated by Meshmagick\n')
                f.write('Meshfile: %s\n' % os.path.abspath(args.infilename))
                f.write('%s\n' % strftime('%c'))
                f.write('meshmagick - version %s\n%s\n' % (__version__, __copyright__))
                f.write('==============================================\n')
                f.write(hs_report)

    
    # ==================================================================================================================
    # WARNING : No more mesh modification should be released below this point and until the end of the main function
    # ==================================================================================================================

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
                        raise IOError('Could not determine a format from input file extension, '
                                      'please specify an input format or an extension')
            else:
                format = os.path.splitext(args.outfilename)[1][1:].lower()

        if verbose:
            print(('Writing %s' % args.outfilename))
        mmio.write_mesh(args.outfilename, mesh.vertices, mesh.faces, format)
        if verbose:
            print('\t-> Done.')

    if verbose:
        print('\n=============================================================')
        print(('Meshmagick - version %s\n%s' % (__version__, __copyright__)))
        print(('Maintainer : %s <%s>' % (__maintainer__, __email__)))
        print('Good Bye!')
        print('=============================================================')


if __name__ == '__main__':
    main()
