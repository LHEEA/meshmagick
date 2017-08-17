#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

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
            print('\t-> %u %s been detected' % (nb_hole, word))

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
        print("\n\t-> %u %s been added successfully\n" % (nb_bound, verb))

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

