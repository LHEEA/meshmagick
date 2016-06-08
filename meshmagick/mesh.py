#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
This module is part of meshmagick software. It provides functions for mesh clipping purpose.
"""

import meshmagick as mm # TODO: voir si oblige d'importer meshmagick ici ...
import MMviewer
import numpy as np
import math
import copy
import vtk
from itertools import count
import sys # Retirer

# TODO: voir si on ne peut pas mettre ces fonctions dans un module dedie ?
def _rodrigues(thetax, thetay):
    """
    Computes the rotation matrix corresponding to angles thetax and thetay using the Olinde-Rodrigues formula

    Parameters
    ----------
    thetax : float
        Angle around Ox axe (rad)
    thetay
        Angle around Oy axe (rad)
    Returns
    -------
    R : ndarray
        Rotation matrix
    """

    theta = math.sqrt(thetax*thetax + thetay*thetay)
    if theta == 0.:
        nx = ny = 0.
        ctheta = 1.
        stheta = 0.
    else:
        nx, ny = thetax/theta, thetay/theta
        ctheta = math.cos(theta)
        stheta = math.sin(theta)
    nxny = nx*ny

    # Olinde Rodrigues formulae
    # FIXME: S'assurer qu'on a effectivement pas de Ctheta devant le I3 !! et repercuter sur l'hydrostatique
    R = ctheta*np.eye(3) \
       + (1-ctheta) * np.array([[nx*nx, nxny, 0.],
                                [nxny, ny*ny, 0.],
                                [0., 0., 0.]]) \
       + stheta * np.array([[0., 0.,  ny],
                            [0., 0., -nx],
                            [-ny, nx, 0.]])
    return R

def _cardan(phi, theta):
    """
    Computes the rotation matrix corresponding to angles phi (roll) and theta (pitch) using Cardan angles convention

    Parameters
    ----------
    phi : float
        Roll angle (rad)
    theta : float
        Pitch angel (rad)

    Returns
    -------
    R : ndarray
        Rotation matrix
    """

    # Rotation matrix
    cphi = math.cos(phi)
    sphi = math.sin(phi)
    ctheta = math.cos(theta)
    stheta = math.sin(theta)

    Re0 = np.zeros((3, 3), dtype=float)
    Re0[0] = [ctheta, 0., -stheta]
    Re0[1] = [sphi*stheta, cphi, sphi*ctheta]
    Re0[2] = [cphi*stheta, -sphi, cphi*ctheta]

    return Re0

def _get_rotation_matrix(thetax, thetay, atype='fixed'):
    """
    Computes rotation matrix using different angle conventions

    Parameters
    ----------
    thetax : float
        Angle around x (rad)
    thetay : float
        Angle around y (rad)
    atype : {'fixed', 'cardan'}, optional
        Angle convention to use. Default to 'fixed' (fixed axes)

    Returns
    -------
    R : ndarray
        Rotation matrix

    """
    if atype == 'fixed':
        R = _rodrigues(thetax, thetay)
    elif atype == 'cardan':
        R = _cardan(thetax, thetay)
    else:
        raise AttributeError, 'Unknown angle convention: %s' % atype

    return R


# Classes
class Plane(object): # TODO: placer cette classe dans un module a part --> utilise dans meshmagick aussi...
    """
    Class to handle planes for clipping purposes
    """
    def __init__(self, normal=[0., 0., 1.], c=0.):
        """
        Plane constructor

        Parameters
        ----------
        normal : array_like, optional
            Normal of the plane. Default to [0., 0., 1.]
        c : float, optional
            Plane scalar parameter. Default to 0.

        Returns
        -------
            Plane object
        """

        normal = np.asarray(normal, dtype=np.float)

        self._normal = normal / np.linalg.norm(normal)
        self.c = float(c)

        # Storing rotation matrix (redundant !) to speedup computations
        # Shall be update in methods !!! --> using decorator ?
        thetax, thetay = self.get_normal_orientation_wrt_z()
        self._rot = _get_rotation_matrix(thetax, thetay)

        return

    @property
    def normal(self):
        return self._normal

    @normal.setter
    def normal(self, value):
        self._normal = np.asarray(value, dtype=np.float)

    def rotate_normal(self, thetax, thetay):
        """
        Rotates the current plane normal by fixed angles thetax and thetay.

        Parameters
        ----------
        thetax : float
            Angle of rotation around Ox (rad)
        thetay : float
            Angle of rotation around Oy (rad)

        """
        R = _get_rotation_matrix(thetax, thetay)
        self.normal = np.dot(R, self.normal)

        # updating self._rot
        self._rot = np.dot(R, self._rot)
        return

    def set_normal_from_angles(self, thetax, thetay):
        """
        Set the normal orientation given angles thetax and thetay.

        Parameters
        ----------
        thetax : float
            Angle around Ox (rad)
        thetay : float
            Anglearound Oy (rad)

        Remarks
        -------
        Equations are:
         $$\theta=\sqrt{\theta_x^2 + \theta_y^2}$$
         $$\sin{\theta} = \sqrt{n_x^2 + n_y^2}$$
         $$\theta_x = -\frac{\theta}{\sin{\theta}} n_y$$
         $$\theta_y =  \frac{\theta}{\sin{\theta}} n_x$$
         $$n_z = \cos{\theta}$$

        """

        theta = math.sqrt(thetax*thetax + thetay*thetay)

        if theta == 0.:
            self.normal[:] = [0., 0., 1.]

            # updating self._rot
            self._rot = np.eye(3)
        else:
            stheta_theta = math.sin(theta) / theta
            ctheta = math.cos(theta)
            self.normal[:] = np.array([stheta_theta*thetay,
                                       -stheta_theta*thetax,
                                       ctheta])
            # Updating self._rot
            self._rot = _get_rotation_matrix(thetax, thetay)
        return

    def get_normal_orientation_wrt_z(self):
        """
        Returns the angles thetax and thetay giving the orientation of the plane normal

        Returns
        -------
        thetax : float
            Angle around Ox (rad)
        thetay : float
            Angle around Oy (rad)
        """

        nx, ny, nz = self.normal
        stheta = math.sqrt(nx*nx + ny*ny)
        ctheta = nz
        if stheta == 0.:
            if nz == 1.:
                thetax = thetay = 0.
            elif nz == -1.:
                thetax = math.pi
                thetay = 0.
        else:
            theta = math.atan2(stheta, ctheta)
            theta_stheta = theta / stheta

            thetax = -theta_stheta * ny
            thetay = theta_stheta * nx

        return thetax, thetay

    def update_plane(self, c, thetax, thetay):
        """
        Updates the plane parameters (normal and scalar parameter) given scalar and angles.

        Parameters
        ----------
        c : float
            Plane scalar parameter (m)
        thetax : float
            Normal angle around Ox (rad)
        thetay : float
            Normal angle around Oy (rad)

        Remark
        ------
        Formula is (reformulate...):
        $$\theta = \sqrt{\theta_x^2 + \thetay^2}$$
        $$n_f = R'R z$$
        $$c_f = c*\cos{\theta} + c'$$

        """
        self.rotate_normal(thetax, thetay)
        ctheta = math.cos(math.sqrt(thetax*thetax + thetay*thetay))
        self.c = self.c * ctheta + c
        return

    def get_point_dist_wrt_plane(self, points):
        """
        Return the orthogonal distance of points with respect to the plane

        Parameters
        ----------
        points : ndarray
            Array of points coordinates

        Returns
        -------
        dist : ndarray
            Array of distances of points with respect to the plane
        """
        return np.dot(points, self.normal) - self.c

    def flip_normal(self):
        """
        Flips the Normal of the plane
        """
        self.normal *= -1
        thetax, thetay = self.get_normal_orientation_wrt_z()
        self._rot = _get_rotation_matrix(thetax, thetay)
        return

    def coord_in_plane(self, points):
        """
        Return the coordinates of points in the frame of the plane

        Parameters
        ----------
        points : ndarray
            Array of points coordinates

        Returns
        -------
        output : ndarray
            Array of points coordinates in the frame of the plane
        """
        # TODO: verifier effectivement que si on prend des points se trouvant dans le plan, leurs coordonnees dans le
        #  plan n'ont pas de composante z
        return -self.c*self.normal + np.transpose(np.dot(self._rot, points.T))

    def get_edge_intersection(self, P0, P1):
        """
        Returns the coordinates of the intersection point between the plane and the edge P0P1.

        Parameters
        ----------
        P0 : ndarray
            Coordinates of point P0
        P1 : ndarray
            Coordinates of point P1

        Returns
        -------
        I : ndarray
            Coordinates of intersection point
        """
        P0n = np.dot(P0, self.normal)
        P1n = np.dot(P1, self.normal)
        t = (P0n-self.c) / (P0n-P1n)
        if t<0. or t>1.:
            raise RuntimeError, 'Intersection is outside the edge'
        return (1-t)*P0 + t*P1

    def orthogonal_projection_on_plane(self, points):
        """
        Returns the coordinates of the orthogonal projection of points

        Parameters
        ----------
        points : ndarray
            Coordinates of the points to be projected

        Returns
        -------
        projected_points : ndarray
            Coordinates of the projection points
        """
        # TODO: passer en vectoriel
        projected_points = np.zeros_like(points)
        for point, projected_point in zip(points, projected_points):
            projected_point[:] = point - self.get_point_dist_wrt_plane(point) * self.normal

        return projected_points

class cached_property(object):
    def __init__(self, func):
        self.__doc__ = getattr(func, '__doc__')
        self.func = func

    def __get__(self, obj, cls):
        try:
            return obj._cached_properties[self.func.__name__]
        except:
            # print 'Computing %s and caching it' % self.func.__name__
            value = self.func(obj)
            try:
                obj._cached_properties[self.func.__name__] = value # FIXME: les valeurs sont enregistrees deux fois...
            except AttributeError:
                obj._cached_properties = {self.func.__name__: value}
            return value

class invalidate_cache(object):
    def __init__(self, func):
        self.__doc__ = getattr(func, '__doc__')
        self.func = func

    def __call__(self, cls, *args):
        self.func(cls, *args)
        # print 'Invalidation of the cache'
        try:
            cls._cached_properties.clear()
        except:
            cls._cached_properties = dict()

# TODO: Mettre la classe mesh dans

class Mesh(object):
    _ids = count(0)
    def __init__(self, vertices, faces, name=None):
        """

        Parameters
        ----------
        vertices : ndarray
            (nv x 3) Array of mesh vertices coordinates. Each line is a vertex.
        faces : ndarray
            Arrays of mesh connectivities for faces.

        Returns
        -------

        """

        self.V = vertices
        self.F = faces
        self._id = self._ids.next()
        if not name:
            self._name = 'mesh_%u' % self._id
        else:
            self._name = name

        self._verbose = False

    def __str__(self):
        str_repr = """
        --------------------------------------------
        \tMESH NAME : %s
        --------------------------------------------

        Number of vertices: %u
        Number of faces:    %u

        Number of triangles:   %u
        Number of quadrangles: %u

        xmin = %f\txmax = %f
        ymin = %f\tymax = %f
        zmin = %f\tzmax = %f


        """ % (self.name,
               self.nb_vertices,
               self.nb_faces,
               self.nb_triangles,
               self.nb_quadrangles,
               self.V[:, 0].min(),
               self.V[:, 0].max(),
               self.V[:, 1].min(),
               self.V[:, 1].max(),
               self.V[:, 2].min(),
               self.V[:, 2].max(),
               )
        return str_repr

    def print_quality(self):
        # This function is reproduced from
        # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/Verdict/Testing/Python/MeshQuality.py

        quality = vtk.vtkMeshQuality()
        quality.SetInput(self.vtk_polydata)

        def DumpQualityStats(iq, arrayname):
            an = iq.GetOutput().GetFieldData().GetArray(arrayname)
            cardinality = an.GetComponent(0, 4)
            range = list()
            range.append(an.GetComponent(0, 0))
            range.append(an.GetComponent(0, 2))
            average = an.GetComponent(0, 1)
            stdDev = math.sqrt(math.fabs(an.GetComponent(0, 3)))
            outStr = '%s%g%s%g\n%s%g%s%g' % (
                '    range: ', range[0], '  -  ', range[1],
                '    average: ', average, '  , standard deviation: ', stdDev)
            return outStr

        # Here we define the various mesh types and labels for output.
        meshTypes = [
            ['Triangle', 'Triangle',
             [['QualityMeasureToArea', ' Area Ratio:'],
              ['QualityMeasureToEdgeRatio', ' Edge Ratio:'],
              ['QualityMeasureToAspectRatio', ' Aspect Ratio:'],
              ['QualityMeasureToRadiusRatio', ' Radius Ratio:'],
              ['QualityMeasureToAspectFrobenius', ' Frobenius Norm:'],
              ['QualityMeasureToMinAngle', ' Minimal Angle:']
              ]
             ],

            ['Quad', 'Quadrilateral',
             [['QualityMeasureToArea', ' Area Ratio:'],
              ['QualityMeasureToEdgeRatio', ' Edge Ratio:'],
              ['QualityMeasureToAspectRatio', ' Aspect Ratio:'],
              ['QualityMeasureToRadiusRatio', ' Radius Ratio:'],
              ['QualityMeasureToMedAspectFrobenius',
               ' Average Frobenius Norm:'],
              ['QualityMeasureToMaxAspectFrobenius',
               ' Maximal Frobenius Norm:'],
              ['QualityMeasureToMinAngle', ' Minimal Angle:']
              ]
             ]
        ]

        if self.vtk_polydata.GetNumberOfCells() > 0:
            res = ''
            for meshType in meshTypes:
                res += '\n%s%s' % (meshType[1], ' quality of the mesh ')
                quality.Update()
                an = quality.GetOutput().GetFieldData().GetArray('Mesh ' + meshType[1] + ' Quality')
                cardinality = an.GetComponent(0, 4)

                res = ''.join((res, '(%u elements):\n' % (cardinality)))

                # res += '('+str(cardinality) +meshType[1]+'):\n'

                for measure in meshType[2]:
                    eval('quality.Set' + meshType[0] + measure[0] + '()')
                    quality.Update()
                    res += '\n%s\n%s' % (measure[1],
                                         DumpQualityStats(quality,
                                        'Mesh ' + meshType[1] + ' Quality'))
                res += '\n'

        info = """\n\nDefinition of the different quality measures is given
        in the verdict library manual :
        http://www.vtk.org/Wiki/images/6/6b/VerdictManual-revA.pdf\n"""

        res += info
        print res
        return

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        self._verbose = bool(value)

    def verbose_on(self):
        self._verbose = True

    def verbose_off(self):
        self._verbose = False

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = str(value)
        return

    @property
    def nb_vertices(self):
        return self._V.shape[0]

    @property
    def nb_faces(self):
        return self._F.shape[0]

    @property
    def V(self):
        # print 'getting V'
        return self._V.copy()

    @property
    def F(self):
        # print 'getting F'
        return self._F.copy()

    @V.setter
    @invalidate_cache
    def V(self, value):
        # print 'setting V'
        self._V = np.asarray(value, dtype=np.float).copy()
        self._V.setflags(write=False)
        return

    @F.setter
    @invalidate_cache
    def F(self, value):
        # print 'setting F'
        self._F = np.asarray(value, dtype=np.int).copy()
        self._F.setflags(write=False)
        return

    @cached_property
    def faces_properties(self):
        # areas, normals, centers = mm.get_all_faces_properties(self._V, self._F)
        nf = self.nb_faces

        # triangle_mask = F[:, 0] == F[:, -1]
        # nb_triangles = np.sum(triangle_mask)
        # quads_mask = np.invert(triangle_mask)
        # nb_quads = nf - nb_triangles

        areas = np.zeros(nf, dtype=np.float)
        normals = np.zeros((nf, 3), dtype=np.float)
        centers = np.zeros((nf, 3), dtype=np.float)

        # Collectively dealing with triangles
        # triangles = F[triangle_mask]
        triangles_id = self.triangles_ids
        triangles = self._F[triangles_id]

        triangles_normals = np.cross(self._V[triangles[:, 1]] - self._V[triangles[:, 0]],
                                     self._V[triangles[:, 2]] - self._V[triangles[:, 0]])
        triangles_areas = np.linalg.norm(triangles_normals, axis=1)
        normals[triangles_id] = triangles_normals / np.array(([triangles_areas, ] * 3)).T
        areas[triangles_id] = triangles_areas / 2.
        centers[triangles_id] = np.sum(self._V[triangles[:, :3]], axis=1) / 3.

        # Collectively dealing with quads
        quads_id = self.quadrangles_ids
        quads = self._F[quads_id]
        # quads = F[quads_mask]

        quads_normals = np.cross(self._V[quads[:, 2]] - self._V[quads[:, 0]],
                                 self._V[quads[:, 3]] - self._V[quads[:, 1]])
        normals[quads_id] = quads_normals / np.array(([np.linalg.norm(quads_normals, axis=1), ] * 3)).T

        a1 = np.linalg.norm(np.cross(self._V[quads[:, 1]] - self._V[quads[:, 0]],
                                     self._V[quads[:, 2]] - self._V[quads[:, 0]]), axis=1) * 0.5
        a2 = np.linalg.norm(np.cross(self._V[quads[:, 3]] - self._V[quads[:, 0]],
                                     self._V[quads[:, 2]] - self._V[quads[:, 0]]), axis=1) * 0.5
        areas[quads_id] = a1 + a2

        C1 = np.sum(self._V[quads[:, :3]], axis=1) / 3.
        C2 = (np.sum(self._V[quads[:, 2:4]], axis=1) + self._V[quads[:, 0]]) / 3.

        centers[quads_id] = (np.array(([a1, ] * 3)).T * C1 + np.array(([a2, ] * 3)).T * C2)
        centers[quads_id] /= np.array(([areas[quads_id], ] * 3)).T

        faces_properties = {'areas': areas,
                            'normals': normals,
                            'centers': centers}

        return faces_properties

    @property
    def faces_areas(self):
        return self.faces_properties['areas']

    @property
    def faces_centers(self):
        return self.faces_properties['centers']

    @property
    def faces_normals(self):
        return self.faces_properties['normals']

    @cached_property
    def triangles_quadrangles_ids(self):
        triangle_mask = (self._F[:, 0] == self._F[:, -1])
        quadrangles_mask = np.invert(triangle_mask)
        triangles_quadrangles = {'triangles_ids': np.where(triangle_mask)[0],
                                 'quadrangles_ids': np.where(quadrangles_mask)[0]}

        return triangles_quadrangles

    @property
    def triangles_ids(self):
        return self.triangles_quadrangles_ids['triangles_ids']

    @property
    def nb_triangles(self):
        return len(self.triangles_ids)

    @property
    def quadrangles_ids(self):
        return self.triangles_quadrangles_ids['quadrangles_ids']

    @property
    def nb_quadrangles(self):
        return len(self.quadrangles_ids)

    def is_triangle(self, face_id):
        if self._F[face_id, 0] == self._F[face_id, -1]:
            return True
        else:
            return False

    def get_face(self, face_id):
        if self.is_triangle(face_id):
            return self._F[face_id, :3]
        else:
            return self._F[face_id]

    def extract_faces(self, id_faces_to_extract, return_id_extracted_vertices=False):
        nv = self.nb_vertices

        # Determination of the vertices to keep
        Vmask = np.zeros(nv, dtype=bool)
        Vmask[self._F[id_faces_to_extract].flatten()] = True
        idV = np.arange(nv)[Vmask]

        # Building up the vertex array
        V_extracted = self._V[idV]
        newID_V = np.arange(nv)
        newID_V[idV] = np.arange(len(idV))

        F_extracted = self._F[id_faces_to_extract]
        F_extracted = newID_V[F_extracted.flatten()].reshape((len(id_faces_to_extract), 4))

        extracted_mesh = Mesh(V_extracted, F_extracted)
        extracted_mesh._verbose = self._verbose

        if return_id_extracted_vertices:
            return extracted_mesh, idV
        else:
            return extracted_mesh

    @cached_property
    def vtk_polydata(self):

        # Create a vtkPoints object and store the points in it
        points = vtk.vtkPoints()
        for point in self._V:
            points.InsertNextPoint(point)

        # Create a vtkCellArray to store faces
        faces = vtk.vtkCellArray()
        for face_ids in self._F:
            if face_ids[0] == face_ids[-1]:
                # Triangle
                curface = face_ids[:3]
                vtk_face = vtk.vtkTriangle()
            else:
                # Quadrangle
                curface = face_ids[:4]
                vtk_face = vtk.vtkQuad()

            for idx, id in enumerate(curface):
                vtk_face.GetPointIds().SetId(idx, id)

            faces.InsertNextCell(vtk_face)

        vtk_polydata = vtk.vtkPolyData()
        vtk_polydata.SetPoints(points)
        vtk_polydata.SetPolys(faces)

        return vtk_polydata

    def show(self):
        self.viewer = MMviewer.MMViewer()
        self.viewer.add_polydata(self.vtk_polydata)
        self.viewer.show()
        self.viewer.finalize()

    @cached_property
    def _connectivity(self):

        nv = self.nb_vertices
        nf = self.nb_faces

        mesh_closed = True

        # Building connectivities

        # Establishing VV and VF connectivities
        VV = dict([(i, set()) for i in xrange(nv)])
        VF = dict([(i, set()) for i in xrange(nv)])
        for (iface, face) in enumerate(self._F):
            if face[0] == face[-1]:
                face_w = face[:3]
            else:
                face_w = face
            for (index, iV) in enumerate(face_w):
                VF[iV].add(iface)
                VV[face_w[index - 1]].add(iV)
                VV[iV].add(face_w[index - 1])

        # Connectivity FF
        boundary_edges = dict()

        FF = dict([(i, set()) for i in xrange(nf)])
        for ivertex in xrange(nv):
            S1 = VF[ivertex]
            for iadjV in VV[ivertex]:
                S2 = VF[iadjV]
                I = list(S1 & S2)
                if len(I) == 2:
                    FF[I[0]].add(I[1])
                    FF[I[1]].add(I[0])

                elif len(I) == 1:
                    boundary_face = self._F[I[0]]

                    if boundary_face[0] == boundary_face[-1]:
                        boundary_face = boundary_face[:3]
                    ids = np.where((boundary_face == ivertex) + (boundary_face == iadjV))[0]

                    if ids[1] != ids[0]+1:
                        iV_orig, iV_target = boundary_face[ids]
                    else:
                        iV_target, iV_orig = boundary_face[ids]

                    boundary_edges[iV_orig] = iV_target
                else:
                    raise RuntimeError, 'Unexpected error while computing mesh connectivities'

        # Computing boundaries
        boundaries = list()
        while True:
            try:
                boundary = list()
                iV0_init, iV1 = boundary_edges.popitem()
                boundary.append(iV0_init)
                boundary.append(iV1)
                iV0 = iV1

                while True:
                    try:
                        iV1 = boundary_edges.pop(iV0)
                        boundary.append(iV1)
                        iV0 = iV1
                    except KeyError:
                        if boundary[0] != boundary[-1]:
                            print 'Boundary is not closed !!!'
                        else:
                            boundaries.append(boundary)
                        break
            except KeyError:
                break

        connectivity = {'VV': VV,
                        'VF': VF,
                        'FF': FF,
                        'boundaries': boundaries}
        return connectivity

    @property
    def VV(self):
        return self._connectivity['VV']

    @property
    def VF(self):
        return self._connectivity['VF']

    @property
    def FF(self):
        return self._connectivity['FF']

    @property
    def boundaries(self):
        return self._connectivity['boundaries']

    @property
    def nb_boundaries(self):
        return len(self._connectivity['boundaries'])

    def is_mesh_closed(self):
        return len(self._connectivity['boundaries']) == 0

    def is_mesh_conformal(self):
        # FIXME: bugge
        tol = 1e-7

        boundaries = self._connectivity['boundaries']
        polygons_areas = np.zeros(len(boundaries), dtype=np.float)

        conformal = True

        for boundary in boundaries:
            boundary_vertices = self._V[boundary]
            # Si les trois projections (Oxy, Oxz, Oyz) de boundary sont des courbes dans ces plans, alors la courbe
            # est colapsee
            # Projecting on Oxy
            plane = Plane(normal=[0., 0., 1.])
            proj0 = plane.orthogonal_projection_on_plane(boundary_vertices)
            # Projecting on Oyz
            plane.normal = [1., 0., 0.]
            proj1 = plane.orthogonal_projection_on_plane(boundary_vertices)
            # Projecting on Oxz
            plane.normal = [0., 1., 0.]
            proj2 = plane.orthogonal_projection_on_plane(boundary_vertices)

            # Compputing areas of curves
            x = proj0[:, 0]
            y = proj0[:, 1]
            a0 = ((np.roll(y, -1)-y) * (np.roll(x, -1)+x)).sum()

            y = proj1[:, 0]
            z = proj1[:, 1]
            a1 = ((np.roll(y, -1)-y) * (np.roll(z, -1)+z)).sum()

            x = proj2[:, 0]
            z = proj2[:, 1]
            a2 = ((np.roll(x, -1)-x) * (np.roll(z, -1)+z)).sum()

            if math.fabs(a0) < tol and math.fabs(a1) < tol and math.fabs(a2) < tol:
                conformal = False

            return conformal

    def rotate_x(self, thetax):
        self.rotate([thetax, 0., 0.])
        return

    def rotate_y(self, thetay):
        self.rotate([0., thetay, 0.])
        return

    def rotate_z(self, thetaz):
        self.rotate([0., 0., thetaz])

    # @invalidate_cache
    def rotate(self, angles):
        angles = np.asarray(angles, dtype=np.float)
        theta = np.linalg.norm(angles)
        if theta == 0.:
            return

        ctheta = math.cos(theta)
        stheta = math.sin(theta)

        nx, ny, nz = angles/theta
        nxny = nx*ny
        nxnz = nx*nz
        nynz = ny*nz
        nx2 = nx*nx
        ny2 = ny*ny
        nz2 = nz*nz

        R = ctheta*np.eye(3) \
          + (1-ctheta) * np.array([[nx2, nxny, nxnz],
                                   [nxny, ny2, nynz],
                                   [nxnz, nynz, nz2]]) \
          + stheta * np.array([[0., -nz, ny],
                               [nz, 0., -nx],
                               [-ny, nx, 0.]])

        self.V = np.transpose(np.dot(R, self._V.copy().T))

        return

    # @invalidate_cache
    def translate_x(self, tx):
        V = self._V.copy()
        V[:, 0] += tx
        self.V = V
        return

    # @invalidate_cache
    def translate_y(self, ty):
        V = self._V.copy()
        V[:, 1] += ty
        self.V = V
        return

    # @invalidate_cache
    def translate_z(self, tz):
        V = self._V.copy()
        V[:, 2] += tz
        self.V = V
        return

    # @invalidate_cache
    def translate(self, t):
        # t = np.asarray(t, dtype=np.float)
        tx, ty, tz = t
        V = self._V.copy()
        V[:, 0] += tx
        V[:, 1] += ty
        V[:, 2] += tz
        self.V = V
        return

    def scale(self, alpha):
        V = self._V.copy()
        V *= alpha
        self.V = V
        return

    def flip_normals(self):
        F = self._F.copy()
        self.F = np.fliplr(F)
        return

    def __add__(self, mesh):
        V = np.concatenate((self.V, mesh.V), axis=0)
        F = np.concatenate((self.F, mesh.F+self.nb_vertices), axis=0)
        new_mesh = Mesh(V, F, name='_'.join([self.name, mesh.name]))
        # new_mesh.merge_duplicates()
        new_mesh._verbose = self._verbose or mesh._verbose
        return new_mesh

    def copy(self):
        return copy.deepcopy(self)

    def merge_duplicates(self, tol=1e-8, return_index=False):
        # TODO: voir ou mettre l'implementation de la fonction merge_duplicates
        output = mm.merge_duplicates(self.V, self.F, verbose=False, tol=tol, return_index=return_index)
        if return_index:
            V, F, newID = output
        else:
            V, F = output
        if self._verbose:
            print "* Merging duplicate vertices"
            delta_n = self.nb_vertices - V.shape[0]
            if delta_n > 0:
                print "\t--> %u vertices have been merged" % delta_n
            else:
                print "\t--> No duplicate vertices have been found"
        self.V, self.F = V, F
        if return_index:
            return newID
        else:
            return

    def heal_normals(self):

        # TODO: return the different groups of a mesh in case it is made of several unrelated groups

        nv = self.nb_vertices
        nf = self.nb_faces
        F = self.F

        # Building connectivities
        VV = self.VV
        VF = self.VF
        FF = self.FF
        boundaries = self.boundaries

        if len(boundaries) > 0:
            mesh_closed = False
        else:
            mesh_closed = True

        # Flooding the mesh to find inconsistent normals
        type_cell = np.zeros(nf, dtype=np.int32)
        type_cell[:] = 4
        # triangles_mask = F[:, 0] == F[:, -1]
        type_cell[self.triangles_ids] = 3

        FVis = np.zeros(nf, dtype=bool)
        FVis[0] = True
        stack = [0]
        nb_reversed = 0
        while 1:
            if len(stack) == 0:
                not_FVis = np.logical_not(FVis)
                if np.any(not_FVis):
                    iface = np.where(np.logical_not(FVis))[0][0]
                    stack.append(iface)
                    FVis[iface] = True
                else:
                    break

            iface = stack.pop()
            face = F[iface]
            S1 = set(face)

            for iadjF in FF[iface]:
                if FVis[iadjF]:
                    continue
                FVis[iadjF] = True
                # Removing the other pointer
                FF[iadjF].remove(iface)  # So as it won't go from iadjF to iface in the future

                # Shared vertices
                adjface = F[iadjF]
                S2 = set(adjface)
                # try:
                common_vertices = list(S1 & S2)
                if len(common_vertices) == 2:
                    iV1, iV2 = common_vertices
                else:
                    print 'WARNING: faces %u and %u have more than 2 vertices in common !' % (iface, iadjF)
                    continue

                # Checking normal consistency
                face_ref = np.roll(face[:type_cell[iface]], -np.where(face == iV1)[0][0])
                adj_face_ref = np.roll(adjface[:type_cell[iadjF]], -np.where(adjface == iV1)[0][0])

                if face_ref[1] == iV2:
                    i = 1
                else:
                    i = -1

                if adj_face_ref[i] == iV2:
                    # Reversing normal
                    nb_reversed += 1
                    F[iadjF] = np.flipud(F[iadjF])

                # Appending to the stack
                stack.append(iadjF)

        if self._verbose:
            print "* Healing normals to make them consistent and if possible outward"
            if nb_reversed > 0:
                print '\t--> %u faces have been reversed to make normals consistent across the mesh' % (nb_reversed)
            else:
                print "\t--> Normals orientations are consistent"

        self.F = F

        # Checking if the normals are outward
        if mesh_closed:
            zmax = np.max(self._V[:, 2])

            areas = self.faces_areas
            normals = self.faces_normals
            centers = self.faces_centers
            # areas, normals, centers = get_all_faces_properties(V, F)

            hs = (np.array([(centers[:, 2] - zmax) * areas, ] * 3).T * normals).sum(axis=0)

            tol = 1e-9
            if math.fabs(hs[0]) > tol or math.fabs(hs[1]) > tol:
                if self._verbose:
                    print "\t--> WARNING: the mesh does not seem watertight althought marked as closed..."

            if hs[2] < 0:
                flipped = True
                self.flip_normals()
            else:
                flipped = False

            if self._verbose and flipped:
                print '\t--> Every normals have been reversed to be outward'


        else:
            if self._verbose:
                print "\t--> Mesh is not closed, meshmagick cannot test if the normals are outward"

        return

    def remove_unused_vertices(self):

        nv = self.nb_vertices
        V, F = self.V, self.F

        usedV = np.zeros(nv, dtype=np.bool)
        usedV[sum(map(list, F), [])] = True
        nb_usedV = sum(usedV)

        if nb_usedV < nv:
            newID_V = np.arange(nv)
            newID_V[usedV] = np.arange(nb_usedV)
            F = newID_V[F]
            V = V[usedV]

        self.V, self.F = V, F

        if self._verbose:
            print "* Removing unused vertices in the mesh:"
            if nb_usedV < nv:
                unusedV = np.where(np.logical_not(usedV))[0]
                vlist_str = '[' + ', '.join(str(iV) for iV in unusedV) + ']'
                print "\t--> %u unused vertices have been removed" % (nv - nb_usedV)
            else:
                print "\t--> No unused vertices"

        return

    def heal_triangles(self):

        F = self.F

        quads = F[:, 0] != F[:, -1]
        nquads_init = sum(quads)

        F[quads] = np.roll(F[quads], 1, axis=1)
        quads = F[:, 0] != F[:, -1]

        F[quads] = np.roll(F[quads], 1, axis=1)
        quads = F[:, 0] != F[:, -1]

        F[quads] = np.roll(F[quads], 1, axis=1)
        quads = F[:, 0] != F[:, -1]
        nquads_final = sum(quads)

        self.F = F

        if self._verbose:
            print "* Ensuring consistent definition of triangles:"
            if nquads_final < nquads_init:
                print "\t--> %u triangles were described the wrong way and have been corrected" % (
                nquads_init - nquads_final)
            else:
                print "\t--> Triangle description is consistent"

        return

    def remove_degenerated_faces(self, rtol=1e-5):

        areas = self.faces_areas
        area_threshold = areas.mean() * rtol

        # Detecting faces that have null area
        F = self.F[np.logical_not(areas < area_threshold)]
        if self._verbose:
            nb_removed = self.nb_faces - F.shape[0]
            print '* Removing degenerated faces'
            if nb_removed > 0:
                print '\t-->%u degenerated faces have been removed' % nb_removed
            else:
                print '\t--> No degenerated faces'

        self.F = F
        return

    def heal_mesh(self):
        self.remove_unused_vertices()
        self.remove_degenerated_faces()
        self.merge_duplicates()
        self.heal_triangles()
        self.heal_normals()
        return

    def triangulate_quadrangles(self):
        # TODO: Ensure the best quality aspect ratio of generated triangles

        # Defining both triangles id lists to be generated from quadrangles
        T1 = (0, 1, 2)
        T2 = (0, 2, 3)

        F = self.F

        # Triangulation
        new_faces = F[self.quadrangles_ids].copy()
        new_faces[:, :3] = new_faces[:, T1]
        new_faces[:, -1] = new_faces[:, 0]

        F[self.quadrangles_ids, :3] = F[:, T2][self.quadrangles_ids]
        F[self.quadrangles_ids, -1] = F[self.quadrangles_ids, 0]

        F = np.concatenate((F, new_faces))

        if self._verbose:
            print '\nTriangulating quadrangles'
            if self.nb_quadrangles != 0:
                print '\t-->{:d} quadrangles have been split in triangles'.format(self.nb_quadrangles)

        self.F = F

        return F

    def symmetrize(self, plane):

        # Symmetrizing the nodes
        V, F = self.V, self.F

        # normal = plane.normal / np.dot(plane.normal, plane.normal)
        V = np.concatenate((V, V - 2 * np.outer(np.dot(V, plane.normal) - plane.c, plane.normal)))
        F = np.concatenate((F, np.fliplr(F.copy() + self.nb_vertices)))

        self.V, self.F = V, F
        verbose = self.verbose
        self.verbose_off()
        self.merge_duplicates()
        self.verbose = verbose
        return

    def _partition_wrt_plane(self, plane, tol=1e-4):

        distances = plane.get_point_dist_wrt_plane(self._V)

        vertices_above_mask = distances > tol
        vertices_on_mask = np.fabs(distances) < tol
        vertices_below_mask = distances < -tol

        nb_vertices_above = vertices_above_mask[self._F].sum(axis=1)
        nb_vertices_on = vertices_on_mask[self._F].sum(axis=1)
        nb_vertices_below = vertices_below_mask[self._F].sum(axis=1)

        above_faces_mask = nb_vertices_above == 4
        below_faces_mask = nb_vertices_below == 4

        crown_faces_mask = np.logical_and(np.logical_not(above_faces_mask), np.logical_not(below_faces_mask))

        partition = {
            'distances': distances,
            'vertices_above_mask': vertices_above_mask,
            'vertices_on_mask': vertices_on_mask,
            'vertices_below_mask': vertices_below_mask,
            'above_faces_ids': np.where(above_faces_mask)[0],
            'crown_faces_ids': np.where(crown_faces_mask)[0],
            'below_faces_ids': np.where(below_faces_mask)[0]
        }

        return partition

    # def get_crown_faces_ids(self, plane):
    #     return self._partition_wrt_plane(plane)['crown_faces_ids']


    def clip(self, plane, return_boundaries=False, assert_closed_boundaries=False):

        nv = self.nb_vertices

        partition = self._partition_wrt_plane(plane)

        below_mesh = self.extract_faces(partition['below_faces_ids'])
        clipped_crown_mesh, idV = self.extract_faces(partition['crown_faces_ids'], return_id_extracted_vertices=True)
        crown_vertices = clipped_crown_mesh.V

        vertices_above_mask = partition['vertices_above_mask'][idV]
        vertices_on_mask = partition['vertices_on_mask'][idV]
        vertices_below_mask = partition['vertices_below_mask'][idV]

        # TODO: Vertices pre-projection
        # vertices_on = partition['vertices_on']
        # V[vertices_on] = plane.orthogonal_projection_on_plane(V[vertices_on])
        # pos[vertices_on] = 0.

        # Init
        crown_faces = list()
        direct_boundary_edges = dict()
        inv_boundary_edges = dict()
        intersections = list()

        nI = clipped_crown_mesh.nb_vertices

        for face_id in xrange(clipped_crown_mesh.nb_faces):

            face = clipped_crown_mesh.get_face(face_id)

            # # Determining the type of face clipping
            v_above_face = np.where(vertices_above_mask[face])[0]
            v_on_face = np.where(vertices_on_mask[face])[0]
            v_below_face = np.where(vertices_below_mask[face])[0]

            nb_above = len(v_above_face)
            nb_on = len(v_on_face)
            nb_below = len(v_below_face)

            face_type = str(nb_above) + str(nb_on) + str(nb_below)

            if face_type == '202':  # Done
                #    0*-----*3
                #     |     |
                # ----o-----o----
                #     |     |
                #    1*-----*2
                if v_above_face[1] == v_above_face[0] + 1:
                    face = np.roll(face, -v_above_face[1])
                P0, P1, P2, P3 = crown_vertices[face]
                Ileft = plane.get_edge_intersection(P0, P1)
                Iright = plane.get_edge_intersection(P2, P3)
                intersections += [Ileft, Iright]
                boundary_edge = [nI, nI + 1]
                crown_faces.append([nI, face[1], face[2], nI + 1])
                nI += 2

            elif face_type == '301':  # Done
                #      *2
                #     / \
                #    /   \
                #   /     \
                # 3*       *1
                #   \     /
                # ---o---o---
                #     \ /
                #      *0
                face = np.roll(face, -v_below_face[0])
                P0, P1, P3 = crown_vertices[face[[0, 1, 3]]]
                Ileft = plane.get_edge_intersection(P0, P3)
                Iright = plane.get_edge_intersection(P0, P1)
                intersections += [Ileft, Iright]
                boundary_edge = [nI, nI + 1]
                crown_faces.append([nI, face[0], nI + 1, nI])
                nI += 2

            elif face_type == '103':  # Done
                #      *0
                #     / \
                # ---o---o---
                #   /     \
                # 1* - - - *3
                #   \     /
                #    \   /
                #     \ /
                #      *2
                face = np.roll(face, -v_above_face[0])
                P0, P1, P3 = crown_vertices[face[[0, 1, 3]]]
                Ileft = plane.get_edge_intersection(P0, P1)
                Iright = plane.get_edge_intersection(P0, P3)
                intersections += [Ileft, Iright]
                boundary_edge = [nI, nI + 1]
                crown_faces.append([nI, face[1], face[3], nI + 1])
                crown_faces.append([face[1], face[2], face[3], face[1]])
                nI += 2

            elif face_type == '102':  # Done
                #      *O
                #     / \
                # ---o---o---
                #   /     \
                # 1*-------*2
                face = np.roll(face, -v_above_face[0])
                P0, P1, P2 = crown_vertices[face]
                Ileft = plane.get_edge_intersection(P0, P1)
                Iright = plane.get_edge_intersection(P0, P2)
                intersections += [Ileft, Iright]
                boundary_edge = [nI, nI + 1]
                crown_faces.append([nI, face[1], face[2], nI + 1])
                nI += 2

            elif face_type == '201':  # done
                #  2*-------*1
                #    \     /
                #  ---o---o---
                #      \ /
                #       *0
                face = np.roll(face, -v_below_face[0])
                P0, P1, P2 = crown_vertices[face]
                Ileft = plane.get_edge_intersection(P0, P2)
                Iright = plane.get_edge_intersection(P0, P1)
                intersections += [Ileft, Iright]
                boundary_edge = [nI, nI + 1]
                crown_faces.append([nI, face[0], nI + 1, nI])
                nI += 2

            elif face_type == '211':  # Done
                #        *3                   *1
                #       / \                  / \
                #      /   *2       or     2*   \
                #    0/   /                  \   \0
                # ---*---o---              ---o---*---
                #     \ /                      \ /
                #      *1                       *3
                #

                face = np.roll(face, -v_on_face[0])
                if pos[face[1]] < 0.:
                    P1, P2 = crown_vertices[face[[1, 2]]]
                    Iright = plane.get_edge_intersection(P1, P2)
                    intersections.append(Iright)
                    boundary_edge = [face[0], nI]
                    crown_faces.append([face[0], face[1], nI, face[0]])
                else:
                    P2, P3 = crown_vertices[face[[2, 3]]]
                    Ileft = plane.get_edge_intersection(P2, P3)
                    intersections.append(Ileft)
                    boundary_edge = [nI, face[0]]
                    crown_faces.append([nI, face[3], face[0], nI])
                nI += 1

            elif face_type == '112':  # Done
                #       *3                     *1
                #      / \                    / \
                #  ---*---o---      or    ---o---*---
                #     0\   \                /   /0
                #       \   *2            2*   /
                #        \ /                \ /
                #         *1                 *3
                face = np.roll(face, -v_on_face[0])
                if pos[face[1]] < 0.:
                    P2, P3 = crown_vertices[face[[2, 3]]]
                    Iright = plane.get_edge_intersection(P2, P3)
                    intersections.append(Iright)
                    boundary_edge = [face[0], nI]
                    crown_faces.append([face[0], face[1], face[2], nI])
                else:
                    P1, P2 = crown_vertices[face[[1, 2]]]
                    Ileft = plane.get_edge_intersection(P1, P2)
                    intersections.append(Ileft)
                    boundary_edge = [nI, face[0]]
                    crown_faces.append([nI, face[2], face[3], face[0]])
                nI += 1

            elif face_type == '013':  # Done
                # -----*-----
                #     / \
                #    /   \
                #   *     *
                #    \   /
                #     \ /
                #      *
                boundary_edge = None
                crown_faces.append(list(face))

            elif face_type == '210' or face_type == '310':  # Done
                #   *-------*               *
                #    \ 210 /               / \ 310
                #     \   /               *   *
                #      \ /                 \ /
                #   ----*----           ----*----
                boundary_edge = None

            elif face_type == '111':  # Done
                #        *2              *1
                #       /|               |\
                #      / |               | \
                #  ---*--o---    or   ---o--*---
                #     0\ |               | /0
                #       \|               |/
                #        *1              *2
                face = np.roll(face, -v_on_face[0])
                P1, P2 = crown_vertices[face[[1, 2]]]
                if pos[face[1]] < 0.:
                    Iright = plane.get_edge_intersection(P1, P2)
                    intersections.append(Iright)
                    boundary_edge = [face[0], nI]
                    crown_faces.append([face[0], face[1], nI, face[0]])
                else:
                    Ileft = plane.get_edge_intersection(P1, P2)
                    intersections.append(Ileft)
                    boundary_edge = [nI, face[0]]
                    crown_faces.append([nI, face[2], face[0], nI])
                nI += 1

            elif face_type == '120':  # Done
                #         *O
                #        / \
                #       /   \
                #     1/     \2
                # ----*-------*----
                face = np.roll(face, -v_above_face[0])
                boundary_edge = [face[1], face[2]]

            elif face_type == '021':  # Done
                #  ----*-------*----
                #      2\     /1
                #        \   /
                #         \ /
                #          *0
                face = np.roll(face, -v_below_face[0])
                boundary_edge = [face[2], face[1]]
                face = list(face)
                face.append(face[0])
                crown_faces.append(face)

            elif face_type == '022':
                # ----*-----*----
                #    0|     |3
                #     |     |
                #    1*-----*2
                if v_on_face[1] == v_on_face[0] + 1:
                    face = np.roll(face, -v_on_face[1])
                boundary_edge = [face[0], face[3]]
                crown_faces.append(list(face))

            elif face_type == '012':  # Done
                #   ------*------
                #        / \
                #       /   \
                #      /     \
                #     *-------*
                boundary_edge = None
                face = list(face)
                face.append(face[0])
                crown_faces.append(face)

            elif face_type == '220':  # Done
                #    0*-----*3
                #     |     |
                #    1|     |2
                # ----*-----*----
                if v_above_face[1] == v_above_face[0] + 1:
                    face = np.roll(face, -v_above_face[1])
                boundary_edge = [face[1], face[2]]

            elif face_type == '121':  # Done
                #       *0
                #      / \
                #     /   \
                # ---*-----*---
                #    1\   /3
                #      \ /
                #       *2
                face = np.roll(face, -v_above_face[0])
                boundary_edge = [face[1], face[3]]
                crown_faces.append([face[1], face[2], face[3], face[1]])

            elif face_type == '300' or face_type == '400':
                #       *               *-----*
                #      / \              |     |
                #     /300\       or    | 400 |
                #    *-----*            *-----*
                # ____________       ______________
                boundary_edge = None

            elif face_type == '003':
                #  -----------
                #       *
                #      / \
                #     /   \
                #    *-----*
                boundary_edge = None
                face = list(face)
                face.append(face[0])
                crown_faces.append(face)

            elif face_type == '004':
                #  ---------------
                #      *-----*
                #      |     |
                #      |     |
                #      *-----*
                boundary_edge = None
                crown_faces.append(list(face))

            # Building boundary connectivity
            if boundary_edge is not None:
                direct_boundary_edges[boundary_edge[0]] = boundary_edge[1]
                inv_boundary_edges[boundary_edge[1]] = boundary_edge[0]

        crown_vertices = np.concatenate((crown_vertices, intersections))

        clipped_crown_mesh = Mesh(crown_vertices, crown_faces)

        newID = clipped_crown_mesh.merge_duplicates(return_index=True)

        immersed_mesh = clipped_crown_mesh + below_mesh

        if return_boundaries:
            # Updating dictionaries
            direct_boundary_edges = dict(
                zip(newID[direct_boundary_edges.keys()], newID[direct_boundary_edges.values()]))
            inv_boundary_edges = dict(zip(newID[inv_boundary_edges.keys()], newID[inv_boundary_edges.values()]))

            # Ordering boundary edges in continuous lines
            closed_polygons = list()
            open_lines = list()
            while True:
                try:
                    line = list()
                    V0_init, V1 = direct_boundary_edges.popitem()
                    line.append(V0_init)
                    line.append(V1)
                    V0 = V1

                    while True:
                        try:
                            V1 = direct_boundary_edges.pop(V0)
                            line.append(V1)
                            V0 = V1
                        except KeyError:
                            if line[0] != line[-1]:
                                # Trying to find an other queue
                                queue = list()
                                V0 = V0_init
                                while True:
                                    try:
                                        V1 = inv_boundary_edges[V0]
                                        direct_boundary_edges.pop(V1)
                                        queue.append(V1)
                                        V0 = V1
                                    except:
                                        queue.reverse()
                                        line = queue + line
                                        open_lines.append(line)
                                        break
                            else:
                                closed_polygons.append(line)
                            break

                except:
                    # TODO: retirer les deux lines suivantes
                    if self._verbose:
                        print "%u closed polygon\n%u open curve" % (len(closed_polygons), len(open_lines))
                    # print open_lines[0]

                    if assert_closed_boundaries:
                        if len(open_lines) > 0:
                            # mm.write_VTP('mesh_clip.vtp', crown_vertices, crown_faces)
                            raise RuntimeError, 'Open intersection curve found'

                    break
            boundaries = {
                'closed_polygons': closed_polygons,
                'open_lines': open_lines
            }
            return immersed_mesh, boundaries

        else:  # return_boundaries = False
            # mm.show(crown_vertices, crown_faces)
            return immersed_mesh


# def clip_mesh_against_plane(V, F):
#     pass

if __name__ == '__main__':

    V, F = mm.load_VTP('SEAREV.vtp')
    plane = Plane()
    # V, F = mm.clip_by_plane(V, F, plane)
    mesh = Mesh(V, F)
    print mesh

    # mesh.verbose_on()
    #
    # clipped_mesh, boundaries = mesh.clip(plane, return_boundaries=True)



    sys.exit(0)

#####################################

    plane = Plane()

    # clip(V, F, plane)
    # sys.exit(0)

    # V, F = mm.symmetrize(V, F, plane)
    # plane.normal = np.array([0, 1, 0])
    # V, F = mm.symmetrize(V, F, plane)

    # Rotation aleatoire
    import hydrostatics as hs

    # R = np.load('buggy_rotation_meshmagick.npy')
    # V = np.transpose(np.dot(R, V.T))
    # normal = [0, 0, 1]
    # c = 0
    # plane = Plane(normal, c)
    # crown_face_ids, below_face_ids = partition_mesh(V, F, plane)
    # Vclip, F_crown, boundary_polygons, boundary_lines = _clip(V, F[crown_face_ids], plane, return_boundaries=True)
    # sys.exit(0)

    iter = 0
    while True:
        iter += 1
        print '\n', iter

        Vc = V.copy()
        thetax, thetay = np.random.rand(2) * 2*math.pi

        R = hs._get_rotation_matrix(thetax, thetay)
        Vc = np.transpose(np.dot(R, Vc.T))

        # mm.show(Vc, F)
        normal = [0, 0, 1]
        c = 0
        plane = Plane(normal, c)

        crown_face_ids, below_face_ids = split_mesh(Vc, F, plane)

        try:
            Vclip, F_crown, boundary_polygons, boundary_lines = _clip(Vc, F[crown_face_ids], plane, return_boundaries=True)
            # mm.show(Vclip, F_crown)
        except:
            np.save('buggy_rotation_meshmagick', R)
            mm.write_VTP('all.vtp', Vc, F)
            mm.write_VTP('ring.vtp', Vc, F[crown_face_ids])
            print R
            print 'state saved'
            break
