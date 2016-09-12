#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
This module is part of meshmagick software. It provides functions for mesh clipping purpose.
"""

from tools import merge_duplicate_rows
import MMviewer
import numpy as np
import math
import copy
import vtk
from itertools import count
from warnings import warn
import sys # TODO: Retirer


# TODO: les points doivent etre des objects nodes...
# TODO: On doit pouvoir specifier des objets frame

# TODO: voir si on ne peut pas mettre ces fonctions dans un module dedie --> module rotation !!!
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

def _get_axis_angle_from_rotation_matrix(rot):
    warn('Fonction _get_axis_angle_from_rotation_matrix a verifier !!!')
    theta = math.acos((np.trace(rot)-1.)*0.5)
    direction = (1./(2.*math.sin(theta))) * np.array([rot[2, 1]-rot[1, 2],
                                                      rot[0, 2]-rot[2, 0],
                                                      rot[1, 0]-rot[0, 1]])
    return theta, direction



# Classes
class Plane(object): # TODO: placer cette classe dans un module a part (surface) --> utilise dans meshmagick aussi...
    """
    Class to handle planes for clipping purposes
    """
    def __init__(self, normal=[0., 0., 1.], scalar=0.):
        """
        Plane constructor

        Parameters
        ----------
        normal : array_like, optional
            Normal of the plane. Default to [0., 0., 1.]
        scalar : float, optional
            Plane scalar parameter. Default to 0.

        Returns
        -------
            Plane object
        """

        normal = np.asarray(normal, dtype=np.float)

        self._normal = normal / np.linalg.norm(normal)
        self._scalar = float(scalar)

        # Storing rotation matrix (redundant !) to speedup computations
        # Shall be _update in methods !!! --> using decorator ?
        thetax, thetay = self.get_normal_orientation_wrt_z()
        self._rot = _get_rotation_matrix(thetax, thetay)

        return

    @property
    def normal(self):
        """

        Returns
        -------
        ndarray
            The 3x1 plane's normal coordinates
        """
        return self._normal

    @normal.setter
    def normal(self, value):
        value = np.asarray(value, dtype=np.float)
        self._normal = value / np.linalg.norm(value)
        return

    @property
    def c(self):
        """

        Returns
        -------
        float
            The plane scalar
        """
        return self._scalar

    @c.setter
    def c(self, value):
        self._scalar = float(value)
        return

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

        Notes
        -----
        Equations are:
        .. math::
         \theta=\sqrt{\theta_x^2 + \theta_y^2}\\
         \sin{\theta} = \sqrt{n_x^2 + n_y^2}\\
         \theta_x = -\frac{\theta}{\sin{\theta}} n_y\\
         \theta_y =  \frac{\theta}{\sin{\theta}} n_x\\
         n_z = \cos{\theta}

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

    def set_plane_parameters(self, scalar, thetax, thetay):
        """
        Updates the plane parameters (normal and scalar parameter) given scalar and angles.

        Parameters
        ----------
        scalar : float
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
        self._scalar = self._scalar * ctheta + scalar
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
        return np.dot(points, self._normal) - self._scalar

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
        return -self._scalar * self.normal + np.transpose(np.dot(self._rot, points.T))

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
        t = (P0n - self._scalar) / (P0n - P1n)
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
    
    def get_origin(self):
        return self.c * self.normal

# class cached_property(object):
#     def __init__(self, func):
#         self.__doc__ = getattr(func, '__doc__')
#         self.func = func
#
#     def __get__(self, obj, cls):
#         try:
#             return obj._cached_properties[self.func.__name__]
#         except:
#             # print 'Computing %s and caching it' % self.func.__name__
#             value = self.func(obj)
#             try:
#                 obj._cached_properties[self.func.__name__] = value # FIXME: les valeurs sont enregistrees deux fois...
#             except AttributeError:
#                 obj._cached_properties = {self.func.__name__: value}
#             return value
#
#
# class invalidate_cache(object):
#     def __init__(self, func):
#         self.__doc__ = getattr(func, '__doc__')
#         self.func = func
#
#
#     def __call__(self, cls, *args):
#         self.func(cls, *args)
#         # print 'Invalidation of the cache'
#         try:
#             cls._cached_properties.clear()
#         except:
#             cls._cached_properties = dict()


class _3DPointsArray(np.ndarray):
    def __new__(cls, points):
        obj = np.asarray(points).view(cls)
        cls.x = property(fget=lambda cls: cls[:, 0])
        cls.y = property(fget=lambda cls: cls[:, 1])
        cls.z = property(fget=lambda cls: cls[:, 2])
        return obj


class Mesh(object):
    _ids = count(0)
    def __init__(self, vertices, faces, name=None):
        """

        Parameters
        ----------
        vertices : ndarray
            (nv x 3) Array of mesh _vertices coordinates. Each line is a vertex.
        faces : ndarray
            Arrays of mesh connectivities for _faces.

        Returns
        -------

        """

        self.__internals__ = dict()
        
        assert np.array(vertices).shape[1] == 3
        assert np.array(faces).shape[1] == 4
        
        self._vertices = np.array(vertices, dtype=np.float)
        self._faces = np.array(faces, dtype=np.int)
        self._id = self._ids.next()

        if not name:
            self._name = 'mesh_%u' % self._id
        else:
            self._name = str(name)

        self._verbose = False

    def __str__(self):
        str_repr = """
        --------------------------------------------
        \tMESH NAME : %s
        --------------------------------------------

        Number of _vertices: %u
        Number of _faces:    %u

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
               self._vertices[:, 0].min(),
               self._vertices[:, 0].max(),
               self._vertices[:, 1].min(),
               self._vertices[:, 1].max(),
               self._vertices[:, 2].min(),
               self._vertices[:, 2].max(),
               )
        return str_repr

    def print_quality(self):
        # This function is reproduced from
        # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/Verdict/Testing/Python/MeshQuality.py

        quality = vtk.vtkMeshQuality()
        quality.SetInput(self._vtk_polydata)

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

        if self._vtk_polydata.GetNumberOfCells() > 0:
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
        return self._vertices.shape[0]

    @property
    def nb_faces(self):
        return self._faces.shape[0]

    @property
    def vertices(self):
        return self._vertices

    @property
    def faces(self):
        return self._faces

    @vertices.setter
    def vertices(self, value):
        self._vertices = np.asarray(value, dtype=np.float).copy()
        self._vertices.setflags(write=False)
        self.__internals__.clear()
        return

    @faces.setter
    def faces(self, value):
        self._faces = np.asarray(value, dtype=np.int).copy()
        self._faces.setflags(write=False)
        self.__internals__.clear()
        return

    def _faces_properties(self):
        # faces_areas, faces_normals, faces_centers = mm.get_all_faces_properties(self._vertices, self._faces)
        nf = self.nb_faces

        # triangle_mask = _faces[:, 0] == _faces[:, -1]
        # nb_triangles = np.sum(triangle_mask)
        # quads_mask = np.invert(triangle_mask)
        # nb_quads = nf - nb_triangles

        faces_areas = np.zeros(nf, dtype=np.float)
        faces_normals = np.zeros((nf, 3), dtype=np.float)
        faces_centers = np.zeros((nf, 3), dtype=np.float)

        # Collectively dealing with triangles
        # triangles = _faces[triangle_mask]
        triangles_id = self.triangles_ids
        triangles = self._faces[triangles_id]

        triangles_normals = np.cross(self._vertices[triangles[:, 1]] - self._vertices[triangles[:, 0]],
                                     self._vertices[triangles[:, 2]] - self._vertices[triangles[:, 0]])
        triangles_areas = np.linalg.norm(triangles_normals, axis=1)
        faces_normals[triangles_id] = triangles_normals / np.array(([triangles_areas, ] * 3)).T
        faces_areas[triangles_id] = triangles_areas / 2.
        faces_centers[triangles_id] = np.sum(self._vertices[triangles[:, :3]], axis=1) / 3.

        # Collectively dealing with quads
        quads_id = self.quadrangles_ids
        quads = self._faces[quads_id]
        # quads = _faces[quads_mask]

        quads_normals = np.cross(self._vertices[quads[:, 2]] - self._vertices[quads[:, 0]],
                                 self._vertices[quads[:, 3]] - self._vertices[quads[:, 1]])
        faces_normals[quads_id] = quads_normals / np.array(([np.linalg.norm(quads_normals, axis=1), ] * 3)).T

        a1 = np.linalg.norm(np.cross(self._vertices[quads[:, 1]] - self._vertices[quads[:, 0]],
                                     self._vertices[quads[:, 2]] - self._vertices[quads[:, 0]]), axis=1) * 0.5
        a2 = np.linalg.norm(np.cross(self._vertices[quads[:, 3]] - self._vertices[quads[:, 0]],
                                     self._vertices[quads[:, 2]] - self._vertices[quads[:, 0]]), axis=1) * 0.5
        faces_areas[quads_id] = a1 + a2
        
        C1 = np.sum(self._vertices[quads[:, :3]], axis=1) / 3.
        C2 = (np.sum(self._vertices[quads[:, 2:4]], axis=1) + self._vertices[quads[:, 0]]) / 3.

        faces_centers[quads_id] = (np.array(([a1, ] * 3)).T * C1 + np.array(([a2, ] * 3)).T * C2)
        faces_centers[quads_id] /= np.array(([faces_areas[quads_id], ] * 3)).T

        faces_properties = {'faces_areas': faces_areas,
                            'faces_normals': faces_normals,
                            'faces_centers': faces_centers}

        self.__internals__.update(faces_properties)

        return

    def _has_faces_properties(self):
        return self.__internals__.has_key('faces_areas')

    def _remove_faces_properties(self):
        del self.__internals__['faces_areas']
        del self.__internals__['faces_centers']
        del self.__internals__['faces_normals']
        return

    @property
    def faces_areas(self):
        if not self.__internals__.has_key('faces_areas'):
            self._faces_properties()
        return self.__internals__['faces_areas']

    @property
    def faces_centers(self):
        if not self.__internals__.has_key('faces_centers'):
            self._faces_properties()
        return self.__internals__['faces_centers']

    @property
    def faces_normals(self):
        if not self.__internals__.has_key('faces_normals'):
            self._faces_properties()
        return self.__internals__['faces_normals']

    def _triangles_quadrangles(self):
        triangle_mask = (self._faces[:, 0] == self._faces[:, -1])
        quadrangles_mask = np.invert(triangle_mask)
        triangles_quadrangles = {'triangles_ids': np.where(triangle_mask)[0],
                                 'quadrangles_ids': np.where(quadrangles_mask)[0]}
        self.__internals__.update(triangles_quadrangles)
        return

    def _has_triangles_quadrangles(self):
        return self.__internals__.has_key('triangles_ids')

    def _remove_triangles_quadrangles(self):
        del self.__internals__['triangles_ids']
        del self.__internals__['quadrangles_ids']
        return

    @property
    def triangles_ids(self):
        if not self.__internals__.has_key('triangles_ids'):
            self._triangles_quadrangles()
        return self.__internals__['triangles_ids']

    @property
    def nb_triangles(self):
        if not self.__internals__.has_key('triangles_ids'):
            self._triangles_quadrangles()
        return len(self.__internals__['triangles_ids'])

    @property
    def quadrangles_ids(self):
        if not self.__internals__.has_key('triangles_ids'):
            self._triangles_quadrangles()
        return self.__internals__['quadrangles_ids']

    @property
    def nb_quadrangles(self):
        if not self.__internals__.has_key('triangles_ids'):
            self._triangles_quadrangles()
        return len(self.__internals__['quadrangles_ids'])

    def is_triangle(self, face_id):
        return self._faces[face_id, 0] == self._faces[face_id, -1]

    def get_face(self, face_id):
        if self.is_triangle(face_id):
            return self._faces[face_id, :3]
        else:
            return self._faces[face_id]

    def extract_faces(self, id_faces_to_extract, return_index=False):
        """
        Extracts a new mesh from a selection of _faces ids

        Parameters
        ----------
        id_faces_to_extract : ndarray

        Returns
        -------
        extracted_mesh : Mesh
        """
        nv = self.nb_vertices

        # Determination of the _vertices to keep
        Vmask = np.zeros(nv, dtype=bool)
        Vmask[self._faces[id_faces_to_extract].flatten()] = True
        idV = np.arange(nv)[Vmask]

        # Building up the vertex array
        V_extracted = self._vertices[idV]
        newID_V = np.arange(nv)
        newID_V[idV] = np.arange(len(idV))

        F_extracted = self._faces[id_faces_to_extract]
        F_extracted = newID_V[F_extracted.flatten()].reshape((len(id_faces_to_extract), 4))

        extracted_mesh = Mesh(V_extracted, F_extracted)
        extracted_mesh._verbose = self._verbose

        extracted_mesh.name = 'mesh_extracted_from_%s' % self.name

        if return_index:
            return extracted_mesh, idV
        else:
            return extracted_mesh

    def _vtk_polydata(self):
        # TODO: placer cette methode dans MMviewer !!
        # Create a vtkPoints object and store the points in it
        points = vtk.vtkPoints()
        for point in self._vertices:
            points.InsertNextPoint(point)

        # Create a vtkCellArray to store _faces
        faces = vtk.vtkCellArray()
        for face_ids in self._faces:
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
        vtk_polydata = self._vtk_polydata()
        self.viewer = MMviewer.MMViewer()
        self.viewer.add_polydata(vtk_polydata)
        self.viewer.show()
        self.viewer.finalize()

    def _connectivity(self):

        nv = self.nb_vertices
        nf = self.nb_faces

        mesh_closed = True

        # Building connectivities

        # Establishing VV and VF connectivities
        VV = dict([(i, set()) for i in xrange(nv)])
        VF = dict([(i, set()) for i in xrange(nv)])
        for (iface, face) in enumerate(self._faces):
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
                    boundary_face = self._faces[I[0]]

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
        # TODO: calculer des boundaries fermees et ouvertes (closed_boundaries et open_boundaries) et mettre dans dict
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
        self.__internals__.update(connectivity)

        return

    def _has_connectivity(self):
        return self.__internals__.has_key('VV')

    def _remove_connectivity(self):
        del self.__internals__['VV']
        del self.__internals__['VF']
        del self.__internals__['FF']
        del self.__internals__['boundaries']
        return

    @property
    def VV(self):
        if not self.__internals__.has_key('VV'):
            self._connectivity()
        return self.__internals__['VV']

    @property
    def VF(self):
        if not self.__internals__.has_key('VF'):
            self._connectivity()
        return self.__internals__['VF']

    @property
    def FF(self):
        if not self.__internals__.has_key('FF'):
            self._connectivity()
        return self.__internals__['FF']

    @property
    def boundaries(self):
        if not self.__internals__.has_key('boundaries'):
            self._connectivity()
        return self.__internals__['boundaries']

    @property
    def nb_boundaries(self):
        if not self.__internals__.has_key('boundaries'):
            self._connectivity()
        return len(self.__internals__['boundaries'])
    
    @property
    def axis_aligned_bbox(self):
        x, y, z = self._vertices.T
        return (x.min(), x.max(),
                y.min(), y.max(),
                z.min(), z.max())
    
    @property
    def squared_axis_aligned_bbox(self):
        xmin, xmax, ymin, ymax, zmin, zmax = self.axis_aligned_bbox
        (x0, y0, z0) = np.array([xmin+xmax, ymin+ymax, zmin+zmax]) * 0.5
        d = (np.array([xmax-xmin, ymax-ymin, zmax-zmin]) * 0.5).max()
        return (x0-d, x0+d, y0-d, y0+d, z0-d, z0+d)
    
    def is_mesh_closed(self):
        if not self.__internals__.has_key('boundaries'):
            self._connectivity()
        return len(self.__internals__['boundaries']) == 0

    def is_mesh_conformal(self):

        warn('This method is not stable yet !! Use with caution')
        # FIXME: bugge
        tol = 1e-7

        boundaries = self._connectivity['boundaries']
        polygons_areas = np.zeros(len(boundaries), dtype=np.float)

        conformal = True

        for boundary in boundaries:
            boundary_vertices = self._vertices[boundary]
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
        return self.rotate([thetax, 0., 0.])

    def rotate_y(self, thetay):
        return self.rotate([0., thetay, 0.])

    def rotate_z(self, thetaz):
        return self.rotate([0., 0., thetaz])

    def rotate(self, angles):
        angles = np.asarray(angles, dtype=np.float)
        theta = np.linalg.norm(angles)
        if theta == 0.:
            return np.eye(3)

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

        self._vertices = np.transpose(np.dot(R, self._vertices.copy().T))

        # Updating _faces properties if any
        if self._has_faces_properties():
            # Rotating normals and centers too
            normals = self.__internals__['faces_normals']
            centers = self.__internals__['faces_centers']
            self.__internals__['faces_normals'] = np.transpose(np.dot(R, normals.T))
            self.__internals__['faces_centers'] = np.transpose(np.dot(R, centers.T))

        return R

    # @invalidate_cache
    def translate_x(self, tx):
        V = self._vertices
        V[:, 0] += tx
        self._vertices = V

        # Updating properties if any
        if self._has_faces_properties():
            centers = self.__internals__['faces_centers']
            centers[:, 0] += tx
            self.__internals__['faces_centers'] = centers

        return

    def translate_y(self, ty):
        V = self._vertices.copy()
        V[:, 1] += ty
        self._vertices = V

        # Updating properties if any
        if self._has_faces_properties():
            centers = self.__internals__['faces_centers']
            centers[:, 1] += ty
            self.__internals__['faces_centers'] = centers
        return

    def translate_z(self, tz):
        V = self._vertices.copy()
        V[:, 2] += tz
        self._vertices = V

        # Updating properties if any
        if self._has_faces_properties():
            centers = self.__internals__['faces_centers']
            centers[:, 2] += tz
            self.__internals__['faces_centers'] = centers
        return

    def translate(self, t):
        # t = np.asarray(t, dtype=np.float)
        tx, ty, tz = t
        V = self._vertices.copy() # FIXME: why doing a copy ???
        V[:, 0] += tx
        V[:, 1] += ty
        V[:, 2] += tz
        self._vertices = V

        # Updating properties if any
        if self._has_faces_properties():
            centers = self.__internals__['faces_centers']
            centers[:, 0] += tx
            centers[:, 1] += ty
            centers[:, 2] += tz
            self.__internals__['faces_centers'] = centers
        return

    def scale(self, alpha):
        V = self._vertices.copy()
        V *= alpha
        self._vertices = V

        if self._has_faces_properties():
            self._remove_faces_properties()

        return

    def flip_normals(self):
        F = self._faces.copy()
        self._faces = np.fliplr(F)

        if self._has_faces_properties():
            self.__internals__['faces_normals'] *= -1

        return

    def __add__(self, mesh_to_add):
        V = np.concatenate((self._vertices, mesh_to_add._vertices), axis=0)
        F = np.concatenate((self._faces, mesh_to_add._faces + self.nb_vertices), axis=0)
        new_mesh = Mesh(V, F, name='_'.join([self.name, mesh_to_add.name]))
        # new_mesh.merge_duplicates()
        new_mesh._verbose = self._verbose or mesh_to_add._verbose

        # TODO: exporter les pptes de facette si elles sont presentes dans les 2 maillages
        # if hasattr(self, '_cached_properties') and hasattr(mesh_to_add, '_cached_properties'):
        #     if self._cached_properties.has_key('_faces_properties') and mesh_to_add._cached_properties.has_key('_faces_properties'):
        #         new_mesh._cached_properties['_faces_properties']['areas'] = \
        #             np.concatenate(self.faces_areas, mesh_to_add.faces_areas)
        #         new_mesh._cached_properties['_faces_properties']['centers'] = \
        #             np.concatenate(self.faces_centers, mesh_to_add.faces_centers)
        #         new_mesh._cached_properties['_faces_properties']['normals'] = \
        #             np.concatenate(self.faces_normals, mesh_to_add.faces_normals)

        return new_mesh

    def copy(self):
        return copy.deepcopy(self)

    # def merge_duplicates(self, tol=1e-8, return_index=False):
    #     # TODO: voir ou mettre l'implementation de la fonction merge_duplicates
    #     output = mm.merge_duplicates(self._vertices, self._faces, verbose=False, tol=tol, return_index=return_index)
    #     if return_index:
    #         V, F, newID = output
    #     else:
    #         V, F = output
    #     if self._verbose:
    #         print "* Merging duplicate vertices"
    #         delta_n = self.nb_vertices - V.shape[0]
    #         if delta_n > 0:
    #             print "\t--> %u vertices have been merged" % delta_n
    #         else:
    #             print "\t--> No duplicate vertices have been found"
    #     self._vertices, self._faces = V, F
    #
    #     if self._has_connectivity():
    #         self._remove_connectivity()
    #
    #     if return_index:
    #         return newID
    #     else:
    #         return

    def merge_duplicates(self, decimals=8, return_index=False):
        uniq, newID = merge_duplicate_rows(self._vertices, decimals=decimals, return_index=True)

        nv_init = self.nb_vertices

        # Updating mesh data
        self._vertices = uniq
        self._faces = newID[self._faces] # Faces vertices ids are updated here

        nv_final = self.nb_vertices

        if self._verbose:
            print "* Merging duplicate vertices that are close to %u decimals..." % decimals
            delta_n = nv_init - nv_final
            if delta_n == 0:
                print "\t--> No duplicate vertices have been found"
            else:
                print "\t--> Initial number of vertices : %u" % nv_init
                print "\t--> Final number of vertices   : %u" % nv_final
                print "\t--> %u vertices have been merged\n" % delta_n

        if self._has_connectivity():
            self._remove_connectivity()

        if return_index:
            return newID
        else:
            return

    def heal_normals(self):

        # TODO: return the different groups of a mesh in case it is made of several unrelated groups

        nv = self.nb_vertices
        nf = self.nb_faces
        F = self._faces

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
        # triangles_mask = _faces[:, 0] == _faces[:, -1]
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

                # Shared _vertices
                adjface = F[iadjF]
                S2 = set(adjface)
                # try:
                common_vertices = list(S1 & S2)
                if len(common_vertices) == 2:
                    iV1, iV2 = common_vertices
                else:
                    print 'WARNING: _faces %u and %u have more than 2 _vertices in common !' % (iface, iadjF)
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
                print '\t--> %u _faces have been reversed to make normals consistent across the mesh' % (nb_reversed)
            else:
                print "\t--> Normals orientations are consistent"

        self._faces = F

        # Checking if the normals are outward
        if mesh_closed:
            zmax = np.max(self._vertices[:, 2])

            areas = self.faces_areas
            normals = self.faces_normals
            centers = self.faces_centers
            # areas, normals, centers = get_all_faces_properties(_vertices, _faces)

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

        if self._has_faces_properties():
            self._remove_faces_properties()

        return

    def remove_unused_vertices(self, return_index=False):
        # TODO: implementer return_index !!
        nv = self.nb_vertices
        V, F = self._vertices, self._faces

        usedV = np.zeros(nv, dtype=np.bool)
        usedV[sum(map(list, F), [])] = True
        nb_usedV = sum(usedV)

        if nb_usedV < nv:
            newID_V = np.arange(nv)
            newID_V[usedV] = np.arange(nb_usedV)
            F = newID_V[F]
            V = V[usedV]

        self._vertices, self._faces = V, F

        if self._verbose:
            print "* Removing unused _vertices in the mesh:"
            if nb_usedV < nv:
                unusedV = np.where(np.logical_not(usedV))[0]
                vlist_str = '[' + ', '.join(str(iV) for iV in unusedV) + ']'
                print "\t--> %u unused _vertices have been removed" % (nv - nb_usedV)
            else:
                print "\t--> No unused _vertices"

        if self._has_connectivity():
            self._remove_connectivity()

        return

    def heal_triangles(self):

        F = self._faces

        quads = F[:, 0] != F[:, -1]
        nquads_init = sum(quads)

        F[quads] = np.roll(F[quads], 1, axis=1)
        quads = F[:, 0] != F[:, -1]

        F[quads] = np.roll(F[quads], 1, axis=1)
        quads = F[:, 0] != F[:, -1]

        F[quads] = np.roll(F[quads], 1, axis=1)
        quads = F[:, 0] != F[:, -1]
        nquads_final = sum(quads)

        self._faces = F

        if self._verbose:
            print "* Ensuring consistent definition of triangles:"
            if nquads_final < nquads_init:
                print "\t--> %u triangles were described the wrong way and have been corrected" % (
                nquads_init - nquads_final)
            else:
                print "\t--> Triangle description is consistent"

        return

    def remove_degenerated_faces(self, rtol=1e-5):
        # TODO: implementer un retour d'index des _faces extraites
        areas = self.faces_areas
        area_threshold = areas.mean() * rtol

        # Detecting _faces that have null area
        F = self._faces[np.logical_not(areas < area_threshold)]
        if self._verbose:
            nb_removed = self.nb_faces - F.shape[0]
            print '* Removing degenerated _faces'
            if nb_removed > 0:
                print '\t-->%u degenerated _faces have been removed' % nb_removed
            else:
                print '\t--> No degenerated _faces'

        self._faces = F
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

        F = self._faces

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

        self.__internals__.clear()

        self._faces = F

        return F

    def symmetrize(self, plane):

        # Symmetrizing the nodes
        V, F = self._vertices, self._faces

        V = np.concatenate((V, V - 2 * np.outer(np.dot(V, plane.normal) - plane.c, plane.normal)))
        F = np.concatenate((F, np.fliplr(F.copy() + self.nb_vertices)))

        self._vertices, self._faces = V, F
        verbose = self.verbose
        self.verbose_off()
        self.merge_duplicates()
        self.verbose = verbose
        return
    
    def mirror(self, plane):
        self._vertices = self._vertices - 2 * np.outer(np.dot(self._vertices, plane.normal) - plane.c, plane.normal)
        self.flip_normals()
        return
    
    def _compute_faces_integrals(self, sum_faces_contrib=False): # TODO: implementer le sum_surface_contrib

        # TODO: faire la machinerie pour pour
        surface_integrals = np.zeros((15, self.nb_faces), dtype=np.float)

        # First triangles
        if self.nb_triangles > 0:
            triangles_ids = self.triangles_ids
            # print self._faces[triangles_ids][:, :3].shape
            triangles_vertices = self._vertices[self._faces[triangles_ids][:, :3]] # Remettre le 3
            surface_integrals[:, triangles_ids] = self._compute_triangles_integrals(triangles_vertices)

        # Now quadrangles by splitting them up
        if self.nb_quadrangles > 0:
            quadrangles_ids = self.quadrangles_ids
            quadrangles = self._faces[quadrangles_ids]

            # First pass
            surface_integrals[:, quadrangles_ids] = \
                self._compute_triangles_integrals(self._vertices[quadrangles[:, (0, 1, 2)]])

            # Second pass
            surface_integrals[:, quadrangles_ids] += \
                self._compute_triangles_integrals(self._vertices[quadrangles[:, (0, 2, 3)]])

        # names = ('sint_x', 'sint_y', 'sint_z',
        #          'sint_yz', 'sint_xz', 'sint_xy',
        #          'sint_x2', 'sint_y2', 'sint_z2',
        #          'sint_x3', 'sint_y3', 'sint_z3'
        #          'sint_x2y', 'sint_y2z', 'sint_z2x')
        #
        # s_int = dict(zip(names, surface_integrals))
        # self.__internals__._update_hydrostatic_properties(s_int)

        self.__internals__['surface_integrals'] = surface_integrals

        return

    def has_surface_integrals(self):
        return self.__internals__.has_key('surface_integrals')

    def get_surface_integrals(self):
        return self.__internals__['surface_integrals']

    def _compute_volume(self):
        if not self.has_surface_integrals():
            self._compute_faces_integrals()

        normals = self.faces_normals
        sigma_0_2 = self.__internals__['surface_integrals'][:3]

        return (normals.T * sigma_0_2).sum() / 3.

    @property
    def volume(self):
        return self._compute_volume()
    
    def _edges_stats(self):
        pass
    
    @property
    def min_edge_length(self):
        pass
    
    @staticmethod
    def _compute_triangles_integrals(triangles_vertices, sum_faces_contrib=False):
        """
        Notes
        -----
        triangles_vertices doit decrire par dimension croissante du general au particulier :
        dimension 0 : informations sur chaque facette -- triangles_vertices[0] -> facette 0)
        dimension 1 : informations sur chaque vertex de la facette -- triangles_vertices[0, 1] -> vertex 1 de la facette 0
        dimension 2 : information sur chacune des coordonnées des vertex -- triangles_vertices[0, 1, 2] -> coordonnee z du vertex 1 de la facette 0
        """


        s_int = np.zeros((15, triangles_vertices.shape[0]), dtype=np.float)

        P0, P1, P2 = map(_3DPointsArray, np.rollaxis(triangles_vertices, 1, 0))

        t0 = P0 + P1
        f1 = t0 + P2
        t1 = P0 * P0
        t2 = t1 + P1*t0
        f2 = t2 + P2*f1
        f3 = P0*t1 + P1*t2 + P2*f2
        g0 = f2 + P0 * (f1+P0)
        g1 = f2 + P1 * (f1+P1)
        g2 = f2 + P2 * (f1+P2)

        e1 = P1 - P0
        e2 = P2 - P0

        delta = np.linalg.norm(np.cross(e1, e2), axis=1)

        s_int[0:3] = np.einsum('i, ij -> ji', delta, f1) / 6.

        s_int[3] = delta * (6.*P0.y*P0.z + 3*(P1.y*P1.z + P2.y*P2.z) - P0.y*f1[:, 2] - P0.z*f1[:, 1]) / 12.
        s_int[4] = delta * (6.*P0.x*P0.z + 3*(P1.x*P1.z + P2.x*P2.z) - P0.x*f1[:, 2] - P0.z*f1[:, 0]) / 12.
        s_int[5] = delta * (6.*P0.x*P0.y + 3*(P1.x*P1.y + P2.x*P2.y) - P0.x*f1[:, 1] - P0.y*f1[:, 0]) / 12.

        s_int[6:9] = np.einsum('i, ij -> ji', delta, f2) / 12.
        s_int[9:12] = np.einsum('i, ij -> ji', delta, f3) / 20.

        # Ne pas oublier le delta
        s_int[12] = delta * ( P0.y*g0[:, 0] + P1.y*g1[:, 0] + P2.y*g2[:, 0]) / 60.
        s_int[13] = delta * ( P0.z*g0[:, 1] + P1.z*g1[:, 1] + P2.z*g2[:, 1]) / 60.
        s_int[14] = delta * ( P0.x*g0[:, 2] + P1.x*g1[:, 2] + P2.x*g2[:, 2]) / 60.

        return s_int
    
    def quick_save(self, filename=None):
        if not filename:
            filename = 'quick_save.vtp'
        
        if not filename.endswith('.vtp'):
            filename += '.vtp'
        try:
            import mmio
            mmio.write_VTP(filename, self.vertices, self.faces)
            print 'File %s written' % filename
        except:
            print 'mmio module not found'


if __name__ == '__main__':

    import mmio
    V, F = mmio.load_VTP('SEAREV.vtp')
    mymesh = Mesh(V, F)

    mymesh._compute_faces_integrals()

    print mymesh.volume

    sys.exit(0)

