#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
This module is part of meshmagick software. It provides functions for mesh clipping purpose.
"""

import meshmagick as mm # TODO: voir si oblige d'importer meshmagick ici ...
import MMviewer
import numpy as np
import math
import vtk
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
        self.c = c

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
            print 'Computing %s and caching it' % self.func.__name__
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
        print 'Invalidation of the cache'
        try:
            cls._cached_properties.clear()
        except:
            cls._cached_properties = dict()


# TODO: cette classe devra se trouver dans meshmagick...

class Mesh(object):

    def __init__(self, vertices, faces):
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

    # TODO: implementer la fonction directement dans la classe mais la splitter pour chaque propriete areas,normals,
    # centers...

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

    def extract_faces(self, id_faces_to_extract):
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

        return Mesh(V_extracted, F_extracted)

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






class MeshClipper(object):
    def __init__(self, mesh=None, clipping_surface=None):
        self._mesh = mesh
        self._clipping_surface = clipping_surface

    @property
    def mesh(self):
        return self._mesh

    @property
    def clipping_surface(self):
        return self._clipping_surface

    @mesh.setter
    @invalidate_cache
    def mesh(self, obj):
        self.mesh = obj
        return

    @clipping_surface.setter
    @invalidate_cache
    def clipping_surface(self, obj):
        self.clipping_surface = obj
        return

    # def partition_mesh(self):
    #     #TODO: split_mesh doit devenir une methode de classe MeshClipper
    #     print split_mesh(self.mesh.V, self.mesh.F, self.clipping_surface)




def split_mesh(V, F, plane):
    """
    Splits the mesh into an intersection crowm and a wet part.

    Parameters
    ----------
    V : ndarray
        Array of mesh vertices coordinates
    F : ndarray
        Array of mesh connectivities
    plane : Plane
        Clipping plane. Vertices that are at the opposite from the normal pointing are considered under

    Returns
    -------
    crown_face_ids : ndarray
        Array of face ids that are intersecting the plane (the intersection crown)
    below_face_ids : ndarray
        Array of face ids that are totally under the plane
    """
    nf = F.shape[0]

    # Triangles and quadrangles masks
    triangle_mask = F[:, 0] == F[:, -1]
    quad_mask = np.invert(triangle_mask)

    triangles = np.where(triangle_mask)[0]
    quadrangles = np.where(quad_mask)[0]

    n_tri = triangle_mask.sum()
    n_quad =  nf - n_tri

    # Construction des half-edges
    # TODO: Ne pas le refaire a chaque fois... Mettre en pre-processing d'un classe
    nhe = 3*n_tri + 4*n_quad
    HE_F = np.zeros(nhe, dtype=np.int)
    HE = np.zeros((nhe, 2), dtype=np.int)

    for itri, triangle in enumerate(F[triangle_mask, :3]):
        istart = 3*itri
        istop = istart+3
        HE[istart:istop, :] = [(triangle[i-1], triangle[i]) for i in xrange(3)]
        HE_F[istart:istop] = triangles[itri]

    itri_offset = 3*n_tri
    for iquad, quadrangle in enumerate(F[quad_mask, :]):
        istart = 4*iquad + itri_offset
        istop = istart+4
        HE[istart:istop, :] = [(quadrangle[i-1], quadrangle[i]) for i in xrange(4)]
        HE_F[istart:istop] = quadrangles[iquad]


    # Tout ce qui au-dessus doit etre mis dans une classe Mesh

    # Position of points with respect to the clipping plane
    pos = plane.get_point_dist_wrt_plane(V)

    v_below = pos < 0.

    nb_below_by_tri =  v_below[F[triangles, :3]].sum(axis=1)
    nb_below_by_quad = v_below[F[quadrangles]].sum(axis=1)

    clipping_triangles = triangles[np.logical_and(nb_below_by_tri > 0, nb_below_by_tri < 3)]
    clipping_quadrangles = quadrangles[np.where(np.logical_and(nb_below_by_quad > 0, nb_below_by_quad < 4))[0]]
    crown_face_ids = np.concatenate((clipping_triangles, clipping_quadrangles))

    below_triangles = triangles[nb_below_by_tri == 3]
    below_quadrangles = quadrangles[nb_below_by_quad == 4]
    below_face_ids = np.concatenate((below_triangles, below_quadrangles))

    return crown_face_ids, below_face_ids


def _clip(V, F, plane, tol=1e-4, pos=None, return_boundaries=False, assert_closed_boundaries=True):
    """
    Returns a clipped mesh by a plane

    Parameters
    ----------
    V : ndarray
        Array of mesh vertices coordinates
    F : ndarray
        Array of mesh connectivities
    plane : Plane
        The clipping plane
    tol : float, optional
        Absolute tolerance for the vertices to be projected on the plane. Default to 1e-4.
    pos : ndarray, optional
        Array of distance of vertices with respect to the plane. Default to None trigs its computations. Useful if
        these distances have already been calculated.
    return_boundaries : bool
        Wether to return lists of intersection boundaries. Default to False
    assert_closed_boundaries : bool
        If the mesh is closed, then every intersection curves must be a closed polygon. If True (default),
        then an open boundary line will raise a Runtime Exception. Only used if return_boundaries is set to True.

    Returns
    -------
    Vclip : ndarray
        Array of vertices coordinates of the clipped mesh
    F_clip : ndarray
        Array of connectivities of the clipped mesh
    boundary_polygons : list, optional
        List of closed boundary polygons, oriented anticlockwise
    boundary_lines : list, optional
        List of open boundary lines, oriented anticlockwise

    Notes
    -----
    The algorithms classifies each face to clip by counting the number of vertices that are above, on and below the
    plane. Begin on the plane is to be at a vicinity of tol to the plane. These vertices are automatically projected
    on the plane, orthogonally. The three numbers are combined into a hash string and used as an identifier to deal
    with each possible case specifically so as to bring total robustness of algorithm as well as quich calculations.
    """

    # TODO: voir a recuperer ces positions en entree... (seulement pour la crown !
    Vclip = V.copy()
    nv = Vclip.shape[0]

    if pos is None:
        pos = plane.get_point_dist_wrt_plane(Vclip) # TODO: voir a recuperer pos en entree

    # Vertices pre-projection
    v_on = np.where(np.fabs(pos) <= tol)[0]
    V[v_on] = plane.orthogonal_projection_on_plane(V[v_on])
    pos[v_on] = 0.

    on_mask = np.zeros(nv, dtype=np.bool)
    on_mask[v_on] = True

    # Classifying vertices
    above_mask = pos > tol
    below_mask = pos < -tol

    v_above = np.where(above_mask)[0]
    # v_on = np.where(on_mask)[0]
    v_below = np.where(below_mask)[0]

    # Init
    F_clip = list()
    direct_boundary_edges = dict()
    inv_boundary_edges = dict()
    intersections = list()
    nI = nv

    for iface, face in enumerate(F):

        if face[0] == face[-1]:
            face = face[:3]

        nv_face = len(face)

        face_pos = pos[face]

        # Determining the type of face clipping
        v_above_face = np.where(above_mask[face])[0]
        v_on_face = np.where(on_mask[face])[0]
        v_below_face = np.where(below_mask[face])[0]

        nb_above = len(v_above_face)
        nb_on = len(v_on_face)
        nb_below = len(v_below_face)

        face_type = str(nb_above)+str(nb_on)+str(nb_below)

        if face_type == '202': # Done
            #    0*-----*3
            #     |     |
            # ----o-----o----
            #     |     |
            #    1*-----*2
            if v_above_face[1] == v_above_face[0]+1:
                face = np.roll(face, -v_above_face[1])
            P0, P1, P2, P3 = Vclip[face]
            Ileft = plane.get_edge_intersection(P0, P1)
            Iright = plane.get_edge_intersection(P2, P3)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI, face[1], face[2], nI+1])
            nI += 2

        elif face_type == '301':#Done
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
            P0, P1, P3 = Vclip[face[[0, 1, 3]]]
            Ileft = plane.get_edge_intersection(P0, P3)
            Iright = plane.get_edge_intersection(P0, P1)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI, face[0], nI+1, nI])
            nI +=2

        elif face_type == '103': # Done
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
            P0, P1, P3 = Vclip[face[[0, 1, 3]]]
            Ileft = plane.get_edge_intersection(P0, P1)
            Iright = plane.get_edge_intersection(P0, P3)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI, face[1], face[3], nI+1])
            F_clip.append([face[1], face[2], face[3], face[1]])
            nI += 2

        elif face_type == '102': #Done
            #      *O
            #     / \
            # ---o---o---
            #   /     \
            # 1*-------*2
            face = np.roll(face, -v_above_face[0])
            P0, P1, P2 = Vclip[face]
            Ileft = plane.get_edge_intersection(P0, P1)
            Iright = plane.get_edge_intersection(P0, P2)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI,  face[1], face[2], nI+1])
            nI += 2

        elif face_type == '201':#done
            #  2*-------*1
            #    \     /
            #  ---o---o---
            #      \ /
            #       *0
            face = np.roll(face, -v_below_face[0])
            P0, P1, P2 = Vclip[face]
            Ileft = plane.get_edge_intersection(P0, P2)
            Iright = plane.get_edge_intersection(P0, P1)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI, face[0], nI+1, nI])
            nI += 2

        elif face_type == '211': # Done
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
                P1, P2 = Vclip[face[[1, 2]]]
                Iright = plane.get_edge_intersection(P1, P2)
                intersections.append(Iright)
                boundary_edge = [face[0], nI]
                F_clip.append([face[0], face[1], nI, face[0]])
            else:
                P2, P3 = Vclip[face[[2, 3]]]
                Ileft = plane.get_edge_intersection(P2, P3)
                intersections.append(Ileft)
                boundary_edge = [nI, face[0]]
                F_clip.append([nI, face[3], face[0], nI])
            nI += 1

        elif face_type == '112': # Done
            #       *3                     *1
            #      / \                    / \
            #  ---*---o---      or    ---o---*---
            #     0\   \                /   /0
            #       \   *2            2*   /
            #        \ /                \ /
            #         *1                 *3
            face = np.roll(face, -v_on_face[0])
            if pos[face[1]] < 0.:
                P2, P3 = Vclip[face[[2, 3]]]
                Iright = plane.get_edge_intersection(P2, P3)
                intersections.append(Iright)
                boundary_edge = [face[0], nI]
                F_clip.append([face[0], face[1], face[2], nI])
            else:
                P1, P2 = Vclip[face[[1, 2]]]
                Ileft = plane.get_edge_intersection(P1, P2)
                intersections.append(Ileft)
                boundary_edge = [nI, face[0]]
                F_clip.append([nI, face[2], face[3], face[0]])
            nI += 1

        elif face_type == '013': # Done
            # -----*-----
            #     / \
            #    /   \
            #   *     *
            #    \   /
            #     \ /
            #      *
            boundary_edge = None
            F_clip.append(list(face))

        elif face_type == '210' or face_type == '310': #Done
            #   *-------*               *
            #    \ 210 /               / \ 310
            #     \   /               *   *
            #      \ /                 \ /
            #   ----*----           ----*----
            boundary_edge = None

        elif face_type == '111': #Done
            #        *2              *1
            #       /|               |\
            #      / |               | \
            #  ---*--o---    or   ---o--*---
            #     0\ |               | /0
            #       \|               |/
            #        *1              *2
            face = np.roll(face, -v_on_face[0])
            P1, P2 = Vclip[face[[1, 2]]]
            if pos[face[1]] < 0.:
                Iright = plane.get_edge_intersection(P1, P2)
                intersections.append(Iright)
                boundary_edge = [face[0], nI]
                F_clip.append([face[0], face[1], nI, face[0]])
            else:
                Ileft = plane.get_edge_intersection(P1, P2)
                intersections.append(Ileft)
                boundary_edge = [nI, face[0]]
                F_clip.append([nI, face[2], face[0], nI])
            nI += 1

        elif face_type == '120': # Done
            #         *O
            #        / \
            #       /   \
            #     1/     \2
            # ----*-------*----
            face = np.roll(face, -v_above_face[0])
            boundary_edge = [face[1], face[2]]

        elif face_type == '021': # Done
            #  ----*-------*----
            #      2\     /1
            #        \   /
            #         \ /
            #          *0
            face = np.roll(face, -v_below_face[0])
            boundary_edge = [face[2], face[1]]
            face = list(face)
            face.append(face[0])
            F_clip.append(face)

        elif face_type == '022':
            # ----*-----*----
            #    0|     |3
            #     |     |
            #    1*-----*2
            if v_on_face[1] == v_on_face[0]+1:
                face = np.roll(face, -v_on_face[1])
            boundary_edge = [face[0], face[3]]
            F_clip.append(list(face))

        elif face_type == '012': # Done
            #   ------*------
            #        / \
            #       /   \
            #      /     \
            #     *-------*
            boundary_edge = None
            face = list(face)
            face.append(face[0])
            F_clip.append(face)

        elif face_type == '220': # Done
            #    0*-----*3
            #     |     |
            #    1|     |2
            # ----*-----*----
            if v_above_face[1] == v_above_face[0]+1:
                face = np.roll(face, -v_above_face[1])
            boundary_edge = [face[1], face[2]]

        elif face_type == '121': # Done
            #       *0
            #      / \
            #     /   \
            # ---*-----*---
            #    1\   /3
            #      \ /
            #       *2
            face = np.roll(face, -v_above_face[0])
            boundary_edge = [face[1], face[3]]
            F_clip.append([face[1], face[2], face[3], face[1]])

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
            F_clip.append(face)

        elif face_type == '004':
            #  ---------------
            #      *-----*
            #      |     |
            #      |     |
            #      *-----*
            boundary_edge = None
            F_clip.append(list(face))

        # Building boundary connectivity
        if boundary_edge is not None:
            direct_boundary_edges[boundary_edge[0]] = boundary_edge[1]
            inv_boundary_edges[boundary_edge[1]] = boundary_edge[0]

    # It is now necessary to merge intersection vertices in order to build the intersection polygons
    intersections, newID = mm.merge_duplicates(np.asarray(intersections), return_index=True, tol=0.1*tol)

    newID = np.concatenate((np.arange(nv), newID + nv))
    F_clip = newID[F_clip]

    Vclip = np.concatenate((Vclip, intersections))

    if return_boundaries:
        # Updating dictionaries
        direct_boundary_edges = dict(zip(newID[direct_boundary_edges.keys()], newID[direct_boundary_edges.values()]))
        inv_boundary_edges = dict(zip(newID[inv_boundary_edges.keys()], newID[inv_boundary_edges.values()]))

        # Ordering boundary edges in continuous lines
        boundary_polygons = list()
        boundary_lines = list()
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
                                    boundary_lines.append(line)
                                    break
                        else:
                            boundary_polygons.append(line)

                        break

            except:
                # TODO: retirer les deux lines suivantes
                print "%u closed polygon\n%u open curve" % (len(boundary_polygons), len(boundary_lines))
                # print boundary_lines[0]

                if assert_closed_boundaries:
                    if len(boundary_lines) > 0:
                        mm.write_VTP('mesh_clip.vtp', Vclip, F_clip)
                        raise RuntimeError, 'Open intersection curve found'

                break

        return Vclip, F_clip, boundary_polygons, boundary_lines

    else: # return_boundaries = False
        # mm.show(Vclip, F_clip)
        return Vclip, F_clip

# def clip_mesh_against_plane(V, F):
#     pass

if __name__ == '__main__':

    V, F = mm.load_VTP('SEAREV.vtp')

    V, F = mm.clip_by_plane(V, F, mm.Plane())

    mesh = Mesh(V, F)

    print mesh.is_mesh_conformal()
    # print mesh._connectivity

    # mesh_extract = mesh.extract_faces([1, 2, 3])
    # mm.show(mesh_extract.V, mesh_extract.F)


    # plane = Plane()
    # clipper = MeshClipper(mesh=mesh, clipping_surface=plane)

    # clipper.partition_mesh()



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
    # crown_face_ids, below_face_ids = split_mesh(V, F, plane)
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
