#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
This module is part of meshmagick software. It provides functions for mesh clipping purpose.
"""

import meshmagick as mm # TODO: voir si oblige d'importer meshmagick ici ...
import numpy as np
import math
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

        self.normal = normal / np.linalg.norm(normal)
        self.c = c

        # Storing rotation matrix (redundant !) to speedup computations
        # Shall be update in methods !!! --> using decorator ?
        thetax, thetay = self.get_normal_orientation_wrt_z()
        self._rot = _get_rotation_matrix(thetax, thetay)

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

    def orthogonal_projection_on_plane(self, point):
        """
        Returns the coordinates of the orthogonal projection of point

        Parameters
        ----------
        point : ndarray
            Coordinates of the point to be projected

        Returns
        -------
        I : ndarray
            Coordinates of the projection point
        """
        # TODO: passer en vectoriel

        return point - self.get_point_dist_wrt_plane(point) * self.normal


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

        # Determining the type of face
        above_mask = face_pos > tol
        on_mask = np.fabs(face_pos) <= tol
        below_mask = face_pos < -tol

        v_above = np.where(face_pos > tol)[0]
        v_on = np.where(np.fabs(face_pos) <= tol)[0]
        v_below = np.where(face_pos < -tol)[0]

        nb_above = above_mask.sum()
        nb_on = on_mask.sum()
        nb_below = below_mask.sum()

        face_type = str(nb_above)+str(nb_on)+str(nb_below)

        if face_type == '202': # Done
            if v_above[1] == v_above[0]+1:
                face = np.roll(face, -v_above[1])
            P0, P1, P2, P3 = Vclip[face]
            Ileft = plane.get_edge_intersection(P0, P1)
            Iright = plane.get_edge_intersection(P2, P3)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI, face[1], face[2], nI+1])
            nI += 2

        elif face_type == '301':#Done
            face = np.roll(face, -v_below[0])
            P0, P1, P3 = Vclip[face[[0, 1, 3]]]
            Ileft = plane.get_edge_intersection(P0, P3)
            Iright = plane.get_edge_intersection(P0, P1)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI, face[0], nI+1, nI])
            nI +=2

        elif face_type == '103': # Done
            face = np.roll(face, -v_above[0])
            P0, P1, P3 = Vclip[face[[0, 1, 3]]]
            Ileft = plane.get_edge_intersection(P0, P1)
            Iright = plane.get_edge_intersection(P0, P3)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI, face[1], face[3], nI+1])
            F_clip.append([face[1], face[2], face[3], face[1]])
            nI += 2

        elif face_type == '102': #Done
            face = np.roll(face, -v_above[0])
            P0, P1, P2 = Vclip[face]
            Ileft = plane.get_edge_intersection(P0, P1)
            Iright = plane.get_edge_intersection(P0, P2)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI,  face[1], face[2], nI+1])
            nI += 2

        elif face_type == '201':#done
            face = np.roll(face, -v_below[0])
            P0, P1, P2 = Vclip[face]
            Ileft = plane.get_edge_intersection(P0, P2)
            Iright = plane.get_edge_intersection(P0, P1)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_clip.append([nI, face[0], nI+1, nI])
            nI += 2

        elif face_type == '211': # Done
            face = np.roll(face, -v_on[0])
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
            Vclip[face[0]] = plane.orthogonal_projection_on_plane(Vclip[face[0]])
            nI += 1

        elif face_type == '112': # Done
            face = np.roll(face, -v_on[0])
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
            Vclip[face[0]] = plane.orthogonal_projection_on_plane(Vclip[face[0]])
            nI += 1

        elif face_type == '013': # Done
            Vclip[face[v_on[0]]] = plane.orthogonal_projection_on_plane(Vclip[face[v_on[0]]])
            boundary_edge = []
            F_clip.append(list(face))

        elif face_type == '210' or face_type == '310': #Done
            Vclip[face[v_on[0]]] = plane.orthogonal_projection_on_plane(Vclip[face[v_on[0]]])
            boundary_edge = []

        elif face_type == '111': #Done
            face = np.roll(face, -v_on[0])
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
            Vclip[face[0]] = plane.orthogonal_projection_on_plane(Vclip[face[0]])
            nI += 1

        elif face_type == '120': # Done
            face = np.roll(face, -v_above[0])
            Vclip[face[1]] = plane.orthogonal_projection_on_plane(Vclip[face[1]])
            Vclip[face[2]] = plane.orthogonal_projection_on_plane(Vclip[face[2]])
            boundary_edge = [face[1], face[2]]

        elif face_type == '021': # Done
            face = np.roll(face, -v_below[0])
            Vclip[face[1]] = plane.orthogonal_projection_on_plane(Vclip[face[1]])
            Vclip[face[2]] = plane.orthogonal_projection_on_plane(Vclip[face[2]])
            boundary_edge = [face[2], face[1]]
            face = list(face)
            face.append(face[0])
            F_clip.append(face)

        elif face_type == '012': # Done
            Vclip[face[v_on[0]]] = plane.orthogonal_projection_on_plane(Vclip[face[v_on[0]]])
            boundary_edge = []
            face = list(face)
            face.append(face[0])
            F_clip.append(face)

        elif face_type == '220': # Done
            if v_above[1] == v_above[0]+1:
                face = np.roll(face, -v_above[1])
            boundary_edge = [face[1], face[2]]

        elif face_type == '121': # Done
            face = np.roll(face, -v_above[0])
            Vclip[face[1]] = plane.orthogonal_projection_on_plane(Vclip[face[1]])
            Vclip[face[3]] = plane.orthogonal_projection_on_plane(Vclip[face[3]])
            boundary_edge = [face[1], face[3]]
            F_clip.append([face[1], face[2], face[3], face[1]])

        elif face_type == '300' or face_type == '400':
            boundary_edge = []

        elif face_type == '003':
            boundary_edge = []
            face = list(face)
            face.append(face[0])
            F_clip.append(face)

        elif face_type == '004':
            boundary_edge = []
            F_clip.append(list(face))

        # Building boundary connectivity
        if len(boundary_edge) == 2:
            direct_boundary_edges[boundary_edge[0]] = boundary_edge[1]
            inv_boundary_edges[boundary_edge[1]] = boundary_edge[0]

    # It is now necessary to merge intersection vertices in order to build the intersection polygons
    intersections, newID = mm.merge_duplicates(np.asarray(intersections), return_index=True, tol=tol)

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

                if assert_closed_boundaries:
                    if len(boundary_lines) > 0:
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
    plane = Plane()

    # clip(V, F, plane)
    # sys.exit(0)

    # V, F = mm.symmetrize(V, F, plane)
    # plane.normal = np.array([0, 1, 0])
    # V, F = mm.symmetrize(V, F, plane)

    # Rotation aleatoire
    import hydrostatics as hs

    R = np.load('buggy_rotation_meshmagick.npy')
    V = np.transpose(np.dot(R, V.T))
    normal = [0, 0, 1]
    c = 0
    plane = Plane(normal, c)
    crown_face_ids, below_face_ids = split_mesh(V, F, plane)
    Vclip, F_crown = _clip(V, F[crown_face_ids], plane, return_boundaries=True)
    sys.exit(0)

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
            _clip(Vc, F[crown_face_ids], plane, return_boundaries=True)
            # mm.show(Vclip, F_crown)
        except:
            np.save('buggy_rotation_meshmagick', R)
            mm.write_VTP('all.vtp', Vc, F)
            mm.write_VTP('ring.vtp', Vc, F[crown_face_ids])
            print R
            print 'state saved'
            break
