#!/usr/bin/env python
#  -*- coding: utf-8 -*-

# TODO: voir si oblige d'importer meshmagick ici ...
import meshmagick as mm
import numpy as np
import math
import sys

# TODO: voir si on ne peut pas mettre ces fonctions dans un module dedie ?
def _rodrigues(thetax, thetay):
    """

    Parameters
    ----------
    thetax
    thetay

    Returns
    -------

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

    Parameters
    ----------
    phi
    theta

    Returns
    -------

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

    Parameters
    ----------
    thetax
    thetay
    atype

    Returns
    -------

    """
    if atype == 'fixed':
        R = _rodrigues(thetax, thetay)
    elif atype == 'cardan':
        R = _cardan(thetax, thetay)
    else:
        raise AttributeError, 'Unknown angle convention: %s' % atype

    return R


# Classes
class Plane(object):
    def __init__(self, normal=[0., 0., 1.], c=0.):
        """

        Parameters
        ----------
        normal
        c

        Returns
        -------

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

        Parameters
        ----------
        thetax
        thetay

        Returns
        -------

        """
        R = _get_rotation_matrix(thetax, thetay)
        self.normal = np.dot(R, self.normal)

        # updating self._rot
        self._rot = np.dot(R, self._rot)
        return

    def set_normal_from_angles(self, thetax, thetay):
        """

        Parameters
        ----------
        thetax
        thetay

        Returns
        -------

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

        Returns
        -------

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

        Parameters
        ----------
        thetax
        thetay
        dz

        Returns
        -------

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

    def get_point_pos_wrt_plane(self, points):
        """

        Parameters
        ----------
        points

        Returns
        -------

        """
        return np.dot(points, self.normal) - self.c

    def flip_normal(self):
        """

        Returns
        -------

        """
        self.normal *= -1
        thetax, thetay = self.get_normal_orientation_wrt_z()
        self._rot = _get_rotation_matrix(thetax, thetay)
        return

    def coord_in_plane(self, points):
        """

        Parameters
        ----------
        points

        Returns
        -------

        """
        # TODO: verifier effectivement que si on prend des points se trouvant dans le plan, leurs coordonnees dans le
        #  plan n'ont pas de composante z
        return -self.c*self.normal + np.transpose(np.dot(self._rot, points.T))

    def get_edge_intersection(self, P0, P1):
        """

        Parameters
        ----------
        P0
        P1

        Returns
        -------

        """
        P0n = np.dot(P0, self.normal)
        P1n = np.dot(P1, self.normal)
        t = (P0n-self.c) / (P0n-P1n)
        if t<0. or t>1.:
            raise RuntimeError, 'Intersection is outside the edge'
        return (1-t)*P0 + t*P1

    def orthogonal_projection_on_plane(self, point):
        """

        Parameters
        ----------
        point

        Returns
        -------

        """
        # TODO: passer en vectoriel

        return point - self.get_point_pos_wrt_plane(point) * self.normal


# class Mesh(object):
#     def __init__(self, V, F):
#         self.V = V
#         self.F = F
#         self.nv = V.shape[0]
#         self.nf = F.shape[0]

    # def build_edges(self):
    #
    #     self.F_HE = np.zeros_like(self.F, dtype=np.int)
    #     nhe = 0
    #     for face in self.F:
    #         if face[0] == face[-1]:
    #             # Triangle




def split_mesh(V, F, plane):
    """

    Parameters
    ----------
    V
    F
    plane : Plane
        Clipping plane. Vertices that are at the opposite from the normal pointing are considered under

    Returns
    -------

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
    pos = plane.get_point_pos_wrt_plane(V)

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
# Remarque: crown_face_ids sont les facettes coupees par le plan
# below_face_ids sont les facettes qui sont conservees et ne sont pas impactees par la decoupe (on pourra donc
# reutiliser les donnes d'integrales dessus !


# def __clip210(data):
#     pass
#
# clip_dict = {'210': _clip210,
#              '120': _clip120,
#              '021': _clip021,
#              '012': _clip012,
#              '111': _clip111,
#              '201': _clip201,
#              '102': _clip102,
#              '310': _clip310,
#              '220': _clip220,
#              '211': _clip211,
#              '121': _clip121,
#              '112': _clip112,
#              '022': _clip022,
#              '211': _clip211,
#              '202': _clip202,
#              '103': _clip103,
#              '301': _clip301,
#              '013': _clip013}
def clip_face(face, V, plane, pos=None):
    """

    Parameters
    ----------
    face : ndarray
    V : ndarray
    plane : Plane
    pos : ndarray, optional

    Returns
    -------

    """

    tol = 1e-4 # A placer en argument

    vertices = V[face]

    if pos is None:
        pos = plane.get_point_pos_wrt_plane(vertices)


    # Building edges
    nbv = len(face)
    edges = np.asarray([(face[i-1], face[i]) for i in xrange(nbv)])
    positions = np.asarray([(pos[i-1], pos[i]) for i in xrange(nbv)])

    v_above = np.where(positions > tol)[0]
    v_on = np.where(np.fabs(positions) <= tol)[0]
    v_below = np.where(positions < -tol)[0]

    data = [v_above, v_on, v_below]

    nb_above = v_above.sum()
    nb_on = v_on.sum()
    nb_below = v_below.sum()

    face_type = str(nb_above)+str(nb_on)+str(nb_below)



    return V, new_face, I






def clip_crown(Vinit, F_crown, plane):

    tol = 1e-4
    #
    # # arrays for triangle and quadrangle extraction in case of clipped face having 5 edges
    # triangle = [2, 3, 4, 2]
    # quadrangle = [0, 1, 2, 4]

    # TODO: voir a recuperer ces positions en entree... (seulement pour la crown !
    V = Vinit.copy()
    nv = V.shape[0]

    pos = plane.get_point_pos_wrt_plane(V) # TODO: voir a recuperer pos en entree

    F_crown_clipped = list()
    direct_boundary_edges = dict()
    inv_boundary_edges = dict()

    intersections = list()
    nI = nv

    for iface, face in enumerate(F_crown):

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

        if face_type == '210' or face_type == '310': #Done
            V[face[v_on[0]]] = plane.orthogonal_projection_on_plane(V[face[v_on[0]]])
            boundary_edge = []

        elif face_type == '120': # Done
            face = np.roll(face, -v_above[0])
            V[face[1]] = plane.orthogonal_projection_on_plane(V[face[1]])
            V[face[2]] = plane.orthogonal_projection_on_plane(V[face[2]])
            boundary_edge = [face[1], face[2]]

        elif face_type == '021': # Done
            face = np.roll(face, -v_below[0])
            V[face[1]] = plane.orthogonal_projection_on_plane(V[face[1]])
            V[face[2]] = plane.orthogonal_projection_on_plane(V[face[2]])
            boundary_edge = [face[2], face[1]]
            face = list(face)
            face.append(face[0])
            F_crown_clipped.append(face)

        elif face_type == '012' or face_type == '013': # Done
            V[face[v_on[0]]] = plane.orthogonal_projection_on_plane(V[face[v_on[0]]])
            boundary_edge = []
            face = list(face)
            face.append(face[0])
            F_crown_clipped.append(face)

        elif face_type == '111': #Done
            face = np.roll(face, -v_on[0])
            P1, P2 = V[face[[1, 2]]]
            if pos[face[1]] < 0.:
                Iright = plane.get_edge_intersection(P1, P2)
                intersections.append(Iright)
                boundary_edge = [face[0], nI]
                F_crown_clipped.append([face[0], face[1], nI, face[0]])
            else:
                Ileft = plane.get_edge_intersection(P1, P2)
                intersections.append(Ileft)
                boundary_edge = [nI, face[0]]
                F_crown_clipped.append([nI, face[2], face[0], nI])
            V[face[0]] = plane.orthogonal_projection_on_plane(V[face[0]])
            nI += 1

        elif face_type == '201':#done
            face = np.roll(face, -v_below[0])
            P0, P1, P2 = V[face]
            Ileft = plane.get_edge_intersection(P0, P2)
            Iright = plane.get_edge_intersection(P0, P1)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_crown_clipped.append([nI, face[0], nI+1, nI])
            nI += 2

        elif face_type == '102': #Done
            face = np.roll(face, -v_above[0])
            P0, P1, P2 = V[face]
            Ileft = plane.get_edge_intersection(P0, P1)
            Iright = plane.get_edge_intersection(P0, P2)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_crown_clipped.append([nI,  face[1], face[2], nI+1])
            nI += 2

        elif face_type == '220': # Done
            if v_above[1] == v_above[0]+1:
                face = np.roll(face, -v_above[1])
            boundary_edge = [face[1], face[2]]

        elif face_type == '211': # Done
            face = np.roll(face, -v_on[0])
            if pos[face[1]] < 0.:
                P1, P2 = V[face[[1, 2]]]
                Iright = plane.get_edge_intersection(P1, P2)
                intersections.append(Iright)
                boundary_edge = [face[0], nI]
                F_crown_clipped.append([face[0], face[1], nI, face[0]])
            else:
                P2, P3 = V[face[[2, 3]]]
                Ileft = plane.get_edge_intersection(P2, P3)
                intersections.append(Ileft)
                boundary_edge = [nI, face[0]]
                F_crown_clipped.append([nI, face[3], face[0], nI])
            V[face[0]] = plane.orthogonal_projection_on_plane(V[face[0]])
            nI += 1

        elif face_type == '121': # Done
            face = np.roll(face, -v_above[0])
            V[face[1]] = plane.orthogonal_projection_on_plane(V[face[1]])
            V[face[3]] = plane.orthogonal_projection_on_plane(V[face[3]])
            boundary_edge = [face[1], face[3]]
            F_crown_clipped.append([face[1], face[2], face[3], face[1]])

        elif face_type == '112': # Done
            face = np.roll(face, -v_on[0])
            if pos[face[1]] < 0.:
                P2, P3 = V[face[[2, 3]]]
                Iright = plane.get_edge_intersection(P2, P3)
                intersections.append(Iright)
                boundary_edge = [face[0], nI]
                F_crown_clipped.append([face[0], face[1], face[2], nI])
            else:
                P1, P2 = V[face[[1, 2]]]
                Ileft = plane.get_edge_intersection(P1, P2)
                intersections.append(Ileft)
                boundary_edge = [nI, face[0]]
                F_crown_clipped.append([nI, face[2], face[3], face[0]])
            V[face[0]] = plane.orthogonal_projection_on_plane(V[face[0]])
            nI += 1

        elif face_type == '202': # Done
            if v_above[1] == v_above[0]+1:
                face = np.roll(face, -v_above[1])
            P0, P1, P2, P3 = V[face]
            Ileft = plane.get_edge_intersection(P0, P1)
            Iright = plane.get_edge_intersection(P2, P3)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_crown_clipped.append([nI, face[1], face[2], nI+1])
            nI += 2

        elif face_type == '103': # Done
            face = np.roll(face, -v_above[0])
            P0, P1, P3 = V[face[[0, 1, 3]]]
            Ileft = plane.get_edge_intersection(P0, P1)
            Iright = plane.get_edge_intersection(P0, P3)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_crown_clipped.append([nI, face[1], face[3], nI+1])
            F_crown_clipped.append([face[1], face[2], face[3], face[1]])
            nI += 2

        elif face_type == '301':#Done
            face = np.roll(face, -v_below[0])
            P0, P1, P3 = V[face[[0, 1, 3]]]
            Ileft = plane.get_edge_intersection(P0, P3)
            Iright = plane.get_edge_intersection(P0, P1)
            intersections += [Ileft, Iright]
            boundary_edge = [nI, nI+1]
            F_crown_clipped.append([nI, face[0], nI+1, nI])
            nI +=2

        if len(boundary_edge) == 2:
            direct_boundary_edges[boundary_edge[0]] = boundary_edge[1]
            inv_boundary_edges[boundary_edge[1]] = boundary_edge[0]


    # It is now necessary to merge intersection vertices in order to build the intersection polygons
    intersections, newID = mm.merge_duplicates(np.asarray(intersections), return_index=True, tol=tol)

    newID = np.concatenate((np.arange(nv), newID + nv))
    F_crown_clipped = newID[F_crown_clipped]

    V = np.concatenate((V, intersections))

    # Updating dictionaries
    direct_boundary_edges = dict(zip(newID[direct_boundary_edges.keys()], newID[direct_boundary_edges.values()]))
    inv_boundary_edges = dict(zip(newID[inv_boundary_edges.keys()], newID[inv_boundary_edges.values()]))

    # mm.show(V, F_crown_clipped)
    # mm.write_VTP('ring_clipped.vtp', V, F_crown_clipped)


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
            print "%u closed polygon(s) found" % len(closed_polygons)
            print closed_polygons
            print "%u open line(s) found" % len(open_lines)
            for line in open_lines:
                print line
            if len(open_lines) > 0: # TODO: A retirer
                raise RuntimeError, 'no closed loop'

            break


    return V, F_crown_clipped



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

    # R = np.load('save.npy')
    # V = np.transpose(np.dot(R, V.T))
    # normal = [0, 0, 1]
    # c = 0
    # plane = Plane(normal, c)
    # crown_face_ids, below_face_ids = split_mesh(V, F, plane)
    # Vclip, F_crown = clip_crown(V, F[crown_face_ids], plane)
    # sys.exit(0)

    while True:
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
            clip_crown(Vc, F[crown_face_ids], plane)
            # mm.show(Vclip, F_crown)
        except:
            np.save('save', R)
            mm.write_VTP('all.vtp', Vc, F)
            mm.write_VTP('ring.vtp', Vc, F[crown_face_ids])
            print R
            print 'state saved'
            break
    # mm.show(Vclip, F_crown)
