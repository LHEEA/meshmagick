#!/usr/bin/env python
#  -*- coding: utf-8 -*-

# TODO: voir si oblige d'importer meshmagick ici ...
import meshmagick as mm
import numpy as np
import math

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

    def orthogonal_projection_on_plane(self, points):
        """

        Parameters
        ----------
        points

        Returns
        -------

        """
        pos = self.get_point_pos_wrt_plane(points)

        return points - np.einsum('i, j', pos, self.normal)


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


def clip_crown(Vinit, F_crown, plane):

    tol = 1e-4

    # arrays for triangle and quadrangle extraction in case of clipped face having 5 edges
    triangle = [2, 3, 4, 2]
    quadrangle = [0, 1, 2, 4]

    # TODO: voir a recuperer ces positions en entree... (seulement pour la crown !
    V = Vinit.copy()
    nv = V.shape[0]

    pos = plane.get_point_pos_wrt_plane(V)

    # Projecting vertices that are very close to the plane
    # close_vertices = np.where(np.fabs(pos) < 1e-4)[0]
    # V[close_vertices] = plane.orthogonal_projection_on_plane(V[close_vertices])

    F_crown_clipped = list()
    direct_boundary_edges = dict()
    inv_boundary_edges = dict()

    intersections = list()
    nI = nv

    for iface, face in enumerate(F_crown):
        # print face

        # TODO: reposer sur les edges generes dans une classe...
        if face[0] == face[-1]:
            face = face[:3]

        edges = np.array([(face[i-1], face[i]) for i in xrange(len(face))])

        if np.all(pos[face] <= 0.):
            # Every face vertex is below the plane. The face has to be kept intact as it does not intersect the plane.
            F_crown_clipped.append(face)

            # Extracting eventual edge already on the plane and to be added to the boundary edges
            for edge in edges:
                if np.all(pos[edge] == 0.):
                    direct_boundary_edges[edge[1]] = edge[0]
                    inv_boundary_edges[edge[0]] = edge[1]

            continue

        else:
            # The face has at least one vertex lying above the plane. We then have to clip the face.
            new_face = list()
            for edge in edges:
                pos0, pos1 = pos[edge]

                if pos0*pos1 < 0.:
                    # Edge has to be clipped
                    if pos0 > 0.:
                        if pos0 < tol:
                            ileft = edge[0]
                            new_face.append(edge[0])
                        elif -pos1 < tol:
                            ileft = edge[1]
                            continue
                        else:
                            P0, P1 = V[edge]
                            I = plane.get_edge_intersection(P0, P1)
                            intersections.append(I)
                            ileft = nI
                            new_face.append(nI)
                            nI+=1
                    else:
                        if -pos0 < tol:
                            iright = edge[0]
                            continue
                        elif pos1 < tol:
                            iright = edge[1]
                            new_face.append(edge[0])
                            new_face.append(edge[1])
                        else:
                            P0, P1 = V[edge]
                            I = plane.get_edge_intersection(P0, P1)
                            intersections.append(I)
                            new_face.append(edge[0])
                            iright = nI
                            new_face.append(nI)
                            nI += 1
                else:

                    # Edge is entirely at one side of the plane
                    if pos0 > 0.:
                        # Edge is rejected
                        continue
                    else:
                        # Edge is kept
                        new_face.append(edge[0])
                        if -pos1 < tol:
                            iright = edge[1]
                            new_face.append(edge[1])

            if len(new_face) == 0:
                # The face has been rejected
                continue

            # Boundary edge
            direct_boundary_edges[ileft] = iright
            inv_boundary_edges[iright] = ileft

            new_face_len = len(new_face)
            if new_face_len == 3:
                # Triangles
                new_face.append(new_face[0])
            elif new_face_len == 5:
                # We obtained a polygon with 5 edges --> has to be split
                new_face = np.asarray(new_face)

                newV = np.where(new_face >= nv)[0]
                print iface
                if newV[1] != newV[0]+1:
                    new_face = np.roll(new_face, 1)
                else:
                    new_face = np.roll(new_face, -newV[0])

                F_crown_clipped.append(new_face[quadrangle]) # A ajouter !
                new_face = new_face[triangle]

                pass # TODO: voir si utile

            F_crown_clipped.append(new_face)

    # It is now necessary to merge intersection vertices in order to build the intersection polygons
    intersections, newID = mm.merge_duplicates(np.asarray(intersections), return_index=True, tol=1e-4)

    newID = np.concatenate((np.arange(nv), newID + nv))
    # Updating new vertices IDs from merged vertices
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
            if len(open_lines) > 0:
                raise RuntimeError, 'no closed loop'

            break


    return V, F_crown_clipped



if __name__ == '__main__':

    V, F = mm.load_VTP('SEAREV.vtp')
    plane = Plane()

    # V, F = mm.symmetrize(V, F, plane)
    # plane.normal = np.array([0, 1, 0])
    # V, F = mm.symmetrize(V, F, plane)

    # Rotation aleatoire
    import hydrostatics as hs
    import sys

    R = np.load('save.npy')
    print R
    V = np.transpose(np.dot(R, V.T))
    normal = [0, 0, 1]
    c = 0
    plane = Plane(normal, c)
    crown_face_ids, below_face_ids = split_mesh(V, F, plane)
    Vclip, F_crown = clip_crown(V, F[crown_face_ids], plane)
    sys.exit(0)

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
            Vclip, F_crown = clip_crown(Vc, F[crown_face_ids], plane)
            # mm.show(Vclip, F_crown)
        except:
            np.save('save', R)
            mm.write_VTP('all.vtp', Vc, F)
            mm.write_VTP('ring.vtp', Vc, F[crown_face_ids])
            print R
            print 'state saved'
            break
    # mm.show(Vclip, F_crown)
