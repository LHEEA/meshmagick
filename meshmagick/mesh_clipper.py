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


def clip_crown(V, F_crown, plane):

    F_crown_clipped = list()
    for iface, face in enumerate(F_crown):
        pos = plane.get_point_pos_wrt_plane(V[face])
        if np.all(pos <=0.):
            # Face has to be kept like that as it does not intersect the plane

            # TODO: construire un edge du polygone
            continue
        else:
            if face[0] == face[-1]:
                nb = 3
            else:
                nb = 4
            edges = np.array([(face[i-1], face[i]) for i in xrange(nb)])





if __name__ == '__main__':

    V, F = mm.load_VTP('Cylinder.vtp')
    plane = Plane()
    V, F = mm.symmetrize(V, F, plane)
    plane.normal = np.array([0, 1, 0])
    V, F = mm.symmetrize(V, F, plane)

    V[:, 2] += 0.5

    normal = [0, 0, 1]
    c = 0
    plane = Plane(normal, c)

    crown_face_ids, below_face_ids = split_mesh(V, F, plane)

    # mm.show(V, F[crown_face_ids])

    clip_crown(V, F[crown_face_ids], plane)
