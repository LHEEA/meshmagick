#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PYTHON_ARGCOMPLTETE_OK

# Python module to manipulate 2D meshes for hydrodynamics purposes

"""
This module contains utility function to manipulate, load, save and
convert surface mesh files used by the hydrodynamics community.
Two numpy arrays are manipulated in this module : V and F.
V is the array of nodes coordinates. It is an array of shape (nv, 3) where
nv is the number of nodes in the mesh.
F is the array of cell connectivities. It is an array of shape (nf, 4) where
nf is the number of cells in the mesh. Not that it has 4 columns as we consider
flat polygonal cells up to 4 edges (quads). Triangles are obtained by repeating
the first node at the end of the cell node ID list.

IMPORTANT NOTE:
IDs of vertices are internally idexed from 0 in meshmagick. However, several mesh
file format use indexing starting at 1. This different convention might be transparent
to user and 1-indexing may not be present outside the I/O functions
"""

import os, sys
import numpy as np
import numpy as np
import math
from datetime import datetime
from numpy.core.multiarray import dtype

__year__ = datetime.now().year

__author__     = "Francois Rongere"
__copyright__  = "Copyright 2014-%u, Ecole Centrale de Nantes" % __year__
__credits__    = "Francois Rongere"
__licence__    = "CeCILL"
__version__    = "1.0"
__maintainer__ = "Francois Rongere"
__email__      = "Francois.Rongere@ec-nantes.fr"
__status__     = "Development"


real_str = r'[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[Ee][+-]?\d+)?' # Regex for floats

# TODO: for vtk usage, replace vtkUnstructuredGrid by vtkPolyData --> .vtp files...

# MESHMAGICK CLASSES
# ------------------

class Plane:
    """Class to manipulate planes.

    Planes are used for symmetrizing and clipping meshes. The internal representation
    of a plane is the equation <N.X> = c where X is a point belonging to the plane,
    N is the normal to the plane and c is the normal distance of the plane with respect
    to the coordinate system origin.

    Example:
        We can define a plane by giving a normal and a scalar.

        >>> my_plane = Plane([1, 0, 0], 0)

        will define the plane 0yz with a normal pointing towards the positive x-axis.

    """
    def __init__(self, normal=np.array([0., 0., 1.]), c=0.):
        """Plane constructor

        Parameters:
            normal: ndarray
                A numpy array of three elements defining the normal to the plane.
            c: float
                A float defining the normal distance of the plane with respect to
                the coordinate system origin.
        """
        self.normal = normal / np.linalg.norm(normal) # Ensuring a unit normal
        self.c = c

        phi, theta = self._get_angles_from_normal()
        self.Re0 = self._get_rotation_matrix(phi, theta)


    def set_position(self, z=0., phi=0., theta=0.):
        """Set the position of the plane from z and phi, theta angles (given in radian)
        instead of the normal and scalar convention."""

        self.Re0 = self._get_rotation_matrix(phi, theta)
        self.normal = self.Re0[2]
        self.c = z

        return 1


    def get_position(self):
        """Returns the position of the plane as a [z, phi, theta] numpy array
        """
        # FIXME : on ne garde plus les angles en attribut
        phi, theta = self._get_angles_from_normal()
        return np.array([self.c, phi, theta], dtype=np.float)


    def update(self, deta):
        """Update the position of the plane with respect to its current position"""

        # Computing the rotation matrix between current and final position
        Rpp_p = self._get_rotation_matrix(deta[1], deta[2])

        self.Re0 = np.dot(Rpp_p, self.Re0)
        self.normal = self.Re0[2]
        self.c = self.c * Rpp_p[2, 2] + deta[0]

        return 1


    def flip(self):
        """Flip the normal orientation"""
        self.normal = -self.normal
        #TODO : faire la mise a jour des infos d'angle !!


    def _get_angles_from_normal(self):
        """Internal method returning the orientation of the plane normal
        with respect to the coordinate system"""

        u, v, w = self.normal

        phi = math.asin(-v)
        theta = math.atan2(u, w)

        return phi, theta


    def _get_rotation_matrix(self, phi, theta):
        """Internal method returning the rotation matrix associated to the
        orientation of the plane's normal orientation angles"""

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


    def coord_in_plane(self, vertices):
        """Returns the coordinates of vertices in the local plane coordinate system"""

        # FIXME : ne fonctionne pas si on envoie un seul vertex !
        if vertices.ndim == 1: # Case where only one vertex is given
            new_vertices = np.dot(self.Re0, vertices)
            new_vertices[2] -= self.c
        else:
            new_vertices = np.array([np.dot(self.Re0, vertices[i]) for i in xrange(vertices.shape[0])], dtype=float)
            new_vertices[:, 2] -= self.c
        return new_vertices


class Mesh:
    """Class to manipulate surface meshes

    """
    def __init__(self, V, F, verbose=False):

        # Verifications on input arrays
        if V.shape[1] != 3:
            raise RuntimeError, 'V must be a nv x 3 array'
        if F.shape[1] != 4:
            raise RuntimeError, 'F must be a nf x 4 array'

        self.verbose = verbose

        # Ensuring proper type (ndarrays)
        V = np.asarray(V,   dtype=np.float)
        F = np.asarray(F, dtype=np.int)

        # Cleaning data
        V, F = clean_mesh(V, F, verbose=verbose)

        self.vertices = V

        self.nv = V.shape[0]
        self.nf = F.shape[0]

        # Partitionning faces array with triangles first and quadrangles last

        # NOTE: Here we modify the order of faces. See in the future if it is a problem.
        # It is of practical interest when it comes to generate the half-edge data structure
        # But it could be possible to order faces like that temporarily, perform fast generation
        # of the half-edge data structure and then get back to the initial faces ordering...
        triangle_mask = F[:, 0] == F[:, -1]
        self.nb_tri = sum(triangle_mask)
        self.nb_quad = self.nf - self.nb_tri

        self.faces = np.zeros((self.nf, 4), dtype=np.int)
        self.faces[:self.nb_tri] = F[triangle_mask]
        self.faces[self.nb_tri:] = F[np.logical_not(triangle_mask)]

        # Connectivity
        self.VV = None
        self.VF = None
        self.FF = None

        # properties
        self.areas   = None
        self.normals = None
        self.centers = None

        # Half-edge data structure
        self.has_HE_connectivity = False
        self.nhe    = None # Number of half-edges in the mesh
        self.F_1HE  = None # One half-edge of the face
        self.HE_F   = None # The face of the half-edge
        self.HE_nHE = None # The next half-edge of the current half-edge
        self.HE_pHE = None # The previous half-edge of the current half-edge
        self.HE_tHE = None # The twin half-edge of the current half-edge
        self.HE_iV  = None # The incident vertex of the current half-edge
        self.HE_tV  = None # The target vertex of the current half-edge
        self.V_1iHE = None # One half-edge having V as incident vertex

        self.HE_edge  = None # Stores the edge that is associated with the half-edge
        self.nb_edges = None
        self.edges  = None   # Dictionary that maps edges to halfedges couples
        self.is_boundary_edge = None
        self.is_boundaryV = None # flag telling if a vertex is a boundary vertex

        # TODO : regarder la lib pyACVD

    def generate_HE_connectivity(self):

        if self.verbose:
            print "\nGenerating half-edge data structure..."

        nbHE_tri = 3*self.nb_tri
        nbHE_quad = 4*self.nb_quad
        self.nhe = nbHE_tri + nbHE_quad

        # Filling HE_F
        self.HE_F = np.concatenate((
            np.repeat(np.arange(self.nb_tri), 3),
            np.repeat(np.arange(self.nb_tri, self.nf), 4)
        ))

        # Filling HE_iV
        self.HE_iV = np.zeros(self.nhe, dtype=np.int)
        self.HE_iV[:nbHE_tri] = np.roll(self.faces[:self.nb_tri, :3], 1, axis=1).flatten()
        self.HE_iV[nbHE_tri:] = np.roll(self.faces[self.nb_tri:, :4], 1, axis=1).flatten()

        # Filling HE_tV
        self.HE_tV = np.zeros(self.nhe, dtype=np.int)
        self.HE_tV[:nbHE_tri] = self.faces[:self.nb_tri, :3].flatten()
        self.HE_tV[nbHE_tri:] = self.faces[self.nb_tri:, :4].flatten()

        # Filling HE_nHE
        self.HE_nHE = np.zeros(self.nhe, dtype=np.int)
        self.HE_nHE[:nbHE_tri] = np.roll(np.arange(nbHE_tri).reshape((self.nb_tri, 3)), -1, axis=1).flatten()
        self.HE_nHE[nbHE_tri:] = np.roll(np.arange(nbHE_tri, self.nhe).reshape((self.nb_quad, 4)), -1, axis=1).flatten()

        # Filling HE_pHE
        self.HE_pHE = np.zeros(self.nhe, dtype=np.int)
        self.HE_pHE[:nbHE_tri] = np.roll(np.arange(nbHE_tri).reshape((self.nb_tri, 3)), 1, axis=1).flatten()
        self.HE_pHE[nbHE_tri:] = np.roll(np.arange(nbHE_tri, self.nhe).reshape((self.nb_quad, 4)), 1, axis=1).flatten()

        # Filling F_1HE
        self.F_1HE = np.concatenate((
            np.arange(nbHE_tri, step=3),
            np.arange(nbHE_tri, self.nhe, step=4)
        ))

        # Filling V_1iHE
        self.V_1iHE = np.zeros(self.nv, dtype=np.int)
        self.V_1iHE[self.faces[:self.nb_tri, :3]] = np.roll(np.arange(nbHE_tri).reshape((self.nb_tri, 3)), -1, axis=1)
        self.V_1iHE[self.faces[self.nb_tri:, :4]] = np.roll(np.arange(nbHE_tri, self.nhe).reshape((self.nb_quad, 4)),
                                                            -1, axis=1)

        # Filling HE_tHE
        edge_dst = dict([(i, dict()) for i in xrange(self.nv)])
        for iHE in xrange(self.nhe):
            iVmin, iVmax = sorted([self.HE_iV[iHE], self.HE_tV[iHE]])
            if edge_dst[iVmin].has_key(iVmax):
                edge_dst[iVmin][iVmax].append(iHE)
            else:
                edge_dst[iVmin][iVmax] = [iHE,]

        self.HE_tHE   = np.zeros(self.nhe, dtype=np.int)

        self.is_boundaryV = np.zeros(self.nv, dtype=np.bool)
        self.HE_edge = np.zeros(self.nhe, dtype=np.int)

        nb_edges = -1
        self.edges = []
        self.is_boundary_edge = []
        # nb_edges = 0
        for iV1 in edge_dst.keys():
            for iV2 in edge_dst[iV1].keys():
                helist = edge_dst[iV1][iV2]
                n_helist = len(helist)

                # Populating edges
                self.edges.append(helist[0])
                nb_edges += 1

                if n_helist == 2:
                    self.HE_tHE[helist[0]] = helist[1]
                    self.HE_tHE[helist[1]] = helist[0]
                    self.HE_edge[helist[0]] = nb_edges
                    self.HE_edge[helist[1]] = nb_edges
                    self.is_boundary_edge.append(False)
                elif n_helist >2:
                    print "WARNING: The mesh is not manifold !!!"
                else:
                    # The half-edge is on a boundary
                    self.is_boundaryV[iV1] = True
                    self.is_boundaryV[iV2] = True
                    self.HE_edge[helist[0]] = nb_edges
                    self.is_boundary_edge.append(True)

        self.edges = np.asarray(self.edges, dtype=np.int)
        self.is_boundary_edge = np.asarray(self.is_boundary_edge, dtype=np.bool)
        self.nb_edges = nb_edges+1

        # Generating half-edges for boundaries
        boundary_HE = self.edges[self.is_boundary_edge]
        nb_boundary_HE = len(boundary_HE)
        HE_F = self.HE_F[boundary_HE]
        HE_tHE = boundary_HE
        self.HE_tHE[boundary_HE] = np.arange(self.nhe, self.nhe+nb_boundary_HE)
        HE_iV  = self.HE_tV[boundary_HE]
        HE_tV  = self.HE_iV[boundary_HE]

        HE_nHE = self.HE_tHE[self.HE_pHE[self.HE_tHE[self.HE_pHE[boundary_HE]]]]
        HE_pHE = self.HE_tHE[self.HE_nHE[self.HE_tHE[self.HE_nHE[boundary_HE]]]]

        HE_edge = self.HE_edge[boundary_HE]

        self.nhe += nb_boundary_HE
        self.HE_F = np.concatenate((self.HE_F, HE_F))
        self.HE_tHE = np.concatenate((self.HE_tHE, HE_tHE))
        self.HE_iV = np.concatenate((self.HE_iV, HE_iV))
        self.HE_tV = np.concatenate((self.HE_tV, HE_tV))
        self.HE_nHE = np.concatenate((self.HE_nHE, HE_nHE))
        self.HE_pHE = np.concatenate((self.HE_pHE, HE_pHE))
        self.HE_edge = np.concatenate((self.HE_edge, HE_edge))


        self.has_HE_connectivity = True
        if self.verbose:
            print "\t -> Done !\n"

        return 1

    def detect_features(self,
                        thetaf=10, # Threshold for l-strongness in DA of half-edges
                        thetaF=65, # Threshold for u-strongness in DA of edges
                        thetaD=40, # Threshold for u-strongness in AD of vertices (sharpness)
                        thetat=20, # Threshold for u-strongness in OSTA of edges
                        thetaT=40, # Threshold for u-strongness in TA of vertices
                        thetae=25, # Threshold for e-strongness in DA of edges
                        thetak=50, # Threshold for obsurity of curves
                        k=5,       # Minimum length of curves for not being obscure
                        verbose=False,
                        debug=False):

        if verbose:
            print "\nDetecting features of the mesh"

        la = np.linalg
        epsilon = math.tan(thetaf*math.pi/180/2)**2
        pi = math.pi
        acos = math.acos

        # generating faces properties
        if self.areas is None:
            self.generate_faces_properties()

        if not self.has_HE_connectivity:
            self.generate_HE_connectivity()


        # Computing Dihedral Angle (DA) of each edge
        dihedral_angle_edge = np.zeros(self.nb_edges, dtype=np.float)
        for iedge, iHE1 in enumerate(self.edges):
            if self.is_boundary_edge[iedge]:
                dihedral_angle_edge[iedge] = pi
            else:
                iface1 = self.HE_F[iHE1]
                n1 = self.normals[iface1]
                iHE2 = self.HE_tHE[iHE1]
                iface2 = self.HE_F[iHE2]
                n2 = self.normals[iface2]
                cosn1n2 = max(-1, min(1, np.dot(n1, n2)))
                dihedral_angle_edge[iedge] = acos(cosn1n2)

        # For each face, computing the covariance matrix, needed to compute
        # the medial quadric of each vertex (giving the ridge direction)
        covF = [np.zeros((3,3), dtype=np.float) for i in xrange(self.nf)]
        for iface, face in enumerate(self.faces):
            normal = self.normals[iface]
            covF[iface] = np.outer(normal, normal)

        # Computing data for each vertex
        angle_defect_V    = np.zeros(self.nv, dtype=np.float)
        medial_quadric    = np.zeros((3,3), dtype=np.float)
        ridge_direction_V = np.zeros((self.nv, 3), dtype=np.float)
        sharp_corner      = np.zeros(self.nv, dtype=np.bool)
        ambiguous_vertex  = np.zeros(self.nv, dtype=np.bool)
        incident_HE_list  = []

        OSTA = np.zeros(self.nhe, dtype=np.float)

        for iV, vertex in enumerate(self.vertices):

            # Getting list of incident half-edges
            iHE_list = self.get_incident_HE_from_V(iV)

            # Building lists of incident half-edges for each vertex, once for all
            incident_HE_list.append(iHE_list)

            # Computing angle defect
            HE1 = self.vertices[self.HE_tV[iHE_list[-1]]] - vertex
            HE1 /= la.norm(HE1)
            for iHE in iHE_list:
                HE2 = self.vertices[self.HE_tV[iHE]] - vertex
                HE2 /= la.norm(HE2)
                # print acos(np.dot(HE1, HE2))*180/pi
                angle_defect_V[iV] -= acos(np.dot(HE1, HE2))
                HE1 = HE2

            angle_defect_V[iV] += 2*pi

            # Is the vertex u-strong in AD (Angle Defect) <=> is a sharp vertex ?
            if angle_defect_V[iV] > thetaD*pi/180:
                # print iV, ' is a sharp corner'
                sharp_corner[iV] = True

            # Computing the ridge direction
            # -----------------------------
            if self.is_boundaryV[iV]:
                # ridge direction is here undefined, we evaluate it differently
                i1, i2 = list(np.where(self.is_boundary_edge[self.HE_edge[iHE_list]])[0])
                iHE1 = iHE_list[i1]
                iHE2 = iHE_list[i2]
                HE1 = vertex - self.vertices[self.HE_tV[iHE1]]
                HE1 /= la.norm(HE1)
                HE2 = self.vertices[self.HE_tV[iHE2]] - vertex
                HE2 /= la.norm(HE2)
                ridge_dir = HE1 + HE2
                ridge_dir /= la.norm(ridge_dir)
                ridge_direction_V[iV] = ridge_dir
                # Computing the one-sided turning angle of incident half-edges
                for iHE in iHE_list:
                    HE = self.vertices[self.HE_tV[iHE]] - vertex
                    HE /= la.norm(HE)
                    cHE = max(-1, min(1, np.dot(HE, ridge_dir)))
                    OSTA[iHE] = acos(cHE)
                continue

            iface_list = self.HE_F[iHE_list]

            adjF_areas = self.areas[iface_list]

            medial_quadric[:] = 0.
            for i, iface in enumerate(iface_list):
                medial_quadric += adjF_areas[i] * covF[iface]

            eigval, eigvec = la.eigh(medial_quadric)
            (lambda3, lambda2, lambda1) = eigval

            if lambda2/lambda1 >= epsilon and lambda3/lambda2 <= 0.7:
                # Ridge direction is defined, we are not on a flat surface
                ridge_direction_V[iV] = eigvec[:, 0] / la.norm(eigvec[:, 0])

                # Computing the one-sided turning angle of each incident half-edge
                for iHE in iHE_list:
                    HE = self.vertices[self.HE_tV[iHE]] - vertex
                    HE /= la.norm(HE)
                    cHE = max(-1, min(1, np.dot(HE, ridge_direction_V[iV])))
                    OSTA[iHE] = acos(cHE)
            else:
                ambiguous_vertex[iV] = True


        # Building u-strongness and e-strongness
        sharp_edge = np.zeros(self.nb_edges, dtype=np.bool)
        e_strong_DA_edge = np.zeros(self.nb_edges, dtype=np.bool)

        sharp_edge[dihedral_angle_edge > thetaF*pi/180] = True
        e_strong_DA_edge[dihedral_angle_edge > thetae*pi/180] = True

        # Each vertex, determining l-strongness in OSTA and DA
        l_strong_OSTA_HE = np.zeros(self.nhe, dtype=np.bool)
        l_strong_DA_HE   = np.zeros(self.nhe, dtype=np.bool)
        for iV in xrange(self.nv):
            if ambiguous_vertex[iV]:
                continue

            # Getting list of incident half-edges
            iHE_list = incident_HE_list[iV]

            # Getting half-edges whose DA is greater than thetaf
            iHE_DA_gt_thetaf = iHE_list[dihedral_angle_edge[self.HE_edge[iHE_list]] > thetaf*pi/180]
            if len(iHE_DA_gt_thetaf) > 0:
                iHE = iHE_DA_gt_thetaf[np.argmin(OSTA[iHE_DA_gt_thetaf])]
                if OSTA[iHE] < thetat*pi/180:
                    l_strong_OSTA_HE[iHE] = True
                iHE = iHE_DA_gt_thetaf[np.argmax(OSTA[iHE_DA_gt_thetaf])]
                if OSTA[iHE] > pi-thetat*pi/180:
                    l_strong_OSTA_HE[iHE] = True

            iHE_OSTA_lt_pi_2 = iHE_list[OSTA[iHE_DA_gt_thetaf] < pi/2]
            if len(iHE_OSTA_lt_pi_2) > 0:
                edges = self.HE_edge[iHE_OSTA_lt_pi_2]
                iedge_maxDA = edges[np.argmax(dihedral_angle_edge[edges])]
                iHE = self.edges[iedge_maxDA]
                if self.HE_iV[iHE] == iV:
                    l_strong_DA_HE[iHE] = True
                else:
                    l_strong_DA_HE[self.HE_tHE[iHE]] = True

            iHE_OSTA_gt_pi_2 = iHE_list[OSTA[iHE_DA_gt_thetaf] > pi/2]
            if len(iHE_OSTA_gt_pi_2) > 0:
                edges = self.HE_edge[iHE_OSTA_gt_pi_2]
                iedge_maxDA = edges[np.argmax(dihedral_angle_edge[edges])]
                iHE = self.edges[iedge_maxDA]
                if self.HE_iV[iHE] == iV:
                    l_strong_DA_HE[iHE] = True
                else:
                    l_strong_DA_HE[self.HE_tHE[iHE]] = True

            # Dealing with border edges
            edges = self.HE_edge[iHE_list]
            iedges = edges[self.is_boundary_edge[edges]]
            l_strong_OSTA_HE[self.edges[iedges]] = True
            l_strong_OSTA_HE[self.HE_tHE[self.edges[iedges]]] = True
            l_strong_DA_HE[self.edges[iedges]] = True
            l_strong_DA_HE[self.HE_tHE[self.edges[iedges]]] = True

        # Determining attachment of half-edges
        attached_HE          = np.zeros(self.nhe, dtype=np.bool)
        strongly_attached_HE = np.zeros(self.nhe, dtype=np.bool)
        strongly_attached_V  = np.zeros(self.nv, dtype=np.bool)
        for iHE in xrange(self.nhe):
            iedge = self.HE_edge[iHE]
            iV = self.HE_iV[iHE] # Incident vertex of the iHE half-edge

            # List of other half-edges incident to iV
            iHE_list = incident_HE_list[iV]
            edges = self.HE_edge[iHE_list] # Associated edge list
            has_sharp_edge = np.any(sharp_edge[edges])

            if dihedral_angle_edge[iedge] > thetaf*pi/180:

                # Attached
                if l_strong_DA_HE[iHE] or \
                    l_strong_OSTA_HE[iHE] or \
                    has_sharp_edge:

                    attached_HE[iHE] = True

                # Strongly attached
                if l_strong_DA_HE[iHE] and l_strong_OSTA_HE[iHE]:
                    strongly_attached_HE[iHE] = True

                if sharp_edge[iedge] or \
                    sharp_corner[iV] or \
                    ambiguous_vertex[iV]:

                    attached_HE[iHE] = True
                    strongly_attached_HE[iHE] = True

                # Strongly attached vertex
                if strongly_attached_HE[iHE]:
                    strongly_attached_V[iV] = True

        # Quasi-strong half-edges
        quasi_strong_HE   = np.zeros(self.nhe, dtype=np.bool)
        for iHE in xrange(self.nhe):
            iV = self.HE_iV[iHE]
            tV = self.HE_tV[iHE]
            if attached_HE[iHE] and \
                strongly_attached_V[iV] and \
                strongly_attached_V[tV]:

                quasi_strong_HE[iHE] = True

        # Quasi-strong edges
        quasi_strong_edge = np.zeros(self.nb_edges, dtype=np.bool)
        for iedge in xrange(self.nb_edges):
            iHE1 = self.edges[iedge]
            iHE2 = self.HE_tHE[iHE1]
            if quasi_strong_HE[iHE1] and quasi_strong_HE[iHE2]:
                # iedge is quasi-strong edge and then, it is a candidate edge
                # Associated half-edges may be added to the ICH list
                quasi_strong_edge[iedge] = True

        # quasi-strong edges are candidate edges
        candidate_edges = list(np.where(quasi_strong_edge)[0])
        candidate_HE  = list(self.edges[candidate_edges])
        candidate_HE += list(self.HE_tHE[candidate_HE])

        # Building ICH (Incident Candidate Half-edge list for each vertex)
        ICH = dict()
        for iV in xrange(self.nv):
            iHE_list = incident_HE_list[iV]
            iedges = self.HE_edge[iHE_list]
            iHE_candidate = iHE_list[quasi_strong_edge[iedges]]
            if len(iHE_candidate) > 0:
                ICH[iV] = list(iHE_candidate)

        # Candidate vertices are ICH.keys()...

        # This is for debug of the algorithm
        if debug:
            _tmp_singleton_HE = []
            _tmp_dangling_HE = []
            _tmp_semi_joint = []
            _tmp_disjoint = []
            _tmp_multi_joint_HE = []

        # Filtration procedure starts here
        end_edge_type     = ['SG', 'DG', 'SJ', 'DJ', 'MJ']
        obscure_edge_type = ['DG', 'SJ', 'DJ']
        while 1:

            obscure_curve_found = False
            # ----------------------------
            # Collecting obscure end-edges
            # ----------------------------
            # print '\n----------------------------'
            # print 'Classifying end edges'
            # print '----------------------------'

            # Those are dangling, semi-joint and disjoint
            obscure_end_HE = []
            end_HE = []
            edge_candidate_type = dict()
            for iV in ICH.keys():
                iHE_list = ICH[iV]

                nb_iHE = len(ICH[iV])
                iedge = self.HE_edge[iHE_list]

                if nb_iHE == 1:
                    iHE = iHE_list[0]
                    edge_candidate_type[self.HE_edge[iHE]] = 'SG'
                    if debug:
                        _tmp_singleton_HE.append(iHE)
                    end_HE.append(iHE)

                    if not sharp_edge[iedge[0]] or not sharp_corner[iV]:
                        obscure_end_HE.append(iHE)
                        if debug:
                            _tmp_dangling_HE.append(iHE)
                        edge_candidate_type[self.HE_edge[iHE]] = 'DG'

                elif nb_iHE == 2:
                    iHE1 = iHE_list[0]
                    iHE2 = iHE_list[1]

                    # Computing the Turning angle of the two candidate helf-edges
                    vertex = self.vertices[iV]
                    HE1 = vertex - self.vertices[self.HE_tV[iHE1]]
                    HE1 /= la.norm(HE1)
                    HE2 = self.vertices[self.HE_tV[iHE2]] - vertex
                    HE2 /= la.norm(HE2)
                    cHE = min(1, max(-1, np.dot(HE1, HE2)))
                    turning_angle = acos(cHE)

                    nb_sharp_edges = len(sharp_edge[iedge])

                    if sharp_corner[iV] or turning_angle > thetaT*pi/180:
                        if nb_iHE != nb_sharp_edges:
                            obscure_end_HE.append(iHE1)
                            obscure_end_HE.append(iHE2)
                            edge_candidate_type[self.HE_edge[iHE1]] = 'SJ'
                            edge_candidate_type[self.HE_edge[iHE2]] = 'SJ'
                            if debug:
                                _tmp_semi_joint.append(iHE1)
                                _tmp_semi_joint.append(iHE2)
                            end_HE.append(iHE1)
                            end_HE.append(iHE2)

                elif nb_iHE > 2:
                    iedge_list = self.HE_edge[iHE_list]
                    nb_sharp_edges = len(sharp_edge[iedge_list])
                    dihedral_angles = dihedral_angle_edge[iedge_list]

                    # acute edges are not border edge that have a DA > 90°
                    acute_edge = np.logical_and(
                        np.logical_not(self.is_boundary_edge[iedge_list]),
                        dihedral_angles > pi/2.
                    )
                    nb_acute_edge = sum(acute_edge)

                    for iHE, iedge in zip(iHE_list, iedge_list):
                        is_disjoint = False

                        if not l_strong_DA_HE[iHE] and not l_strong_OSTA_HE[iHE]:
                            if not sharp_corner[iV] and not ambiguous_vertex[iV]:
                                if not sharp_edge[iedge] or not e_strong_DA_edge[iedge]:
                                    is_disjoint = True

                        elif nb_sharp_edges > 0 and not e_strong_DA_edge[iedge]:
                            is_disjoint = True

                        elif nb_acute_edge > 0 and not sharp_edge[iedge]:
                            is_disjoint = True

                        if is_disjoint:
                            obscure_end_HE.append(iHE)
                            edge_candidate_type[self.HE_edge[iHE]] = 'DJ'
                            if debug:
                                _tmp_disjoint.append(iHE)
                            end_HE.append(iHE)
                        else:
                            edge_candidate_type[self.HE_edge[iHE]] = 'MJ'
                            if debug:
                                _tmp_multi_joint_HE.append(iHE)
                            is_disjoint = False
                            end_HE.append(iHE)

            if debug:
                self._write_HE('candidate_HE.vtp', candidate_HE)
                self._write_HE('end_HE.vtp', end_HE)
                self._write_HE('multi_joint.vtp', _tmp_multi_joint_HE)
                self._write_HE('dangling.vtp', _tmp_dangling_HE)
                self._write_HE('singleton.vtp', _tmp_singleton_HE)
                self._write_HE('semi-joint.vtp', _tmp_semi_joint)
                self._write_HE('disjoint.vtp', _tmp_disjoint)

                sharp_HE = self.edges[np.where(sharp_edge)[0]]
                self._write_HE('sharp_HE.vtp', sharp_HE)
                self._write_HE('obscure_end_HE.vtp', obscure_end_HE)

            if len(obscure_end_HE) == 0:
                break

            # -----------------------------------------------------------
            # Traversing candidate curves starting from obscure end-edges
            # in order to determine obscure curves and remove them from
            # ICH list and candidate half-edges
            # -----------------------------------------------------------
            while len(obscure_end_HE) > 0:
                is_curve_closed = False
                iHE_init = obscure_end_HE.pop()

                iedge_init = self.HE_edge[iHE_init]
                iedge_init_type = edge_candidate_type[iedge_init]

                curve = [iHE_init, ]
                iHE = iHE_init

                while 1:
                    he_list = set(ICH[self.HE_tV[iHE]])

                    # Looking for the next half-edge
                    next_candidate_HE = list(he_list - set([self.HE_tHE[iHE], ]))
                    if len(next_candidate_HE) != 1:
                        iHE_end = iHE
                        iedge_end = self.HE_edge[iHE]
                        iedge_end_type = edge_candidate_type[iedge_end]
                        break

                    iHE = next_candidate_HE[0]

                    curve.append(iHE)

                    if iHE == iHE_init:
                        # closed curve
                        is_curve_closed = True
                        iHE_end = iHE
                        iedge_end = iedge_init
                        iedge_end_type = iedge_init_type
                        break

                    iedge = self.HE_edge[iHE]
                    if edge_candidate_type.has_key(iedge):
                        edge_type = edge_candidate_type[iedge]
                        if edge_type in end_edge_type:
                            # print 'Curve traversal finished with type %s'%edge_type
                            iHE_end = iHE
                            iedge_end = iedge
                            iedge_end_type = edge_candidate_type[iedge_init]
                            break

                # Determining the type of the curve
                is_curve_obscure = False

                if not is_curve_closed:
                    curve_dihedral_angles = dihedral_angle_edge[self.HE_edge[curve]]
                    nb_edge_DA_gt_thetak = sum(curve_dihedral_angles > thetak*pi/180.)

                    if ( (iedge_init_type in end_edge_type) and (iedge_end_type in end_edge_type) ) or \
                            ( iedge_init_type == 'DG' or iedge_end_type == 'DG' ):

                        # Counting the number of edges of the curve whose dihedral angle is greater than thetak
                        if nb_edge_DA_gt_thetak < k:
                            is_curve_obscure = True

                    elif (iedge_init_type in obscure_edge_type) != (iedge_end_type in obscure_edge_type):
                        he_list = incident_HE_list[self.HE_iV[iHE_init]]
                        is_any_sharp_edge = np.any(sharp_edge[self.HE_edge[he_list]])
                        he_list = incident_HE_list[self.HE_tV[iHE_end]]
                        is_any_sharp_edge *= np.any(sharp_edge[self.HE_edge[he_list]])

                        if is_any_sharp_edge and nb_edge_DA_gt_thetak == 0:
                            is_curve_obscure = True
                    else:
                        is_curve_obscure = False

                    if is_curve_obscure:
                        # The curve is obscure
                        obscure_curve_found = True

                        # Removing half-edges of the curve from the ICH connectivity table
                        for iHE in curve:
                            iV = self.HE_iV[iHE]
                            tV = self.HE_tV[iHE]
                            ICH[iV].remove(iHE)
                            ICH[tV].remove(self.HE_tHE[iHE])
                            if len(ICH[iV]) == 0:
                                ICH.pop(iV)
                            candidate_HE.remove(iHE)
                            candidate_HE.remove(self.HE_tHE[iHE])

                        # Removing also end half-edge if it is in the obscure half-edges list
                        if self.HE_tHE[iHE_end] in obscure_end_HE:
                            obscure_end_HE.remove(self.HE_tHE[iHE_end])

            if not obscure_curve_found:
                # We're done, every obscure curve has been removed
                break

        # Building feature curves from ICH list, starting from each end half_edge
        curves = []
        while len(end_HE) > 0:
            is_curve_closed = False
            iHE_init = end_HE.pop()

            curve = [iHE_init, ]
            iHE = iHE_init
            # TODO: mutualiser cette procédure de parcours avec la precedente !!
            while 1:
                he_list = set(ICH[self.HE_tV[iHE]])

                # Looking for the next half-edge
                next_candidate_HE = list(he_list - set([self.HE_tHE[iHE], ]))
                if len(next_candidate_HE) > 1:
                    curves.append(curve)
                    break

                iHE = next_candidate_HE[0]
                curve.append(iHE)

                if iHE == iHE_init:
                    # closed curve
                    is_curve_closed = True # FIXME: utile ? non utilise par la suite... --> a retirer ?
                    # Removing the last half-edge
                    curve.pop()
                    curves.append(curve)
                    break

                iedge = self.HE_edge[iHE]
                if edge_candidate_type.has_key(iedge):
                    edge_type = edge_candidate_type[iedge]
                    if edge_type in end_edge_type:
                        # Curve traversal is finished with an end edge
                        curves.append(curve)
                        break

        # ----------------------
        # Post_processing curves
        # ----------------------

        # Using a flood-fill algorithm to find out surface patches
        # --------------------------------------------------------

        # Here, we only use the remaining candidate half-edges list...
        # Precedent curve traversal is not used...

        visited_faces_mask = np.zeros(self.nf, dtype=np.bool)
        feature_half_edges_mask = np.zeros(self.nhe, dtype=np.bool)
        feature_half_edges_mask[candidate_HE] = True

        surfaces = []
        surface = []

        F_stack = []
        while 1:
            if len(F_stack) == 0: # TODO : voir si le else de while peut servir pour rerentrer dans la boucle

                # Looking for a remaining unvisited face
                unvisited_faces_ids = np.where(np.logical_not(visited_faces_mask))[0]
                if len(surface) > 0:
                    surfaces.append(surface)

                if len(unvisited_faces_ids) == 0:
                    # We're done, every faces have been visited :)
                    break
                else:
                    F_stack = [unvisited_faces_ids[0]]
                    visited_faces_mask[F_stack[0]] = True
                    surface = [F_stack[0]]
                    continue

            iface = F_stack.pop()
            face_he_list = self.get_face_half_edges(iface)
            twin_he_list = self.HE_tHE[face_he_list]

            # TODO: extraire ici les courbes qui sotn touchées dans un set...
            feature_half_edges_list = face_he_list[feature_half_edges_mask[face_he_list]]

            # Looking for half-edges that are not feature and whose twin half-edge does not
            # own to a visited face
            propagation_half_edges = twin_he_list[np.logical_and(
                np.logical_not(feature_half_edges_mask[face_he_list]),
                np.logical_not(visited_faces_mask[self.HE_F[twin_he_list]])
            )]

            propagation_faces = list(self.HE_F[propagation_half_edges])
            if len(propagation_faces) > 0:
                surface += propagation_faces
                F_stack += propagation_faces
                visited_faces_mask[propagation_faces] = True


        for surface in surfaces:
            V, F = extract_faces(self.vertices, self.faces, surface)
            polydata = _build_vtkPolyData(V, F)
            if debug:
                _tmp_viewer.add_polydata(polydata, color=Color(pick_for=polydata).get_rgb())
        if debug:
            _tmp_viewer.show()
            _tmp_viewer.finalize()


        # long_surface = []
        # for i, surface in enumerate(surfaces):
        #     long_surface += surface
        #     V, F = extract_faces(self.vertices, self.faces, long_surface)
        #     write_VTP('surf%u.vtp'%i, V, F)

        if verbose:
            print "\t-> Features detected!"

        return 1

    def _write_HE(self, filename, iHE_list, color=None):
        import vtk

        nhe = len(iHE_list)

        iV = self.HE_iV[iHE_list]
        tV = self.HE_tV[iHE_list]

        half_edges = vtk.vtkPolyData()

        points = vtk.vtkPoints()

        for iV1, iV2 in zip(iV, tV):
            points.InsertNextPoint(self.vertices[iV1])
            points.InsertNextPoint(self.vertices[iV2])

        half_edges.SetPoints(points)

        if color is not None:
            # Building color data array
            colors = vtk.vtkUnsignedCharArray()
            colors.SetNumberOfComponents(3)
            colors.SetName('color')

        lines = vtk.vtkCellArray()
        for iHE in xrange(nhe):
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, 2*iHE)
            line.GetPointIds().SetId(1, 2*iHE+1)
            lines.InsertNextCell(line)
            if color is not None:
                colors.InsertNextTupleValue(color)


        half_edges.SetLines(lines)

        if color is not None:
            half_edges.GetCellData().SetScalars(colors)


        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInput(half_edges)
        writer.Write()

        return 0

    def get_face_half_edges(self, iface):
        iHE_init = self.F_1HE[iface]
        HE_list = [iHE_init, ]
        iHE = self.HE_nHE[iHE_init]
        while iHE != iHE_init:
            HE_list.append(iHE)
            iHE = self.HE_nHE[iHE]
        return np.asarray(HE_list, dtype=np.int)

    def get_adjacent_faces_to_face(self, iface):
        HE_list = self.get_face_half_edges(iface)
        twin_HE_list = self.HE_tHE[HE_list]
        return self.HE_F[twin_HE_list]

    def get_HE_vertices(self, iHE):
        return self.HE_iV[iHE], self.HE_tV[iHE]

    def get_HE_vector(self, iHE):
        return self.vertices[self.HE_tV[iHE]] - self.vertices[self.HE_iV[iHE]]


    def generate_faces_properties(self):
        self.areas, self.normals, self.centers = get_all_faces_properties(self.vertices, self.faces)
        return 1

    def show(self):
        show(self.vertices, self.faces)
        return 1

    def shown(self):
        show(self.vertices, self.faces, normals=True)
        return 1

    def switch_verbosity(self):
        self.verbose = not self.verbose
        return 1

    def generate_connectivity(self):
        self.VV, self.VF, self.FF, self.boundaries = generate_connectivity(self.vertices, self.faces)

    def half_edge_to_mesh_arrays(self):
        if self.F_1HE is None:
            raise RuntimeError, 'No half-edge data structure'

        F = []
        for initHE in self.F_1HE:
            face = [self.HE_iV[initHE]]
            iHE = self.HE_nHE[initHE]
            while iHE != initHE:
                face.append(self.HE_iV[iHE])
                iHE = self.HE_nHE[iHE]
            if len(face) == 3:
                face.append(face[0])
            F.append(face)

        return self.vertices, np.asarray(F, dtype=np.int)

    def get_incident_HE_from_V(self, iV):
        """Returns an ordered list of half-edges that are incident to V"""

        initHE = self.V_1iHE[iV]
        helist = [initHE,]
        iHE = self.HE_tHE[self.HE_pHE[initHE]]
        while iHE != initHE:
            helist.append(iHE)
            iHE = self.HE_tHE[self.HE_pHE[iHE]]
        return np.asarray(helist, dtype=np.int)

    def apply_new_numbering_vertices(self, new_V_ids):
        pass

    def apply_new_numbering_faces(self, new_F_ids):
        pass

    def plus(self):
        """Pour implementer une operation de concatenation de maillage"""
        pass


def clip_by_plane(Vinit, Finit, plane, abs_tol=1e-3, infos=False):
    """clip_by_plane(Vinit, Finit, plane, abs_tol=1e-3, infos=False)

    Performs a mesh clipping by plane

    Parameters:
        Vinit: ndarray
            numpy array of shape (nv, 3) specifying the nodes's coodinates
            of the initial mesh to be clipped. nv is the number of nodes.
        Finit: ndarray
            numpy array of shape (nf, 4) specifying the node connectivity
            for the faces description of the initial mesh to be clipped.
            nf is the number of faces in the mesh. Every face is a line of
            Finit. It specifies 4 nodes ids. In case of a triangle face,
            the convention is to repeat the first node id as the last element.
        plane: Plane
            Plane instance that defines the clipping plane
        abs_tol: float
            tolerance under which a node is considered as belonging to the
            clipping plane
        infos: boolean
            if set to True, the function also returns a dictionary named
            clip_info which embed information on the clipping procedure
            such as the clipping polygons

    Returns:
        clipped_V: ndarray
            numpy array of shape (new_nv, 3) giving the nodes of the clipped
            mesh. new_nv is the number of vertices in the clipped mesh
        clipped_F: ndarray
            numpy array of shape (new_nf, 4) giving the connectivity of faces
            for the clipped mesh. new_nf is the number of faces in the clipped
            mesh
        clip_infos [optional]: dict
            dictionary giving the following informations:
                'FkeptOldID':
                    a numpy array of size new_nf specifying for
                    each kept face the old ID it had in the initial, non clipped
                    mesh
                 'FkeptNewID':
                    a numpy array of size nf specifying for each old face, in case
                    it has been kept in the clipped mesh, its new ID in the clipped
                    mesh
                 'FToUpdateNewID':
                    a numpy array of size nm giving the list of face ids that have
                     been modified in the clipping procedure (following IDs in the
                     clipped mesh). nm is the number of modified faces. Modified
                     faces are those faces that have been intersected by the plane
                     and that potentially have degenerated from triangles to
                     quadrangles or other modification linked to that kind of
                     intersection. It also concerns the newly created faces (mainly
                     those that appeared in qudrangle intersection by plane that
                     gives a polygon with five edges and that need to be subdivided
                     into a triangle plus a quadrangle)
                 'PolygonsNewID':
                    a list of numpy arrays specifying for each clipping polygon
                    (generated by the intersection of the mesh with the clipping
                    plane) the list of node ids of the polygon. The size of this
                    list gives the number of polygons issued from the intersection
    """

    # Working on different arrays
    V = Vinit.copy()
    F = Finit.copy()

    # To store information about extraction
    # TODO : create a class clip_infos ?
    clip_infos = {'FkeptOldID':[], 'FkeptNewID':[], 'FToUpdateNewID':[], 'PolygonsNewID':[]}

    # Necessary to deal with clipping of quadrangle that give a pentagon
    triangle = [2, 3, 4, 2]
    quadrangle = [1, 2, 4, 5]

    # Classification of vertices
    nv = V.shape[0]
    nf = F.shape[0]

    # Getting the position of each vertex with respect to the plane (projected distance)
    positions = np.dot(V, plane.normal)-plane.c

    boundary_v_mask = (np.fabs(positions) <= abs_tol) # Vertices that are already on the boundary

    # Getting the vertices we are sure to keep
    keepV = positions <= abs_tol # Vertices we know we will keep

    # If the mesh is totally at one side of the plane, no need to go further !
    nb_kept_V = np.sum(keepV)
    if nb_kept_V == 0:
        # Mesh is totally above the plane, no intersection, nothing to keep --> error
        raise RuntimeError, 'Mesh is totally above the clipping plane. No cells to keep...'

    # Getting triangles and quads masks
    triangle_mask = F[:, 0] == F[:, -1]
    nb_triangles = np.sum(triangle_mask)
    quad_mask = np.invert(triangle_mask)
    nb_quads = nf-nb_triangles

    # Getting the number of kept vertex by face
    nb_V_kept_by_face = np.zeros(nf, dtype=np.int32)
    nb_V_kept_by_face[triangle_mask] = \
        np.sum(keepV[(F[triangle_mask,:3]).flatten()].reshape((nb_triangles, 3)), axis=1)
    nb_V_kept_by_face[quad_mask] = \
        np.sum(keepV[(F[quad_mask]).flatten()].reshape((nb_quads, 4)), axis=1)

    # Getting the number of vertex below the plane by face
    nb_V_below_by_face = np.zeros(nf, dtype=np.int32)
    V_below_mask = positions < -abs_tol
    nb_V_below_by_face[triangle_mask] = np.sum(V_below_mask[(F[triangle_mask, :3]).flatten()].reshape(nb_triangles,
                                                                                                       3), axis=1)
    nb_V_below_by_face[quad_mask] = np.sum(V_below_mask[(F[quad_mask]).flatten()].reshape(nb_quads, 4), axis=1)

    # Getting the faces that are kept as every of their vertices are kept
    keepF = np.zeros(nf, dtype=bool)
    keepF[np.logical_and(triangle_mask, nb_V_kept_by_face == 3)] = True
    keepF[np.logical_and(quad_mask, nb_V_kept_by_face == 4)] = True

    clip_infos['FkeptOldID'] = np.arange(nf)[keepF] # TODO : voir si on ne peut pas mettre cette ligne dans le bloc
    # infos suivant

    if infos:
        # Getting the boundary faces
        nb_V_on_boundary_by_face = np.zeros(nf, dtype=np.int32)
        nb_V_on_boundary_by_face[triangle_mask] = \
            np.sum(boundary_v_mask[(F[triangle_mask, :3]).flatten()].reshape((nb_triangles, 3)), axis=1)
        nb_V_on_boundary_by_face[quad_mask] = \
            np.sum(boundary_v_mask[(F[quad_mask]).flatten()].reshape((nb_quads, 4)), axis=1)

        # Faces that are at the boundary but that have to be clipped, sharing an edge with the boundary
        boundary_faces_mask = np.zeros(nf, dtype=bool)
        boundary_faces_mask[triangle_mask] =  np.logical_and(nb_V_on_boundary_by_face[triangle_mask] == 2,
                                                             nb_V_kept_by_face[triangle_mask] == 3)
        boundary_faces_mask[quad_mask] =  np.logical_and(nb_V_on_boundary_by_face[quad_mask] == 2,
                                                             nb_V_kept_by_face[quad_mask] == 4)

        # Building the boundary edges that are formed by the boundary_vertices
        # boundary_faces = np.array([i for i in xrange(nf)])[boundary_faces_mask]
        boundary_faces = np.arange(nf)[boundary_faces_mask]
        boundary_edges = {}
        for face in F[boundary_faces]:
            if face[0] == face[-1]:
                face_w = face[:3]
            else:
                face_w = face
            boundary_v_face_mask = boundary_v_mask[face_w]
            for (index, is_V_on_boundary) in enumerate(boundary_v_face_mask):
                if is_V_on_boundary:
                    if boundary_v_face_mask[index-1]:
                        boundary_edges[face_w[index]] = face_w[index-1]
                    else:
                        boundary_edges[face_w[index+1]] = face_w[index]
                    break

    # All the faces are kept, the mesh is totally under the plane
    if nb_kept_V == nv:
        if infos:
            # Detecting if the mesh has intersection polygons
            if boundary_v_mask.sum() > 0:
                # Computing the boundary polygons
                initV = boundary_edges.keys()[0]
                polygons = []
                while len(boundary_edges) > 0:
                    polygon = [initV]
                    closed = False
                    iV = initV
                    while 1:
                        iVtarget = boundary_edges.pop(iV)
                        polygon.append(iVtarget)
                        iV = iVtarget
                        if iVtarget == initV:
                            polygons.append(polygon)
                            if len(boundary_edges) > 0:
                                initV = boundary_edges.keys()[0]
                            break
                clip_infos['PolygonsNewID'] = polygons
            # clip_infos['FkeptOldID'] = np.arange(nf, dtype=np.int32)
            clip_infos['FkeptNewID'] = clip_infos['FkeptOldID'].copy()

            return V, F, clip_infos
        else:
            return V, F

    clipped_mask = np.zeros(nf, dtype=bool)
    clipped_mask[triangle_mask] = np.logical_and(nb_V_kept_by_face[triangle_mask] < 3,
                                                 nb_V_below_by_face[triangle_mask] > 0)
    clipped_mask[quad_mask] = np.logical_and(nb_V_kept_by_face[quad_mask] < 4,
                                                 nb_V_below_by_face[quad_mask] > 0)

    keepF[clipped_mask] = True
    clipped_faces = np.arange(nf)[clipped_mask]

    nb_kept_F = np.sum(keepF)

    # TODO : etablir ici une connectivite des faces a couper afin d'aider a la projection des vertex sur le plan

    # Initializing the mesh clipping
    nb_new_V = 0
    newV = []
    nb_new_F = 0
    newF = []
    edges = dict() # keys are ID of vertices that are above the plane


    # Loop on the faces to clip
    for (iface, face) in enumerate(F[clipped_faces]):
        # face is a copy (not a reference) of the line of F
        clipped_face_id = clipped_faces[iface]

        if triangle_mask[clipped_face_id]:
            nb = 3
        else:
            nb = 4

        pos_lst = list(keepV[face[:nb]])
        face_lst = list(face[:nb])

        for iv in range(nb-1, -1, -1):
            # For loop on vertices
            if pos_lst[iv-1] != pos_lst[iv]: # TODO : Gerer les projections ici !!!!
                # We get an edge
                # TODO : use a switch to activate vertices projections (or doing it outside based on the clip_infos data)

                iV0 = face_lst[iv-1]
                iV1 = face_lst[iv]
                V0 = V[iV0]
                V1 = V[iV1]

                if any(boundary_v_mask[[iV0, iV1]]):
                    # Case where the true vertex is on the boundary --> do not compute any intersection
                    continue

                # Storing the edge and the vertex
                if edges.has_key(iV0):
                    if iV1 not in edges[iV0][0]:
                        # We have to compute the intersection
                        Q = get_edge_intersection_by_plane(plane, V0, V1)
                        nb_new_V += 1
                        newV.append(Q)
                        id_Q = int(nv) + nb_new_V - 1

                        edges[iV0][0].append(iV1)
                        edges[iV0][1].append(id_Q)
                    else:
                        # Intersection has already been calculated
                        id_Q = edges[iV0][1][edges[iV0][0].index(iV1)]
                else:
                    # We have to compute the intersection
                    Q = get_edge_intersection_by_plane(plane, V0, V1)
                    nb_new_V += 1
                    newV.append(Q)
                    id_Q = int(nv) + nb_new_V - 1

                    edges[iV0] = [[iV1], [id_Q]]

                # Here, we know the intersection
                if edges.has_key(iV1):
                    if iV0 not in edges[iV1][0]:
                        edges[iV1][0].append(iV0)
                        edges[iV1][1].append(id_Q)
                else:
                    edges[iV1] = [[iV0], [id_Q]]


                face_lst.insert(iv, id_Q)
                pos_lst.insert(iv, True)

        face_w = np.asarray(face_lst, dtype=np.int32)
        pos = np.asarray(pos_lst, dtype=bool)

        clipped_face = face_w[pos]

        if infos:
            # Storing the boundary edge, making the orientation so that the normals of the final closed polygon will be
            # upward
            for index, ivertex in enumerate(clipped_face):
                if ivertex >= nv:
                    if clipped_face[index-1] >= nv or boundary_v_mask[clipped_face[index-1]]:
                        boundary_edges[ivertex] = clipped_face[index-1]
                    else:
                        if index < len(clipped_face)-1:
                            boundary_edges[clipped_face[index+1]] = ivertex
                        else:
                            boundary_edges[clipped_face[0]] = ivertex
                    break

        if len(clipped_face) == 3: # We get a triangle
            clipped_face = np.append(clipped_face, clipped_face[0])

        if len(clipped_face) == 5: # A quad has degenerated in a pentagon, we have to split it in two faces
            n_roll = np.where(pos==False)[0][0]
            clipped_face = np.roll(face_w, -n_roll)

            nb_new_F += 1
            quad = clipped_face[quadrangle]
            newF.append(quad)

            clipped_face = clipped_face[triangle] # Modified face

        # Updating the initial face with the clipped face
        F[clipped_face_id] = clipped_face

    # Adding new elements to the initial mesh
    # TODO : use np.append(V, ..., axis=0) instead of np.concatenate
    if nb_new_V > 0:
        V = np.concatenate((V, np.asarray(newV, dtype=np.float)))
    if nb_new_F > 0:
        F = np.concatenate((F, np.asarray(newF, dtype=np.int32)))

    extended_nb_V = nv + nb_new_V
    extended_nb_F = nf + nb_new_F
    new_nb_V = nb_kept_V + nb_new_V
    new_nb_F = nb_kept_F + nb_new_F

    # Getting the new IDs of kept faces before concatenation
    if infos:
        newID_F = np.arange(extended_nb_F)
        newID_F[keepF] = np.arange(new_nb_F)

    # extending the masks to comply with the extended mesh
    keepV = np.concatenate((keepV, np.ones(nb_new_V, dtype=bool)))
    keepF = np.concatenate((keepF, np.ones(nb_new_F, dtype=bool)))

    # Extracting the kept mesh
    clipped_V = V[keepV]
    # clipped_F = F[keepF]

    # Upgrading connectivity array with new indexing
    newID_V = np.arange(extended_nb_V)
    newID_V[keepV] = np.arange(new_nb_V)
    clipped_F = newID_V[(F[keepF]).flatten()].reshape((new_nb_F, 4)) # modif ici

    if infos:
        # Grabing faces that have been modified or added
        modifiedF = newID_F[clipped_faces]
        if nb_new_F > 0:
            modifiedF = np.concatenate((modifiedF, np.arange(nb_kept_F, new_nb_F)))
        clip_infos['FToUpdateNewID'] = modifiedF
        clip_infos['FkeptNewID'] = newID_F[clip_infos['FkeptOldID']]

        # Computing the boundary polygons
        initV = boundary_edges.keys()[0]
        polygons = []
        while len(boundary_edges) > 0:
            polygon = [initV]
            closed = False
            iV = initV
            while 1:
                iVtarget = boundary_edges.pop(iV)
                polygon.append(iVtarget)
                iV = iVtarget
                if iVtarget == initV:
                    polygons.append(polygon)
                    if len(boundary_edges) > 0:
                        initV = boundary_edges.keys()[0]
                    break
        # Upgrading with the new connectivity
        for (index, polygon) in enumerate(polygons):
            polygons[index] = newID_V[polygon]
        clip_infos['PolygonsNewID'] = polygons

    if infos:
        return clipped_V, clipped_F, clip_infos
    else:
        return clipped_V, clipped_F


def extract_faces(V, F, idF):
    """
    extract_faces(V, F, idF)

    performs the extraction of a subset of the initial mesh giving the
    ids of faces that we want to extract.

    Parameters:
        V: ndarray
            numpy array that define the coordinates of the initial mesh's
            nodes
        F: ndarray
            numpy array that defines the initial mesh's faces by their node
            connectivity
        idF: ndarray
            numpy array that defines the ids of the faces we want to extract
            from the initial mesh defined by V, F arrays

    Returns:
        Vring: ndarray
            numpy array that defines the nodes coordinates of the extracted mesh
        Fring: ndarray
            numpy array that defines the faces of the extracted mesh by their
            node connectivity
    """

    nv = V.shape[0]
    nf = F.shape[0]

    # Determination of the vertices to keep
    Vmask = np.zeros(nv, dtype=bool)
    Vmask[F[idF].flatten()] = True
    idV = np.arange(nv)[Vmask]

    # Building up the vertex array
    Vring = V[idV]
    newID_V = np.arange(nv)
    newID_V[idV] = np.arange(len(idV))

    Fring = F[idF]
    Fring = newID_V[Fring.flatten()].reshape((len(idF), 4))

    return Vring, Fring


def get_edge_intersection_by_plane(plane, V0, V1):
    """get_edge_intersection_by_plane(plane, V0, V1)

    Computes the intersection of an edge with a plane

    Parameters:
        plane: Plane
            instance of Plane used as a clipping plane
        V0: ndarray
            numpy array of size 3 giving the coordinates of the first
            edge's end point
        V1: ndarray
            numpy array of size 3 giving the coordinates of the second
            edge's end point

    Returns:
        VI: ndarray
            numpy array giving the coordinates of the intersection point
    """
    d0 = np.dot(plane.normal, V0) - plane.c
    d1 = np.dot(plane.normal, V1) - plane.c
    t = d0 / (d0-d1)
    return V0+t*(V1-V0)


def get_face_properties(V):
    """get_face_properties(V)

    Returns the prperties of a face defined by the coordinates of its
    vertices. The properties are the area, the normal and the center of
    the face. This function only deals with one face at once. For a
    boradcast version, see get_all_faces_properties function

    Parameters:
        V: ndarray
            numpy array of shape (nv, 3) defining the coordinates of the
            face's nodes. nv may be 3 of 4 for triangle or quadrangle,
            respectively

    Returns:
        area: float
            area of the face
        normal: ndarray
            coordinates of the face's normal
        center: ndarray
            coordinates of the face's center
    """

    nv = V.shape[0]

    if nv == 3: # triangle
        normal = np.cross(V[1]-V[0], V[2]-V[0])
        area = np.linalg.norm(normal)
        normal /= area
        area /= 2.
        center = np.sum(V, axis=0) / 3.
    else: # quadrangle
        normal = np.cross(V[2]-V[0], V[3]-V[1])
        normal /= np.linalg.norm(normal)
        a1 = np.linalg.norm(np.cross(V[1]-V[0], V[2]-V[0])) / 2.
        a2 = np.linalg.norm(np.cross(V[2]-V[0], V[3]-V[1])) / 2.
        area = a1 + a2
        C1 = np.sum(V[:3], axis=0) / 3.
        C2 = (np.sum(V[2:4], axis=0) + V[0]) / 3.
        center = (a1*C1 + a2*C2) / area

    return area, normal, center


def get_all_faces_properties(V, F):
    """get_all_faces_properties(V, F)

    Computes the properties for all the faces of a mesh defined by its
    V, F arrays. It exploits numpy capabilities for vectorizing computations.
    For a version dealing with only one face, see get_face_properties function.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Returns:
        areas: ndarray
            numpy array of the faces areas
        normals: ndarray
            numpy array of the faces' normals coordinates
        centers:
            numpy arrau of the faces' centers coordinates

    """

    nf = F.shape[0]

    triangle_mask = F[:,0] == F[:,-1]
    nb_triangles = np.sum(triangle_mask)
    quads_mask = np.invert(triangle_mask)
    nb_quads = nf-nb_triangles

    areas = np.zeros(nf, dtype=np.float)
    normals = np.zeros((nf, 3), dtype=np.float)
    centers = np.zeros((nf, 3), dtype=np.float)

    # Collectively dealing with triangles
    triangles = F[triangle_mask]

    triangles_normals = np.cross(V[triangles[:,1]] - V[triangles[:,0]], V[triangles[:,2]] - V[triangles[:,0]])
    triangles_areas = np.linalg.norm(triangles_normals, axis=1)
    normals[triangle_mask] = triangles_normals / np.array(([triangles_areas,]*3)).T
    areas[triangle_mask] = triangles_areas/2.
    centers[triangle_mask] = np.sum(V[triangles[:, :3]], axis=1)/3.

    # Collectively dealing with quads
    quads = F[quads_mask]

    quads_normals = np.cross(V[quads[:,2]] - V[quads[:,0]], V[quads[:,3]] - V[quads[:,1]])
    normals[quads_mask] = quads_normals / np.array(([np.linalg.norm(quads_normals, axis=1),]*3)).T

    a1 = np.linalg.norm(np.cross(V[quads[:,1]] - V[quads[:,0]], V[quads[:,2]] - V[quads[:,0]]), axis=1)/2.
    a2 = np.linalg.norm(np.cross(V[quads[:,3]] - V[quads[:,0]], V[quads[:,2]] - V[quads[:,0]]), axis=1)/2.
    areas[quads_mask] = a1 + a2

    C1 = np.sum(V[quads[:, :3]], axis=1) / 3.
    C2 = (np.sum(V[quads[:, 2:4]], axis=1) + V[quads[:, 0]]) / 3.

    centers[quads_mask] = ( np.array(([a1,]*3)).T * C1 +
          np.array(([a2,]*3)).T * C2 ) /  np.array(([areas[quads_mask],]*3)).T

    return areas, normals, centers


_mult_surf = np.array([1/6., 1/6., 1/6., 1/12., 1/12., 1/12., 1/12., 1/12., 1/12., 1/20., 1/20., 1/20., 1/60., 1/60., 1/60.], dtype=float) # Defines the array coefficient to compute surface integrals efficiently
def _get_surface_integrals(V, F, sum=True):
    """_get_surface_integrals(V, F, sum=True)

    Internal function
    Computes all the faces' integrals that may be used in several computations such
    as inertial properties of meshes. This function is partly based on the work of
    David Eberly:
    ...

    Parameters:
        V: ndarray
            numpy array of the mesh's nodes coordinates
        F: ndarray
            numpy array of the mesh's faces nodes connectivity
        sum[optional]: bool
            if let to True, the results will be summed over all faces.
            Otherwise, no sum will be performed and individual integral
            on faces will be returned

    Return:
        sint: ndarray
            numpy array of shape (nf, 15) that contains different integrals
            over the mesh faces. If sum is False, nf is equal to the number
            of faces in the mesh. Otherwise, nf=1.

            The different integrals that sint contains are:

             sint[0]  = \int x dS
             sint[1]  = \int y dS
             sint[2]  = \int z dS
             sint[3]  = \int yz dS
             sint[4]  = \int xz dS
             sint[5]  = \int xy dS
             sint[6]  = \int x^2 dS
             sint[7]  = \int y^2 dS
             sint[8]  = \int z^2 dS
             sint[9]  = \int x^3 dS
             sint[10] = \int y^3 dS
             sint[11] = \int z^3 dS
             sint[12] = \int x^2y dS
             sint[13] = \int y^2z dS
             sint[14] = \int z^2x dS
    """

    # TODO : put reference to the Eberly's work in the docstring

    nf = F.shape[0]

    if sum:
        sint = np.zeros(15, dtype=float)
    else:
        sint = np.zeros((nf, 15), dtype=float)

    tri1 = [0, 1, 2]
    tri2 = [0, 2, 3]

    cross = np.cross

    sint_tmp = np.zeros(15)
    for (iface, face) in enumerate(F):
        sint_tmp *= 0.
        # sint_tmp = np.zeros(15) # FIXME : Essai, la version precedente serait mieux !

        if face[0] == face[-1]:
            nb = 1
        else:
            nb = 2

        vertices = V[face]

        # Loop on triangles of the face
        for itri in xrange(nb):
            if itri == 0:
                triangle = tri1
            else:
                triangle = tri2

            V0, V1, V2 = vertices[triangle]
            x0, y0, z0 = V0
            x1, y1, z1 = V1
            x2, y2, z2 = V2

            d0, d1, d2 = cross(V1-V0, V2-V0)
            e1_c_e2 = math.sqrt(d0**2 + d1**2 + d2**2)

            temp0 = V0 + V1
            f1 = temp0 + V2
            temp1 = V0*V0
            temp2 = temp1 + V1*temp0
            f2 = temp2 + V2*f1
            f3 = V0*temp1 + V1*temp2 + V2*f2
            g0 = f2 + V0*(f1+V0)
            g1 = f2 + V1*(f1+V1)
            g2 = f2 + V2*(f1+V2)

            yz = z0*(4*y0 - y1 - y2)  - y0*(z1+z2) + 3*(y1*z1 + y2*z2)
            xz = x0*(4*z0 - z1 - z2)  - z0*(x1+x2) + 3*(x1*z1 + x2*z2)
            xy = y0*(4*x0 - x1 - x2)  - x0*(y1+y2) + 3*(x1*y1 + x2*y2)

            # Update integrals
            sint_tmp[0] += d0 * f1[0] # order 1 in vol, x in surf
            sint_tmp[1] += d1 * f1[1] # order 1 in vol, y in surf
            sint_tmp[2] += d2 * f1[2] # order 1 in vol, z in surf

            sint_tmp[3] +=  d0 * yz # order yz in surf
            sint_tmp[4] +=  d1 * xz # order xz in surf
            sint_tmp[5] +=  d2 * xy # order xy in surf

            sint_tmp[6] += d0 * f2[0] # order x in vol, x**2 in surf
            sint_tmp[7] += d1 * f2[1] # order y in vol, y**2 in surf
            sint_tmp[8] += d2 * f2[2] # order z in vol, z**2 in surf
            sint_tmp[9] += d0 * f3[0] # order x**2 in vol, x**3 in surf
            sint_tmp[10] += d1 * f3[1] # order y**2 in vol, y**3 in surf
            sint_tmp[11] += d2 * f3[2] # order z**2 in vol, z**3 in surf
            sint_tmp[12] += d0 * (y0*g0[0] + y1*g1[0] + y2*g2[0]) # order xy in vol, x**2*y in surf
            sint_tmp[13] += d1 * (z0*g0[1] + z1*g1[1] + z2*g2[1]) # order yz in vol, y**2*z in surf
            sint_tmp[14] += d2 * (x0*g0[2] + x1*g1[2] + x2*g2[2]) # order zx in vol, z**2*x in surf

            if sum:
                sint += sint_tmp
            else:
                sint[iface] = sint_tmp

    if sum:
        sint *= _mult_surf
    else:
        sint = np.array([sint[j]*_mult_surf for j in xrange(nf)], dtype=float)

    return sint


def get_mass_cog(V, F, rho=1.):
    """get_mass_cog(V, F, rho=1.)

    Returns the mass and the center of gravity of a mesh

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        rho[optional]: float
            specifies the density of the material enclosed by the mesh

    Returns:

    """
    # TODO: allow to specify as mush as options as in get_inertial_properties... or remove this function !

    return get_inertial_properties(V, F, rho=rho)[:2]


_mult_vol = np.array([1., 1., 1., 1., 1., 1., 1/2., 1/2., 1/2., 1/3., 1/3., 1/3., 1/2., 1/2., 1/2.]) # Defines the array coefficient to compute volume integrals on meshes
def get_inertial_properties(V, F, rho=7500., mass=None, thickness=None, shell=False, verbose=False):
    """get_inertial_properties(V, F, rho=7500., mass=None, thickness=None, shell=False, verbose=False)

    Returns the inertial properties of a mesh. The mesh may be considred as being
    filled with homogeneous material or as being a shell.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        rho[optional]: float
            the density of the material. By default, it is the steel density
            (7500 kg/m**3)
        mass[optional]: float
            the mass of the mesh. If it is specified, it will overwrite
            the density (even if explicitely given)
        thickness[optional]: float
            specifies the thickness of the hull. Used if shell is set to True.
        shell[optional]: bool
            if set to True, the mesh will be considered as a shell
        verbose[optional]: bool
            if set to True, the function will display a report on inertial
            properties of the mesh

    Returns:
        mass: float
            The mass of the mesh (computed or given...)
        cog: ndarray
            The coordinates of the center of gravity of the mesh
        inertia_matrix: ndarray
            The inertia matrix of the mesh expressed at the center
             of gravity

    """
    # TODO : allow to specify a reduction point...
    # The default density rho is that of steel

    tol = 1e-8

    # FIXME : la gestion des options n'est pas claire ici !!! Remettre les choses a plat

    if shell: # The geometry is a shell with homogeneous thickness and density.
        areas, normals = get_all_faces_properties(V, F)[:2]
        St = areas.sum() # Total surface
        # The geometry is considered as being a shell with thickness given
        if mass is None:

            # FIXME : may not work if thickness is not given !!

            # if thickness is None:
            #     # Assuming a standard thickness of 1cm
            #     thickness = 1e-2

            sigma = thickness * rho
            mass = St*sigma

        else:# A mass has been specified
            sigma = mass/St # Surfacic density

            if thickness is None:
                # Computing the equivalent thickness
                thickness = sigma / rho
            else:
                # thickness has been given, overwriting the density of the medium accordingly
                rho_tmp = sigma / thickness
                rho = rho_tmp

        # Getting surface integrals
        sint = _get_surface_integrals(V, F, sum=False)

        normals[normals==0.] = 1. # To avoid division by zero
        # Correcting these integrals by normals
        sint[:, :3] /= normals
        sint[:, 3:6] /= normals
        sint[:, 6:9] /= normals
        sint[:, 9:12] /= normals
        sint[:, 12:15] /= normals

        nu = sint.sum(axis=0)

        cog = np.array([nu[0], nu[1], nu[2]], dtype=np.float) * sigma / mass

        inertia_matrix = np.array([
            [nu[7]+nu[8] ,   -nu[5]   ,   -nu[4]],
            [  -nu[5]    , nu[6]+nu[8],   -nu[3]],
            [  -nu[4]    ,   -nu[3]   , nu[6]+nu[7]]
        ], dtype=np.float) * sigma

        if verbose:
            print '\nPrincipal inertia parameters report:'
            print '------------------------------------\n'
            print 'Total surface         : %f m**2' % St
            print 'Thickness             : %f m' % thickness
            print 'Density               : %f kg/m**3' % rho
            print 'Surface density       : %f kg/m**2' % sigma

    else:
        # The geometry is full
        sint = _get_surface_integrals(V, F)

        # Appliying multipliers
        sint *= _mult_vol

        # We take the mean of 3 possible computations from surface integrals
        vol = (sint[0] + sint[1] + sint[2]) / 3.
        cog = np.array([sint[6]/vol, sint[7]/vol, sint[8]/vol], dtype=float)

        # Inertia matrix is expressed for the moment in O
        # xx = sint[10] + sint[11] - vol*(cog[1]**2 + cog[2]**2)
        # yy = sint[9] + sint[11] - vol*(cog[2]**2 + cog[0]**2)
        # zz = sint[9] + sint[10] - vol*(cog[0]**2 + cog[1]**2)
        # xy = -(sint[12] - vol*cog[0]*cog[1])
        # yz = -(sint[13] - vol*cog[1]*cog[2])
        # xz = -(sint[14] - vol*cog[2]*cog[0])

        # Inertia matrix expressed in cog
        xx = sint[10] + sint[11]
        yy = sint[9] + sint[11]
        zz = sint[9] + sint[10]
        xy = -sint[12]
        yz = -sint[13]
        xz = -sint[14]

        mass = rho * vol
        # The inertia matrix is expressed in
        inertia_matrix = rho * np.array(
            [
                [xx, xy, xz],
                [xy, yy, yz],
                [xz, yz, zz]
            ], dtype=np.float)

    # Cleaning
    cog[np.fabs(cog) < tol] = 0.
    inertia_matrix[np.fabs(inertia_matrix) < tol] = 0.

    if verbose:
        print 'Mass                  : %f kg' % mass
        print 'COG                   : (%f, %f, %f) m' % tuple(cog)
        print 'Inertia matrix in COG : '
        print '\t%E, %E, %E\n\t%E, %E, %E\n\t%E, %E, %E\n' % tuple(inertia_matrix.flatten())

    return mass, cog, inertia_matrix


def transport_inertia_matrix(mass, cog, Ig, point, rot=np.eye(3, dtype=np.float)):
    """transport_inertia_matrix(mass, cog, Ig, point, rot)

    Performs the transport of the inertia matrix of a mesh at an other
    reduction point

    Parameters:
        mass: float
            The mass of the mesh
        cog: ndarray
            The coordinates of the center of gravity
        Ig: ndarray
            The 3x3 inertia matrix of the mesh expressed at cog
        point: ndarray
            The coordinates of the reduction point
        rot: ndarray
            The rotation matrix defining the orientation of the
            new axis system with respect to the current

    Returns:
        Ipoint:
            The 3x3 inertia matrix of the mesh expressed at the
            new reduction point and in a frame rotated by rot with
            respect to the initial coordinate system
    """

    point_cog = cog - point
    Ipoint = rot * Ig * rot.T + \
             mass * (np.eye(3, dtype=float) * np.dot(point_cog, point_cog)
                     - np.outer(point_cog, point_cog))
    return Ipoint


def get_volume(V, F):
    """get_volume(V, F)

    Returns the volume of the mesh

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Returns:
        vol: float
            The volume of the mesh
    """
    return _get_surface_integrals(V, F)[0]


def get_COM(V, F):
    """get_COM(V, F)

    Returns the center of mass (center of gravity) of the mesh

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Returns:
        com: ndarray
            Coordinates of the center of gravity of the mesh
    """
    # FIXME: la sortie est tres etonnante. Ne doit pas faire ce que ca dit !!!

    return _get_surface_integrals(V, F)


def generate_connectivity(V, F, verbose=False):

    nv = V.shape[0]
    nf = F.shape[0]

    mesh_closed = True

    # Building connectivities

    # Establishing VV and VF connectivities
    VV = dict([(i, set()) for i in xrange(nv)])
    VF = dict([(i, set()) for i in xrange(nv)])
    for (iface, face) in enumerate(F):
        if face[0] == face[-1]:
            face_w = face[:3]
        else:
            face_w = face
        for (index, iV) in enumerate(face_w):
            VF[iV].add(iface)
            VV[face_w[index-1]].add(iV)
            VV[iV].add(face_w[index-1])

    # Connectivity FF
    boundary_edges = dict()

    FF = dict([(i, set()) for i in xrange(nf)])
    for ivertex in xrange(nv):
        S1 = VF[ivertex]
        for iadjV in VV[ivertex]:
            S2 = VF[iadjV]
            I = list(S1 & S2)
            if len(I) != 1:
                FF[I[0]].add(I[1])
                FF[I[1]].add(I[0])
            else:
                boundary_face = F[I[0]]
                [iV1, iV2] = boundary_face[np.where((boundary_face==ivertex) + (boundary_face==iadjV))[0]]
                boundary_edges[iV2] = iV1

    # Computing boundaries
    boundaries = []
    while len(boundary_edges) > 0:
        boundary = [boundary_edges.pop(boundary_edges.keys()[0])]
        while boundary_edges.has_key(boundary[-1]):
            boundary.append(boundary_edges.pop(boundary[-1]))
        boundaries.append(boundary)

    return VV, VF, FF, boundaries


def remove_unused_vertices(V, F, verbose=False):
    if verbose:
        print "* Removing unused vertices in the mesh:"
    nv = V.shape[0]
    usedV = np.zeros(nv, dtype=np.bool)
    usedV[sum(map(list, F), [])] = True
    nb_usedV = sum(usedV)
    if nb_usedV < nv:
        newID_V = np.arange(nv)
        newID_V[usedV] = np.arange(nb_usedV)
        F = newID_V[F]
        V = V[usedV]
    if verbose:
        if nb_usedV < nv:
            unusedV = np.where(np.logical_not(usedV))[0]
            vlist_str = '[' + ', '.join(str(iV) for iV in unusedV) + ']'
            print "\t -> The %u following vertices were unused in the mesh and have been removed: \n\t\t%s" \
                  % (nv -nb_usedV, vlist_str)
        else:
            print "\t -> Everything is good :)"
    return V, F


def heal_normals(V, F, verbose=False): # TODO : mettre le flag a 0 en fin d'implementation -> ???
    """heal_normals(V, F, verbose=False)

    Returns the mesh with a consistent normal orientiation. Detects if
    the mesh is closed. If it is closed, it will also makke the normals outward.
    It uses a flood fill algorithm to propagate normal orientations through the
    unstructured mesh so it first build connectivity information between faces
    that are not embedded in the {V, F} mesh representation

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        verbose[optional]:
            if set to True, displays informations on the procedure

    Returns:
        F: ndarray
            numpy array of the faces' nodes connectivities corrected
            to give a consistent orientation of normals

    """
    if verbose:
        print "* Healing normals to make them consistent and if possible outward:"
    # TODO: return the different groups of a mesh in case it is made of several unrelated groups

    nv = V.shape[0]
    nf = F.shape[0]

    # Building connectivities
    VV, VF, FF, boundaries = generate_connectivity(V, F, verbose=verbose)

    if len(boundaries) > 0:
        mesh_closed = False
    else:
        mesh_closed = True

    # Flooding the mesh to find inconsistent normals
    type_cell = np.zeros(nf, dtype=np.int32)
    type_cell[:] = 4
    triangles_mask = F[:,0] == F[:,-1]
    type_cell[triangles_mask] = 3

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
            FF[iadjF].remove(iface) # So as it won't go from iadjF to iface in the future

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

    if verbose:
        if nb_reversed > 0:
            print '\t -> %u faces have been reversed to make normals consistent across the mesh' % (nb_reversed)
        else:
            print "\t -> Everything is good :)"


    # Checking if the normals are outward
    if mesh_closed:
        zmax = np.max(V[:,2])

        areas, normals, centers = get_all_faces_properties(V, F)

        hs = (np.array([(centers[:, 2]-zmax)*areas,]*3).T * normals).sum(axis=0)

        tol = 1e-9
        if math.fabs(hs[0]) > tol or math.fabs(hs[1]) > tol:
            if verbose:
                print "WARNING: the mesh does not seem watertight althought marked as closed..."

        if hs[2] < 0:
            flipped = True
            F = flip_normals(F)
        else:
            flipped = False

        if verbose and flipped:
            print '\t -> Every normals have been reversed to be outward'


    else:
        if verbose:
            #TODO : adding the possibility to plot normals on visualization
            print "\t -> Mesh is not closed, meshmagick cannot test if the normals are outward. Please consider " \
                  "checking it visually (e.g. by using --show option of meshmagick)"

    return F

def clean_mesh(V, F, verbose=False):

    if verbose:
        print "\n-----------------"
        print "Cleaning the mesh"
        print "-----------------"
    # Ensuring a correct triangle definition
    F = reformat_triangles(F, verbose=verbose) # TODO: en faire une methode enforce_triangle_rule

    # Removing unused vertices if any
    V, F = remove_unused_vertices(V, F, verbose=verbose)

    # Removing duplicate vertices
    V, F = merge_duplicates(V, F, verbose=verbose)

    # Healing normals
    F = heal_normals(V, F, verbose=verbose)

    if verbose:
        print "-----------------"
    return V, F

def merge_duplicates(V, F, verbose=False, tol=1e-8):
    """merge_duplicates(V, F, verbose=False, tol=1e-8)

    Returns a new node array where close nodes have been merged into one node (following tol). It also returns
    the connectivity array F with the new node IDs.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        verbose[optional]: bool
            if set to True, displays information on the merge procedure
        tol[optional]: float
            the tolerance used to define nodes that are coincident and
            that have to be merged

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes where
            every node is different
        F: ndarray
            numpy array of the faces' nodes connectivities, accordingly
            to the new node list that has been merged
    """

    # TODO : Set a tolerance option in command line arguments
    if verbose:
        print "* Removing duplicate vertices:"
    nv, nbdim = V.shape

    levels = [0, nv]
    Vtmp = []
    iperm = np.arange(nv)

    for dim in range(nbdim):
        # Sorting the first dimension
        values = V[:, dim].copy()
        if dim > 0:
            values = values[iperm]
        levels_tmp = []
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            if istop-istart > 1:
                level_values = values[istart:istop]
                iperm_view = iperm[istart:istop]

                iperm_tmp = level_values.argsort()

                level_values[:] = level_values[iperm_tmp]
                iperm_view[:] = iperm_view[iperm_tmp]

                levels_tmp.append(istart)
                vref = values[istart]

                for idx in xrange(istart, istop):
                    cur_val = values[idx]
                    if np.abs(cur_val - vref) > tol:
                        levels_tmp.append(idx)
                        vref = cur_val

            else:
                levels_tmp.append(levels[ilevel])
        if len(levels_tmp) == nv:
            # No duplicate vertices
            if verbose:
                print "\t -> No duplicate vertices detected :)"
            break

        levels_tmp.append(nv)
        levels = levels_tmp

    else:
        # Building the new merged node list
        Vtmp = []
        newID = np.arange(nv)
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            Vtmp.append(V[iperm[istart]])
            newID[iperm[range(istart, istop)]] = ilevel
        V = np.array(Vtmp, dtype=float)
        # Applying renumbering to cells
        for cell in F:
            cell[:] = newID[cell]

        if verbose:
            nv_new = V.shape[0]
            print "\t -> Initial number of nodes : {:d}".format(nv)
            print "\t -> New number of nodes     : {:d}".format(nv_new)
            print "\t -> {:d} nodes have been merged".format(nv-nv_new)

    return V, F

def concatenate(V1, F1, V2, F2):
    """
    Concatenates two meshes

    """
    nv1 = V1.shape[0]
    nv2 = V2.shape[0]

    V = np.concatenate((V1, V2), axis=0)
    F = np.concatenate((F1, F2+nv1), axis=0)

    return V, F


# =======================================================================
# MESH LOADERS
# ======================================================================
# Contains here all functions to load meshes from different file formats

def load_mesh(filename, format):
    """load_mesh(filename, format)

    Driver function that loads every mesh file format known by meshmagick
    and returns the node list and the connectivity array

    Parameters:
        filename: str
            name of the meh file on disk
        format: str
            format of the mesh defined in the extension_dict dictionary

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
    """
    os.path.isfile(filename)

    if not extension_dict.has_key(format):
        raise IOError, 'Extension ".%s" is not known' % format

    loader = extension_dict[format][0]

    V, F = loader(filename)

    return V, F


def load_RAD(filename):
    """load_RAD(filename)

    Loads RADIOSS mesh files. This export file format may be chosen in ICEM
    meshing program.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: RAD files have a 1-indexing
    """

    import re

    ifile = open(filename, 'r')
    data = ifile.read()
    ifile.close()

    # node_line = r'\s*\d+(?:\s*' + real_str + '){3}'
    node_line = r'\s*\d+\s*(' + real_str + ')\s*(' + real_str + ')\s*(' + real_str + ')'
    node_section = r'((?:' + node_line + ')+)'

    elem_line = r'^\s*(?:\d+\s+){6}\d+\s*[\r\n]+'
    elem_section = r'((?:' + elem_line + '){3,})'

    pattern_node_line = re.compile(node_line, re.MULTILINE)
    pattern_node_line_group = re.compile(node_line, re.MULTILINE)
    pattern_elem_line = re.compile(elem_line, re.MULTILINE)
    pattern_node_section = re.compile(node_section, re.MULTILINE)
    pattern_elem_section = re.compile(elem_section, re.MULTILINE)

    V = []
    node_section = pattern_node_section.search(data).group(1)
    for node in pattern_node_line.finditer(node_section):
        V.append(map(float, list(node.groups())))
    V = np.asarray(V, dtype=float)

    F = []
    elem_section = pattern_elem_section.search(data).group(1)
    for elem in pattern_elem_line.findall(elem_section):
        F.append(map(int, elem.strip().split()[3:]))
    F = np.asarray(F, dtype=np.int) - 1

    return V, F


def load_HST(filename):
    """load_HST(filename)

    Loads HYDROSTAR (Bureau Veritas (c)) mesh files.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: HST files have a 1-indexing
    """
    ifile = open(filename, 'r')
    data = ifile.read()
    ifile.close()

    import re

    node_line = r'\s*\d+(?:\s+' + real_str + '){3}'
    node_section = r'((?:' + node_line + ')+)'

    elem_line = r'^\s*(?:\d+\s+){3}\d+\s*[\r\n]+'
    elem_section = r'((?:' + elem_line + ')+)'

    pattern_node_line = re.compile(node_line, re.MULTILINE)
    pattern_elem_line = re.compile(elem_line, re.MULTILINE)
    pattern_node_section = re.compile(node_section, re.MULTILINE)
    pattern_elem_section = re.compile(elem_section, re.MULTILINE)

    Vtmp = []
    nv = 0
    for node_section in pattern_node_section.findall(data):
        for node in pattern_node_line.findall(node_section):
            Vtmp.append(map(float, node.split()[1:]))
        nvtmp = len(Vtmp)
        Vtmp = np.asarray(Vtmp, dtype=np.float)
        if nv == 0:
            V = Vtmp.copy()
            nv = nvtmp
        else:
            V = np.concatenate((V, Vtmp))
            nv += nvtmp

    Ftmp = []
    nf = 0
    for elem_section in pattern_elem_section.findall(data):
        for elem in pattern_elem_line.findall(elem_section):
            Ftmp.append(map(int, elem.split()))
        nftmp = len(Ftmp)
        Ftmp = np.asarray(Ftmp, dtype=np.int)
        if nf == 0:
            F = Ftmp.copy()
            nf = nftmp
        else:
            F = np.concatenate((F, Ftmp))
            nf += nftmp

    return V, F-1


def load_DAT(filename):
    """Not implemented.
    Intended to load .DAT files used in DIODORE (PRINCIPIA (c))
    """
    raise NotImplementedError


def load_INP(filename):
    """load_INP(filename)

    Loads DIODORE (PRINCIPIA (c)) configuration file format. It parses
    the .INP file and extract meshes defined in subsequent .DAT files
    using the different informations contained in the .INP file.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: INP/DAT files use a 1-indexing
    """
    import re

    ifile = open(filename, 'r')
    text = ifile.read()
    ifile.close()

    # Retrieving frames into a dictionnary frames
    pattern_FRAME_str = r'^\s*\*FRAME,NAME=(.+)[\r\n]+(.*)'
    pattern_FRAME = re.compile(pattern_FRAME_str, re.MULTILINE)

    frames = {}
    for match in pattern_FRAME.finditer(text):
        framename = match.group(1).strip()
        framevector = re.split(r'[, ]', match.group(2).strip())
        frames[framename] = np.asarray(map(float, framevector))

    # Storing the inp layout into a list of dictionnary
    pattern_NODE_ELEMENTS = re.compile(r'^\s*\*(NODE|ELEMENT),(.*)', re.MULTILINE)
    layout = []
    meshfiles = {}
    for match in pattern_NODE_ELEMENTS.finditer(text):
        fielddict = {}
        fielddict['type'] = match.group(1)
        if fielddict['type'] == 'NODE':
            fielddict['INCREMENT'] = 'NO'
        opts = match.group(2).split(',')
        for opt in opts:
            key, pair = opt.split('=')
            fielddict[key] = pair.strip()

        # Retrieving information on meshfiles and their usage
        file = fielddict['INPUT']
        if file in meshfiles:
            meshfiles[file][fielddict['type'] + '_CALL_INP'] += 1
        else:
            meshfiles[file] = {}
            meshfiles[file]['NODE_CALL_INP'] = 0
            meshfiles[file]['ELEMENT_CALL_INP'] = 0
            meshfiles[file][fielddict['type'] + '_CALL_INP'] += 1

        layout.append(fielddict)

        # RETRIEVING DATA SECTIONS FROM MESHFILES
        # patterns for recognition of sections
    node_line = r'\s*\d+(?:\s+' + real_str + '){3}'
    node_section = r'((?:' + node_line + ')+)'
    elem_line = r'^ +\d+(?: +\d+){3,4}[\r\n]+'  # 3 -> triangle, 4 -> quadrangle
    elem_section = r'((?:' + elem_line + ')+)'
    pattern_node_line = re.compile(node_line, re.MULTILINE)
    pattern_elem_line = re.compile(elem_line, re.MULTILINE)
    pattern_node_section = re.compile(node_section, re.MULTILINE)
    pattern_elem_section = re.compile(elem_section, re.MULTILINE)

    for file in meshfiles:
        try:
            meshfile = open(os.path.join(os.path.dirname(filename), file + '.DAT'), 'r')
        except:
            raise IOError, u'File {0:s} not found'.format(file + '.DAT')
        data = meshfile.read()
        meshfile.close()

        node_section = pattern_node_section.findall(data)
        if len(node_section) > 1:
            raise IOError, """Several NODE sections into a .DAT file is not supported by meshmagick
                              as it is considered as bad practice"""
        node_array = []
        idx_array = []
        for node in pattern_node_line.findall(node_section[0]):
            node = node.split()

            node[0] = int(node[0])
            idx_array.append(node[0])
            node[1:] = map(float, node[1:])
            node_array.append(node[1:])

        meshfiles[file]['NODE_SECTION'] = node_array

        # Detecting renumberings to do
        real_idx = 0
        # renumberings = []
        id_new = - np.ones(max(idx_array) + 1, dtype=np.int)
        # FIXME: cette partie est tres buggee !!!
        for i, idx in enumerate(idx_array):
            id_new[idx] = i+1

        meshfiles[file]['ELEM_SECTIONS'] = []
        for elem_section in pattern_elem_section.findall(data):

            elem_array = []
            for elem in pattern_elem_line.findall(elem_section):
                elem = map(int, elem.split())
                # for node in elem[1:]:
                elem = id_new[elem[1:]].tolist()
                if len(elem) == 3:  # Case of a triangle, we repeat the first node at the last position
                    elem.append(elem[0])

                elem_array.append(map(int, elem))
            meshfiles[file]['ELEM_SECTIONS'].append(elem_array)
        meshfiles[file]['nb_elem_sections'] = len(meshfiles[file]['ELEM_SECTIONS'])

        meshfiles[file]['nb_elem_sections_used'] = 0

    nbNodes = 0
    nbElems = 0
    for field in layout:
        file = field['INPUT']
        if field['type'] == 'NODE':
            nodes = np.asarray(meshfiles[file]['NODE_SECTION'], dtype=np.float)
            # Translation of nodes according to frame option id any
            nodes = translate(nodes, frames[field['FRAME']])  # TODO: s'assurer que frame est une options obligatoire...

            if nbNodes == 0:
                V = nodes.copy()
                nbNodes = V.shape[0]
                increment = False
                continue

            if field['INCREMENT'] == 'NO':
                V[idx, :] = nodes.copy()
                increment = False
            else:
                V = np.concatenate((V, nodes))
                nbNodes = V.shape[0]
                increment = True
        else:  # this is an ELEMENT section
            elem_section = np.asarray(meshfiles[file]['ELEM_SECTIONS'][meshfiles[file]['nb_elem_sections_used']],
                                      dtype=np.int)

            meshfiles[file]['nb_elem_sections_used'] += 1
            if meshfiles[file]['nb_elem_sections_used'] == meshfiles[file]['nb_elem_sections']:
                meshfiles[file]['nb_elem_sections_used'] = 0

            # Updating to new id of nodes
            elems = elem_section
            if increment:
                elems += nbNodes

            if nbElems == 0:
                F = elems.copy()
                nbElems = F.shape[0]
                continue
            else:
                F = np.concatenate((F, elems))
                nbElems = F.shape[0]

    return V, F-1


def load_TEC(filename):
    """load_TEC(filename)

    Loads TECPLOT (Tecplot (c)) mesh files. It relies on the tecplot file
    reader from the VTK library.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: TEC files have a 0-indexing
    """

    from vtk import vtkTecplotReader

    reader = vtkTecplotReader()

    # Importing the mesh from the file
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()

    nv = 0
    nf = 0

    for iblock in range(data.GetNumberOfBlocks()):
        block = data.GetBlock(iblock)
        if block.GetClassName() == 'vtkStructuredGrid':
            continue
        nvblock = block.GetNumberOfPoints()
        nfblock = block.GetNumberOfCells()

        Vtmp = np.zeros((nvblock, 3), dtype=np.float)
        for k in range(nvblock):
            Vtmp[k] = np.array(block.GetPoint(k))

        if nv == 0:
            V = Vtmp
        else:
            V = np.concatenate((V, Vtmp))

        nv += nvblock

        # Facet extraction
        Ftmp = np.zeros((nfblock, 4), dtype=np.int)
        for k in range(nfblock):
            cell = block.GetCell(k)
            nv_facet = cell.GetNumberOfPoints()
            for l in range(nv_facet):
                Ftmp[k][l] = cell.GetPointId(l)
            if nv_facet == 3:
                Ftmp[k][l] = Ftmp[k][0]

        if nf == 0:
            F = Ftmp
        else:
            F = np.concatenate((F, Ftmp))

        nf += nfblock

    return V, F


def load_VTU(filename):
    """load_VTU(filename)

    Loads VTK file format in the new XML format (vtu file extension for
    unstructured meshes). It relies on the reader from the VTK library.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: VTU files have a 0-indexing
    """

    from vtk import vtkXMLUnstructuredGridReader
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    V, F = _dump_vtk(vtk_mesh)
    return V, F


def load_VTP(filename):
    """load_VTP(filename)

    Loads VTK file format in the new XML format (vtp file extension for
    polydata meshes). It relies on the reader from the VTK library.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: VTP files have a 0-indexing
    """

    from vtk import vtkXMLPolyDataReader
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    V, F = _dump_vtk(vtk_mesh)
    return V, F


def load_VTK(filename):
    """load_VTK(filename)

    Loads VTK file format in the legacy format (vtk file extension).
    It relies on the reader from the VTK library.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: VTU files have a 0-indexing
    """
    from vtk import vtkUnstructuredGridReader
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    V, F = _dump_vtk(vtk_mesh)
    return V, F


def _dump_vtk(vtk_mesh):
    """_dump_vtk(vtk_mesh)

    Internal driver function that uses the VTK library to read VTK polydata
    or vtk unstructured grid data structures

    Parameters:
        filename: str
            name of the meh file on disk
        reader: Reader
            the reader to use (new XML format ot legacy vtk format)

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
    """
    # Importing the mesh from the file
    # reader.SetFileName(filename)
    # reader.Update()
    # vtk_mesh = reader.GetOutput()

    nv = vtk_mesh.GetNumberOfPoints()
    V = np.zeros((nv, 3), dtype=np.float)
    for k in range(nv):
        V[k] = np.array(vtk_mesh.GetPoint(k))

    nf = vtk_mesh.GetNumberOfCells()
    F = np.zeros((nf, 4), dtype=np.int)
    for k in range(nf):
        cell = vtk_mesh.GetCell(k)
        nv_facet = cell.GetNumberOfPoints()
        for l in range(nv_facet):
            F[k][l] = cell.GetPointId(l)
        if nv_facet == 3:
            F[k][3] = F[k][0]

    return V, F


def load_STL(filename):
    """load_STL(filename)

    Loads STL file format. It relies on the reader from the VTK library.
    As STL file format maintains a redundant set of vertices for each faces
    of the mesh, it returns a merged list of nodes and connectivity array
    by using the merge_duplicates function.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: STL files have a 0-indexing
    """

    from vtk import vtkSTLReader

    reader = vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()

    data = reader.GetOutputDataObject(0)

    nv = data.GetNumberOfPoints()
    V = np.zeros((nv, 3), dtype=np.float)
    for k in range(nv):
        V[k] = np.array(data.GetPoint(k))
    nf = data.GetNumberOfCells()
    F = np.zeros((nf, 4), dtype=np.int)
    for k in range(nf):
        cell = data.GetCell(k)
        if cell is not None:
            for l in range(3):
                F[k][l] = cell.GetPointId(l)
                F[k][3] = F[k][0]  # always repeating the first node as stl is triangle only

    # Merging duplicates nodes
    V, F = merge_duplicates(V, F)

    return V, F


def load_NAT(filename):
    """load_NAT(filename)

    This function loads natural file format for meshes.

    Format spec :
    -------------------
    xsym    ysym
    n    m
    x1    y1    z1
    .
    .
    .
    xn    yn    zn
    i1    j1    k1    l1
    .
    .
    .
    im    jm    km    lm
    -------------------

    where :
    n : number of nodes
    m : number of cells
    x1 y1 z1 : cartesian coordinates of node 1
    i1 j1 k1 l1 : counterclock wise Ids of nodes for cell 1
    if cell 1 is a triangle, i1==l1


    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: NAT files have a 1-indexing
    """

    ifile = open(filename, 'r')
    xsym, ysym = map(int, ifile.readline().split())
    nv, nf = map(int, ifile.readline().split())

    V = []
    for i in range(nv):
        V.append(map(float, ifile.readline().split()))
    V = np.array(V, dtype=np.float)

    F = []
    for i in range(nf):
        F.append(map(int, ifile.readline().split()))
    F = np.array(F, dtype=np.int)

    ifile.close()
    return V, F-1


def load_GDF(filename):
    """load_GDF(filename)

    Loads WAMIT (Wamit INC. (c)) GDF mesh files. As GDF file format maintains
    a redundant set of vertices for each faces of the mesh, it returns a merged
    list of nodes and connectivity array by using the merge_duplicates function.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: GDF files have a 1-indexing
    """
    ifile = open(filename, 'r')

    ifile.readline()  # skip one header line
    line = ifile.readline().split()
    ulen = line[0]
    grav = line[1]

    line = ifile.readline().split()
    isx = line[0]
    isy = line[1]

    line = ifile.readline().split()
    nf = int(line[0])

    V = np.zeros((4 * nf, 3), dtype=np.float)
    F = np.zeros((nf, 4), dtype=np.int)

    iv = -1
    for icell in range(nf):

        for k in range(4):
            iv += 1
            V[iv, :] = np.array(ifile.readline().split())
            F[icell, k] = iv

    ifile.close()
    V, F = merge_duplicates(V, F, verbose=True)

    return V, F


def load_MAR(filename):
    """load_MAR(filename)

    Loads Nemoh (Ecole Centrale de Nantes) mesh files.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: MAR files have a 1-indexing
    """

    ifile = open(filename, 'r')

    ifile.readline()  # Skipping the first line of the file
    V = []
    while 1:
        line = ifile.readline()
        line = line.split()
        if line[0] == '0':
            break
        V.append(map(float, line[1:]))

    V = np.array(V, dtype=np.float)
    F = []
    while 1:
        line = ifile.readline()
        line = line.split()
        if line[0] == '0':
            break
        F.append(map(int, line))

    F = np.array(F, dtype=np.int)

    ifile.close()

    return V, F-1


# def load_STL2(filename):
#     import re
#
#     ifile = open(filename, 'r')
#     text = ifile.read()
#     ifile.close()
#
#     endl = r'(?:\n|\r|\r\n)'
#     patt_str = r"""
#             ^\s*facet\s+normal(.*)""" + endl + """
#             ^\s*outer\sloop""" + endl + """
#             ^\s*vertex\s+(.*)""" + endl + """
#             ^\s*vertex\s+(.*)""" + endl + """
#             ^\s*vertex\s+(.*)""" + endl + """
#             ^\s*endloop""" + endl + """
#             ^\s*endfacet""" + endl + """
#            """
#     pattern = re.compile(patt_str, re.MULTILINE | re.VERBOSE)
#
#     normal = []
#     V = []
#     for match in pattern.finditer(text):
#         normal.append(map(float, match.group(1).split()))
#         V.append(map(float, match.group(2).split()))
#         V.append(map(float, match.group(3).split()))
#         V.append(map(float, match.group(4).split()))
#
#     V = np.array(V, dtype=float, order='fortran')
#
#     nf = np.size(V, 0) / 3
#     F = np.zeros((nf, 4), dtype=np.int32, order='fortran')
#
#     base = np.array([1, 2, 3, 1])
#     for i in range(nf):
#         F[i, :] = base + 3 * i
#
#     return V, F


def load_MSH(filename):
    """load_MSH(filename)

    Loads GMSH .MSH mesh files. It curretly uses an external module but should rely
    on a reader written for meshmagick (or on meshpy module but this one is hard to
    install on Windows...)

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Note: MSH files have a 0-indexing
    """
    import gmsh

    myMesh = gmsh.Mesh()
    myMesh.read_msh(filename)
    V = np.array(myMesh.Verts, dtype=np.float)

    ntri = myMesh.nElmts.get(2)
    nquad = myMesh.nElmts.get(3)
    if ntri is None:
        ntri = 0
    if nquad is None:
        nquad = 0

    nel = ntri + nquad

    F = np.zeros((nel, 4), dtype=np.int)

    if ntri != 0:
        F[:ntri, :3] = myMesh.Elmts.get(2)[1]
        F[:, 3] = F[:, 0]

    if nquad != 0:
        F[ntri:, :] = myMesh.Elmts.get(3)[1]

    return V, F


#=======================================================================
#                             MESH WRITERS
#=======================================================================

def write_mesh(filename, V, F, format):
    """write_mesh(filename, format)

    Driver function that writes every mesh file format known by meshmagick

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        format: str
            format of the mesh defined in the extension_dict dictionary

    """

    if not extension_dict.has_key(format):
        raise IOError, 'Extension "%s" is not known' % format

    writer = extension_dict[format][1]

    writer(filename, V, F)

    return 1


def write_DAT(filename, V, F):
    """write_DAT(filename, V, F)

    Writes .DAT file format for the DIODORE (PRINCIPA (c)) software.
    It also displays suggestions for inclusion into the .INP configuration
    file.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    import time
    import os

    rootfilename, ext = os.path.splitext(filename)
    filename = rootfilename+ext.upper()
    ofile = open(filename, 'w')

    ofile.write('$\n$ Data for DIODORE input file : {0}\n'.format(rootfilename.upper()))
    ofile.write('$ GENERATED BY MESHMAGICK ON {0}\n$\n'.format(time.strftime('%c')))

    ofile.write('$ NODE\n')
    vertex_block = \
        ''.join(
            (
                '\n'.join(
                    ''.join(
                        (
                            '{:8d}'.format(idx+1),
                            ''.join('{:13.5E}'.format(elt) for elt in node)
                        )
                    ) for (idx, node) in enumerate(V)
                ),

                '\n*RETURN\n'
            )
        )
    ofile.write(vertex_block)

    quad_block = '$\n$ ELEMENT,TYPE=Q4C000,ELSTRUCTURE={0}'.format(rootfilename.upper())
    tri_block  = '$\n$ ELEMENT,TYPE=T3C000,ELSTRUCTURE={0}'.format(rootfilename.upper())
    nq = 0
    nt = 0
    for (idx, cell) in enumerate(F+1):
        if cell[0] != cell[-1]:
            # quadrangle
            nq += 1
            quad_block = ''.join(
                (quad_block,
                 '\n',
                 '{:8d}'.format(idx+1),
                 ''.join('{:8d}'.format(node_id) for node_id in cell)
                )
            )

        else:
            # Triangle
            nt += 1
            tri_block = ''.join(
                (tri_block,
                '\n',
                '{:8d}'.format(idx+1),
                ''.join('{:8d}'.format(node_id) for node_id in cell[:3])
                )
            )

    print '-------------------------------------------------'
    print 'Suggestion for .inp DIODORE input file :'
    print ''
    print '*NODE,INPUT={0},FRAME=???'.format(rootfilename)

    if nq > 0:
        quad_block = ''.join((quad_block, '\n*RETURN\n'))
        ofile.write(quad_block)
        print '*ELEMENT,TYPE=Q4C000,ELSTRUCTURE={0},INPUT={0}'.format(rootfilename)
    if nt > 0:
        tri_block = ''.join((tri_block, '\n*RETURN\n'))
        ofile.write(tri_block)
        print '*ELEMENT,TYPE=T3C000,ELSTRUCTURE={0},INPUT={0}'.format(rootfilename)

    print ''
    print '-------------------------------------------------'
    ofile.close()

    return 1


def write_HST(filename, V, F):
    """write_HST(filename, V, F)

    Writes .HST file format for the HYDROSTAR (Bureau Veritas (c)) software.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    ofile = open(filename, 'w')

    ofile.write(''.join((
        'PROJECT:\n',
        'USERS:   meshmagick\n\n'
        'NBODY   1\n'
        'RHO   1025.0\n'
        'GRAVITY   9.81\n\n'
    )))

    coordinates_block = ''.join((  # block
            'COORDINATES\n',
            '\n'.join(  # line
                ''.join(
                    (
                        '{:10d}'.format(idx+1),  # index
                        ''.join('{:16.6E}'.format(elt) for elt in node)  # node coordinates
                    )
                ) for (idx, node) in enumerate(V)
            ),
            '\nENDCOORDINATES\n\n'
    ))

    ofile.write(coordinates_block)

    cells_coordinates = ''.join((  # block
        'PANEL TYPE 0\n',
        '\n'.join(  # line
            ''.join(
                '{:10d}'.format(node_idx) for node_idx in cell
            ) for cell in F+1
        ),
        '\nENDPANEL\n\n'
    ))

    ofile.write(cells_coordinates)

    ofile.write('ENDFILE\n')

    ofile.close()

    print u'File {0:s} written'.format(filename)


def write_TEC(filename, V, F):
    """write_TEC(filename, V, F)

    Writes .TEC file format for the TECPLOT (Tecplot (c)) visualisation
    software. It relies on the VTK library for its writer.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    ofile = open(filename, 'w')

    nv = V.shape[0]
    nf = F.shape[0]

    ofile.write('TITLE = \" THIS FILE WAS GENERATED BY MESHMAGICK - FICHIER : {} \" \n'.format(filename))

    ofile.write('VARIABLES = \"X\",\"Y\",\"Z\" \n')
    ofile.write('ZONE T=\"MESH\" \n')
    ofile.write('N={nv:10d} ,E={nf:10d} ,F=FEPOINT, ET=QUADRILATERAL\n'.format(nv=nv, nf=nf))

    node_block = '\n'.join( # block
        ''.join(
            ''.join('{:16.6E}'.format(elt) for elt in node)
        ) for node in V
    ) + '\n'
    ofile.write(node_block)

    cells_block = '\n'.join(  # block
        ''.join(
            ''.join('{:10d}'.format(node_id) for node_id in cell)
        ) for cell in F+1
    ) + '\n'
    ofile.write(cells_block)

    ofile.close()

    return 1


def write_VTU(filename, V, F):
    """write_VTU(filename, V, F)

    Writes .vtu file format for the paraview (Kitware (c)) visualisation
    software. It relies on the VTK library for its writer. VTU files use
    the last XML file format of the VTK library.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    from vtk import vtkXMLUnstructuredGridWriter
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetDataModeToAscii()
    writer.SetFileName(filename)

    unstructured_grid = _build_vtkUnstructuredGrid(V, F)
    writer.SetInput(unstructured_grid)
    writer.Write()

    return 1


def write_VTP(filename, V, F):
    """write_VTP(filename, V, F)

    Writes .vtp file format for the paraview (Kitware (c)) visualisation
    software. It relies on the VTK library for its writer. VTP files use
    the last XML file format of the VTK library and correspond to polydata.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    from vtk import vtkXMLPolyDataWriter
    writer = vtkXMLPolyDataWriter()
    writer.SetDataModeToAscii()
    writer.SetFileName(filename)

    polydata = _build_vtkPolyData(V, F)
    writer.SetInput(polydata)
    writer.Write()

    return 1


def write_VTK(filename, V, F):
    """write_VTK(filename, V, F)

    Writes .vtk file format for the paraview (Kitware (c)) visualisation
    software. It relies on the VTK library for its writer. VTK files use
    the legagy ASCII file format of the VTK library.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    from vtk import vtkUnstructuredGridWriter
    writer = vtkUnstructuredGridWriter()
    writer.SetFileName(filename)

    unstructured_grid = _build_vtkUnstructuredGrid(V, F)
    writer.SetInput(unstructured_grid)
    writer.Write()

    return 1


def _write_paraview(filename, V, F, writer):
    """_write_paraview(filename, V, F)

    Internal driver function that writes vtk files to be visualised into
    the Paraview software. It relies on the VTK library.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        writer: Writer
            The writer to be used

    """
    writer.SetFileName(filename)
    vtk_mesh = _build_vtkUnstructuredGrid(V, F)
    writer.SetInput(vtk_mesh)
    writer.Write()

    return 1


def _build_vtkUnstructuredGrid(V, F):
    """_build_vtk_mesh_obj(V, F)

    Internal function that builds a VTK object for manipulation by the VTK library.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Returns: vtkObject
        the vtk object instance
    """
    import vtk

    nv = max(np.shape(V))
    nf = max(np.shape(F))

    vtk_mesh = vtk.vtkUnstructuredGrid()
    vtk_mesh.Allocate(nf, nf)

    # Building the vtkPoints data structure
    vtk_points = vtk.vtkPoints()
    vtk_points.SetNumberOfPoints(nv)
    for idx, vertex in enumerate(V):
        vtk_points.SetPoint(idx, vertex)

    vtk_mesh.SetPoints(vtk_points)  # Storing the points into vtk_mesh

    # Building the vtkCell data structure
    for cell in F:
        if cell[-1] in cell[:-1]:
            vtk_cell = vtk.vtkTriangle()
            nc = 3
        else:
            # #print 'quadrangle'
            vtk_cell = vtk.vtkQuad()
            nc = 4

        for k in range(nc):
            vtk_cell.GetPointIds().SetId(k, cell[k])

        vtk_mesh.InsertNextCell(vtk_cell.GetCellType(), vtk_cell.GetPointIds())
    return vtk_mesh

def _build_vtkPolyData(V, F):
    import vtk

    # Create a vtkPoints object and store the points in it
    points = vtk.vtkPoints()
    for point in V:
        points.InsertNextPoint(point)

    # Create a vtkCellArray to store faces
    faces = vtk.vtkCellArray()
    for face_ids in F:
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

    polyDataMesh = vtk.vtkPolyData()
    polyDataMesh.SetPoints(points)
    polyDataMesh.SetPolys(faces)

    return polyDataMesh

def write_NAT(filename, V, F):
    """write_NAT(filename, V, F)

    Writes .nat file format as defined into the load_NAT function.

    See:
        load_NAT

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    ofile = open(filename, 'w')

    nv = max(np.shape(V))
    nf = max(np.shape(F))

    ofile.write('%6u%6u\n' % (0, 0))  # lire les symmetries dans args...
    ofile.write('%6u%6u\n' % (nv, nf))
    for vertex in V:
        ofile.write('%15.6E%15.6E%15.6E\n' % (vertex[0], vertex[1], vertex[2]))
    for cell in F+1:
        ofile.write('%10u%10u%10u%10u\n' % (cell[0], cell[1], cell[2], cell[3]))

    ofile.close()

    return 1


def write_GDF(filename, V, F):
    """write_GDF(filename, V, F)

    Writes .gdf file format for the WAMIT (Wamit INC. (c)) BEM software.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    nf = max(np.shape(F))

    ofile = open(filename, 'w')

    ofile.write('GDF file generated by meshmagick\n')

    ofile.write('%16.6f%16.6f\n' % (100.0, 9.81))
    ofile.write('%12u%12u\n' % (0, 1))  # TODO : mettre les symetries en argument
    ofile.write('%12u\n' % nf)

    for cell in F:
        for k in range(4):
            Vcur = V[cell[k], :]
            ofile.write('%16.6E%16.6E%16.6E\n' % (Vcur[0], Vcur[1], Vcur[2]))

    ofile.close()

    return 1


def write_MAR(filename, V, F):
    """write_MAR(filename, V, F)

    Writes mesh files to be used with Nemoh BEM software (Ecole Centrale de Nantes)

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    # TODO: detecter la symetrie de plan Oxz

    ofile = open(filename, 'w')

    ofile.write('{0:6d}{1:6d}\n'.format(2, 0))  # TODO : mettre les symetries en argument

    nv = V.shape[0]
    for (idx, vertex) in enumerate(V):
        ofile.write('{0:6d}{1:16.6f}{2:16.6f}{3:16.6f}\n'.format(idx+1, vertex[0], vertex[1], vertex[2]))

    ofile.write('{0:6d}{1:6d}{2:6d}{3:6d}{4:6d}\n'.format(0, 0, 0, 0, 0))

    cell_block = '\n'.join(
        ''.join(u'{0:10d}'.format(elt) for elt in cell)
        for cell in F+1
    ) + '\n'
    ofile.write(cell_block)
    ofile.write('%6u%6u%6u%6u\n' % (0, 0, 0, 0))

    ofile.close()

    print 'WARNING: if you described only one part of the mesh using symmetry for Nemoh, you may manually modify the ' \
          'file header accordingly'

    return 1

def write_RAD(filename, V, F):
    raise NotImplementedError

def write_STL(filename, V, F):
    """write_STL(filename, V, F)

    Writes .stl file format. It relies on the VTK library for its writer.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    # TODO : replace this implementation by using the vtk functionalities

    # Triangulating quads
    F = triangulate_quadrangles(F, verbose=False)

    ofile = open(filename, 'w')

    ofile.write('solid meshmagick\n')

    for facet in F:
        if facet[0] != facet[3]:
            raise RuntimeError, """Only full triangle meshes are accepted in STL files.
              Please consider using the --triangulate-quadrangles option (-tq) to
              perform a prior triangulation of the mesh"""

        # Computing normal
        v0 = V[facet[0], :]
        v1 = V[facet[1], :]
        v2 = V[facet[2], :]

        n = np.cross(v1 - v0, v2 - v0)
        n /= np.linalg.norm(n)

        block_facet = ''.join(['  facet normal ', ''.join('%15.6e' % ni for ni in n) + '\n',
                               '    outer loop\n',
                               '      vertex', ''.join('%15.6e' % Vi for Vi in v0) + '\n',
                               '      vertex', ''.join('%15.6e' % Vi for Vi in v1) + '\n',
                               '      vertex', ''.join('%15.6e' % Vi for Vi in v2) + '\n',
                               '    endloop\n',
                               '  endfacet\n'])
        ofile.write(block_facet)
    ofile.write('endsolid meshmagick\n')
    ofile.close()

    return 1

def write_INP(filename, V, F):
    raise NotImplementedError


def write_MSH(filename, V, F):
    raise NotImplementedError


#=======================================================================
#                         MESH MANIPULATION HELPERS
#=======================================================================
def mesh_quality(V, F):
    """mesh_quality(V, F)

    Generates the mesh quality informations (aspect ratio...) using the verdict
    library that is embeded into the VTK library. It displays a quality report.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    # This function is reproduced from
    # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/Verdict/Testing/Python/MeshQuality.py
    import vtk
    import math

    vtk_mesh = _build_vtkUnstructuredGrid(V, F)
    quality = vtk.vtkMeshQuality()
    quality.SetInput(vtk_mesh)

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

    if vtk_mesh.GetNumberOfCells() > 0:
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
    return vtk_mesh, res


def get_info(V, F):
    """get_info(V, F)

    Displays mesh informations such as the number of nodes, faces, triangles,
    quadrangles... as well as a quality report.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    nv = np.size(V, 0)
    nf = np.size(F, 0)
    print ''
    print 'o--------------------------------------------------o'
    print '|               MESH CHARACTERISTICS               |'  #28
    print '|--------------------------------------------------|'
    print '| Number of nodes  :     %15u           |' % nv
    print '|--------------------------------------------------|'
    print '| Number of facets :     %15u           |' % nf
    print '|--------------------------------------------------|'  #51
    print '|      |          Min        |          Max        |'
    print '|------|---------------------|---------------------|'
    print '|   X  |%21E|%21E|' % (V[:, 0].min(), V[:, 0].max())
    print '|------|---------------------|---------------------|'
    print '|   Y  |%21E|%21E|' % (V[:, 1].min(), V[:, 1].max())
    print '|------|---------------------|---------------------|'
    print '|   Z  |%21E|%21E|' % (V[:, 2].min(), V[:, 2].max())
    print 'o--------------------------------------------------o'
    print ''

    _, res = mesh_quality(V, F)
    print res


def translate(V, P):
    """translate(V, P)

    Applies a 3D translation on a mesh.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        P: ndarray
            coordinates of the translation vector

    Returns:
        V: ndarray
            numpy array of the coordinates of the translated mesh's nodes
    """

    if not isinstance(P, np.ndarray):
        P = np.asarray(P, dtype=float)

    try:
        for i in range(np.size(V, 0)):
            V[i, :] += P
    except:
        raise RuntimeError, 'second argument must be a 3D list or numpy array for the translation'

    return V


def translate_1D(V, t, ddl):
    """translate_1D(V, t, ddl)

    Applies a 1D translation on a mesh, following a specified cartesian
    direction

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        t: float
            magnitude of the translation
        ddl: 'str'
            direction to chosen among, 'x', 'y' and 'z'

    Returns:
        V: ndarray
            numpy array of the coordinates of the translated mesh's nodes

    """
    if ddl == 'x':
        j = 0
    elif ddl == 'y':
        j = 1
    elif ddl == 'z':
        j = 2
    else:
        raise IOError, "ddl should be chosen among ('x', 'y', 'z')"
    V[:, j] += t
    return V


def rotate(V, rot):
    """rotate(V, rot)

    Applies a 3D rotation on a mesh.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        rot: ndarray
            3x3 numpy array specifying the rotation matrix

    Returns:
        V: ndarray
            numpy array of the coordinates of the rotated mesh's nodes
    """
    from math import cos, sin

    if not isinstance(rot, np.ndarray):
        rot = np.asarray(rot, dtype=float)

    R = np.zeros((3, 3), dtype=float, order='f')

    phi = rot[0]
    theta = rot[1]
    psi = rot[2]

    cphi = cos(phi)
    sphi = sin(phi)
    ctheta = cos(theta)
    stheta = sin(theta)
    cpsi = cos(psi)
    spsi = sin(psi)

    R[0, 0] = cpsi * ctheta
    R[0, 1] = -spsi * cphi + cpsi * stheta * sphi
    R[0, 2] = spsi * sphi + cpsi * cphi * stheta
    R[1, 0] = spsi * ctheta
    R[1, 1] = cpsi * cphi + sphi * stheta * spsi
    R[1, 2] = -cpsi * sphi + spsi * cphi * stheta
    R[2, 0] = -stheta
    R[2, 1] = ctheta * sphi
    R[2, 2] = ctheta * cphi

    return np.transpose(np.dot(R, V.T))


def rotate_1D(V, angle, ddl):
    """rotate_1D(V, angle, ddl)

    Applies a rotation on a mesh around a cartesian axis

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        angle: float
            specifies the angle of rotation
        ddl: str
            axis of rotation to be chosen among 'x', 'y' and 'z'

    Returns:
        V: ndarray
            numpy array of the coordinates of the rotated mesh's nodes
    """
    if ddl == 'x':
        j = 0
    elif ddl == 'y':
        j = 1
    elif ddl == 'z':
        j = 2
    else:
        raise IOError, "ddl should be chosen among ('x', 'y', 'z')"

    rotvec = np.zeros(3, dtype=float)
    rotvec[j] = angle
    return rotate(V, rotvec)


def scale(V, alpha):
    """scale(V, alpha)

    Scales a mesh

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        alpha: float
            scaling factor

    Returns:
        V: ndarray
            numpy array of the coordinates of the rotated mesh's nodes

    """
    return alpha * V


def flip_normals(F):
    """flip_normals(F)

    Reverse the normals of a mesh

    Parameters:
        F: ndarray
            numpy array of the faces' nodes connectivities

    Returns:
        F: ndarray
            numpy array of the faces' nodes connectivities with
            reversed normals

    """
    return np.fliplr(F)


def symmetrize(V, F, plane):
    """symmetrize(V, F, plane)

    Symmetrize the mesh against a plane

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        plane: Plane
            The symmetry plane

    Returns:
        V: ndarray
            numpy array of the coordinates of the symmetrised
            mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities for
            the symmetrised mesh

    """

    # Symmetrizing the nodes
    nv = V.shape[0]

    normal = plane.normal/np.dot(plane.normal, plane.normal)
    V = np.concatenate((V, V-2*np.outer(np.dot(V, normal)-plane.c, normal)))
    F = np.concatenate((F, np.fliplr(F.copy()+nv)))

    # TODO: be more fine and detect vertices that are duplicated, avoiding using the general merge_duplicates()
    # function
    return merge_duplicates(V, F, verbose=False)

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

    # FIXME : Do we have to repeat the first point at the last position
    # of polygon ???

    x = point[0]
    y = point[1]

    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in xrange(n):
        p2x, p2y = poly[i]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def generate_lid(V, F, max_area=None, verbose=False):
    """generate_lid(V, F, max_area=None, verbose=False)

    Meshes the lid of a mesh with triangular faces to be used in irregular frequency
    removal in BEM softwares. It clips the mesh againt the plane Oxy, extract the intersection
    polygon and relies on meshpy (that is a wrapper around the TRIANGLE meshing library).
    It is able to deal with moonpools.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        max_area[optional]: float
            The maximum area of triangles to be generated
        verbose[optional]: bool
            If set to True, generates output along the processing

    """

    # TODO: Put the reference of TRIANGLE and meshpy (authors...) in the docstring

    # TODO: remove verbose mode and place it into the main of meshmagick !!!

    try:
        import meshpy.triangle as triangle
    except:
        raise ImportError, 'Meshpy has to be available to use the generate_lid() function'

    if verbose:
        print '\n--------------'
        print 'Lid generation'
        print '--------------\n'

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
        sarea = np.array([points[j][0]*points[j+1][1] - points[j+1][0]*points[j][1] for j in xrange(n-1)],
                         dtype=np.float).sum()
        if sarea < 0.:
            holes.append(polygon)
        else:
            boundaries.append(polygon)

    nb_hole = len(holes)
    nb_bound = len(boundaries)

    hole_dict = dict([(j, []) for j in xrange(nb_bound)])
    if nb_hole > 0:
        if verbose:
            if nb_hole == 1:
                word = 'moonpool has'
            else:
                word = 'moonpools have'
            print '\t-> %u %s been detected' % (nb_hole, word)

        # TODO : getting a point inside the hole polygon

        def pick_point_inside_hole(hole):

            # First testing with the geometric center of the hole
            point = np.array(hole).sum(axis=0)/len(hole)
            if not _is_point_inside_polygon(point, hole):
                # Testing something else
                raise RuntimeError, 'The algorithm should be refined to more complex polygon topologies... up to you ?'

            return point

        # Assigning holes to boundaries
        if nb_bound == 1 and nb_hole == 1:
            # Obvious case
            hole_dict[0].append( ( 0, pick_point_inside_hole(V[holes[0]][:, :2]) ) )
        else:
            # We may do a more elaborate search
            for ihole, hole in enumerate(holes):
                P0 = V[hole[0]][:2]
                # Testing against all boundary polygons
                for ibound, bound in enumerate(boundaries):
                    if _is_point_inside_polygon(P0, V[bound][:, :2]):
                        hole_dict[ibound].append( ( ihole, pick_point_inside_hole(V[hole][:, :2]) ) )
                        break

    def round_trip_connect(start, end):
        return [(j, j+1) for j in xrange(start, end)] + [(end, start)]

    # Meshing every boundaries, taking into account holes
    for ibound, bound in enumerate(boundaries):

        nvp = len(bound)-1

        # Building the loop
        points = map(tuple, list(V[bound][:-1, :2]))

        edges = round_trip_connect(0, nvp-1)

        info = triangle.MeshInfo()

        if len(hole_dict[ibound]) > 0:
            for ihole, point in hole_dict[ibound]:
                hole = holes[ihole]
                points.extend(map(tuple, list(V[hole][:-1, :2])))
                edges.extend(round_trip_connect(nvp, len(points)-1))

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

        # show(V, F)
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
        print "\n\t-> %u %s been added successfully\n" % (nb_bound, verb)

    return V, F

def fill_holes(V, F, verbose=False):

    import vtk

    if verbose:
        print "Filling holes"

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
        print "\t--> Done!"

    return V, F


def remove_degenerated_faces(V, F, rtol=1e-5, verbose=False):

    # First computing the faces properties
    areas, normals, centers = get_all_faces_properties(V, F)

    area_threshold = areas.mean()*rtol

    # Detecting faces that have null area
    tol = 1e-5
    F = F[np.logical_not(areas < area_threshold)]

    return F

def triangulate_quadrangles(F, verbose=False):
    """triangulate_quadrangles(V, F)

    Return a copy of the V, F mesh where all quadrangles have been triangulated
    in order to obtain a triangle only mesh. A quadrangle is split into two triangles
    following a dummy rule (may be enhanced in the future).

    Parameters:
        F: ndarray
            numpy array of the faces' nodes connectivities
        verbose[optional]: bool
            If set to True, displays report

    Returns:
        F: ndarray
            numpy array of the faces' nodes connectivities for the triangulated mesh
    """

    # TODO: Ensure the best quality aspect ratio of generated triangles

    # Defining both triangles id lists to be generated from quadrangles
    T1 = (0, 1, 2)
    T2 = (0, 2, 3)

    quad_mask = F[:, 0] != F[:, -1]

    # Triangulation
    new_faces = F[quad_mask].copy()
    new_faces[:, :3] = new_faces[:, T1]
    new_faces[:, -1] = new_faces[:, 0]

    F[quad_mask, :3] = F[:, T2][quad_mask]
    F[quad_mask, -1] = F[quad_mask, 0]

    F = np.concatenate((F, new_faces))

    if verbose:
        nquad = len(new_faces)
        if nquad != 0:
            print '{:d} quadrangles have been split in triangles'.format(nquad)
            print '\t-> Done.\n'

    return F


def detect_features(V, F, verbose=True):

    mesh = Mesh(V, F, verbose=verbose)
    mesh.detect_features(verbose=verbose)

    return


def reformat_triangles(F, verbose=False):
    """reformat_triangles(F)

    Returns a connectivity array F that ensures that triangles
    are described using the rule that the first node id is equal
    to the last id

    Parameters:
        F: ndarray
            numpy array of the faces' nodes connectivities

    Returns:
        F: ndarray
            numpy array of the faces' nodes connectivities with the right rule
            for triangle description

    """

    if verbose:
        print "* Ensuring consistent definition of triangles:"
    quads = F[:, 0] != F[:, -1]
    nquads_init = sum(quads)

    F[quads] = np.roll(F[quads], 1, axis=1)
    quads = F[:, 0] != F[:, -1]

    F[quads] = np.roll(F[quads], 1, axis=1)
    quads = F[:, 0] != F[:, -1]

    F[quads] = np.roll(F[quads], 1, axis=1)
    quads = F[:, 0] != F[:, -1]
    nquads_final = sum(quads)
    if verbose:
        if nquads_final < nquads_init:
            print "\t -> %u triangles were described the wrong way and have been corrected" % (nquads_init-nquads_final)
        else:
            print "\t -> Everything is good :)"

    return F

def _build_polyline(curve):
    import vtk

    npoints = len(curve)

    points = vtk.vtkPoints()
    for point in curve:
        points.InsertNextPoint(point)

    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(npoints)

    for id in xrange(npoints):
        polyline.GetPointIds().SetId(i, i)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(cells)

    return polydata

def show(V, F, normals=False):
    import MMviewer
    polydata = _build_vtkPolyData(V, F)
    my_viewer = MMviewer.MMViewer()
    if normals:
        my_viewer.normals_on()
    my_viewer.add_polydata(polydata)
    my_viewer.show()
    my_viewer.finalize()

# =======================================================================
#                         COMMAND LINE USAGE
# =======================================================================
extension_dict = { #keyword           reader,   writer
                  'mar':             (load_MAR, write_MAR),
                  'nemoh':           (load_MAR, write_MAR),
                  'wamit':           (load_GDF, write_GDF),
                  'gdf':             (load_GDF, write_GDF),
                  'diodore-inp':     (load_INP, write_INP),
                  'inp':             (load_INP, write_INP),
                  'diodore-dat':     (load_DAT, write_DAT),
                  'hydrostar':       (load_HST, write_HST),
                  'hst':             (load_HST, write_HST),
                  'natural':         (load_NAT, write_NAT),
                  'nat':             (load_NAT, write_NAT),
                  'gmsh':            (load_MSH, write_MSH),
                  'msh':             (load_MSH, write_MSH),
                  'rad':             (load_RAD, write_RAD),
                  'radioss':         (load_RAD, write_RAD),
                  'stl':             (load_STL, write_STL),# FIXME: Verifier que ce n'est pas load_STL2
                  'vtu':             (load_VTU, write_VTU),
                  'vtp':             (load_VTP, write_VTP),
                  'paraview-legacy': (load_VTK, write_VTK),# VTK
                  'vtk':             (load_VTK, write_VTK),
                  'tecplot':         (load_TEC, write_TEC),
                  'tec':             (load_TEC, write_TEC)
                  }

def main():

    import argparse
    import sys

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

                    *-----------*------------*-----------------*----------------------*
                    | File      | R: Reading | Software        | Keywords             |
                    | extension | W: writing |                 |                      |
                    *-----------*------------*-----------------*----------------------*
                    |   .mar    |    R/W     | NEMOH (1)       | nemoh, mar           |
                    |   .gdf    |    R/W     | WAMIT (2)       | wamit, gdf           |
                    |   .inp    |    R       | DIODORE (3)     | diodore-inp, inp     |
                    |   .DAT    |    W       | DIODORE (3)     | diodore-dat          |
                    |   .hst    |    R/W     | HYDROSTAR (4)   | hydrostar, hst       |
                    |   .nat    |    R/W     |    -            | natural, nat         |
                    |   .msh    |    R       | GMSH (5)        | gmsh, msh            |
                    |   .rad    |    R       | RADIOSS         | rad, radioss         |
                    |   .stl    |    R/W     |    -            | stl                  |
                    |   .vtu    |    R/W     | PARAVIEW (6)    | vtu                  |
                    |   .vtp    |    R/W     | PARAVIEW (6)    | vtp                  |
                    |   .vtk    |    R/W     | PARAVIEW (6)    | paraview-legacy, vtk |
                    |   .tec    |    R/W     | TECPLOT (7)     | tecplot, tec         |
                    *---------- *-----------------------------------------------------*

                    By default, Meshmagick uses the filename extensions to choose the
                    appropriate reader/writer. This behaviour might be bypassed using the
                    -ifmt and -ofmt optional arguments. When using these options, keywords
                    defined in the table above must be used as format identifiers.

                    (1) NEMOH is an open source BEM Software for seakeeping developped at
                        Ecole Centrale de Nantes (LHHEA)
                    (2) WAMIT is a BEM Software for seakeeping developped by WAMIT, Inc.
                    (3) DIODORE is a BEM Software for seakeeping developped by PRINCIPIA
                    (4) HYDROSTAR is a BEM Software for seakeeping developped by
                        BUREAU VERITAS
                    (5) GMSH is an open source meshing software developped by C. Geuzaine
                    and J.-F. Remacle
                    (6) PARAVIEW is an open source visualization software developped by
                        Kitware
                    (7) TECPLOT is a visualization software developped by Tecplot


                    """,
        epilog='--  Copyright 2014-2015  -  Francois Rongere  /  Ecole Centrale de Nantes  --',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('infilename', # TODO : voir pour un typ=file pour tester l'existence
                        help='path of the input mesh file in any format')

    parser.add_argument('-o', '--outfilename', type=str,
                        help="""path of the output mesh file. The format of
                         this file is determined from the extension given
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

    parser.add_argument('-t', '--translate',
                        nargs=3, type=float,
                        help="""translates the mesh in 3D
                        Usage -translate tx ty tz""")

    parser.add_argument('-tx', '--translatex', type=float,
                        help="""translates the mesh following the x direction""")

    parser.add_argument('-ty', '--translatey', type=float,
                        help="""translates the mesh following the y direction""")

    parser.add_argument('-tz', '--translatez', type=float,
                        help="""translates the mesh following the z direction""")

    parser.add_argument('-r', '--rotate',
                        nargs=3, type=float,
                        help="""rotates the mesh in 3D""")

    parser.add_argument('-rx', '--rotatex', type=float,
                        help="""rotates the mesh around the x direction""")

    parser.add_argument('-ry', '--rotatey', type=float,
                        help="""rotates the mesh around the y direction""")

    parser.add_argument('-rz', '--rotatez', type=float,
                        help="""rotates the mesh around the z direction""")

    parser.add_argument('-s', '--scale', type=float,
                        help="""scales the mesh. CAUTION : if used along
                         with a translation option, the scaling is done before
                        the translations. The translation magnitude should be set
                        accordingly to the newly scaled mesh.
                        """)

    parser.add_argument('-hn', '--heal_normals', action='store_true',
                        help="""Checks and heals the normals consistency and
                        verify if they are outward.
                        """)

    parser.add_argument('-fn', '--flip-normals', action='store_true',
                        help="""flips the normals of the mesh""")

    parser.add_argument('-p', '--plane', nargs='+', action='append',
                        help="""Defines a plane used by the --clip and --symmetrize options.
                        It can be defined by the floats nx ny nz c where [nx, ny, nz]
                        is a normal vector to the plane and c defines its position
                        following the equation <N|X> = c with X a point belonging
                        to the plane.
                        It can also be defined by a string among [Oxy, Oxz, Oyz, \Oxy, \Oxz, \Oyz]
                        for quick definition. Several planes may be defined on the same command
                        line. Planes with a prepended '\' have normals inverted i.e. if Oxy has its
                        normal pointing upwards, \Oxy has its normal pointing downwards.
                        In that case, the planes are indexed by an integer starting by
                        0 following the order given in the command line.
                        """)

    parser.add_argument('-c', '--clip', nargs='*', action='append',
                        help="""cuts the mesh with a plane. Is no arguments are given, the Oxy plane
                        is used. If an integer is given, it should correspond to a plane defined with
                        the --plane option. If a key string is given, it should be a valid key (see
                        help of --plane option for valid plane keys). A normal and a scalar could
                        also be given for the plane definition just as for the --plane option. Several
                        clipping planes may be defined on the same command line.""")

    parser.add_argument('-md', '--merge-duplicates', nargs='?', const='1e-8', default=None,
                        help="""merges the duplicate nodes in the mesh with the absolute tolerance
                        given as argument (default 1e-8)""")

    parser.add_argument('-tq', '--triangulate_quadrangles', action='store_true',
                        help="""Triangulate all quadrangle faces by a simple splitting procedure.
                        Twho triangles are generated and from both solution, the one with the best
                        aspect ratios is kept. This option may be used in conjunction with a
                        mesh export in a format that only deal with triangular cells like STL format.""")

    parser.add_argument('-sym', '--symmetrize', nargs='*', action='append',
                        help="""Symmetrize the mesh by a plane defined wether by 4 scalars, i.e.
                        the plane normal vector coordinates and a scalar c such as N.X=c is the
                        plane equation (with X a point of the plane) or a string among Oxz, Oxy
                        and Oyz which are shortcuts for planes passing by the origin and whose
                        normals are the reference axes. Default is Oxz if only -y is specified.
                        Be careful that symmetry is applied before any rotation so as the plane
                        equation is defined in the initial frame of reference.""")

    # Arguments concerning the hydrostatics
    # TODO : permettre de specifier un fichier de sortie
    parser.add_argument('-hs', '--hydrostatics', action='store_true',
                        help="""Compute hydrostatics. When used with the --verbose option, hydrostatic
                        data will be shown on the command line.""") # TODO : completer l'aide avec les options

    # parser.add_argument('--mass', default=None, type=float,
    #                     help="""Specifies the mass of the device for hydrostatics calculations.
    #                     When used, an hydrostatic equilibrium is resolved. Depending on the join
    #                     use of --mass, --cog and --zcog options, several behavior may be obtained.
    #
    #                     1) --mass is given a value:
    #
    #                         This options trigs an equilibrium resolution.
    #
    #                         a) No options --cog or --zcog are given:
    #                             Equilibrium is searched in z only, no hydrostatic stiffness
    #                             matrix is computed but the stiffness in heave.
    #                         b) --zcog is given alone:
    #                             Equilibrium is searched in z only, the hydrostatics stiffness
    #                             matrix is given.
    #                         c) --cog is given alone:
    #                             Equilibrium is searched in 6 dof and the hydrostatics stiffness
    #                             matrix is given
    #
    #                         Note that if both options --cog and --zcog are given, --zcog will
    #                         be ignored but taken from --cog and the result will be the same as c).
    #
    #                     2) --mass is used:
    #
    #                         No equilibrium resolution is performed if no mass is specified.
    #                         Mesh position is considered as being the equilibrium position and
    #                         hydrostatics data are generated as-is. Mass of the device is then
    #                         an output and is considered to be the computed displacement * rho_water.
    #
    #                         a) No options --cog or --zcog are given:
    #                             Only the stiffness in heave can be calculated
    #                         b) --zcog is given alone:
    #                             The full stiffness matrix is calculated and the estimated cog
    #                             position is given.
    #                         c) --cog is given alone:
    #                             Idem b) but the residual may be not identically 0. as the mesh may
    #                             not be at equilibrium.
    #
    #                         Note that if both options --cog and --zcog are given, --zcog will
    #                         be ignored but taken from --cog and the result will be the same as c).
    #                     """)

    # parser.add_argument('--cog', nargs=3, default=None, type=float,
    #                     help="""Specifies the center of gravity to be used along with the
    #                     --hydrostatics option. See --mass for more information. See also the
    #                     --inertias option for side effects.
    #                     """)

    parser.add_argument('--zcog', default=None, type=float,
                        help="""Specifies the z position of the center of gravity to be used along
                        the --hydrostatics option.
                        """)

    # parser.add_argument('--zcog', default=None, type=float,
    #                     help="""Specifies the z position of the center of gravity to be used along
    #                     the --hydrostatics option. See --mass for more information. See also the
    #                     --inertias option for side effects.
    #                     """)

    # parser.add_argument('--anim', action='store_true',
    #                     help="""Generates animation files for paraview visualization of the equilibrium
    #                     resolution.
    #                     """)

    parser.add_argument('--rho-water', default=1023., type=float,
                        help="""Specifies the density of salt water. Default is 1023 kg/m**3.
                        """)

    parser.add_argument('-g', '--grav', default=9.81, type=float,
                        help="""Specifies the acceleration of gravity on the earth surface.
                        Default is 9.81 m/s**2.
                        """)

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

    parser.add_argument('--lid', nargs='?', const=1., default=None, type=float,
                        help="""Generate a triangle mesh lid on the mesh clipped by the Oxy plane.
                        """)

    parser.add_argument('--fill-holes', '-fh', action='store_true',
                        help="""Fill little holes by triangulation if any.
                        """)

    parser.add_argument('--show', action='store_true',
                        help="""Shows the input mesh in an interactive window""")

    parser.add_argument('--version', action='version',
                        version='meshmagick - version %s\n%s'%(__version__, __copyright__),
                        help="""Shows the version number and exit""")


    if acok:
        argcomplete.autocomplete(parser)

    # TODO : Utiliser des sous-commandes pour l'utilisation de meshmagick

    args, unknown = parser.parse_known_args()

    if args.quiet:
        verbose = False
    else:
        verbose = True

    if verbose:
        print '\n============================================='
        print 'meshmagick - version %s\n%s' % (__version__, __copyright__)
        print '============================================='

    # LOADING DATA FROM FILE
    if args.input_format is not None:
        format = args.input_format
    else:
        # Format based on extension
        _ ,ext = os.path.splitext(args.infilename)
        format = ext[1:].lower()
        if format == '':
            raise IOError, 'Unable to determine the input file format from its extension. Please specify an input format.'

    # Loading mesh elements from file
    if os.path.isfile(args.infilename):
        V, F = load_mesh(args.infilename, format)
        if verbose:
            print '%s successfully loaded' % args.infilename
    else:
        raise IOError, 'file %s not found'%args.infilename

    # Ensuring triangles are following the right convention (last id = first id)
    F = reformat_triangles(F)

    # Merge duplicate vertices
    if args.merge_duplicates is not None:
        tol = float(args.merge_duplicates)
        if verbose:
            print '\nOPERATION: Merge duplicate nodes'
        V, F = merge_duplicates(V, F, verbose=verbose, tol=tol)
        if verbose:
            print '\t-> Done.'

    # TODO : put that dict at the begining of the main function
    plane_str_list = {'Oxy':[0.,0.,1.],
                      'Oxz':[0.,1.,0.],
                      'Oyz':[1.,0.,0.],
                      '\Oxy':[0.,0.,-1.],
                      '\Oxz':[0.,-1.,0.],
                      '\Oyz':[-1.,0.,0.]}

    # Defining planes
    if args.plane is not None:
        nb_planes = len(args.plane)

        if verbose:
            if nb_planes == 1:
                verb = 'plane has'
            else:
                verb = 'planes have'
            print '\n%u %s been defined' % (nb_planes, verb)

        planes = [Plane() for i in xrange(nb_planes)]
        for (iplane, plane) in enumerate(args.plane):
            if len(plane) == 4:
                # plane is defined by normal and scalar
                try:
                    planes[iplane].normal = np.array(map(float, plane[:3]), dtype=np.float)
                    planes[iplane].c = float(plane[3])
                except:
                    raise AssertionError, 'Defining a plane by normal and scalar requires four scalars'

            elif len(plane) == 1:
                if plane_str_list.has_key(plane[0]):
                    planes[iplane].normal = np.array(plane_str_list[plane[0]], dtype=np.float)
                    planes[iplane].c = 0.
                else:
                    raise AssertionError, '%s key for plane is not known. Choices are [%s].' % (plane[0], ', '.join(plane_str_list.keys()) )
            else:
                raise AssertionError, 'Planes should be defined by a normal and a scalar or by a key to choose among [%s]' % (', '.join(plane_str_list.keys()))



    # Symmetrizing the mesh
    if args.symmetrize is not None:
        nb_sym = len(args.symmetrize)

        if verbose:
            if nb_sym == 1:
                verb = 'plane'
            else:
                verb = 'planes'
            print '\nMesh is being symmetrized by %u %s' % (nb_sym, verb)

        sym_plane = Plane()
        for plane in args.symmetrize:
            if len(plane) == 0:
                # Default symmetry by plane Oxz
                sym_plane.normal = np.array([0., 1., 0.], dtype=np.float)
                sym_plane.c = 0.
            elif len(plane) == 1:
                try:
                    # Plane ID
                    plane_id = int(plane[0])
                    if plane_id < nb_planes:
                        sym_plane = planes[plane_id]
                    else:
                        raise AssertionError, 'Plane with ID %u has not been defined with option --plane' % plane_id
                except:
                    # A key string
                    if plane_str_list.has_key(plane[0]):
                        sym_plane.normal = np.asarray(plane_str_list[plane[0]], dtype=np.float)
                        sym_plane.c = 0.
                    else:
                        raise AssertionError, 'Planes should be defined by a normal and a scalar or by a key to choose among [%s]' % (', '.join(plane_str_list.keys()))
            elif len(plane) == 4:
                # Plane defined by a normal and a scalar
                try:
                    sym_plane.normal = np.array(map(float, plane[:3]), dtype=np.float)
                    sym_plane.c = float(plane[3])
                except:
                    raise AssertionError, 'Defining a plane by normal and scalar requires four scalars'
            else:
                raise AssertionError, 'Unknown mean to define a plane for symmetry'
            V, F = symmetrize(V, F, sym_plane)
            if verbose:
                print '\t-> Done.'

    # Heal normals
    if args.heal_normals:
        if verbose:
            print '\nOPERATION: heal normals'
        F = heal_normals(V, F, verbose=verbose)
        if verbose:
            print '\t-> Done.'

    # Mesh translations
    if args.translate is not None:
        if verbose:
            print '\nOPERATION: Translation by [%f, %f, %f]' % tuple(args.translate)
        V = translate(V, args.translate)
        if verbose:
            print '\t-> Done.'

    if args.translatex is not None:
        if verbose:
            print '\nOPERATION: Translation by %f along X' % args.translatex
        V = translate_1D(V, args.translatex, 'x')
        if verbose:
            print '\t-> Done.'

    if args.translatey is not None:
        if verbose:
            print '\nOPERATION: Translation by %f along Y' % args.translatey
        V = translate_1D(V, args.translatey, 'y')
        if verbose:
            print '\t-> Done.'

    if args.translatez is not None:
        if verbose:
            print '\nOPERATION: Translation by %f along Z' % args.translatez
        V = translate_1D(V, args.translatez, 'z')
        if verbose:
            print '\t-> Done.'

    # Mesh rotations
    # FIXME : supprimer le cast angles et ne prendre que des degres
    if args.rotate is not None:
        if verbose:
            print '\nOPERATION: Rotation by [%f, %f, %f]' % tuple(args.rotate)
        V = rotate(V, args.rotate*math.pi/180.)
        if verbose:
            print '\t-> Done.'

    if args.rotatex is not None:
        if verbose:
            print '\nOPERATION: Rotation by %f around X (Roll)' % args.rotatex
        V = rotate_1D(V, args.rotatex*math.pi/180., 'x')
        if verbose:
            print '\t-> Done.'

    if args.rotatey is not None:
        if verbose:
            print '\nOPERATION: Rotation by %f around Y (Pitch)' % args.rotatey
        V = rotate_1D(V, args.rotatey*math.pi/180., 'y')
        if verbose:
            print '\t-> Done.'

    if args.rotatez is not None:
        if verbose:
            print '\nOPERATION: Rotation by %f around Z (Yaw)' % args.rotatez
        V = rotate_1D(V, args.rotatez*math.pi/180., 'z')
        if verbose:
            print '\t-> Done.'

    if args.scale is not None:
        if verbose:
            print '\nOPERATION: Scaling by %f' % args.scale
        V = scale(V, args.scale)
        if verbose:
            print '\t-> Done.'

    if args.flip_normals:
        if verbose:
            print '\nOPERATION: Flipping normals'
        F = flip_normals(F)
        if verbose:
            print '\t-> Done.'

    if args.triangulate_quadrangles:
        F = triangulate_quadrangles(F, verbose=verbose)

    # Clipping the mesh
    if args.clip is not None:
        clipping_plane = Plane()
        nb_clip = len(args.clip)

        if verbose:
            if nb_clip == 1:
                verb = 'plane'
            else:
                verb = 'planes'
            print '\nMesh is being clipped by %u %s' % (nb_clip, verb)

        for plane in args.clip:
            if len(plane) == 0:
                # Default clipping plane Oxy
                clipping_plane.normal = np.array([0., 0., 1.], dtype=np.float)
                clipping_plane.c = 0.
            elif len(plane) == 1:
                try:
                    # Plane ID
                    plane_id = int(plane[0])
                    if plane_id < nb_planes:
                        clipping_plane = planes[plane_id]
                    else:
                        raise AssertionError, 'Plane with ID %u has not been defined with option --plane' % plane_id
                except:
                    # A key string
                    if plane_str_list.has_key(plane[0]):
                        clipping_plane.normal = np.asarray(plane_str_list[plane[0]], dtype=np.float)
                        clipping_plane.c = 0.
                    else:
                        raise AssertionError, 'Planes should be defined by a normal and a scalar or by a key to choose among [%s]' % (', '.join(plane_str_list.keys()))
            elif len(plane) == 4:
                # Plane defined by a normal and a scalar
                try:
                    clipping_plane.normal = np.array(map(float, plane[:3]), dtype=np.float)
                    clipping_plane.c = float(plane[3])
                except:
                    raise AssertionError, 'Defining a plane by normal and scalar requires four scalars'
            else:
                raise AssertionError, 'Unknown mean to define a plane for clipping'
            V, F = clip_by_plane(V, F, clipping_plane)
            if verbose:
                print '\t-> Done.'

    # Compute principal inertia parameters
    # if args.inertias:
    #     # TODO : completer l'aide avec la logique de cette fonction !!
    #     if verbose:
    #         print '\n------------------'
    #         print 'Computing inertias'
    #         print '------------------'
    #     if args.no_hull:
    #         hull = False
    #     else:
    #         hull = True
    #
    #     mass, cog, inertia_matrix = get_inertial_properties(V, F,
    #                                     rho=args.rho_medium,
    #                                     mass=args.mass,
    #                                     thickness=args.thickness,
    #                                     shell=hull,
    #                                     verbose=verbose)
    #     # Replacing values in command line arguments in the eventuality of hydrostatics computations
    #     args.mass = mass
    #     args.cog = cog
    #     if verbose:
    #         print '\t-> Done.'

    # if args.gz_curves is not None:
    #     if verbose:
    #         print '\n-------------------'
    #         print 'Computing GZ curves'
    #         print '-------------------'
    #
    #     spacing = args.gz_curves
    #     try:
    #         import hydrostatics as hs
    #     except:
    #         raise ImportError, '--gz-curves option relies on the hydrostatics module that can not be found'
    #
    #     # if args.hydrostatics:
    #     #     raise RuntimeError, """GZ computations can not be performed at the same time as a hydrostatics equilibrium
    #     #                            resolution as it needs a full mesh to perform clipping at different angles"""
    #
    #     if args.zcog is None:
    #         raise RuntimeError, 'For the GZ computations, the --zcog option is mandatory'
    #
    #     hsMesh = hs.HydrostaticsMesh(V, F, rho_water=args.rho_water, g=args.grav)
    #     hs.get_GZ_curves(hsMesh, args.zcog,
    #                      spacing=spacing,
    #                      rho_water=args.rho_water,
    #                      g=args.grav,
    #                      verbose=verbose)
    #     if verbose:
    #         print '\t-> Done.'


    # Compute hydrostatics
    # ESSAI
    if args.hydrostatics:
        try:
            import hydrostatics as hs
        except:
            raise ImportError, '--hydrostatics option relies on the hydrostatics module that can not be found'

        if args.zcog is None:
            raise ValueError, 'The hydrostatics option shall be used along with the --zcog option for the hydrostatic stiffness matrix to be computed'

        outputHS = hs.compute_hydrostatics(V, F, args.zcog, verbose=verbose)
        V = outputHS['Vc']
        F = outputHS['Fc']

    # if args.hydrostatics:
    #     try:
    #         import hydrostatics as hs
    #     except:
    #         raise ImportError, '--hydrostatics option relies on the hydrostatics module that can not be found'
    #
    #     if args.anim:
    #         anim=True
    #     else:
    #         anim=False
    #
    #     if args.cog is None:
    #         cog = None
    #     else:
    #         cog = np.asarray(args.cog, dtype=np.float)
    #
    #     # TODO : Revoir la structure afin de ne jouer que sur l'objet !!
    #     hsMesh = hs.HydrostaticsMesh(V, F, rho_water=args.rho_water, g=args.grav)
    #     V, F = hs.get_hydrostatics(hsMesh,
    #                                mass=args.mass,
    #                                zcog=args.zcog,
    #                                cog=cog,
    #                                rho_water=args.rho_water,
    #                                g=args.grav,
    #                                anim=anim,
    #                                verbose=verbose)[:2]


    # Lid generation on a clipped mesh
    if args.lid is not None:
        V, F = generate_lid(V, F, max_area=args.lid, verbose=verbose)

    if args.fill_holes:
        V, F = fill_holes(V, F, verbose=verbose)


    # WARNING : No more mesh modification should be released from this point until the end of the main

    if args.info:
        get_info(V, F)

    if args.show:
        show(V, F)

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
                    if not extension_dict.has_key(format):
                        raise IOError, 'Could not determine a format from input file extension, please specify an input format or an extension'
            else:
                format = os.path.splitext(args.outfilename)[1][1:].lower()

        if verbose:
            print 'Writing %s' % args.outfilename
        write_mesh(args.outfilename, V, F, format)
        if verbose:
            print '\t-> Done.'

    if verbose:
        print '\n============================================================='
        print 'meshmagick - version %s\n%s' % (__version__, __copyright__)
        print 'Maintainer : %s <%s>' % (__maintainer__, __email__)
        print 'Good Bye!'
        print '============================================================='

if __name__ == '__main__':
    main()
