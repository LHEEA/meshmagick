#!/usr/bin/env python
#  -*- coding: utf-8 -*-
import numpy as np
import meshmagick as mm
import math


class MeshHE:
    """
    Class to detect features of mesh using a half-edge data structure
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
        V, F = mm.clean_mesh(V, F, verbose=verbose)

        self.vertices = V

        self.nv = V.shape[0]
        self.nf = F.shape[0]

        # Partitionning _faces array with triangles first and quadrangles last

        # NOTE: Here we modify the order of _faces. See in the future if it is a problem.
        # It is of practical interest when it comes to generate the half-edge data structure
        # But it could be possible to order _faces like that temporarily, perform fast generation
        # of the half-edge data structure and then get back to the initial _faces ordering...
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

        return

    def detect_features(self,
                        thetaf=10, # Threshold for l-strongness in DA of half-edges
                        thetaF=65, # Threshold for u-strongness in DA of edges
                        thetaD=40, # Threshold for u-strongness in AD of _vertices (sharpness)
                        thetat=20, # Threshold for u-strongness in OSTA of edges
                        thetaT=40, # Threshold for u-strongness in TA of _vertices
                        thetae=25, # Threshold for e-strongness in DA of edges
                        thetak=50, # Threshold for obsurity of curves
                        k=5,       # Minimum length of curves for not being obscure
                        verbose=False):

        if verbose:
            print "\nDetecting features of the mesh"

        la = np.linalg
        epsilon = math.tan(math.radians(thetaf)/2)**2
        pi = math.pi
        acos = math.acos

        # generating _faces properties
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

        # Candidate _vertices are ICH.keys()...

        # This is for debug of the algorithm
        # if debug:
        #     _tmp_singleton_HE = []
        #     _tmp_dangling_HE = []
        #     _tmp_semi_joint = []
        #     _tmp_disjoint = []
        #     _tmp_multi_joint_HE = []

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
                            # if debug:
                            #     _tmp_multi_joint_HE.append(iHE)
                            is_disjoint = False
                            end_HE.append(iHE)

            # if debug:
            #     self._write_HE('candidate_HE.vtp', candidate_HE)
            #     self._write_HE('end_HE.vtp', end_HE)
            #     self._write_HE('multi_joint.vtp', _tmp_multi_joint_HE)
            #     self._write_HE('dangling.vtp', _tmp_dangling_HE)
            #     self._write_HE('singleton.vtp', _tmp_singleton_HE)
            #     self._write_HE('semi-joint.vtp', _tmp_semi_joint)
            #     self._write_HE('disjoint.vtp', _tmp_disjoint)
            #
            #     sharp_HE = self.edges[np.where(sharp_edge)[0]]
            #     self._write_HE('sharp_HE.vtp', sharp_HE)
            #     self._write_HE('obscure_end_HE.vtp', obscure_end_HE)

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


        print ICH

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
                    # We're done, every _faces have been visited :)
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
            V, F = mm.extract_faces(self.vertices, self.faces, surface)
            polydata = mm._build_vtkPolyData(V, F)
        #     if debug:
        #         _tmp_viewer.add_polydata(polydata, color=Color(pick_for=polydata).get_rgb())
        # if debug:
        #     _tmp_viewer.show()
        #     _tmp_viewer.finalize()


        # long_surface = []
        # for i, surface in enumerate(surfaces):
        #     long_surface += surface
        #     _vertices, _faces = extract_faces(self._vertices, self._faces, long_surface)
        #     write_VTP('surf%u.vtp'%i, _vertices, _faces)

        if verbose:
            print "\t-> Features detected!"

        return

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

        return

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
        self.areas, self.normals, self.centers = mm.get_all_faces_properties(self.vertices, self.faces)
        return

    def show(self):
        mm.show(self.vertices, self.faces)
        return

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


if __name__ == '__main__':

    V, F = mm.load_VTP('SEAREV.vtp')

    HE_mesh = MeshHE(V, F, verbose=True)

    HE_mesh.detect_features(verbose=True)
