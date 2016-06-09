#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from mesh import *


class MeshClipper(object):
    def __init__(self, source_mesh=None, plane=None, vicinity_tol=1e-4, assert_closed_boundaries=False, verbose=False):
        source_mesh._verbose = verbose
        self._source_mesh = source_mesh
        self._plane = plane

        self._vicinity_tol = vicinity_tol

        self._assert_closed_boundaries = assert_closed_boundaries
        self._verbose = verbose

        self._internal_data = dict()

        self._update()

    @property
    def verbose(self):
        return self._verbose

    def verbose_on(self):
        self._verbose = True

    def verbose_off(self):
        self._verbose = False

    @property
    def assert_closed_boundaries(self):
        return self._assert_closed_boundaries

    def assert_closed_boundaries_on(self):
        self._assert_closed_boundaries = True
        return

    def assert_closed_boundaries_off(self):
        self._assert_closed_boundaries = False
        return

    @property
    def vicinity_tol(self):
        return self._vicinity_tol

    @vicinity_tol.setter
    def vicinity_tol(self, value):
        self._internal_data.clear()
        self._vicinity_tol = float(value)
        self._update()

    @property
    def source_mesh(self):
        return self._source_mesh

    @source_mesh.setter
    def source_mesh(self, value):
        self._source_mesh = value
        self._internal_data.clear()
        self._update()
        return

    @property
    def plane(self):
        return self._plane

    @plane.setter
    def plane(self, value):
        self._internal_data.clear()
        self._plane = value
        self._update()
        return

    def _update(self):
        self._vertices_positions_wrt_plane()
        self._partition_mesh()
        self.clip()
        return

    def _vertices_positions_wrt_plane(self):
        vertices_distances = self._plane.get_point_dist_wrt_plane(self._source_mesh._V)

        vertices_positions = {'vertices_distances': vertices_distances,
                              'vertices_above_mask': vertices_distances > self._vicinity_tol,
                              'vertices_on_mask': np.fabs(vertices_distances) < self._vicinity_tol,
                              'vertices_below_mask': vertices_distances < -self._vicinity_tol
                              }
        self._internal_data.update(vertices_positions)
        return

    def _partition_mesh(self):
        vertices_distances = self._internal_data['vertices_distances']
        vertices_above_mask = self._internal_data['vertices_above_mask']
        vertices_on_mask = self._internal_data['vertices_on_mask']
        vertices_below_mask = self._internal_data['vertices_below_mask']

        source_mesh_faces = self.source_mesh.F

        nb_vertices_above = vertices_above_mask[source_mesh_faces].sum(axis=1)
        nb_vertices_below = vertices_below_mask[source_mesh_faces].sum(axis=1)

        # Simple criteria ensuring that faces are totally above or below the plane (4 vertices at the same side)
        # Works for both triangles and quadrangles
        above_faces_mask = nb_vertices_above == 4
        below_faces_mask = nb_vertices_below == 4
        crown_faces_mask = np.logical_and(np.logical_not(above_faces_mask), np.logical_not(below_faces_mask))

        above_faces_ids = np.where(above_faces_mask)[0]
        below_faces_ids = np.where(below_faces_mask)[0]
        crown_faces_ids = np.where(crown_faces_mask)[0]

        partition = dict()
        def generate_mesh(faces_ids, key):
            new_mesh, ids = self._source_mesh.extract_faces(faces_ids, return_index=True)
            new_mesh.name = key
            partition['_'.join((key, 'vertices_ids'))] = ids
            partition['_'.join((key, 'vertices_distances'))] = vertices_distances[ids]
            partition['_'.join((key, 'above_vertices_mask'))] = vertices_above_mask[ids]
            partition['_'.join((key, 'on_vertices_mask'))] = vertices_on_mask[ids]
            partition['_'.join((key, 'below_vertices_mask'))] = vertices_below_mask[ids]
            partition[key] = new_mesh

        # Generating partition meshes
        generate_mesh(above_faces_ids, 'upper_mesh')
        generate_mesh(crown_faces_ids, 'crown_mesh')
        generate_mesh(below_faces_ids, 'lower_mesh')
        partition['above_faces_ids'] = above_faces_ids
        partition['crown_faces_ids'] = crown_faces_ids
        partition['below_faces_ids'] = below_faces_ids

        self._internal_data.update(partition)

        return


    @property
    def lower_mesh(self):
        return self._internal_data['lower_mesh']

    @property
    def crown_mesh(self):
        return self._internal_data['crown_mesh']

    @property
    def upper_mesh(self):
        return self._internal_data['upper_mesh']

    @property
    def closed_polygons(self):
        return self._internal_data['closed_polygons']

    @property
    def nb_closed_polygons(self):
        return len(self._internal_data['closed_polygons'])

    @property
    def open_lines(self):
        return self._internal_data['open_lines']

    @property
    def nb_open_lines(self):
        return len(self._internal_data['open_lines'])

    @property
    def clipped_crown_mesh(self):
        return self._internal_data['clipped_crown_mesh']

    # @cached_property
    def clip_crown_by_plane(self):

        crown_mesh = self.crown_mesh
        vertices = crown_mesh.V
        vertices_on_mask = self._internal_data['crown_mesh_on_vertices_mask']

        # TODO: Vertices pre-projection
        # vertices_on = partition['vertices_on']
        # V[vertices_on] = plane.orthogonal_projection_on_plane(V[vertices_on])
        # pos[vertices_on] = 0.

        vertices_above_mask = self._internal_data['crown_mesh_above_vertices_mask']
        vertices_below_mask = self._internal_data['crown_mesh_below_vertices_mask']

        vertices_distances = self._internal_data['crown_mesh_vertices_distances']

        # Init
        crown_faces = list()
        direct_boundary_edges = dict()
        inv_boundary_edges = dict()
        intersections = list()

        nI = crown_mesh.nb_vertices

        for face_id in xrange(crown_mesh.nb_faces):

            face = crown_mesh.get_face(face_id)

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
                P0, P1, P2, P3 = vertices[face]
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
                P0, P1, P3 = vertices[face[[0, 1, 3]]]
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
                P0, P1, P3 = vertices[face[[0, 1, 3]]]
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
                P0, P1, P2 = vertices[face]
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
                P0, P1, P2 = vertices[face]
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
                if vertices_distances[face[1]] < 0.:
                    P1, P2 = vertices[face[[1, 2]]]
                    Iright = plane.get_edge_intersection(P1, P2)
                    intersections.append(Iright)
                    boundary_edge = [face[0], nI]
                    crown_faces.append([face[0], face[1], nI, face[0]])
                else:
                    P2, P3 = vertices[face[[2, 3]]]
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
                if vertices_distances[face[1]] < 0.:
                    P2, P3 = vertices[face[[2, 3]]]
                    Iright = plane.get_edge_intersection(P2, P3)
                    intersections.append(Iright)
                    boundary_edge = [face[0], nI]
                    crown_faces.append([face[0], face[1], face[2], nI])
                else:
                    P1, P2 = vertices[face[[1, 2]]]
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
                P1, P2 = vertices[face[[1, 2]]]
                if vertices_distances[face[1]] < 0.:
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

        vertices = np.concatenate((vertices, intersections))

        clipped_crown_mesh = Mesh(vertices, crown_faces)

        # TODO: faire un merge uniquement sur la liste instersections et non sur tout le maillage clipped_crown
        newID = clipped_crown_mesh.merge_duplicates(return_index=True)

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

                if self._assert_closed_boundaries:
                    if len(open_lines) > 0:
                        # mm.write_VTP('mesh_clip.vtp', vertices, crown_faces)
                        raise RuntimeError, 'Open intersection curve found'

                break
        boundaries = {
            'closed_polygons': closed_polygons,
            'open_lines': open_lines
        }

        output = {'clipped_crown_mesh': clipped_crown_mesh,
                  'closed_polygons': closed_polygons,
                  'open_lines': open_lines}

        self._internal_data.update(output)
        return

    @property
    def clipped_mesh(self):
        return self._internal_data['clipped_mesh']


    def clip(self):
        self.clip_crown_by_plane()
        self._internal_data['clipped_mesh'] = self.lower_mesh + self.clipped_crown_mesh
        return

if __name__ == '__main__':

    import mmio

    V, F = mmio.load_VTP('SEAREV.vtp')
    mymesh = Mesh(V, F)

    plane = Plane()

    clipper = MeshClipper(mymesh, plane)


    # lower_mesh = clipper.lower_mesh
    # upper_mesh = clipper.upper_mesh
    # crown_mesh = clipper.crown_mesh
    #
    #
    # mmio.write_VTP('lower.vtp', lower_mesh.V, lower_mesh.F)
    # mmio.write_VTP('crown.vtp', crown_mesh.V, crown_mesh.F)
    # mmio.write_VTP('upper.vtp', upper_mesh.V, upper_mesh.F)
    #
    # print lower_mesh.faces_areas.sum() + clipped_crown_mesh.faces_areas.sum()
    # # print
    # print clipped_mesh.faces_areas.sum()

    for iter in xrange(100):
        print iter
        thetax, thetay = np.random.rand(2)*2*math.pi
        plane.rotate_normal(thetax, thetay)
        clipper.plane = plane
        # clipper.clipped_mesh.show()
        # print clipper.lower_mesh
        # mmio.write_VTP('mesh_%u.vtp'%iter, clipper.clipped_mesh.V, clipper.clipped_mesh.F)


    # print clipped_mesh.faces_areas.sum()
    # clipped_mesh = clipper.clip()
    # clipped_mesh.show()
    #
    # plane.normal = [1, 1, 1]
    # clipper.plane = plane
    #
    # clipped_mesh = clipper.clip()
    # clipped_mesh.show()
