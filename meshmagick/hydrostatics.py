__author__ = 'frongere'

import meshmagick as mm
import numpy as np


class HydrostaticsMesh:
    def __init__(self, V, F, rho_water=1023.):
        self.V = V
        self.F = F
        self.rho_water = rho_water

        # TODO : ne conserver ici que les operation d'initialisation, on mettra les operations de decoupe et de
        # calcul des pptes dans un fonction d'update qui se chargera de faire que les infos sur le maillage coupe
        # soient bien a jour.

        props = mm.get_all_faces_properties(V, F)
        self.areas, self.normals, self.centers = props

        # Defining the first plane
        self.plane = mm.Plane()

        self.clipped_V, self.clipped_F, polygons, props = \
            mm.clip_by_plane(V, F, self.plane, get_polygon=True, props=props)

        if len(polygons) == 0:
            raise RuntimeError, 'could not compute any intersection polygon'

        self.clipped_areas, self.clipped_normals, self.clipped_centers = props

        # Projecting the intersection polygons to the plane coordinates
        self.intersection_vertices = []
        for polygon in polygons:
            self.intersection_vertices.append(self.plane.coord_in_plane(self.clipped_V[polygon]))

        # Updating volume integrals
        self.volint = mm._get_volume_integrals(self.clipped_V, self.clipped_F)
        self.sfint = self._get_surface_integrals()

        # Computing the flottation surface area
        self.sf = self.sfint[0]
        # Computing the immersed volume
        self.Vw = self.get_vw()

        # Computing the volume integrals




        self.Vw = self.get_displacement()

    def get_vw(self):

        R11 = self.plane.Re0[0,0]
        R21 = self.plane.Re0[1,0]
        Vw = self.volint[0] + self.plane.normal[0] * (R11*self.sfint[1] + R21*self.sfint[2] +
                                                      self.plane.e*self.plane.normal[0]*self.sf)

        return Vw

    def get_buoyancy_center(self):
        pass

    def get_displacement(self):
        return self.rho_water * self.Vw

    # def clip_by_plane(self, plane):
    #     mm.clip_by_plane(self.V, self.F, plane, get_polygon=True)

    def get_wetted_surface(self):
        return np.sum(self.clipped_areas)

    def _get_surface_integrals(self):

        mult = np.array([1/2., 1/6., -1/6., 1/12., 1/12., -1/12.], dtype=float) # FIXME : a completer

        sint = np.zeros(6, dtype=float)

        for ring_vertices in self.intersection_vertices:
            nv = len(ring_vertices)-1
            x = ring_vertices[:, 0]
            y = ring_vertices[:, 1]

            # int(1)
            sint[0] += np.array([ y[j+1] * x[j] - y[j] * x[j+1]  for j in xrange(nv) ], dtype=float).sum()
            # int(x)
            sint[1] += np.array([ ((x[j]+x[j+1])**2 - x[j]*x[j+1]) * (y[j+1]-y[j]) for j in xrange(nv) ],
                                dtype=float).sum()
            # int(y)
            sint[2] += np.array([ ((y[j]+y[j+1])**2 - y[j]*y[j+1]) * (x[j+1]-x[j]) for j in xrange(nv) ],
                                dtype=float).sum()
            # int(xy)
            sint[3] += np.array(
                [(y[j+1]-y[j]) * ( (y[j+1]-y[j]) * (2*x[j]*x[j+1]-x[j]**2) + 2*y[j]*(x[j]**2+x[j+1]**2) +
                                    2*x[j]*x[j+1]*y[j+1]) for j in xrange(nv)],
            dtype=float).sum()
            # int(x**2)
            sint[4] += np.array([ (y[j+1]-y[j]) * (x[j]**3 + x[j]**2*x[j+1] + x[j]*x[j+1]**2 + x[j+1]**3)
                                 for j in xrange(nv)], dtype=float).sum()
            # int(y**2)
            sint[5] += np.array([ (x[j+1]-x[j]) * (y[j]**3 + y[j]**2*y[j+1] + y[j]*y[j+1]**2 + y[j+1]**3)
                                 for j in xrange(nv)], dtype=float).sum()

        sint *= mult

        return sint





def get_hydrostatics(V, F, mass=None, CG=None):
    """Computes the hydrostatics of the mesh and return the clipped mesh.

        Computes the hydrostatics properties of the mesh. Depending on the information given, the equilibrium is
        computed iteratively.
        1) If none of the mass and the center of gravity position are given,
        1) If only the mass of the body is given, the mesh position will be adjusted to comply with the """
    pass
