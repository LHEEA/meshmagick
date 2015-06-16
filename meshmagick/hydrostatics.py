
__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-2015, Ecole Centrale de Nantes"
__credits__ = "Francois Rongere"
__licence__ = "CeCILL"
__version__ = "0.3.1"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@ec-nantes.fr"
__status__ = "Development"

import meshmagick as mm
import numpy as np

mult_sf = np.array([1/2., 1/6., -1/6., 1/12., 1/12., -1/12.], dtype=float)

class HydrostaticsMesh:
    def __init__(self, V, F, rho_water=1023., g=9.81):
        self.V = V
        self.F = F
        self.rho_water = rho_water
        self.g = g

        # Defining protected attributes
        self._cV = None
        self._cF = None
        self._c_areas = None
        self._c_normals = None
        self._c_centers = None
        self._boundary_vertices = None
        self._sfint = np.zeros(6, dtype=float)
        self._sf = 0.
        self._vw = 0.
        self._cw = np.zeros(3, dtype=np.float)

        # Computing initial mesh properties
        self.areas, self.normals, self.centers = mm.get_all_faces_properties(V, F) # FIXME : pas utile a priori...

        # Computing once the volume integrals on faces of the initial mesh to be used by the clipped mesh
        self._surfint = mm._get_surface_integrals(V, F, sum=False)

        # Defining the clipping plane Oxy and updating hydrostatics
        self._plane = mm.Plane()
        self.update([0., 0., 0.])

        self._has_plane_changed = True # FIXME : utile ?


    def update(self, eta):

        # Updating the clipping plane position
        self._plane.set_position(z=eta[0], phi=eta[1], theta=eta[2])
        # TODO : ajouter une fonction permettant de recuperer la position du plan

        # Clipping the mesh by the plane
        self._cV, self._cF, clip_infos = mm.clip_by_plane(self.V, self.F, self._plane, infos=True)

        # Testing if the mesh presents intersections and storing the clipped mesh properties
        if len(clip_infos['PolygonsNewID']) == 0:
            raise RuntimeError, 'could not compute any intersection polygon'

        # TODO : mettre les updates dans des methodes
        # Extracting a mesh composed by only the faces that have to be updated
        V_update, F_update = mm.extract_faces(self._cV, self._cF, clip_infos['FToUpdateNewID'])

        # Updating faces properties for the clipped mesh
        self._update_faces_properties(V_update, F_update, clip_infos) # FIXME : useless !!!

        # Updating surface integrals for underwater faces of the clipped mesh
        self._update_surfint(V_update, F_update, clip_infos)

        # Projecting the boundary polygons into the frame of the clipping plane
        self._boundary_vertices = []
        for polygon in clip_infos['PolygonsNewID']:
            self._boundary_vertices.append(self._plane.coord_in_plane(self._cV[polygon]))

        # Computing surface integrals for the floating plane
        self._sfint = self._get_floating_surface_integrals()

        # Area of the flotation surface
        self._sf = self._sfint[0]

        # Computing the immersed volume
        self._vw = self.get_vw()

        # Computing the center of buoyancy
        self._cw = self.get_buoyancy_center()

        return 1



    def _update_surfint(self, V_update, F_update, clip_infos):
        """Extraction of volume integrals from the initial mesh to the clipped mesh"""
        # On a besoin ici des informations sur l'extraction du maillage par rapport au maillage initial. Il faut donc
        #  sortir les infos d'extraction, tant au niveau des facettes conservees. Pour les facettes crees ou
        # modifiees, il convient de relancer un calcul d'integrales de volume.
        up_surfint = mm._get_surface_integrals(V_update, F_update)

        c_surfint = np.zeros((self._cF.shape[0], 12), dtype=np.float)
        c_surfint[clip_infos['FkeptNewID']] = self._surfint[clip_infos['FkeptOldID']]
        c_surfint[clip_infos['FToUpdateNewID']] = up_surfint

        self._c_surfint = c_surfint.sum(axis=0)

        return

    def get_hydrostatic_stiffness_matrix(self, cog):

        tol = 1e-9

        z_0p = self._plane.Re0[:, 2]

        # z of the buoyancy center in the frame of the flotation plane
        z_c = np.dot(z_0p, self._cw) - self._plane.e # FIXME : devrait etre fait directement dans Plane

        # z of the center of gravity in the frame of the flotation plane
        z_g = np.dot(z_0p, cog) - self._plane.e

        corr = self._vw * (z_c-z_g)

        k33 = self._sfint[0]
        k34 = self._sfint[2]
        k35 = -self._sfint[1]
        k44 = self._sfint[5] + corr
        k45 = -self._sfint[3]
        k55 = self._sfint[4] + corr

        Khs = self.rho_water * self.g * \
            np.array([
                [k33, k34, k35],
                [k34, k44, k45],
                [k35, k45, k55]
            ], dtype=np.float)
        Khs[Khs<tol] = 0.
        return Khs

    def _update_faces_properties(self, V_update, F_update, clip_infos):

        up_areas, up_normals, up_centers = mm.get_all_faces_properties(V_update, F_update)

        # Collectively updating properties of wetted mesh
        nf = self._cF.shape[0]
        self._c_areas = np.zeros(nf, dtype=float)
        self._c_areas[clip_infos['FkeptNewID']] = self.areas[clip_infos['FkeptOldID']]
        self._c_areas[clip_infos['FToUpdateNewID']] = up_areas

        self._c_normals = np.zeros((nf, 3), dtype=float)
        self._c_normals[clip_infos['FkeptNewID']] = self.normals[clip_infos['FkeptOldID']]
        self._c_normals[clip_infos['FToUpdateNewID']] = up_normals

        self._c_centers = np.zeros((nf, 3), dtype=float)
        self._c_centers[clip_infos['FkeptNewID']] = self.centers[clip_infos['FkeptOldID']]
        self._c_centers[clip_infos['FToUpdateNewID']] = up_centers

        return

    def get_vw(self):

        r13 = self._plane.Re0[0, 2]
        r23 = self._plane.Re0[1, 2]
        vw = self._c_surfint[2] + self._plane.normal[2] * (r13*self._sfint[1] + r23*self._sfint[2] +
                                                      self._plane.e*self._plane.normal[2]*self._sf)
        return vw

    def get_buoyancy_center(self):

        tol = 1e-9

        R11 = self._plane.Re0[0, 0]
        R21 = self._plane.Re0[1, 0]
        R12 = self._plane.Re0[0, 1]
        R22 = self._plane.Re0[1, 1]
        R13 = self._plane.Re0[0, 2]
        R23 = self._plane.Re0[1, 2]

        s1 = self._sfint[1]
        s2 = self._sfint[2]
        s3 = self._sfint[3]
        s4 = self._sfint[4]
        s5 = self._sfint[5]

        (up, vp, wp) = self._plane.normal
        e = self._plane.e
        e2 = e*e

        cw = np.zeros(3, dtype=np.float)
        cw[0] = self._c_surfint[3] + up * (R11**2*s4 + R21**2*s5 + e2*up**2*self._sf +
                                          2*(R11*R21*s3 + e*up*(R11*s1+R21*s2)))
        cw[1] = self._c_surfint[4] + vp * (R12**2*s4 + R22**2*s5 + e2*vp**2*self._sf +
                                          2*(R12*R22*s3 + e*vp*(R12*s1+R22*s2)))
        cw[2] = self._c_surfint[5] + wp * (R13**2*s4 + R23**2*s5 + e2*wp**2*self._sf +
                                          2*(R13*R23*s3 + e*wp*(R13*s1+R23*s2)))

        cw /= (2*self._vw)
        cw[np.fabs(cw)<tol] = 0.
        return cw

    def get_displacement(self):
        # This function should not be used in loops for performance reasons, please inline the code
        return self.rho_water * self._vw

    def get_wet_surface(self):
        return np.sum(self._c_areas)

    def _get_floating_surface_integrals(self):


        sint = np.zeros(6, dtype=float)

        for ring_vertices in self._boundary_vertices:
            nv = len(ring_vertices)-1
            x = ring_vertices[:, 0]
            y = ring_vertices[:, 1]

            # int(1)
            sint[0] += np.array([ y[j+1] * x[j] - y[j] * x[j+1]  for j in xrange(nv) ], dtype=float).sum()
            # int(x)
            sint[1] += np.array([ ((x[j]+x[j+1])**2 - x[j]*x[j+1]) * (y[j+1]-y[j]) for j in xrange(nv) ],dtype=float).sum()
            # int(y)
            sint[2] += np.array([ ((y[j]+y[j+1])**2 - y[j]*y[j+1]) * (x[j+1]-x[j]) for j in xrange(nv) ],dtype=float).sum()
            # int(xy)
            sint[3] += np.array(
                [(y[j+1]-y[j]) * ( (y[j+1]-y[j]) * (2*x[j]*x[j+1]-x[j]**2) + 2*y[j]*(x[j]**2+x[j+1]**2) +
                                    2*x[j]*x[j+1]*y[j+1] ) for j in xrange(nv)],dtype=float).sum()
            # int(x**2)
            sint[4] += np.array([(y[j+1]-y[j]) * (x[j]**3 + x[j]**2*x[j+1] + x[j]*x[j+1]**2 + x[j+1]**3)
                                 for j in xrange(nv)], dtype=float).sum()
            # int(y**2)
            sint[5] += np.array([(x[j+1]-x[j]) * (y[j]**3 + y[j]**2*y[j+1] + y[j]*y[j+1]**2 + y[j+1]**3)
                                 for j in xrange(nv)], dtype=float).sum()

        sint *= mult_sf

        return sint



def get_hydrostatics(V, F, mass=None, cog=None, zcog=None, rho_water=1023, g=9.81):
    """Computes the hydrostatics of the mesh and return the clipped mesh.

        Computes the hydrostatics properties of the mesh. Depending on the information given, the equilibrium is
        computed iteratively.
        1) If none of the mass and the center of gravity position are given,
        1) If only the mass of the body is given, the mesh position will be adjusted to comply with the """

    hsMesh = HydrostaticsMesh(V, F, rho_water=rho_water, g=g)

    # cog = mm.get_inertial_properties(V, F)[1]
    cog = np.zeros(3, dtype=np.float)

    Khs = hsMesh.get_hydrostatic_stiffness_matrix(cog)



    return hsMesh._cV, hsMesh._cF
