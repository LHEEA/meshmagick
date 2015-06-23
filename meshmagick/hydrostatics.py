
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
import math
import warnings

import sys

mult_sf = np.array([1/2., 1/6., -1/6., 1/24., 1/12., -1/12.], dtype=float)

class HydrostaticsMesh:
    def __init__(self, V, F, rho_water=1023., g=9.81):
        self.V = V
        self.F = F
        self.rho_water = rho_water # FIXME : rho_water est en doublon dans les fonctions de module
        self.g = g

        # Defining protected attributes
        self._boundary_vertices = None
        self._sfint = np.zeros(6, dtype=float)
        self._sf = 0.
        self._vw = 0.
        self._cw = np.zeros(3, dtype=np.float)

        # Computing once the volume integrals on all the faces of the initial mesh to speedup subsequent operations
        self._surfint = mm._get_surface_integrals(V, F, sum=False)

        # Computing the characteristic dimension in z
        # TODO : voir si on ne peut pas affiner ce critere et le rendre plus general
        height = V[:, 2].max() - V[:, 2].min()
        self._zmax = height * 1e-1
        self._dz = height * 1e-3
        # TODO : Tuner plus precisement ce critere !
        # TODO : ne pas permettre lors des iterations d'appliqier des corrections plus faibles que ce critere

        # Defining the clipping plane Oxy and updating hydrostatics
        self._plane = mm.Plane()
        self.update(np.zeros(3))

        vw0 = self._vw
        self.update([self._dz, 0., 0.])
        vw1 = self._vw
        self._abs_tol_vol = math.fabs(vw1-vw0)

        # Back to the initial state
        self.update(np.zeros(3))




    def update(self, eta, rel=True):

        # Updating the clipping plane position
        if rel:
            self._plane.update(eta)
        else:
            self._plane.set_position(z=eta[0], phi=eta[1], theta=eta[2])

        # Clipping the mesh by the plane
        self._cV, self._cF, clip_infos = mm.clip_by_plane(self.V, self.F, self._plane, infos=True)

        # Testing if the mesh presents intersections and storing the clipped mesh properties
        if len(clip_infos['PolygonsNewID']) == 0:
            raise RuntimeError, 'could not compute any intersection polygon'

        # TODO : mettre les updates dans des methodes
        # Extracting a mesh composed by only the faces that have to be updated
        V_update, F_update = mm.extract_faces(self._cV, self._cF, clip_infos['FToUpdateNewID'])

        # Updating surface integrals for underwater faces of the clipped mesh
        self._update_surfint(clip_infos)

        # Projecting the boundary polygons into the frame of the clipping plane
        self._boundary_vertices = []
        for polygon in clip_infos['PolygonsNewID']:
            self._boundary_vertices.append(self._plane.coord_in_plane(self._cV[polygon]))

        # Computing surface integrals for the floating plane
        self._sfint = self._get_floating_surface_integrals()

        # Area of the flotation surface
        self._sf = self._sfint[0]

        # Computing the immersed volume
        self._vw = self.get_immersed_volume()

        # Computing the center of buoyancy
        self._cw = self.get_buoyancy_center()

        return 1



    def _update_surfint(self, clip_infos):
        """Extraction of volume integrals from the initial mesh to the clipped mesh"""
        # On a besoin ici des informations sur l'extraction du maillage par rapport au maillage initial. Il faut donc
        #  sortir les infos d'extraction, tant au niveau des facettes conservees. Pour les facettes crees ou
        # modifiees, il convient de relancer un calcul d'integrales de volume.

        V_update, F_update = mm.extract_faces(self._cV, self._cF, clip_infos['FToUpdateNewID'])

        # Essai :
        self._c_surfint = mm._get_surface_integrals(V_update, F_update, sum=False).sum(axis=0) + \
                          self._surfint[clip_infos['FkeptOldID']].sum(axis=0)

        return

    def get_hydrostatic_stiffness_matrix(self, cog):

        # Warning, this function works in the frame of the clipped mesh !!

        tol = 1e-8

        z_0p = self._plane.Re0[:, 2]

        # z of the buoyancy center in the frame of the flotation plane
        z_c = np.dot(z_0p, self._cw) - self._plane.c # FIXME : devrait etre fait directement dans Plane

        # z of the center of gravity in the frame of the flotation plane
        z_g = np.dot(z_0p, cog) - self._plane.c

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

        # if (Khs < 0.).any():
        #     # FIXME : A retirer et voir pourquoi on a des valeurs negatives parfois !!!
        #     warnings.warn('Some coefficients of the stiffness matrix are negative, correction', RuntimeWarning)
        #     # raise RuntimeWarning, 'Some coefficients of the stiffness matrix are negative, this should not happen'
        Khs = np.fabs(Khs)

        Khs[Khs < tol] = 0.
        return Khs

    # def get_generalized_position(self):
    #     return self._plane.c

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

    def get_immersed_volume(self):

        r13 = self._plane.Re0[0, 2]
        r23 = self._plane.Re0[1, 2]
        vw = self._c_surfint[2] + self._plane.normal[2] * (r13*self._sfint[1] + r23*self._sfint[2] +
                                                      self._plane.c*self._plane.normal[2]*self._sf)
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
        e = self._plane.c
        e2 = e*e

        cw = np.zeros(3, dtype=np.float)
        cw[0] = self._c_surfint[6] + up * (R11**2*s4 + R21**2*s5 + e2*up**2*self._sf +
                                          2*(R11*R21*s3 + e*up*(R11*s1+R21*s2)))
        cw[1] = self._c_surfint[7] + vp * (R12**2*s4 + R22**2*s5 + e2*vp**2*self._sf +
                                          2*(R12*R22*s3 + e*vp*(R12*s1+R22*s2)))
        cw[2] = self._c_surfint[8] + wp * (R13**2*s4 + R23**2*s5 + e2*wp**2*self._sf +
                                          2*(R13*R23*s3 + e*wp*(R13*s1+R23*s2)))

        cw /= (2*self._vw)
        cw[np.fabs(cw)<tol] = 0.
        return cw

    def get_flotation_center(self):
        tol = 1e-9
        Cf = np.asarray([self._sfint[1], self._sfint[2], 0.], dtype=np.float) / self._sf
        Cf[np.fabs(Cf) < tol] = 0.
        return Cf

    def get_displacement(self):
        # This function should not be used in loops for performance reasons, please inline the code
        return self.rho_water * self._vw

    def get_metacentric_radius(self):
        rhox = self._sfint[5] / self._vw
        rhoy = self._sfint[4] / self._vw
        return np.asarray([rhox, rhoy], dtype=np.float)

    def get_wet_surface(self):
        return mm.get_all_faces_properties(self._cV, self._cF)[0].sum()

    def _get_floating_surface_integrals(self):

        sint = np.zeros(6, dtype=float)

        for ring_vertices in self._boundary_vertices:
            nv = len(ring_vertices)-1

            iter = xrange(nv)

            x = ring_vertices[:, 0]
            y = ring_vertices[:, 1]

            # Precomputing some patterns for every vertices
            xjj_xj = np.array([ x[j+1]-x[j] for j in iter], dtype=np.float)
            yjj_yj = np.array([ y[j+1]-y[j] for j in iter], dtype=np.float)
            xjpxjj = np.array([ x[j]+x[j+1] for j in iter], dtype=np.float)
            yjpyjj = np.array([ y[j]+y[j+1] for j in iter], dtype=np.float)
            xjxjj = np.array([ x[j]*x[j+1] for j in iter], dtype=np.float)
            yjyjj = np.array([ y[j]*y[j+1] for j in iter], dtype=np.float)
            xj2 = np.append(np.array([ x[j]*x[j] for j in iter], dtype=np.float), x[0]*x[0])
            yj2 = np.append(np.array([ y[j]*y[j] for j in iter], dtype=np.float), y[0]*y[0])


            # int(1)
            sint[0] += np.array([ xjpxjj[j] * yjj_yj[j] for j in iter ], dtype=np.float).sum()

            # int(x)
            sint[1] += np.array([ (xj2[j] + xjxjj[j] + xj2[j+1])*yjj_yj[j] for j in iter], dtype=np.float).sum()

            # int(y)
            sint[2] += np.array([ (yj2[j] + yjyjj[j] + yj2[j+1])*xjj_xj[j] for j in iter], dtype=np.float).sum()

            # int(xy)
            sint[3] += np.array([ (xj2[j]*(2*y[j]+yjpyjj[j])
                                + xj2[j+1]*(2*y[j+1]+yjpyjj[j])
                                + 2*xjxjj[j]*yjpyjj[j]) * yjj_yj[j] for j in iter], dtype=np.float).sum()

            # int(x**2)
            sint[4] += np.array([ (xj2[j]+xj2[j+1]) * xjpxjj[j] * yjj_yj[j] for j in iter], dtype=np.float).sum()

            # int(y**2)
            sint[5] += np.array([ (yj2[j]+yj2[j+1]) * yjpyjj[j] * xjj_xj[j] for j in iter], dtype=np.float).sum()

        sint *= mult_sf

        return sint

# ======================================================================================================================

def print_hysdrostatics_report(hs_data):

    # TODO : Ajouter les metacentres --> infos sur la stab
    # TODO : mettre un ordre d'apparition fixe !

    hs_text = {
        'disp' : 'Displacement (m**3):\n\t%E\n',
        'Cw'   : 'Center of Buoyancy (m):\n\t%E, %E, %E\n',
        'Sf'   : 'Waterplane area (m**2):\n\t%E\n',
        'mass' : 'Mass (kg):\n\t%E\n',
        'res'  : 'Residual (kg, Nm, Nm):\n\t%E, %E, %E\n',
        'cog'  : 'Gravity center (m):\n\t%E, %E, %E\n',
        'K33'  : 'Heave stiffness (N/m):\n\t%E\n',
        'Khs'  : 'Hydrostatic Stiffness matrix:\n'
                 '\t%E, %E, %E\n'
                 '\t%E, %E, %E\n'
                 '\t%E, %E, %E\n',
        'draft': 'Draft (m):\n\t%E\n',
        'Ws'   : 'Wetted surface:\n\t%E\n',
        'Cf'   : 'Center of flotation (m):\n\t%E, %E, %E\n',
        'meta_radius' : 'Initial metacentric radius (m):\n\tx: %f\n\ty: %f\n'
    }

    print '\nHydrostatic Report'
    print '------------------\n'
    for key in hs_text:
        if hs_data.has_key(key):
            repl = hs_data[key]
            if isinstance(repl, np.ndarray):
                repl = tuple(repl.flatten())
            print hs_text[key] % repl

    return 1



def _get_residual(rho_water, g, vw, cw, mass, cog):

    rgvw = rho_water*g*vw
    mg = mass*g

    res = np.array([
         rgvw - mg,
         rgvw * cw[1] - mg * cog[1],
        -rgvw * cw[0] + mg * cog[0]
    ], dtype=np.float)
    return res

def get_hydrostatics(hsMesh, mass=None, cog=None, zcog=None, rho_water=1023, g=9.81, anim=False, verbose=False):
    """Computes the hydrostatics of the mesh and return the clipped mesh.

        Computes the hydrostatics properties of the mesh. Depending on the information given, the equilibrium is
        computed iteratively.
        1) If none of the mass and the center of gravity position are given,
        1) If only the mass of the body is given, the mesh position will be adjusted to comply with the """

    # TODO : recuperer le deplacement total pour verifier que la masse fournie est consistante

    # TODO : decouper les deux cas principaux en deux fonctions

    # TODO : mettre ces procedures en methodes de la classe HydrostaticMesh
    # Instantiation of the hydrostatic mesh object
    # hsMesh = HydrostaticsMesh(V, F, rho_water=rho_water, g=g)

    hs_data = dict()

    #==================================================================
    if mass is None: # No equilibrium is performed if mass is not given
    #==================================================================
        if verbose:
            print 'Computation of hydrostatics with the given mesh position'

        disp = hsMesh._vw           # displacement
        Cw = hsMesh._cw             # Center of buoyancy
        Sf = hsMesh._sf             # Area of the flotation plane
        mass = rho_water * disp     # Mass of the device
        cV = hsMesh._cV             # Vertices of the mesh
        cF = hsMesh._cF             # Faces of the mesh

        # Choosing wether we return a stiffness in heave only or a stiffness matrix
        if cog is None:
            if zcog is None:
                # Return only the stiffness in heave
                hs_data['K33'] = rho_water*g*Sf
            else:
                # Computing the stiffness matrix
                hs_data['Khs'] = hsMesh.get_hydrostatic_stiffness_matrix(np.array([0., 0., zcog], dtype=np.float))
                hs_data['cog'] = Cw.copy()
                hs_data['cog'][2] = zcog
        else:
            # Computing the stiffness matrix with the cog given
            hs_data['Khs'] = hsMesh.get_hydrostatic_stiffness_matrix(cog)
            hs_data['res'] = _get_residual(rho_water, g, disp, Cw, mass, cog)
            hs_data['cog'] = cog

        hs_data['disp'] = hsMesh._vw
        hs_data['Cw'] = hsMesh._cw
        hs_data['Sf'] = hsMesh._sf
        hs_data['mass'] = mass
        hs_data['draft'] = cV[:,2].min()
        hs_data['Ws'] = hsMesh.get_wet_surface()

        hs_data['Cf'] = hsMesh.get_flotation_center()

        hs_data['meta_radius'] = hsMesh.get_metacentric_radius()

    #========================================================================
    else: # mass is given explicitly, iterative resolution of the equilibrium
    #========================================================================

        maxiter = 100
        rg = rho_water * g
        mg = mass * g
        niter = 0

        zmax = hsMesh._zmax
        abs_tol_pos = rho_water * hsMesh._abs_tol_vol

        if anim:
            # Removing all files eq*.vtu
            import os, glob
            for eqx in glob.glob('eq*.vtu'):
                os.remove(eqx)

            filename = 'eq0.vtu'
            mm.write_VTU(filename, hsMesh._cV, hsMesh._cF)

        #-----------------------------------------------------
        if cog is None: # Equilibrium resolution in heave only
        #-----------------------------------------------------
            if verbose:
                print "Equilibrium resolution knowing only mass --> iterations on z only"

            # FIXME : il faut que la mise en equilibre se fasse suivant la normale du plan !!!

            res = 0.
            while 1:
                # Iteration loop

                if niter == maxiter:
                    status = 0
                    break

                res_old = res
                res = rho_water * hsMesh._vw - mass # residual

                if verbose:
                    print 'Iteration %u:' % niter
                    print '\t-> Residual = %E (kg)' % res
                    print '\t-> Target mass: %E (kg); Current: %E (kg)' % (mass, rho_water*hsMesh._vw)

                # Convergence criteria
                if math.fabs(res) < abs_tol_pos:
                    status = 1
                    break

                niter += 1
                stiffness = rg * hsMesh._sf # K33 : stiffness in heave
                dz = g * res/stiffness

                # Checking for a sign modification in the residual
                # TODO : play on the sign of the correction instead of the sign of the residual...
                # TODO : harmoniser ce critere avec celui utilise dans la version 6dof
                if res*res_old < 0.:
                    if math.fabs(res) > math.fabs(res_old):
                        reduc = 1/4.
                    else:
                        reduc = 1/2.
                    zmax *= reduc

                if math.fabs(dz) > zmax:
                    dz = math.copysign(zmax, dz)

                zcur = hsMesh._plane.c
                hsMesh.update([-dz, 0., 0.]) # The - sign is here to make the plane move, not the mesh

                if anim:
                    filename = 'eq%u.vtu'%niter
                    mm.write_VTU(filename, hsMesh._cV, hsMesh._cF)

            hs_data['res'] = np.array([res, 0., 0.], dtype=np.float)

            # Moving the mesh
            # FIXME : il faut que la mise en equilibre se fasse suivant la normale du plan !!!
            cV = hsMesh._plane.coord_in_plane(hsMesh._cV)
            cF = hsMesh._cF

            hs_data['disp'] = hsMesh._vw
            hs_data['Cw'] = hsMesh._cw # FIXME : Cw doit etre fourni dans le nouveau repere
            hs_data['Sf'] = hsMesh._sf
            hs_data['mass'] = rho_water * hsMesh._vw
            hs_data['draft'] = cV[:, 2].min()

            if verbose:
                if status == 1:
                    print "\nEquilibrium found in %u iterations" % niter
                else:
                    print "\nEquilibrium approached but the mesh is not refined enough to reach convergence"
                print '\nZ translation on the initial mesh : %f (m)' % (-zcur)

            if zcog is None:
                hs_data['K33'] = rg*hsMesh._sf
                # hs_data['cog'] = np.array([0., 0., zcog])
            else:
                hs_data['Khs'] = hsMesh.get_hydrostatic_stiffness_matrix(np.array([0., 0., zcog], dtype=np.float))
                hs_data['cog'] = hsMesh._cw.copy()
                hs_data['cog'][2] = zcog # FIXME : cog doit etre fourni dans le nouveau repere

        # ---------------------------------------------------------
        else: # cog has been specified, 6dof equilibrium resolution
        # ---------------------------------------------------------
            if verbose:
                print "Equilibrium resolution knowing mass and center of gravity --> iterations on 6dof"

            deta = np.zeros(3, dtype=np.float)

            while 1:
                # Iteration loop

                if niter == maxiter:
                    status = 0
                    break

                # Projecting cog and buoyancy center in the clipping plane frame
                # TODO : ajouter une methode dans hsMesh pour ne pas appeler directement les methodes de Plane...
                cog_e, cw_e = hsMesh._plane.coord_in_plane(np.array([cog, hsMesh._cw]))

                # TODO : faire que _get_residual n'ait pas autant de parametres d'entree
                res = _get_residual(rho_water, g, hsMesh._vw, cw_e, mass, cog_e) # expressed in the plane frame

                if verbose:
                    print 'Iteration %u: ' % niter
                    print '\t-> Residual (N, Nm, Nm) = %E, %E, %E' % tuple(res.flatten())
                    print '\t-> Target mass: %E (kg); Current: %E (kg)' % (mass, rho_water*hsMesh._vw)

                    print '\t-> Relative x position of cog and B : %E (m)' % (math.fabs(cog_e[0]-cw_e[0]))
                    print '\t-> Relative y position of cog and B : %E (m)' % (math.fabs(cog_e[1]-cw_e[1]))

                # Convergence criteria
                if (np.fabs(res) < np.ones(3)*abs_tol_pos).all(): # TODO : travailler sur ce critere (differencier pos et rot)
                    status=1
                    break

                niter += 1
                Khs = hsMesh.get_hydrostatic_stiffness_matrix(cog)

                deta_prev = deta
                deta = np.linalg.solve(Khs, res)

                deta_sign = deta * deta_prev

                # TODO : vectoriser la partie suivante
                if deta_sign[0] < 0.: # Change sign
                    zmax = min([math.fabs(deta[0]-deta_prev[0])/2., zmax/2.])

                if math.fabs(deta[0]) > zmax:
                    deta[0] = math.copysign(zmax, deta[0])

                # Updating the plane position
                hsMesh.update(-deta) # The - sign makes the plane move, not the mesh...

                if anim:
                    filename = 'eq%u.vtu'%niter
                    mm.write_VTU(filename, hsMesh._cV, hsMesh._cF)

            hs_data['res'] = res

            # Moving the mesh
            cV = hsMesh._cV
            cF = hsMesh._cF

            cV = hsMesh._plane.coord_in_plane(cV)

            hs_data['disp'] = hsMesh._vw
            hs_data['Cw'] = hsMesh._plane.coord_in_plane(cw_e)
            hs_data['Sf'] = hsMesh._sf
            hs_data['mass'] = rho_water * hsMesh._vw
            hs_data['draft'] = cV[:, 2].min()
            hs_data['res'] = res
            hs_data['Khs'] = Khs
            hs_data['cog'] = cog_e


            if verbose:
                if status == 1:
                    print "\nEquilibrium found in %u iterations" % niter
                else:
                    print "\nEquilibrium approached but the mesh is not refined enough to reach convergence"
                z, phi, theta = hsMesh._plane.get_position()
                print '\nTransformation of the initial mesh '
                print 'z (m)       : %f' % (-z)
                print 'phi (deg)   : %f' % (-phi*180./math.pi)
                print 'theta (deg) : %f' % (-theta*180./math.pi)

    if verbose:
        print_hysdrostatics_report(hs_data)

    # TODO : renvoyer egalement les infos hydrostatiques sous forme de dictionnaire -> ou alors sortir un fichier !
    return cV, cF, hs_data


def get_GZ_curves(hsMesh, zcog, spacing=2., rho_water=1023, g=9.81, verbose=False):
    # TODO : mettre rho_water et g en variable de module avec les valeurs par defaut !
    # TODO : verifier que hsMesh est bien du type HydrostaticMesh

    if not isinstance(hsMesh, HydrostaticsMesh):
        raise RuntimeError, 'hsMesh argument must be an instance of HydrostaticMesh class'
    # Computing hydrostatics for the initial mesh
    # hsMesh = HydrostaticsMesh(V, F, rho_water=rho_water, g=g)
    cV, cF, hs_data = get_hydrostatics(hsMesh, zcog=zcog)

    Cw0 = hs_data['Cw']
    Vw0 = hs_data['disp']
    mass = rho_water*Vw0

    G = np.array([Cw0[0], Cw0[1], zcog], dtype=np.float)

    angles = np.arange(0., 180.+spacing, spacing) * math.pi/180.
    # TODO : voir pour parametrer les bornes

    # Computing the GZ curve in phi
    for (index, phi) in enumerate(angles):
        print 'phi: %f (deg)' % (phi*180/math.pi)
        # Setting the plane
        # hsMesh._plane.set_position(phi=phi)
        eta = np.array([0., phi, 0.], dtype=np.float)
        hsMesh.update(eta, rel=False)
        cV, cF, hs_data = get_hydrostatics(hsMesh, mass=mass, verbose=False)
        # mm.show(cV, cF)
        Cwj = hs_data['Cw']

        filename = 'eq%u.vtu'%index
        mm.write_VTU(filename, cV, cF)




    def get_area_curve(V, F, dir):

        raise NotImplementedError
