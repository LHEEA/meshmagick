
__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-2015, Ecole Centrale de Nantes"
__credits__ = "Francois Rongere"
__licence__ = "CeCILL"
__version__ = "1.0"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@ec-nantes.fr"
__status__ = "Development"

import meshmagick as mm
import numpy as np
import math

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

        # Defining the clipping plane Oxy and updating hydrostatics
        self._plane = mm.Plane()
        self.update(np.zeros(3))


    def _init_tolerances(self):

        # TODO : voir si on ne peut pas affiner ce critere et le rendre plus general
        height = self.V[:, 2].max() - self.V[:, 2].min()
        self._zmax = height * 1e-1
        self._dz = height * 1e-3
        # TODO : Tuner plus precisement ce critere !
        # TODO : ne pas permettre lors des iterations d'appliqier des corrections plus faibles que ce critere

        vw0 = self._vw
        self.update([self._dz, 0., 0.])
        vw1 = self._vw
        self._abs_tol_vol = math.fabs(vw1-vw0)

        # Back to the initial state
        self.update(np.zeros(3))


    def update(self, eta, rel=True):

        # TODO : Make mm.clip_by_plane function to return a code to manage void intersections

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
        # V_update, F_update = mm.extract_faces(self._cV, self._cF, clip_infos['FToUpdateNewID'])

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

    print '\n------------------'
    print 'Hydrostatic Report'
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
            print '\n--------------------------------------------------------'
            print 'Computation of hydrostatics with the given mesh position'
            print '--------------------------------------------------------'

        disp = hsMesh._vw             # displacement
        Cw   = hsMesh._cw             # Center of buoyancy
        Sf   = hsMesh._sf             # Area of the flotation plane
        mass = rho_water * disp       # Mass of the device
        cV   = hsMesh._cV             # Vertices of the mesh
        cF   = hsMesh._cF             # Faces of the mesh

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

        hs_data['disp']  = hsMesh._vw
        hs_data['Cw']    = hsMesh._cw
        hs_data['Sf']    = hsMesh._sf
        hs_data['mass']  = mass
        hs_data['draft'] = cV[:,2].min()
        hs_data['Ws']    = hsMesh.get_wet_surface()

        hs_data['Cf'] = hsMesh.get_flotation_center()

        hs_data['meta_radius'] = hsMesh.get_metacentric_radius()

    #========================================================================
    else: # mass is given explicitly, iterative resolution of the equilibrium
    #========================================================================

        # Initialization of tolerances
        hsMesh._init_tolerances()

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
                print '\n----------------------------------------------------'
                print 'Hydrostatic equilibrium resolution knowing only mass'
                print '----------------------------------------------------'

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
                    print '\t-> Target mass: %E (kg); Current: %E (kg)\n' % (mass, rho_water*hsMesh._vw)

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
                print '\n---------------------------------------------------------'
                print "Equilibrium resolution knowing mass and center of gravity"
                print '---------------------------------------------------------'

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
    hs_data = get_hydrostatics(hsMesh, zcog=zcog)[2]

    Cw0 = hs_data['Cw']
    Vw0 = hs_data['disp']
    mass = hs_data['mass']
    cog = hs_data['cog']
    a = cog[2] - Cw0[2]

    angles = np.arange(spacing, 180.+spacing, spacing)
    # TODO : voir pour parametrer les bornes
    # TODO : permettre de balayer egalement les angles negatifs

    # TODO : n'accepter qu'une direction a la fois...

    GZ_phi = np.zeros(angles.shape[0]+1, dtype=np.float)

    # Computing the GZ curve in phi
    for (index, phi) in enumerate(angles * math.pi/180.):
        if verbose:
            print 'phi: %f (deg)' % (phi*180/math.pi)
        # Setting the plane
        eta = np.array([0., phi, 0.], dtype=np.float)
        hsMesh.update(eta, rel=False)
        # Ensuring iso-displacement
        hs_data = get_hydrostatics(hsMesh, mass=mass, verbose=False)[2]
        if verbose:
            print 'Target disp : %E ; current disp : %E' % (Vw0, hs_data['disp'])
        Cwj = hs_data['Cw']

        # Computing the transverse metacentric point relative to angle phi
        t = (cog[1] - Cwj[1]) / hsMesh._plane.normal[1]
        hz = Cwj[2] + t*hsMesh._plane.normal[2]

        # Metacentric height:
        h =  hz-Cwj[2]

        GZ_phi[index+1] = (h-a)*math.sin(phi)


    GZ_theta = np.zeros(angles.shape[0]+1, dtype=np.float)
    for (index, theta) in enumerate(angles * math.pi/180.):
        if verbose:
            print 'theta: %f (deg)' % (angles[index])
        # Setting the plane
        eta = np.array([0., 0., theta], dtype=np.float)
        hsMesh.update(eta, rel=False)
        # Ensuring iso-displacement
        cV, cF, hs_data = get_hydrostatics(hsMesh, mass=mass, verbose=False)
        if verbose:
            print 'Target disp : %E ; obtained disp : %E' % (Vw0, hs_data['disp'])
        filename = 'eq%u.vtu'%(index+1)
        mm.write_VTU(filename, cV, cF)

        Cwj = hs_data['Cw']

        # Computing the transverse metacentric point relative to angle phi
        t = (cog[0] - Cwj[0]) / hsMesh._plane.normal[0]
        hz = Cwj[2] + t*hsMesh._plane.normal[2]

        # Metacentric height:
        h =  hz-Cwj[2]

        GZ_theta[index+1] = (h-a)*math.sin(theta)

    try:
        import matplotlib.pyplot as plt
        plt.figure(1)
        plt.subplot(211)
        plt.plot(np.arange(0., 180.+spacing, spacing), GZ_phi)
        plt.xlabel('Roll angle (deg)')
        plt.ylabel('GZ (m)')
        plt.grid()

        plt.subplot(212)
        plt.plot(np.arange(0., 180.+spacing, spacing), GZ_theta)
        plt.xlabel('Pitch angle (deg)')
        plt.ylabel('GZ (m)')
        plt.grid()
        plt.show()
    except:
        print 'No visualization of GZ curves from meshmagick as matplotlib is not available'


    return 1

    def get_area_curve(V, F, dir):

        raise NotImplementedError





#### END OLD CODE






def show(V, F, CG=None, B=None):

    # TODO: afficher polygones...
    import MMviewer
    import vtk
    my_viewer = MMviewer.MMViewer()

    # Adding mesh
    polydata_mesh = mm._build_vtkPolyData(V, F)
    my_viewer.add_polydata(polydata_mesh)

    # Flotation plane
    plane = vtk.vtkPlaneSource()

    # plane.SetNormal(0, 0, 1)
    xmin, ymin = V[:, :2].min(axis=0)
    xmax, ymax = V[:, :2].max(axis=0)
    plane.SetOrigin(xmin, ymax, 0)
    plane.SetPoint1(xmin, ymin, 0)
    plane.SetPoint2(xmax, ymax, 0)

    my_viewer.add_polydata(plane.GetOutput(), color=[0.1, 0.9, 0.7])

    # Adding O
    pO = vtk.vtkPoints()
    vO = vtk.vtkCellArray()

    iO = pO.InsertNextPoint([0, 0, 0])
    vO.InsertNextCell(1)
    vO.InsertCellPoint(iO)

    pdO = vtk.vtkPolyData()
    pdO.SetPoints(pO)
    pdO.SetVerts(vO)

    my_viewer.add_polydata(pdO, color=[0, 0, 0])

    # Adding CG
    if CG is not None:
        pCG = vtk.vtkPoints()
        vCG = vtk.vtkCellArray()

        iCG = pCG.InsertNextPoint(CG)
        vCG.InsertNextCell(1)
        vCG.InsertCellPoint(iCG)

        pdCG = vtk.vtkPolyData()
        pdCG.SetPoints(pCG)
        pdCG.SetVerts(vCG)

        my_viewer.add_polydata(pdCG, color=[1, 0, 0])

        # Ploting also a line between O and CG
        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(CG)

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, 0)
        line.GetPointIds().SetId(1, 1)

        lines = vtk.vtkCellArray()
        lines.InsertNextCell(line)

        lines_pd = vtk.vtkPolyData()
        lines_pd.SetPoints(points)
        lines_pd.SetLines(lines)

        my_viewer.add_polydata(lines_pd, color=[0, 0, 0])

    # Adding B
    if B is not None:
        pB = vtk.vtkPoints()
        vB = vtk.vtkCellArray()

        iB = pB.InsertNextPoint(B)
        vB.InsertNextCell(1)
        vB.InsertCellPoint(iB)

        pdB = vtk.vtkPolyData()
        pdB.SetPoints(pB)
        pdB.SetVerts(vB)

        my_viewer.add_polydata(pdB, color=[0, 1, 0])

        # Ploting also a line between O and B
        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(B)

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, 0)
        line.GetPointIds().SetId(1, 1)

        lines = vtk.vtkCellArray()
        lines.InsertNextCell(line)

        lines_pd = vtk.vtkPolyData()
        lines_pd.SetPoints(points)
        lines_pd.SetLines(lines)

        my_viewer.add_polydata(lines_pd, color=[0, 0, 0])




    # # Display intersection polygons
    # if polygons is not None:
    #     for polygon_ids in polygons:
    #         np = polygon_ids.size
    #         points_coord = V[polygon_ids]
    #         points = vtk.vtkPoints()
    #         for point in points_coord:
    #             points.InsertNextPoint(point)
    #         polygon = vtk.vtkPolygon()
    #         polygon.GetPointIds().SetNumberOfIds(np)
    #         for i in xrange(np):
    #             polygon.GetPointIds().SetId(i, i)
    #         cell = vtk.vtkCellArray()
    #         cell.InsertNextCell(polygon)
    #
    #         polygon_pd = vtk.vtkPolyData()
    #         polygon_pd.SetPoints(points)
    #         polygon_pd.SetPolys(cell)
    #
    #         my_viewer.add_polydata(polygon_pd, color=[1, 0, 0])


    # Adding corner annotation
    ca = vtk.vtkCornerAnnotation()
    ca.SetLinearFontScaleFactor(2)
    ca.SetNonlinearFontScaleFactor(1)
    ca.SetMaximumFontSize(20)
    labels = "O: Black; CG: Red; B: Green\nTo see points, press 'w'"
    ca.SetText(2, labels)
    ca.GetTextProperty().SetColor(0., 0., 0.)
    my_viewer.renderer.AddViewProp(ca)


    my_viewer.show()
    my_viewer.finalize()





def _get_rotation_matrix(thetax, thetay):

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
    R = ctheta*np.eye(3) \
       + (1-ctheta) * np.array([[nx*nx, nxny, 0.],
                                [nxny, ny*ny, 0.],
                                [0., 0., 0.]]) \
       + stheta * np.array([[0., 0.,  ny],
                            [0., 0., -nx],
                            [-ny, nx, 0.]])
    return R



def compute_equilibrium(V, F, disp, CG, rho_water=1023, grav=9.81,
                        reltol=1e-2, itermax=100, max_nb_relaunch=10, theta_relax=2, z_relax=0.1,
                        verbose=False, anim=False):
    """

    Parameters
    ----------
    V : ndarray
        Array of mesh vertices
    F :
        Array of mesh connectivities
    disp : float
        Mass displacement of the floater (in tons)
    CG : ndarray
        Position of the center of gravity
    rho_water : float, optional
        Density of water (Default 1023 kg/m**3)
    grav : float, optional
        Acceleration of gravity (default 9.81 m/s**2)
    abstol : float, optional
        Absolute tolerance
    itermax : int, optional
        Maximum number of iterations
    verbose : bool, optional
        If True, prints results on screen. Default is False.

    Returns
    -------

    """
    max_nb_relaunch = 10

    # Relaxation parameters
    z_relax = 0.1 # in meters
    theta_relax_x = theta_relax * math.pi/180.
    theta_relax_y = theta_relax * math.pi/180.


    rhog = rho_water*grav
    mg = disp*grav*1e3 # Conversion of displacement in kg

    dz = 0.
    thetax = 0.
    thetay = 0.
    iter = 0

    R = np.eye(3, 3)
    dz_update = 0.

    Vc = V.copy()
    Fc = F.copy()
    CGc = CG.copy()

    # origin = np.zeros(3)

    nb_relaunch = 0

    while True:
        print '\n', iter
        if iter == (nb_relaunch+1)*itermax:
            # Max iterations reach

            if nb_relaunch < max_nb_relaunch:
                nb_relaunch += 1
                # Random on the position of the body
                print 'Max iteration reached: relaunching number %u with random orientation' % nb_relaunch
                thetax, thetay = np.random.rand(2)*2*math.pi
                R = _get_rotation_matrix(thetax, thetay)
                Vc = np.transpose(np.dot(R, Vc.T))
                CGc = np.dot(R, CGc)
                # origin = np.dot(R, origin)
                dz = thetax = thetay = 0.
            else:
                code = 0
                break

        # Transformation
        R = _get_rotation_matrix(thetax, thetay)

        # Applying transformation to the mesh
        Vc[:, -1] += dz # Z translation
        Vc = np.transpose(np.dot(R, Vc.T)) # rotation

        if anim:
            mm.write_VTP('mesh%u.vtp'%iter, Vc, Fc)

        # Applying transformation to the center of gravity
        CGc[-1] += dz
        CGc = np.dot(R, CGc)
        xg, yg, zg = CGc

        # origin[-1] += dz
        # origin = np.dot(R, origin)

        # Computing hydrostatics
        hs_output = compute_hydrostatics(Vc, Fc, zg, rho_water=rho_water, grav=grav, verbose=False)

        # Computing characteristic length
        # TODO: L et l doivent etre calcules a partir du polygone d'intersection !!!
        xmin, xmax = hs_output['Sf_x_lim']
        ymin, ymax = hs_output['Sf_y_lim']
        L, l = xmax-xmin, ymax-ymin

        xb, yb, zb = hs_output['B']
        KH = hs_output['KH'][2:5, 2:5]
        if KH[1, 1] < 0.:
            xstable = False
            print 'Unstable in roll'
        else:
            xstable = True

        if KH[2, 2] < 0.:
            ystable = False
            print 'Unstable in pitch'
        else:
            ystable = True

        Vw = hs_output['Vw']

        # Correcting KH
        # corr = grav*zg * (rho_water*Vw - disp)
        # KH[3, 3] += corr
        # KH[4, 4] += corr

        # Residual to compensate
        rhogV = rhog*Vw

        res = np.array([rhogV - mg,
                        rhogV*yb - mg*yg,
                       -rhogV*xb + mg*xg])

        # Convergence criteria
        scale = np.array([mg, mg*l, mg*L])

        if np.all(np.fabs(res/scale) < reltol):
            # Convergence at an equilibrium
            if xstable and ystable:
                # Stable equilibrium
                code = 1
                break
            else:
                # code = 2
                # if not ystable:
                #     print '\nConvergence reach at an unstable configuration in PITCH at iteration %u' % iter
                #     print '\t--> Keep going iterations with an opposite orientation\n'
                #     # Choose a new position at 180 deg in pitch
                #     R = _get_rotation_matrix(0., math.pi)
                #     Vc = np.transpose(np.dot(R, Vc.T)) # rotation
                #     CGc = np.dot(R, CGc)
                #
                # if not xstable:
                #     print '\nConvergence reach at an unstable configuration in ROLL at iteration %u' % iter
                #     print '\t--> Keep going iterations with an opposite orientation\n'
                #     # Choose a new position at 180 deg in roll
                #     R = _get_rotation_matrix(math.pi, 0.)
                #     Vc = np.transpose(np.dot(R, Vc.T)) # rotation
                #     CGc = np.dot(R, CGc)

                iter += 1
                continue

        dz_old = dz
        thetax_old = thetax
        thetay_old = thetay

        # ESSAI d'utilisation de KH diagonale seulement
        # KH = np.diag(np.diag(KH))

        # Computing correction
        dz, thetax, thetay = np.linalg.solve(KH, res)

        # Adaptativity
        # if dz*dz_old < 0.:
        #     z_relax /= 2
        #
        # if thetax*thetax_old < 0.:
        #     theta_relax_x /= 2
        #
        # if thetay*thetay_old < 0.:
        #     theta_relax_y /= 2

        # Relaxation
        if math.fabs(dz) > z_relax:
            dz = math.copysign(z_relax, dz)

        if math.fabs(thetax) > theta_relax_x:
            thetax = math.copysign(theta_relax_x, thetax)

        if math.fabs(thetay) > theta_relax_y:
            thetay = math.copysign(theta_relax_y, thetay)
        # print dz, thetax, thetay

        # if z_relax < 1e-4 and theta_relax_x < 1e-3 and theta_relax_y < 1e-3:
        #     print "Can't converge to a solution better than %f %%" % (np.max(np.fabs(res/scale))*100)
        #     break

        iter += 1

    # Zeroing xcog and ycog
    Vc[:, 0] -= CGc[0]
    Vc[:, 1] -= CGc[1]
    CGc[0] = CGc[1] = 0.

    # origin[0] -= CGc[0]
    # origin[1] -= CGc[1]

    if verbose:
        if code == 0:
            # Max iterations
            print 'No convergence after %u iterations' % itermax
        elif code == 1:
            print 'Convergence reached after %u iterations at %f %% of the displacement.' % (iter, reltol*100)
        elif code == 2:
            print 'Convergence reached but at an unstable configuration'

    # print origin + CGc
    # print CG

    return Vc, Fc, CGc



def _get_Sf_Vw(V, F):
    """
    Computes only the flotation surface area and the immersed volume of the mesh

    Parameters
    ----------
    V : ndarray
        Array of the mesh vertices
    F : ndarray
        Array of the mesh connectivities

    Returns
    -------
    Vc : ndarray
    Fc : ndarray
    Sf : float
    Vw : float
    """
    Vc, Fc, clip_infos = mm.clip_by_plane(V, F, mm.Plane(), infos=True)

    areas, normals, centers = mm.get_all_faces_properties(Vc, Fc)

    Vw = (areas*(normals*centers).sum(axis=1)).sum()/3. # Formule approchee mais moyennee

    Sf = 0.
    polygons = clip_infos['PolygonsNewID']
    for polygon in polygons:
        polyverts = Vc [polygon]

        x = polyverts[:, 0]
        y = polyverts[:, 1]

        Sf += ((np.roll(y, -1)-y) * (np.roll(x, -1)+x)).sum()
    Sf *= 0.5

    return Vc, Fc, Sf, Vw


def set_displacement(V, F, disp, rho_water=1023, grav=9.81, abs_tol= 1., itermax=25, verbose=False):
    """
    Displaces mesh at a prescribed displacement and returns

    Parameters
    ----------
    V : ndarray
        Array of the mesh vertices
    F : ndarray
        Array of the mesh connectivities
    disp : float
        Mass displacement of the hull (in tons)
    rho_water : float
        Density of water (in kg/m**3)
    grav : float
        Acceleration of gravity (in m/s**2)
    abs_tol : float
        absolute tolerance on the volume (in m**3)
    itermax : int
        Maximum number of iterations
    verbose : bool
        False by default. If True, a convergence report is printed on the screen
    Returns
    -------
        dz : float
            The quantity necessary to
        Vc : ndarray
        Fc : ndarray
    """

    dz = 0.
    iter = 0
    rho = rho_water*1e-3 # To convert automatically the displacement in kg in calculations

    while True:
        if iter == itermax:
            if verbose:
                print 'No convergence of the displacement after %u iterations' % itermax
            break

        Vc = V.copy()
        Fc = F.copy()

        # Translating the mesh
        mm.translate_1D(Vc, dz, 'z') # TODO: le faire directement

        # Getting hydrostatics
        # TODO: ecrire une fonction privee permettant de ne calculer que Sf et Vw ainsi que Vc et Fc
        # On pourra alors ne plus utiliser grav dans cette fonction !
        Vc, Fc, Sf, Vw = _get_Sf_Vw(Vc, Fc)
        # hs_output = compute_hydrostatics(Vc, Fc, 0, rho_water=rho_water, grav=grav, verbose=False)
        # Sf = hs_output['Sf']
        # Vw = hs_output['Vw']

        dV = Vw - disp/rho
        if math.fabs(dV) < abs_tol:
            break

        iter += 1

        # TODO: mettre une relaxation sur la MAJ de dz...
        # Updating dz
        dz += dV / Sf

    if verbose:
        print 'Displacement obtained after %u iterations' % iter
        print 'The mesh has been displaced by %f m along z to reach a displacement of %f tons' % (dz, rho_water*Vw*1e-3)

    return dz, Vc, Fc


def compute_hydrostatics(V, F, zg, rho_water=1023, grav=9.81, verbose=False):
    """
    Computes the hydrostatics properties of a mesh.

    Parameters
    ----------
    V : ndarray
        Array of the mesh vertices
    F : ndarray
        Array of the mesh connectivities
    zg : float
        Vertical position of the center of gravity
    rho_water : float
        Density of water (default: 2013 kg/m**3)
    grav : float
        Gravity acceleration (default: 9.81 m/s**2)
    verbose : bool, optional
        False by default. If True, a hydrostatic report is printed on screen.

    Returns
    -------
    output : dict
        Dictionary containing the hydrosatic properties and whose keys are:
            'Sw': (float) Wet surface area (in m**2)
            'Vw': (float) Immersed volume (in m**3)
            'disp': (float) Displacement (in Tons)
            'B': (ndarray) Coordinates of the buoyancy center (in m)
            'F': (ndarray) Coordinates of the center of flotation (in m)
            'Sf': (float) Area of the flotation surface (in m**2)
            'r': (float) Transverse metacentric radius (in m)
            'R': (float) Longitudinal metacentrix radius (in m)
            'GMx': (float) Transverse metacentric height (in m)
            'GMy': (float) Longitudinal metacentric height (in m)
            'KH': (ndarray) Hydrostatic stiffness matrix
            'Vc': (ndarray) Array of hydrostatic mesh vertices
            'Fc': (ndarray) Array of hydrostatic mesh connectivities
    """

    eps = 1e-4 # For zeroing tiny coefficients in the hydrostatic stiffness matrix

    # Clipping the mesh by the Oxy plane
    plane = mm.Plane() # Oxy plane
    try:
        Vc, Fc, clip_infos = mm.clip_by_plane(V, F, plane, infos=True)
    except:
        show(V, F)
        raise Exception, 'Hydrostatic module only work with watertight hull. Please consider using the --sym option.'

    # Calculs des pptes des facettes
    areas, normals, centers = mm.get_all_faces_properties(Vc, Fc)

    # Calcul surface mouillee
    Sw = areas.sum()

    # Calcul volume de carene
    Vw = (areas*(normals*centers).sum(axis=1)).sum()/3. # Formule approchee mais moyennee

    # Buoyancy center calculation
    xb = (areas * normals[:, 1] * centers[:, 1] * centers[:, 0]).sum() / Vw
    yb = (areas * normals[:, 2] * centers[:, 2] * centers[:, 1]).sum() / Vw
    zb = (areas * normals[:, 1] * centers[:, 1] * centers[:, 2]).sum() / Vw

    # Computing quantities from intersection polygons
    sigma0 = 0. # \int_{Sf} dS = Sf
    sigma1 = 0. # \int_{Sf} x dS
    sigma2 = 0. # \int_{Sf} y dS
    sigma3 = 0. # \int_{Sf} xy dS
    sigma4 = 0. # \int_{Sf} x^2 dS
    sigma5 = 0. # \int_{Sf} y^2 dS

    xmin = []
    xmax = []
    ymin = []
    ymax = []

    polygons = clip_infos['PolygonsNewID']
    for polygon in polygons:
        polyverts = Vc[polygon]

        # TODO: voir si on conserve ce test...
        if np.any(np.fabs(polyverts[:, 2]) > 1e-3):
            print 'The intersection polygon is not on the plane z=0'

        xi, yi = polyverts[0, :2]
        for (xii, yii) in polyverts[1:, :2]:

            dx = xii - xi
            dy = yii - yi
            px = xi + xii
            py = yi + yii
            a = xi*xi + xii*xii

            sigma0 += dy*px
            sigma1 += dy * (px*px-xi*xii)
            sigma2 += dx * (py*py-yi*yii)
            sigma3 += dy * ( py*a + 2*px*(xi*yi + xii*yii) )
            sigma4 += dy * a * px
            sigma5 += dx * (yi*yi + yii*yii) * py

            xi, yi = xii, yii

        xmin.append(polyverts[:, 0].min())
        xmax.append(polyverts[:, 0].max())
        ymin.append(polyverts[:, 1].min())
        ymax.append(polyverts[:, 1].max())

    Sf_x_lim = [min(xmin), max(xmax)]
    Sf_y_lim = [min(ymin), max(ymax)]

    sigma0 /= 2
    sigma1 /= 6
    sigma2 /= -6
    sigma3 /= 24
    sigma4 /= 12
    sigma5 /= -12

    # Flotation surface
    Sf = sigma0

    rhog = rho_water * grav

    # Stiffness matrix coefficients that do not depend on the position of the gravity center
    S33 = rhog * Sf
    S34 = rhog * sigma2
    S35 = -rhog * sigma1
    S45 = -rhog * sigma3

    # Metacentric radius (Bouguer formulae)
    r = sigma5 / Vw # Around Ox
    R = sigma4 / Vw # Around Oy

    # Metacentric height
    a = zg - zb # BG
    GMx = r - a
    GMy = R - a
    # print 'GMx=%f; GMy=%f' % (GMx, GMy)

    # Stiffness matrix coefficients that depend on the position of the gravity center
    S44 = rhog*Vw * GMx
    S55 = rhog*Vw * GMy

    # Assembling matrix
    KH = np.zeros((6, 6))
    KH[2, 2] = S33
    KH[3, 3] = S44
    KH[4, 4] = S55
    KH[2, 3] = S34
    KH[3, 2] = S34
    KH[2, 4] = S35
    KH[4, 2] = S35
    KH[3, 4] = S45
    KH[4, 3] = S45

    # Zeroing tiny coefficients
    KH[np.fabs(KH) < eps] = 0.

    # Flotation center F:
    xF = -S35/S33
    yF =  S34/S33

    #Displacement
    disp = rho_water * Vw * 1e-3 # in tons

    if verbose:
        # Data for DNV standards
        GM_min = 0.15

        print '\nWet surface = %f (m**2)\n' % Sw
        print 'Immersed volume = %f (m**3)\n' % Vw
        print 'Displacement = %f (tons)\n' % disp
        print 'Buoyancy center (m): xb=%f, yb=%f, zb=%f\n' % (xb, yb, zb)
        print 'Flottation surface = %f (m**2)\n' % Sf
        print 'Flotation center (m): xf=%f, yf=%f\n' % (xF, yF)
        print 'Transverse metacentric radius = %f (m)\n' % r
        print 'Longitudinal metacentric radius = %f (m)\n' % R

        print 'Transverse metacentric height GMx = %f (m)' % GMx
        if GMx < 0.:
            print '\t --> Unstable in roll !'
            print '\t     To be stable, you should have at least zg < %f (m)' % (r+zb)
            print '\t     DNV Standards say : zg < %f (m) to get GMx > %f m\n' % (r+zb-GM_min, GM_min)
        else:
            print '\t --> Stable in roll\n'

        print 'Longitudinal metacentric height GMy = %f (m)' % GMy
        if GMy < 0.:
            print '\t --> Unstable in pitch !'
            print '\t     To be stable, you should have at least zg < %f (m)' % (R+zb)
            print '\t     DNV Standards say : zg < %f (m) to get GMy > %f m\n' % (R+zb-GM_min, GM_min)
        else:
            print '\t --> Stable in pitch\n'

        print 'Hydrostatic stiffness matrix:'
        for line in KH:
            print '%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E' % (line[0], line[1], line[2], line[3], line[4], line[5])

    # Output data
    output = dict()
    output['Sw'] = Sw
    output['Vw'] = Vw
    output['disp'] = disp
    output['B'] = np.array([xb, yb, zb], dtype=np.float)
    output['F'] = np.array([xF, yF, 0.], dtype=np.float)
    output['Sf'] = Sf
    output['r'] = r
    output['R'] = R
    output['GMx'] = GMx
    output['GMy'] = GMy
    output['KH'] = KH
    output['Sf_x_lim'] = Sf_x_lim
    output['Sf_y_lim'] = Sf_y_lim
    output['Vc'] = Vc
    output['Fc'] = Fc
    output['polygons'] = polygons

    return output


if __name__ == '__main__':

    # The following code are only for testng purpose
    import pickle

    V, F = mm.load_VTP('meshmagick/tests/data/SEAREV.vtp')
    
    hs_output = compute_hydrostatics(V, F, 0)
    
    pickle.dump(hs_output, open('hs_ouput_r177.p', 'w'))
    
