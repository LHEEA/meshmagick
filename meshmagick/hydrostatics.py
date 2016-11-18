
__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-2015, Ecole Centrale de Nantes"
__credits__ = "Francois Rongere"
__licence__ = "CeCILL"
__version__ = "1.0"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@ec-nantes.fr"
__status__ = "Development"

import numpy as np
import math
from warnings import warn

import mmio
from mesh_clipper import MeshClipper

# Data for DNV standards
GM_min = 0.15

# TODO: pass this code to a version that only rotates the free surface plane and not the mesh...
# TODO: throw the transformation needed to get the equilibrium from initial position after a successful equilibrium computation
# TODO: make the mesh not to diverge from principal axis

class AdditionalForce(object):
    def __init__(self, point=[0, 0, 0], value=[0, 0, 0], name='', mode='relative'):
        
        assert len(point) == 3
        assert len(value) == 3
        assert mode in ('relative', 'absolute')
        
        self.point = np.asarray(point, dtype=np.float)
        self.value = np.asarray(value, dtype=np.float)
        self.name = str(name)
        self.mode = mode
    
    def __str__(self):
        str_repr = "Force: %s\n\tPoint: %s\n\tValue: %s\n\tMode: %s" % (self.name, self.point, self.value, self.mode)
        return str_repr
    
    def update(self, dz=0., rot=np.eye(3, dtype=np.float)):
        self.point[2] += dz
        self.point = np.dot(rot, self.point)
        if self.mode == 'relative':
            self.value = np.dot(rot, self.value)
        return
    
    def update_xy(self, dx, dy):
        self.point[0] += dx
        self.point[1] += dy
        return
    
    @property
    def hs_force(self):
        Fz = self.value[2]
        moment = np.cross(self.point, self.value)
        Mx, My = moment[:2]
        return np.array([Fz, Mx, My], dtype=np.float)
    
    def reset(self): # TODO: a appeler depuis le reset de Hydrostatics
        raise NotImplementedError
    

class Hydrostatics(object):
    """
    Class to perform hydrostatic computations on meshes.
    
    Rules:
    ------
    *** At class instantiation, no modification on the mesh position is done and hydrostatic properties are accessible directly for this position. The mesh may not be at equilibrium thought and residual force balance between hydrostatic force and gravity force may be not null. To trig the equilibrium computations, we have to call the equilibrate() method.
    
       Note that it is useless to call any method to get the hydrostatic properties up to date as it is done internally and automatically at each modification.
    
    *** Intrinsic mesh properties are
    
        - _mass
        - cog
        - zg
        - rho_water
        - gravity
        
        Setting one of these property by its setter will trig a search for a new mesh position that fullfill the hydrostatic equilibrium.
        Note that zg is mandatory to get the hydrostatic stiffness matrix. By default, be carefull that it is set to zero at instantiation.
    
    *** Modification of the displacement will produce the mesh to move along z so that
    it reaches the correct displacement. The _mass is then set accordingly and equals the displacement. Note that this is different from modifying the _mass and
    
    
    
    
    
    
    *** Setting the buoyancy center will trig a search for a position of the gravity center that fullfill this hydrostatic property, while keeping the correct displacement. NOT IMPLEMENTED, THIS IS AN ADVANCED FEATURE
    
    """
    def __init__(self, mesh, CG=np.zeros(3, dtype=np.float), zg=0., mass=None, rho_water=1023, grav=9.81,
                 verbose=False, animate=False):
        
        # super(HSmesh, self).__init__(vertices, faces, name=None)
        self.backup = dict()
        
        # Ajout
        # init_mesh = mesh.Mesh(vertices, faces, name=name)
        self.backup['init_mesh'] = mesh
        self.mesh = mesh.copy()
        # FIN
        
        CG = np.array(CG, dtype=np.float)
        assert CG.shape[0] == 3
        self.backup['gravity_center'] = CG
        self._gravity_center = CG
        self._rho_water = rho_water
        self._gravity = grav
        
        self._rhog = rho_water * grav
        
        # Solver parameters
        self._solver_parameters = {'reltol': 1e-2,
                                   'itermax': 100,
                                   'max_nb_restart': 10,
                                   'theta_relax': 2,
                                   'z_relax': 0.1,
                                   'stop_at_unstable': False}
        
        
        # TODO: ajouter le calcul du tirant d'eau et d'air
        self.hs_data = dict()
        self._update_hydrostatic_properties()
        
        if mass:
            self.mass = mass * 1000. # conversion in kg
        else:
            self._mass = self.hs_data['disp_mass']
        
        self._mg = self._mass * self._gravity
        
        self._verbose = verbose
        
        self.animate = animate
        
        self._rotation = np.eye(3, dtype=np.float)
        
        self.additional_forces = []
        
    
    @property
    def verbose(self):
        return self._verbose
    
    def verbose_on(self):
        self._verbose = True
        return
    
    def verbose_off(self):
        self._verbose = False
        return
    
    @property
    def gravity(self):
        return self._gravity
    
    @gravity.setter
    def gravity(self, value):
        self._gravity = float(value)
        self._rhog = self._rho_water * value
        self._mg = self._mass * value
        self._update_hydrostatic_properties()
        return
    
    @property
    def rho_water(self):
        return self._rho_water
    
    @rho_water.setter
    def rho_water(self, value):
        self._rho_water = float(value)
        self._rhog = value * self._gravity
        self._update_hydrostatic_properties()
        return
    
    @property
    def mass(self):
        return self._mass / 1000.
    
    @mass.setter
    def mass(self, value):
        
        self._mass = float(value) * 1000. # Conversion from tons to kg
        self._mg = self._mass * self._gravity # SI units
        
        if self.is_sinking():
            raise ValueError('%s is sinking as it is too heavy.' % self.mesh.name)
        
        return
    
    def _max_displacement(self):
        return self._rho_water * self.mesh._compute_volume() # in kg
    
    def is_sinking(self):
        if self._mass > self._max_displacement():
            return True
        else:
            return False
    
    @property
    def gravity_center(self):
        return self._gravity_center
    
    @gravity_center.setter
    def gravity_center(self, value, update=True):
        value = np.asarray(value, dtype=np.float)
        self.backup['gravity_center'] = value
        assert value.shape[0] == 3
        self._gravity_center = value
        if update:
            self._update_hydrostatic_properties()
        return
    
    @property
    def zg(self):
        return self._gravity_center[-1]
    
    @zg.setter
    def zg(self, value, update=True):
        self._gravity_center[-1] = float(value)
        if update:
            self._update_hydrostatic_properties()
        return
    
    @property
    def wetted_surface_area(self):
        return self.hs_data['Sw']
    
    @property
    def displacement_volume(self):
        return self.hs_data['disp_volume']
    
    @property
    def displacement(self):
        return self.hs_data['disp_mass'] / 1000.
    
    @property
    def buoyancy_center(self):
        return self.hs_data['buoy_center']
    
    @property
    def flotation_surface_area(self):
        return self.hs_data['Aw']
    
    @property
    def flotation_center(self):
        return self.hs_data['flot_center']
    
    @property
    def transversal_metacentric_radius(self):
        return self.hs_data['r']
    
    @property
    def longitudinal_metacentric_radius(self):
        return self.hs_data['R']
    
    @property
    def transversal_metacentric_height(self):
        return self.hs_data['GMx']
    
    @property
    def longitudinal_metacentric_height(self):
        return self.hs_data['GMy']
    
    @property
    def hydrostatic_stiffness_matrix(self):
        return self.hs_data['KH']

    @property
    def hydrostatic_mesh(self):
        return self.hs_data['hs_mesh']
    
    def reset(self):
        
        # TODO: Utiliser plutot la rotation generale pour retrouver le maillage initial
        
        self.mesh = self.backup['init_mesh']
        self._gravity_center = self.backup['gravity_center']
        
        self._update_hydrostatic_properties()
        
        self._rotation = np.eye(3, dtype=np.float)
        return
    
    def is_stable_in_roll(self):
        return self.transversal_metacentric_height > 0.
        
    def is_stable_in_pitch(self):
        return self.longitudinal_metacentric_height > 0.
    
    def isstable(self):
        return self.is_stable_in_pitch() and self.is_stable_in_roll()
    
    def is_at_equilibrium(self):
        
        residual = self.residual
        
        mg = self._mg
        B = self.hs_data['B']
        Lpp = self.hs_data['Lpp']
        scale = np.array([mg, mg * B, mg * Lpp])
        
        if np.all(np.fabs(residual/scale) < self._solver_parameters['reltol']):
            return True
        else:
            return False
    
    @property
    def delta_fz(self):
        fz, _, _ = self.residual
        return fz
    
    @property
    def delta_mx(self):
        _, mx, _ = self.residual
        return mx
    
    @property
    def delta_my(self):
        _, _, my = self.residual
        return my
    
    # TODO: create a Hydrostatic stiffness matrix class
    @property
    def S33(self):
        return self.hs_data['KH'][0, 0]
    
    @property
    def S34(self):
        return self.hs_data['KH'][0, 1]
    
    @property
    def S35(self):
        return self.hs_data['KH'][0, 2]
    
    @property
    def S44(self):
        return self.hs_data['KH'][1, 1]
    
    @property
    def S45(self):
        return self.hs_data['KH'][1, 2]
    
    @property
    def S55(self):
        return self.hs_data['KH'][2, 2]
    
    @property
    def reltol(self):
        return self._solver_parameters['reltol']
    
    @reltol.setter
    def reltol(self, value):
        value = float(value)
        assert value > 0. and value < 1
        self._solver_parameters['reltol'] = value
        return
    
    @property
    def theta_relax(self):
        return self._solver_parameters['theta_relax']
    
    @theta_relax.setter
    def theta_relax(self, value):
        value = float(value)
        assert value > 0. and value < 10
        self._solver_parameters['theta_relax'] = value
        return
        
    @property
    def z_relax(self):
        return self._solver_parameters['z_relax']
    
    @z_relax.setter
    def z_relax(self, value):
        value = float(value)
        assert value > 0. and value < 1.
        self._solver_parameters['z_relax'] = value
        return
    
    @property
    def max_iterations(self):
        return self._solver_parameters['itermax']
    
    @max_iterations.setter
    def max_iterations(self, value):
        value = int(value)
        assert value > 0
        self._solver_parameters['itermax'] = value
        return
    
    @property
    def max_restart(self):
        return self._solver_parameters['max_nb_restart']
    
    @max_restart.setter
    def max_restart(self, value):
        value = int(value)
        assert value >= 0
        self._solver_parameters['max_nb_restart'] = value
        return
    
    @property
    def allow_unstable(self, value):
        return self._solver_parameters['stop_at_unstable']
    
    @allow_unstable.setter
    def allow_unstable(self, value):
        value = bool(value)
        self._solver_parameters['stop_at_unstable'] = value
        return
    
    def allow_unstable_on(self):
        self.allow_unstable = True
        return
    
    def allow_unstable_off(self):
        self.allow_unstable = False
        return
    
    def rotate(self, angles, update=True):
        # R = super(HSmesh, self).rotate(angles)
        R = self.mesh.rotate(angles)
        self._gravity_center = np.dot(R, self._gravity_center)
        if update:
            self._update_hydrostatic_properties()
        self._rotation = np.dot(R, self._rotation)
        return R
    
    def rotate_x(self, thetax, update=True):
        return self.rotate([thetax, 0., 0.], update=update)
        
    def rotate_y(self, thetay, update=True):
        return self.rotate([0., thetay, 0.], update=update)
    
    def rotate_z(self, thetaz, update=True):
        return self.rotate([0., 0., thetaz], update=update)
    
    def translate_x(self, tx, update=True):
        self.mesh.translate_x(tx)
        # super(HSmesh, self).translate_x(tx)
        self._gravity_center[0] += tx
        if update:
            self._update_hydrostatic_properties()
        return
        
    def translate_y(self, ty, update=True):
        self.mesh.translate_y(ty)
        self._gravity_center[1] += ty
        if update:
            self._update_hydrostatic_properties()
        return
        
    def translate_z(self, tz, update=True):
        self.mesh.translate_z(tz)
        self._gravity_center[2] += tz
        if update:
            self._update_hydrostatic_properties()
        return
        
    def translate(self, t, update=True):
        super(Hydrostatics, self).translate(t)
        self._gravity_center += t
        if update:
            self._update_hydrostatic_properties()
        return
     
    def _update_hydrostatic_properties(self):
        """
        Computes the hydrostatics properties of a mesh.

        Parameters
        ----------
        mymesh : Mesh
            Mesh object for hydrostatic computations
        zg : float
            Vertical position of the center of gravity
        rho_water : float
            Density of water (default: 2013 kg/m**3)
        grav : float
            Gravity acceleration (default: 9.81 m/s**2)
        verbose : bool, optional
            False by default. If True, the hydrostatic report is printed on screen.

        Returns
        -------
        output : dict
            Dictionary containing the hydrosatic properties and whose keys are:
                'Sw': (float) Wet surface area (in m**2)
                'disp_volume': (float) Immersed volume (in m**3)
                'disp_mass': (float) Displacement (in Tons)
                'B': (ndarray) Coordinates of the buoyancy center (in m)
                'F': (ndarray) Coordinates of the center of flotation (in m)
                'Aw': (float) Area of the flotation surface (in m**2)
                'r': (float) Transverse metacentric radius (in m)
                'R': (float) Longitudinal metacentrix radius (in m)
                'GMx': (float) Transverse metacentric height (in m)
                'GMy': (float) Longitudinal metacentric height (in m)
                'KH': (ndarray) Hydrostatic stiffness matrix

            """
    
        eps = 1e-4  # For zeroing tiny coefficients in the hydrostatic stiffness matrix
    
        # Clipping the mesh by the Oxy plane
        # TODO: ne pas recreer une instance de MeshClipper a fois mais la stocker et mettre a jour son maillage ainsi
        #  que ses donnees

        clipper = MeshClipper(self.mesh, assert_closed_boundaries=True, verbose=False)
        clipped_mesh = clipper.clipped_mesh  # TODO: enregistrer le maillage coupe !!!
        
        # Retrieving faces properties for the clipped mesh
        areas = clipped_mesh.faces_areas
        normals = clipped_mesh.faces_normals
        centers = clipped_mesh.faces_centers
        
        if np.any(areas == 0.): # TODO: bloc a retirer
            print 'probleme de facette'
            self.mesh.quick_save()
            raise Exception
        
        # Storing clipped mesh
        self.hs_data['hs_mesh'] = clipped_mesh
        
        # Wetted surface area
        Sw = areas.sum()
    
        # TODO: utiliser des formules analytiques et non approchees comme celles-ci !
        # Volume displacement
        disp_volume = (areas * (normals * centers).sum(axis=1)).sum() / 3.  # Formule approchee mais moyennee
    
        # Buoyancy center calculation
        xb = (areas * normals[:, 1] * centers[:, 1] * centers[:, 0]).sum() / disp_volume
        yb = (areas * normals[:, 2] * centers[:, 2] * centers[:, 1]).sum() / disp_volume
        zb = (areas * normals[:, 1] * centers[:, 1] * centers[:, 2]).sum() / disp_volume
    
        # Computing quantities from intersection polygons
        sigma0 = 0.  # \int_{Aw} dS = Aw
        sigma1 = 0.  # \int_{Aw} x dS
        sigma2 = 0.  # \int_{Aw} y dS
        sigma3 = 0.  # \int_{Aw} xy dS
        sigma4 = 0.  # \int_{Aw} x^2 dS
        sigma5 = 0.  # \int_{Aw} y^2 dS
    
        xmin = []
        xmax = []
        ymin = []
        ymax = []
    
        polygons = clipper.closed_polygons
        for polygon in polygons:
            polyverts = clipper.clipped_crown_mesh.vertices[polygon]
        
            # TODO: voir si on conserve ce test...
            if np.any(np.fabs(polyverts[:, 2]) > 1e-3):
                print 'The intersection polygon is not on the plane z=0'
        
            xi, yi = polyverts[0, :2]
            for (xii, yii) in polyverts[1:, :2]:
                dx = xii - xi
                dy = yii - yi
                px = xi + xii
                py = yi + yii
                a = xi * xi + xii * xii
            
                sigma0 += dy * px
                sigma1 += dy * (px * px - xi * xii)
                sigma2 += dx * (py * py - yi * yii)
                sigma3 += dy * (py * a + 2 * px * (xi * yi + xii * yii))
                sigma4 += dy * a * px
                sigma5 += dx * (yi * yi + yii * yii) * py
            
                xi, yi = xii, yii
        
            xmin.append(polyverts[:, 0].min())
            xmax.append(polyverts[:, 0].max())
            ymin.append(polyverts[:, 1].min())
            ymax.append(polyverts[:, 1].max())
    
        minx, maxx = [min(xmin), max(xmax)]
        miny, maxy = [min(ymin), max(ymax)]
    
        sigma0 /= 2
        sigma1 /= 6
        sigma2 /= -6
        sigma3 /= 24
        sigma4 /= 12
        sigma5 /= -12
    
        # Flotation surface
        Aw = sigma0
    
        # Stiffness matrix coefficients that do not depend on the position of the gravity center
        rhog = self._rhog
        S33 = rhog * Aw
        S34 = rhog * sigma2
        S35 = -rhog * sigma1
        S45 = -rhog * sigma3
    
        # Metacentric radius (Bouguer formulae)
        r = sigma5 / disp_volume  # Around Ox
        R = sigma4 / disp_volume  # Around Oy
    
        # Metacentric height
        a = self.zg - zb  # BG
        GMx = r - a
        GMy = R - a
    
        # Stiffness matrix coefficients that depend on the position of the gravity center
        S44 = rhog * disp_volume * GMx
        S55 = rhog * disp_volume * GMy
    
        # Assembling stiffness matrix
        KH = np.array([[S33, S34, S35],
                       [S34, S44, S45],
                       [S35, S45, S55]], dtype=np.float)
        
        # Zeroing tiny coefficients
        KH[np.fabs(KH) < eps] = 0.
    
        # Flotation center F:
        xF = -S35 / S33
        yF = S34 / S33
    
        # Storing data
        self.hs_data['Sw'] = Sw
        self.hs_data['disp_volume'] = disp_volume
        self.hs_data['disp_mass'] = self._rho_water * disp_volume
        self.hs_data['buoy_center'] = np.array([xb, yb, zb], dtype=np.float)
        self.hs_data['flot_center'] = np.array([xF, yF, 0.], dtype=np.float)
        self.hs_data['Aw'] = Aw
        self.hs_data['r'] = r
        self.hs_data['R'] = R
        self.hs_data['GMx'] = GMx
        self.hs_data['GMy'] = GMy
        self.hs_data['KH'] = KH
        self.hs_data['Lpp'] = maxx - minx
        self.hs_data['B'] = maxy - miny
    
        return
    
    def add_force(self, force):
        assert isinstance(force, AdditionalForce)
        self.additional_forces.append(force)
        return
    
    @property
    def scale(self):
        mg = self._mg
        B = self.hs_data['B']
        Lpp = self.hs_data['Lpp']
        scale = np.array([mg, mg * B, mg * Lpp])
        return scale

    @property
    def residual(self):
        rhogV = self._rhog * self.hs_data['disp_volume']
        mg = self._mg
    
        xb, yb, _ = self.hs_data['buoy_center']
        xg, yg, _ = self._gravity_center
    
        residual = np.array([rhogV - mg,
                             rhogV * yb - mg * yg,
                             -rhogV * xb + mg * xg])
        
        for force in self.additional_forces:
            residual += force.hs_force
        # # Essai
        # Fz = hs.force['force'][2]
        # M = np.cross(hs.force['point'], hs.force['force'])
        # Mx = M[0]
        # My = M[1]
        #
        # residual += [Fz, Mx, My]
        #
        # # FIN
        
        return residual
    
    def set_displacement(self, disp):
        """
        Displaces mesh at a prescribed displacement and returns

        Parameters
        ----------
        mymesh : Mesh
            Mesh to be used for computations
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
        
        if self.verbose:
            print "\nComplying with a displacement of %.3f tons" % disp
            print "----------------------------------------------"
        
        self.mass = disp
        
        itermax = self._solver_parameters['itermax']
        reltol = self._solver_parameters['reltol']
        z_relax = self._solver_parameters['z_relax']
        
        dz = 0.
        iter = 0

        while True:
            # print iter
            if iter == itermax:
                if self.verbose:
                    print '\t-> No convergence of the displacement after %u iterations' % itermax
                break
            
            # Translating the mesh
            self.mesh.translate_z(dz)
            self._gravity_center[2] += dz
            
            for force in self.additional_forces:
                force.update(dz=dz)
            
            self._update_hydrostatic_properties()
    
            residual = self.delta_fz
            if math.fabs(residual / self._mg) < reltol:
                if self.verbose:
                    print '\t-> Convergence obtained after %u iterations' % iter
                break
            
            dz = residual / (self._rhog * self.flotation_surface_area)
            if math.fabs(dz) > z_relax:
                dz = math.copysign(z_relax, dz)
            iter += 1
        
        # self.mass = self.hs_data['disp_mass']
        
        return

    def equilibrate(self, init_disp=True):
        
        # Initial displacement equilibrium
        if init_disp:
            if self.verbose:
                print "First placing the mesh at the target displacement"
                print "-------------------------------------------------"
            self.set_displacement(self.mass)
        
        if self.verbose:
            print '\nComputing equilibrium from initial condition.'
            print '----------------------------------------------'
        
        # Retrieving solver parameters
        z_relax = self._solver_parameters['z_relax']
        theta_relax = self._solver_parameters['theta_relax']
        theta_relax_x = math.radians(theta_relax)
        theta_relax_y = math.radians(theta_relax)
        
        itermax = self._solver_parameters['itermax']
        reltol = self._solver_parameters['reltol']
        max_nb_restart = self._solver_parameters['max_nb_restart']
        
        rhog = self._rhog
        mg = self._mg
        
        # Initialization
        dz = thetax = thetay = 0.
        iter = nb_restart = 0
        unstable_config = False
    
        while True:
            # print '\n', iter
            if unstable_config or ( iter == (nb_restart + 1) * (itermax-1) ):
                if self.verbose:
                    if unstable_config:
                        print 'Unstable equilibrium reached.'
                        print '\t-> Keep searching a stable configuration by random restart number %u.' % (nb_restart+1)
                    else:
                        print 'Failed to find an equilibrium configuration with these initial conditions in %u ' \
                              'iterations.' % itermax
                        print '\t-> Keep searching by random restart number %u.' % (nb_restart+1)
                
                unstable_config = False
                
                # Max iterations reached
                if nb_restart < max_nb_restart-1:
                    nb_restart += 1
                    # Random on the position of the body
                    # print 'Max iteration reached: restarting for the %uth time with random orientation' % nb_restart
                    thetax, thetay = np.random.rand(2) * math.pi
                    R = self.mesh.rotate([thetax, thetay, 0.])
                    self._gravity_center = np.dot(R, self._gravity_center)
                    self._rotation = np.dot(R, self._rotation)
                    self._update_hydrostatic_properties()
                    
                    for force in self.additional_forces:
                        force.update(rot=R)
                    
                    dz = thetax = thetay = 0.
                    
                else:
                    # Max number of restart allowed. Failed to find an equilibrium configuration.
                    code = 0
                    break
        
            # Applying transformation to the mesh
            self.mesh.translate_z(dz)
            self._gravity_center[2] += dz
            
            R = self.mesh.rotate([thetax, thetay, 0.])
            self._gravity_center = np.dot(R, self._gravity_center)
            self._rotation = np.dot(R, self._rotation)
            
            # Updating force data
            for force in self.additional_forces:
                force.update(dz=dz, rot=R)
            
            self._update_hydrostatic_properties()

            # if self.animate:
            #     mmio.write_VTP('mesh%u.vtp' % iter, mymesh_c.vertices, mymesh_c.faces)
        
            residual = self.residual
            scale = self.scale
            
            if np.all(np.fabs(self.residual / scale) < reltol):
            #     # Convergence at an equilibrium
                if self.isstable():
                    # Stable equilibrium
                    code = 1
                    break
                else:
                    if self._solver_parameters['stop_at_unstable']:
                        # Unstable configuration reached
                        code = 2
                        break
                    # TODO: mettre une option permettant de choisir si on s'arrete a des configs instables
                    else:
                        # We force the restart with an other random orientation as initial condition
                        unstable_config = True
                        iter += 1
                        continue
                    
            # Computing correction
            # TODO: voir pour une resolution ne faisant pas intervenir np.linalg ?
            dz, thetax, thetay = np.linalg.solve(self.hs_data['KH'], residual)
        
            # Relaxation
            if math.fabs(dz) > z_relax:
                dz = math.copysign(z_relax, dz)
        
            if math.fabs(thetax) > theta_relax_x:
                thetax = math.copysign(theta_relax_x, thetax)
        
            if math.fabs(thetay) > theta_relax_y:
                thetay = math.copysign(theta_relax_y, thetay)
        
            iter += 1
    
        # Zeroing xcog and ycog
        self.mesh.translate([-self._gravity_center[0], -self._gravity_center[1], 0.])
        self.hs_data['buoy_center'][:2] -= self._gravity_center[:2]
        self._gravity_center[:2] = 0.
        for force in self.additional_forces:
            force.update_xy(-self._gravity_center[0], -self._gravity_center[1])
        
        if self.verbose:
            if code == 0:
                print "\t-> Maximum number of restart reached. Failed to find an equilibrum position."
            elif code == 1:
                print "Stable equilibrium reached after %u iterations and %u random restart" % (iter, nb_restart)
            elif code == 2:
                print 'Unstable equilibrium reached after %u iterations and %u random restart' % (iter, nb_restart)
        
        return code

    def get_hydrostatic_report(self):
        
        msg = '\n'
        
        msg += ('Gravity acceleration = %.2f m/s**2\n' % self.gravity)
        msg += ('Density of water = %.2f (kg/m**3)\n' % self.rho_water)
        msg += ('Wetted surface area = %.2f (m**2)\n' % self.wetted_surface_area)
        msg += ('Displacement volume = %.2f (m**3)\n' % self.displacement_volume)
        msg += ('Displacement = %.3f (tons)\n' % (self.displacement))
        xb, yb, zb = self.buoyancy_center
        xg, yg, zg = self.gravity_center
        msg += ('Buoyancy center (m): xb=%.3f, yb=%.3f, zb=%.3f\n' % (xb, yb, zb))
        msg += ('Gravity center (m):  xg=%.3f, yg=%.3f, zg=%.3f\n' % (xg, yg, zg))
        msg += ('Flottation surface = %f (m**2)\n' % self.flotation_surface_area)
        xF, yF, _ = self.flotation_center
        msg += ('Flotation center (m): xf=%.3f, yf=%.3f\n' % (xF, yF))
        r = self.transversal_metacentric_radius
        R = self.longitudinal_metacentric_radius
        GMx = self.transversal_metacentric_height
        GMy = self.longitudinal_metacentric_height
        msg += ('Transverse metacentric radius = %.3f (m)\n' % r)
        msg += ('Longitudinal metacentric radius = %.3f (m)\n' % R)
        msg += ('Transverse metacentric height GMx = %.3f (m)\n' % GMx)
        if GMx < 0.:
            msg += ('\t --> Unstable in roll !\n')
            msg += ('\t     To be stable, you should have at least zg < %.3f (m)\n' % (r + zb))
            msg += ('\t     DNV Standards say : zg < %.3f (m) to get GMx > %.3f m\n' % (r + zb - GM_min, GM_min))
        else:
            msg += ('\t --> Stable in roll\n')
    
        msg += ('Longitudinal metacentric height GMy = %.3f (m)\n' % GMy)
        if GMy < 0.:
            msg += ('\t --> Unstable in pitch !\n')
            msg += ('\t     To be stable, you should have at least zg < %.3f (m)\n' % (R + zb))
            msg += ('\t     DNV Standards say : zg < %.3f (m) to get GMy > %.3f m\n' % (R + zb - GM_min, GM_min))
        else:
            msg += ('\t --> Stable in pitch\n')
    
        msg += ('\nHydrostatic stiffness matrix:\n')
        KH = self.hydrostatic_stiffness_matrix
        for row in KH:
            msg += ('%.4E\t%.4E\t%.4E\n' % (row[0], row[1], row[2]))
        
        residual = self.residual
        msg += ('\nResidual:\n')
        msg += ('Delta Fz = %.3f N\n' % residual[0])
        msg += ('Delta Mx = %.3f Nm\n' % residual[1])
        msg += ('Delta My = %.3f Nm\n' % residual[2])
        
        rel_res = residual / self.scale
        msg += ('\nRelative residual:\n')
        msg += ('Delta Fz = %E\n' % rel_res[0])
        msg += ('Delta Mx = %E\n' % rel_res[1])
        msg += ('Delta My = %E\n' % rel_res[2])
        msg += ('Relative tolerance of the solver: %.1E\n' % self.reltol)
        
        return msg
    
    def show(self):
        # TODO: Ce n'est pas ce module qui doit savoir utiliser vtk !!!
        
        import MMviewer
        import vtk
        
        vtk_polydata = self.mesh._vtk_polydata()
        self.viewer = MMviewer.MMViewer()
        self.viewer.add_polydata(vtk_polydata)
        
        # Removing the following for the flotation plane as the base class does it now
        # # Adding the flotation plane
        self.viewer.plane_on()
        
        # # TODO: Use the add_plane method of the viewer...
        # plane = vtk.vtkPlaneSource()
        # (xmin, xmax, ymin, ymax, _, _) = self.axis_aligned_bbox
        # dx = 0.1 * (xmax-xmin)
        # dy = 0.1 * (ymax-ymin)
        #
        # plane.SetOrigin(xmin-dx, ymax+dy, 0)
        # plane.SetPoint1(xmin-dx, ymin-dy, 0)
        # plane.SetPoint2(xmax+dx, ymax+dy, 0)
        # plane.Update()
        # # self.viewer.add_polydata(plane.GetOutput(), color=[0.1, 0.9, 0.7])
        # self.viewer.add_polydata(plane.GetOutput(), color=[0., 102./255, 204./255])
        
        # Adding a point for the origin
        pO = vtk.vtkPoints()
        vO = vtk.vtkCellArray()

        iO = pO.InsertNextPoint([0, 0, 0])
        vO.InsertNextCell(1)
        vO.InsertCellPoint(iO)

        pdO = vtk.vtkPolyData()
        pdO.SetPoints(pO)
        pdO.SetVerts(vO)

        self.viewer.add_polydata(pdO, color=[0, 0, 0])
        
        # Adding the center of gravity
        pCG = vtk.vtkPoints()
        vCG = vtk.vtkCellArray()

        iCG = pCG.InsertNextPoint(self.gravity_center)
        vCG.InsertNextCell(1)
        vCG.InsertCellPoint(iCG)

        pdCG = vtk.vtkPolyData()
        pdCG.SetPoints(pCG)
        pdCG.SetVerts(vCG)

        self.viewer.add_polydata(pdCG, color=[1, 0, 0])

        # Ploting also a line between O and CG
        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(self.gravity_center)

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, 0)
        line.GetPointIds().SetId(1, 1)

        lines = vtk.vtkCellArray()
        lines.InsertNextCell(line)

        lines_pd = vtk.vtkPolyData()
        lines_pd.SetPoints(points)
        lines_pd.SetLines(lines)

        self.viewer.add_polydata(lines_pd, color=[0, 0, 0])
        
        # Adding the buoyancy center
        pB = vtk.vtkPoints()
        vB = vtk.vtkCellArray()

        iB = pB.InsertNextPoint(self.buoyancy_center)
        vB.InsertNextCell(1)
        vB.InsertCellPoint(iB)

        pdB = vtk.vtkPolyData()
        pdB.SetPoints(pB)
        pdB.SetVerts(vB)

        self.viewer.add_polydata(pdB, color=[0, 1, 0])

        # Ploting also a line between O and B
        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(self.buoyancy_center)

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, 0)
        line.GetPointIds().SetId(1, 1)

        lines = vtk.vtkCellArray()
        lines.InsertNextCell(line)

        lines_pd = vtk.vtkPolyData()
        lines_pd.SetPoints(points)
        lines_pd.SetLines(lines)

        self.viewer.add_polydata(lines_pd, color=[0, 0, 0])
        
        # Adding corner annotation
        ca = vtk.vtkCornerAnnotation()
        ca.SetLinearFontScaleFactor(2)
        ca.SetNonlinearFontScaleFactor(1)
        ca.SetMaximumFontSize(20)
        labels = "Origin: Black; Gravity Center: Red; Buoyancy Center: Green\nTo see points, press 'w'"
        ca.SetText(2, labels)
        ca.GetTextProperty().SetColor(0., 0., 0.)
        self.viewer.renderer.AddViewProp(ca)
        
        # Showing the viewer
        self.viewer.show()
        self.viewer.finalize()
        
        return


if __name__ == '__main__':

    # The following code are only for testing purpose
    import mmio
    import sys
    import mesh
    
    vertices, faces = mmio.load_VTP('meshmagick/tests/data/SEAREV.vtp')
    
    searev = mesh.Mesh(vertices, faces, name='SEAREV')
    
    hs = Hydrostatics(searev)
    hs.verbose_on()
    
    # hs.mass = 2300
    hs.gravity_center = [0, 0, -2]
    hs.allow_unstable_off()
    hs.reltol = 1e-3
    
    F = [0, 0, -50e3*9.81]
    force1 = AdditionalForce(point=[0, 0, 5], value=[0, 1e7, 0], mode='relative')
    # force2 = AdditionalForce(point=[0, 15, -5], value=[0, 0, -50e3*9.81], mode='absolute')
    # force = AdditionalForce(point=[0, 0, 3], value=[0, 1e7, 0], mode='absolute')
    hs.add_force(force1)
    # hs.add_force(force2)
    
    
    hs.equilibrate()
    hs.show()
    
    print hs.get_hydrostatic_report()
    
    print 'Done !'
    
    

