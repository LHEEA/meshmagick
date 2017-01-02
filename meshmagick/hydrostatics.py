#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""This module allows to perform hydrostatics computations on meshes"""

import numpy as np
import math

from mesh_clipper import MeshClipper

__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-2015, Ecole Centrale de Nantes"
__credits__ = "Francois Rongere"
__licence__ = "CeCILL"
__version__ = "1.0"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@ec-nantes.fr"
__status__ = "Development"

# Data for DNV standards
GM_MIN = 0.15


# TODO: pass this code to a version that only rotates the free surface plane and not the mesh...
# TODO: throw the transformation needed to get the equilibrium from initial position after a successful equilibrium computation
# TODO: make the mesh not to diverge from principal axis


class Force(object):
    """Class to handle a linear force (resultant).
    
    Parameters
    ----------
    point : array_like, optional
        The force application point. Default is (0, 0, 0)
    value : array_like, optional
        The value of the force. Default is (0, 0, 0)
    name : str, optional
        The name of the force
    mode : str, ['relative', 'absolute']
        The mode of the force whether its direction follows the body motion ('relative') or remains fixed with
        respect to the reference frame ('absolute'). Default is 'relative'.
    
    """
    def __init__(self, point=(0, 0, 0), value=(0, 0, 0), name='', mode='relative'):
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
        """Updates the vertical position and orientation of the force.
        
        Parameters
        ----------
        dz : float, optional
            Update in the vertical position of the application point. Default is 0.
        rot : ndarray
            The 3x3 rotation matrix for the force direction. Default is identity.
        """
        
        self.point[2] += dz
        self.point = np.dot(rot, self.point)
        if self.mode == 'relative':
            self.value = np.dot(rot, self.value)

    def update_xy(self, dx, dy):
        """Updates the horizontal position of the application point of the force.
        
        Parameters
        ----------
        dx : float, optional
            Update in the x position of the application point. Defaults is 0.
        dy : float, optional
            Update in the y position of the application point. Defaults is 0.
        """
        
        self.point[0] += dx
        self.point[1] += dy
        return

    @property
    def hs_force(self):
        """Returns the relevant force vector for hydrostatics computations.
        
        The return array's components are the vertical force, x moment and y moment around the application point.
        
        Returns
        -------
        ndarray
        """
        fz = self.value[2]
        moment = np.cross(self.point, self.value)
        mx, my = moment[:2]
        return np.array([fz, mx, my], dtype=np.float)

    def reset(self):  # TODO: a appeler depuis le reset de Hydrostatics
        raise NotImplementedError


class Hydrostatics(object):
    # TODO: refactor this docstring
    """Class to perform hydrostatic computations on meshes.
    
    Parameters
    ----------
    working_mesh : Mesh
        The mesh we want to work with for hydrostatics computations
    cog : array_like, optional
        The mesh's center of gravity coordinates. Default is the origin (0, 0, 0)
    mass : float, optional
        The mesh's mass. Default is 0. Unit is Ton.
    rho_water : float, optional
        The density of water (in kg/m**3). Default is that of salt water (1023 kg//m**3)
    grav : float, optional
        The acceleration of gravity. Default is 9.81 m/s**2.
    
    
    Warnings
    --------
    * At class instantiation, no modification on the mesh position is done and hydrostatic properties are accessible
      directly for this position. The mesh may not be at equilibrium thought and residual force balance between
      hydrostatic force and gravity force may be not null. To trig the equilibrium computations, we have to call the
      equilibrate() method.
    
    * Note that it is useless to call any method to get the hydrostatic properties up to date as it is done internally
      and automatically at each modification.
      
    * Mass unit is the ton !
    """

    def __init__(self, working_mesh, cog=[0., 0., 0.], mass=None, rho_water=1023, grav=9.81,
                 verbose=False, animate=False):

        self.backup = dict()

        self.backup['init_mesh'] = working_mesh.copy()
        self.mesh = working_mesh.copy()

        cog = np.array(cog, dtype=np.float)
        assert cog.shape[0] == 3
        self.backup['gravity_center'] = cog.copy()
        self._gravity_center = cog.copy()
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
            self.mass = mass * 1000.  # conversion in kg
        else:
            self._mass = self.hs_data['disp_mass']

        self._mg = self._mass * self._gravity

        self._verbose = verbose

        self.animate = animate

        self._rotation = np.eye(3, dtype=np.float)

        self.additional_forces = []

    @property
    def verbose(self):
        """Get the verbosity"""
        return self._verbose

    def verbose_on(self):
        """Switch ON the verbosity."""
        self._verbose = True

    def verbose_off(self):
        """Switch OFF the verbosity."""
        self._verbose = False

    @property
    def gravity(self):
        """Get the gravity acceleration"""
        return self._gravity

    @gravity.setter
    def gravity(self, value):
        """Set the gravity acceleration"""
        self._gravity = float(value)
        self._rhog = self._rho_water * value
        self._mg = self._mass * value
        self._update_hydrostatic_properties()

    @property
    def rho_water(self):
        """Get the water density"""
        return self._rho_water

    @rho_water.setter
    def rho_water(self, value):
        """Set the water density"""
        self._rho_water = float(value)
        self._rhog = value * self._gravity
        self._update_hydrostatic_properties()

    @property
    def mass(self):
        """Get the mass (in tons)"""
        return self._mass / 1000.

    @mass.setter
    def mass(self, value):
        """Set the mass. It must be given in tons."""
        self._mass = float(value) * 1000.  # Conversion from tons to kg
        self._mg = self._mass * self._gravity  # SI units

        if self.is_sinking():
            raise ValueError('%s is sinking as it is too heavy.' % self.mesh.name)

    def _max_displacement(self):
        return self._rho_water * self.mesh._compute_volume()  # in kg

    def is_sinking(self):
        """Returns whether the mesh is sinking with the current mass.
        
        Returns
        -------
        bool
        """
        if self._mass > self._max_displacement():
            return True
        else:
            return False

    @property
    def gravity_center(self):
        """Get the gravity center position
        
        Returns
        -------
        ndarray
        """
        return self._gravity_center

    @gravity_center.setter
    def gravity_center(self, value):
        """Set the gravity center position"""
        value = np.asarray(value, dtype=np.float)
        self.backup['gravity_center'] = value
        assert value.shape[0] == 3
        self._gravity_center = value
        self._update_hydrostatic_properties()

    @property
    def zg(self):
        """Get the gravity center vertical position"""
        return self._gravity_center[-1]

    @zg.setter
    def zg(self, value):
        """Set the gravity center vertical position"""
        self._gravity_center[-1] = float(value)
        self._update_hydrostatic_properties()

    @property
    def wet_surface_area(self):
        """Get the area of the underwater surface (m**2)"""
        return self.hs_data['wet_surface_area']

    @property
    def displacement_volume(self):
        """Get the volume displacement (m**2)"""
        return self.hs_data['disp_volume']

    @property
    def displacement(self):
        """Get the mass displacement (tons)"""
        return self.hs_data['disp_mass'] / 1000.

    @property
    def buoyancy_center(self):
        """Get the position of the buoyancy center"""
        return self.hs_data['buoy_center']

    @property
    def flotation_surface_area(self):
        """Get the area of the flotation plane"""
        return self.hs_data['waterplane_area']

    @property
    def flotation_center(self):
        """Get the position of the center of the flotation plane"""
        return self.hs_data['flotation_center']

    @property
    def transversal_metacentric_radius(self):
        """Get the transversal metacentric radius"""
        return self.hs_data['transversal_metacentric_radius']

    @property
    def longitudinal_metacentric_radius(self):
        """Get the longitudinal metacentric radius"""
        return self.hs_data['longitudinal_metacentric_radius']

    @property
    def transversal_metacentric_height(self):
        """Get the transversal metacentric height (GMx)"""
        return self.hs_data['gm_x']

    @property
    def longitudinal_metacentric_height(self):
        """Get the longitudinal metacentric height (GMy)"""
        return self.hs_data['gm_y']

    @property
    def hydrostatic_stiffness_matrix(self):
        """Get the hydrostatic stiffness matrix"""
        return self.hs_data['stiffness_matrix']

    @property
    def hydrostatic_mesh(self):
        """Get the underwater part of the mesh"""
        return self.hs_data['clipper'].clipped_mesh

    def reset(self):
        """Reset hydrostatics with respect to the initial mesh"""
        # TODO: Utiliser plutot la rotation generale pour retrouver le maillage initial

        self.mesh = self.backup['init_mesh'].copy()
        self._gravity_center = self.backup['gravity_center'].copy()
        self._reinit_clipper()

        self._update_hydrostatic_properties()

        self._rotation = np.eye(3, dtype=np.float)

    def is_stable_in_roll(self):
        """Returns whether the mesh is stable in roll (GMx positive)
        
        Returns
        -------
        bool
        """
        return self.transversal_metacentric_height > 0.

    def is_stable_in_pitch(self):
        """Returns whether the mesh is stable in pitch (GMy positive)

        Returns
        -------
        bool
        """
        return self.longitudinal_metacentric_height > 0.

    def isstable(self):
        """Returns whether the mesh is both stable in roll and in pitch

        Returns
        -------
        bool
        """
        return self.is_stable_in_pitch() and self.is_stable_in_roll()

    def is_at_equilibrium(self):
        """Returns whether the mesh is actually in an equilibrium configuration
        
        Returns
        -------
        bool
        """
        residual = self.residual

        mg = self._mg
        breadth = self.hs_data['breadth']
        lwl = self.hs_data['lwl']
        scale = np.array([mg, mg * breadth, mg * lwl])

        if np.all(np.fabs(residual / scale) < self._solver_parameters['reltol']):
            return True
        else:
            return False

    @property
    def delta_fz(self):
        """The residual vertical force
        
        A non-zero value indicates that the mesh is not at hydrostatic equilibrium
        
        Returns
        -------
        float
        """
        fz, _, _ = self.residual
        return fz

    @property
    def delta_mx(self):
        """The residual moment around x
        
        A non-zero value indicates that the mesh is not at hydrostatic equilibrium
        
        Returns
        -------
        float
        """
        _, mx, _ = self.residual
        return mx

    @property
    def delta_my(self):
        """The residual moment around x

        A non-zero value indicates that the mesh is not at hydrostatic equilibrium
        
        Returns
        -------
        float
        """
        _, _, my = self.residual
        return my

    # TODO: create a Hydrostatic stiffness matrix class which should be a symmetric array class
    @property
    def S33(self):
        """The S33 heave-heave stiffness coefficient"""
        return self.hs_data['stiffness_matrix'][0, 0]

    @property
    def S34(self):
        """The S34 heave-roll stiffness coefficient"""
        return self.hs_data['stiffness_matrix'][0, 1]

    @property
    def S35(self):
        """The S35 heave-pitch stiffness coefficient"""
        return self.hs_data['stiffness_matrix'][0, 2]

    @property
    def S44(self):
        """The S44 roll-roll stiffness coefficient"""
        return self.hs_data['stiffness_matrix'][1, 1]

    @property
    def S45(self):
        """The S45 roll-pitch stiffness coefficient"""
        return self.hs_data['stiffness_matrix'][1, 2]

    @property
    def S55(self):
        """The S55 pitch-pitch stiffness coefficient"""
        return self.hs_data['stiffness_matrix'][2, 2]

    @property
    def reltol(self):
        """The relative tolerance for hydrostatic equilibrium solver"""
        return self._solver_parameters['reltol']

    @reltol.setter
    def reltol(self, value):
        """Set the relative tolerance for hydrostatic equilibrium solver"""
        value = float(value)
        assert 0. < value < 1.
        self._solver_parameters['reltol'] = value

    @property
    def theta_relax(self):
        """The angle relaxation value of the hydrostatic solver"""
        return self._solver_parameters['theta_relax']

    @theta_relax.setter
    def theta_relax(self, value):
        """Set the angle relaxation value of the hydrostatic solver"""
        value = float(value)
        assert 0. < value < 10.
        self._solver_parameters['theta_relax'] = value

    @property
    def z_relax(self):
        """The distance relaxation value of the hydrostatic solver"""
        return self._solver_parameters['z_relax']

    @z_relax.setter
    def z_relax(self, value):
        """Set the distance relaxation value of the hydrostatic solver"""
        value = float(value)
        assert 0. < value < 1.
        self._solver_parameters['z_relax'] = value

    @property
    def max_iterations(self):
        """The maximum number of iterations to find a hydrostatic equilibrium starting from an initial configuration"""
        return self._solver_parameters['itermax']

    @max_iterations.setter
    def max_iterations(self, value):
        """Set the maximum number of iterations to find a hydrostatic equilibrium starting from an initial
        configuration"""
        value = int(value)
        assert value > 0.
        self._solver_parameters['itermax'] = value

    @property
    def max_restart(self):
        """The maximum number of random restart for hydrostatic equilibrium computations"""
        return self._solver_parameters['max_nb_restart']

    @max_restart.setter
    def max_restart(self, value):
        """Set the maximum number of random restart for hydrostatic equilibrium computations"""
        value = int(value)
        assert value >= 0.
        self._solver_parameters['max_nb_restart'] = value

    @property
    def allow_unstable(self):
        """Whether unstable equilibrium configurations allowed"""
        return self._solver_parameters['stop_at_unstable']

    def allow_unstable_on(self):
        """Switches ON the allowance of unstable equilibrium configurations finding"""
        self._solver_parameters['stop_at_unstable'] = True

    def allow_unstable_off(self):
        """Switches OFF the allowance of unstable equilibrium configurations finding"""
        self._solver_parameters['stop_at_unstable'] = False
        return

    def _reinit_clipper(self):
        try:
            del self.hs_data['clipper']
        except KeyError:
            pass
        return

    def _update_hydrostatic_properties(self):
        """Updates the hydrostatics properties of the mesh.
        """

        eps = 1e-4  # For zeroing tiny coefficients in the hydrostatic stiffness matrix

        # Clipping the mesh by the Oxy plane
        # TODO: ne pas recreer une instance de MeshClipper a fois mais la stocker et mettre a jour son maillage ainsi
        #  que ses donnees

        # FIXME: on ne devrait recouper le maillage que si ce dernier a ete modifie !!!
        try:
            clipper = self.hs_data['clipper']
            clipped_mesh = clipper.clipped_mesh
        except KeyError:
            clipper = MeshClipper(self.mesh, assert_closed_boundaries=True, verbose=False)
            self.hs_data['clipper'] = clipper
            clipped_mesh = clipper.clipped_mesh

        # Retrieving faces properties for the clipped mesh
        areas = clipped_mesh.faces_areas
        # normals = clipped_mesh.faces_normals
        # centers = clipped_mesh.faces_centers

        if np.any(np.fabs(areas) < 1e-10):  # TODO: bloc a retirer
            print 'probleme de facette'
            self.mesh.quick_save()
            raise Exception

        wet_surface_area = areas.sum()

        inertia = clipped_mesh.eval_plain_mesh_inertias(rho_medium=self.rho_water)
        xb, yb, zb = inertia.gravity_center
        disp_volume = inertia.mass / self.rho_water

        # Computing quantities from intersection polygons
        sigma0 = 0.  # \iint_{waterplane_area} dS = waterplane_area
        sigma1 = 0.  # \iint_{waterplane_area} x dS
        sigma2 = 0.  # \iint_{waterplane_area} y dS
        sigma3 = 0.  # \iint_{waterplane_area} xy dS
        sigma4 = 0.  # \iint_{waterplane_area} x^2 dS
        sigma5 = 0.  # \iint_{waterplane_area} y^2 dS

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
        waterplane_area = sigma0

        # Stiffness matrix coefficients that do not depend on the position of the gravity center
        rhog = self._rhog
        s33 = rhog * waterplane_area
        s34 = rhog * sigma2
        s35 = -rhog * sigma1
        s45 = -rhog * sigma3

        # Metacentric radius (Bouguer formulae)
        transversal_metacentric_radius = sigma5 / disp_volume  # Around Ox
        longitudinal_metacentric_radius = sigma4 / disp_volume  # Around Oy

        # Metacentric height
        a = self.zg - zb  # BG
        gm_x = transversal_metacentric_radius - a
        gm_y = longitudinal_metacentric_radius - a

        # Stiffness matrix coefficients that depend on the position of the gravity center
        s44 = rhog * disp_volume * gm_x
        s55 = rhog * disp_volume * gm_y

        # Assembling stiffness matrix
        stiffness_matrix = np.array([[s33, s34, s35],
                                     [s34, s44, s45],
                                     [s35, s45, s55]], dtype=np.float)

        # Zeroing tiny coefficients
        stiffness_matrix[np.fabs(stiffness_matrix) < eps] = 0.

        # Flotation center F:
        x_f = -s35 / s33
        y_f = s34 / s33

        xmin, xmax, ymin, ymax, zmin, zmax = clipped_mesh.axis_aligned_bbox

        # Storing data
        self.hs_data['wet_surface_area'] = wet_surface_area
        self.hs_data['disp_volume'] = disp_volume
        self.hs_data['disp_mass'] = self._rho_water * disp_volume
        self.hs_data['buoy_center'] = np.array([xb, yb, zb], dtype=np.float)
        self.hs_data['flotation_center'] = np.array([x_f, y_f, 0.], dtype=np.float)
        self.hs_data['waterplane_area'] = waterplane_area
        self.hs_data['transversal_metacentric_radius'] = transversal_metacentric_radius
        self.hs_data['longitudinal_metacentric_radius'] = longitudinal_metacentric_radius
        self.hs_data['gm_x'] = gm_x
        self.hs_data['gm_y'] = gm_y
        self.hs_data['stiffness_matrix'] = stiffness_matrix
        self.hs_data['lwl'] = maxx - minx
        self.hs_data['los'] = xmax - xmin
        self.hs_data['bos'] = ymax - ymin
        self.hs_data['draught'] = math.fabs(zmin)
        self.hs_data['fp'] = maxx
        self.hs_data['breadth'] = maxy - miny

        # TODO: we should better store the inertia object !
        inertia.shift_at_cog()
        self.hs_data['Ixx'] = inertia.xx
        self.hs_data['Iyy'] = inertia.yy
        self.hs_data['Izz'] = inertia.zz
        self.hs_data['Ixy'] = inertia.xy
        self.hs_data['Ixz'] = inertia.xz
        self.hs_data['Iyz'] = inertia.yz

        return

    def get_gravity_force(self):
        """Returns the gravity force applied on the body
        
        Returns
        -------
        Force
        """
        return Force(point=self._gravity_center, value=[0, 0, -self._mass * self._gravity], mode='absolute')

    def get_buoyancy_force(self):
        """Returns the buoyancy force applied on the body
        
        Returns
        -------
        Force
        """
        value = [0, 0, self.rho_water * self._gravity * self.displacement_volume]
        return Force(point=self.buoyancy_center, value=value, mode='absolute')

    def add_force(self, force):
        """Add a custom force applying on the body
        
        Parameters
        ----------
        force : Force
            External force
        """
        # TODO: allow to remove those forces...
        assert isinstance(force, Force)
        self.additional_forces.append(force)

    @property
    def _scale(self):
        """The scaling factor for relative convergence"""
        mg = self._mg
        breadth = self.hs_data['breadth']
        lwl = self.hs_data['lwl']
        scale = np.array([mg, mg * breadth, mg * lwl])
        return scale

    @property
    def residual(self):
        """The residual force resulting from the balance between gravity, buoyancy and external forces
        
        It is composed of the vertical force and horizontal moments only.
        
        Returns
        -------
        ndarray
        """
        rhog_v = self._rhog * self.hs_data['disp_volume']
        mg = self._mg

        xb, yb, _ = self.hs_data['buoy_center']
        xg, yg, _ = self._gravity_center

        residual = np.array([rhog_v - mg,
                             rhog_v * yb - mg * yg,
                             -rhog_v * xb + mg * xg])

        # Accounting for additional forces in the static equilibrium
        for force in self.additional_forces:
            residual += force.hs_force

        return residual

    def set_displacement(self, disp):
        """
        Displaces the mesh at a prescribed displacement

        Parameters
        ----------
        disp : float
            Mass displacement of the hull (in tons)
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

            self._reinit_clipper()
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

    def equilibrate(self, init_disp=True):
        """Performs 3D equilibrium search.
        
        Parameters
        ----------
        init_disp : bool, optional
            Flag to indicate if the mesh has to be first placed at its displacement. Default is True.
        Returns
        -------
        int
            A code indicating the state of the solver at the end of the computations
            
        Notes
        -----
        
        The return code can have the following values:
        
        * 0 : Failed to find an equilibrium position
        * 1 : A stable equilibrium configuration has been reached
        * 2 : An unstable equilibrium configuration has been reached
        """
        
        linear_solver = np.linalg.solve
        
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
            if unstable_config or (iter == (nb_restart + 1) * (itermax - 1)):
                if self.verbose:
                    if unstable_config:
                        print 'Unstable equilibrium reached.'
                        print '\t-> Keep searching a stable configuration by random restart number %u.' % (
                        nb_restart + 1)
                    else:
                        print 'Failed to find an equilibrium configuration with these initial conditions in %u ' \
                              'iterations.' % itermax
                        print '\t-> Keep searching by random restart number %u.' % (nb_restart + 1)

                unstable_config = False

                # Max iterations reached
                if nb_restart < max_nb_restart - 1:
                    nb_restart += 1
                    # Random on the position of the body
                    thetax, thetay = np.random.rand(2) * math.pi
                    rot_matrix = self.mesh.rotate([thetax, thetay, 0.])
                    self._gravity_center = np.dot(rot_matrix, self._gravity_center)
                    self._rotation = np.dot(rot_matrix, self._rotation)
                    self._reinit_clipper()
                    self._update_hydrostatic_properties()

                    for force in self.additional_forces:
                        force.update(rot=rot_matrix)

                    dz = thetax = thetay = 0.

                else:
                    # Max number of restart allowed. Failed to find an equilibrium configuration.
                    code = 0
                    break

            # Applying transformation to the mesh
            self.mesh.translate_z(dz)
            self._gravity_center[2] += dz

            rot_matrix = self.mesh.rotate([thetax, thetay, 0.])
            self._gravity_center = np.dot(rot_matrix, self._gravity_center)
            self._rotation = np.dot(rot_matrix, self._rotation)

            # Updating force data
            for force in self.additional_forces:
                force.update(dz=dz, rot=rot_matrix)

            self._reinit_clipper()
            self._update_hydrostatic_properties()

            # TODO: animation may be trigged here

            residual = self.residual
            scale = self._scale

            if np.all(np.fabs(self.residual / scale) < reltol):
                # Convergence at an equilibrium
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
            stiffness_matrix = self.hs_data['stiffness_matrix']
            dz, thetax, thetay = linear_solver(stiffness_matrix, residual)

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

        for force in self.additional_forces:
            force.update_xy(-self._gravity_center[0], -self._gravity_center[1])

        self._gravity_center[:2] = 0.

        if self.verbose:
            if code == 0:
                print "\t-> Maximum number of restart reached. Failed to find an equilibrum position."
            elif code == 1:
                print "Stable equilibrium reached after %u iterations and %u random restart" % (iter, nb_restart)
            elif code == 2:
                print 'Unstable equilibrium reached after %u iterations and %u random restart' % (iter, nb_restart)

        return code

    def get_hydrostatic_report(self):
        """Returns a hydrostatic report for the current configuration
        
        Returns
        -------
        str
        """
        
        def hspace():
            return '\n'

        def build_line(text, data, precision=3, dtype='f'):
            # TODO: ajouter unit
            textwidth = 40
            try:
                line = '\t{:-<{textwidth}}>  {:< .{precision}{dtype}}\n'.format(str(text).upper(), data,
                                                                                precision=precision,
                                                                                textwidth=textwidth,
                                                                                dtype=dtype
                                                                                )
            except ValueError:
                if isinstance(data, np.ndarray):
                    if data.ndim == 1:
                        data_str = ''.join(['{:< 10.{precision}{dtype}}'.format(val, precision=precision, dtype=dtype)
                                            for val in data])
                        line = '\t{:-<{textwidth}}>  {}\n'.format(str(text).upper(), data_str, textwidth=textwidth)

                        # else:
                        #     print data

            return line

        # TODO: ajouter la reference au point de calcul de la matrice raideur

        msg = '\n'
        # title = 'Hydrostatic report ({0})\n
        # \tGenerated by meshmagick on {1}'.format(self.mesh.name, strftime('%c')).upper()
        # msg += header(title)

        msg += hspace()
        msg += build_line('Gravity acceleration (M/S**2)', self.gravity, precision=2)
        msg += build_line('Density of water (kg/M**3)', self.rho_water, precision=1)

        msg += hspace()
        msg += build_line('Waterplane area (M**2)', self.flotation_surface_area, precision=1)
        msg += build_line('Waterplane center (M)', self.flotation_center[:2], precision=3)
        msg += build_line('Wet area (M**2)', self.wet_surface_area, precision=1)
        msg += build_line('Displacement volume (M**3)', self.displacement_volume, precision=3)
        msg += build_line('Displacement mass (tons)', self.displacement, precision=3)
        msg += build_line('Buoyancy center (M)', self.buoyancy_center, precision=3)
        msg += build_line('Center of gravity (M)', self.gravity_center, precision=3)

        msg += hspace()
        msg += build_line('Draught (M)', self.hs_data['draught'], precision=3)  # TODO
        msg += build_line('Length overall submerged (M)', self.hs_data['los'], precision=2)
        msg += build_line('Breadth overall submerged (M)', self.hs_data['bos'], precision=2)
        msg += build_line('Length at Waterline LWL (M)', self.hs_data['lwl'], precision=2)
        msg += build_line('Forward perpendicular (M)', self.hs_data['fp'], precision=2)

        msg += hspace()
        msg += build_line('Transversal metacentric radius (M)', self.transversal_metacentric_radius, precision=3)
        msg += build_line('Transversal metacentric height GMt (M)', self.transversal_metacentric_height, precision=3)
        msg += build_line('Longitudinal metacentric radius (M)', self.longitudinal_metacentric_radius, precision=3)
        msg += build_line('Longitudinal metacentric height GMl (M)', self.longitudinal_metacentric_height, precision=3)

        msg += hspace()
        msg += '\tHYDROSTATIC STIFFNESS COEFFICIENTS:\n'
        msg += build_line('K33 (N/M)', self.S33, precision=4, dtype='E')
        msg += build_line('K34 (N)', self.S34, precision=4, dtype='E')
        msg += build_line('K35 (N)', self.S35, precision=4, dtype='E')
        msg += build_line('K44 (N.M)', self.S44, precision=4, dtype='E')
        msg += build_line('K45 (N.M)', self.S45, precision=4, dtype='E')
        msg += build_line('K55 (N.M)', self.S55, precision=4, dtype='E')

        # Il faut faire une correction avec le plan de la flottaison de certains coeffs
        msg += hspace()
        msg += '\tINERTIAS:\n'
        msg += build_line('Ixx', self.hs_data['Ixx'], precision=3, dtype='E')
        msg += build_line('Ixy', self.hs_data['Ixy'], precision=3, dtype='E')
        msg += build_line('Ixz', self.hs_data['Ixz'], precision=3, dtype='E')
        msg += build_line('Iyy', self.hs_data['Iyy'], precision=3, dtype='E')
        msg += build_line('Iyz', self.hs_data['Iyz'], precision=3, dtype='E')
        msg += build_line('Izz', self.hs_data['Izz'], precision=3, dtype='E')

        msg += hspace()
        msg += '\tRESIDUALS:\n'
        msg += build_line('Absolute', self.residual, precision=3, dtype='E')
        msg += build_line('Relative', self.residual / self._scale, precision=3, dtype='E')
        msg += build_line('Relative tolerance', self.reltol, precision=1, dtype='E')
        # residual = self.residual
        # msg += ('\nResidual:\n')
        # msg += ('Delta Fz = %.3f N\n' % residual[0])
        # msg += ('Delta Mx = %.3f Nm\n' % residual[1])
        # msg += ('Delta My = %.3f Nm\n' % residual[2])
        #
        # rel_res = residual / self._scale
        # msg += ('\nRelative residual:\n')
        # msg += ('Delta Fz = %E\n' % rel_res[0])
        # msg += ('Delta Mx = %E\n' % rel_res[1])
        # msg += ('Delta My = %E\n' % rel_res[2])
        # msg += ('Relative tolerance of the solver: %.1E\n' % self.reltol)

        return msg

    # FIXME: la methode show ne devrait pas faire appel explicitement a des fonctions vtk...
    # Tout devrait etre gere dans MMViewer
    def show(self):
        """Displays the mesh in the meshmagick viewer with additional graphics proper to hydrostatics"""
        # TODO: Ce n'est pas ce module qui doit savoir utiliser vtk !!!

        import MMviewer
        import vtk

        vtk_polydata = self.mesh._vtk_polydata()
        self.viewer = MMviewer.MMViewer()
        self.viewer.add_polydata(vtk_polydata)

        pd_orig = self.viewer.add_point([0, 0, 0], color=[0, 0, 0])

        pd_cog = self.viewer.add_point(self.gravity_center, color=[1, 0, 0])

        pd_orig_cog = self.viewer.add_line([0, 0, 0], self.gravity_center, color=[1, 0, 0])

        pd_buoy = self.viewer.add_point(self.buoyancy_center, color=[0, 1, 0])

        pd_orig_buoy = self.viewer.add_line([0, 0, 0], self.buoyancy_center, color=[0, 1, 0])

        scale = self.mass * 1000

        # Adding glyph for gravity force
        gforce = self.get_gravity_force()
        self.viewer.add_vector(gforce.point, gforce.value, scale=scale, color=[1, 0, 0])

        # Adding glyph for buoyancy force
        bforce = self.get_buoyancy_force()
        self.viewer.add_vector(bforce.point, bforce.value, scale=scale, color=[0, 1, 0])

        # Adding glyph for additional forces
        for force in self.additional_forces:
            self.viewer.add_point(force.point, color=[0, 0, 1])
            self.viewer.add_vector(force.point, force.value, scale=scale, color=[0, 0, 1])

        # TODO: mettre la suite dans une methode de MMViewer
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
