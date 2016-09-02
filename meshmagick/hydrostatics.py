
__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-2015, Ecole Centrale de Nantes"
__credits__ = "Francois Rongere"
__licence__ = "CeCILL"
__version__ = "1.0"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@ec-nantes.fr"
__status__ = "Development"

# import meshmagick as mm
import mesh
import mmio
from mesh_clipper import MeshClipper

import numpy as np
import math


class HSmesh(mesh.Mesh):
    """
    Class to perform hydrostatic computations on meshes.
    
    Rules:
    ------
    *** At class instantiation, no modification on the mesh position is done and hydrostatic properties are accessible directly for this position. The mesh may not be at equilibrium thought and residual force balance between hydrostatic force and gravity force may be not null. To trig the equilibrium computations, we have to call the equilibrate() method.
    
       Note that it is useless to call any method to get the hydrostatic properties up to date as it is done internally and automatically at each modification.
    
    *** Intrinsic mesh properties are
    
        - mass
        - cog
        - zg
        - rho_water
        - gravity
        
        Setting one of these property by its setter will trig a search for a new mesh position that fullfill the hydrostatic equilibrium.
        Note that zg is mandatory to get the hydrostatic stiffness matrix. By default, be carefull that it is set to zero at instantiation.
    
    *** Modification of the displacement will produce the mesh to move along z so that
    it reaches the correct displacement. The mass is then set accordingly and equals the displacement. Note that this is different from modifying the mass and
    
    
    
    
    
    
    *** Setting the buoyancy center will trig a search for a position of the gravity center that fullfill this hydrostatic property, while keeping the correct displacement. NOT IMPLEMENTED, THIS IS AN ADVANCED FEATURE
    
    """
    def __init__(self, vertices, faces, name=None,
                 CG=np.zeros(3, dtype=np.float), zg=0., mass=None, rho_water=1023, grav=9.81,
                 animate=False, mass_unit='tons'):
        
        super(HSmesh, self).__init__(vertices, faces, name=None)
        self.__internals__['backup'] = dict()
        self.__internals__['backup']['initial_vertices'] = vertices
        self.__internals__['backup']['initial_faces'] = faces
        
        assert mass_unit in ['kg', 'tons']
        
        self.__internals__['mass_unit'] = mass_unit # TODO: l'utiliser...
        
        self._gravity_center = CG
        self._rho_water = rho_water
        self._gravity = grav
        
        self.__internals__['rhog'] = rho_water*grav
        
        self.__internals__['hs_data'] = dict()
        self.__internals__['hs_data']['Sw'] = 0.
        self.__internals__['hs_data']['Vw'] = 0.
        self.__internals__['hs_data']['disp'] = 0.
        self.__internals__['hs_data']['B'] = np.zeros(3, dtype=np.float)
        self.__internals__['hs_data']['Sf'] = 0.
        self.__internals__['hs_data']['F'] = np.zeros(3, dtype=np.float)
        self.__internals__['hs_data']['r'] = 0.
        self.__internals__['hs_data']['R'] = 0.
        self.__internals__['hs_data']['GMx'] = 0.
        self.__internals__['hs_data']['GMy'] = 0.
        self.__internals__['hs_data']['KH'] = np.zeros((3, 3), dtype=np.float)
        
        # TODO: ajouter le calcul du tirant d'eau et d'air
        
        self._update_hydrostatic_properties()
        
        if mass:
            self._mass = mass
        else:
            self._mass = self.__internals__['hs_data']['disp']
    
    @property
    def gravity(self):
        return self._gravity
    
    @gravity.setter
    def gravity(self, value):
        self._gravity = value
        self._update_hydrostatic_properties()
        return
    
    @property
    def rho_water(self):
        return self._rho_water
    
    @rho_water.setter
    def rho_water(self, value):
        self._rho_water = value
        self._update_hydrostatic_properties()
        return
    
    @property
    def zg(self):
        return self._gravity_center[-1]
    
    @zg.setter
    def zg(self, value):
        self._gravity_center[-1] = value
        self._update_hydrostatic_properties()
        return
    
    @property
    def mass_unit(self):
        return self.__internals__['mass_unit']
    
    def mass_in_kg(self):
        self.__internals__['mass_unit'] = 'kg'
    
    def mass_in_tons(self):
        self.__internals__['mass_unit'] = 'tons'
    
    
    
    @property
    def wetted_surface_area(self):
        return self.__internals__['hs_data']['Sw']
    
    @property
    def immersed_volume(self):
        return self.__internals__['hs_data']['Vw']
    
    @property
    def displacement(self):
        if self.mass_unit == 'tons':
            return self.__internals__['hs_data']['disp'] / 1000.
        else:
            return self.__internals__['hs_data']['disp']
    
    @property
    def buoyancy_center(self):
        return self.__internals__['hs_data']['B']
    
    @property
    def flotation_surface_area(self):
        return self.__internals__['hs_data']['Sf']
    
    @property
    def flotation_center(self):
        return self.__internals__['hs_data']['F']
    
    @property
    def transversal_metacentric_radius(self):
        return self.__internals__['hs_data']['r']
    
    @property
    def longitudinal_metacentric_radius(self):
        return self.__internals__['hs_data']['R']
    
    @property
    def transversal_metacentric_height(self):
        return self.__internals__['hs_data']['GMx']
    
    @property
    def longitudinal_metacentric_height(self):
        return self.__internals__['hs_data']['GMy']
    
    @property
    def hydrostatic_stiffness_matrix(self):
        return self.__internals__['hs_data']['KH']
    
    def reset(self):
        vertices = self.__internals__['initial_vertices']
        faces = self.__internals__['initial_faces']
        self.vertices = vertices
        self.faces = faces
        
        self._update_hydrostatic_properties()
        return
    
    def _update_hydrostatic_properties(self): # TODO: voir si on rend la methode privee...
        "Computation of hydrostatics"
        data = compute_hydrostatics(self, self.zg, rho_water=self._rho_water, grav=self._gravity)
        # print hs_output
        self.__internals__['hs_data']['Sw'] = data['Sw']
        self.__internals__['hs_data']['Vw'] = data['Vw']
        self.__internals__['hs_data']['disp'] = data['disp']
        self.__internals__['hs_data']['B'] = data['B']
        self.__internals__['hs_data']['Sf'] = data['Sf']
        self.__internals__['hs_data']['F'] = data['F']
        self.__internals__['hs_data']['r'] = data['r']
        self.__internals__['hs_data']['R'] = data['R']
        self.__internals__['hs_data']['GMx'] = data['GMx']
        self.__internals__['hs_data']['GMy'] = data['GMy']
        self.__internals__['hs_data']['KH'] = data['KH']
        
        return
    
    @property
    def residual(self):
        rhogV = self.__internals__['rhog'] * self.__internals__['hs_data']['Vw']
        xb, yb, _ = self.__internals__['hs_data']['B']
        mg = self._mass * self._gravity
        
        residual = np.array([rhogV - mg,
                             rhogV*yb - mg*yg,
                            -rhogV*xb + mg*xg])
        
        return residual
    
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
        
        













# TODO: voir a mettre ca dans la classe mesh
def show(mymesh, CG=None, B=None):

    # TODO: afficher polygones...
    import MMviewer
    import vtk
    my_viewer = MMviewer.MMViewer()

    # Adding mesh
    polydata_mesh = mm._build_vtkPolyData(vertices, faces)
    my_viewer.add_polydata(polydata_mesh)

    # Flotation plane
    plane = vtk.vtkPlaneSource()

    # plane.SetNormal(0, 0, 1)
    xmin, ymin = vertices[:, :2].min(axis=0)
    xmax, ymax = vertices[:, :2].max(axis=0)
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
    #         points_coord = _vertices[polygon_ids]
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

def compute_equilibrium(mymesh, disp, CG, rho_water=1023, grav=9.81,
                        reltol=1e-2, itermax=100, max_nb_relaunch=10, theta_relax=2, z_relax=0.1,
                        verbose=False, anim=False):
    """

    Parameters
    ----------
    mymesh : Mesh
        mesh to be placed at hydrostatic equilibrium
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
    # max_nb_relaunch = 10

    # Relaxation parameters
    # z_relax = 0.1 # in meters
    theta_relax_x = math.radians(theta_relax)
    theta_relax_y = math.radians(theta_relax)


    rhog = rho_water*grav
    mg = disp*grav*1e3 # Conversion of displacement in kg

    dz = 0.
    thetax = 0.
    thetay = 0.
    iter = 0

    # rot = np.eye(3, 3)
    dz_update = 0.

    # Vc = vertices.copy()
    # Fc = faces.copy()
    mymesh_c = mymesh.copy()
    # Vc = mymesh.vertices.copy()
    # Fc = mymesh.faces.copy()
    CGc = np.asarray(CG, dtype=np.float).copy()

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
                rot = _get_rotation_matrix(thetax, thetay)
                Vc = np.transpose(np.dot(rot, Vc.T))
                CGc = np.dot(rot, CGc)
                # origin = np.dot(rot, origin)
                dz = thetax = thetay = 0.
            else:
                code = 0
                break
        
        # Transformation
        # rot = mesh._rodrigues(thetax, thetay)

        # Applying transformation to the mesh
        mymesh_c.translate_z(dz)
        rot = mymesh_c.rotate([thetax, thetay, 0.])
        # Vc = np.transpose(np.dot(rot, Vc.T)) # rotation

        if anim: # FIXME: IO have been moved to mmio module...
            mmio.write_VTP('mesh%u.vtp'%iter, mymesh_c.vertices, mymesh_c.faces)

        # Applying transformation to the center of gravity
        CGc[-1] += dz
        CGc = np.dot(rot, CGc)
        xg, yg, zg = CGc

        # origin[-1] += dz
        # origin = np.dot(rot, origin)

        # Computing hydrostatics
        hs_output = compute_hydrostatics(mymesh_c, zg, rho_water=rho_water, grav=grav, verbose=False)

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
                #     rot = _get_rotation_matrix(0., math.pi)
                #     Vc = np.transpose(np.dot(rot, Vc.T)) # rotation
                #     CGc = np.dot(rot, CGc)
                #
                # if not xstable:
                #     print '\nConvergence reach at an unstable configuration in ROLL at iteration %u' % iter
                #     print '\t--> Keep going iterations with an opposite orientation\n'
                #     # Choose a new position at 180 deg in roll
                #     rot = _get_rotation_matrix(math.pi, 0.)
                #     Vc = np.transpose(np.dot(rot, Vc.T)) # rotation
                #     CGc = np.dot(rot, CGc)

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
    mymesh_c.translate([-CGc[0], -CGc[1], 0.])
    # Vc[:, 0] -= CGc[0]
    # Vc[:, 1] -= CGc[1]
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

    return mymesh_c, CGc

def _get_Sf_Vw(mymesh):
    """
    Computes only the flotation surface area and the immersed volume of the mesh

    Parameters
    ----------
    vertices : ndarray
        Array of the mesh _vertices
    faces : ndarray
        Array of the mesh connectivities

    Returns
    -------
    Vc : ndarray
    Fc : ndarray
    Sf : float
    Vw : float
    """
    Vc, Fc, clip_infos = mm.clip_by_plane(vertices, faces, mm.Plane(), infos=True)

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

def set_displacement(mymesh, disp, rho_water=1023, grav=9.81, abs_tol= 1., itermax=25, verbose=False):
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

    dz = 0.
    iter = 0
    rho = rho_water*1e-3 # To convert automatically the displacement in kg in calculations

    while True:
        if iter == itermax:
            if verbose:
                print 'No convergence of the displacement after %u iterations' % itermax
            break

        Vc = vertices.copy()
        Fc = faces.copy()

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

def compute_hydrostatics(mymesh, zg, rho_water=1023, grav=9.81, verbose=False):
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

    """

    eps = 1e-4 # For zeroing tiny coefficients in the hydrostatic stiffness matrix

    # TODO: initialiser un clipper

    # Clipping the mesh by the Oxy plane
    plane = mesh.Plane([0., 0., 1.], 0.) # Oxy plane
    # try:
        # Vc, Fc, clip_infos = mm.clip_by_plane(_vertices, _faces, plane, infos=True)

    clipper = MeshClipper(mymesh, plane, assert_closed_boundaries=True, verbose=False)

    clipped_mesh = clipper.clipped_mesh
        # clipped_mesh, boundaries = mymesh.clip(plane, return_boundaries=True, assert_closed_boundaries=True)
    # except:
    #     show(mesh._vertices, mesh._faces)
    #     raise Exception, 'Hydrostatic module only work with watertight hull. Please consider using the --sym option.'

    # TODO: retenir les pptes du maillage initial
    # Calculs des pptes des facettes
    areas = clipped_mesh.faces_areas
    normals = clipped_mesh.faces_normals
    # centers = clipped_mesh.faces_centers
    centers = clipped_mesh.faces_centers

    # areas, normals, centers = mm.get_all_faces_properties(Vc, Fc)

    # Calcul surface mouillee
    Sw = areas.sum()


    # TODO: utiliser des formules analytiques et non approchees comme celles-ci !
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
    # output['Vc'] = Vc
    # output['Fc'] = Fc
    # output['polygons'] = polygons

    return output


if __name__ == '__main__':

    # The following code are only for testing purpose
    import mmio
    # import pickle

    vertices, faces = mmio.load_VTP('meshmagick/tests/data/SEAREV.vtp')
    
    searev = HSmesh(vertices, faces, name='SEAREV', mass_unit='tons')
    
    
    # print searev.immersed_volume
    # print searev.displacement
    # print searev.hydrostatic_stiffness_matrix

    # searev = mesh.Mesh(vertices, faces)
    #
    # searevc, CGc = compute_equilibrium(searev, 1500, [0, 45, 0], verbose=True)
    #
    # searevc.show()
    # hs_output = compute_hydrostatics(searevc, CGc[-1], verbose=True)
    
    
    
    
    # dz_list = [-5, -2, 0, 2, 4]
    #
    # for dz in dz_list:
    #     searev_cp = searev.copy()
    #     searev_cp.translate_z(dz)
    #     hs_output = compute_hydrostatics(searev_cp, 0)
    #
    #     hs_output_old = pickle.load(open('hs_output_%s.p'%dz, 'r'))
    #     print hs_output, hs_output_old
    
    
    

    # vertices, faces = mmio.load_VTP('Cylinder.vtp')
    # cylinder = mesh.Mesh(vertices, faces)
    #
    # # SEAREV :
    # disp = 2000
    # CG = np.array([-1, 0, -3])


    # vertices, faces, CGc = compute_equilibrium(searev, disp, CG, rho_water=1023, grav=9.81,
    #                                 reltol=1e-2, itermax=100, max_nb_relaunch=15, theta_relax=1, z_relax=0.1,
    #                                 verbose=True, anim=False)


    # hs_output = compute_hydrostatics(cylinder, CG[-1], verbose=True)
    # print '\nCenter of gravity new_location: ', CG
    # print '\nBuoyancy center               : ', hs_output['B']

    # show(vertices, faces, CG=CG, B=hs_output['B'])


