
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
    
