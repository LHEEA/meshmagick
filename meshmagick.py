#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python module to manipulate 2D meshes for hydrodynamics purposes

"""
This module contains utility function to manipulate, load, save and
convert mesh files
Two numpy arrays are manipulated in this module : V and F.
V is the array of nodes coordinates. It is an array of shape (nv, 3) where
nv is the number of nodes in the mesh.
F is the array of cell connectivities. It is an array of shape (nf, 4) where
nf is the number of cells in the mesh. Not that it has 4 columns as we consider
flat polygonal cells up to 4 edges (quads). Triangles are obtained by repeating
the first node at the end of the cell node ID list.

ATTENTION : IDs of vertices should always start at 1, which is the usual case in
file formats. However, for VTK, it starts at 0. To avoid confusion, we should
always start at 1 and do the appropriated conversions on F.
"""

import os, sys
import numpy as np

__author__ = "Francois Rongere"
__copyright__ = "Copyright 2014-2015, Ecole Centrale de Nantes"
__credits__ = ["Francois Rongere"]
__licence__ = "CeCILL"
__version__ = "0.1"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@ec-nantes.fr"
__status__ = "Development"

real_str = r'[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[Ee][+-]?\d+)?'


class HalfEdge:
    def __init__(self, Vorig, Vtarget):
        self.orig = Vorig
        self.target = Vtarget

    def reverse(self):
        self.orig, self.target = self.target, self.orig
        self.next, self.prev = self.prev, self.next


class Cell:
    def __init__(self, nodesIdx):
        if len(nodesIdx) != 4:
            raise ValueError, "The node list should have 4 elements"

        if nodesIdx.__class__ is not np.ndarray:
            nodesIdx = np.asarray(nodesIdx, dtype=np.int32, order='F')

        if nodesIdx[0] == nodesIdx[-1]:
            self.type = 3
        else:
            self.type = 4

        self.nodesIdx = nodesIdx[:self.type]

    def __iter__(self):
        self.__index = 0
        return self

    def next(self):
        if self.__index <= self.type - 1:
            Vorig = self.nodesIdx[self.__index]
            if self.__index == self.type - 1:
                Vtarget = self.nodesIdx[0]
            else:
                Vtarget = self.nodesIdx[self.__index + 1]
            curHE = HalfEdge(Vorig, Vtarget)
            curHE.__index = self.__index
            self.__index += 1

            return curHE
        else:  # Out of bounds
            raise StopIteration

    def area(self):
        return 0.0


class Mesh:
    def __init__(self, V, F):
        self.V = V
        self.F = F - 1
        self.nf = F.shape[0]
        self.nv = V.shape[0]

        V_oneOutgoingHE = [-1 for i in range(self.nv)]

        cells = []
        for cell in self.F:
            cell = Cell(cell)
            cells.append(cell)
            cell.idx = len(cells) - 1
        self.cells = cells

        self.nhe = 0
        self.HE = []
        edges = [[[], []] for j in range(self.nv - 1)]
        for cell in self.cells:
            for curHE in cell:
                self.nhe += 1
                curHE.idx = self.nhe - 1
                curHE.facet = cell.idx

                if V_oneOutgoingHE[curHE.orig] == -1:
                    V_oneOutgoingHE[curHE.orig] = curHE.idx

                if curHE._Cell__index == cell.type - 1:
                    curHE.prev = curHE.idx - 1
                    curHE.next = curHE.idx - cell.type + 1
                elif curHE._Cell__index == 0:
                    curHE.prev = curHE.idx + cell.type - 1
                    curHE.next = curHE.idx + 1
                else:
                    curHE.prev = curHE.idx - 1
                    curHE.next = curHE.idx + 1

                if curHE.orig < curHE.target:
                    idx = curHE.orig
                    V = curHE.target
                else:
                    idx = curHE.target
                    V = curHE.orig

                try:
                    ipos = edges[idx][0].index(V)
                    curHE.opposite = edges[idx][1][ipos]
                    self.HE[curHE.opposite].opposite = curHE.idx
                except:
                    edges[idx][0].append(V)
                    edges[idx][1].append(curHE.idx)

                self.HE.append(curHE)
            cell.oneHE = curHE.idx

        for he in self.HE:
            if hasattr(he, 'opposite'):
                he.border = False
            else:
                he.border = True
                print "(%u,%u) est une bordure !" % (he.orig, he.target)

        return


def merge_duplicates(V, F, verbose=False, tol=1e-8):
    """
    Returns a new node array where close nodes have been merged into one node (following tol). It also returns
    the connectivity array F with the new node IDs.
    :param V:
    :param F:
    :param verbose:
    :param tol:
    :return:
    """

    # TODO : Set a tolerance option in command line arguments
    nv, nbdim = V.shape

    levels = [0, nv]
    Vtmp = []
    iperm = np.array([i for i in xrange(nv)])

    for dim in range(nbdim):
        # Sorting the first dimension
        values = V[:, dim].copy()
        if dim > 0:
            values = values[iperm]
        levels_tmp = []
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            if istop-istart > 1:
                level_values = values[istart:istop]
                iperm_view = iperm[istart:istop]

                iperm_tmp = level_values.argsort()

                level_values[:] = level_values[iperm_tmp]
                iperm_view[:] = iperm_view[iperm_tmp]

                levels_tmp.append(istart)
                vref = values[istart]

                for idx in xrange(istart, istop):
                    cur_val = values[idx]
                    if np.abs(cur_val - vref) > tol:
                        levels_tmp.append(idx)
                        vref = cur_val

            else:
                levels_tmp.append(levels[ilevel])
        if len(levels_tmp) == nv:
            # No duplicate vertices
            if verbose:  # TODO : verify it with SEAREV mesh
                print "The mesh has no duplicate vertices"
            break

        levels_tmp.append(nv)
        levels = levels_tmp

    else:
        # Building the new merged node list
        Vtmp = []
        newID = np.array([i for i in xrange(nv)])
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            Vtmp.append(V[iperm[istart]])
            newID[iperm[range(istart, istop)]] = ilevel
        V = np.array(Vtmp, dtype=float, order='F')
        # Applying renumbering to cells
        for cell in F:
            cell[:] = newID[cell-1]+1

        if verbose:
            nv_new = V.shape[0]
            print "Initial number of nodes : {:d}".format(nv)
            print "New number of nodes     : {:d}".format(nv_new)
            print "{:d} nodes have been merged".format(nv-nv_new)

    return V, F


# =======================================================================
# MESH LOADERS
#=======================================================================
# Contains here all functions to load meshes from different file formats

def load_mesh(filename):
    """
    Function to load every known mesh file format
    """
    _, ext = os.path.splitext(filename)
    ext = ext.lower()
    if ext == '.vtk' or ext == '.vtu':
        V, F = load_VTK(filename)
    elif ext == '.gdf':
        V, F = load_GDF(filename)
    elif ext == '.mar':
        V, F = load_MAR(filename)
    elif ext == '.nat':
        V, F = load_NAT(filename)
    elif ext == '.stl':
        V, F = load_STL(filename)
    elif ext == '.msh':
        V, F = load_MSH(filename)
    elif ext == '.inp':
        V, F = load_INP(filename)
    elif ext == '.dat':
        raise NotImplementedError, "Not implemented"
        # V, F = load_DAT(filename)
    elif ext == '.tec':
        V, F = load_TEC(filename)
    elif ext == '.hst':
        V, F = load_HST(filename)
    elif ext == '.rad':
        V, F = load_RAD(filename)
    else:
        raise RuntimeError, 'extension %s is not recognized' % ext

    return V, F


def load_RAD(filename):
    """
    Loads RADIOSS files
    :param filename:
    :return:
    """

    import re

    ifile = open(filename, 'r')
    data = ifile.read()
    ifile.close()

    # node_line = r'\s*\d+(?:\s*' + real_str + '){3}'
    node_line = r'\s*\d+\s*(' + real_str + ')\s*(' + real_str + ')\s*(' + real_str + ')'
    node_section = r'((?:' + node_line + ')+)'

    elem_line = r'^\s*(?:\d+\s+){6}\d+\s*[\r\n]+'
    elem_section = r'((?:' + elem_line + '){3,})'

    pattern_node_line = re.compile(node_line, re.MULTILINE)
    pattern_node_line_group = re.compile(node_line, re.MULTILINE)
    pattern_elem_line = re.compile(elem_line, re.MULTILINE)
    pattern_node_section = re.compile(node_section, re.MULTILINE)
    pattern_elem_section = re.compile(elem_section, re.MULTILINE)

    V = []
    node_section = pattern_node_section.search(data).group(1)
    for node in pattern_node_line.finditer(node_section):
        V.append(map(float, list(node.groups())))
    V = np.asarray(V, dtype=float, order='F')

    F = []
    elem_section = pattern_elem_section.search(data).group(1)
    for elem in pattern_elem_line.findall(elem_section):
        F.append(map(int, elem.strip().split()[3:]))
    F = np.asarray(F, dtype=np.int32, order='F')

    return V, F


def load_HST(filename):
    """
    This function loads data from HYDROSTAR software.
    :param filename:
    :return:
    """
    ifile = open(filename, 'r')
    data = ifile.read()
    ifile.close()

    import re

    project = re.search(r'PROJECT\s*:\s*(.*)\s*', data).group(1).strip()

    node_line = r'\s*\d+(?:\s+' + real_str + '){3}'
    node_section = r'((?:' + node_line + ')+)'

    elem_line = r'^\s*(?:\d+\s+){3}\d+\s*[\r\n]+'
    elem_section = r'((?:' + elem_line + ')+)'

    pattern_node_line = re.compile(node_line, re.MULTILINE)
    pattern_elem_line = re.compile(elem_line, re.MULTILINE)
    pattern_node_section = re.compile(node_section, re.MULTILINE)
    pattern_elem_section = re.compile(elem_section, re.MULTILINE)

    Vtmp = []
    nv = 0
    for node_section in pattern_node_section.findall(data):
        for node in pattern_node_line.findall(node_section):
            Vtmp.append(map(float, node.split()[1:]))
        nvtmp = len(Vtmp)
        Vtmp = np.asarray(Vtmp, dtype=float, order='F')
        if nv == 0:
            V = Vtmp.copy()
            nv = nvtmp
        else:
            V = np.concatenate((V, Vtmp))
            nv += nvtmp

    Ftmp = []
    nf = 0
    for elem_section in pattern_elem_section.findall(data):
        for elem in pattern_elem_line.findall(elem_section):
            Ftmp.append(map(int, elem.split()))
        nftmp = len(Ftmp)
        Ftmp = np.asarray(Ftmp, dtype=np.int32, order='F')
        if nf == 0:
            F = Ftmp.copy()
            nf = nftmp
        else:
            F = np.concatenate((F, Ftemp))
            nf += nftmp

    return V, F


def load_INP(filename):
    """
    This function loads data from DIODORE (PRINCIPIA) INP file format.

    """
    import re

    ifile = open(filename, 'r')
    text = ifile.read()
    ifile.close()

    # End lines for every system
    # endl = r'(?:\n|\r|\r\n)'



    # Retrieving frames into a dictionnary frames
    pattern_FRAME_str = r'^\s*\*FRAME,NAME=(.+)[\r\n]+(.*)'
    pattern_FRAME = re.compile(pattern_FRAME_str, re.MULTILINE)

    frames = {}
    for match in pattern_FRAME.finditer(text):
        framename = match.group(1).strip()
        framevector = re.split(r'[, ]', match.group(2).strip())
        frames[framename] = np.asarray(map(float, framevector), order='F')

    # Storing the inp layout into a list of dictionnary
    pattern_NODE_ELEMENTS = re.compile(r'^\s*\*(NODE|ELEMENT),(.*)', re.MULTILINE)
    layout = []
    meshfiles = {}
    for match in pattern_NODE_ELEMENTS.finditer(text):
        fielddict = {}
        fielddict['type'] = match.group(1)
        if fielddict['type'] == 'NODE':
            fielddict['INCREMENT'] = 'NO'
        opts = match.group(2).split(',')
        for opt in opts:
            key, pair = opt.split('=')
            fielddict[key] = pair.strip()

        # Retrieving information on meshfiles and their usage
        file = fielddict['INPUT']
        if file in meshfiles:
            meshfiles[file][fielddict['type'] + '_CALL_INP'] += 1
        else:
            meshfiles[file] = {}
            meshfiles[file]['NODE_CALL_INP'] = 0
            meshfiles[file]['ELEMENT_CALL_INP'] = 0
            meshfiles[file][fielddict['type'] + '_CALL_INP'] += 1

        layout.append(fielddict)

        # RETRIEVING DATA SECTIONS FROM MESHFILES
        # patterns for recognition of sections
    node_line = r'\s*\d+(?:\s+' + real_str + '){3}'
    node_section = r'((?:' + node_line + ')+)'
    elem_line = r'^ +\d+(?: +\d+){3,4}[\r\n]+'  # 3 -> triangle, 4 -> quadrangle
    elem_section = r'((?:' + elem_line + ')+)'
    pattern_node_line = re.compile(node_line, re.MULTILINE)
    pattern_elem_line = re.compile(elem_line, re.MULTILINE)
    pattern_node_section = re.compile(node_section, re.MULTILINE)
    pattern_elem_section = re.compile(elem_section, re.MULTILINE)

    for file in meshfiles:
        try:
            meshfile = open(file + '.DAT', 'r')
        except:
            raise IOError, u'File {0:s} not found'.format(file + '.DAT')
        data = meshfile.read()
        meshfile.close()

        node_section = pattern_node_section.findall(data)
        if len(node_section) > 1:
            raise IOError, """Several NODE sections into a .DAT file is not supported by meshmagick
                              as it is considered as bad practice"""
        node_array = []
        idx_array = []
        for node in pattern_node_line.findall(node_section[0]):
            node = node.split()

            node[0] = int(node[0])
            idx_array.append(node[0])
            node[1:] = map(float, node[1:])
            node_array.append(node[1:])

        meshfiles[file]['NODE_SECTION'] = node_array

        # Detecting renumberings to do
        real_idx = 0
        # renumberings = []
        id_new = - np.ones(max(idx_array) + 1, dtype=int)
        for idx in idx_array:
            real_idx += 1
            if real_idx != idx:  # Node number and line number in the array are not consistant...
                id_new[idx] = real_idx
                # renumberings.append( (r'(?<!\n)(\s*)' + str(node[0]) + '(\s+)', str(idx))  )

        meshfiles[file]['ELEM_SECTIONS'] = []
        for elem_section in pattern_elem_section.findall(data):

            # Renumbering if any
            # for renumber in renumberings:
            #     elem_section = re.sub(renumber[0], '\g<1>'+renumber[1]+'\g<2>', elem_section)

            elem_array = []
            for elem in pattern_elem_line.findall(elem_section):
                elem = map(int, elem.split())
                # for node in elem[1:]:
                elem = id_new[elem[1:]].tolist()
                if len(elem) == 3:  # Case of a triangle, we repeat the first node at the last position
                    elem.append(elem[0])

                elem_array.append(map(int, elem))
            meshfiles[file]['ELEM_SECTIONS'].append(elem_array)
        meshfiles[file]['nb_elem_sections'] = len(meshfiles[file]['ELEM_SECTIONS'])

        meshfiles[file]['nb_elem_sections_used'] = 0

    nbNodes = 0
    nbElems = 0
    for field in layout:
        file = field['INPUT']
        if field['type'] == 'NODE':
            nodes = np.asarray(meshfiles[file]['NODE_SECTION'])
            # Translation of nodes according to frame option id any
            nodes = translate(nodes, frames[field['FRAME']])  # TODO: s'assurer que frame est une options obligatoire...

            if nbNodes == 0:
                V = nodes.copy(order='F')
                nbNodes = V.shape[0]
                increment = False
                continue

            if field['INCREMENT'] == 'NO':
                V[idx, :] = nodes.copy(order='F')
                increment = False
            else:
                V = np.concatenate((V, nodes))
                nbNodes = V.shape[0]
                increment = True
        else:  # this is an ELEMENT section
            elem_section = np.asarray(meshfiles[file]['ELEM_SECTIONS'][meshfiles[file]['nb_elem_sections_used']])

            meshfiles[file]['nb_elem_sections_used'] += 1
            if meshfiles[file]['nb_elem_sections_used'] == meshfiles[file]['nb_elem_sections']:
                meshfiles[file]['nb_elem_sections_used'] = 0

            # Updating to new id of nodes
            elems = elem_section
            if increment:
                elems += nbNodes

            if nbElems == 0:
                F = elems.copy(order='F')
                nbElems = F.shape[0]
                continue
            else:
                F = np.concatenate((F, elems))
                nbElems = F.shape[0]

    return V, F


def load_TEC(filename):
    """
    This function loads data from XML and legacy VTK file format.
    At that time, only unstructured meshes are supported.

    Usage:
        V, F = load_VTK(filename)
    """

    from vtk import vtkTecplotReader

    reader = vtkTecplotReader()

    # Importing the mesh from the file
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()

    nv = 0
    nf = 0

    for iblock in range(data.GetNumberOfBlocks()):
        block = data.GetBlock(iblock)
        if block.GetClassName() == 'vtkStructuredGrid':
            continue
        nvblock = block.GetNumberOfPoints()
        nfblock = block.GetNumberOfCells()

        Vtmp = np.zeros((nvblock, 3), dtype=float, order='F')
        for k in range(nvblock):
            Vtmp[k] = np.array(block.GetPoint(k))

        if nv == 0:
            V = Vtmp
        else:
            V = np.concatenate((V, Vtmp))

        nv += nvblock

        # Facet extraction
        Ftmp = np.zeros((nfblock, 4), dtype=np.int32, order='F')
        for k in range(nfblock):
            cell = block.GetCell(k)
            nv_facet = cell.GetNumberOfPoints()
            for l in range(nv_facet):
                Ftmp[k][l] = cell.GetPointId(l)
            if nv_facet == 3:
                Ftmp[k][l] = Ftmp[k][0]

        if nf == 0:
            F = Ftmp
        else:
            F = np.concatenate((F, Ftmp))

        nf += nfblock

    F += 1
    return V, F


def load_VTK(filename):
    """
    This function loads data from XML and legacy VTK file format.
    At that time, only unstructured meshes are supported.

    Usage:
        V, F = load_VTK(filename)
    """
    ext = os.path.splitext(filename)[1]

    if ext == '.vtu':  # XML file format for unstructured grid
        from vtk import vtkXMLUnstructuredGridReader

        reader = vtkXMLUnstructuredGridReader()
    elif ext == '.vtk':  # Legacy file format for unstructured grid
        from vtk import vtkUnstructuredGridReader

        reader = vtkUnstructuredGridReader()

    else:
        raise RuntimeError('Unknown file type %s ' % filename)

    # Importing the mesh from the file
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    nv = vtk_mesh.GetNumberOfPoints()
    V = np.zeros((nv, 3), dtype=float, order='fortran')
    for k in range(nv):
        V[k] = np.array(vtk_mesh.GetPoint(k))

    nf = vtk_mesh.GetNumberOfCells()
    F = np.zeros((nf, 4), dtype=np.int32, order='fortran')
    for k in range(nf):
        cell = vtk_mesh.GetCell(k)
        nv_facet = cell.GetNumberOfPoints()
        for l in range(nv_facet):
            F[k][l] = cell.GetPointId(l)
        if nv_facet == 3:
            F[k][3] = F[k][0]

    F += 1
    return V, F


def load_STL(filename):
    """
    This function reads an STL file to extract the mesh
    :param filename:
    :return:
    """

    from vtk import vtkSTLReader

    reader = vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()

    data = reader.GetOutputDataObject(0)

    nv = data.GetNumberOfPoints()
    V = np.zeros((nv, 3), dtype=float, order='F')
    for k in range(nv):
        V[k] = np.array(data.GetPoint(k))
    nf = data.GetNumberOfCells()
    F = np.zeros((nf, 4), dtype=np.int32, order='F')
    for k in range(nf):
        cell = data.GetCell(k)
        if cell is not None:
            for l in range(3):
                F[k][l] = cell.GetPointId(l)
                F[k][3] = F[k][0]  # always repeating the first node as stl is triangle only
    F += 1

    merge_duplicates(V, F)

    return V, F


def load_NAT(filename):
    """
    This function loads natural file format for meshes.

    Format spec :
    -------------------
    xsym    ysym
    n    m
    x1    y1    z1
    .
    .
    .
    xn    yn    zn
    i1    j1    k1    l1
    .
    .
    .
    im    jm    km    lm
    -------------------

    where :
    n : number of nodes
    m : number of cells
    x1 y1 z1 : cartesian coordinates of node 1
    i1 j1 k1 l1 : counterclock wise Ids of nodes for cell 1
    if cell 1 is a triangle, i1==l1
    """

    ifile = open(filename, 'r')
    xsym, ysym = map(int, ifile.readline().split())
    nv, nf = map(int, ifile.readline().split())

    V = []
    for i in range(nv):
        V.append(map(float, ifile.readline().split()))
    V = np.array(V, dtype=float, order='fortran')

    F = []
    for i in range(nf):
        F.append(map(int, ifile.readline().split()))
    F = np.array(F, dtype=np.int32, order='fortran')

    ifile.close()
    return V, F


def load_GDF(filename):
    """
    This function loads GDF files from WAMIT mesh file format.
    """
    ifile = open(filename, 'r')

    ifile.readline()  # skip one header line
    line = ifile.readline().split()
    ulen = line[0]
    grav = line[1]

    line = ifile.readline().split()
    isx = line[0]
    isy = line[1]

    nf = int(ifile.readline())

    V = np.zeros((4 * nf, 3), dtype=float, order='fortran')
    F = np.zeros((nf, 4), dtype=np.int32, order='fortran')

    iv = -1
    for icell in range(nf):

        for k in range(4):
            iv += 1
            V[iv, :] = np.array(ifile.readline().split())
            F[icell, k] = iv + 1

    ifile.close()
    V, F = merge_duplicates(V, F, verbose=True)

    return V, F


def load_MAR(filename):
    """
    This function loads .mar files in memory.
    """

    ifile = open(filename, 'r')

    ifile.readline()  # Skipping the first line of the file
    V = []
    while 1:
        line = ifile.readline()
        line = line.split()
        if line[0] == '0':
            break
        V.append(map(float, line[1:]))

    V = np.array(V, dtype=float, order='fortran')
    F = []
    while 1:
        line = ifile.readline()
        line = line.split()
        if line[0] == '0':
            break
        F.append(map(int, line))

    F = np.array(F, dtype=np.int32, order='fortran')

    ifile.close()

    return V, F


def load_STL2(filename):
    import re

    ifile = open(filename, 'r')
    text = ifile.read()
    ifile.close()

    endl = r'(?:\n|\r|\r\n)'
    patt_str = r"""
            ^\s*facet\s+normal(.*)""" + endl + """
            ^\s*outer\sloop""" + endl + """
            ^\s*vertex\s+(.*)""" + endl + """
            ^\s*vertex\s+(.*)""" + endl + """
            ^\s*vertex\s+(.*)""" + endl + """
            ^\s*endloop""" + endl + """
            ^\s*endfacet""" + endl + """
           """
    pattern = re.compile(patt_str, re.MULTILINE | re.VERBOSE)

    normal = []
    V = []
    for match in pattern.finditer(text):
        normal.append(map(float, match.group(1).split()))
        V.append(map(float, match.group(2).split()))
        V.append(map(float, match.group(3).split()))
        V.append(map(float, match.group(4).split()))

    V = np.array(V, dtype=float, order='fortran')

    nf = np.size(V, 0) / 3
    F = np.zeros((nf, 4), dtype=np.int32, order='fortran')

    base = np.array([1, 2, 3, 1])
    for i in range(nf):
        F[i, :] = base + 3 * i

    return V, F


def load_MSH(filename):
    import gmsh

    myMesh = gmsh.Mesh()
    myMesh.read_msh(filename)
    V = np.array(myMesh.Verts, dtype=float, order='fortran')

    ntri = myMesh.nElmts.get(2)
    nquad = myMesh.nElmts.get(3)
    if ntri is None:
        ntri = 0
    if nquad is None:
        nquad = 0

    nel = ntri + nquad

    F = np.zeros((nel, 4), dtype=np.int32, order='fortran')

    if ntri != 0:
        F[:ntri, :3] = myMesh.Elmts.get(2)[1]
        F[:, 3] = F[:, 0]

    if nquad != 0:
        F[ntri:, :] = myMesh.Elmts.get(3)[1]

    F += 1
    return V, F


#=======================================================================
#                             MESH WRITERS
#=======================================================================

def write_mesh(filename, V, F):
    """
    This function writes mesh data into filename following its extension
    """
    _, ext = os.path.splitext(filename)

    ext = ext.lower()
    if ext == '.vtk' or ext == '.vtu':
        write_VTK(filename, V, F)
    elif ext == '.gdf':
        write_GDF(filename, V, F)
    elif ext == '.mar':
        write_MAR(filename, V, F)
    elif ext == '.nat':
        write_NAT(filename, V, F)
    elif ext == '.stl':
        write_STL(filename, V, F)
    elif ext == '.msh':
        raise NotImplementedError, 'MSH file writer not implemented'
    elif ext == '.inp':
        raise NotImplementedError, 'INP file writer not implemented'
    elif ext == '.dat':
        write_DAT(filename, V, F)
    elif ext == '.tec':
        write_TEC(filename, V, F)
    elif ext == '.hst':
        write_HST(filename, V, F)
    else:
        raise RuntimeError, 'extension %s is not recognized' % ext


def write_DAT(filename, V, F):
    """
    Writes DAT files for DIODORE
    :param filename:
    :param V:
    :param F:
    :return:
    """

    import time
    import os

    rootfilename, ext = os.path.splitext(filename)
    filename = rootfilename+ext.upper()
    ofile = open(filename, 'w')

    ofile.write('$\n$ Data for DIODORE input file : {0}\n'.format(rootfilename.upper()))
    ofile.write('$ GENERATED BY MESHMAGICK ON {0}\n$\n'.format(time.strftime('%c')))

    ofile.write('$ NODE\n')
    vertex_block = \
        ''.join(
            (
                '\n'.join(
                    ''.join(
                        (
                            '{:8d}'.format(idx+1),
                            ''.join('{:13.5E}'.format(elt) for elt in node)
                        )
                    ) for (idx, node) in enumerate(V)
                ),

                '\n*RETURN\n'
            )
        )
    ofile.write(vertex_block)

    quad_block = '$\n$ ELEMENT,TYPE=Q4C000,ELSTRUCTURE={0}'.format(rootfilename.upper())
    tri_block  = '$\n$ ELEMENT,TYPE=T3C000,ELSTRUCTURE={0}'.format(rootfilename.upper())
    for (idx, cell) in enumerate(F):
        if cell[3] != cell[2]:
            # quadrangle
            quad_block = ''.join(
                (quad_block,
                 '\n',
                 '{:8d}'.format(idx+1),
                 ''.join('{:8d}'.format(node_id) for node_id in cell)
                )
            )


        else:
            # Triangle
            tri_block = ''.join(
                (tri_block,
                '\n',
                '{:8d}'.format(idx+1),
                ''.join('{:8d}'.format(node_id) for node_id in cell[:3])
                )
            )

    print '-------------------------------------------------'
    print 'Suggestion for .inp DIODORE input file :'
    print ''
    print '*NODE,INPUT={0},FRAME=???'.format(rootfilename)
    if quad_block != '':
        quad_block = ''.join((quad_block, '\n*RETURN\n'))
        ofile.write(quad_block)
        print '*ELEMENT,TYPE=Q4C000,ELSTRUCTURE={0},INPUT={0}'.format(rootfilename)
    if tri_block != '':
        tri_block = ''.join((tri_block, '\n*RETURN\n'))
        ofile.write(tri_block)
        print '*ELEMENT,TYPE=T3C000,ELSTRUCTURE={0},INPUT={0}'.format(rootfilename)

    print ''
    print '-------------------------------------------------'
    ofile.close()

    print 'File %s written' % filename

    return


def write_HST(filename, V, F):
    """
    This function writes mesh into a HST file format
    :param filename:
    :param V:
    :param F:
    :return:
    """

    ofile = open(filename, 'w')

    ofile.write(''.join((
        'PROJECT:\n',
        'USERS:   meshmagick\n\n'
        'NBODY   1\n'
        'RHO   1025.0\n'
        'GRAVITY   9.81\n\n'
    )))

    coordinates_block = ''.join((  # block
            'COORDINATES\n',
            '\n'.join(  # line
                ''.join(
                    (
                        '{:10d}'.format(idx+1),  # index
                        ''.join('{:16.6E}'.format(elt) for elt in node)  # node coordinates
                    )
                ) for (idx, node) in enumerate(V)
            ),
            '\nENDCOORDINATES\n\n'
    ))

    ofile.write(coordinates_block)

    cells_coordinates = ''.join((  # block
        'PANEL TYPE 0\n',
        '\n'.join(  # line
            ''.join(
                '{:10d}'.format(node_idx) for node_idx in cell
            ) for cell in F
        ),
        '\nENDPANEL\n\n'
    ))

    ofile.write(cells_coordinates)

    ofile.write('ENDFILE\n')

    ofile.close()

    print u'File {0:s} written'.format(filename)


def write_TEC(filename, V, F):
    """
    This function writes data in a tecplot file

    :param filename:
    :param V:
    :param F:
    :return:
    """
    ofile = open(filename, 'w')

    nv = V.shape[0]
    nf = F.shape[0]

    ofile.write('TITLE = \" THIS FILE WAS GENERATED BY MESHMAGICK - FICHIER : {} \" \n'.format(filename))

    ofile.write('VARIABLES = \"X\",\"Y\",\"Z\" \n')
    ofile.write('ZONE T=\"MESH\" \n')
    ofile.write('N={nv:10d} ,E={nf:10d} ,F=FEPOINT, ET=QUADRILATERAL\n'.format(nv=nv, nf=nf))

    node_block = '\n'.join( # block
        ''.join(
            ''.join('{:16.6E}'.format(elt) for elt in node)
        ) for node in V
    ) + '\n'
    ofile.write(node_block)

    cells_block = '\n'.join(  # block
        ''.join(
            ''.join('{:10d}'.format(node_id) for node_id in cell)
        ) for cell in F
    ) + '\n'
    ofile.write(cells_block)

    ofile.close()
    print 'File %s written' % filename


def write_VTK(filename, V, F):
    """ This function writes data in a VTK XML file.
    Currently, it only support writing unstructured grids
    """

    vtk_mesh = build_vtk_mesh_obj(V, F)
    # Writing it to the file
    _, ext = os.path.splitext(filename)
    if ext == '.vtk':
        from vtk import vtkUnstructuredGridWriter

        writer = vtkUnstructuredGridWriter()
    elif ext == '.vtu':
        from vtk import vtkXMLUnstructuredGridWriter

        writer = vtkXMLUnstructuredGridWriter()
        writer.SetDataModeToAscii()
    else:
        raise RuntimeError, 'Unknown file extension %s' % ext

    writer.SetFileName(filename)
    writer.SetInput(vtk_mesh)
    writer.Write()
    print 'File %s written' % filename


def build_vtk_mesh_obj(V, F):
    import vtk

    nv = max(np.shape(V))
    nf = max(np.shape(F))

    vtk_mesh = vtk.vtkUnstructuredGrid()
    vtk_mesh.Allocate(nf, nf)

    # Building the vtkPoints data structure
    vtk_points = vtk.vtkPoints()
    vtk_points.SetNumberOfPoints(nv)
    idx = -1
    for vertex in V:
        idx += 1
        vtk_points.SetPoint(idx, vertex)

    vtk_mesh.SetPoints(vtk_points)  # Storing the points into vtk_mesh

    # Building the vtkCell data structure
    F = F - 1
    for cell in F:
        if cell[-1] in cell[:-1]:
            vtk_cell = vtk.vtkTriangle()
            nc = 3
        else:
            # #print 'quadrangle'
            vtk_cell = vtk.vtkQuad()
            nc = 4

        for k in range(nc):
            vtk_cell.GetPointIds().SetId(k, cell[k])

        vtk_mesh.InsertNextCell(vtk_cell.GetCellType(), vtk_cell.GetPointIds())
    return vtk_mesh


def write_NAT(filename, V, F, *args):
    """
    This function writes mesh to file
    """
    ofile = open(filename, 'w')

    nv = max(np.shape(V))
    nf = max(np.shape(F))

    ofile.write('%6u%6u\n' % (0, 0))  # lire les symmetries dans args...
    ofile.write('%6u%6u\n' % (nv, nf))
    for vertex in V:
        ofile.write('%15.6E%15.6E%15.6E\n' % (vertex[0], vertex[1], vertex[2]))
    for cell in F:
        ofile.write('%10u%10u%10u%10u\n' % (cell[0], cell[1], cell[2], cell[3]))

    ofile.close()
    print 'File %s written' % filename


def write_GDF(filename, V, F, *args):
    """
    This function writes mesh data into a GDF file for Wamit computations
    """

    nf = max(np.shape(F))

    ofile = open(filename, 'w')

    ofile.write('GDF file generated by meshmagick\n')

    ofile.write('%16.6f%16.6f\n' % (100.0, 9.81))
    ofile.write('%12u%12u\n' % (0, 1))  # TODO : mettre les symetries en argument
    ofile.write('%12u\n' % nf)

    for cell in F:
        for k in range(4):
            Vcur = V[cell[k] - 1, :]
            ofile.write('%16.6E%16.6E%16.6E\n' % (Vcur[0], Vcur[1], Vcur[2]))

    ofile.close()
    print 'File %s written' % filename


def write_MAR(filename, V, F, *args):
    ofile = open(filename, 'w')

    ofile.write('{0:6d}{1:6d}\n'.format(2, 0))  # TODO : mettre les symetries en argument

    nv = V.shape[0]
    for (idx, vertex) in enumerate(V):
        ofile.write('{0:6d}{1:16.6f}{2:16.6f}{3:16.6f}\n'.format(idx+1, vertex[0], vertex[1], vertex[2]))

    ofile.write('{0:6d}{1:6d}{2:6d}{3:6d}{4:6d}\n'.format(0, 0, 0, 0, 0))

    cell_block = '\n'.join(
        ''.join(u'{0:10d}'.format(elt) for elt in cell)
        for cell in F
    ) + '\n'
    ofile.write(cell_block)
    ofile.write('%6u%6u%6u%6u\n' % (0, 0, 0, 0))

    ofile.close()
    print 'File %s written' % filename


def write_STL(filename, V, F):
    """
    :type filename: str
    """

    # TODO : replace this implementation by using the vtk functionalities

    ofile = open(filename, 'w')

    ofile.write('solid meshmagick\n')
    F -= 1  # STL format specifications tells that numerotation starts at 0

    for facet in F:
        if facet[0] != facet[3]:
            raise RuntimeError, 'Only full triangle meshes are accepted in STL files'

        # Computing normal
        v0 = V[facet[0], :]
        v1 = V[facet[1], :]
        v2 = V[facet[2], :]

        n = np.cross(v1 - v0, v2 - v0)
        n /= np.linalg.norm(n)

        block_facet = ''.join(['  facet normal ', ''.join('%15.6e' % ni for ni in n) + '\n',
                               '    outer loop\n',
                               '      vertex', ''.join('%15.6e' % Vi for Vi in v0) + '\n',
                               '      vertex', ''.join('%15.6e' % Vi for Vi in v1) + '\n',
                               '      vertex', ''.join('%15.6e' % Vi for Vi in v2) + '\n',
                               '    endloop\n',
                               '  endfacet\n'])
        ofile.write(block_facet)
    ofile.write('endsolid meshmagick\n')
    ofile.close()

    print 'File %s written' % filename


#=======================================================================
#                         MESH MANIPULATION HELPERS
#=======================================================================
def mesh_quality(V, F):
    # This function is reproduced from
    # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/Verdict/Testing/Python/MeshQuality.py
    import vtk
    import math

    vtk_mesh = build_vtk_mesh_obj(V, F)
    quality = vtk.vtkMeshQuality()
    quality.SetInput(vtk_mesh)

    def DumpQualityStats(iq, arrayname):
        an = iq.GetOutput().GetFieldData().GetArray(arrayname)
        cardinality = an.GetComponent(0, 4)
        range = list()
        range.append(an.GetComponent(0, 0))
        range.append(an.GetComponent(0, 2))
        average = an.GetComponent(0, 1)
        stdDev = math.sqrt(math.fabs(an.GetComponent(0, 3)))
        outStr = '%s%g%s%g\n%s%g%s%g' % (
                '    range: ', range[0], '  -  ', range[1],
                '    average: ', average, '  , standard deviation: ', stdDev)
        return outStr

    # Here we define the various mesh types and labels for output.
    meshTypes = [
                ['Triangle', 'Triangle',
                 [['QualityMeasureToArea', ' Area Ratio:'],
                  ['QualityMeasureToEdgeRatio', ' Edge Ratio:'],
                  ['QualityMeasureToAspectRatio', ' Aspect Ratio:'],
                  ['QualityMeasureToRadiusRatio', ' Radius Ratio:'],
                  ['QualityMeasureToAspectFrobenius', ' Frobenius Norm:'],
                  ['QualityMeasureToMinAngle', ' Minimal Angle:']
                 ]
                 ],

                ['Quad', 'Quadrilateral',
                 [['QualityMeasureToArea', ' Area Ratio:'],
                  ['QualityMeasureToEdgeRatio', ' Edge Ratio:'],
                  ['QualityMeasureToAspectRatio', ' Aspect Ratio:'],
                  ['QualityMeasureToRadiusRatio', ' Radius Ratio:'],
                  ['QualityMeasureToMedAspectFrobenius',
                  ' Average Frobenius Norm:'],
                  ['QualityMeasureToMaxAspectFrobenius',
                  ' Maximal Frobenius Norm:'],
                  ['QualityMeasureToMinAngle', ' Minimal Angle:']
                 ]
                ]
                ]

    if vtk_mesh.GetNumberOfCells() > 0:
        res = ''
        for meshType in meshTypes:
            res += '\n%s%s' % (meshType[1], ' quality of the mesh ')
            quality.Update()
            an = quality.GetOutput().GetFieldData().GetArray('Mesh ' + meshType[1] + ' Quality')
            cardinality = an.GetComponent(0, 4)

            res = ''.join((res, '(%u elements):\n' % (cardinality)))

            # res += '('+str(cardinality) +meshType[1]+'):\n'

            for measure in meshType[2]:
                eval('quality.Set' + meshType[0] + measure[0] + '()')
                quality.Update()
                res += '\n%s\n%s' % (measure[1],
                        DumpQualityStats(quality,
                                 'Mesh ' + meshType[1] + ' Quality'))
            res += '\n'

    info = """\n\nDefinition of the different quality measures is given
in the verdict library manual :
http://www.vtk.org/Wiki/images/6/6b/VerdictManual-revA.pdf\n"""
    res += info
    return vtk_mesh, res


def get_info(V, F):
    nv = np.size(V, 0)
    nf = np.size(F, 0)
    print ''
    print 'o--------------------------------------------------o'
    print '|               MESH CHARACTERISTICS               |'  #28
    print '|--------------------------------------------------|'
    print '| Number of nodes  :     %15u           |' % nv
    print '|--------------------------------------------------|'
    print '| Number of facets :     %15u           |' % nf
    print '|--------------------------------------------------|'  #51
    print '|      |          Min        |          Max        |'
    print '|------|---------------------|---------------------|'
    print '|   X  |%21E|%21E|' % (V[:, 0].min(), V[:, 0].max())
    print '|------|---------------------|---------------------|'
    print '|   Y  |%21E|%21E|' % (V[:, 1].min(), V[:, 1].max())
    print '|------|---------------------|---------------------|'
    print '|   Z  |%21E|%21E|' % (V[:, 2].min(), V[:, 2].max())
    print 'o--------------------------------------------------o'
    print ''

    _, res = mesh_quality(V, F)
    print res


def translate(V, P):
    """
    Translates an array V with respect to translation vector P.
    :param V:
    :param P:
    :return:
    """

    if not isinstance(P, np.ndarray):
        P = np.asarray(P, dtype=float)

    try:
        for i in range(np.size(V, 0)):
            V[i, :] += P
    except:
        raise RuntimeError, 'second argument must be a 3D list or numpy array for the translation'

    return V


def translate_1D(V, t, ddl):
    if ddl == 'x':
        j = 0
    elif ddl == 'y':
        j = 1
    elif ddl == 'z':
        j = 2
    else:
        raise IOError, "ddl should be chosen among ('x', 'y', 'z')"
    V[:, j] += t
    return V


def rotate(V, rot):
    from math import cos, sin

    if not isinstance(rot, np.ndarray):
        rot = np.asarray(rot, dtype=float)

    R = np.zeros((3, 3), dtype=float, order='f')

    phi = rot[0]
    theta = rot[1]
    psi = rot[2]

    cphi = cos(phi)
    sphi = sin(phi)
    ctheta = cos(theta)
    stheta = sin(theta)
    cpsi = cos(psi)
    spsi = sin(psi)

    R[0, 0] = cpsi * ctheta
    R[0, 1] = -spsi * cphi + cpsi * stheta * sphi
    R[0, 2] = spsi * sphi + cpsi * cphi * stheta
    R[1, 0] = spsi * ctheta
    R[1, 1] = cpsi * cphi + sphi * stheta * spsi
    R[1, 2] = -cpsi * sphi + spsi * cphi * stheta
    R[2, 0] = -stheta
    R[2, 1] = ctheta * sphi
    R[2, 2] = ctheta * cphi

    return np.transpose(np.dot(R, V.T))


def rotate_1D(V, rot, ddl):
    if ddl == 'x':
        j = 0
    elif ddl == 'y':
        j = 1
    elif ddl == 'z':
        j = 2
    else:
        raise IOError, "ddl should be chosen among ('x', 'y', 'z')"

    rotvec = np.zeros(3, dtype=float)
    rotvec[j] = rot
    return rotate(V, rotvec)


def scale(V, alpha):
    return alpha * V


def flip_normals(F):
    return np.fliplr(F)


def _fix_python_windows_install():
    """
    Fix a bug in the install of python into Windows. It modifies a registry key in order
    to be able to use python scripts in the DOS prompt with command line arguments.
    In the initial install, the system is unable to catch command line arguments.
    :return:
    """

    try:
        import _winreg as wr
    except:
        import sys
        raise ImportError, "This function only concerns Windows " \
                                "environment ! your environment is {:s}".format(sys.platform)

    try:
        key = wr.OpenKey(HKEY_CLASSES_ROOT, r'py_auto_file\shell\open\command', 0, KEY_ALL_ACCESS)
    except:
        raise OSError, "The key is not present on you system, " \
                       "please check you have properly set the .py file association"

    value = wr.QueryValue(key, '')

    if value.find('%*') != -1:
        print "The fix has already been done."
        wr.CloseKey(key)
        return

    value = ' '.join((value, '%*'))

    try:
        wr.SetValueEx(key, '', 0, REG_SZ, value)
        wr.CloseKey(key)
    except:
        raise OSError, """The function failed to change the registry key. You will have to do it by yourself.
                        Please append %* to the default value of the following registry key :
                        HKEY_CLASSES_ROOT\py_auto_file\shell\open\command\
                        """


    return

# =======================================================================
#                         COMMAND LINE USAGE
# =======================================================================
if __name__ == '__main__':
    import argparse
    import sys
	
    try:
        import argcomplete

        acok = True
    except:
        acok = False

    parser = argparse.ArgumentParser(
        description="""  --  MESHMAGICK --
                    A python module and a command line utility to manipulate meshes from different format used in
                    hydrodynamics as well as for visualization.

                    The formats currently supported by meshmagick are :

                    *---------------*-------------*-----------------*-----------------------*
                    | File          | R: Reading  | Software        | Keywords              |
                    | extension     | W: writing  |                 |                       |
                    *---------------*-------------*-----------------*-----------------------*
                    |     .mar      |    R/W      | NEMOH (1)       | nemoh, mar            |
                    |     .gdf      |    R/W      | WAMIT (2)       | wamit, gdf            |
                    |     .inp      |    R        | DIODORE (3)     | diodore-inp, inp      |
                    |     .DAT      |    R/W      | DIODORE (3)     | diodore-dat           |
                    |     .hst      |    R/W      | HYDROSTAR (4)   | hydrostar, hst        |
                    |     .nat      |    R/W      |    -            | natural, nat          |
                    |     .msh      |    R        | GMSH (5)        | gmsh, msh             |
                    |     .stl      |    R/W      |    -            | stl                   |
                    |     .vtu      |    R/W      | PARAVIEW (6)    | paraview, vtu         |
                    |     .vtk      |    R/W      | PARAVIEW (6)    | paraview-legacy, vtk  |
                    |     .tec      |    R/W      | TECPLOT (7)     | tecplot, tec          |
                    *---------------*-------------------------------------------------------*

                    By default, Meshmagick uses the filename extensions to choose the appropriate reader/writer.
                    This behaviour might be bypassed using the -ifmt and -ofmt optional arguments. When using these
                    options, keywords defined in the table above must be used as format identifiers.

                    (1) NEMOH is an open source BEM Software for seakeeping developped at Ecole Centrale de Nantes (LHHEA)
                    (2) WAMIT is a BEM Software for seakeeping developped by WAMIT, Inc.
                    (3) DIODORE is a BEM Software for seakeeping developped by PRINCIPIA
                    (4) HYDROSTAR is a BEM Software for seakeeping developped by BUREAU VERITAS
                    (5) GMSH is an open source meshing software developped by C. Geuzaine and J.-F. Remacle
                    (6) PARAVIEW is an open source visualization software developped by Kitware
                    (7) TECPLOT is a visualization software developped by Tecplot


                    """,
        epilog="""             --  Copyright 2014-2015  --
        --  Francois Rongere  /  Ecole Centrale de Nantes  --""",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('infilename',
                        help='path of the input mesh file in any format')
    parser.add_argument('-o', '--outfilename',
                        help='path of the output mesh file. The format of ' +
                             'this file is determined from the extension given')
    parser.add_argument('-ifmt', '--input-format',
                        help="""Input format. Meshmagick will read the input file considering the INPUT_FORMAT rather than using the extension""")
    parser.add_argument('-ofmt', '--output-format',
                        help="""Output format. Meshmagick will write the output file considering the OUTPUT_FORMAT rather than using the extension""")
    # parser.add_argument('-fmt', '--format', type=str,
    #                     help="""specifies the output file format and use the base name
    #                     of the infilename instead of the output file. File format has
    #                     to be chosen among VTK, VTU, MAR, GDF, STL, NAT, MHE. Default is VTU.
    #                     If -o option has been done, it has no effect""")
    parser.add_argument('-v', '--verbose',
                        help="""make the program give more informations on the computations""",
                        action='store_true')
    parser.add_argument('-i', '--info',
                        help="""extract informations on the mesh on the standard output""",
                        action='store_true')
    # parser.add_argument('-b', '--binary',
    #                     help="""specifies wether the output file will be in a binary
    #                     format. It is only used when used with --convert option and
    #                     with .vtu output format""",
    #                     action='store_true')

    parser.add_argument('-t', '--translate',
                        nargs=3, type=float,
                        help="""translates the mesh in 3D
                        Usage -translate tx ty tz""")
    parser.add_argument('-tx', '--translatex',
                        nargs=1, type=float,
                        help="""translates the mesh following the x direction""")
    parser.add_argument('-ty', '--translatey',
                        nargs=1, type=float,
                        help="""translates the mesh following the y direction""")
    parser.add_argument('-tz', '--translatez',
                        nargs=1, type=float,
                        help="""translates the mesh following the z direction""")

    parser.add_argument('-r', '--rotate',
                        nargs=3, type=str,
                        help="""rotates the mesh in 3D""")
    parser.add_argument('-rx', '--rotatex',
                        nargs=1, type=str,
                        help="""rotates the mesh around the x direction""")
    parser.add_argument('-ry', '--rotatey',
                        nargs=1, type=str,
                        help="""rotates the mesh around the y direction""")
    parser.add_argument('-rz', '--rotatez',
                        nargs=1, type=str,
                        help="""rotates the mesh around the z direction""")
    # parser.add_argument('-u', '--unit',
    #                     type=str, choices=['rad', 'deg'],
    #                     default='deg',
    #                     help="""sets the unit for rotations. Default is deg""")
    parser.add_argument('-s', '--scale',
                        type=float,
                        help="""scales the mesh. CAUTION : if used along
                         with a translation option, the scaling is done before
                        the translations. The translation magnitude should be set
                        accordingly to the newly scaled mesh.""")
    parser.add_argument('--flip-normals', action='store_true',
                        help="""flips the normals of the mesh""")
    # parser.add_argument('--cut', action='store_true',
    #                     help="""cuts the mesh with the plane defined with the option
    #                     --plane""")
    # parser.add_argument('--symmetrize', type=str,
    #                     nargs='*', default=None, choices=['x', 'y', 'z', 'p'],
    #                     help="""performs a symmetry of the mesh.
    #                     If x argument is given, then a symmetry with respect to plane
    #                     0yz etc... If p or no argument is given, then the --plane option
    #                     is used.""")
    # parser.add_argument('-p', '--plane', nargs=4, type=float,
    #                     help="""Defines a plane used by the --cut option.
    #                     It is followed by the floats nx ny nz c where [nx, ny, nz]
    #                     is a normal vector to the plane and c defines its position
    #                     following the equation <N|X> = c with X a point belonging
    #                     to the plane""")
    parser.add_argument('-m', '--merge-duplicates', nargs='?', const='1e-8', default='1e-8',
                        help="""merges the duplicate nodes in the mesh with the absolute tolerance
                        given as argument (default 1e-8)""")


    # parser.add_argument('--renumber', action='store_true',
    #                     help="""renumbers the cells and nodes of the mesh so as
    #                     to optimize cache efficiency in algorithms. Uses the Sloan
    #                     method""")
    # parser.add_argument('--optimize', action='store_true',
    #                     help="""optimizes the mesh. Same as --merge-duplicates
    #                     and --renumber used together""")
    # parser.add_argument('--fix-windows', action='store_true',
		# 				help="""Fix the python installation to be able to run
		# 				python scripts with command line arguments. It is not a
		# 				meshmagick fix but a windows fix. It modifies a registry key
		# 				of your Windows installation. It has to be run once, althought
		# 				the key is checked to ensure that the key has not already been
		# 				modified.""")

    if acok:
        argcomplete.autocomplete(parser)

    args, unknown = parser.parse_known_args()

    tol = float(args.merge_duplicates)

    # if args.fix_windows:
		# _fix_python_windows_install()
		# sys.exit(1)
    

    extension_dict = ('vtk', 'vtu', 'gdf', 'mar', 'nat', 'stl', 'msh', 'inp', 'dat', 'tec', 'hst', 'rad')

    write_file = False  # switch to decide if data should be written to outfilename

    _, ext = os.path.splitext(args.infilename)
    if ext[1:].lower() not in extension_dict:
        raise IOError, 'Extension "%s" is not known' % ext

    if args.outfilename is not None:
        _, ext = os.path.splitext(args.outfilename)
        write_file = True
        if ext[1:].lower() not in extension_dict:
            raise IOError, 'Extension "%s" is not known' % ext
    # else:
        # if args.format is not None:
        #     write_file = True
        #     if args.format.lower() not in extension_dict:
        #         raise IOError, 'Extension ".%s" is not known' % args.format
        #     root, _ = os.path.splitext(args.infilename)
        #     args.outfilename = root + '.' + args.format
        # else:
        #     args.outfilename = args.infilename




    # Importing data from file
    V, F = load_mesh(args.infilename)

    # myMesh = Mesh(V, F)

    # Dealing with different options
    # if args.optimize:
    #     args.merge_duplicates = True
    #     args.renumber = True

    if args.merge_duplicates:
        try:
            from pymesh import pymesh

            pymesh.set_mesh_data(V, F)
            pymesh.cleandata_py()
            V = np.copy(pymesh.vv)
            F = np.copy(pymesh.ff)
            pymesh.free_mesh_data()
            # TODO: implement the setting of tolerance externally into pymesh
        except:
            V, F = merge_duplicates(V, F, verbose=args.verbose, tol=tol)

    # if args.renumber:
    #     raise NotImplementedError, "Renumbering is not implemented yet into meshmagick"
        # try:
        #     from pymesh import pymesh
        # except:
        #     raise ImportError, 'pymesh module is not available, unable to use the --renumber option'
        #     ##raise NotImplementedError, '--renumber option is not implemented yet'

    # if args.binary:
    #     raise NotImplementedError, 'Not implemented yet'

    if args.translate is not None:
        V = translate(V, args.translate)
        write_file = True
    if args.translatex is not None:
        V = translate_1D(V, args.translatex, 'x')
        write_file = True
    if args.translatey is not None:
        V = translate_1D(V, args.translatey, 'y')
        write_file = True
    if args.translatez is not None:
        V = translate_1D(V, args.translatez, 'z')
        write_file = True

    # def cast_angles(angles, unit):
    #     from math import pi
    #
    #     if unit == 'deg':
    #         try:
    #             angles = map(float, angles)
    #         except:
    #             raise IOError, 'Bad input for rotation arguments'
    #         angles = np.array(angles) * pi / 180.
    #     else:
    #         def evalpi(elt):
    #             from math import pi
    #
    #             elttemp = elt.replace('pi', '1.')
    #             try:
    #                 elttemp = eval(elttemp)
    #             except:
    #                 raise IOError, 'Bad input %s' % elt
    #             return float(eval(elt))
    #
    #         angles = map(evalpi, angles)
    #     return angles

    if args.rotate is not None:
        args.rotate = cast_angles(args.rotate, args.unit)
        V = rotate(V, args.rotate)

    if args.rotatex is not None:
        args.rotatex = cast_angles(args.rotatex, args.unit)
        V = rotate_1D(V, args.rotatex[0], 'x')
    if args.rotatey is not None:
        args.rotatey = cast_angles(args.rotatey, args.unit)
        V = rotate_1D(V, args.rotatey[0], 'y')
    if args.rotatez is not None:
        args.rotatez = cast_angles(args.rotatez, args.unit)
        V = rotate_1D(V, args.rotatez[0], 'z')

    if args.scale is not None:
        V = scale(V, args.scale)

    # if args.symmetrize is not None:
    #     raise NotImplementedError, "Symmetrization is not implemented yet into meshmagick"
        # try:
        #     from pymesh import pymesh
        # except:
        #     raise ImportError, 'pymesh module is not available, unable to use the --symmetrize option'

        # if 'p' in args.symmetrize:
        #     if len(args.symmetrize) > 1:
        #         raise RuntimeError, 'When using p argument of --symmetrize option, only one argument should be provided'
        #     else:
        #         # Using plane to symmetrize
        #         if args.plane is None:
        #             raise RuntimeError, 'When using --symmetrize p option argument, you must also provide the --plane option'
        #         pass
        #         # Idee ici on se debrouille pour tourner le maillage afin que le plan de --plane soit
        #         # le plan yOz pour s'appuer sur la fct Fortran de pymesh, on symmetrize puis on retransforme
        #         # le maillage pour le remettre dans la bonne direction.
        #     raise NotImplementedError, '--symmetrize option is not implemented yet for argument p'
        # else:
        #     sym = np.zeros(3, dtype=np.int32)
        #     for plane in args.symmetrize:
        #         if plane == 'x':
        #             sym[0] = 1
        #         elif plane == 'y':
        #             sym[1] = 1
        #         else:
        #             sym[2] = 1
        #     pymesh.set_mesh_data(V, F)
        #     pymesh.symmetrizemesh_py(sym)
        #     V = np.copy(pymesh.vv)
        #     F = np.copy(pymesh.ff)
        #     pymesh.free_mesh_data()

    # if args.cut:
    #     raise NotImplementedError, '--cut option is not implemented yet'

    if args.flip_normals:
        F = flip_normals(F)

    if args.info:
        get_info(V, F)

    if write_file:
        write_mesh(args.outfilename, V, F)

