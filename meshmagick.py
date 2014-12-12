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
__copyright__ = "Copyright 2014, Ecole Centrale de Nantes"
__credits__ = ["Francois Rongere"]
__licence__ = "CeCILL"
__version__ = "0.1"
__maintainer__ = "Francois Rongere"
__email__ = "Francois.Rongere@ec-nantes.fr"
__status__ = "Development"

real_str = r'[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[Ee][+-]?\d+)?'


# class HalfEdge:
#     def __init__(self, Vorig, Vtarget):
#         self.orig = Vorig
#         self.target = Vtarget
#     def reverse(self):
#         self.orig, self.target = self.target, self.orig
#         self.next, self.prev = self.prev, self.next
#
#
# class Cell:
#     def __init__(self, nodesIdx):
#         if len(nodesIdx) != 4:
#             raise ValueError, "The node list should have 4 elements"
#
#         if nodesIdx.__class__ is not np.ndarray:
#             nodesIdx = np.asarray(nodesIdx, dtype=np.int32, order='F')
#
#         if nodesIdx[0] == nodesIdx[-1]:
#             self.type = 3
#         else:
#             self.type = 4
#
#         self.nodesIdx = nodesIdx[:self.type]
#
#     def __iter__(self):
#         self.__index = 0
#         return self
#
#     def next(self):
#         if self.__index <= self.type - 1:
#             Vorig = self.nodesIdx[self.__index]
#             if self.__index == self.type - 1:
#                 Vtarget = self.nodesIdx[0]
#             else:
#                 Vtarget = self.nodesIdx[self.__index + 1]
#             curHE = HalfEdge(Vorig, Vtarget)
#             curHE.__index = self.__index
#             self.__index += 1
#
#             return curHE
#         else:  # Out of bounds
#             raise StopIteration
#
#     def area(self):
#
#
#
#         return 0.0


class Mesh:
    def __init__(self, V, F):
        # Storing vertices
        if V.shape[0] > 3:
            self.Vcoord = V.T
        else:
            self.Vcoord = V
        self.nv = self.Vcoord.shape[1]

        # if 0 not in F:
        #     F -= 1

        if F.shape[0] == 4:
            F = F.T

        self.nf = F.shape[0]

        def store_edges_informations(iV1, iV2, ihe, edges):
            # Storing edge informations
            iVmin = iV1
            iVmax = iV2
            if iVmin > iVmax:
                iVmin, iVmax = iVmax, iVmin

            if not edges.has_key(iVmin):
                edges[iVmin] = {}
                edges[iVmin][iVmax] = [ihe]
            else:
                if not edges[iVmin].has_key(iVmax):
                    edges[iVmin][iVmax] = [ihe]
                else:
                    edges[iVmin][iVmax].append(ihe)
            return edges

        # Building half-edge data structure
        edges = {}
        self.V_oneOutgoingHE = [-1 for i in range(self.nv)]
        self.F_1HE = [-1 for i in range(self.nf)]
        self.HE_targetV = []
        self.HE_F = []
        self.HE_nextHE = []
        self.HE_prevHE = []
        HE_origV = [] # FIXME : Utile ?

        nhe = 0
        icell = -1
        for cell in F-1:

            icell += 1
            istop = 3
            if cell[0] == cell[-1]:  # Triangle
                istop -= 1

            # Iterating over the half-edges of the cell
            for idx in range(istop):
                # New half-edge
                ihe = nhe
                nhe += 1
                HE_origV.append(cell[idx])
                self.HE_targetV.append(cell[idx+1])

                if self.V_oneOutgoingHE[cell[idx]] == -1:
                    self.V_oneOutgoingHE[cell[idx]] = ihe

                if idx == 0:
                    self.HE_prevHE.append(ihe+istop)
                    self.HE_nextHE.append(ihe+1)
                else:
                    self.HE_prevHE.append(ihe-1)
                    self.HE_nextHE.append(ihe+1)

                self.HE_F.append(icell)

                # Storing edge informations
                edges = store_edges_informations(cell[idx], cell[idx+1], ihe, edges)

            # New half-edge
            ihe = nhe
            nhe += 1
            HE_origV.append(cell[idx])
            self.HE_targetV.append(cell[0])
            self.HE_prevHE.append(ihe-1)
            self.HE_nextHE.append(ihe-istop)
            self.HE_F.append(icell)

            # Storing edge informations
            edges = store_edges_informations(cell[istop], cell[0], ihe, edges)

            self.F_1HE[icell] = ihe

        self.nhe = nhe

        # Getting connectivity of twin half-edges
        self.HE_twinHE = [-1 for i in range(nhe)]
        self.HE_boundary = [False for i in range(nhe)]
        boundaryF = []
        for iVmin in edges.keys():
            for iVmax in edges[iVmin].keys():

                heList = edges[iVmin][iVmax]
                helen = len(heList)
                if helen == 0:
                    raise "Erreur!!!"
                elif helen == 1:
                    # This is a boundary half-edge
                    ihe = heList[0]
                    self.HE_boundary[ihe] = True
                    boundaryF.append(self.HE_F[ihe])
                elif helen == 2:
                    # We have a twin
                    self.HE_twinHE[heList[0]] = heList[1]
                    self.HE_twinHE[heList[1]] = heList[0]
                else:
                    # More than 2 half-edges are following the same edge...
                    print "More than 2 half-edges are following the same edge"
                    print "Facets :"
                    for ihe in heList:
                        print self.HE_F[ihe]
        print boundaryF
        return


def merge_duplicates(V, F, verbose=False):
    nv, nbdim = V.shape

    tol = 1e-8

    blocks = [np.array([j for j in range(nv)])]
    Vtmp = []
    iperm = np.zeros((nv), dtype=np.int32)
    for dim in range(nbdim):
        # Sorting the first dimension
        newBlocks = []
        for block in blocks:
            col = V[block, dim]
            # Sorting elements in ascending order along dimension dim
            indices = np.argsort(col)  # Indices are here relative to the block array
            diff = np.abs(np.diff(col[indices])) < tol
            if dim > 0:
                indices = block[indices]  # Making the indices relative to the global node list

            # Forming blocks of contiguous values
            newBlock = [indices[0]]
            for idx in range(len(block)-1):  # -1 as size of diff is size of vector -1
                if diff[idx]:  # elements at index idx and idx+1 are equal
                    newBlock.append(indices[idx+1])
                else:
                    # We add an additional block if it has more than one element
                    if len(newBlock) > 1:
                        newBlock = np.array(newBlock)
                        newBlock.sort()
                        newBlocks.append(newBlock)
                        if dim == nbdim-1:
                            Vtmp.append(V[newBlock[0], :])
                            iperm[newBlock] = len(Vtmp)-1
                    else:
                        Vtmp.append(V[newBlock[0], :])  # Keeping the vertex as it is a singleton
                        iperm[newBlock] = len(Vtmp)-1
                    newBlock = [indices[idx+1]]
            if len(newBlock) > 1:
                newBlock = np.array(newBlock)
                newBlock.sort()
                newBlocks.append(np.array(newBlock))
                if dim == nbdim-1:
                    Vtmp.append(V[newBlock[0], :])
                    iperm[newBlock] = len(Vtmp)-1
            else:
                Vtmp.append(V[newBlock[0], :])  # Keeping the vertex as it is a singleton
                iperm[newBlock] = len(Vtmp)-1

        # Termination condition if every block in newBlocks is a singleton
        if len(newBlocks) == 0:
            if verbose:
                print "No duplicate nodes"
            return V, F

        blocks = newBlocks

    # New array of vertices
    V = np.array(Vtmp, dtype=float, order='F')
    nvNew = V.shape[0]


    # Updating the connectivities
    nf = F.shape[0]
    for idx in range(nf):
        F[idx, :] = iperm[F[idx, :]-1]+1

    if verbose :
        print 'Initial number of vertices : %u' % nv
        print 'New number of vertices     : %u' % nvNew
        print "%u nodes have been merged" % (nv-nvNew)

    return V, F


def merge_cells(F):

    # Sorting the connectivities in order to make the comparison efficient
    nf = F.shape[0]
    F_backup = F.copy()
    F.sort()
    indices = F[:,0].argsort(axis=0)
    # print F[indices, :]
    F = F[indices]
    precell = F[-1,:]
    # print precell
    # print '--'
    for idx in xrange(nf-2, -1, -1):
        cell = F[idx,:]
        if np.all(cell == precell):
            print "Duplicate facet !!"
        precell = cell

    return F_backup


# =======================================================================
#                             MESH LOADERS
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

    V, F = merge_duplicates(V, F)
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

    ifile.readline() # skip one header line
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

    V, F = merge_duplicates(V, F)

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
    elif ext == '.tec':
        write_TEC(filename, V, F)
    elif ext == '.hst':
        write_HST(filename, V, F)
    else:
        raise RuntimeError, 'extension %s is not recognized' % ext


def write_HST(filename, V, F):
    """
    This function writes mesh into a HST file format
    :param filename:
    :param V:
    :param F:
    :return:
    """

    ls = os.linesep

    ofile = open(filename, 'w')

    ofile.write('PROJECT:' + ls)
    ofile.write('USERS:   meshmagick' + ls + ls)

    ofile.write('NBODY   1' + ls)
    ofile.write('RHO   1025.0' + ls)
    ofile.write('GRAVITY   9.81' + ls + ls)

    ofile.write('COORDINATES' + ls)
    idx = 0
    line = '%10u%16.6E%16.6E%16.6E' + ls
    for node in V:
        idx += 1
        ofile.write(line % (idx, node[0], node[1], node[2]))
    ofile.write('ENDCOORDINATES' + ls + ls)

    ofile.write('PANEL TYPE 0' + ls)
    line = '%10u%10u%10u%10u' + ls
    for elem in F:
        ofile.write(line % (elem[0], elem[1], elem[2], elem[3]))
    ofile.write('ENDPANEL' + ls + ls)

    ofile.write('ENDFILE' + ls)

    ofile.close()
    print 'File %s written' % filename


def write_TEC(filename, V, F):
    """
    This function writes data in a tecplot file

    :param filename:
    :param V:
    :param F:
    :return:
    """
    ofile = open(filename, 'w')

    nv = max(np.shape(V))
    nf = max(np.shape(F))

    ofile.write('TITLE = \" THIS FILE WAS GENERATED BY MESHMAGICK - FICHIER : %s \" \n' % filename)

    ofile.write('VARIABLES = \"X\",\"Y\",\"Z\" \n')
    ofile.write('ZONE T=\"MESH\" \n')
    ofile.write('N=%10u ,E=%10u ,F=FEPOINT, ET=QUADRILATERAL\n' % (nv, nf))

    for node in V:
        ofile.write('%16.6E%16.6E%16.6E\n' % (node[0], node[1], node[2]))
    for elem in F:
        ofile.write('%10u%10u%10u%10u\n' % (elem[0], elem[1], elem[2], elem[3]))

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
        if cell[0] == cell[-1]:
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
    # nf = max(np.shape(F))

    ofile = open(filename, 'w')

    ofile.write('%6u%6u\n' % (2, 0))  # TODO : mettre les symetries en argument
    idx = 0
    for vertex in V:
        idx += 1
        ofile.write('%6u%12.5f%12.5f%12.5f\n' % (idx, vertex[0], vertex[1],
                                                 vertex[2]))
    ofile.write('%6u%6u%6u%6u%6u\n' % (0, 0, 0, 0, 0))

    idx = 0
    for cell in F:
        idx += 1
        ofile.write('%6u%6u%6u%6u\n' % (cell[0], cell[1], cell[2], cell[3]))
    ofile.write('%6u%6u%6u%6u\n' % (0, 0, 0, 0))
    ofile.write('%6u\n' % 0)

    ofile.close()
    print 'File %s written' % filename


def write_STL(filename, V, F):
    # TODO : implement a STL writer based on vtk lib
    ofile = open(filename, 'w')

    ofile.write('solid meshmagick' + os.linesep)
    F -= 1
    for facet in F:
        if facet[0] != facet[3]:
            raise RuntimeError, 'Only full triangle meshes are accepted in STL files'

        # Computing normal
        V0 = V[facet[0], :]
        V1 = V[facet[1], :]
        V2 = V[facet[2], :]

        n = np.cross(V1 - V0, V2 - V0)
        n /= np.linalg.norm(n)

        ofile.write('  facet normal %15.6e%15.6e%15.6e' % (n[0], n[1], n[2]))
        ofile.write(os.linesep)
        ofile.write('    outer loop' + os.linesep)
        ofile.write('      vertex %15.6e%15.6e%15.6e' % (V0[0], V0[1], V0[2]))
        ofile.write(os.linesep)
        ofile.write('      vertex %15.6e%15.6e%15.6e' % (V1[0], V1[1], V1[2]))
        ofile.write(os.linesep)
        ofile.write('      vertex %15.6e%15.6e%15.6e' % (V2[0], V2[1], V2[2]))
        ofile.write(os.linesep)
        ofile.write('    endloop' + os.linesep)
        ofile.write('  endfacet' + os.linesep)

    ofile.write('endsolid meshmagick' + os.linesep)
    ofile.close()
    print 'File %s written' % filename


#=======================================================================
#                         MESH MANIPULATIONS
#=======================================================================
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


# =======================================================================
#                         COMMAND LINE USAGE
# =======================================================================
if __name__ == '__main__':
    import argparse

    try:
        import argcomplete

        acok = True
    except:
        acok = False

    parser = argparse.ArgumentParser(
        description="""A python module to manipulate meshes from different format used in hydrodynamics as well as for
                    visualization of these meshes into paraview and tecplot.

                    The formats actually supported are (R: reading; W: writing):

                         - mar      (R/W) : the format used by BEM NEMOH software (Ecole Centrale de Nantes)
                         - gdf      (R/W) : the format used by BEM WAMIT software (WAMIT, Inc.)
                         - inp      (R)   : the format used by BEM DIODORE software (PRINCIPIA)
                         - hst      (R/W) : the format used by BEM HYDROSTAR software (BV)
                         - nat      (R/W) : a natural format to store unstructured 2D meshes

                         - msh      (R)   : the file format from the GMSH mesher (C. Geuzaine, J.-F. Remacle)
                         - stl      (R/W) : broadly used file format for meshes, mostly generated from CAO softwares
                         - vtk, vtu (R/W) : file format for visualization in Paraview software (Kitware) (vtk is the old release and vtu is the new release)
                         - tec      (R/W) : file format for visualization in Tecplot

                    """,
        epilog="""Copyright 2014  - Francois Rongere\nEcole Centrale de Nantes""",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('infilename',
                        help='path of the input mesh file in any format')
    parser.add_argument('-o', '--outfilename',
                        help='path of the output mesh file. The format of ' +
                             'this file is determined from the extension given')

    parser.add_argument('-fmt', '--format', type=str,
                        help="""specifies the output file format and use the base name
                        of the infilename instead of the output file. File format has
                        to be chosen among VTK, VTU, MAR, GDF, STL, NAT, MHE. Default is VTU.
                        If -o option has been done, it has no effect""")
    parser.add_argument('-v', '--verbose',
                        help="""make the program give more informations on the computations""",
                        action='store_true')
    parser.add_argument('-i', '--info',
                        help="""extract informations on the mesh on the standard output""",
                        action='store_true')
    parser.add_argument('-b', '--binary',
                        help="""specifies wether the output file will be in a binary
                        format. It is only used when used with --convert option and
                        with .vtu output format""",
                        action='store_true')

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
    parser.add_argument('-u', '--unit',
                        type=str, choices=['rad', 'deg'],
                        default='deg',
                        help="""sets the unit for rotations. Default is deg""")
    parser.add_argument('-s', '--scale',
                        type=float,
                        help="""scales the mesh""")
    parser.add_argument('--flip-normals', action='store_true',
                        help="""flips the normals of the mesh""")
    parser.add_argument('--cut', action='store_true',
                        help="""cuts the mesh with the plane defined with the option
                        --plane""")
    parser.add_argument('--symmetrize', type=str,
                        nargs='*', default=None, choices=['x', 'y', 'z', 'p'],
                        help="""performs a symmetry of the mesh.
                        If x argument is given, then a symmetry with respect to plane
                        0yz etc... If p or no argument is given, then the --plane option
                        is used.""")
    parser.add_argument('-p', '--plane', nargs=4, type=float,
                        help="""Defines a plane used by the --cut option.
                        It is followed by the floats nx ny nz c where [nx, ny, nz]
                        is a normal vector to the plane and c defines its position
                        following the equation <N|X> = c with X a point belonging
                        to the plane""")
    parser.add_argument('--merge-duplicates', action='store_true',
                        help="""merges the duplicate nodes in the mesh""")
    parser.add_argument('--renumber', action='store_true',
                        help="""renumbers the cells and nodes of the mesh so as
                        to optimize cache efficiency in algorithms. Uses the Sloan
                        method""")
    parser.add_argument('--optimize', action='store_true',
                        help="""optimizes the mesh. Same as --merge-duplicates
                        and --renumber used together""")

    if acok:
        argcomplete.autocomplete(parser)

    args, unknown = parser.parse_known_args()

    extension_dict = ('vtk', 'vtu', 'gdf', 'mar', 'nat', 'stl',
                      'msh', 'inp', 'tec', 'hst', 'rad')

    write_file = False  # switch to decide if data should be written to outfilename

    _, ext = os.path.splitext(args.infilename)
    if ext[1:].lower() not in extension_dict:
        raise IOError, 'Extension "%s" is not known' % ext

    if args.outfilename is not None:
        _, ext = os.path.splitext(args.outfilename)
        write_file = True
        if ext[1:].lower() not in extension_dict:
            # Unknown extension
            raise IOError, 'Extension "%s" is not known' % ext
    else:
        if args.format is not None:
            write_file = True
            if args.format.lower() not in extension_dict:
                raise IOError, 'Extension ".%s" is not known' % args.format
            root, _ = os.path.splitext(args.infilename)
            args.outfilename = root + '.' + args.format
        else:
            args.outfilename = args.infilename

    # IMPORTING DATA FROM FILE
    try:
        V, F = load_mesh(args.infilename)
    except:
        raise IOError, "Can't open %s" % args.infilename

    # TESTING
    V, F = merge_duplicates(V, F, verbose=True)
    # F = merge_cells(F)
    myMesh = Mesh(V, F)
    # FIN TESTING

    # Dealing with different options
    if args.optimize:
        args.merge_duplicates = True
        args.renumber = True

    if args.merge_duplicates:
        try:
            from pymesh import pymesh
            pymesh.set_mesh_data(V, F)
            pymesh.cleandata_py()
            V = np.copy(pymesh.vv)
            F = np.copy(pymesh.ff)
            pymesh.free_mesh_data()
        except:
            V, F = merge_duplicates(V, F, args.verbose)

    if args.renumber:
        raise NotImplementedError, "Renumbering is not implemented yet"
        # try:
        #     from pymesh import pymesh
        # except:
        #     raise ImportError, 'pymesh module is not available, unable to use the --renumber option'
        #     ##raise NotImplementedError, '--renumber option is not implemented yet'

    if args.binary:
        raise NotImplementedError, 'Not implemented yet'

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

    def cast_angles(angles, unit):
        from math import pi

        if unit == 'deg':
            try:
                angles = map(float, angles)
            except:
                raise IOError, 'Bad input for rotation arguments'
            angles = np.array(angles) * pi / 180.
        else:
            def evalpi(elt):
                from math import pi

                elttemp = elt.replace('pi', '1.')
                try:
                    elttemp = eval(elttemp)
                except:
                    raise IOError, 'Bad input %s' % elt
                return float(eval(elt))

            angles = map(evalpi, angles)
        return angles

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

    if args.symmetrize is not None:
        raise NotImplementedError, "Symmetrization is not implemented yet into meshmagick"
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

    if args.cut:
        raise NotImplementedError, '--cut option is not implemented yet'

    if args.flip_normals:
        F = flip_normals(F)

    if args.info:
        get_info(V, F)

    if write_file:
        write_mesh(args.outfilename, V, F)

