#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import os
import numpy as np

# =======================================================================
# MESH LOADERS
# ======================================================================
# Contains here all functions to load meshes from different file formats

def load_mesh(filename, format):
    """load_mesh(filename, format)

    Driver function that loads every mesh file format known by meshmagick
    and returns the node list and the connectivity array

    Parameters:
        filename: str
            name of the meh file on disk
        format: str
            format of the mesh defined in the extension_dict dictionary

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities
    """
    os.path.isfile(filename)

    if not extension_dict.has_key(format):
        raise IOError, 'Extension ".%s" is not known' % format

    loader = extension_dict[format][0]

    V, F = loader(filename)

    return V, F

def load_RAD(filename):
    """load_RAD(filename)

    Loads RADIOSS mesh files. This export file format may be chosen in ICEM
    meshing program.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: RAD files have a 1-indexing
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
    V = np.asarray(V, dtype=float)

    F = []
    elem_section = pattern_elem_section.search(data).group(1)
    for elem in pattern_elem_line.findall(elem_section):
        F.append(map(int, elem.strip().split()[3:]))
    F = np.asarray(F, dtype=np.int) - 1

    return V, F

def load_HST(filename):
    """load_HST(filename)

    Loads HYDROSTAR (Bureau Veritas (c)) mesh files.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: HST files have a 1-indexing
    """
    ifile = open(filename, 'r')
    data = ifile.read()
    ifile.close()

    import re

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
        Vtmp = np.asarray(Vtmp, dtype=np.float)
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
        Ftmp = np.asarray(Ftmp, dtype=np.int)
        if nf == 0:
            F = Ftmp.copy()
            nf = nftmp
        else:
            F = np.concatenate((F, Ftmp))
            nf += nftmp

    return V, F-1

def load_DAT(filename):
    """Not implemented.
    Intended to load .DAT files used in DIODORE (PRINCIPIA (c))
    """
    raise NotImplementedError

def load_INP(filename):
    """load_INP(filename)

    Loads DIODORE (PRINCIPIA (c)) configuration file format. It parses
    the .INP file and extract meshes defined in subsequent .DAT files
    using the different informations contained in the .INP file.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: INP/DAT files use a 1-indexing
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
        frames[framename] = np.asarray(map(float, framevector))

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
            meshfile = open(os.path.join(os.path.dirname(filename), file + '.DAT'), 'r')
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
        id_new = - np.ones(max(idx_array) + 1, dtype=np.int)
        # FIXME: cette partie est tres buggee !!!
        for i, idx in enumerate(idx_array):
            id_new[idx] = i+1

        meshfiles[file]['ELEM_SECTIONS'] = []
        for elem_section in pattern_elem_section.findall(data):

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
            nodes = np.asarray(meshfiles[file]['NODE_SECTION'], dtype=np.float)
            # Translation of nodes according to frame option id any
            nodes = translate(nodes, frames[field['FRAME']])  # TODO: s'assurer que frame est une options obligatoire...

            if nbNodes == 0:
                V = nodes.copy()
                nbNodes = V.shape[0]
                increment = False
                continue

            if field['INCREMENT'] == 'NO':
                V[idx, :] = nodes.copy()
                increment = False
            else:
                V = np.concatenate((V, nodes))
                nbNodes = V.shape[0]
                increment = True
        else:  # this is an ELEMENT section
            elem_section = np.asarray(meshfiles[file]['ELEM_SECTIONS'][meshfiles[file]['nb_elem_sections_used']],
                                      dtype=np.int)

            meshfiles[file]['nb_elem_sections_used'] += 1
            if meshfiles[file]['nb_elem_sections_used'] == meshfiles[file]['nb_elem_sections']:
                meshfiles[file]['nb_elem_sections_used'] = 0

            # Updating to new id of nodes
            elems = elem_section
            if increment:
                elems += nbNodes

            if nbElems == 0:
                F = elems.copy()
                nbElems = F.shape[0]
                continue
            else:
                F = np.concatenate((F, elems))
                nbElems = F.shape[0]

    return V, F-1

def load_TEC(filename):
    """load_TEC(filename)

    Loads TECPLOT (Tecplot (c)) mesh files. It relies on the tecplot file
    reader from the VTK library.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: TEC files have a 0-indexing
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

        Vtmp = np.zeros((nvblock, 3), dtype=np.float)
        for k in range(nvblock):
            Vtmp[k] = np.array(block.GetPoint(k))

        if nv == 0:
            V = Vtmp
        else:
            V = np.concatenate((V, Vtmp))

        nv += nvblock

        # Facet extraction
        Ftmp = np.zeros((nfblock, 4), dtype=np.int)
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

    return V, F

def load_VTU(filename):
    """load_VTU(filename)

    Loads VTK file format in the new XML format (vtu file extension for
    unstructured meshes). It relies on the reader from the VTK library.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: VTU files have a 0-indexing
    """

    from vtk import vtkXMLUnstructuredGridReader
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    V, F = _dump_vtk(vtk_mesh)
    return V, F

def load_VTP(filename):
    """load_VTP(filename)

    Loads VTK file format in the new XML format (vtp file extension for
    polydata meshes). It relies on the reader from the VTK library.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: VTP files have a 0-indexing
    """

    from vtk import vtkXMLPolyDataReader
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    V, F = _dump_vtk(vtk_mesh)
    return V, F

def load_VTK(filename):
    """load_VTK(filename)

    Loads VTK file format in the legacy format (vtk file extension).
    It relies on the reader from the VTK library.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: VTU files have a 0-indexing
    """
    from vtk import vtkUnstructuredGridReader
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    V, F = _dump_vtk(vtk_mesh)
    return V, F

def _dump_vtk(vtk_mesh):
    """_dump_vtk(vtk_mesh)

    Internal driver function that uses the VTK library to read VTK polydata
    or vtk unstructured grid data structures

    Parameters:
        filename: str
            name of the meh file on disk
        reader: Reader
            the reader to use (new XML format ot legacy vtk format)

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities
    """
    # Importing the mesh from the file
    # reader.SetFileName(filename)
    # reader.Update()
    # vtk_mesh = reader.GetOutput()

    nv = vtk_mesh.GetNumberOfPoints()
    V = np.zeros((nv, 3), dtype=np.float)
    for k in range(nv):
        V[k] = np.array(vtk_mesh.GetPoint(k))

    nf = vtk_mesh.GetNumberOfCells()
    F = np.zeros((nf, 4), dtype=np.int)
    for k in range(nf):
        cell = vtk_mesh.GetCell(k)
        nv_facet = cell.GetNumberOfPoints()
        for l in range(nv_facet):
            F[k][l] = cell.GetPointId(l)
        if nv_facet == 3:
            F[k][3] = F[k][0]

    return V, F

def load_STL(filename):
    """load_STL(filename)

    Loads STL file format. It relies on the reader from the VTK library.
    As STL file format maintains a redundant set of vertices for each faces
    of the mesh, it returns a merged list of nodes and connectivity array
    by using the merge_duplicates function.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: STL files have a 0-indexing
    """

    from vtk import vtkSTLReader

    reader = vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()

    data = reader.GetOutputDataObject(0)

    nv = data.GetNumberOfPoints()
    V = np.zeros((nv, 3), dtype=np.float)
    for k in range(nv):
        V[k] = np.array(data.GetPoint(k))
    nf = data.GetNumberOfCells()
    F = np.zeros((nf, 4), dtype=np.int)
    for k in range(nf):
        cell = data.GetCell(k)
        if cell is not None:
            for l in range(3):
                F[k][l] = cell.GetPointId(l)
                F[k][3] = F[k][0]  # always repeating the first node as stl is triangle only

    # Merging duplicates nodes
    V, F = merge_duplicates(V, F)

    return V, F

def load_NAT(filename):
    """load_NAT(filename)

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


    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: NAT files have a 1-indexing
    """

    ifile = open(filename, 'r')
    xsym, ysym = map(int, ifile.readline().split())
    nv, nf = map(int, ifile.readline().split())

    V = []
    for i in range(nv):
        V.append(map(float, ifile.readline().split()))
    V = np.array(V, dtype=np.float)

    F = []
    for i in range(nf):
        F.append(map(int, ifile.readline().split()))
    F = np.array(F, dtype=np.int)

    ifile.close()
    return V, F-1

def load_GDF(filename):
    """load_GDF(filename)

    Loads WAMIT (Wamit INC. (c)) GDF mesh files. As GDF file format maintains
    a redundant set of vertices for each faces of the mesh, it returns a merged
    list of nodes and connectivity array by using the merge_duplicates function.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: GDF files have a 1-indexing
    """
    ifile = open(filename, 'r')

    ifile.readline()  # skip one header line
    line = ifile.readline().split()
    ulen = line[0]
    grav = line[1]

    line = ifile.readline().split()
    isx = line[0]
    isy = line[1]

    line = ifile.readline().split()
    nf = int(line[0])

    V = np.zeros((4 * nf, 3), dtype=np.float)
    F = np.zeros((nf, 4), dtype=np.int)

    iv = -1
    for icell in range(nf):

        for k in range(4):
            iv += 1
            V[iv, :] = np.array(ifile.readline().split())
            F[icell, k] = iv

    ifile.close()
    V, F = merge_duplicates(V, F, verbose=True)

    return V, F

def load_MAR(filename):
    """load_MAR(filename)

    Loads Nemoh (Ecole Centrale de Nantes) mesh files.

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: MAR files have a 1-indexing
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

    V = np.array(V, dtype=np.float)
    F = []
    while 1:
        line = ifile.readline()
        line = line.split()
        if line[0] == '0':
            break
        F.append(map(int, line))

    F = np.array(F, dtype=np.int)

    ifile.close()

    return V, F-1

def load_MSH(filename):
    """load_MSH(filename)

    Loads GMSH .MSH mesh files. It curretly uses an external module but should rely
    on a reader written for meshmagick (or on meshpy module but this one is hard to
    install on Windows...)

    Parameters:
        filename: str
            name of the meh file on disk

    Returns:
        vertices: ndarray
            numpy array of the coordinates of the mesh's nodes
        faces: ndarray
            numpy array of the faces' nodes connectivities

    Note: MSH files have a 0-indexing
    """
    import gmsh

    myMesh = gmsh.Mesh()
    myMesh.read_msh(filename)
    V = np.array(myMesh.Verts, dtype=np.float)

    ntri = myMesh.nElmts.get(2)
    nquad = myMesh.nElmts.get(3)
    if ntri is None:
        ntri = 0
    if nquad is None:
        nquad = 0

    nel = ntri + nquad

    F = np.zeros((nel, 4), dtype=np.int)

    if ntri != 0:
        F[:ntri, :3] = myMesh.Elmts.get(2)[1]
        F[:, 3] = F[:, 0]

    if nquad != 0:
        F[ntri:, :] = myMesh.Elmts.get(3)[1]

    return V, F
# def load_STL2(filename):
#     import re
#
#     ifile = open(filename, 'r')
#     text = ifile.read()
#     ifile.close()
#
#     endl = r'(?:\n|\r|\r\n)'
#     patt_str = r"""
#             ^\s*facet\s+normal(.*)""" + endl + """
#             ^\s*outer\sloop""" + endl + """
#             ^\s*vertex\s+(.*)""" + endl + """
#             ^\s*vertex\s+(.*)""" + endl + """
#             ^\s*vertex\s+(.*)""" + endl + """
#             ^\s*endloop""" + endl + """
#             ^\s*endfacet""" + endl + """
#            """
#     pattern = re.compile(patt_str, re.MULTILINE | re.VERBOSE)
#
#     normal = []
#     vertices = []
#     for match in pattern.finditer(text):
#         normal.append(map(float, match.group(1).split()))
#         vertices.append(map(float, match.group(2).split()))
#         vertices.append(map(float, match.group(3).split()))
#         vertices.append(map(float, match.group(4).split()))
#
#     vertices = np.array(vertices, dtype=float, order='fortran')
#
#     nf = np.size(vertices, 0) / 3
#     faces = np.zeros((nf, 4), dtype=np.int32, order='fortran')
#
#     base = np.array([1, 2, 3, 1])
#     for i in range(nf):
#         faces[i, :] = base + 3 * i
#
#     return vertices, faces





#=======================================================================
#                             MESH WRITERS
#=======================================================================

def write_mesh(filename, V, F, format):
    """write_mesh(filename, format)

    Driver function that writes every mesh file format known by meshmagick

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        format: str
            format of the mesh defined in the extension_dict dictionary

    """

    if not extension_dict.has_key(format):
        raise IOError, 'Extension "%s" is not known' % format

    writer = extension_dict[format][1]

    writer(filename, V, F)

    return 1

def write_DAT(filename, V, F):
    """write_DAT(filename, vertices, faces)

    Writes .DAT file format for the DIODORE (PRINCIPA (c)) software.
    It also displays suggestions for inclusion into the .INP configuration
    file.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

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
    nq = 0
    nt = 0
    for (idx, cell) in enumerate(F+1):
        if cell[0] != cell[-1]:
            # quadrangle
            nq += 1
            quad_block = ''.join(
                (quad_block,
                 '\n',
                 '{:8d}'.format(idx+1),
                 ''.join('{:8d}'.format(node_id) for node_id in cell)
                )
            )

        else:
            # Triangle
            nt += 1
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

    if nq > 0:
        quad_block = ''.join((quad_block, '\n*RETURN\n'))
        ofile.write(quad_block)
        print '*ELEMENT,TYPE=Q4C000,ELSTRUCTURE={0},INPUT={0}'.format(rootfilename)
    if nt > 0:
        tri_block = ''.join((tri_block, '\n*RETURN\n'))
        ofile.write(tri_block)
        print '*ELEMENT,TYPE=T3C000,ELSTRUCTURE={0},INPUT={0}'.format(rootfilename)

    print ''
    print '-------------------------------------------------'
    ofile.close()

    return 1

def write_HST(filename, V, F):
    """write_HST(filename, vertices, faces)

    Writes .HST file format for the HYDROSTAR (Bureau Veritas (c)) software.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

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
            ) for cell in F+1
        ),
        '\nENDPANEL\n\n'
    ))

    ofile.write(cells_coordinates)

    ofile.write('ENDFILE\n')

    ofile.close()

    print u'File {0:s} written'.format(filename)

def write_TEC(filename, V, F):
    """write_TEC(filename, vertices, faces)

    Writes .TEC file format for the TECPLOT (Tecplot (c)) visualisation
    software. It relies on the VTK library for its writer.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    ofile = open(filename, 'w')

    nv = V.shape[0]
    nf = F.shape[0]

    ofile.write('TITLE = \" THIS FILE WAS GENERATED BY MESHMAGICK - FICHIER : {} \" \n'.format(filename))

    ofile.write('VARIABLES = \"X\",\"Y\",\"Z\" \n')
    ofile.write('ZONE T=\"MESH\" \n')
    ofile.write('N={nv:10d} ,E={nf:10d} ,faces=FEPOINT, ET=QUADRILATERAL\n'.format(nv=nv, nf=nf))

    node_block = '\n'.join( # block
        ''.join(
            ''.join('{:16.6E}'.format(elt) for elt in node)
        ) for node in V
    ) + '\n'
    ofile.write(node_block)

    cells_block = '\n'.join(  # block
        ''.join(
            ''.join('{:10d}'.format(node_id) for node_id in cell)
        ) for cell in F+1
    ) + '\n'
    ofile.write(cells_block)

    ofile.close()

    return 1

def write_VTU(filename, V, F):
    """write_VTU(filename, vertices, faces)

    Writes .vtu file format for the paraview (Kitware (c)) visualisation
    software. It relies on the VTK library for its writer. VTU files use
    the last XML file format of the VTK library.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    from vtk import vtkXMLUnstructuredGridWriter
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetDataModeToAscii()
    writer.SetFileName(filename)

    unstructured_grid = _build_vtkUnstructuredGrid(V, F)
    writer.SetInput(unstructured_grid)
    writer.Write()

    return 1

def write_VTP(filename, V, F):
    """write_VTP(filename, vertices, faces)

    Writes .vtp file format for the paraview (Kitware (c)) visualisation
    software. It relies on the VTK library for its writer. VTP files use
    the last XML file format of the VTK library and correspond to polydata.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    from vtk import vtkXMLPolyDataWriter
    writer = vtkXMLPolyDataWriter()
    writer.SetDataModeToAscii()
    writer.SetFileName(filename)

    polydata = _build_vtkPolyData(V, F)
    writer.SetInput(polydata)
    writer.Write()

    return 1

def write_VTK(filename, V, F):
    """write_VTK(filename, vertices, faces)

    Writes .vtk file format for the paraview (Kitware (c)) visualisation
    software. It relies on the VTK library for its writer. VTK files use
    the legagy ASCII file format of the VTK library.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    from vtk import vtkUnstructuredGridWriter
    writer = vtkUnstructuredGridWriter()
    writer.SetFileName(filename)

    unstructured_grid = _build_vtkUnstructuredGrid(V, F)
    writer.SetInput(unstructured_grid)
    writer.Write()

    return 1

def _write_paraview(filename, V, F, writer):
    """_write_paraview(filename, vertices, faces)

    Internal driver function that writes vtk files to be visualised into
    the Paraview software. It relies on the VTK library.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities
        writer: Writer
            The writer to be used

    """
    writer.SetFileName(filename)
    vtk_mesh = _build_vtkUnstructuredGrid(V, F)
    writer.SetInput(vtk_mesh)
    writer.Write()

    return 1

def _build_vtkUnstructuredGrid(V, F):
    """_build_vtk_mesh_obj(vertices, faces)

    Internal function that builds a VTK object for manipulation by the VTK library.

    Parameters:
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    Returns: vtkObject
        the vtk object instance
    """
    import vtk

    nv = max(np.shape(V))
    nf = max(np.shape(F))

    vtk_mesh = vtk.vtkUnstructuredGrid()
    vtk_mesh.Allocate(nf, nf)

    # Building the vtkPoints data structure
    vtk_points = vtk.vtkPoints()
    vtk_points.SetNumberOfPoints(nv)
    for idx, vertex in enumerate(V):
        vtk_points.SetPoint(idx, vertex)

    vtk_mesh.SetPoints(vtk_points)  # Storing the points into vtk_mesh

    # Building the vtkCell data structure
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

def _build_vtkPolyData(V, F):
    import vtk

    # Create a vtkPoints object and store the points in it
    points = vtk.vtkPoints()
    for point in V:
        points.InsertNextPoint(point)

    # Create a vtkCellArray to store faces
    faces = vtk.vtkCellArray()
    for face_ids in F:
        if face_ids[0] == face_ids[-1]:
            # Triangle
            curface = face_ids[:3]
            vtk_face = vtk.vtkTriangle()
        else:
            # Quadrangle
            curface = face_ids[:4]
            vtk_face = vtk.vtkQuad()

        for idx, id in enumerate(curface):
            vtk_face.GetPointIds().SetId(idx, id)

        faces.InsertNextCell(vtk_face)

    polyDataMesh = vtk.vtkPolyData()
    polyDataMesh.SetPoints(points)
    polyDataMesh.SetPolys(faces)

    return polyDataMesh

def write_NAT(filename, V, F):
    """write_NAT(filename, vertices, faces)

    Writes .nat file format as defined into the load_NAT function.

    See:
        load_NAT

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """
    ofile = open(filename, 'w')

    nv = max(np.shape(V))
    nf = max(np.shape(F))

    ofile.write('%6u%6u\n' % (0, 0))  # lire les symmetries dans args...
    ofile.write('%6u%6u\n' % (nv, nf))
    for vertex in V:
        ofile.write('%15.6E%15.6E%15.6E\n' % (vertex[0], vertex[1], vertex[2]))
    for cell in F+1:
        ofile.write('%10u%10u%10u%10u\n' % (cell[0], cell[1], cell[2], cell[3]))

    ofile.close()

    return 1

def write_GDF(filename, V, F):
    """write_GDF(filename, vertices, faces)

    Writes .gdf file format for the WAMIT (Wamit INC. (c)) BEM software.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    nf = max(np.shape(F))

    ofile = open(filename, 'w')

    ofile.write('GDF file generated by meshmagick\n')

    ofile.write('%16.6f%16.6f\n' % (100.0, 9.81))
    ofile.write('%12u%12u\n' % (0, 1))  # TODO : mettre les symetries en argument
    ofile.write('%12u\n' % nf)

    for cell in F:
        for k in range(4):
            Vcur = V[cell[k], :]
            ofile.write('%16.6E%16.6E%16.6E\n' % (Vcur[0], Vcur[1], Vcur[2]))

    ofile.close()

    return 1

def write_MAR(filename, V, F):
    """write_MAR(filename, vertices, faces)

    Writes mesh files to be used with Nemoh BEM software (Ecole Centrale de Nantes)

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    # TODO: detecter la symetrie de plan Oxz

    ofile = open(filename, 'w')

    ofile.write('{0:6d}{1:6d}\n'.format(2, 0))  # TODO : mettre les symetries en argument

    nv = V.shape[0]
    for (idx, vertex) in enumerate(V):
        ofile.write('{0:6d}{1:16.6f}{2:16.6f}{3:16.6f}\n'.format(idx+1, vertex[0], vertex[1], vertex[2]))

    ofile.write('{0:6d}{1:6d}{2:6d}{3:6d}{4:6d}\n'.format(0, 0, 0, 0, 0))

    cell_block = '\n'.join(
        ''.join(u'{0:10d}'.format(elt) for elt in cell)
        for cell in F+1
    ) + '\n'
    ofile.write(cell_block)
    ofile.write('%6u%6u%6u%6u\n' % (0, 0, 0, 0))

    ofile.close()

    print 'WARNING: if you described only one part of the mesh using symmetry for Nemoh, you may manually modify the ' \
          'file header accordingly'

    return 1

def write_RAD(filename, V, F):
    raise NotImplementedError

def write_STL(filename, V, F):
    """write_STL(filename, vertices, faces)

    Writes .stl file format. It relies on the VTK library for its writer.

    Parameters:
        filename: str
            name of the mesh file to be written on disk
        V: ndarray
            numpy array of the coordinates of the mesh's nodes
        F: ndarray
            numpy array of the faces' nodes connectivities

    """

    # TODO : replace this implementation by using the vtk functionalities

    # Triangulating quads
    F = triangulate_quadrangles(F, verbose=False)

    ofile = open(filename, 'w')

    ofile.write('solid meshmagick\n')

    for facet in F:
        if facet[0] != facet[3]:
            raise RuntimeError, """Only full triangle meshes are accepted in STL files.
              Please consider using the --triangulate-quadrangles option (-tq) to
              perform a prior triangulation of the mesh"""

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

    return 1

def write_INP(filename, V, F):
    raise NotImplementedError

def write_MSH(filename, V, F):
    raise NotImplementedError


def know_extension(ext):
    return extension_dict.has_key(ext)

extension_dict = {  # keyword           reader,   writer
    'mar': (load_MAR, write_MAR),
    'nemoh': (load_MAR, write_MAR),
    'wamit': (load_GDF, write_GDF),
    'gdf': (load_GDF, write_GDF),
    'diodore-inp': (load_INP, write_INP),
    'inp': (load_INP, write_INP),
    'diodore-dat': (load_DAT, write_DAT),
    'hydrostar': (load_HST, write_HST),
    'hst': (load_HST, write_HST),
    'natural': (load_NAT, write_NAT),
    'nat': (load_NAT, write_NAT),
    'gmsh': (load_MSH, write_MSH),
    'msh': (load_MSH, write_MSH),
    'rad': (load_RAD, write_RAD),
    'radioss': (load_RAD, write_RAD),
    'stl': (load_STL, write_STL),  # FIXME: Verifier que ce n'est pas load_STL2
    'vtu': (load_VTU, write_VTU),
    'vtp': (load_VTP, write_VTP),
    'paraview-legacy': (load_VTK, write_VTK),  # VTK
    'vtk': (load_VTK, write_VTK),
    'tecplot': (load_TEC, write_TEC),
    'tec': (load_TEC, write_TEC)
}
