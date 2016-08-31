#!/usr/bin/env python
#  -*- coding: utf-8 -*-

# import pytest
from mmio import *
from mesh import Mesh


# Loaders
# ------------------

def test_load_GDF():
    load_GDF('tests/data/coque.gdf')

def test_load_HST():
    load_HST('tests/data/SEAREV.hst')

def test_load_MAR():
    load_MAR('tests/data/SEAREV.mar')

def test_load_TEC():
    load_TEC('tests/data/SEAREV.tec')

def test_load_INP():
    load_INP('tests/data/Buoys.inp')

def test_load_VTU():
    load_VTU('tests/data/SEAREV.vtu')

def test_load_VTP():
    load_VTP('tests/data/SEAREV.vtp')

def test_load_MSH():
    load_MSH('tests/data/SEAREV.msh')

def test_load_STL():
    load_STL('tests/data/SEAREV.stl')
    
# def test_load_write_all():
#     vertices, faces = load_VTP('tests/data/SEAREV.vtp')
#
#     for (loader, writer) in extension_dict.values():
#         try:
#             writer('tests/data/temp.dat', vertices, faces)
#             loader('tests/data/temp.dat')
#         except NotImplementedError:
#             print 'NotImplemented'
#             pass
#         else:
#             print 'error'
#             pass
