#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from meshmagick.mmio import *
import os

# Loaders
# ------------------

def test_load_GDF():
    load_GDF('meshmagick/tests/data/coque.gdf')

def test_load_HST():
    load_HST('meshmagick/tests/data/SEAREV.hst')

def test_load_MAR():
    load_MAR('meshmagick/tests/data/SEAREV.mar')

def test_load_TEC():
    load_TEC('meshmagick/tests/data/SEAREV.tec')

def test_load_INP():
    load_INP('meshmagick/tests/data/Buoys.inp')

def test_load_VTU():
    load_VTU('meshmagick/tests/data/SEAREV.vtu')

def test_load_VTP():
    load_VTP('meshmagick/tests/data/SEAREV.vtp')

def test_load_MSH():
    load_MSH('meshmagick/tests/data/SEAREV.msh')

def test_load_STL():
    load_STL('meshmagick/tests/data/SEAREV.stl')

def test_all_io():
    vertices, faces = load_VTP('meshmagick/tests/data/SEAREV.vtp')
    for (loader, writer) in extension_dict.values():
        try:
            print writer
            writer('meshfile', vertices, faces)
            can_try_to_load = True
        except NotImplementedError:
            can_try_to_load = False
        
        if can_try_to_load:
            try:
                print loader
                loader('meshfile')
            except NotImplementedError:
                pass
    
    os.remove('meshfile')
        


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
