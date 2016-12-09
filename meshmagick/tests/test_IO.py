#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from meshmagick.mmio import *
import os
from hashlib import md5

# Loaders
# ------------------

def test_load_INP():
    vertices, faces = load_INP('meshmagick/tests/data/Buoys.inp')
    mv = md5()
    mv.update(vertices)
    assert mv.digest() == '\xea\x97dv`\xf7F\xb1\xf85\x1e}\xd4F\x0fi'
    mf = md5()
    mf.update(faces)
    assert mf.digest() == '\x1e\xf7\xca\x98\xec\x19\x86\xb2\x1b\xc8\x87$\xa7\xf4\xb2\xbb'
    return
    

def test_load_MSH():
    vertices, faces = load_MSH('meshmagick/tests/data/SEAREV.msh')
    mv = md5()
    mv.update(vertices)
    assert mv.digest() == 'p\x97\xc0\xeb\x87\x08?\xd8\xbb\x9d\x9e\x88\x14I\x03\xfe'
    mf = md5()
    mf.update(faces)
    assert mf.digest() == 'x[p\xf4\x8b\x86*H\xf5\x06\xa8\xf2\xc7\x85\x98\xa7'
    return


def test_load_WRL():
    vertices, faces = load_WRL('meshmagick/tests/data/submarine.wrl')
    mv = md5()
    mv.update(vertices)
    assert mv.digest() == 'H\xbf\xc5=;-\xe0\xab\xa1\x8b\x04A\xe6\x99U\x96'
    mf = md5()
    mf.update(faces)
    assert mf.digest() == '\x15|\xda1\x8f\xe7]\xb733v\xf6\xdd.\x05\x0f'
    return


def test_load_MED():
    vertices, faces = load_MED('meshmagick/tests/data/barge.med')
    mv = md5()
    mv.update(vertices)
    assert mv.digest() =='\xea{p\xb9\xf6m\xe2z\x83C\xf9\rJ=\x0f\xba'
    mf = md5()
    mf.update(faces)
    assert mf.digest() =='\xb2\xa93\xf2\xf9>\xa4{jJ\xd8UbD\xe1\x88'
    return


def test_all_io():
    vertices, faces = load_VTP('meshmagick/tests/data/SEAREV.vtp')
    mv = md5()
    mv.update(vertices)
    hash_vertices = mv.digest()
    
    mf = md5()
    mf.update(faces)
    hash_faces = mv.digest()
    
    for (loader, writer) in extension_dict.values():
        try:
            writer('meshfile', vertices, faces)
            can_try_to_load = True
        except NotImplementedError:
            can_try_to_load = False
            print writer
        
        if can_try_to_load:
            try:
                loader('meshfile')
            except NotImplementedError:
                pass
    
    os.remove('meshfile')
        
