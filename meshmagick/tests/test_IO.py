#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from meshmagick.mmio import *
import os


def test_all_io():
    vertices, faces = load_VTP('meshmagick/tests/data/SEAREV.vtp')
    
    for (loader, writer) in extension_dict.values():
        try:
            writer('meshfile', vertices, faces)
            can_try_to_load = True
        except NotImplementedError:
            can_try_to_load = False
        
        if can_try_to_load:
            try:
                loader('meshfile')
            except NotImplementedError:
                pass
    
    os.remove('meshfile')
