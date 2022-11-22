#!/usr/bin/env python
#  -*- coding: utf-8 -*-
from pathlib import Path

from meshmagick.mmio import *
from meshmagick.inertia import sphere
import os


def test_all_io():
    vertices, faces = load_VTP('meshmagick/tests/data/SEAREV.vtp')
    
    for (loader, writer) in list(extension_dict.values()):
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


def test_load_gdf_compressed():
    
    body = sphere(10)
    vertices, faces = load_VTP('meshmagick/tests/data/SEAREV.vtp')
    body_path = Path("temp_mesh.gdf")
    write_GDF(body_path, vertices, faces, ulen=1, gravity=9.81, isx=0, isy=0)

    gdf_vertices, gdf_faces = load_GDF(str(body_path))
    gdf_compressed_vertices, gdf_compressed_faces = load_GDF_compressed(str(body_path))
    
    np.testing.assert_allclose(
        gdf_vertices[gdf_faces], 
        gdf_compressed_vertices[gdf_compressed_faces]
        )

    body_path.unlink()
    