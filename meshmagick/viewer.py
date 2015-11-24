#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PYTHON_ARGCOMPLTETE_OK

"""
This module is part of meshmagick. It implements a viewer based on vtk
"""

import vtk

class mmViewer:
    def __init__(self):

        # Building renderer
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.7706, 0.8165, 1.0)

        # Building render window
        self.renWin = vtk.vtkRenderWindow()
        self.renWin.SetSize(1024, 769)
        self.renWin.SetWindowName("Meshmagick viewer")
        self.renWin.AddRenderer(self.renderer)

        # Building interactor
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.renWin)
        self.iren.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

        # Building axes view
        axes = vtk.vtkAxesActor()
        widget = vtk.vtkOrientationMarkerWidget()
        widget.SetOrientationMarker(axes)
        self.widget = widget

        # Building command annotations
        command_text = "left mouse : rotate\n" + \
                        "right mouse : zoom\n" + \
                        "middle mouse : pan\n" + \
                        "ctrl+left mouse : spin\n" + \
                        "f : focus on the mouse cursor\n" + \
                        "r : reset view\n" + \
                        "s : surface representation\n" + \
                        "w : wire representation\n" + \
                        "e : quit"

        corner_annotation = vtk.vtkCornerAnnotation()
        corner_annotation.SetLinearFontScaleFactor(2)
        corner_annotation.SetNonlinearFontScaleFactor(1)
        corner_annotation.SetMaximumFontSize(20)
        corner_annotation.SetText(3, command_text)
        corner_annotation.GetTextProperty().SetColor(0., 0., 0.)
        self.renderer.AddViewProp(corner_annotation)

        self.normals = False

    def normals_on(self):
        self.normals = True

    def normals_off(self):
        self.normals = False

    def add_polydata(self, polydata):
        if not isinstance(polydata, vtk.vtkPolyData):
            raise AssertionError, 'polydata must be a vtkPolyData object'

        # Building mapper
        mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.SetInput(polydata)
        else:
            mapper.SetInputData(polydata)

        # Building actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # Properties setting
        actor.GetProperty().SetColor(1, 1, 0)
        actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetEdgeColor(0, 0, 0)
        actor.GetProperty().SetLineWidth(1)

        self.renderer.AddActor(actor)

        if self.normals:
            normals = vtk.vtkPolyDataNormals()
            normals.ConsistencyOff()
            # normals.ComputePointNormalsOn()
            normals.ComputeCellNormalsOn()
            if vtk.VTK_MAJOR_VERSION <= 5:
                normals.SetInput(polydata)
            else:
                normals.SetInputData(polydata)
            # normals.SplittingOff()
            # normals.AutoOrientNormalsOn()
            normals.Update()

            normals_at_centers = vtk.vtkCellCenters()
            normals_at_centers.SetInputConnection(normals.GetOutputPort())

            normals_mapper = vtk.vtkPolyDataMapper()
            if vtk.VTK_MAJOR_VERSION <= 5:
                normals_output = normals.GetOutput()
                normals_mapper.SetInput(normals_output)
            else:
                normals_mapper.SetInputData(normals.GetOutput())
            normals_actor = vtk.vtkActor()
            normals_actor.SetMapper(normals_mapper)

            arrows = vtk.vtkArrowSource()
            arrows.SetTipResolution(16)
            arrows.SetTipLength(0.3)
            arrows.SetTipRadius(0.1)

            glyph = vtk.vtkGlyph3D()
            glyph.SetSourceConnection(arrows.GetOutputPort())
            glyph.SetInputConnection(normals_at_centers.GetOutputPort())
            glyph.SetVectorModeToUseNormal()
            glyph.Update()

            glyph_mapper = vtk.vtkPolyDataMapper()
            glyph_mapper.SetInputConnection(glyph.GetOutputPort())

            glyph_actor = vtk.vtkActor()
            glyph_actor.SetMapper(glyph_mapper)

            self.renderer.AddActor(glyph_actor)

    def show(self):

        self.widget.SetInteractor(self.iren)
        self.widget.SetEnabled(1)
        self.widget.InteractiveOn()

        self.iren.Initialize()
        self.renWin.Render()
        self.iren.Start()

