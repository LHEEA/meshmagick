#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module is part of meshmagick. It implements a viewer based on vtk
"""

import vtk
from os import getcwd
from datetime import datetime
from math import degrees
import numpy as np
import re

from multiprocessing import Process


try:
    import zmq
    has_zmq = True
except ImportError:
    has_zmq = False


__year__ = datetime.now().year

# TODO: See links below for interactive update of actors
# https://stackoverflow.com/questions/31075569/vtk-rotate-actor-programmatically-while-vtkrenderwindowinteractor-is-active
# https://stackoverflow.com/questions/32417197/vtk-update-vtkpolydata-in-renderwindow-without-blocking-interaction

def animation(mymesh, computations_port=6000):
    if not has_zmq:
        raise ImportError("zeromq module must be installed to use this functionality")
    
    visualization_port = computations_port + 1 # TODO: voir si c'est malin...
    
    # Starting the broker in a process
    Process(target=_visualization_broker, args=(computations_port, visualization_port)).start()
    
    # Starting the visualization window
    Process(target=_visualization_window, args=(mymesh, visualization_port)).start()
    
    
    return computations_port
    


def _visualization_broker(computations_port, visualization_port):
    from Queue import Queue, Empty
    
    context = zmq.Context()
    
    frontend = context.socket(zmq.REP)
    frontend.bind("tcp://127.0.0.1:%u" % visualization_port)
    
    backend = context.socket(zmq.SUB)
    backend.bind("tcp://127.0.0.1:%u" % computations_port)
    backend.setsockopt(zmq.SUBSCRIBE, 'viz')
    
    poller = zmq.Poller()
    poller.register(frontend, zmq.POLLIN)
    poller.register(backend, zmq.POLLIN)
    
    data_queue = Queue()
    
    while True:
        try:
            socks = dict(poller.poll())
        except KeyboardInterrupt:
            break
            
        # REP/REQ pattern for visualization frontend
        if frontend in socks:
            msg = frontend.recv()
            try:
                reply = data_queue.get(block=False)
            except Empty:
                reply = 'Empty'
            frontend.send(reply)
            
        # PUB/SUB pattern for computations backend
        if backend in socks:
            msg = ' '.join(backend.recv().split()[1:])
            data_queue.put(msg)
    
    
    


# def _visualization_window(mymesh, port):
#     print 'Creation visu dans un autre processus'
#     polydata = mymesh._vtk_polydata()
#
#     viewer = MMViewer()
#
#     transform = vtk.vtkTransform()
#     transform.Identity()
#
#     transform_filter = vtk.vtkTransformPolyDataFilter()
#     transform_filter.SetInputData(polydata) # Attention compatibilite versions VTK (VTK_MAJOR_VERSION)
#     transform_filter.SetTransform(transform)
#
#     mapper = vtk.vtkPolyDataMapper()
#     mapper.SetInputConnection(transform_filter.GetOutputPort()) # Voir si on est obliges de faire ca...
#
#     actor = vtk.vtkActor()
#     actor.SetMapper(mapper)
#
#     viewer.renderer.AddActor(actor)
#
#     cb = vtkTimerCallback(port, transform)
#
#     viewer.render_window_interactor.Initialize()
#     viewer.render_window_interactor.CreateRepeatingTimer(100)
#     viewer.render_window_interactor.AddObserver('TimerEvent', cb.execute)
#
#     viewer.show()
#     viewer.finalize()
#     return

def _visualization_window(mymesh, port):
    polydata = mymesh._vtk_polydata()
    viewer = MMViewer()
    pd_id = viewer.add_polydata(polydata)
    
    cb = vtkTimerCallback(port, viewer.transforms[pd_id])
    
    viewer.render_window_interactor.Initialize()
    viewer.render_window_interactor.CreateRepeatingTimer(100)
    viewer.render_window_interactor.AddObserver('TimerEvent', cb.execute)
    
    viewer.show()
    viewer.finalize()
    # TODO: envoyer un signal au broker pour qu'il se stoppe aussi ... --> PUBLISH
    return
    


class vtkTimerCallback(object):
    def __init__(self, port, transform):
        context = zmq.Context()
        socket = context.socket(zmq.REQ)
        socket.connect("tcp://127.0.0.1:%u" % port)
        
        self.socket = socket
        
        self.transform = transform
        
        
        
        
    
    def execute(self, obj, event):
        # print 'waiting message'
        
        self.socket.send('')
        msg = self.socket.recv()
        
        if msg == 'Empty':
            return
        
        # print 'VISU has received %s' % msg
        topic, tx, ty, tz = msg.split()
        # print tx, ty, tz
        # theta = np.array([tx, ty, 0.], dtype=np.float)
        #
        # angle = np.linalg.norm(theta)
        # ux, uy, uz = theta / angle
        # angle = degrees(angle)
        #
        # self.transform.RotateWXYZ(angle, ux, uy, uz)
        
        self.transform.Translate(map(float, [tx, ty, tz]))
        
        obj.GetRenderWindow().Render()
        return


class MMViewer:
    def __init__(self):

        # Building renderer
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.7706, 0.8165, 1.0)

        # Building render window
        self.render_window = vtk.vtkRenderWindow()
        self.render_window.SetSize(1024, 768)
        self.render_window.SetWindowName("Meshmagick viewer")
        self.render_window.AddRenderer(self.renderer)

        # Building interactor
        self.render_window_interactor = vtk.vtkRenderWindowInteractor()
        self.render_window_interactor.SetRenderWindow(self.render_window)
        self.render_window_interactor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
        self.render_window_interactor.AddObserver('KeyPressEvent', self.on_key_press, 0.0)

        # Building axes view
        axes = vtk.vtkAxesActor()
        widget = vtk.vtkOrientationMarkerWidget()
        widget.SetOrientationMarker(axes)
        self.widget = widget

        self.widget.SetInteractor(self.render_window_interactor)
        self.widget.SetEnabled(1)
        self.widget.InteractiveOn()

        # Building command annotations
        command_text = "left mouse : rotate\n" + \
                        "right mouse : zoom\n" + \
                        "middle mouse : pan\n" + \
                        "ctrl+left mouse : spin\n" + \
                        "n : (un)show normals\n" + \
                        "b : (un)show axes box\n" + \
                        "f : focus on the mouse cursor\n" + \
                        "r : reset view\n" + \
                        "s : surface representation\n" + \
                        "w : wire representation\n" + \
                        "h : (un)show Oxy plane\n" + \
                        "x : save\n" + \
                        "c : screenshot\n" + \
                        "q : quit"

        corner_annotation = vtk.vtkCornerAnnotation()
        corner_annotation.SetLinearFontScaleFactor(2)
        corner_annotation.SetNonlinearFontScaleFactor(1)
        corner_annotation.SetMaximumFontSize(20)
        corner_annotation.SetText(3, command_text)
        corner_annotation.GetTextProperty().SetColor(0., 0., 0.)
        self.renderer.AddViewProp(corner_annotation)
        
        # Building copyright text annotations
        copyright_text = "Meshmagick Viewer\nCopyright 2014-%u, Ecole Centrale de Nantes" % __year__

        copyright_annotation = vtk.vtkCornerAnnotation()
        copyright_annotation.SetLinearFontScaleFactor(1)
        copyright_annotation.SetNonlinearFontScaleFactor(1)
        copyright_annotation.SetMaximumFontSize(12)
        copyright_annotation.SetText(1, copyright_text)
        copyright_annotation.GetTextProperty().SetColor(0., 0., 0.)
        self.renderer.AddViewProp(copyright_annotation)
        
        self.normals = []
        self.axes = []
        self.oxy_plane = None

        self.polydatas = list()
        self.transforms = list()
        self.actors = list()
        # self.hiden = dict() # TODO: A terminer -> cf methode self.hide()
    
    def show_normals(self):
        for polydata in self.polydatas:
            normals = vtk.vtkPolyDataNormals()
            normals.ConsistencyOff()
            # normals.ComputePointNormalsOn()
            normals.ComputeCellNormalsOn()
            if vtk.VTK_MAJOR_VERSION <= 5:
                normals.SetInput(polydata)
            else:
                normals.SetInputData(polydata)
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
            arrows.SetTipLength(0.5)
            arrows.SetTipRadius(0.1)

            glyph = vtk.vtkGlyph3D()
            glyph.SetSourceConnection(arrows.GetOutputPort())
            glyph.SetInputConnection(normals_at_centers.GetOutputPort())
            glyph.SetVectorModeToUseNormal()
            glyph.SetScaleFactor(1) # FIXME: may be too big ...
            # glyph.SetVectorModeToUseNormal()
            # glyph.SetVectorModeToUseVector()
            # glyph.SetScaleModeToDataScalingOff()
            glyph.Update()

            glyph_mapper = vtk.vtkPolyDataMapper()
            glyph_mapper.SetInputConnection(glyph.GetOutputPort())

            glyph_actor = vtk.vtkActor()
            glyph_actor.SetMapper(glyph_mapper)

            self.renderer.AddActor(glyph_actor)
            self.normals.append(glyph_actor)
            return
        
    def normals_on(self):
        self.normals = True

    def normals_off(self):
        self.normals = False
    
    def plane_on(self):
        pd = self.polydatas[0] # FIXME: pas robuste
        
        plane = vtk.vtkPlaneSource()
        (xmin, xmax, ymin, ymax, _, _) = pd.GetBounds()
        
        dx = 0.1 * (xmax - xmin)
        dy = 0.1 * (ymax - ymin)

        plane.SetOrigin(xmin - dx, ymax + dy, 0)
        plane.SetPoint1(xmin - dx, ymin - dy, 0)
        plane.SetPoint2(xmax + dx, ymax + dy, 0)
        plane.Update()
        polydata = plane.GetOutput()
        
        mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.SetInput(polydata)
        else:
            mapper.SetInputData(polydata)
            
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        
        color = [0., 102. / 255, 204. / 255]
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetEdgeColor(0, 0, 0)
        actor.GetProperty().SetLineWidth(1)
        
        self.renderer.AddActor(actor)
        self.renderer.Modified()
        self.oxy_plane = actor
        
        
        return
    
    def add_plane(self, center, normal):
        plane = vtk.vtkPlaneSource()
        plane.SetCenter(center)
        plane.SetNormal(normal)

        mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.SetInput(plane.GetOutput())
        else:
            mapper.SetInputData(plane.GetOutput())

    def add_polydata(self, polydata, color=[1, 1, 0]):
        if not isinstance(polydata, vtk.vtkPolyData):
            raise TypeError, 'polydata must be a vtkPolyData object'

        pd_id = len(self.polydatas)
        
        transform = vtk.vtkTransform()
        transform.Identity()

        transform_filter = vtk.vtkTransformPolyDataFilter()
        if vtk.VTK_MAJOR_VERSION <= 5:
            transform_filter.SetInput(polydata)
        else:
            transform_filter.SetInputData(polydata)
            
        transform_filter.SetTransform(transform)
        
        # Building mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(transform_filter.GetOutputPort())
        # if vtk.VTK_MAJOR_VERSION <= 5:
        #     mapper.SetInput(polydata)
        # else:
        #     mapper.SetInputData(polydata)

        # Building actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # Properties setting
        actor.GetProperty().SetColor(color)
        actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetEdgeColor(0, 0, 0)
        actor.GetProperty().SetLineWidth(1)
        actor.GetProperty().SetPointSize(10)

        self.polydatas.append(polydata)
        self.transforms.append(transform)
        self.actors.append(actor)
        
        self.renderer.AddActor(actor)
        self.renderer.Modified()
        
        return pd_id

    def show_axes(self):

        tprop = vtk.vtkTextProperty()
        tprop.SetColor(0., 0., 0.)
        tprop.ShadowOn()

        axes = vtk.vtkCubeAxesActor2D()
        if vtk.VTK_MAJOR_VERSION <= 5:
            axes.SetInput(self.polydatas[0])
        else:
            axes.SetInputData(self.polydatas[0])

        axes.SetCamera(self.renderer.GetActiveCamera())
        axes.SetLabelFormat("%6.4g")
        axes.SetFlyModeToOuterEdges()
        axes.SetFontFactor(0.8)
        axes.SetAxisTitleTextProperty(tprop)
        axes.SetAxisLabelTextProperty(tprop)
        # axes.DrawGridLinesOn()

        self.renderer.AddViewProp(axes)
        self.axes.append(axes)

    def show(self):
        self.renderer.ResetCamera()
        self.render_window.Render()
        self.render_window_interactor.Start()

    def hide(self, id):
        if id > len(self.polydatas):
            print "No mesh with id %u" % id
            return
        
        self.hiden

    def save(self):
        from vtk import vtkXMLPolyDataWriter

        writer = vtkXMLPolyDataWriter()
        writer.SetDataModeToAscii()
        writer.SetFileName('mmviewer_save.vtp')

        for polydata in self.polydatas:
            if vtk.VTK_MAJOR_VERSION <= 5:
                writer.SetInput(polydata)
            else:
                writer.SetInputData(polydata)
        writer.Write()

        print "File 'mmviewer_save.vtp' written in %s" % getcwd()
        return

    def screenshot(self):
        w2if = vtk.vtkWindowToImageFilter()
        w2if.SetInput(self.render_window)
        w2if.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName("screenshot.png")
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(w2if.GetOutput())
        else:
            writer.SetInputData(w2if.GetOutput())
        writer.Write()

        print "File 'screenshot.png' written in %s" % getcwd()
        return

    def finalize(self):
        del self.render_window
        del self.render_window_interactor

    def on_key_press(self, obj, event):
        key = obj.GetKeySym()
        
        if key == 'n':
            if self.normals:
                # self.normals = False
                for actor in self.normals:
                    self.renderer.RemoveActor(actor)
                self.renderer.Render()
                self.normals = []
            else:
                self.show_normals()
                self.renderer.Render()
        elif key == 'b':
            if self.axes:
                for axis in self.axes:
                    self.renderer.RemoveActor(axis)
                self.axes = []
            else:
                self.show_axes()

        elif key == 'e' or key == 'q':
            self.render_window_interactor.GetRenderWindow().Finalize()
            self.render_window_interactor.TerminateApp()

        elif key == 'x':
            self.save()

        elif key == 'c':
            # pass
            self.screenshot()
        
        elif key == 'h':
            if self.oxy_plane:
                self.renderer.RemoveActor(self.oxy_plane)
                self.oxy_plane = None
            else:
                self.plane_on()
        
        
