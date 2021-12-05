#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import vtk
import math
import time
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from integration import SolveEvolution

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

import warnings


class PlanetSatelliteSystem:
    def __init__(self, semimajor_axis, eccentricity, inclination):
        self.setParams(semimajor_axis, eccentricity, inclination)

    def setParams(self, semimajor_axis, eccentricity, inclination):
        self.semimajor_axis = semimajor_axis
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.semiminor_axis = self.getSemiminorAxis()
        self.focus = self.getFocus()

    def getFocus(self):
        return self.eccentricity * self.semimajor_axis

    def getSemiminorAxis(self):
        return self.semimajor_axis * math.sqrt(1 - math.pow(self.eccentricity, 2))

    def makePlanet(self):
        # Create planet
        planet = vtk.vtkSphereSource()
        planet.SetThetaResolution(30)
        planet.SetPhiResolution(12)
        planet.SetCenter(0, 0, 0)
        planet.SetRadius(2)

        # Create a planet mapper
        planet_mapper = vtk.vtkPolyDataMapper()
        planet_mapper.SetInputConnection(planet.GetOutputPort())

        # Create a planet actor
        planet_actor = vtk.vtkActor()
        planet_actor.SetMapper(planet_mapper)
        return planet_actor

    def makePlanetPlane(self):
        # create source
        resolution = 17
        planet_plane = vtk.vtkPlaneSource()
        planet_plane.SetCenter(15.0, 5.0, 3.0)
        planet_plane.SetNormal(0, 0, 1)
        planet_plane.SetOrigin(-15.0, -15.0, 0.0)
        planet_plane.SetPoint1(15.0, -15.0, 0.0)
        planet_plane.SetPoint2(-15.0, 15.0, 0.0)
        planet_plane.SetXResolution(resolution)
        planet_plane.SetYResolution(resolution)
        transform = vtk.vtkTransform()
        transform.Translate(-5.0, -5.0, 0.0)
        # planet_plane.SetUserTransform(transform)
        planet_plane.Update()

        # mapper
        planet_plane_mapper = vtk.vtkPolyDataMapper()
        planet_plane_mapper.SetInput(planet_plane.GetOutput())

        # actor
        planet_plane_actor = vtk.vtkActor()
        planet_plane_actor.SetMapper(planet_plane_mapper)
        planet_plane_actor.GetProperty().SetOpacity(.40)
        # Color
        planet_plane_actor.GetProperty().SetColor(0.5, 0, 0)

        return planet_plane_actor

    def makeSatellite(self):
        # Create satellite
        satellite = vtk.vtkSphereSource()
        satellite.SetThetaResolution(30)
        satellite.SetPhiResolution(12)
        satellite.SetCenter(0, 0, 0)
        satellite.SetRadius(1)

        # transform = vtk.vtkTransform()
        # # transform->RotateWXYZ(double angle, double x, double y, double z);
        # transform.RotateWXYZ(45, 0, 1, 1)
        #
        # transformFilter = vtk.vtkTransformPolyDataFilter()
        #
        # transformFilter.SetTransform(transform)
        # transformFilter.SetInputConnection(satellite.GetOutputPort())
        # transformFilter.Update()


        # Create a satellite mapper
        satellite_mapper = vtk.vtkPolyDataMapper()
        satellite_mapper.SetInputConnection(satellite.GetOutputPort())
        # satellite_mapper.SetInputConnection(transformFilter.GetOutputPort())
        # Create a satellite actor
        satellite_actor = vtk.vtkActor()
        satellite_actor.SetMapper(satellite_mapper)
        # satellite_actor.SetPosition(9, 0, 0)
        # satellite_actor.RotateWXYZ(180, 0, 0, 1)
        # satellite_actor.RotateZ(90)
        # satellite_actor.RotateY(-180)
        # satellite_actor.RotateX(-45)
        return satellite_actor

    def makeOrbit(self):
        # Create orbit ellipse
        a = self.semimajor_axis
        b = self.semiminor_axis

        # vtkPoints represents 3D points. The data model for vtkPoints is an array of
        # vx-vy-vz triplets accessible by (point or cell) id.
        # Добавляем точки для нашего эллипса
        cont = 0
        i = 0
        orbitPoints = vtk.vtkPoints()
        while i <= (2 * math.pi):
            orbitPoints.InsertPoint(cont, a * math.cos(i), b * math.sin(i), 0)  # math.pi*math.cos(i)
            i += (math.pi / 20)
            cont += 1

        orbitPoints.InsertPoint(++cont, self.focus, 0, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, self.focus, 1.0, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, self.focus, -1.0, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, self.focus, 0, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, 0, 0, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, 0.5, 0.5, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, -0.5, -0.5, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, 0, 0, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, -0.5, 0.5, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, 0.5, -0.5, 0)
        cont += 1
        orbitPoints.InsertPoint(cont, 0, 0, 0)
        cont += 1
        # orbitPoints.InsertPoint(cont, 0, 0, 0)
        orbitPoints.InsertPoint(cont, a * math.cos(math.pi), 0, 0)
        cont += 1

        # vtkCellArray is a supporting object that explicitly represents cell connectivity.
        # The cell array structure is a raw integer list of the form:
        # (n,id1,id2,...,idn, n,id1,id2,...,idn, ...) where n is the number of points in
        # the cell, and id is a zero-offset index into an associated point list.
        # Создаем список точек
        orbit_lines = vtk.vtkCellArray()
        orbit_lines.InsertNextCell(cont)
        k = 0
        while k < cont:
            orbit_lines.InsertCellPoint(k)
            k += 1

        # vtkPolyData is a data object that is a concrete implementation of vtkDataSet.
        # vtkPolyData represents a geometric structure consisting of vertices, lines,
        # polygons, and/or triangle strips
        orbit_polyData = vtk.vtkPolyData()
        orbit_polyData.SetPoints(orbitPoints)
        orbit_polyData.SetLines(orbit_lines)
        # orbit_polyData.SetCenter(9,0,0)

        transform = vtk.vtkTransform()
        # transform->RotateWXYZ(double angle, double x, double y, double z);
        transform.RotateWXYZ(-self.inclination * 180 / math.pi, 0, 1, 0)
        transform.Translate(self.focus, 0, 0.0)
        transformFilter = vtk.vtkTransformPolyDataFilter()

        transformFilter.SetTransform(transform)
        if vtk.VTK_MAJOR_VERSION <= 5:
            transformFilter.SetInputConnection(orbit_polyData.GetProducerPort())  # orbit_polyData.GetProducerPort()
        else:
            transformFilter.SetInputConnection(orbit_polyData)
        transformFilter.Update()

        # Create a mapper
        # vtkPolyDataMapper is a class that maps polygonal data (i.e., vtkPolyData)
        # to graphics primitives
        orbit_polyDataMapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            orbit_polyDataMapper.SetInputConnection(transformFilter.GetOutputPort())  # orbit_polyData.GetProducerPort()
        else:
            orbit_polyDataMapper.SetInputData(orbit_polyData)
            orbit_polyDataMapper.Update()

        # Create an actor
        # Create an actor to represent the polygon. The actor orchestrates rendering of
        # the mapper's graphics primitives. An actor also refers to properties via a
        # vtkProperty instance, and includes an internal transformation matrix. We
        # set this actor's mapper to be polygonMapper which we created above.
        orbit_actor = vtk.vtkActor()
        orbit_actor.SetMapper(orbit_polyDataMapper)
        orbit_actor.SetPosition(0, 0, 0)
        # Opacity
        # orbit_actor.GetProperty().SetOpacity(1)
        # Width
        orbit_actor.GetProperty().SetLineWidth(2.0)
        # Color
        orbit_actor.GetProperty().SetColor(0, 0.9, 1)
        return orbit_actor

    @staticmethod
    def position(a, ec, E):
        # a=semimajor axis, ec=eccentricity, E=eccentric anomaly
        # x,y = coordinates of the satellite with respect to the planet
        coords = {}
        coords['x'] = a * (math.cos(E) - ec)
        coords['y'] = a * math.sqrt(1.0 - ec * ec) * math.sin(E)
        return coords

    @staticmethod
    def EccAnom(ec, m, dp):
        # arguments:
        # ec=eccentricity, m=mean anomaly,
        # dp=number of decimal places

        maxIter = 30
        i = 0
        delta = math.pow(10, -dp)
        k = 0.85
        E = m + math.copysign(1, math.sin(m)) * k * ec
        F = E - ec * math.sin(E) - m

        while ((abs(F) > delta) and (i < maxIter)):
            E = E - F / (1.0 - ec * math.cos(E))
            F = E - ec * math.sin(E) - m
            i = i + 1

        return round(E * math.pow(10, dp)) / math.pow(10, dp)


class vtkTimerCallback():
    def __init__(self):
        self.timer_count = 0
        self.time = 0
        self.start_time = time.time()

    def execute(self, obj, event):
        # print(self.timer_count)
        new_time = time.time() - self.start_time

        m = new_time * self.pss.n * 86400 * math.pow(10, self.pss.speed)
        a = self.pss.semimajor_axis
        ecc = self.pss.eccentricity

        E = self.pss.EccAnom(ecc, m, 4)
        coords = self.pss.position(a, ecc, E)

        angle = self.pss.inclination  # * math.pi / 180.0
        x1 = coords['x']  # + self.pss.semimajor_axis + self.pss.focus
        x = x1 * math.cos(angle)
        y = coords['y']
        z = x1 * math.sin(angle)
        # print(angle)

        # self.actor.SetPosition(coords['x'], coords['y'], 0)
        self.actor.SetPosition(-x, -y, -z)

        iren = obj
        iren.GetRenderWindow().Render()
        self.timer_count += 1


def Polygon():
    # Setup four points
    points = vtk.vtkPoints()
    points.InsertNextPoint(0.0, 0.0, 0.0)
    points.InsertNextPoint(1.0, 0.0, 0.0)
    points.InsertNextPoint(1.0, 1.0, 0.0)
    points.InsertNextPoint(0.0, 1.0, 0.0)

    # Create the polygon
    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(4)  # make a quad
    polygon.GetPointIds().SetId(0, 0)
    polygon.GetPointIds().SetId(1, 1)
    polygon.GetPointIds().SetId(2, 2)
    polygon.GetPointIds().SetId(3, 3)

    # Add the polygon to a list of polygons
    polygons = vtk.vtkCellArray()
    polygons.InsertNextCell(polygon)

    # Create a PolyData
    polygonPolyData = vtk.vtkPolyData()
    polygonPolyData.SetPoints(points)
    polygonPolyData.SetPolys(polygons)

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        mapper.SetInput(polygonPolyData)
    else:
        mapper.SetInputData(polygonPolyData)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    return actor


class Widget3D(object):
    xRotationChanged = QtCore.pyqtSignal(int)
    yRotationChanged = QtCore.pyqtSignal(int)
    zRotationChanged = QtCore.pyqtSignal(int)
    timePositionChanged = QtCore.pyqtSignal(int)

    def connect(self, MainWindow, ui):

        self.ui = ui
        self.setSliders(ui)
        self.setParams(ui)
        self.setTimeSlider()

        frame = ui.frame3D
        self.vl = QtGui.QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(frame)
        self.vl.addWidget(self.vtkWidget)
        frame.setLayout(self.vl)

        # Create the Renderer. A renderer is like a
        # viewport. It is part or all of a window on the screen and it is
        # responsible for drawing the actors it has.  We also set the
        # background color here.
        self.ren = vtk.vtkRenderer()
        ##self.ren.SetBackground(0.1, 0.2, 0.4)
        self.ren.SetBackground(1, 1,1)

        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)

        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()

        interactorStyle = vtk.vtkInteractorStyleImage()
        self.iren.SetInteractorStyle(interactorStyle)

        self.pss = PlanetSatelliteSystem(10, 0.1, 0)
        self.orbit = self.pss.makeOrbit()
        self.ren.AddActor(self.orbit)
        self.setGraphic()
        self.goIntegrate()

        self.setButtons()
        self.setMenu()

        # Assign actors to the Renderer.
        self.ren.AddActor(self.pss.makePlanet())

        sattelite = self.pss.makeSatellite()
        self.ren.AddActor(sattelite)
        # self.ren.AddActor(Polygon())
        self.ren.AddActor(self.pss.makePlanetPlane())

        # self.ren.Render()
        #
        # self.ren.AddActor(self.orbit)
        # self.pss.setParams(10, self.ecc, self.incl)

        # self.ren.RemoveActor(self.orbit)
        # self.orbit = self.pss.makeOrbit()


        transform = vtk.vtkTransform()
        transform.Translate(6.0, 6.0, 0.0)

        axes = vtk.vtkAxesActor()
        #  The axes are positioned with a user transform
        axes.SetUserTransform(transform)

        # properties of the axes labels can be set as follows
        # this sets the x axis label to red
        # axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1,0,0);
        axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1, 0, 0)
        # the actual text of the axis label can be changed:
        # axes->SetXAxisLabelText("test");
        # self.ren.AddActor(axes)

        # Automatically set up the camera based on the visible actors.
        # The camera will reposition itself to view the center point of the actors,
        # and move along its initial view plane normal
        # (i.e., vector defined from camera position to focal point) so that all of the
        # actors can be seen.
        self.ren.ResetCamera()

        self.camera = vtk.vtkCamera()
        # self.camera.SetPosition(80, 0, 0)
        # self.camera.SetFocalPoint(0, 1, 0)
        # self.camera.Azimuth(90)
        # self.camera.SetViewUp(0,0,0)

        # Create a renderer, render window, and interactor
        self.ren.SetActiveCamera(self.camera)
        self.camera.Elevation(-60)
        self.ren.ResetCamera()

        # Sign up to receive TimerEvent
        self.cb = vtkTimerCallback()
        self.cb.actor = sattelite
        self.cb.pss = self.pss
        self.iren.AddObserver('TimerEvent', self.cb.execute)
        # self.timerId = self.iren.CreateRepeatingTimer(100000)
        # self.iren.InvokeEvent('TimerEvent', self.cb.execute)

        # self.frame.setLayout(self.vl)
        # ww.addWidget(self.frame)
        # self.widget3D = QtGui.QWidget(self.centralwidget)

        # MainWindow.show()
        # self.show()
        # self.iren.Initialize()
        # self.iren.Start()

        MainWindow.ren = self.ren
        MainWindow.iren = self.iren
        self.MainWindow = MainWindow
        # MainWindow.w2D = self.w2D

    def setSliders(self, ui):
        self.xRot = 180
        self.yRot = 30
        ui.ySlider.setValue(30)
        ui.ySlider.setRange(1, 179)
        ui.xSlider.valueChanged.connect(self.setXRotation)
        ui.ySlider.valueChanged.connect(self.setYRotation)

    def setXRotation(self, angle):
        if angle != self.xRot:
            self.camera.Azimuth(angle - self.xRot)
            self.xRot = angle

    def setYRotation(self, angle):
        if angle != self.yRot:
            self.camera.Elevation(angle - self.yRot)
            self.yRot = angle

    def setParams(self, ui):
        self.param_a = 384.4
        self.param_i = 5.15
        self.param_e = 0.0549
        self.param_mp = 5.9736
        self.param_r0 = 6.378
        self.param_Tp = 23.93419
        self.param_ms = 734.9

        self.param_speed = 0
        self.param_time_integration = 10

        ui.edit_a.setValue(self.param_a)
        ui.edit_e.setValue(self.param_e)
        ui.edit_i.setValue(self.param_i)
        ui.edit_mp.setValue(self.param_mp)
        ui.edit_r0.setValue(self.param_r0)
        ui.edit_Tp.setValue(self.param_Tp)
        ui.edit_ms.setValue(self.param_ms)

        ui.edit_speed.setValue(self.param_speed)
        ui.edit_time_integration.setValue(self.param_time_integration)

        ui.edit_a.valueChanged.connect(self.setA)
        ui.edit_e.valueChanged.connect(self.setE)
        ui.edit_i.valueChanged.connect(self.setI)
        ui.edit_mp.valueChanged.connect(self.setMp)
        ui.edit_r0.valueChanged.connect(self.setR0)
        ui.edit_Tp.valueChanged.connect(self.setTp)
        ui.edit_ms.valueChanged.connect(self.setMs)

        ui.edit_speed.valueChanged.connect(self.setSpeed)
        ui.edit_time_integration.valueChanged.connect(self.setTimeIntegration)

    def setA(self, value):
        if value != self.param_a:
            self.param_a = value

    def setE(self, value):
        if value != self.param_e:
            self.param_e = value

    def setI(self, value):
        if value != self.param_i:
            self.param_i = value

    def setMp(self, value):
        if value != self.param_mp:
            self.param_mp = value

    def setR0(self, value):
        if value != self.param_r0:
            self.param_r0 = value

    def setTp(self, value):
        if value != self.param_Tp:
            self.param_Tp = value

    def setMs(self, value):
        if value != self.param_ms:
            self.param_ms = value

    def setSpeed(self, value):
        if value != self.param_speed:
            self.param_speed = value
            self.pss.speed = self.param_speed

    def setTimeIntegration(self, value):
        if value != self.param_time_integration:
            self.param_time_integration = value

    def goIntegrate(self):
        a = self.param_a * (10 ** 6)
        mp = self.param_mp * (10 ** 24)
        r0 = self.param_r0 * (10 ** 6)
        ms = self.param_ms * (10 ** 20)
        self.solve_eq = SolveEvolution(ecc=self.param_e, semi_a=a, incl=self.param_i, mp=mp, r0=r0, Tp=self.param_Tp,
                                       ms=ms, time=self.param_time_integration)
        with warnings.catch_warnings():
            try:
                self.integr_result = self.solve_eq.integrate()
            except:
                self.ui.msgTextBox.setText("Ошибка интегрирования, попробуйте изменить начальные условия или время")
        #
        self.ui.timeSlider.setValue(0)
        self.setEvolutionParams(0)
        self.plot(0)

    def setTimeSlider(self):
        self.ui.timeSlider.valueChanged.connect(self.setEvolutionParams)

    def setEvolutionParams(self, value):
        self.timeValue = value
        self.ui.evol_e.setText(str(self.integr_result[2][value]))
        self.ui.evol_i.setText(str(self.integr_result[3][value] * 180.0 / math.pi))
        self.ui.evol_a.setText(str(self.integr_result[4][value]))
        self.plot(value)
        self.ecc = self.integr_result[2][value]
        self.incl = self.integr_result[3][value]
        self.pss.n = self.integr_result[5][value]
        self.pss.speed = self.param_speed
        self.pss.setParams(10, self.ecc, self.incl)
        self.ren.RemoveActor(self.orbit)
        self.orbit = self.pss.makeOrbit()
        self.ren.AddActor(self.orbit)

    def setGraphic(self):
        frame_e = self.ui.tab_e
        # a figure instance to plot on
        self.figure_e = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas_e = FigureCanvas(self.figure_e)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        # self.toolbar = NavigationToolbar(self.canvas, frame)

        # set the layout
        layout_e = QtGui.QVBoxLayout()
        # layout.addWidget(self.toolbar)
        layout_e.addWidget(self.canvas_e)
        frame_e.setLayout(layout_e)

        frame_i = self.ui.tab_i
        self.figure_i = plt.figure()
        self.canvas_i = FigureCanvas(self.figure_i)
        layout_i = QtGui.QVBoxLayout()
        layout_i.addWidget(self.canvas_i)
        frame_i.setLayout(layout_i)

        frame_a = self.ui.tab_a
        self.figure_a = plt.figure()
        self.canvas_a = FigureCanvas(self.figure_a)
        layout_a = QtGui.QVBoxLayout()
        layout_a.addWidget(self.canvas_a)
        frame_a.setLayout(layout_a)
        #
        # frame_en0 = self.ui.tab_en0
        # self.figure_en0 = plt.figure()
        # self.canvas_en0 = FigureCanvas(self.figure_en0)
        # layout_en0 = QtGui.QVBoxLayout()
        # layout_en0.addWidget(self.canvas_en0)
        # frame_en0.setLayout(layout_en0)
        #
        # frame_in0 = self.ui.tab_in0
        # self.figure_in0 = plt.figure()
        # self.canvas_in0 = FigureCanvas(self.figure_in0)
        # layout_in0 = QtGui.QVBoxLayout()
        # layout_in0.addWidget(self.canvas_in0)
        # frame_in0.setLayout(layout_in0)

        self.canvas_list = [
            [self.canvas_e, self.figure_e, 0, 2],
            [self.canvas_i, self.figure_i, 0, 3],
            [self.canvas_a, self.figure_a, 0, 4],
            # [self.canvas_en0, self.figure_en0, 2, 1],
            # [self.canvas_in0, self.figure_in0, 3, 1]
        ]

        # self.ui.graphicTabWidget.setCurrentIndex(0)
        self.currentTab = 0
        # self.ui.graphicTabWidget.currentChanged.connect(self.onGraphicTabChange)

    def setButtons(self):
        self.ui.solveButton.clicked.connect(self.goIntegrate)
        self.ui.startButton.clicked.connect(self.startAnimation)
        self.ui.stopButton.clicked.connect(self.stopAnimation)

    def setMenu(self):
        self.ui.action_save.triggered.connect(self.saveParamsToFile)
        self.ui.action_load.triggered.connect(self.loadParamsFromFile)

    def saveParamsToFile(self):
        """
        Function that will save params to file.
        """
        # save to TESTING.txt in current working directory
        savePath = os.path.join(os.getcwd(),
                                'params', 'params.json')

        # allow user to override location if requested
        savePath = QtGui.QFileDialog.getSaveFileName(self.MainWindow,
                                                     'Save params to file',
                                                     savePath,
                                                     '*.json')

        # if just saving, or user didn't cancel, make and save file
        if len(savePath) > 0:
            with open(savePath, 'w') as theFile:
                data = self.getParams()
                json.dump(data, theFile, sort_keys=False, indent=4, ensure_ascii=False)

    def getParams(self):
        params = {
            'a': self.param_a,
            'i': self.param_i,
            'e': self.param_e,
            'mp': self.param_mp,
            'r0': self.param_r0,
            'Tp': self.param_Tp,
            'ms': self.param_ms,
            'speed': self.param_speed,
            'time_integration': self.param_time_integration
        }

        return params

    def setLoadedParams(self, params):

        if params['a']:
            self.ui.edit_a.setValue(params['a'])

        if params['e']:
            self.ui.edit_e.setValue(params['e'])

        if params['i']:
            self.ui.edit_i.setValue(params['i'])

        if params['mp']:
            self.ui.edit_mp.setValue(params['mp'])

        if params['r0']:
            self.ui.edit_r0.setValue(params['r0'])

        if params['Tp']:
            self.ui.edit_Tp.setValue(params['Tp'])

        if params['ms']:
            self.ui.edit_ms.setValue(params['ms'])

        if params['speed']:
            self.ui.edit_speed.setValue(params['speed'])

        if params['time_integration']:
            self.ui.edit_time_integration.setValue(params['time_integration'])


    def loadParamsFromFile(self):
        """
        Function that will load params from file.
        """
        # save to TESTING.txt in current working directory
        loadPath = os.path.join(os.getcwd(),
                                'params')
        # allow user to override location if requested
        loadPath = QtGui.QFileDialog.getOpenFileName(self.MainWindow,
                                                     'Load params from file',
                                                     loadPath,
                                                     '*.json')

        # if just saving, or user didn't cancel, make and save file
        if len(loadPath) > 0:
            with open(loadPath, 'r') as theFile:
                data = json.load(theFile)
                self.setLoadedParams(data)

    def startAnimation(self):
        self.cb.start_time = time.time()
        self.timerId = self.iren.CreateRepeatingTimer(1000)

    def stopAnimation(self):
        if self.timerId:
            self.iren.DestroyTimer(self.timerId)

    def plot(self, value):
        ''' plot canvas in current tab '''
        # self.drawCanvas(self.canvas_list[self.currentTab][0],
        #                 self.canvas_list[self.currentTab][1],
        #                 self.integr_result[self.canvas_list[self.currentTab][2]],
        #                 self.integr_result[self.canvas_list[self.currentTab][3]], value)
        self.drawCanvas(self.canvas_e, self.figure_e, self.integr_result[0], self.integr_result[2], value)
        self.drawCanvas(self.canvas_i, self.figure_i, self.integr_result[0], self.integr_result[3], value)
        self.drawCanvas(self.canvas_a, self.figure_a, self.integr_result[0], self.integr_result[4], value)
        # self.drawCanvas(self.canvas_en0, self.figure_en0, self.integr_result[2], self.integr_result[1], value)
        # self.drawCanvas(self.canvas_in0, self.figure_in0, self.integr_result[3], self.integr_result[1], value)

    def drawCanvas(self, canvas, figure, t, x, value):
        # create an axis
        ax = figure.add_subplot(111)
        # discards the old graph
        ax.hold(False)
        # plot data
        ax.plot(t, x, '-')
        ax.hold(True)
        ax.plot(t[value], x[value], '*')
        canvas.draw()

    def onGraphicTabChange(self):
        self.currentTab = self.ui.graphicTabWidget.currentIndex()
        self.plot(self.timeValue)
