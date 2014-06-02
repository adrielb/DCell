import vtk
import vtk.util.colors as colors
import numpy as np
import pandas as pd
import igraph as ig
from numpy.random import randn

# edge / vert / fiber import {{{
FILE_OFFSET = 10000
MAX_TIME = 0

def readEdgeList( time_index ):
    return pd.read_csv( 'edgelist.'+str(time_index + FILE_OFFSET)+'.0.csv' )

def readVertPositionList( time_index ):
    # > - endian, f - float, 8 - bytes
    dt = np.dtype( [('x','>f8'), ('y','>f8'), ('z','>f8')] )
    filename = 'verts.'+str(time_index + FILE_OFFSET)+'.Real64'
    file = open( filename, 'rb' ) 
    # header is 2 * sizeof(int) for petsc vecs
    file.read(2*4)
    return np.fromfile( file, dt)

def edgeListToFibers( edgelist, verts ):
    v0 = edgelist['vPO0']
    v1 = edgelist['vPO1']

    g = ig.Graph()
    g.add_vertices(verts.shape[0])

    for i in xrange(edgelist.shape[0]):
        g.add_edge( v0.iat[i], v1.iat[i] )

    # end points of fibers have vertex degree = 1
    vs = g.vs.select( _degree = 1 )
    endpoints = vs.indices
    
    # for each end point
    #   extract fiber
    #   remove its other endpoint from list
    fibers = []
    while len(endpoints) != 0:
        e = endpoints.pop()
        fiber = g.subcomponent(e)
        otherend = fiber[-1]
        endpoints.remove( otherend )
        fibers.append( fiber )

    return fibers

def read_fibers( time_index ):
    edgelist = readEdgeList( time_index )
    verts = readVertPositionList( time_index )
    fibers = edgeListToFibers( edgelist, verts )
    return (edgelist, verts, fibers)

''' Alternative methods for extracting fibers from an edgelist
clusters = g.clusters()
subgraphs =clusters.subgraphs()
g.get_shortest_paths( 4, to=2 )
g.degree()

# append fiber ID column
# initially ID = -1 meaning not associated with any fiber
edgelist['fID'] = pd.Series( -np.ones(numEdges),index = edgelist.index)
'''
# }}}

MAX_TIME = 99
time_index = 0
(edgelist, verts, fibers) = read_fibers( time_index );

# text actor {{{
textActor = vtk.vtkTextActor()
textActor.SetTextScaleModeToProp()
textActor.SetInput( "Dynamics Test: " + str(time_index) )

tprop = textActor.GetTextProperty()
tprop.SetFontSize( 19 )
tprop.BoldOn()
tprop.SetFontFamilyToArial()
tprop.SetColor( 0, 1, 0)
#}}}

# polydata Actor {{{
edgelist = pd.read_csv( 'edgelist.10000.0.csv' )
numEdges = edgelist.shape[0]
dt = np.dtype( [('x','>f8'), ('y','>f8'), ('z','>f8')] )
verts = np.fromfile( 'verts.10000.Real64', dt)

v0 = edgelist['vPO0']
v1 = edgelist['vPO1']

linePoints = vtk.vtkPoints()
linePoints.SetNumberOfPoints( verts.size )
for i in range(verts.size):
    linePoints.InsertPoint(i, *verts[i] )

cells = vtk.vtkCellArray()
for i in xrange(numEdges): 
    aLine = vtk.vtkLine()
    aLine.GetPointIds().SetId(0, v0.iloc[i])
    aLine.GetPointIds().SetId(1, v1.iloc[i])
    cells.InsertNextCell( aLine )

polydata = vtk.vtkPolyData()
polydata.SetPoints( linePoints )
polydata.SetLines(cells)

tube = vtk.vtkTubeFilter()
tube.SetInputData( polydata )
tube.SetRadius( 0.1 )
tube.SetNumberOfSides( 64 )
tube.CappingOn()

#smoother = vtk.vtkSmoothPolyDataFilter()
#smoother.SetInputConnection( tube.GetOutputPort () )
#smoother.SetNumberOfIterations( 1 )

polyDataMapper = vtk.vtkPolyDataMapper()
polyDataMapper.SetInputConnection( tube.GetOutputPort( ))

polydataActor = vtk.vtkActor()
polydataActor.SetMapper( polyDataMapper )
# }}}

# polyline actor {{{
polylinePoints = vtk.vtkPoints()
cells = vtk.vtkCellArray()
polydata = vtk.vtkPolyData()

def UpdateData( time_index ):
    print time_index

    textActor.SetInput( "Dynamics Test:" + str(time_index) )
    
    (edgelist, verts, fibers) = read_fibers( time_index );

    polylinePoints.SetNumberOfPoints( verts.size )
    for i,X in enumerate(verts):
        polylinePoints.InsertPoint(i, X[0], X[1], X[2] )

    cells.SetNumberOfCells(0)
    for f in fibers:
        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds( len(f) )
        for i,v in enumerate(f):
            polyline.GetPointIds().SetId( i, v )
        cells.InsertNextCell( polyline )

    polydata.Modified()

polydata.SetPoints( polylinePoints )
polydata.SetLines(cells)

tube = vtk.vtkTubeFilter()
tube.SetInputData( polydata )
tube.SetRadius( 0.01 )
tube.SetNumberOfSides( 64 )
tube.CappingOn()

#smoother = vtk.vtkSmoothPolyDataFilter()
#smoother.SetInputConnection( tube.GetOutputPort () )
#smoother.SetNumberOfIterations( 1 )

polyDataMapper = vtk.vtkPolyDataMapper()
polyDataMapper.SetInputConnection( tube.GetOutputPort( ))

polylineActor = vtk.vtkActor()
polylineActor.SetMapper( polyDataMapper )

polylineExport = vtk.vtkOBJExporter()
# }}}

# ball actors {{{
balls = vtk.vtkSphereSource()
balls.SetRadius(0.2)
balls.SetPhiResolution(10)
balls.SetThetaResolution(10)

glyphPoints = vtk.vtkGlyph3D()
glyphPoints.SetInputData( polydata )
glyphPoints.SetSourceConnection( balls.GetOutputPort() )

glyphMapper = vtk.vtkPolyDataMapper()
glyphMapper.SetInputConnection( glyphPoints.GetOutputPort() )

glyphActor = vtk.vtkActor()
glyphActor.SetMapper(glyphMapper)
glyphActor.GetProperty().SetDiffuseColor(colors.dark_orange)
glyphActor.GetProperty().SetSpecular(0.3)
glyphActor.GetProperty().SetSpecularPower(30)
#}}}

# line actor {{{
edgelist = pd.read_csv( 'edgelist.10000.0.csv' )
numEdges = edgelist.shape[0]
dt = np.dtype( [('x','>f8'), ('y','>f8'), ('z','>f8')] )
verts = np.fromfile( 'verts.10000.Real64', dt)

v0 = edgelist['vPO0']
v1 = edgelist['vPO1']

linePoints = vtk.vtkPoints()
linePoints.SetNumberOfPoints( verts.size )
for i in range(verts.size):
    linePoints.InsertPoint(i, *verts[i] )

aLineGrid = vtk.vtkUnstructuredGrid()
aLineGrid.Allocate(numEdges, 1)
aLineGrid.SetPoints(linePoints)

for i in xrange(numEdges): 
    aLine = vtk.vtkLine()
    aLine.GetPointIds().SetId(0, v0.iloc[i])
    aLine.GetPointIds().SetId(1, v1.iloc[i])
    aLineGrid.InsertNextCell(aLine.GetCellType(), aLine.GetPointIds())

aLineMapper = vtk.vtkDataSetMapper()
aLineMapper.SetInputData(aLineGrid)
aLineActor = vtk.vtkActor()
aLineActor.SetMapper(aLineMapper)
aLineActor.GetProperty().SetDiffuseColor(.2, 1, 1)
# }}}


# Renderer, RenderWindow, RenderWindowInteractor {{{
ren = vtk.vtkRenderer()
ren.SetBackground( 0, 0, 0)

camera = ren.GetActiveCamera()
#camera.SetViewUp( 0, 0, 1)
camera.SetViewAngle( 85 )
camera.SetPosition( 0.5, 0.5, 1)
camera.SetFocalPoint( 0.5, 0.5, 0)
camera.ParallelProjectionOn()

renWin = vtk.vtkRenderWindow()
renWin.AddRenderer( ren )

iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow( renWin )
# }}}

# Add actors
ren.AddActor2D( textActor )
#ren.AddActor( aLineActor )
#ren.AddActor(polydataActor)
ren.AddActor( polylineActor )
#ren.AddActor( glyphActor)

def shifttime( delta ):
    global time_index
    time_index = time_index + delta
    time_index = np.clip( time_index, 0, MAX_TIME )
    UpdateData(time_index)
    renWin.Render()

def Keypress(obj, event):
    key = obj.GetKeySym()
    if key == "j":
        shifttime(1)
    elif key == "k":
        shifttime(-1)


iren.AddObserver("KeyPressEvent", Keypress)
iren.Initialize()
iren.Start()

# modifying source points and updating renWin: Examples/Annotation/Python/labeledMesh.py

def save():
    iren.Disable()

    win2image = vtk.vtkWindowToImageFilter()
    win2image.SetInput( renWin )
    writer = vtk.vtkPNGWriter()
    writer.SetInputConnection( win2image.GetOutputPort() )
    writer.SetFileName( "" )
    writer.Write()

