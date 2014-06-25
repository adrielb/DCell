import sys
import vtk
import glob
import vtk.util.colors as colors
import numpy as np
import pandas as pd
import igraph as ig

# edge / vert / fiber import {{{
FILE_OFFSET = 10000
time_index = 0

def GetMaxTime( ):
    files = glob.glob( 'verts.*.Real64' )
    return len(files) - 1

MAX_TIME = GetMaxTime()

def readEdgeList( time_index ):
    files = glob.glob( 'edgelist.'+str(time_index + FILE_OFFSET)+'.*.csv' )
    files.sort()
    lists = [ pd.read_csv( f ) for f in files ]
    return pd.concat( lists, ignore_index=True ).astype(np.int32)

def readVertPositionList( time_index ):
    # > - endian, f - float, 8 - bytes
    dt = np.dtype( [('x','>f8'), ('y','>f8'), ('z','>f8')] )
    filename = 'verts.'+str(time_index + FILE_OFFSET)+'.Real64'
    file = open( filename, 'rb' ) 
    # header is 2 * sizeof(int) for petsc vecs
    file.read(2*4)
    return np.fromfile( file, dt)

def edgeListToFibers( edgelist, verts ):
    etype = edgelist['etype'] == 1
    v0 = edgelist[etype]['vPO0']
    v1 = edgelist[etype]['vPO1']

    g = ig.Graph()
    g.add_vertices(verts.shape[0])

    for i in xrange(v0.shape[0]):
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

class Time_IndexActor(vtk.vtkTextActor): # {{{
    def __init__(self):
        #vtk.vtkTextActor.__init__(self)
        #super(Time_IndexActor, self).__init__()
        #self.textActor = vtk.vtkTextActor()
        #textActor.SetTextScaleModeToProp()

        tprop = self.GetTextProperty()
        tprop.SetFontSize( 19 )
        tprop.BoldOn()
        tprop.SetFontFamilyToArial()
        tprop.SetColor( 0, 1, 0)

    def UpdateTimeIndex( self, time_index ):
        self.SetInput( "time_index = " + str(time_index) )
#}}}

class VertexActor(vtk.vtkActor): # {{{
    def __init__(self):
        self.polyPoints = vtk.vtkPoints()
        
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints( self.polyPoints )

        self.balls = vtk.vtkSphereSource()
        self.balls.SetRadius(0.02)
        self.balls.SetPhiResolution(10)
        self.balls.SetThetaResolution(10)

        self.glyphPoints = vtk.vtkGlyph3D()
        self.glyphPoints.SetInputData( self.polydata )
        self.glyphPoints.SetSourceConnection( self.balls.GetOutputPort() )

        self.glyphMapper = vtk.vtkPolyDataMapper()
        self.glyphMapper.SetInputConnection( self.glyphPoints.GetOutputPort() )

        self.SetMapper(self.glyphMapper)
        self.GetProperty().SetDiffuseColor(colors.dark_orange)
        self.GetProperty().SetSpecular(0.3)
        self.GetProperty().SetSpecularPower(30)


    def UpdateTimeIndex( self, time_index ):
        self.verts = readVertPositionList( time_index )
        self.polyPoints.SetNumberOfPoints( self.verts.size )
        for i, X in enumerate(self.verts):
            self.polyPoints.InsertPoint( i, *X )
        self.polydata.Modified()
# }}}

class FiberActor( vtk.vtkActor ): #{{{
    def __init__(self):
        self.polylinePoints = vtk.vtkPoints()
        self.cells          = vtk.vtkCellArray()
        self.polydata       = vtk.vtkPolyData()
        self.tube           = vtk.vtkTubeFilter()
        self.polyDataMapper = vtk.vtkPolyDataMapper()

        self.polydata.SetPoints( self.polylinePoints )
        self.polydata.SetLines( self.cells )

#smoother = vtk.vtkSmoothPolyDataFilter()
#smoother.SetInputConnection( tube.GetOutputPort () )
#smoother.SetNumberOfIterations( 1 )

        self.tube.SetInputData( self.polydata )
        self.tube.SetRadius( 0.01 )
        self.tube.SetNumberOfSides( 64 )
        self.tube.CappingOn()

        self.polyDataMapper.SetInputConnection( self.tube.GetOutputPort( ))

        self.SetMapper( self.polyDataMapper )


    def UpdateTimeIndex( self, time_index ):
        (edgelist, verts, fibers) = read_fibers( time_index );

        self.polylinePoints.SetNumberOfPoints( verts.size )
        for i,X in enumerate(verts):
            self.polylinePoints.InsertPoint(i, X[0], X[1], X[2] )

        self.cells.Reset()
        for f in fibers:
            self.polyline = vtk.vtkPolyLine()
            self.polyline.GetPointIds().SetNumberOfIds( len(f) )
            for i,v in enumerate(f):
                self.polyline.GetPointIds().SetId( i, v )
            self.cells.InsertNextCell( self.polyline )

        self.polydata.Modified()
#}}}

#polylineExport = vtk.vtkOBJExporter()


#factGraphics = vtk.vtkGraphicsFactory()
#factGraphics.SetUseMesaClasses(1);


#ren = vtk.vtkMesaRenderer()
#renWin = vtk.vtkXMesaRenderWindow()
#renWin.OffScreenRenderingOn()

class MyRenderer(vtk.vtkRenderer): # {{{
    def __init__(self):
        # No cmd line time_index, use on screen rendering
        self.onScreen = len(sys.argv) == 1
        self.SetBackground( 0, 0, 0)

        self.renWin = vtk.vtkRenderWindow()
        self.renWin.AddRenderer( self )

        if self.onScreen:
            self.renWin.FullScreenOn()
            self.iren = vtk.vtkRenderWindowInteractor()
            self.iren.SetRenderWindow( self.renWin )
            self.iren.AddObserver("KeyPressEvent", self.Keypress_ShiftTime)
        else:
            self.renWin.OffScreenRenderingOn()
            self.renWin.SetSize( 1280, 720 )
            
        self.axesActor = vtk.vtkCubeAxesActor();
        self.axesActor.SetFlyModeToStaticTriad()
        #actorAxes->SetCamera(camOrtho);
        #double b[6] = {2,16,2,16,3,8};
        #actorAxes->SetBounds(b);
        self.axesCamera = vtk.vtkCamera()
        self.axesCamera.SetPosition( 7, 7, 7)
        self.axesCamera.SetFocalPoint( 0, 0, 0)
        self.axesActor.SetCamera( self.axesCamera )

        self.AddActor( self.axesActor )

    def SetUpdateData( self, func ):
        self.UpdateData = func

    def ShiftTimeIndex( self, delta ):
        global time_index
        time_index = time_index + delta
        time_index = np.clip( time_index, 0, MAX_TIME )
        self.UpdateData(time_index)
        self.renWin.Render()

    def OrthogonalView( self, axis):
        orthoHeight = 5
        eps = 1e-4
        pos = [eps,eps,eps]
        pos[axis] = orthoHeight
        self.GetActiveCamera().SetPosition(*pos) 
        self.GetActiveCamera().SetRoll(0) 
        self.renWin.Render()

    def Keypress_ShiftTime(self, obj, event):
        A = 5
        key = obj.GetKeySym()
        #print 'key = ' + key
        if key == "j":
            self.ShiftTimeIndex(1)
        elif key == "k":
            self.ShiftTimeIndex(-1)
        elif key == "h":
            self.ShiftTimeIndex(-10)
        elif key == "l":
            self.ShiftTimeIndex(10)
        elif key == "semicolon":
            self.ShiftTimeIndex(-999999)
        elif key == "apostrophe":
            self.ShiftTimeIndex( 999999) 

        elif key == "u":
            self.OrthogonalView( 0 )
        elif key == "i":
            self.OrthogonalView( 1 )
        elif key == "o":
            self.OrthogonalView( 2 )
        elif key == "7":
            self.GetActiveCamera().SetPosition(-A,0,0) 
        elif key == "8":
            self.GetActiveCamera().SetPosition(0,-A,0) 
        elif key == "9":
            self.GetActiveCamera().SetPosition(0,0,-A) 

    def Draw( self ):
        global time_index
        if self.onScreen:
            self.UpdateData( time_index )
            self.iren.Initialize()
            self.iren.Start()
        else:
            #pass
            time_index = int(sys.argv[1])
            self.UpdateData( time_index )
            #self.Render()
            self.renWin.Render()

            windowToImageFilter = vtk.vtkWindowToImageFilter()
            windowToImageFilter.SetInput(self.renWin)
            windowToImageFilter.Update()
            
            writer = vtk.vtkPNGWriter()
            writer.SetFileName("output."+str(FILE_OFFSET+time_index)+".png")
            writer.SetInputConnection(windowToImageFilter.GetOutputPort())
            writer.Write()


# }}}
