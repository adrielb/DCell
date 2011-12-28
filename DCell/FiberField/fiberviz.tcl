package require vtk
package require vtkinteraction
package require vtktesting

# Generate some random points
#
vtkMath math
vtkPoints points
set len 500
set x1 0.0
set x2 0.0
set x3 0.0
for {set i 0} {$i<$len} {incr i 1} {
    set x1 [expr $x1+[math Random -1 1]]
    set x2 [expr $x2+[math Random -1 1]]
    set x3 [expr $x3+[math Random -1 1]]
    eval points InsertPoint $i $x1 $x2 $x3
}

vtkPolyLine polyline
    [polyline GetPointIds] SetNumberOfIds $len
    for {set i 0} {$i<$len} {incr i 1} {
        [polyline GetPointIds] SetId $i $i
    }
vtkCellArray cells
  cells InsertNextCell polyline
vtkPolyData profile
    profile SetPoints points
    profile SetLines cells
vtkTubeFilter tubes
    tubes SetInput profile
    tubes SetRadius 1
    tubes SetNumberOfSides 64
    tubes CappingOn
vtkSmoothPolyDataFilter smoother
    smoother SetInputConnection [tubes GetOutputPort]
    smoother SetNumberOfIterations 10
vtkPolyDataMapper mapEdges
    mapEdges SetInputConnection [smoother GetOutputPort]
vtkActor edgeActor
    edgeActor SetMapper mapEdges
eval [edgeActor GetProperty] SetColor $peacock
    [edgeActor GetProperty] SetSpecularColor 1 1 1
    [edgeActor GetProperty] SetSpecular 0.3
    [edgeActor GetProperty] SetSpecularPower 20
    [edgeActor GetProperty] SetAmbient 0.2
    [edgeActor GetProperty] SetDiffuse 0.8


# Create graphics objects
# Create the rendering window, renderer, and interactive renderer
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# Add the actors to the renderer, set the background and size
ren1 AddActor edgeActor
ren1 SetBackground 1 1 1
renWin SetSize 500 500

# render the image
#
iren AddObserver UserEvent {wm deiconify .vtkInteract}
ren1 ResetCamera
[ren1 GetActiveCamera] Zoom 1.5
iren Initialize

# prevent the tk window from showing up then start the event loop
wm withdraw .
