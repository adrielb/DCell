from DCellVTK import *

def UpdateData( time_index ):
    #print time_index
    tiActor.UpdateTimeIndex( time_index )
    vertActor.UpdateTimeIndex( time_index )
    fiberActor.UpdateTimeIndex( time_index )

vertActor = VertexActor()
vertActor.balls.SetRadius(0.02)
        
fiberActor = FiberActor()
fiberActor.tube.SetRadius( 0.01 )

tiActor = Time_IndexActor()
tiActor.GetTextProperty().SetFontSize( 18 )

UpdateData( time_index )
    

ren = MyRenderer(False)
ren.SetUpdateData( UpdateData )
ren.AddActor( vertActor )
ren.AddActor( fiberActor )
ren.AddActor( tiActor )


#ren.ResetCamera()
camera = ren.GetActiveCamera()
#camera.SetViewUp( 0, 0, 1)
#camera.SetViewAngle( 2)
camera.SetPosition( -6, 5, -4)
camera.SetFocalPoint( 0, 0, 0)
#camera.ParallelProjectionOn()

ren.Draw()


