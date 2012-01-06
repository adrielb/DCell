#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkVolume16Reader.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkOutlineFilter.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkPolyDataNormals.h"
#include "vtkContourFilter.h"
#include "vtkRIBExporter.h"
#include "vtkImageReader2.h"
#include "vtkImageDataGeometryFilter.h"
#include "vtkWarpScalar.h"
#include "vtkDataSetMapper.h"
#include "vtkCubeAxesActor2D.h"

int main (int argc, char **argv)
{
  // Create the renderer, the render window, and the interactor. The renderer
  // draws into the render window, the interactor enables mouse- and 
  // keyboard-based interaction with the data within the render window.
  //
  vtkRenderer *aRenderer = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(aRenderer);
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

  
  vtkImageReader2 *reader = vtkImageReader2::New(); 
    reader->SetFileName("/home/abergman/Research/DCell/LevelSetMethod/result.0.Real64");
    reader->SetDataScalarTypeToDouble();
    reader->SetDataExtent( 0, 64, 0, 64, 0, 0);
  
  vtkImageDataGeometryFilter *geom = vtkImageDataGeometryFilter::New();
    geom->SetInputConnection( reader->GetOutputPort() );
  
  vtkWarpScalar *warp = vtkWarpScalar::New();
    warp->SetInputConnection( geom->GetOutputPort() );
    warp->SetNormal( 0, 0, 1);
    warp->UseNormalOn();
  
  vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
     normals->SetInputConnection( warp->GetOutputPort() );
     
  vtkDataSetMapper *mapper = vtkDataSetMapper::New(); 
    mapper->SetInputConnection( normals->GetOutputPort() );
    
  vtkActor *actor = vtkActor::New();
    actor->SetMapper( mapper );
  // It is convenient to create an initial view of the data. The FocalPoint
  // and Position form a vector direction. Later on64 (ResetCamera() method)
  // this vector is used to position the camera to look at the data in
  // this direction.
  vtkCamera *aCamera = vtkCamera::New();
    aCamera->SetViewUp (0, 0, 1);
    aCamera->SetPosition (64, 64, 128);
    aCamera->SetFocalPoint (0, 0, 0);
    aCamera->ComputeViewPlaneNormal();
    aCamera->SetViewAngle( 90 );

  vtkCubeAxesActor2D *axes = vtkCubeAxesActor2D::New();
    axes->SetInput( normals->GetOutput() );
    axes->SetCamera( aCamera );
    axes->SetLabelFormat( "%6.4g" );
    axes->SetFlyModeToOuterEdges();
    axes->SetFontFactor( 5.8 );
  aRenderer->AddViewProp( axes ); 
  
  // Actors are added to the renderer. An initial camera view is created.
  // The Dolly() method moves the camera towards the FocalPoint,
  // thereby enlarging the image.
  aRenderer->AddActor( actor );
  aRenderer->SetActiveCamera(aCamera);
  aRenderer->ResetCamera ();
  aCamera->Dolly(1.5);

  // Set a background color for the renderer and set the size of the
  // render window (expressed in pixels).
  aRenderer->SetBackground(0,0,0);
  renWin->SetSize(640, 480);

  // Note that when camera movement occurs (as it does in the Dolly()
  // method), the clipping planes often need adjusting. Clipping planes
  // consist of two planes: near and far along the view direction. The 
  // near plane clips out objects in front of the plane; the far plane
  // clips out objects behind the plane. This way only what is drawn
  // between the planes is actually rendered.
  aRenderer->ResetCameraClippingRange ();
  char dir[128];
  
  while(1)
  {
    for( int i = 0; i < 1000; i+=10)
    {
      sprintf(dir,"/home/abergman/Research/DCell/temp/result.%d.Real64",i);
      reader->SetFileName(dir);
      renWin->Render();
    }
  }

  // Initialize the event loop and then start it.
  iren->Initialize();
  iren->Start(); 
  
  vtkRIBExporter *rib = vtkRIBExporter::New();
    rib->SetFilePrefix("myrib");
    rib->SetRenderWindow(renWin);
    rib->Write();
    rib->Delete();
  
  // It is important to delete all objects created previously to prevent
  // memory leaks. In this case, since the program is on its way to
  // exiting, it is not so important. But in applications it is
  // essential.
  aCamera->Delete();
  iren->Delete();
  renWin->Delete();
  aRenderer->Delete();

  return 0;
}