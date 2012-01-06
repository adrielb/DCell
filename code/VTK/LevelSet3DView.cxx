#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkImageReader2.h"
#include "vtkImageShiftScale.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeRayCastMIPFunction.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkVolume.h"
#include "vtkWindowToImageFilter.h"
#include "vtkJPEGWriter.h"
#include "vtkCamera.h"
#include "vtkImageExtractComponents.h"
#include "vtkGraphicsFactory.h"
#include "vtkImagingFactory.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkDepthSortPolyData.h"

#define OFFSCREEN
#undef OFFSCREEN

#ifdef OFFSCREEN
#include "vtkMesaRenderer.h"
#include "vtkXMesaRenderWindow.h"
#endif

/* TODO: need GUI for:
 *        * Camera angle
 *        * Color transfer function
 *        * Opaque transfer function
 */

int main(int argc, char **argv)
{
  int d1=100, d2=d1, d3=d1, ddd = d1*d2*d3;
  char *inputDataName= argv[1],
       *outputImageName = argv[2];

  /* 
   * Indicate that we want to use Mesa Offscreen
   */
#ifdef OFFSCREEN  
  vtkGraphicsFactory *factGraphics = vtkGraphicsFactory::New();
    factGraphics->SetUseMesaClasses(1);
    factGraphics->Delete();

  vtkImagingFactory *factImage = vtkImagingFactory::New();
    factImage->SetUseMesaClasses(1);
    factImage->Delete();
#endif
    
  vtkImageReader2 *reader = vtkImageReader2::New();
    reader->SetFileName(inputDataName);
    reader->SetDataScalarTypeToDouble();
    reader->SetDataExtent(0,d1-1,0,d2-1,0,d3-1);
    reader->SetDataSpacing(1,1,1);
    reader->SetFileDimensionality(3);
    reader->SetHeaderSize(2*sizeof(int)); //For reading Petsc Vecs
      
  vtkCamera *aCamera = vtkCamera::New();
    aCamera->SetViewUp(0, 0, 1);
    aCamera->SetPosition(d1/3,d2,3*d3);
    aCamera->SetFocalPoint(d1/2, d2/2, d3/2);
    aCamera->ComputeViewPlaneNormal();
    aCamera->SetViewAngle( 45 );
  
  vtkMarchingCubes *contours = vtkMarchingCubes::New();
    contours->SetInputConnection( reader->GetOutputPort() );
    //    contours->GenerateValues(5, .1, 1);
    contours->SetValue(0,0.);
    
  vtkDepthSortPolyData *depthSort = vtkDepthSortPolyData::New();
    depthSort->SetInputConnection( contours->GetOutputPort());
    depthSort->SetDirectionToBackToFront();
    depthSort->SetCamera(aCamera);
    
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection( depthSort->GetOutputPort() );
    mapper->ScalarVisibilityOff();
  
  vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetOpacity(0.4);
    actor->GetProperty()->SetColor(0,1,0); 

#ifdef OFFSCREEN
/* 
 * OFF-SCREEN MESA RENDERING
 */
  vtkMesaRenderer *ren = vtkMesaRenderer::New();
  
  vtkXMesaRenderWindow *renWin = vtkXMesaRenderWindow::New();
    renWin->OffScreenRenderingOn();
#else
/*
 * ON-SCREEN RENDERING
 */
  vtkRenderer *ren = vtkRenderer::New();
  
  vtkRenderWindow *renWin = vtkRenderWindow::New();    
#endif

    ren->AddActor( actor );
    ren->SetActiveCamera(aCamera);
    ren->SetBackground(1,1,1);  // white  
    ren->SetBackground(0,0,0);  // black

    renWin->AddRenderer(ren);
    renWin->SetSize(1280, 720);    // HD
    renWin->SetSize(800, 600);   // Projector full-screen    
//    renWin->Render();

#ifndef OFFSCREEN 
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    iren->Initialize();
    iren->Start();
#endif
    
  vtkWindowToImageFilter *w2if = vtkWindowToImageFilter::New();
    w2if->SetInput(renWin);
    
  vtkJPEGWriter *writer = vtkJPEGWriter::New();
    writer->SetInputConnection(w2if->GetOutputPort());
    writer->SetFileName(outputImageName);
    writer->Write();
}
