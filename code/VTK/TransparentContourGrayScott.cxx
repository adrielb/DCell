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
//#undef OFFSCREEN

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
  int d1=100, d2=100, d3=25, ddd = d1*d2*d3, max=65000;
  char *inputDataName= argv[1],
       *outputImageName = argv[2];

  /* 
   * Indicate that we want to use Mesa Offscreen
   */
  vtkGraphicsFactory *factGraphics = vtkGraphicsFactory::New();
    factGraphics->SetUseMesaClasses(1);
    factGraphics->Delete();

  vtkImagingFactory *factImage = vtkImagingFactory::New();
    factImage->SetUseMesaClasses(1);
    factImage->Delete();

  vtkImageReader2 *reader = vtkImageReader2::New();
    reader->SetFileName(inputDataName);
    reader->SetDataScalarTypeToDouble();
    reader->SetDataExtent(0,d1-1,0,d2-1,0,d3-1);
    reader->SetDataSpacing(1,1,1);
    reader->SetFileDimensionality(3);
    reader->SetHeaderSize(2*sizeof(int)); //For reading Petsc Vecs
    reader->SetNumberOfScalarComponents(2);
  
  vtkImageExtractComponents *extractImage = vtkImageExtractComponents::New(); 
    extractImage->SetInputConnection( reader->GetOutputPort() );
    extractImage->SetComponents( 1 );
  
  vtkImageShiftScale *shiftScale = vtkImageShiftScale::New(); 
    shiftScale->SetInputConnection( extractImage->GetOutputPort() );
    shiftScale->SetShift( 0 );
    shiftScale->SetScale( max ); // max / (double)ddd
    shiftScale->SetOutputScalarTypeToUnsignedShort();
      
  vtkCamera *aCamera = vtkCamera::New();
    aCamera->SetViewUp(0, 0, 1);
    aCamera->SetPosition(d1/3,d2,3*d3);
    aCamera->SetFocalPoint(d1/2, d2/2, d3/2);
    aCamera->ComputeViewPlaneNormal();
    aCamera->SetViewAngle( 45 );
  
  vtkMarchingCubes *contours = vtkMarchingCubes::New();
    contours->SetInputConnection( extractImage->GetOutputPort() );
    //    contours->GenerateValues(5, .1, 1);
    contours->SetValue(0,.1);
    
  vtkDepthSortPolyData *depthSort = vtkDepthSortPolyData::New();
    depthSort->SetInputConnection( contours->GetOutputPort());
    depthSort->SetDirectionToBackToFront();
    depthSort->SetCamera(aCamera);
    
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection( depthSort->GetOutputPort() );
    mapper->ScalarVisibilityOff();
  
  vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);
    //    actor->GetProperty()->SetOpacity(0.1);
    actor->GetProperty()->SetColor(0,1,0);

#ifdef OFFSCREEN
/* 
 * OFF-SCREEN MESA RENDERING
 */
  vtkMesaRenderer *ren = vtkMesaRenderer::New();
  
  vtkXMesaRenderWindow *renWin = vtkXMesaRenderWindow::New();
    renWin->AddRenderer(ren);
    renWin->OffScreenRenderingOn();
#else
/*
 * ON-SCREEN RENDERING
 */
  vtkRenderer *ren = vtkRenderer::New();
  
  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren);
#endif

    ren->AddActor( actor );
    ren->SetActiveCamera(aCamera);
    ren->SetBackground(1,1,1);  // white  
    ren->SetBackground(0,0,0);  // black
    
    renWin->SetSize(1280, 720);    // HD
    renWin->SetSize(800, 600);   // Projector full-screen    
    renWin->Render();

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
