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
#include "vtkMesaRenderer.h"
#include "vtkXMesaRenderWindow.h"
#include "vtkGraphicsFactory.h"
#include "vtkImagingFactory.h"

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
  
  vtkMesaRenderer *ren = vtkMesaRenderer::New();
    ren->SetBackground(0,0,0);  // black

  /*
   * On-screen 
  vtkRenderer *ren = vtkRenderer::New();
    ren->SetBackground(1,1,1);  // white  
    ren->SetBackground(0,0,0);  // black
  */    
    
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
  
  vtkPiecewiseFunction *opacityTransfer = vtkPiecewiseFunction::New(); 
    opacityTransfer->AddPoint(  5000, 0.0);
    opacityTransfer->AddPoint( 30000, 1.0);

  vtkPiecewiseFunction *grayTransfer = vtkPiecewiseFunction::New();
    grayTransfer->AddPoint(0, 0.0);
    grayTransfer->AddPoint(max, 1.0);
    
  vtkColorTransferFunction *colorTransfer = vtkColorTransferFunction::New(); 
    colorTransfer->AddRGBPoint( 1000, 1.0, 1.0, 1.0);
    colorTransfer->AddRGBPoint( 1001, 0.0, 1.0, 0.0);

  vtkVolumeProperty *volumeProperty = vtkVolumeProperty::New();
    volumeProperty->SetColor( colorTransfer );
    volumeProperty->SetScalarOpacity(opacityTransfer);
    volumeProperty->ShadeOn();
    volumeProperty->SetInterpolationTypeToLinear();
//    volumeProperty->SetInterpolationTypeToNearest();
    
  vtkVolumeRayCastCompositeFunction *compositeFunction = vtkVolumeRayCastCompositeFunction::New();
    compositeFunction->SetCompositeMethodToClassifyFirst();
//    compositeFunction->SetCompositeMethodToInterpolateFirst();
    
  vtkVolumeRayCastMIPFunction *mipFunction = vtkVolumeRayCastMIPFunction::New();
   
  vtkVolumeRayCastMapper *volumeMapper = vtkVolumeRayCastMapper::New();
    volumeMapper->SetInputConnection(shiftScale->GetOutputPort());  
    volumeMapper->SetVolumeRayCastFunction( mipFunction );
    volumeMapper->SetSampleDistance(.1);
    
  vtkVolume *volume = vtkVolume::New();
    volume->SetMapper(volumeMapper);
    volume->SetProperty(volumeProperty);
    
  vtkCamera *aCamera = ren->GetActiveCamera();
    aCamera->SetViewUp(0, 0, 1);
    aCamera->SetPosition(d1/3,d2,3*d3);
    aCamera->SetFocalPoint(d1/2, d2/2, d3/2);
    aCamera->ComputeViewPlaneNormal();
    aCamera->SetViewAngle( 45 );
  
  
    ren->AddVolume(volume);

/* 
 * OFF-SCREEN MESA RENDERING
 */
  vtkXMesaRenderWindow *renWin = vtkXMesaRenderWindow::New();
    renWin->AddRenderer(ren);
    renWin->OffScreenRenderingOn();
    renWin->SetSize(1280, 720);    // HD
/* */
    
    
/*
 * ON-SCREEN RENDERING
 *
  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren);
//    renWin->SetSize(1280, 720);    // HD
    renWin->SetSize(800, 600);   // Projector full-screen


  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    iren->Initialize();
    iren->Start();
  */
    
    renWin->Render();
    
  vtkWindowToImageFilter *w2if = vtkWindowToImageFilter::New();
    w2if->SetInput(renWin);
    
  vtkJPEGWriter *writer = vtkJPEGWriter::New();
    writer->SetInputConnection(w2if->GetOutputPort());
    writer->SetFileName(outputImageName);
    writer->Write();
}
