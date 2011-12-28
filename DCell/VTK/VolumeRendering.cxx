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

int main(int argc, char **argv)
{
  int d1=300, d2=d1, d3=d2, ddd = d1*d2*d3, max=65000;
  char *inputDataName= argv[1],
       *outputImageName = argv[2];

  /* Indicate that we want to use Mesa Offscreen
  vtkGraphicsFactory *factGraphics = vtkGraphicsFactory::New();
    factGraphics->SetUseMesaClasses(1);
    factGraphics->SetOffScreenOnlyMode(1);

  vtkImagingFactory *factImage = vtkImagingFactory::New();
    factImage->SetUseMesaClasses(1);
  */

  vtkRenderer *aRenderer = vtkRenderer::New();
    aRenderer->SetBackground(0,0,0);  // black
    aRenderer->SetBackground(1,1,1);  // white
    
  vtkImageReader2 *reader = vtkImageReader2::New();
    reader->SetFileName(inputDataName);
    reader->SetDataScalarTypeToDouble();
//    reader->SetFileName("/home/abergman/Research/DCell/temp/mat.Real64");
    reader->SetDataByteOrderToBigEndian();
    reader->SetHeaderSize(2*sizeof(int)); //For reading Petsc Vecs
    reader->SetDataExtent(0,d1-1,0,d2-1,0,d3-1);
    reader->SetDataSpacing(1,1,1);
    reader->SetFileDimensionality(3);
    reader->SetNumberOfScalarComponents(2);

  vtkImageExtractComponents *extractImage = vtkImageExtractComponents::New();
    extractImage->SetInputConnection( reader->GetOutputPort() );
    extractImage->SetComponents( 0 );
  
  vtkImageShiftScale *shiftScale = vtkImageShiftScale::New(); 
    shiftScale->SetInputConnection( extractImage->GetOutputPort() );
    shiftScale->SetShift( 0 );
    shiftScale->SetScale( 60000 ); // max / (double)ddd
    shiftScale->SetOutputScalarTypeToUnsignedShort();
  
  vtkPiecewiseFunction *opacityTransfer = vtkPiecewiseFunction::New(); 
    opacityTransfer->AddPoint( 10, 0.0);
    opacityTransfer->AddPoint( 400, 1.0);

  vtkPiecewiseFunction *grayTransfer = vtkPiecewiseFunction::New();
    grayTransfer->AddPoint(0, 0.0);
    grayTransfer->AddPoint(max, 1.0);
    
  vtkColorTransferFunction *colorTransfer = vtkColorTransferFunction::New(); 
    colorTransfer->AddRGBPoint( 0  , 1.0, 1.0, 1.0);
    colorTransfer->AddRGBPoint( 1.1, 0.0, 0.0, 1.0);

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
    volumeMapper->SetVolumeRayCastFunction( compositeFunction );
    volumeMapper->SetSampleDistance(1);
    volumeMapper->SetInputConnection(shiftScale->GetOutputPort());
    
  vtkVolume *volume = vtkVolume::New();
    volume->SetMapper(volumeMapper);
    volume->SetProperty(volumeProperty);
    
  vtkCamera *aCamera = aRenderer->GetActiveCamera();
    aCamera->SetViewUp(0, 0, 1);
    aCamera->SetPosition(2*d1,2*d2,2*d3);
    aCamera->SetFocalPoint(d1/2, d2/2, d3/2);
    aCamera->ComputeViewPlaneNormal();
    aCamera->SetViewAngle( 45 );
  
  aRenderer->AddVolume(volume);
 
  vtkRenderWindow *renWin = vtkRenderWindow::New();
//    renWin->SetOffScreenRendering(1);
    renWin->AddRenderer(aRenderer);
//    renWin->SetSize(1280, 720);    // HD
    renWin->SetSize(800, 600);     // SD

// ON SCREEN
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    iren->Initialize();
    iren->Start();

    renWin->Render();


  vtkWindowToImageFilter *w2if = vtkWindowToImageFilter::New();
    w2if->SetInput(renWin);
    
  vtkJPEGWriter *writer = vtkJPEGWriter::New();
    writer->SetInputConnection(w2if->GetOutputPort());
    writer->SetFileName(outputImageName);
    writer->Write();
}
