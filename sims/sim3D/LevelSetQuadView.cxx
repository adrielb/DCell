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
#include "vtkPNGWriter.h"
#include "vtkCamera.h"
#include "vtkImageExtractComponents.h"
#include "vtkGraphicsFactory.h"
#include "vtkImagingFactory.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkDepthSortPolyData.h"
#include "vtkCubeAxesActor.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkSliderWidget.h"
#include "vtkCommand.h"
#include "vtkInteractorStyleSwitch.h"
#include <sstream>

/* TODO: need GUI for:
 *        * Camera angle
 *        * Color transfer function
 *        * Opaque transfer function
 */

using namespace std;

void ReadSizeFile( string filename, vector<int>& vec ) {
  int temp;
  ifstream ifs( filename.c_str() );
  while( ifs ) {
    vec.push_back(0);
    ifs >> temp;
    vec.push_back(temp-1);
  }
  cout << filename << "  " << vec.size() << endl;
}

void ReadPosFile( string filename, double dx, vector<double>& vec ) {
  double temp;
  ifstream ifs( filename.c_str() );
  while( ifs ) {
    ifs >> temp;
    vec.push_back(dx*temp);
  }
  cout << filename << "  " << vec.size()  << endl;
}

void ReadGrooves( string inputDataDir, double dx, vtkActor **grooves ) {
  int idx = 1;
  double grey = 0.5;
  string pre = "grooves";
  string inputFileName = inputDataDir + "/" + pre + ".10001.Real64";
  cout << inputFileName << endl;
  vector<int> size;
  vector<double> pos;
  ReadSizeFile( inputDataDir + "/" + pre + ".size", size );
  ReadPosFile(  inputDataDir + "/" + pre + ".pos", dx, pos );
  vtkImageReader2 *reader = vtkImageReader2::New();
    reader->SetFileName( inputFileName.c_str() );
    reader->SetDataScalarTypeToDouble();
    reader->SetDataExtent(&size[6*idx]);
    reader->SetDataOrigin(&pos[ 3*idx]);
    reader->SetDataSpacing(dx,dx,dx);
    reader->SetFileDimensionality(3);
    reader->SwapBytesOn();
    reader->SetHeaderSize(2*sizeof(int)); //For reading Petsc Vecs
  vtkMarchingCubes *contours = vtkMarchingCubes::New();
    contours->SetInputConnection( reader->GetOutputPort() );
    contours->SetValue(0,0.);
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection( contours->GetOutputPort() );
    mapper->ScalarVisibilityOff();
  vtkActor *actorGrooves = vtkActor::New();
    actorGrooves->SetMapper(mapper);
    actorGrooves->GetProperty()->SetColor(grey,grey,grey);

  *grooves = actorGrooves;
}

class vtkSliderCallback : public vtkCommand {
public:
  static vtkSliderCallback *New() {
    return new vtkSliderCallback;
  }

  virtual void Execute(vtkObject *caller, unsigned long, void*) {
    vtkSliderWidget *sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
    int idx = static_cast<int>(static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue());
    UpdateFileIndex( idx );
  }

  void UpdateFileIndex( int idx ) {
    inputFileName << inputDataDir << "/mycell." << (idx+10000) << ".Real64";
    cout << inputFileName.str() << endl;
    reader->SetFileName( inputFileName.str().c_str() );
    reader->SetDataExtent( &phisize[6*idx] );
    reader->SetDataOrigin( &phipos [3*idx] );
    inputFileName.clear();
    inputFileName.str("");
  }

  vtkSliderCallback():reader(0) {}
  vtkImageReader2 *reader;
  vector<int> phisize;
  vector<double> phipos;
  string inputDataDir;

private:
  ostringstream inputFileName;
};

// /data/dir dx 1999
int main(int argc, char **argv)
{
  float d1=20, d2=20, d3=9;
  string inputDataDir = argv[1];
  float dx = atof( argv[2] );
  char *fileIdx = argv[3];
  int   endIdx = atoi( fileIdx ) - 10000;

  string inputFileName = inputDataDir + "/mycell." + fileIdx + ".Real64";
  string outputImageName = inputDataDir + "/g." + fileIdx + ".png";
  cout << inputFileName << endl;
  cout << outputImageName << endl;
  cout << "idx = " << endIdx << endl;

  vector<int> phisize;
  vector<double> phipos;
  ReadSizeFile( inputDataDir + "/mycell.size", phisize );
  ReadPosFile(  inputDataDir + "/mycell.pos", dx, phipos );

  vtkImageReader2 *reader = vtkImageReader2::New();
    reader->SetFileName( inputFileName.c_str() );
    reader->SetDataScalarTypeToDouble();
    reader->SetDataExtent(&phisize[endIdx*6]);
    reader->SetDataOrigin(&phipos[3*endIdx]);
    reader->SetDataSpacing(dx,dx,dx);
    reader->SetFileDimensionality(3);
    reader->SwapBytesOn();
    reader->SetHeaderSize(2*sizeof(int)); //For reading Petsc Vecs
  
  vtkMarchingCubes *contours = vtkMarchingCubes::New();
    contours->SetInputConnection( reader->GetOutputPort() );
    contours->SetValue(0,0.);
/* only for transparency
  vtkDepthSortPolyData *depthSort = vtkDepthSortPolyData::New();
    depthSort->SetInputConnection( contours->GetOutputPort());
    depthSort->SetDirectionToBackToFront();
    depthSort->SetCamera(aCamera);
*/
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection( contours->GetOutputPort() );
    mapper->ScalarVisibilityOff();
  
  vtkCamera *camOrtho = vtkCamera::New();
    camOrtho->SetViewUp(0, 0, 1);
    camOrtho->ComputeViewPlaneNormal();
    camOrtho->SetViewAngle( 55 );
    camOrtho->SetPosition(1,1,d3-1);
    camOrtho->SetFocalPoint(d1/2,d2/2,0);

  vtkCamera *camOrtho2 = vtkCamera::New();
    camOrtho2->SetViewUp(0, 0, 1);
    camOrtho2->ComputeViewPlaneNormal();
    camOrtho2->SetViewAngle( 55 );
    camOrtho2->SetPosition(d1-1,d2-1,d3-1);
    camOrtho2->SetFocalPoint(d1/2,d2/2,0);

  vtkCamera *camCut = vtkCamera::New();
    camCut = vtkCamera::New();
    camCut->SetPosition(d1/2,0,5);
    camCut->SetFocalPoint(d1/2,d2/2,0);
    camCut->SetClippingRange(10,1000);

  vtkCamera *camTop = vtkCamera::New();
    double zpos = 800;
    camTop->SetPosition(  d1/2, d2/2, zpos);
    camTop->SetFocalPoint(d1/2, d2/2, 0);
    camTop->SetClippingRange(zpos-0.8,1000);
    camTop->ParallelProjectionOn();
    camTop->Zoom(0.1);

  vtkCamera *cameras[4] = {camTop, camCut, camOrtho, camOrtho2 };

  vtkActor *actorMyCell = vtkActor::New();
    actorMyCell->SetMapper(mapper);
    actorMyCell->GetProperty()->SetColor(0,1,0);
  vtkActor *actorGrooves = vtkActor::New();
    ReadGrooves(inputDataDir, dx, &actorGrooves );
  vtkCubeAxesActor *actorAxes = vtkCubeAxesActor::New();
    actorAxes->SetCamera(camOrtho);
    double b[6] = {2,16,2,16,3,8};
    actorAxes->SetBounds(b);

  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->SetSize(800, 600);   // Projector full-screen

  double xmins[4] = {0.0, 0.5, 0.0, 0.5};
  double xmaxs[4] = {0.5, 1.0, 0.5, 1.0};
  double ymins[4] = {0.0, 0.0, 0.5, 0.5};
  double ymaxs[4] = {0.5, 0.5, 1.0, 1.0};
  for(unsigned i = 0; i < 4; i++) {
    vtkRenderer *ren = vtkRenderer::New();
      ren->SetViewport(xmins[i],ymins[i],xmaxs[i],ymaxs[i]);
      ren->AddActor( actorAxes );
      ren->AddActor( actorMyCell );
      ren->AddActor( actorGrooves );
      ren->SetActiveCamera(cameras[i]);
      ren->SetBackground(0,0,0);  // black
      renWin->AddRenderer(ren);
  }

  vtkInteractorStyleSwitch *istyle = vtkInteractorStyleSwitch::New();
      istyle->SetCurrentStyleToTrackballCamera();
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    iren->SetInteractorStyle(istyle);
    renWin->Render();

  vtkSliderRepresentation2D *sliderRep = vtkSliderRepresentation2D::New();
    sliderRep->SetMinimumValue(0);
    sliderRep->SetMaximumValue(endIdx);
    sliderRep->SetValue(endIdx);
    sliderRep->SetTitleText("Index");
    sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint1Coordinate()->SetValue(0.01 , 0.03);
    sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint2Coordinate()->SetValue(0.45, 0.03);

  vtkSliderWidget *sliderWidget = vtkSliderWidget::New();
    sliderWidget->SetInteractor(iren);
    sliderWidget->SetRepresentation(sliderRep);
    sliderWidget->SetAnimationModeToAnimate();
    sliderWidget->EnabledOn();

  vtkSliderCallback *callback = vtkSliderCallback::New();
    callback->reader = reader;
    callback->inputDataDir = inputDataDir;
    callback->phipos = phipos;
    callback->phisize = phisize;
    sliderWidget->AddObserver( vtkCommand::InteractionEvent, callback);

  iren->Initialize();
  renWin->Render();
  iren->Start();

  vtkWindowToImageFilter *w2if = vtkWindowToImageFilter::New();
    w2if->SetInput(renWin);

  vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputConnection(w2if->GetOutputPort());
    writer->SetFileName(outputImageName.c_str());
    writer->Write();
}
