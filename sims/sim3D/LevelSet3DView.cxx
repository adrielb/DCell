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
#include "vtkGlobFileNames.h"
#include <sstream>
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkObjectFactory.h"
#include "vtkClipPolyData.h"
#include "vtkPlane.h"
#include "vtkLight.h"
#include "vtkImplicitPlaneWidget.h"
using namespace std;

void ReadGrooves( string inputDataDir, double dx, vtkActor **grooves );

void ReadSizeFile( string filename, vector<int>& vec ) {
  int temp;
  ifstream ifs( filename.c_str() );
  vec.clear();
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
  vec.clear();
  while( ifs ) {
    ifs >> temp;
    vec.push_back(dx*temp);
  }
  cout << filename << "  " << vec.size()  << endl;
}

class vtkSliderCallback : public vtkCommand {
public:
  vtkSliderCallback():reader(0) {
    globFiles = vtkGlobFileNames::New();
  }

  static vtkSliderCallback *New() {
    return new vtkSliderCallback;
  }

  virtual void Execute(vtkObject *caller, unsigned long, void*) {
    vtkSliderWidget *sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
    int idx = static_cast<int>(static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue());
    UpdateFileIndex( idx );
  }

  const char* MakeFileName( string prefix, int idx, string suffix ) {
    osFileName.clear();
    osFileName.str("");
    osFileName << inputDataDir << "/" << prefix << "." << (idx+10000) << "." << suffix;
    cout << osFileName.str() << endl;
    return osFileName.str().c_str();
  }

  void UpdateFileIndex( int idx ) {
    const char* inputFileName = MakeFileName( "mycell", idx, "Real64" );
    reader->SetFileName( inputFileName );
    reader->SetDataExtent( &phisize[6*idx] );
    reader->SetDataOrigin( &phipos [3*idx] );
  }

  void SetDataDir( string datadir ) {
    inputDataDir = datadir;
    UpdateEndIndex();
  }

  void UpdateEndIndex() {
    globFiles->Reset();
    globFiles->SetDirectory( inputDataDir.c_str() );
    globFiles->AddFileNames( "mycell.*.Real64" );
    int endIdx = globFiles->GetNumberOfFileNames() - 1;
    cout << "Last File: " << globFiles->GetNthFileName(endIdx) << endl;
    ReadSizeFile( inputDataDir + "/mycell.size", phisize );
    ReadPosFile(  inputDataDir + "/mycell.pos", dx, phipos );
    sliderRep->SetMaximumValue(endIdx);
    sliderRep->SetValue(endIdx);
    UpdateFileIndex( endIdx );
  }

  void MakeMovie() {
    iren->Disable();
    int endIdx = sliderRep->GetMaximumValue();
    for (int idx = 0; idx < endIdx; ++idx) {
      UpdateFileIndex(idx);
      renWin->Render();
      const char* outputFileName = MakeFileName("g", idx, "png");
      writer->SetFileName(outputFileName);
      writer->Write();
    }
  }

  void ShiftFrame( int delta ) {
    int idx = sliderRep->GetValue() + delta;
    int max = sliderRep->GetMaximumValue();
    int min = sliderRep->GetMinimumValue();
    if( idx < min ) idx = min;
    if( idx > max ) idx = max;
    sliderRep->SetValue( idx );
    UpdateFileIndex( idx );
    renWin->Render();
  }

  vtkImageReader2 *reader;
  vtkImageWriter *writer;
  vtkSliderRepresentation2D *sliderRep;
  vtkRenderWindowInteractor *iren;
  vtkRenderWindow *renWin;
  double dx;

private:
  ostringstream osFileName;
  string inputDataDir;
  vtkGlobFileNames *globFiles;
  vector<int> phisize;
  vector<double> phipos;
};

class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  vtkSliderCallback *callback;

  static KeyPressInteractorStyle* New();
  vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

  virtual void OnKeyPress()
  {
    vtkRenderWindowInteractor *rwi = this->Interactor;
    std::string key = rwi->GetKeySym();
//    std::cout << "Pressed " << key << std::endl;

    if(key == "Left")  callback->ShiftFrame(-1);
    if(key == "Right") callback->ShiftFrame( 1);
    if(key == "Down")  callback->ShiftFrame(-10);
    if(key == "Up")    callback->ShiftFrame( 10);
    if(key == "u")     callback->UpdateEndIndex();
    if(key == "m")     callback->MakeMovie();

    // Forward events
    vtkInteractorStyleTrackballCamera::OnKeyPress();
  }
};
vtkStandardNewMacro(KeyPressInteractorStyle);

// /data/dir dx 1999
int main(int argc, char **argv)
{
  float d1=20, d2=20, d3=9;
  string inputDataDir = argv[1];
  float dx = atof( argv[2] );

  vtkSliderRepresentation2D *sliderRep = vtkSliderRepresentation2D::New();
  vtkImageReader2 *reader = vtkImageReader2::New();

  vtkSliderCallback *callback = vtkSliderCallback::New();
    callback->dx = dx;
    callback->sliderRep = sliderRep;
    callback->reader = reader;
    callback->SetDataDir( inputDataDir );

    reader->SetDataScalarTypeToDouble();
    reader->SetDataSpacing(dx,dx,dx);
    reader->SetFileDimensionality(3);
    reader->SwapBytesOn();
    reader->SetHeaderSize(2*sizeof(int)); //For reading Petsc Vecs
  vtkMarchingCubes *contour = vtkMarchingCubes::New();
    contour->SetInputConnection( reader->GetOutputPort() );
    contour->SetValue(0,0.);
  vtkPlane *clipPlane = vtkPlane::New();
    clipPlane->SetOrigin(d1/2,d2/2,1);
    clipPlane->SetNormal(0,-1,0);
//    clipPlane->SetNormal(0,0,-1);
  vtkClipPolyData *clip = vtkClipPolyData::New();
    clip->SetInputConnection(contour->GetOutputPort());
    clip->SetClipFunction( clipPlane );
  vtkImplicitPlaneWidget *planeWidget = vtkImplicitPlaneWidget::New();
    planeWidget->SetInput(clip->GetOutput());
    planeWidget->PlaceWidget();
    planeWidget->SetPlaceFactor(10);
    planeWidget->AddObserver(  vtkCommand::InteractionEvent, callback);
    planeWidget->GetPlane( clipPlane );

/* only for transparency
  vtkDepthSortPolyData *depthSort = vtkDepthSortPolyData::New();
    depthSort->SetInputConnection( contours->GetOutputPort());
    depthSort->SetDirectionToBackToFront();
    depthSort->SetCamera(aCamera);
*/
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection( clip->GetOutputPort() );
    mapper->ScalarVisibilityOff();

  vtkCamera *camOrtho = vtkCamera::New();
    camOrtho->SetViewUp(0, 0, 1);
    camOrtho->ComputeViewPlaneNormal();
    camOrtho->SetViewAngle( 85 );
    double r = 7;
    camOrtho->SetPosition(0.5*r+d1/2, 0.5*r+d2/2, 0.71*r);
    camOrtho->SetFocalPoint(d1/2,d2/2,0);

  vtkActor *actorMyCell = vtkActor::New();
    actorMyCell->SetMapper(mapper);
    actorMyCell->GetProperty()->SetColor(0,1,0);
  vtkActor *actorGrooves = vtkActor::New();
    ReadGrooves(inputDataDir, dx, &actorGrooves );
  vtkCubeAxesActor *actorAxes = vtkCubeAxesActor::New();
    actorAxes->SetCamera(camOrtho);
    double b[6] = {2,16,2,16,3,8};
    actorAxes->SetBounds(b);

  vtkRenderer *ren = vtkRenderer::New();
    ren->AddActor( actorAxes );
    ren->AddActor( actorMyCell );
    ren->AddActor( actorGrooves );
    ren->SetActiveCamera(camOrtho);
    ren->SetBackground(0,0,0);  // black
  vtkLight *light = vtkLight::New();
    light->SetPosition(0,0,d3);
    ren->AddLight( light );
  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->SetSize(800, 600);   // Projector full-screen
    renWin->AddRenderer(ren);
    callback->renWin = renWin;
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    callback->iren = iren;
    planeWidget->SetInteractor(iren);
    renWin->Render();


    sliderRep->SetMinimumValue(0);
    sliderRep->SetTitleText("Index");
    sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint1Coordinate()->SetValue(0.01 ,.1);
    sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint2Coordinate()->SetValue(0.99, .1);

  vtkSliderWidget *sliderWidget = vtkSliderWidget::New();
    sliderWidget->SetInteractor(iren);
    sliderWidget->SetRepresentation(sliderRep);
    sliderWidget->SetAnimationModeToAnimate();
    sliderWidget->EnabledOn();

  KeyPressInteractorStyle *style = KeyPressInteractorStyle::New();
    style->callback = callback;
    style->SetCurrentRenderer( ren );
    iren->SetInteractorStyle( style );
    sliderWidget->AddObserver( vtkCommand::InteractionEvent, callback);

  vtkWindowToImageFilter *w2if = vtkWindowToImageFilter::New();
    w2if->SetInput(renWin);
  vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputConnection(w2if->GetOutputPort());
    callback->writer = writer;

  iren->Initialize();
  renWin->Render();
  iren->Start();
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
  vtkMarchingCubes *contour = vtkMarchingCubes::New();
    contour->SetInputConnection( reader->GetOutputPort() );
    contour->SetValue(0,0.);
  vtkPlane *clipPlane = vtkPlane::New();
    clipPlane->SetOrigin(10,10,1);
    clipPlane->SetNormal(0,0,-1);
  vtkClipPolyData *clip = vtkClipPolyData::New();
    clip->SetInputConnection(contour->GetOutputPort());
    clip->SetClipFunction( clipPlane );
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection( clip->GetOutputPort() );
    mapper->ScalarVisibilityOff();
  vtkActor *actorGrooves = vtkActor::New();
    actorGrooves->SetMapper(mapper);
    actorGrooves->GetProperty()->SetColor(grey,grey,grey);

  *grooves = actorGrooves;
}
