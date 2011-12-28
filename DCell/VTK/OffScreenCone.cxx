#include "vtkXMesaRenderWindow.h"
#include "vtkMesaRenderer.h"
#include "vtkConeSource.h"
#include "vtkMesaPolyDataMapper.h"
#include "vtkMesaActor.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

int main(int argc, char **argv)
{
  vtkXMesaRenderWindow *renWin = vtkXMesaRenderWindow::New();
    renWin->OffScreenRenderingOn();
    
  vtkMesaRenderer *ren = vtkMesaRenderer::New();
    renWin->AddRenderer( ren );
    
  vtkConeSource *cone = vtkConeSource::New();

  vtkMesaPolyDataMapper *mp = vtkMesaPolyDataMapper::New();
    mp->SetInputConnection( cone->GetOutputPort() );
  
  vtkMesaActor *actor = vtkMesaActor::New();
    actor->SetMapper(mp);

  ren->AddActor(actor);
  renWin->Render();
  
  vtkWindowToImageFilter *w2if = vtkWindowToImageFilter::New();
    w2if->SetInput(renWin);

  vtkPNGWriter *wr = vtkPNGWriter::New();
    wr->SetInputConnection( w2if->GetOutputPort() );
    wr->SetFileName( "MangledMesaTest.png" );
    wr->Write();
}