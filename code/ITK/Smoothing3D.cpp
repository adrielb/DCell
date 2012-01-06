#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkJoinSeriesImageFilter.h"
#include <iostream>
 
int main( int argc, char ** argv)
{
  typedef double      PixelType;
  const unsigned int  Dimension3 = 3;
  const char * inputFilename  = argv[1];
  const char * outputFilename = argv[2];
  const unsigned int numberOfIterations = atoi( argv[3] );
  const double       timeStep = 0.0625; //Recommended for 3D images

  typedef itk::Image< PixelType, Dimension3 > ImageType3D;
  
  typedef itk::ImageFileReader< ImageType3D >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFilename  );
  
  typedef itk::CurvatureFlowImageFilter< ImageType3D, ImageType3D >  FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetTimeStep( timeStep );
  filter->SetInput( reader->GetOutput() );
    
  typedef itk::ImageFileWriter< ImageType3D >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFilename );
  writer->SetInput( filter->GetOutput() );
  
  try { 
    writer->Update(); 
  } catch( itk::ExceptionObject & err ) { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
