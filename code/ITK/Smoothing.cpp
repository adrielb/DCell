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
  const unsigned int  Dimension2 = 2;
  const unsigned int  Dimension3 = 3;
  const char * inputFilename  = argv[1];
  const char * outputFilename = argv[2];
  const unsigned int numberOfIterations = atoi( argv[3] );
  const double       timeStep = 0.125; //Recommended for 2D images

  typedef itk::Image< PixelType, Dimension2 > ImageType2D;
  typedef itk::Image< PixelType, Dimension3 > ImageType3D;
  
  typedef itk::ImageFileReader< ImageType3D >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFilename  );
  reader->Update();
  //TODO: since joinseriesimagefilter has problems, try TileImageFilter or PasteImageFilter
  typedef itk::JoinSeriesImageFilter< ImageType2D, ImageType3D > JoinSeriesType;
  JoinSeriesType::Pointer joinseries = JoinSeriesType::New();
  
  ImageType3D::RegionType inputRegion = reader->GetOutput()->GetLargestPossibleRegion();
  ImageType3D::SizeType size = inputRegion.GetSize();
  unsigned int max = size[2];
  size[2] = 0;
  ImageType3D::IndexType start = inputRegion.GetIndex();
  ImageType3D::RegionType desiredRegion;
  desiredRegion.SetSize( size );
  for (unsigned int i = 0; i < max; i++)
  {
//    std::cout << "i: " << i << std::endl;
    /* TODO: Only first image in stack is filtered and then repeated
     * for every time point. Allocating new filters for every slice 
     * solves this but is inefficient, and calling delete causes a seg fault.
     * How to allocate only once?
     */  
    typedef itk::ExtractImageFilter< ImageType3D, ImageType2D > ExtractFilterType;
    ExtractFilterType::Pointer extractfilter = ExtractFilterType::New();
    extractfilter->SetInput( reader->GetOutput() );
    start[2] = i;
    desiredRegion.SetIndex( start );
//    std::cout << i << ": " << desiredRegion << std::endl;
    extractfilter->SetExtractionRegion( desiredRegion );
    extractfilter->Update();
    typedef itk::CurvatureFlowImageFilter< ImageType2D, ImageType2D >  FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetNumberOfIterations( numberOfIterations );
    filter->SetTimeStep( timeStep );
    filter->SetInput( extractfilter->GetOutput() );
    filter->Update();
    joinseries->PushBackInput( filter->GetOutput() );
  }
  
  typedef itk::ImageFileWriter< ImageType3D >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFilename );
  writer->SetInput( joinseries->GetOutput() );

  try { 
    writer->Update(); 
  } catch( itk::ExceptionObject & err ) { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
