#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRawImageIO.h"
#include "itkImage.h"
#include "itkRescaleIntensityImageFilter.h"

int main( int argc, char ** argv)
{
  typedef unsigned short CascadePixelType;
  typedef unsigned char   OutputPixelType;
  const   unsigned int    Dimension = 2;
  typedef itk::Image< CascadePixelType, Dimension > InputImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::ImageFileReader< InputImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[1] );
  
  typedef itk::RawImageIO< CascadePixelType, Dimension> MyRawType;
    MyRawType::Pointer raw = MyRawType::New();
    raw->SetDimensions(0, 512);
    raw->SetDimensions(1, 512);
    raw->SetByteOrderToLittleEndian();
    reader->SetImageIO( raw );

  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > RescaleFilterType;
    RescaleFilterType::Pointer filter = RescaleFilterType::New();
//    filter->SetOutputMinimum(0);
//    filter->SetOutputMaximum(255);
    filter->SetInput( reader->GetOutput() );
  
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[2] );
    writer->SetInput( filter->GetOutput());
  
  try { 
    writer->Update(); 
  } catch( itk::ExceptionObject & err ) { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}