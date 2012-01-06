#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRawImageIO.h"
#include "itkImage.h"
#include "itkConnectedComponentImageFilter.h"

int main( int argc, char ** argv)
{
  const   unsigned int    Dimension = 2;
  typedef unsigned char InputPixelType;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;
  typedef unsigned char OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::ImageFileReader< InputImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[1] );
  
  typedef itk::RawImageIO< InputPixelType, Dimension> MyRawType;
    MyRawType::Pointer raw = MyRawType::New();
    raw->SetDimensions(0, 421);
    raw->SetDimensions(1, 451);
    raw->SetByteOrderToLittleEndian();
    reader->SetImageIO( raw );

  typedef itk::ConnectedComponentImageFilter< InputImageType, OutputImageType > ConnectedFilterType;
    ConnectedFilterType::Pointer filter = ConnectedFilterType::New();
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
