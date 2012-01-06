#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRawImageIO.h"
#include "itkCastImageFilter.h"
#include "itkVnlFFTRealToComplexConjugateImageFilter.h"
#include "itkVnlFFTComplexConjugateToRealImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkMultiplyImageFilter.h"

int main( int argc, char * argv [] )
{
  if( argc < 4 )
  {
    std::cerr << "Usage: " << argv[0] << " inputScalarImage  inputMaskImage";
    std::cerr << " outputFilteredImage" << std::endl;
  }
// Import image  
  const   unsigned int   timepoints = 128; //atoi( argv[4] );
  const   unsigned int   Dimension = 3;
  typedef unsigned short CascadePixelType;
  typedef itk::Image< CascadePixelType, Dimension > InputImageType;

  typedef itk::RawImageIO< CascadePixelType, Dimension> MyRawType;
    MyRawType::Pointer raw = MyRawType::New();
    raw->SetDimensions(0, 512);
    raw->SetDimensions(1, 512);
    raw->SetDimensions(2, timepoints);
    raw->SetFileDimensionality(Dimension);
    raw->SetByteOrderToLittleEndian();

  typedef itk::ImageFileReader< InputImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[1] );
    reader->SetImageIO( raw );
  
  typedef float FloatPixelType;
  typedef itk::Image< FloatPixelType, Dimension > FloatImageType;
  typedef itk::CastImageFilter< InputImageType, FloatImageType > CastFilterType;
    CastFilterType::Pointer filter = CastFilterType::New();
    filter->SetInput( reader->GetOutput() );
  
  typedef itk::VnlFFTRealToComplexConjugateImageFilter< 
      FloatPixelType, Dimension >  FFTFilterType;
    FFTFilterType::Pointer fftFilter = FFTFilterType::New();
    fftFilter->SetInput( filter->GetOutput() );
  typedef FFTFilterType::OutputImageType    SpectralImageType;
  
// Read in filter 
  typedef itk::RawImageIO< FloatPixelType, Dimension> RawFloatType;
    RawFloatType::Pointer rawfloat = RawFloatType::New();
    rawfloat->SetDimensions(0, 512);
    rawfloat->SetDimensions(1, 512);
    rawfloat->SetDimensions(2, timepoints);
    rawfloat->SetFileDimensionality(Dimension);
    rawfloat->SetByteOrderToLittleEndian();
  typedef itk::ImageFileReader< FloatImageType > FloatReaderType;
    FloatReaderType::Pointer floatreader = FloatReaderType::New();
    floatreader->SetFileName( argv[2] );
    floatreader->SetImageIO( rawfloat );
  
  typedef itk::MultiplyImageFilter< SpectralImageType, FloatImageType, SpectralImageType> MultFilterType;
    MultFilterType::Pointer multFilter = MultFilterType::New();
    multFilter->SetInput1( fftFilter->GetOutput() );
    multFilter->SetInput2( floatreader->GetOutput() );
  
  typedef itk::VnlFFTComplexConjugateToRealImageFilter< 
      FloatPixelType, Dimension >  IFFTFilterType;
    IFFTFilterType::Pointer fftInverseFilter = IFFTFilterType::New();
    fftInverseFilter->SetInput( multFilter->GetOutput() );

// Write filtered image
  typedef itk::RawImageIO< FloatPixelType, Dimension> FloatIOType;
    FloatIOType::Pointer floatio = FloatIOType::New();
    floatio->SetDimensions(0, 512);
    floatio->SetDimensions(1, 512);
    floatio->SetDimensions(2, timepoints);
    floatio->SetFileDimensionality(Dimension);
    floatio->SetByteOrderToLittleEndian();
  typedef itk::ImageFileWriter< FloatImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetImageIO( floatio );
  
  try
    {
    writer->SetFileName( argv[3] );
    writer->SetInput( fftInverseFilter->GetOutput() );
//    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error: " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

/* Debugging: obtaining real/imaginary parts */ 
  typedef itk::ComplexToRealImageFilter<
      SpectralImageType, FloatImageType > RealFilterType;
    RealFilterType::Pointer realFilter = RealFilterType::New();
    realFilter->SetInput( fftFilter->GetOutput() );

  typedef itk::ComplexToImaginaryImageFilter<
        SpectralImageType, FloatImageType > ImaginaryFilterType;
      ImaginaryFilterType::Pointer imaginaryFilter = ImaginaryFilterType::New();
      imaginaryFilter->SetInput( fftFilter->GetOutput() );
    
  try
    {
    writer->SetFileName( "import.Real32" );
    writer->SetInput( filter->GetOutput() );
    writer->Update();
    
    writer->SetFileName( "real.Real32" );
    writer->SetInput( realFilter->GetOutput() );
    writer->Update();
        
    writer->SetFileName( "imaginary.Real32" );
    writer->SetInput( imaginaryFilter->GetOutput() );
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error: " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
