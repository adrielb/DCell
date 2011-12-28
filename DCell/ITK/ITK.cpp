//============================================================================
// Name        : ITK.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRawImageIO.h"
#include "itkImage.h"
#include "itkFixedArray.h"
#include "itkScalarToArrayCastImageFilter.h"
#include <iostream>

#include "itkRescaleIntensityImageFilter.h"
#include "itkMRFImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkMinimumDecisionRule.h"
#include "itkImageClassifierBase.h"

/*
 * importing data from buffer
 * KDTree
 * Curvature smoothing
 * itk::MultiResolutionPyramidImageFilter
 * 
 */

int main()
{
  const unsigned int numberOfIterations = 50;
  const double       smoothingFactor    = 3;
  const unsigned int numberOfClasses    = 3;
  const double       means[3]           = {5156,10180,21197};
  
  const char * inputFilename  = "/home/abergman/workspace/ITK/image.a.1.gray";
  const char * outputFilename = "/home/abergman/workspace/ITK/image.a.1.tif";
  typedef unsigned short      PixelType;
  const unsigned int          Dimension = 2;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFilename  );
  typedef itk::RawImageIO<PixelType, Dimension> MyRawType;
  MyRawType::Pointer raw = MyRawType::New();
  raw->SetDimensions(0, 512);
  raw->SetDimensions(1, 512);
  raw->SetByteOrderToBigEndian();
  reader->SetImageIO( raw );
  
  typedef unsigned char       LabelPixelType;
  typedef itk::Image<LabelPixelType, Dimension > LabelImageType;
  typedef itk::FixedArray<LabelPixelType,1>  ArrayPixelType;
  typedef itk::Image< ArrayPixelType, Dimension > ArrayImageType;
  typedef itk::ScalarToArrayCastImageFilter< 
                     ImageType, ArrayImageType > ScalarToArrayFilterType;
  ScalarToArrayFilterType::Pointer 
    scalarToArrayFilter = ScalarToArrayFilterType::New();
  scalarToArrayFilter->SetInput( reader->GetOutput() );
  
  typedef itk::ImageClassifierBase< 
    ArrayImageType,LabelImageType > SupervisedClassifierType;
  SupervisedClassifierType::Pointer classifier = 
    SupervisedClassifierType::New();
  typedef itk::MinimumDecisionRule DecisionRuleType;
  DecisionRuleType::Pointer  classifierDecisionRule = DecisionRuleType::New();
  classifier->SetDecisionRule( classifierDecisionRule.GetPointer() );
  typedef itk::Statistics::DistanceToCentroidMembershipFunction< 
    ArrayPixelType > MembershipFunctionType;
  typedef MembershipFunctionType::Pointer MembershipFunctionPointer;
  double meanDistance = 0;
  vnl_vector<double> centroid(1); 
  for( unsigned int i=0; i < numberOfClasses; i++ )
    {
    MembershipFunctionPointer membershipFunction = 
      MembershipFunctionType::New();
    centroid[0] = means[i]; 
    membershipFunction->SetCentroid( centroid );
    classifier->AddMembershipFunction( membershipFunction );
    meanDistance += static_cast< double > (centroid[0]);
    }
  meanDistance /= numberOfClasses;
  std::vector< double > weights;
  weights.push_back(1.5);
  weights.push_back(2.0);
  weights.push_back(1.5);
  weights.push_back(2.0);
  weights.push_back(0.0); // This is the central pixel
  weights.push_back(2.0);
  weights.push_back(1.5);
  weights.push_back(2.0);
  weights.push_back(1.5);
  double totalWeight = 0;
  for(std::vector< double >::const_iterator wcIt = weights.begin(); 
    wcIt != weights.end(); ++wcIt )
    {
    totalWeight += *wcIt;
    }
  for(std::vector< double >::iterator wIt = weights.begin(); 
    wIt != weights.end(); wIt++ )
    {
    *wIt = static_cast< double > ( (*wIt) * meanDistance / (2 * totalWeight));
    }

  typedef itk::MRFImageFilter< ArrayImageType, LabelImageType > MRFFilterType;
  MRFFilterType::Pointer mrfFilter = MRFFilterType::New();
  mrfFilter->SetNumberOfClasses( numberOfClasses );
  mrfFilter->SetMaximumNumberOfIterations( numberOfIterations );
  mrfFilter->SetErrorTolerance( 1e-7 );
  mrfFilter->SetSmoothingFactor( smoothingFactor );
  mrfFilter->SetNeighborhoodRadius( 1 );
  mrfFilter->SetMRFNeighborhoodWeight( weights );
  mrfFilter->SetClassifier( classifier );
  mrfFilter->SetInput( scalarToArrayFilter->GetOutput() );
  
  typedef itk::ImageFileWriter< LabelImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFilename );
  writer->SetInput( mrfFilter->GetOutput() );
  try { 
    writer->Update(); 
  } catch( itk::ExceptionObject & err ) { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
  }
  
  std::cout << "Number of Iterations : ";
  std::cout << mrfFilter->GetNumberOfIterations() << std::endl;
  std::cout << "Stop condition: " << std::endl;
  std::cout << "  (1) Maximum number of iterations " << std::endl;
  std::cout << "  (2) Error tolerance:  "  << std::endl;
  std::cout << mrfFilter->GetStopCondition() << std::endl;

  return 0;
}
