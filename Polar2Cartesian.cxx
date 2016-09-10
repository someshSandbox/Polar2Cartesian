#include "itkPolarToCartesianTransform.h"
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkResampleImageFilter.h>
 
// Dimension 
const unsigned int Dimension2D = 2;
// Anatomic pixel type
typedef unsigned int AnatomicPixelType;
typedef itk::Image< AnatomicPixelType, Dimension2D > AnatomicImageType;

/* Write Image to File */
template< typename TImage>
void WriteImage(std::string filename, typename TImage::Pointer image) {

    typedef itk::ImageFileWriter<TImage> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(image);
    writer->Update();
}

/* Create a polar image */
void CreatePolarImage(AnatomicImageType::PointType origin, 
                            double maxTheta, 
                            double thetaSpacing,
                            double radius,
                            double radiusSpacing,
                            void* pImg)
{

   AnatomicImageType::Pointer imageSource = AnatomicImageType::New();
        AnatomicImageType::IndexType start;
        start[0] = 0;
        start[1] = 0;
        AnatomicImageType::SizeType size;
        size[0] = (radius - origin[0]) / radiusSpacing;
        size[1] = (maxTheta - origin[1]) / thetaSpacing;
        AnatomicImageType::RegionType region;
        region.SetSize(size);
        region.SetIndex(start);
        AnatomicImageType::SpacingType spacingImage;
        spacingImage[0] = radiusSpacing;
        spacingImage[1] = thetaSpacing;

        imageSource->SetSpacing(spacingImage);
        imageSource->SetOrigin(origin);
        imageSource->SetRegions(region);
        imageSource->Allocate();

        itk::ImageRegionIteratorWithIndex<AnatomicImageType > it(imageSource, imageSource->GetLargestPossibleRegion().GetSize());
        it.GoToBegin();
        while (!it.IsAtEnd()) 
        {          
          AnatomicImageType::IndexType index = it.GetIndex();
          it.Set(index[0]);
          ++it;
        }

        AnatomicImageType::Pointer *outimage = (AnatomicImageType::Pointer *)pImg;
        *outimage = imageSource;


   } 


/* Create a cartesian image */
void CreateCartesianImage(AnatomicImageType::PointType origin, 
                            double maxTheta, 
                            double radius,
                            double XSpacing,                            
                            double YSpacing,
                            void* pImg)
{

   AnatomicImageType::Pointer imageSource = AnatomicImageType::New();
        AnatomicImageType::IndexType start;
        start[0] = 0;
        start[1] = 0;
        AnatomicImageType::SizeType size;
        size[0] = (2.0 * radius) / XSpacing;
        size[1] = (2.0 * radius) / YSpacing;
        AnatomicImageType::RegionType region;
        region.SetSize(size);
        region.SetIndex(start);
        AnatomicImageType::SpacingType spacingImage;
        spacingImage[0] = XSpacing;
        spacingImage[1] = YSpacing;

        AnatomicImageType::PointType originImage;
        originImage[0] = origin[0] - radius;
        originImage[1] = origin[1] - radius;
        imageSource->SetSpacing(spacingImage);
        imageSource->SetOrigin(originImage);
        imageSource->SetRegions(region);
        imageSource->Allocate();

        itk::ImageRegionIteratorWithIndex<AnatomicImageType > it(imageSource, imageSource->GetLargestPossibleRegion().GetSize());
        it.GoToBegin();
        while (!it.IsAtEnd()) 
        {          
          AnatomicImageType::IndexType index = it.GetIndex();
          AnatomicImageType::PointType point;
          imageSource->TransformIndexToPhysicalPoint(index, point); 
          double r = sqrt(point[0]*point[0] + point[1]*point[1]);
    		  // Theta
    		  double radianToDegrees = 180.0 / vnl_math::pi ;
    		  double theta = (std::atan2(point[1], point[0])) * radianToDegrees;
    		  if(theta < 0)
    		  	theta += 360;

          // Set pixels inside given radius to R, and outside radius to 0
    		  if((r < radius) && ( theta <= maxTheta))
                it.Set(r);
              else
            	it.Set(0);

              ++it;
            }

        AnatomicImageType::Pointer *outimage = (AnatomicImageType::Pointer *)pImg;
        *outimage = imageSource;


   } 

/* Convert a cartesian image to polar image */
template<typename TImage>
void CartesianToPolar(const typename TImage::Pointer sourceImage, void* destImage) {

    typedef TImage ImageType;
    typename ImageType::Pointer *outImage = (typename ImageType::Pointer *)destImage;

    // Resize
    typename ImageType::SpacingType outputSpacing;
    outputSpacing.Fill(1);

    typename ImageType::SizeType inputSize = sourceImage->GetLargestPossibleRegion().GetSize() ;

    typename ImageType::SizeType outputSize;
    outputSize[0] = inputSize[0] / 2;;
    outputSize[1] = 360;
    typename ImageType::PointType outputOrigin;
    outputOrigin.Fill(0);

    /** Transform typedef. */
    typedef itk::PolarToCartesianTransform< double, Dimension2D > PolarToCartesianTransformType;
  	PolarToCartesianTransformType::Pointer cylindricalTransform = PolarToCartesianTransformType::New();
    cylindricalTransform->SetForwardPolarToCartesian();
    
    /** Interpolator typedef. */
    typedef itk::LinearInterpolateImageFunction< ImageType, double > DefaultInterpolatorType;
    typedef typename DefaultInterpolatorType::Pointer InterpolatorPointer;
    InterpolatorPointer imInterpolator;
    imInterpolator = DefaultInterpolatorType::New();
    // Apply resample filter
    typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
    typename ResampleFilterType::Pointer filterResample = ResampleFilterType::New();
    filterResample->SetInput(sourceImage);
    filterResample->SetTransform(cylindricalTransform);
    filterResample->SetInterpolator(imInterpolator);
    filterResample->SetSize(outputSize);
    filterResample->SetOutputSpacing(outputSpacing);
    filterResample->SetOutputOrigin(outputOrigin);
    filterResample->SetDefaultPixelValue(0);
    filterResample->Update();
    sourceImage->DisconnectPipeline();
    *outImage = filterResample->GetOutput();

}


/* Convert a polar image to cartesian image */
template<typename TImage>
void PolarToCartesian(const typename TImage::Pointer sourceImage, void* destImage) {

    typedef TImage ImageType;
    typename ImageType::Pointer *outImage = (typename ImageType::Pointer *)destImage;

    // Resize
    typename ImageType::SpacingType outputSpacing;
    outputSpacing.Fill(1);
    typename ImageType::SizeType outputSize;
    typename ImageType::SizeType inputSize = sourceImage->GetLargestPossibleRegion().GetSize() ;
    sourceImage->Print(std::cout);
    outputSize[0] = 2*inputSize[0];
    outputSize[1] = 2*inputSize[0];
    typename ImageType::PointType outputOrigin;
    outputOrigin.Fill(0);
    outputOrigin[0] = -(double)(inputSize[0]);
    outputOrigin[1] = -(double)(inputSize[0]);

        /** Transform typedef. */
    typedef itk::PolarToCartesianTransform< double, Dimension2D > PolarToCartesianTransformType;
    PolarToCartesianTransformType::Pointer cylindricalTransform = PolarToCartesianTransformType::New();
    cylindricalTransform->SetForwardCartesianToPolar();    /** Interpolator typedef. */
    typedef itk::LinearInterpolateImageFunction< ImageType, double > DefaultInterpolatorType;
    typedef typename DefaultInterpolatorType::Pointer InterpolatorPointer;
    InterpolatorPointer imInterpolator;
    imInterpolator = DefaultInterpolatorType::New();
    // Apply resample filter
    typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
    typename ResampleFilterType::Pointer filterResample = ResampleFilterType::New();
    filterResample->SetInput(sourceImage);
    filterResample->SetTransform(cylindricalTransform);
    filterResample->SetInterpolator(imInterpolator);
    filterResample->SetSize(outputSize);
    filterResample->SetOutputSpacing(outputSpacing);
    filterResample->SetOutputOrigin(outputOrigin);
    filterResample->SetDefaultPixelValue(0);
    filterResample->Update();
    sourceImage->DisconnectPipeline();
    *outImage = filterResample->GetOutput();
}


int main(int, char*[])
{

  itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);

  AnatomicImageType::PointType origin;
  origin.Fill(0);

  double maxTheta = 360;
  double radius = 20;
  double xSpacing = 1.0;
  double ySpacing = 1.0;
  double thetaSpacing = 1.0;
  double radiusSpacing = 1.0;

  AnatomicImageType::Pointer polarImage = NULL;
  void* pImg = &polarImage;
  CreatePolarImage(origin,
                         maxTheta,
                         thetaSpacing, 
                         radius,
                         radiusSpacing, 
                         pImg);
  WriteImage<AnatomicImageType>("Polar.nii", polarImage);

  AnatomicImageType::Pointer polar2CartImage = NULL;
  pImg = &polar2CartImage;
  PolarToCartesian<AnatomicImageType>(polarImage, pImg);
  WriteImage<AnatomicImageType>("Polar2Cart.nii", polar2CartImage); 


  AnatomicImageType::Pointer cartImage = NULL;
  pImg = &cartImage;
  CreateCartesianImage( origin, 
                        maxTheta, 
                        radius,
                        xSpacing,                            
                        ySpacing,
                        pImg);
  WriteImage<AnatomicImageType>("Cart.nii", cartImage); 


  AnatomicImageType::Pointer cart2PolarImage = NULL;
  pImg = &cart2PolarImage;
  CartesianToPolar<AnatomicImageType>(cartImage, pImg);
  WriteImage<AnatomicImageType>("Cart2Polar.nii", cart2PolarImage); 


  return EXIT_SUCCESS;
}