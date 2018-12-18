#include "itkChangeInformationImageFilter.h"
#include <itkRegionOfInterestImageFilter.h>
#include "itkFlipImageFilter.h"
#include "itkFixedArray.h"
#include "itkAffineTransform.h"
#include "itkCenteredEuler3DTransform.h"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkRecursiveGaussianImageFilter.h"
#include <itkImageRegionIterator.h>
#include <itkNiftiImageIO.h>

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"

#include "itkSimpleFilterWatcher.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkAffineTransform.h"

#include <itkExtractImageFilter.h>
#include "itkImageMaskSpatialObject.h"
#include "itkCastImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkImportImageFilter.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"

#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"

#include <itkPoint.h>
#include <itkImageSliceIteratorWithIndex.h>
#include <itkImageIteratorWithIndex.h>
#include <itkPointSetToImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkPolygonSpatialObject.h>
#include <itkImageFileWriter.h>
#include <gdcmTypes.h>  
#include <gdcmReader.h>
#include <gdcmSmartPointer.h>
#include <gdcmAttribute.h>
#include <gdcmSequenceOfItems.h>
#include <gdcmFileMetaInformation.h>

namespace itkTypes {
	typedef itk::Image< float, 3 > ImageType;
	typedef itk::Image< float, 2 > ImageType2D;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImportImageFilter<float, 3> ImportFilterType;
	typedef itk::ImportImageFilter<float, 2> ImportFilterType2D;
	typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
	typedef itk::ImageMaskSpatialObject<3>      ImageMaskSpatialObject;
	typedef ImageMaskSpatialObject::ImageType   ImageTypeMask;
	typedef ImageTypeMask::RegionType               RegionType;
	typedef itk::ImageFileReader< ImageTypeMask >   ReaderTypeMask;
	typedef itk::Image<unsigned char, 3> charImage;
	typedef itk::Image<int, 3> intImage;
	typedef itk::CastImageFilter<ImageType, intImage> CastFilterType;
	typedef itk::CastImageFilter<ImageType, charImage> CastFilterTypeChar;
	typedef ImageMaskSpatialObject::ImageType ImageTypeSpatial;
	typedef ImageTypeSpatial::RegionType RegionType;
	//adjust spacing of mask to image
	typedef itk::ChangeInformationImageFilter< ImageType > ChangeInfoFilter;
	typedef itk::LinearInterpolateImageFunction<ImageType, float> LinearInterpolatorType;

	typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
	typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> T_Interpolator;
	typedef itk::ResampleImageFilter<ImageType, ImageType>    ResampleFilterType;
	typedef itk::AffineTransform <double, 3> TransformType;
	typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;

	//typedef to read RTstruct
	typedef itk::Point< double, 3 > PointType;

	typedef itk::PolygonSpatialObject<2> PolygonType;
	typedef PolygonType::Pointer PolygonPointer;
	typedef itk::SpatialObjectPoint<2> PolygonPointType;

	typedef itk::Image< float, 2 >   ImageSliceType;

	typedef itk::GroupSpatialObject<2> GroupType;
	typedef itk::SpatialObjectToImageFilter<GroupType, ImageSliceType> SpatialObjectToImageFilterType;
	typedef itk::ImageFileWriter< ImageType > WriterType;
	typedef itk::ExtractImageFilter< ImageType, ImageSliceType> FilterTypeRTS;

}
