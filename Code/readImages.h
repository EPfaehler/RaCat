#ifndef _READIMAGES_H_
#define _READIMAGES_H_
#include "itkGDCMImageIO.h"
#include "itkMetaDataObject.h"
#include "readInFeatureSelection.h"
#include "morphologicalFeatures.h"
#include "itkTypes.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <bitset>
#include <sstream>

using namespace itkTypes;
/*! \file */
//    define the types we will use later for reading the image and getting a specific region

ImageType::Pointer readImage(string imageName);
RegionType getBoundingBoxMask(ImageType *mask);
ImageType::Pointer getImageMasked(ImageType *image, RegionType boundingBoxRegion);
ImageType::Pointer getMaskNewSpacing(ImageType *imageFiltered, ImageType *maskFiltered);
ImageType::Pointer smoothImage(ImageType *image, float kernel);
vector<int> getImageSizeInterpolated(ImageType *imageFiltered, ImageType::SizeType imageSize, double (&outputSpacing)[3], ConfigFile config);
//change values inside mask to 1
ImageType::Pointer maskValues2One(ImageType *originalMask);


#include "readImages.cpp"

#endif