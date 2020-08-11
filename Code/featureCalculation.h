#ifndef _FEATURECALCULATION_H_
#define _FEATURECALCULATION_H_
#include <boost/algorithm/string/replace.hpp>
#include <itkTestingComparisonImageFilter.h>
#include "readDicom.h"
#include "readPrj.h"
#include "image.h"
#include "softwareParameters.h"
#include "itkImageFileWriter.h"
#include <itkImageRegionConstIterator.h>
#include "itkBinaryThresholdImageFilter.h"
#include <boost/bind.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <iostream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "processing.h"
#include "dispersityFeatures.h"
/*! \file */


void readImageAndMask(ConfigFile config);
void calculateFeatures(ImageType *imageFiltered, ImageType *maskNewSpacing, ConfigFile config);
void writeImageData2Log(ConfigFile config);
ImageType::Pointer flipNII(ImageType::Pointer mask);
int getNrVoxels(ImageType *mask, float threshold);

#include "featureCalculation.cpp"
#endif