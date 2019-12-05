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

#include "string"
/*! \file */


ImageType::Pointer thresholdMask(ImageType *mask, float threshold);
void prepareDataForFeatureCalculation(ConfigFile config);
void calculateFeaturesForConfig(ImageType *imageFiltered, ImageType *maskNewSpacing, ConfigFile config);
ImageType::Pointer readVoiFilePET(string prjPath, string voiPath, ImageType *image, ConfigFile config, unsigned int(&dim)[3], float voxelSize[3]);
void writeImageData2Log(ConfigFile config);
void writeExactVolume(float volume, ConfigFile config);
void calculatePETmetrics(ImageType::Pointer image, ImageType::Pointer mask, int volume, ConfigFile config);
void writePETmetrics(float value, string nameVariable, ConfigFile config);
float getOriginalVolume(ImageType::Pointer mask);
#include "featureCalculation.cpp"
#endif