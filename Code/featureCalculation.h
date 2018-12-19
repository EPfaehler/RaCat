#ifndef _FEATURECALCULATION_H_
#define _FEATURECALCULATION_H_

#include "readDicom.h"
#include "readPrj.h"
#include "image.h"
#include "softwareParameters.h"

/*! \file */



void prepareDataForFeatureCalculation(ConfigFile config);
void calculateFeaturesForConfig(ImageType *imageFiltered, ImageType *maskNewSpacing, ConfigFile config);
ImageType::Pointer readVoiFilePET(string prjPath, string voiPath, ImageType *image, ConfigFile config, unsigned int(&dim)[3], float voxelSize[3]);
void writeImageData2Log(ConfigFile config);

#include "featureCalculation.cpp"
#endif