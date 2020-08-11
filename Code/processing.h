#ifndef _PROCESSING_H_
#define _PROCESSING_H_
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
#include "intensityHistogram.h"
#include "morphologicalFeatures.h"
#include "localIntensityFeatures.h"

#include "GLCMFeatures2DVMRG.h"
#include "GLCMFeatures2DMRG.h"
#include "GLCMFeatures2DAVG.h"
#include "GLCMFeatures2DDMRG.h"
#include "GLCMFeatures3DMRG.h"
#include "GLCMFeatures3DAVG.h"

#include "GLRLMFeatures2DVMRG.h"
#include "GLRLMFeatures2DDMRG.h"
#include "GLRLMFeatures2DMRG.h"
#include "GLRLMFeatures2DAVG.h"
#include "GLRLMFeatures3D.h"
#include "GLRLMFeatures3DAVG.h"

#include "GLSZMFeatures2D.h"
#include "GLSZMFeatures2DAVG.h"
#include "GLSZMFeatures3D.h"

#include "NGLDMFeatures2DMRG.h"
#include "NGLDMFeatures2DAVG.h"
#include "NGLDMFeatures3D.h"

#include "NGTDM2DMRG.h"
#include "NGTDM2DAVG.h"
#include "NGTDM3D.h"

#include "GLDZMFeatures2DMRG.h"
#include "GLDZMFeatures3D.h"
#include "GLDZMFeatures2DAVG.h"

#include "intensityVolumeFeatures.h"
#include "localIntensityFeatures.h"
#include "dispersityFeatures.h"
#include "string"

void fillCSVwithNANs(ConfigFile config);
ImageType::Pointer thresholdMask(ImageType *mask, float threshold);
void writePETmetrics(float value, string nameVariable, ConfigFile config);
void writeExactVolume(float volume, ConfigFile config);
void calculatePETmetrics(ImageType::Pointer image, ImageType::Pointer mask, int volume, ConfigFile config);
void storePreInterpolationFeatures(ImageType::Pointer image, ImageType::Pointer mask, ConfigFile config);
float getOriginalVolume(ImageType::Pointer mask);
#include "processing.cpp"
#endif