#ifndef READINFEATURESELECTION_H_INCLUDED
#define READINFEATURESELECTION_H_INCLUDED

#include "configFlags.h"
#include "getNeighborhoodMatrices.h"


//void readInFeatureSelection(EFoobar::Flags &featureFlags, string featureSelectionPath);
void CalculateRelFeatures(Image<float, 3> imageAttr, ImageType *image, ImageType *mask, ConfigFile config);
void calculateRelFeaturesDiscretized(Image<float, 3> imageAttr, vector<double> spacing, ConfigFile config);
void writeLogFile(string logFileName, std::string &text);
#include "readInFeatureSelection.cpp"

#endif // readInFeatureSelection_H_INCLUDED
