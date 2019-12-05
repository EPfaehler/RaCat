#ifndef GETNEIGHBORHOODMATRICES_H_INCLUDED
#define GETNEIGHBORHOODMATRICES_H_INCLUDED

#include "matrixFunctions.h"
#include "itkMaskImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkImportImageFilter.h"
#include "itkCropImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "image.h"
#include "NGLDMFeatures2DMRG.h"
#include "itkConvolutionImageFilter.h"
#include "helpFunctions.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include "amp.h"
//#include "readPrj.h"
using namespace Concurrency;
using namespace std;

template<typename T>
vector<T> getNeighborhood3D(boost::multi_array<T, 3> inputMatrix, int *indexOfElement, vector<float> spacing, ConfigFile config);
template<typename T>
void getNeighborhoodMatrix3D(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm3D, boost::multi_array<float, 2> &ngldm3D, vector<double> spacing, ConfigFile config);
template<typename T>
void getNeighborhoodMatrix3DNGTDM(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm3D, boost::multi_array<T, 3> nrNeighborMatrix, boost::multi_array<T, 3> sumMatrix, vector<double> spacing, ConfigFile config);
template<typename T>
void getNeighborhoodMatrix2DNGTDM(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm2D, boost::multi_array<T, 3> nrNeighborMatrix, boost::multi_array<T, 3> sumMatrix);
template<typename T>
void getNeighborhoodMatrix2D(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm2D, vector<double> spacing, ConfigFile config); 
template<typename T>
void fillConvolutionalVector(vector<T> &convolutionalVector, int is3D);
template<typename T>
ImageType::Pointer constructEmptyNewImage(Image<T, 3> imageAttr, boost::multi_array<T, 3> actualMatrix, int is3D);
template<typename T>
vector<T> getNeighborhood(boost::multi_array<T, 3> inputMatrix, boost::multi_array<T, 3> &ngtdm2D, boost::multi_array<T, 3> &ngldm2D,int *indexOfElement, vector<double> spacing, ConfigFile config); 
template<typename T>
ImageType::Pointer convolutionImage(ImageType::Pointer image, ImageType::Pointer kernel);
template<typename T>
ImageType::Pointer generateImageForConvolutionNGTDM(Image<T, 3> imageAttr);
template<typename T>
ImageType::Pointer generateMaskForConvolutionNGTDM(Image<T, 3> imageAttr, ConfigFile config);
template<typename T>
void getNeighborhoodMatrix2DNGLDM(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngldm2D, vector<float> spacing, ConfigFile config);
template<typename T>
void getNeighborhoodMatrix3D_convolution(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm3D, vector<float> spacing, ConfigFile config);
template<typename T>
void getNGLDMatrix3D_convolution(Image<T, 3> imageAttr,  boost::multi_array<float, 2> &ngldm3D, vector<float> spacing, ConfigFile config);
template<typename T>
ImageType::Pointer constructNewImage(ImageType::Pointer actImage, boost::multi_array<T, 3> actualMatrix);
template<typename T>
vector<T> getNeighborhood3D_convolution(boost::multi_array<T, 3> inputMatrix, int *indexOfElement, vector<float> spacing, ConfigFile config);
#include "getNeighborhoodMatrices.cpp"
#endif
