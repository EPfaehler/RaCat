#ifndef GETNEIGHBORHOODMATRICES_H_INCLUDED
#define GETNEIGHBORHOODMATRICES_H_INCLUDED

#include "matrixFunctions.h"
#include "image.h"
#include "NGLDMFeatures2DMRG.h"
#include "helpFunctions.h"
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;
template<typename T>
vector<T> getNeighborhood3D(boost::multi_array<T, 3> inputMatrix, int *indexOfElement, vector<double> spacing, ConfigFile config);
template<typename T>
void getNeighborhoodMatrix3D(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm3D, boost::multi_array<double, 2> &ngldm3D, vector<double> spacing, ConfigFile config);
template<typename T>
void getNeighborhoodMatrix3DNGTDM(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm3D, boost::multi_array<T, 3> nrNeighborMatrix, boost::multi_array<T, 3> sumMatrix, vector<double> spacing, ConfigFile config);
template<typename T>
void getNeighborhoodMatrix2DNGTDM(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm2D, boost::multi_array<T, 3> nrNeighborMatrix, boost::multi_array<T, 3> sumMatrix);
template<typename T>
void getNeighborhoodMatrix2D(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm2D, boost::multi_array<T, 3> &ngldm2D, vector<double> spacing, ConfigFile config); 
template<typename T>
vector<T> getNeighborhood(boost::multi_array<T, 3> inputMatrix, boost::multi_array<T, 3> &ngtdm2D, boost::multi_array<T, 3> &ngldm2D,int *indexOfElement, vector<double> spacing, ConfigFile config); 

#include "getNeighborhoodMatrices.cpp"
#endif
