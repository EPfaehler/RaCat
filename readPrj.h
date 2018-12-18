#ifndef _READPRJ_H_
#define _READPRJ_H_

#include "readImages.h"
#include "featureCalculation.h"
#include "itkTypes.h"
//functions to read accurate files
ImageType::Pointer converArray2Image(float *imageArray, unsigned int* dim, float *voxelSize);
//read the file into array and stores it as ITK image
ImageType::Pointer readPrjFilePET(string prjPath, string imageType, float smoothingKernel, unsigned int(&dim)[3], float (&voxelSize)[3]);
//get image dimensions etc from the prj file
void getImageDimension(ifstream &inFile, unsigned int(&dim)[3]);
void getVoxelSize(ifstream &inFile, float(&voxelSize)[3]);
void getImageValues(ifstream &inFile, vector<float> &imageArray);
void getHalfLifeScale(ifstream &inFile, float &halfLife, float &volumeScale, float &frameScale);
#include "readPrj.cpp"

#endif