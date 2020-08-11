#ifndef _READDICOM_H_
#define _READDICOM_H_
#include "readImages.h"
#include "itkSpatialObject.h"

namespace gdcm { class Reader; }

ImageType::Pointer readDicom(string prjPath);
ImageType::Pointer resetImage(ImageType *mask);
void trim(std::string &str);
void reset2DImage(ImageSliceType::Pointer imageSlice);
void mergeImages(ImageSliceType::Pointer tempSlice, ImageType::Pointer finalImage, int iRequiredSlice);
void insertRegion(GroupType * group, ImageType::Pointer finalImage, int iSlice);
void writeFile(ImageType::Pointer mask, std::string outputFileName);
ImageType::Pointer readRTstruct(string voiPath, ImageType *imageDicom, ImageType *mask);
#include "readDicom.cpp"

#endif