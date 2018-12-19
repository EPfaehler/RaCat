#ifndef IMAGE_INCLUDED
#define IMAGE_INCLUDED


#include <iostream>
#include "math.h"

#include <string>
#include <algorithm>
#include <vector>

#include "boost/multi_array.hpp"
#include "itkMesh.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSimplexMesh.h"
#include "itkSimplexMeshVolumeCalculator.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"

#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"

#include "itkTypes.h"

#include <itkImageFileWriter.h>

#include "readConfigFile.h"
/*! \file */

/*!
The class image has as attributes the important properties of an image: \n
the image intensity values of the VOI are stored in a boost::multi_array \n
the different grey levels of the image are stored in one vector \n

Furthermore it contains, the functions to interpolate and discretize an image and store the discretized and interpolated
image values again in a boost::multi_array.

*/


//TODO
//here in the methods often i use the calculation of min und max
//I could calculate it once and use it
//I could save these values as an attribute of the image to avoid double calculations
using namespace std;

using namespace itkTypes;


template <class T, size_t R = 3>
class Image {

private:
	int nrBins;
	double binWidth;
	//factor to convert the image from kBq/ml to 
	double SUVfactor;

	int maxValueInMask;

public:
	//vector containing all different grey levels
	vector<T> diffGreyLevels;
	//create an image matrix
	boost::multi_array<T, R> imageMatrix;
	boost::multi_array<T, R> imageMatrixOriginal;
	boost::multi_array<T, R> imageMatrixLocalInt;
	boost::multi_array<T, R> imageMatrixIVH;
	//vector where every intensity value of a voxel inside the VOI is stored
	vector<T> vectorOfMatrixElements;
	//to discretize a resegmented image with a fixed number of bins, I need the original min and max value
	vector<T> vectorOfMatrixElementsOriginal;
	int nrRows;
	int nrCols;
	int nrDepth;
	T minGreyLevel;
	T maxGreyLevel;

	typename ImageType::Pointer image;
	typename ImageType::Pointer mask;


	//constructor
	Image(unsigned int row, unsigned int col, unsigned int depth);
	~Image();
	//interpolate the image mask
	void getInterpolatedImageMask(ConfigFile config, ImageType *image, ImageType *mask);
	
	//resample the image to the desired outputSpacing and outputSize
	ImageType::Pointer getResampledImage(ImageType *originalImage, double *outputSpacing, itk::Size<3> outputSize, string interpolationMethod);
	//assign the values of the image that are lying inside the mask to a boost multi array
	boost::multi_array<T, R> get3Dimage(ImageType *image, ImageType *mask, ConfigFile config);
	boost::multi_array<T, R> get3DimageResegmented(ImageType *image, ImageType *mask, ConfigFile config);
	boost::multi_array<T, R> get3DimageLocalInt(ImageType *image, ImageType *mask);
	//determine the different grey levels of the image
	vector<T> getGreyLevels();
	//save all values which are not NAN in one vector
	vector<T> getVectorOfMatrixElementsNotNAN(boost::multi_array<T, R> inputMatrix);
	//methods for discretization
	void discretizationFixedWidth(boost::multi_array<T, R> &inputMatrix, double intervalWidth, ConfigFile config);
	void discretizationFixedBinNr(boost::multi_array<T, R> &inputMatrix, vector<T> elementVector, double binNr);
	//convert the image from kBq/ml to SUV
	void calculateSUV(boost::multi_array<T, R> &inputMatrix, ConfigFile config);
	void calculateSUL(boost::multi_array<T, R> &inputMatrix, ConfigFile config);
	//assign values to the image attributes
	void getImageAttributes(ImageType *filteredImage, ImageType *mask, ConfigFile configName);
	void getImageAttributesDiscretized(ImageType *filteredImage, ImageType *maskFilter, ConfigFile configName);
	//with this function I determine, which value is inside the mask (100 or 1?)
	//that is needed to create the mesh correctly
	int getValueInMask(ImageType::Pointer mask);
};



//every time when an image object is created, the image matrix is also created
//and the config file is read to read in the image information
template<class T, size_t R>
Image<T, R>::Image(unsigned int row, unsigned int col, unsigned int depth) : imageMatrix(boost::extents[row][col][depth]), imageMatrixOriginal(boost::extents[row][col][depth]), imageMatrixLocalInt(boost::extents[row][col][depth]), imageMatrixIVH(boost::extents[row][col][depth]) {
	diffGreyLevels.clear();
	vectorOfMatrixElements.clear();
}

template<class T, size_t R>
Image<T, R>::~Image() {
}

/*!
\brief getInterpolatedImageMask

The image and the mask are interpolated using the nearest neighbor algorithm using up- or downsampling,
what was set by the user \n
In a later version, the user should be able to choose the interpolation method
*/

template<class T, size_t R>
void Image<T, R>::getInterpolatedImageMask(ConfigFile config, ImageType *imageFiltered, ImageType *maskFiltered) {

	const typename ImageType::SpacingType& inputSpacing = imageFiltered->GetSpacing();
	//downsampling
	double outputSpacing[3];
	// Fetch original image size
	const typename ImageType::RegionType& inputRegion = imageFiltered->GetLargestPossibleRegion();
	const typename ImageType::SizeType& inputSize = inputRegion.GetSize();
	unsigned int oldWidth = inputSize[0];
	unsigned int oldHeight = inputSize[1];
	unsigned int oldDepth = inputSize[2];
	unsigned int newWidth;
	unsigned int newHeight;
	unsigned int newDepth;
	if (config.useDownSampling == 1 ) {
		double minimum = inputSpacing[0];
		if (inputSpacing[1]<minimum) {
			minimum = inputSpacing[1];
		}
		if (inputSpacing[2] < minimum) {
			minimum = inputSpacing[2];
		}
		outputSpacing[0] = minimum;
		outputSpacing[1] = minimum;
		outputSpacing[2] = minimum;
		newWidth = (double)oldWidth * inputSpacing[0] / minimum;
		newHeight = (double)oldHeight * inputSpacing[1] / minimum;
		newDepth = (double)oldDepth * inputSpacing[2] / minimum;
	}
	else if (config.useUpSampling == 1) {
		int maximum = inputSpacing[0];
		if (inputSpacing[1]>maximum) {
			maximum = inputSpacing[1];
		}
		if (inputSpacing[2] > maximum) {
			maximum = inputSpacing[2];
		}
		outputSpacing[0] = maximum;
		outputSpacing[1] = maximum;
		outputSpacing[2] = maximum;
		newWidth = (double)oldWidth * inputSpacing[0] / maximum;
		newHeight = (double)oldHeight * inputSpacing[1] / maximum;
		newDepth = (double)oldDepth * inputSpacing[2] / maximum;
	}
	else{
		outputSpacing[0] = double(2);
		outputSpacing[1] = double(2);
		//outputSpacing[2] = inputSpacing[2];
		outputSpacing[2] = double(2);
		newWidth = ceil((double)oldWidth * inputSpacing[0] / double(2));
		newHeight = ceil((double)oldHeight * inputSpacing[1] / double(2));
		//newDepth = (double)oldDepth;
		newDepth = ceil((double)oldDepth * inputSpacing[2] / double(2));
	}
	
	if (config.imageType == "PET") {
		if (config.useReSegmentation == 1 ||config.excludeOutliers ==1) {
			imageMatrixOriginal = get3Dimage(image, mask, config);
			vectorOfMatrixElementsOriginal = getVectorOfMatrixElementsNotNAN(imageMatrixOriginal);
			imageMatrix = get3DimageResegmented(image, mask, config);
			imageMatrixLocalInt = get3DimageLocalInt(image, mask);

		}
		else {
			imageMatrix = get3Dimage(image, mask, config);
			imageMatrixLocalInt = get3DimageLocalInt(image, mask);

		}
	}
	else {
		if (config.useReSegmentation == 1 || config.excludeOutliers == 1) {
			imageMatrixOriginal = get3Dimage(image, mask, config);
			vectorOfMatrixElementsOriginal = getVectorOfMatrixElementsNotNAN(imageMatrixOriginal);
			imageMatrix = get3DimageResegmented(image, mask, config);
			imageMatrixLocalInt = get3DimageLocalInt(image, mask);
		}
		else {
			imageMatrix = get3Dimage(image, mask, config);
			imageMatrixLocalInt = get3DimageLocalInt(image, mask);
		}
	}


}


/*!
\brief getResampledImage
resample the image to the desired outputSize with desired spacing
@param[in] originalImage
@param[in] double outputSpacing: the spacing of the resampled image
@param[in] itk::size: outputSize: the size of the resampled image
*/
template<class T, size_t R>
ImageType::Pointer Image<T, R>::getResampledImage(ImageType *originalImage, double *outputSpacing, itk::Size<3> outputSize, string interpolationMethod) {
	//use trilinear interpolation
	typename LinearInterpolatorType::Pointer linearInterpolator = LinearInterpolatorType::New();
	typename ImageType::DirectionType direction;
	typename TransformType::Pointer transform = TransformType::New();
	typename InterpolatorType::Pointer interpolator;
	typename ResampleFilterType::Pointer resizeFilterImage = ResampleFilterType::New();
	typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType > GaussianFilterType;
	//if interpolation method is trilinear
	if (interpolationMethod == "Linear" || interpolationMethod == "linear") {
		interpolator= InterpolatorType::New();
		resizeFilterImage->SetInterpolator(interpolator);
	}
	//if interpolation method is spline
	if (interpolationMethod == "Spline" || interpolationMethod == "spline") {
		using splineInterpolatorType = itk::BSplineInterpolateImageFunction<ImageType, double >;
		splineInterpolatorType::Pointer interpolator = splineInterpolatorType::New();
		resizeFilterImage->SetInterpolator(interpolator);
	}
	if (interpolationMethod == "NearestNeighbor" || interpolationMethod == "nearestNeighbor") {
		using nnInterpolatorType = itk::NearestNeighborInterpolateImageFunction<ImageType, double >;
		nnInterpolatorType::Pointer interpolator = nnInterpolatorType::New();
		resizeFilterImage->SetInterpolator(interpolator);
	}
	
	
	resizeFilterImage->SetDefaultPixelValue(0);
	// Set the output spacing as specified on the command line
	resizeFilterImage->SetOutputSpacing(outputSpacing);
	//align voxel centers in interpolation step (and not the origin of the images) to be in compliance with IBSI
	ImageType::PointType newOrigin = originalImage->GetOrigin();
	const typename ImageType::SpacingType& inputSpacing = originalImage->GetSpacing();
	const typename ImageType::RegionType& inputRegion = originalImage->GetLargestPossibleRegion();
	const typename ImageType::SizeType& inputSize = inputRegion.GetSize();
	newOrigin[0] += 0.5*(inputSize[0] - 1)*inputSpacing[0] - 0.5*(outputSize[0] - 1)*outputSpacing[0];
	newOrigin[1] += 0.5*(inputSize[1] - 1)*inputSpacing[1] - 0.5*(outputSize[1] - 1)*outputSpacing[1];
	newOrigin[2] += 0.5*(inputSize[2] - 1)*inputSpacing[2] - 0.5*(outputSize[2] - 1)*outputSpacing[2];
	resizeFilterImage->SetOutputOrigin(newOrigin);
	resizeFilterImage->SetOutputDirection(originalImage->GetDirection());
	// Set the computed size
	resizeFilterImage->SetSize(outputSize);
	// Specify the input for the resamplersitk::Size<3> outputSize
	resizeFilterImage->SetInput(originalImage);
	resizeFilterImage->Update();

	ImageType::Pointer newImage = resizeFilterImage->GetOutput();
	return resizeFilterImage->GetOutput();
	
}



/*!
In the function getValueInMask the value inside the mask is detected. \n
It can vary from 1 - 100. The value inside the mask is needed in order to create a mesh from this mask.
@param[in]: ImageType::Pointer mask, the mask of the VOI
@param[out]: value inside mask (1-100)
*/
template <class T, size_t R>
int Image<T, R>::getValueInMask(ImageType::Pointer mask) {
	int maskValue = 0;
	int row = 0;
	int actPosition = 0;
	const typename ImageType::RegionType region = mask->GetBufferedRegion();
	const typename ImageType::SizeType imageSize = region.GetSize();
	while (actPosition < imageSize[0] * imageSize[1] * imageSize[2]) {
		if (mask->GetBufferPointer()[actPosition] > 0) {
			if (mask->GetBufferPointer()[actPosition] > maskValue) {
				maskValue = mask->GetBufferPointer()[actPosition];
			}
		}
		actPosition += 1;
	}
	return maskValue;
}

/*!
\brief get3Dimage
@param[in] ImageType image: original image
@param[in] ImageType mask: original mask
the image values are assigned to a boost::multi_array image matrix, if they are lying inside the mask
*/
template <class T, size_t R>
boost::multi_array<T, R> Image<T, R>::get3Dimage(ImageType *image, ImageType *mask, ConfigFile config) {
	//compare the size of image and mask and give an error if they are not the sam
	const typename ImageType::RegionType region = image->GetBufferedRegion();

	const typename ImageType::SizeType imageSize = region.GetSize();
	const typename ImageType::RegionType regionMask = mask->GetBufferedRegion();
	const typename ImageType::SizeType maskSize = regionMask.GetSize();

	if (maskSize[0] != imageSize[0] || maskSize[1] != imageSize[1] || maskSize[2] != imageSize[2]) {

		std::cout << "Image size and mask size are not corresponding!";
	}
	maxValueInMask = getValueInMask(mask);
	//store the values of the image in an multi dimensional array
	boost::multi_array<T, R> A(boost::extents[imageSize[0]][imageSize[1]][imageSize[2]]);
	int row = 0;
	int depth = 0;
	int col = 0;
	int actPosition = 0;
	int count = 0;
	T max = -1000000;
	T min= 1000000;
	while (actPosition < imageSize[0] * imageSize[1] * imageSize[2]) {
		if (row < imageSize[0] && depth < imageSize[2]) {
			if (mask->GetBufferPointer()[actPosition] >= config.threshold*maxValueInMask) {
				A[row][col][depth] = image->GetBufferPointer()[actPosition];
				count++;
				if (A[row][col][depth] > max) {
					max = A[row][col][depth];
				}
				if (A[row][col][depth] < min) {
					min = A[row][col][depth];
				}
				//if image is CT image, the image values have to be rounded to the closest integer after interpolation
				if (config.imageType == "CT") {
					if (config.useSampling2mm == 1 || config.useDownSampling == 1 || config.useUpSampling == 1) {
						A[row][col][depth] = round(A[row][col][depth]);	
					}
				}
			}
			else {
				A[row][col][depth] = NAN;
			}
			row = row + 1;
		}

		if (row == imageSize[0] && col < imageSize[1]) {
			col = col + 1;
			row = 0;
		}
		if (col == imageSize[1] && depth < imageSize[2]) {
			col = 0;
			depth = depth + 1;
		}
		actPosition += 1;
	}
	//std::cout << "nrVoxels"<<count<<" "<<max<<" "<<min << std::endl;
	return A;
}


/*!
\brief get3DimageLocalInt
@param[in] ImageType image: original image
@param[in] ImageType mask: original mask
the image values are assigned to a boost::multi_array image matrix, if they are lying in the bounding box. \n
For the calculation of the local intensity features, not only the image values inside the mask are required, but also the 
values lying outside the mask.
*/
template <class T, size_t R>
boost::multi_array<T, R> Image<T, R>::get3DimageLocalInt(ImageType *image, ImageType *mask) {
	//compare the size of image and mask and give an error if they are not the sam
	const typename ImageType::RegionType region = image->GetBufferedRegion();

	const typename ImageType::SizeType imageSize = region.GetSize();
	const typename ImageType::RegionType regionMask = mask->GetBufferedRegion();
	const typename ImageType::SizeType maskSize = regionMask.GetSize();

	if (maskSize[0] != imageSize[0] || maskSize[1] != imageSize[1] || maskSize[2] != imageSize[2]) {

		std::cout << "Image size and mask size are not corresponding!";
	}
	//store the values of the image in an multi dimensional array
	boost::multi_array<T, R> A(boost::extents[imageSize[0]][imageSize[1]][imageSize[2]]);
	int row = 0;
	int depth = 0;
	int col = 0;
	int actPosition = 0;
	while (actPosition < imageSize[0] * imageSize[1] * imageSize[2]) {
		if (row < imageSize[0] && depth < imageSize[2]) {
			A[row][col][depth] = image->GetBufferPointer()[actPosition];
			row = row + 1;
		}
		if (row == imageSize[0] && col < imageSize[1]) {
			col = col + 1;
			row = 0;
		}
		if (col == imageSize[1] && depth < imageSize[2]) {
			col = 0;
			depth = depth + 1;
		}
		actPosition += 1;
	}
	return A;
}


/*!
\brief get3DimageResegmented
@param[in] ImageType image: original image
@param[in] ImageType mask: original mask
the image values are assigned to a boost::multi_array image matrix, if they are lying inside the mask. Intensity values below or above the set 
resegmentation values are deleted from the mask. 
*/
template <class T, size_t R>
boost::multi_array<T, R> Image<T, R>::get3DimageResegmented(ImageType *image, ImageType *mask, ConfigFile config) {
	//compare the size of image and mask and give an error if they are not the sam
	const typename ImageType::RegionType region = image->GetBufferedRegion();

	const typename ImageType::SizeType imageSize = region.GetSize();
	const typename ImageType::RegionType regionMask = mask->GetBufferedRegion();
	const typename ImageType::SizeType maskSize = regionMask.GetSize();
	T max;
	T min;
	if (config.excludeOutliers == 1) {
		vector<T> vectorElementsThreshold;
		for (int i = 0; i < vectorOfMatrixElementsOriginal.size(); i++) {
			if (vectorOfMatrixElementsOriginal[i] >= config.minValueReSeg && vectorOfMatrixElementsOriginal[i] <= config.maxValueReSeg) {
				vectorElementsThreshold.push_back(vectorOfMatrixElementsOriginal[i]);
			}
		}
		
		//calculate the mean value of all elements
		T sum = std::accumulate(vectorElementsThreshold.begin(), vectorElementsThreshold.end(), 0.0);
		T mean = sum / vectorElementsThreshold.size();
		//calculate std
		transform(vectorElementsThreshold.begin(), vectorElementsThreshold.end(), vectorElementsThreshold.begin(), bind2nd(minus<T>(), mean));
		T sq_sum = std::inner_product(vectorElementsThreshold.begin(), vectorElementsThreshold.end(), vectorElementsThreshold.begin(), 0.0);;
		
		T stdev = std::sqrt(sq_sum / (vectorElementsThreshold.size()-1));
		max = mean + 3 * stdev;
		min = mean - 3 * stdev;  
		if (config.useReSegmentation == 1) {
			//apply both outlier methods in parallel
			if (max < config.maxValueReSeg) {
				config.maxValueReSeg = max;
			}
			if (min > config.minValueReSeg) {
				config.minValueReSeg = min;
			}
		}
		else {
			config.maxValueReSeg = max;
			config.minValueReSeg = min;
		}
	}
	if (maskSize[0] != imageSize[0] || maskSize[1] != imageSize[1] || maskSize[2] != imageSize[2]) {
		std::cout << "Image size and mask size are not corresponding!";
	}

	//store the values of the image in an multi dimensional array
	boost::multi_array<T, R> A(boost::extents[imageSize[0]][imageSize[1]][imageSize[2]]);
	int row = 0;
	int depth = 0;
	int col = 0;
	int actPosition = 0;
	int count = 0;
	while (actPosition < imageSize[0] * imageSize[1] * imageSize[2]) {
		if (row < imageSize[0] && depth < imageSize[2]) {
			//get resegmented mask
			if (mask->GetBufferPointer()[actPosition] >= config.threshold*maxValueInMask && image->GetBufferPointer()[actPosition] <= config.maxValueReSeg && image->GetBufferPointer()[actPosition] >= config.minValueReSeg) {
				A[row][col][depth] = image->GetBufferPointer()[actPosition];
				count++;
				//if image is CT image, the image values have to be rounded to the closest integer after interpolation
				if (config.imageType == "CT") {
					if (config.useSampling2mm == 1 || config.useDownSampling == 1 || config.useUpSampling == 1) {
						A[row][col][depth] = round(A[row][col][depth]);
					}
				}
			}
			else {
				A[row][col][depth] = NAN;
			}
			row = row + 1;
		}
		if (row == imageSize[0] && col < imageSize[1]) {
			col = col + 1;
			row = 0;
		}
		if (col == imageSize[1] && depth < imageSize[2]) {
			col = 0;
			depth = depth + 1;
		}
		actPosition += 1;
	}
	//std::cout <<"morphMask"<< count << std::endl;
	return A;
}


/*!
\brief getGreyLevels()
the fuction returns a vector, containing the different grey levels of the matrix
*/
template <class T, size_t R>
vector<T> Image<T, R>::getGreyLevels() {
	vector<T> temGreyLevels;
	//get Grey Levels which are not NAN
	vectorOfMatrixElements = getVectorOfMatrixElementsNotNAN(imageMatrix);
	std::copy(vectorOfMatrixElements.begin(), vectorOfMatrixElements.end(), back_inserter(temGreyLevels));
	//sort the vector and extract every element exactly once
	std::sort(temGreyLevels.begin(), temGreyLevels.end());
#ifdef _WIN32
	auto it = std::unique(temGreyLevels.begin(), temGreyLevels.end());
	temGreyLevels.resize(std::distance(temGreyLevels.begin(), it));
#else
	typename vector<T>::iterator it;
	it = std::unique(temGreyLevels.begin(), temGreyLevels.end());
	temGreyLevels.resize(std::distance(temGreyLevels.begin(), it));
#endif
	diffGreyLevels = temGreyLevels;
	return temGreyLevels;
}

/*!
\brief getVectorOfMatrixElementsNotNAN
In this function, all matrix elements, which are contained in the mask, are stored in a vector
*/
template<class T, size_t R>
vector<T> Image<T, R>::getVectorOfMatrixElementsNotNAN(boost::multi_array<T, R> inputMatrix) {
	vectorOfMatrixElements.clear();
	vector<T> elementVector;
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row < inputMatrix.shape()[0]; row++) {
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				if (!std::isnan(inputMatrix[row][col][depth])) {
					elementVector.push_back(inputMatrix[row][col][depth]);
				}
			}
		}
	}
	return elementVector;
}



/*!
\brief discretizationFixedWidth
@param inputMatrix
@param intervalWidth

The intensity values of the VOI are discretized using a fixed bin size \n
The bin size is given by the user in the .config file
*/

template <class T, size_t R>
void Image<T, R>::discretizationFixedWidth(boost::multi_array<T, R> &inputMatrix, double intervalWidth,ConfigFile config) {
	T minimumValue;
	T maximumValue;
	if (config.useReSegmentation == 1 && config.excludeOutliers==0) {
		minimumValue = config.minValueReSeg;
	}
	
	else {
		//get the minimum intensity value of the VOI
		minimumValue = *min_element(vectorOfMatrixElements.begin(), vectorOfMatrixElements.end());
	}
	
	//subtract the minimumValue from every pixel
	transform(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(),
		inputMatrix.origin(), std::bind2nd(minus<T>(), minimumValue));

	//divide the outcome by the interval width
	transform(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(),
		inputMatrix.origin(), std::bind2nd(std::divides<T>(), intervalWidth));
	//set the elements which are 0 to 1 (the elements which are now 0 were the minimum value before subtracting
	//replace(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), 0, 1);
	//round the values
	for (int i = 0; i<inputMatrix.shape()[0]; i++) {
		for (int j = 0; j<inputMatrix.shape()[1]; j++) {
			for (int k = 0; k<inputMatrix.shape()[2]; k++) {
				inputMatrix[i][j][k] = floor(inputMatrix[i][j][k])+1;
			}
		}
	}
}

/*!
\brief discretizationFixedBinNr
@param inputMatrix
@param binNr
Discretizes the intensity values inside the VOI to a fixed number of bins \n
The number of bins can be set by the user in the .config file
*/
template <class T, size_t R>
void Image<T, R>::discretizationFixedBinNr(boost::multi_array<T, R> &inputMatrix, vector<T> elementVector, double binNr) {
	float minimumValue = *min_element(elementVector.begin(), elementVector.end());
	//get the maximum intensity value
	float maximumValue = *max_element(elementVector.begin(), elementVector.end());
	if (minimumValue == 0 && minimumValue == maximumValue) {
		std::cout << "error in calculating discretization, VOI contains only 0" << std::endl;

		exit(0);
	}
	else if (minimumValue > 0 && minimumValue == maximumValue) {
		std::cout << "error in calculating discretization, VOI contains only one intensity value, minimum value is set to 0" << std::endl;
		minimumValue = 0;
	}
	//subtract minimum value from every matrix element
	transform(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(),
		inputMatrix.origin(), std::bind2nd(minus<T>(), minimumValue));
	//get the range
	
	float range = (maximumValue - minimumValue) / binNr;
	//divide every element of the matrix by the range
	transform(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(),
		inputMatrix.origin(), std::bind2nd(std::divides<T>(), range));
	replace(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), 0, 1);
	//do the rounding
	for (int i = 0; i<inputMatrix.shape()[0]; i++) {
		for (int j = 0; j<inputMatrix.shape()[1]; j++) {
			for (int k = 0; k<inputMatrix.shape()[2]; k++) {
				inputMatrix[i][j][k] = ceil(inputMatrix[i][j][k]) ;
				


				
			}
		}
	}
}

/*!
\brief calculateSUV
this function is only called if the image is a PET-image \n
it calculates the SUV-values from the original image using the information given by the user
*/
template<class T, size_t R>
void Image<T, R>::calculateSUV(boost::multi_array<T, R> &inputMatrix, ConfigFile config) {
	std::cout << "Convert intensity values to SUV" << std::endl;
	double correctionParam;
	if (config.correctionParam != 0) {
		std::cout << "The scaling factor is set. All image values will be multiplied by this parameter." << std::endl;
		correctionParam = config.correctionParam;
	}
	else {
		correctionParam = config.patientWeight / (config.initActivity * 1000);
	}
	for (int i = 0; i<inputMatrix.shape()[0]; i++) {
		for (int j = 0; j<inputMatrix.shape()[1]; j++) {
			for (int k = 0; k<inputMatrix.shape()[2]; k++) {
				if (!std::isnan(inputMatrix[i][j][k])) {
					inputMatrix[i][j][k] = inputMatrix[i][j][k] * correctionParam;
				}
			}
		}
	}
}



/*!
\brief calculateSUL
this function is only called if the image is a PET-image \n
it calculates the SUL-values from the original image using the information given by the user.
*/
template<class T, size_t R>
void Image<T, R>::calculateSUL(boost::multi_array<T, R> &inputMatrix, ConfigFile config) {
	std::cout << "Convert intensity values to SUL" << std::endl;
	double correctionParam;

	if (config.correctionParam != 0) {
		std::cout << "The scale factor parameter is set. All image values will be multiplied by this parameter." << std::endl;
		correctionParam = config.correctionParam;
	}
	else {
		if (config.malePatient == 1) {
			correctionParam = (1.1*config.patientWeight - 128 * (config.patientWeight / config.patientHeight)*(config.patientWeight / config.patientHeight)) / (config.initActivity * 1000);
		}
		else if (config.malePatient == 0) {
			correctionParam = (1.07*config.patientWeight - 148 * (config.patientWeight / config.patientHeight)*(config.patientWeight / config.patientHeight)) / (config.initActivity * 1000);
		}
	}
	for (int i = 0; i<inputMatrix.shape()[0]; i++) {
		for (int j = 0; j<inputMatrix.shape()[1]; j++) {
			for (int k = 0; k<inputMatrix.shape()[2]; k++) {
				inputMatrix[i][j][k] = inputMatrix[i][j][k] * correctionParam;
			}
		}
	}
}

/*!
\brief getImageAttributes
@param filteredImage
@param maskFilter
The function getImageAttributes takes filtered image and mask and stores the image information in a multi-dimensional array (using the get3Dimage-function) \n
It stores the elements which are inside the mask in an array - using the getVectorElementsNotNAN-function (in order to calculate the statistical features) \n
*/
template<class T, size_t R>
void Image<T, R>::getImageAttributes(ImageType *filteredImage, ImageType *maskFilter, ConfigFile config) {
	image = filteredImage;
	mask = maskFilter;
	//get the values inside the mask image
	maxValueInMask = getValueInMask(mask);
	//get image size
	const typename ImageType::RegionType region = filteredImage->GetBufferedRegion();
	const typename ImageType::SizeType imageSize = region.GetSize();
	nrRows = imageSize[0];
	nrCols = imageSize[1];
	nrDepth = imageSize[2];
	const typename ImageType::SpacingType& inputSpacing = filteredImage->GetSpacing();
	//std::cout << "spacing" << inputSpacing[0] << " " << inputSpacing[1] << " " << inputSpacing[2] << std::endl;
	
	if (config.useReSegmentation == 1 || config.excludeOutliers == 1) {
		imageMatrixOriginal = get3Dimage(image, mask, config);
		imageMatrix = get3DimageResegmented(image, mask, config);
		vectorOfMatrixElementsOriginal = getVectorOfMatrixElementsNotNAN(imageMatrixOriginal);
		imageMatrixLocalInt = get3DimageLocalInt(image, mask);
	}
	else {
		imageMatrix = get3Dimage(filteredImage, maskFilter, config);
		imageMatrixLocalInt = get3DimageLocalInt(filteredImage, maskFilter);
	}
	
	if (config.imageType == "PET") {
		if (config.useSUV == 1) {
			calculateSUV(imageMatrix, config);
			calculateSUV(imageMatrixLocalInt, config);

		}
		else if (config.useSUL == 1) {
			calculateSUL(imageMatrix, config);
			calculateSUL(imageMatrixLocalInt, config);
		}
	}
	vectorOfMatrixElements = getVectorOfMatrixElementsNotNAN(imageMatrix);
	diffGreyLevels = getGreyLevels();
}

/*!
\brief getImageAttributesDIscretized
@param filteredImage
@param maskFilter
The function getImageAttributesDiscretized is called before the textural features are calculated: \n
It discretizes the image matrix using the discretization method set by the user - using the discretization function \n
It stores the values of the discretized image in a matrix (using the get3Dimage-function) \n
It stores the elements of this matrix which are inside the mask in an array (using getVectorOfMatrixElementsNotNAN) \n
*/

template<class T, size_t R>
void Image<T, R>::getImageAttributesDiscretized(ImageType *filteredImage, ImageType *maskFilter, ConfigFile config) {
	image = filteredImage;
	mask = maskFilter;
	maxValueInMask = getValueInMask(mask);
	const typename ImageType::RegionType region = filteredImage->GetBufferedRegion();
	const typename ImageType::SizeType imageSize = region.GetSize();
	const typename ImageType::SpacingType& inputSpacing = filteredImage->GetSpacing();

	if (config.useReSegmentation == 1 || config.excludeOutliers ==1) {
		imageMatrixOriginal = get3Dimage(image, mask, config);
		vectorOfMatrixElementsOriginal = getVectorOfMatrixElementsNotNAN(imageMatrixOriginal);
		imageMatrix = get3DimageResegmented(image, mask, config);
	}
	else {
		imageMatrix = get3Dimage(filteredImage, maskFilter, config);
	}
	

	if (config.imageType == "PET") {
		if (config.useSUV == 1) {
			calculateSUV(imageMatrix, config);
		}
		else if (config.useSUL == 1) {
			calculateSUL(imageMatrix, config);

		}
	}
	vectorOfMatrixElements = getVectorOfMatrixElementsNotNAN(imageMatrix);
	if (config.useFixedBinWidthIVH == 1 && config.discretizeIVHSeparated == 1) {
		imageMatrixIVH = imageMatrix;
		discretizationFixedWidth(imageMatrixIVH, config.binWidthIVH, config);
	}
	//if the IVH values have to be calculated separately (can be set by the user)
	else if (config.useFixedNrBinsIVH == 1 && config.discretizeIVHSeparated == 1 && config.useReSegmentation == 1) {
		imageMatrixIVH = imageMatrix;
		discretizationFixedBinNr(imageMatrixIVH, vectorOfMatrixElementsOriginal, config.nrBinsIVH);
	}
	else if (config.useFixedNrBinsIVH == 1 && config.discretizeIVHSeparated == 1 && config.useReSegmentation == 0) {
		imageMatrixIVH = imageMatrix;
		discretizationFixedBinNr(imageMatrixIVH, vectorOfMatrixElements, config.nrBinsIVH);
	}
	if (config.useFixedBinWidth == 1) {
		discretizationFixedWidth(imageMatrix, config.binWidth, config);
	}

	if (config.useFixedNrBins == 1) {
		discretizationFixedBinNr(imageMatrix, vectorOfMatrixElements, config.nrBins);
	}
	diffGreyLevels = getGreyLevels();
}


#endif // IMAGE_INCLUDED