#ifndef MORPHOLOGICALFEATURES_H_INCLUDED
#define MORPHOLOGICALFEATURES_H_INCLUDED

/*! \file */


#include <cmath>
#include "matrixFunctions.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#ifdef _WIN32
#else
#include "itkVTKPolyDataReader.h"
#endif
#include "itkMesh.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkSimplexMesh.h"
#include "itkSimplexMeshVolumeCalculator.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"

#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkConnectedComponentImageFilter.h"

#include "itkImageMaskSpatialObject.h"
#include <boost/math/special_functions/legendre.hpp>
#include "itkCastImageFilter.h"

#include "image.h"
#include "itkTypes.h"

using namespace itkTypes;
template <class T, size_t R>
class MorphologicalFeatures {

private:

	ImageType::Pointer image;
	ImageType::Pointer mask;
	typedef itk::DefaultStaticMeshTraits<float, 3, 3, float, float, float> MeshTrait;

	//in order to calculate the surface we have to convert our mask in a mesh
	typedef itk::Mesh<float, 3, MeshTrait> MeshType;
	typedef itk::BinaryMask3DMeshSource<ImageType, MeshType> MeshSourceType;
	typedef unsigned char PixelType;
	//to calculate the surface we have to iterate over all cells of the mesh
	typedef MeshType::CellType CellType;
	typedef MeshType::CellsContainer::ConstIterator Celliterator;

#ifdef _WIN32
#else
	//to get the coordinates of the points of our mesh, we need a VTKPolyDataReader
	typedef itk::VTKPolyDataReader< MeshType > ReaderTypeVTK;
	typedef ReaderTypeVTK::PointType PointType;
#endif
	//in order to calculate volume and flatness etc we need a shape label object
	typedef itk::ConnectedComponentImageFilter <intImage, intImage> ConnectedComponentFilterType;

	typedef itk::LabelImageToShapeLabelMapFilter<intImage> LabelImageToShapeLabelMapFilterType;
	typedef int LabelType;
	typedef itk::ShapeLabelObject<LabelType, R> ShapeLabelObjectType;
	typedef itk::LabelMap<ShapeLabelObjectType> LabelMapType;
	//the attributes of the class are the radiomics features
	float volume;
	float appVolume;
	float surface;
	float surface2volumeRatio;
	float compactness1;
	float compactness2;
	float sphericalDisproportion;
	float sphericity;
	float asphericity;
	float majorAxisLength;
	float minorAxisLength;
	float leastAxisLength;
	float maximumDiameter;
	float flatness;
	float elongation;
	T integratredInt;
	T centerOfMassShift;
	T meanValue;
	T absMean;
	float moransI;
	float gearysC;

	int nrPixels;

	float volDensityAEE;
	float areaDensityAEE;

	float volDensityOMBB;
	float areaDensityOMBB;

	float volDensityAABB;
	float areaDensityAABB;

	float volDensityMEE;
	float areaDensityMEE;

	vector<vector<float> > vec;
	vector<float> coordinates;
	itk::Vector< double, 3> principalMoments;
	typename ImageType::PointType origin;

	double imageSpacingX;
	double imageSpacingY;
	double imageSpacingZ;

	//two help functions that are needed in order to calculate some of the features
	//calculate the euclidean distance
	float calculateEuclideanDistance(int rowDiff, int colDiff, int depthDiff);
	//convert the rows, cols and depths to cartesian coordinates
	vector<float>  convertToCoordinates(int row, int col, int depth);
	//this function subsamples the image 
	//for big tumors this has to be done before MoransI and GearysC can be calculated
	ImageType::Pointer subsampleImage(ImageType::Pointer image, int factor);

	//get all features that can be calculated by using the itk::label object
	void getLabelObjectFeatures(ImageType::Pointer mask);

	void getBoundingBoxValues(ImageType::Pointer mask);
	void calculateVADensity(float &volDensity, float &areaDensity, itk::Size<R> regionSize);
	void calculateApproximateVolume(boost::multi_array<T, R> inputMatrix, vector<T> vectorOfMatrElements);
	void calculateSurface2Volume();
	void calculateCompactness1();
	void calculateCompactness2();
	void calculateSphericalDisproportion();
	void calculateSphericity();
	void calculateAsphericity();

	void calculateMajorAxisLength();
	void calculateMinorAxisLength();
	void calculateLeastAxisLength();

	void calculateElongation();
	void calculateFlatness();
	//calculate MoransI calculates MoransI, GearysC and the maximum 3D diameter (because all these features need to
	//do two for loops)
	void calculateMoransI(boost::multi_array<T, R> inputMatrix);
	void calculateIntegratedIntensity(vector<T> vectorOfMatrixElements);
	void calculateCentreOfMassShift(boost::multi_array<T, R> inputMatrix, vector<T> vectorOfMatrElements);
	int binomialCoefficient(int n, int k);
	float legendrePolynom(float x, int exponent);
	void calculateVolDensityAEE();
	void calculateAreaDensityAEE();

	void calculateVolDensityMEE();
	void calculateAreaDensityMEE();

	void defineMorphologicalFeatures(vector<string> &features);
	void extractMorphologicalData(vector<T> &morphData, MorphologicalFeatures<T, R> morphFeatures);

public:
	//constructor
	MorphologicalFeatures() {
	}
	//destructor
	~MorphologicalFeatures() {
	}
	vector<T> cartCoordinateROI;
	vector<T> voxVolumeROI;
	const double pi = 3.141592653589793238463;
	boost::multi_array<vector<T>, R> coordinatesMatrix;
	void calculateAllMorphologicalFeatures(MorphologicalFeatures<T, R> &morphFeatures, Image<float, 3> imageAttr, ConfigFile config);
	void writeCSVFileMorphological(MorphologicalFeatures<T, R> morph, string outputFolder);
	void writeOneFileMorphological(MorphologicalFeatures<T, R> morph, string outputFolder);


};

/*!
In the function calculateEuclideanDistance the euclidean distance is calculated. \n
As input the function gets differences of row, columns and depth values and the euclidean distance is calculated
using this values.\n
@param[in]: int rowDiff, int colDiff, int depthDiff
@param[out]: float dist, Euclidean distance of voxels
*/

template <class T, size_t R>
float  MorphologicalFeatures<T, R>::calculateEuclideanDistance(int rowDiff, int colDiff, int depthDiff) {
	float dist;
	dist = sqrt(pow(rowDiff * imageSpacingX, 2) + pow(colDiff*imageSpacingY, 2) + pow(depthDiff*imageSpacingZ, 2));
	//dist = sqrt(pow(rowDiff, 2) + pow(colDiff, 2) + pow(depthDiff, 2));
	return dist;
}





/*!
In the function convertToCoordinates we convert the row, col and depth value of a voxel in a matrix to image coordinates
@param[in] int row, int col, int depth: position of voxel in matrix
@param[out] vector<float> coordinates: image coordinates of these points
*/
template <class T, size_t R>
vector<float>  MorphologicalFeatures<T, R>::convertToCoordinates(int row, int col, int depth) {
	float xcoord;
	float ycoord;
	float zcoord;
	origin = mask->GetOrigin();
	vector<float> newCoordinates;
	xcoord = origin[0] + row*imageSpacingX + double(imageSpacingX) / 2;
	ycoord = origin[1] + col*imageSpacingY + double(imageSpacingY) / 2;
	zcoord = origin[2] + depth * imageSpacingZ + double(imageSpacingZ) / 2;
	newCoordinates.push_back(xcoord);
	newCoordinates.push_back(ycoord);
	newCoordinates.push_back(zcoord);
	return newCoordinates;
}

/*!
In order ot accelerate the calculations of MoransI and GearysC for big VOIs, the image has to be downsampled
@param[in] ImageType::Pointer image
@param[in] int factor: factor for downsampling
@param[out] ImageType::Pointer downsampled image
*/
template<class T, size_t R>
ImageType::Pointer MorphologicalFeatures<T, R>::subsampleImage(ImageType::Pointer imageR, int factor) {
	typedef itk::Image< T, R > ImageType;
	const typename ImageType::SpacingType& inputSpacing = imageR->GetSpacing();
	const typename ImageType::RegionType& inputRegion = imageR->GetLargestPossibleRegion();
	const typename ImageType::SizeType& inputSize = inputRegion.GetSize();
	unsigned int oldWidth = inputSize[0];
	unsigned int oldHeight = inputSize[1];
	unsigned int oldDepth = inputSize[2];
	unsigned int newWidth;
	unsigned int newHeight;
	unsigned int newDepth;
	newWidth = (double)oldWidth * inputSpacing[0] / (double)factor;
	newHeight = (double)oldHeight * inputSpacing[1] / (double)factor;
	newDepth = (double)oldDepth * inputSpacing[2] / (double)factor;
	double outputSpacing[3];
	outputSpacing[0] = inputSpacing[0] * factor;
	outputSpacing[1] = inputSpacing[1] * factor;
	outputSpacing[2] = inputSpacing[2] * factor;
	imageSpacingX = outputSpacing[0];
	imageSpacingY = outputSpacing[1];
	imageSpacingZ = outputSpacing[2];
	itk::Size<3> outputSize = { { newWidth, newHeight, newDepth } };
	Image<float, 3> imageMask(outputSize[0], outputSize[1], outputSize[2]);
	image = imageMask.getResampledImage(imageR, outputSpacing, outputSize, "Spline");
	return image;
}

/*!
In the function getLabelObjectFeatures a set of morphological features is calculated using
the label image to shape label map filter
*/
template<class T, size_t R>
void MorphologicalFeatures<T, R>::getLabelObjectFeatures(ImageType::Pointer mask) {
	
	//With the cast filter type function, the mask is converted to an int-image
	//because otherwise the Label image to shape label map filter is not working
	typename CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(mask);
	//With the ConnectedComponentFilter from the itk library the connected components of the mask are calculated
	typename ConnectedComponentFilterType::Pointer connectedComponentImageFilter = ConnectedComponentFilterType::New();
	connectedComponentImageFilter->SetInput(castFilter->GetOutput());
	connectedComponentImageFilter->Update();
	//With the label image to shape label map filter the mask is converted to a labeled image
	typename LabelImageToShapeLabelMapFilterType::Pointer labelImageToShapeLabelMapFilter = LabelImageToShapeLabelMapFilterType::New();
	labelImageToShapeLabelMapFilter->SetComputeOrientedBoundingBox(true);
	labelImageToShapeLabelMapFilter->SetComputeFeretDiameter(true);

	labelImageToShapeLabelMapFilter->SetInput(connectedComponentImageFilter->GetOutput());
	labelImageToShapeLabelMapFilter->SetComputePerimeter(true);
	labelImageToShapeLabelMapFilter->Update();

	LabelMapType *labelMap = labelImageToShapeLabelMapFilter->GetOutput();
	//For every connected component a labelObject is created
	//Because we are only interested in the object with the label one, we check the label number
	volume = 0;
	surface = 0;
	nrPixels = 0;
	maximumDiameter = 0;
	for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); n++) {
		ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
		int labelNr = labelObject->GetLabel();
		if (labelNr >0) {
#ifdef _WIN32
			surface += labelObject->GetPerimeter();
			volume += labelObject->GetPhysicalSize();
			maximumDiameter += labelObject->GetFeretDiameter();
			if (labelNr == 1) {
				principalMoments = labelObject->GetPrincipalMoments();
			}
			else {
				principalMoments += labelObject->GetPrincipalMoments();
			}
			nrPixels += labelObject->GetNumberOfPixels();
			const itk::Vector<double, R> bbSize = labelObject->GetOrientedBoundingBoxSize();
			itk::Size<R> bbSizeITK = { { bbSize[0], bbSize[1], bbSize[2] } };
			calculateVADensity(volDensityOMBB, areaDensityOMBB, bbSizeITK);
#else
			volume += labelObject->GetPhysicalSize();

			principalMoments = labelObject->GetBinaryPrincipalMoments();
			nrPixels += labelObject->GetSize();
#endif

		}
	}
}


/*!
In the function getBoundingBoxValues the bounding box region is extracted from the image\n
From this region, volume and surface are extracted in order to calculate volume and area density.
*/
template<class T, size_t R>
void MorphologicalFeatures<T, R>::getBoundingBoxValues(ImageType::Pointer mask) {
	typedef itk::ImageMaskSpatialObject<R> ImageMaskSpatialObject;
	typedef ImageMaskSpatialObject::ImageType ImageTypeSpatial;
	typedef ImageTypeSpatial::RegionType RegionType;
	ImageMaskSpatialObject::Pointer maskSO = ImageMaskSpatialObject::New();
	CastFilterTypeChar::Pointer castFilterChar = CastFilterTypeChar::New();
	castFilterChar->SetInput(mask);
	charImage *castImage = castFilterChar->GetOutput();
	castImage->Update();
	maskSO->SetImage(castImage);
	maskSO->Update();
	
	RegionType boundingBoxRegion = maskSO->GetAxisAlignedBoundingBoxRegion();
	itk::Size<R> bbSize = boundingBoxRegion.GetSize();
	calculateVADensity(volDensityAABB, areaDensityAABB, bbSize);
}


/*!
In the function calculateVADensity the volume and the area density of a bounding box are calculated, given the size of the bounding box.
*/
template<class T, size_t R>
void MorphologicalFeatures<T, R>::calculateVADensity(float &volDensity, float &areaDensity, itk::Size<R> regionSize) {
	volDensity = regionSize[0]* imageSpacingX* regionSize[1] * imageSpacingY * regionSize[2] * imageSpacingZ;
	volDensity = volume / volDensity;
	areaDensity = 2*(regionSize[0] * imageSpacingX * regionSize[1] * imageSpacingY + regionSize[1] * imageSpacingY * regionSize[2] * imageSpacingZ + regionSize[0] * imageSpacingX * regionSize[2] * imageSpacingZ);
	areaDensity = surface/ areaDensity;
}

/*!
To calculate the integrated intensity, we calculate the mean value of the VOI and mulitplicate
the value by the volume
*/
template<class T, size_t R>
void MorphologicalFeatures<T, R>::calculateIntegratedIntensity(vector<T> vectorOfMatrixElements) {
	typedef boost::accumulators::features <tag::mean> Features;
	typedef accumulator_set <T, Features> Accumulator;
	Accumulator acc;
	for_each(vectorOfMatrixElements.begin(), vectorOfMatrixElements.end(), boost::bind<void>(boost::ref(acc), _1));
	meanValue = mean(acc);
	integratredInt = meanValue*volume;
}

/*!
To calculate the volume, we count the voxels present in the image
The amount of voxels is multiplied by the volume of one voxel
*/
template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateApproximateVolume(boost::multi_array<T, R> inputMatrix, vector<T> vectorOfMatrElement) {
	int nrVoxels = 0;
	for (int i = 0; i < inputMatrix.shape()[0]; i++) {
		for (int j = 0; j < inputMatrix.shape()[1]; j++) {
			for (int k = 0; k < inputMatrix.shape()[2]; k++) {
				if (!isnan(inputMatrix[i][j][k])) {
					nrVoxels += 1;
				}
			}
		}
	}
	appVolume = nrVoxels*imageSpacingX*imageSpacingY*imageSpacingZ;
}



template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateSurface2Volume() {
	surface2volumeRatio = surface / appVolume;
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateCompactness1() {
	compactness1 = appVolume / (sqrt(pi)*pow(surface, 1.5));
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateCompactness2() {
	compactness2 = 36 * pi*pow(appVolume, 2) / pow(surface, 3);
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateSphericalDisproportion() {
	sphericalDisproportion = pow(compactness2, 0.333333);
	sphericalDisproportion = 1 / sphericalDisproportion;
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateSphericity() {
	sphericity = pow(compactness2, 0.333333);
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateAsphericity() {
	asphericity = pow(1 / (36 * pi)*(pow(surface, 3) / pow(appVolume, 2)), 0.33333) - 1;
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateMajorAxisLength() {
	majorAxisLength = 4 * sqrt(principalMoments[2]);
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateMinorAxisLength() {
	minorAxisLength = 4 * sqrt(principalMoments[1]);
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateLeastAxisLength() {
	leastAxisLength = 4 * sqrt(principalMoments[0]);
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateElongation() {
	elongation = sqrt(minorAxisLength / majorAxisLength);
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateFlatness() {
	flatness = sqrt(leastAxisLength / majorAxisLength);
}


/*!
In the function calculateCentreOfMassShift, the center of mass shift is calculated, taking
as input the original matrix and a vector containing the matrix elements
*/
template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateCentreOfMassShift(boost::multi_array<T, R> inputMatrix, vector<T> vectorOfMatrElements) {
	//calculate the sum of all matrix elements
	T sumMatrElement = accumulate(vectorOfMatrElements.begin(), vectorOfMatrElements.end(), 0);
	vector<float> tempVector;
	//total number of voxels of the VOI
	int nrVoxels = boost::size(vectorOfMatrElements);
	double centerX = 0;
	double centerY = 0;
	double centerZ = 0;
	double x = 0;
	double y = 0;
	double z = 0;
	//calculate for every matrix element the center of the element in cartesian coordinates
	//as well as the grey level centre
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row < inputMatrix.shape()[0]; row++) {
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				if (!std::isnan(inputMatrix[row][col][depth])) {
					//get the cartesian coordinates of the center of the actual voxel
					tempVector = convertToCoordinates(row, col, depth);
					//get the grey level multiplied by the center
					centerX += inputMatrix[row][col][depth] * tempVector[0] / sumMatrElement;
					centerY += inputMatrix[row][col][depth] * tempVector[1] / sumMatrElement;
					centerZ += inputMatrix[row][col][depth] * tempVector[2] / sumMatrElement;

					x += tempVector[0] / nrPixels;
					y += tempVector[1] / nrPixels;
					z += tempVector[2] / nrPixels;
				}
			}
		}
	}
	centerOfMassShift = sqrt(pow((x - centerX), 2) + pow((y - centerY), 2)
		+ pow((z - centerZ), 2));
}



/*!
In the function calculateMoransI, MoransI and Gearys C are calculated. \n
If the VOI is very big, the image is downsampled before this calculation to avoid very long computational time.\n
Still, for big VOIs, the calculation of these two features can take some time.
*/
template <class T, size_t R>
void  MorphologicalFeatures<T, R>::calculateMoransI(boost::multi_array<T, R> inputMatrix) {
	T actualElement;
	int rowDiff;
	int colDiff;
	int depthDiff;
	int nrElements = 0;
	float dist;
	T actMean;
	T sum = 0;
	T sumGeary = 0;
	T sumDivident = 0;
	float sumDistances = 0;
	//sum over all matrix elements
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row < inputMatrix.shape()[0]; row++) {
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				//actal matrix element, include it in the calculation only if it is not NAN
				actualElement = inputMatrix[row][col][depth];
				if (!std::isnan(actualElement)) {

					nrElements += 1;
					actMean = (actualElement - meanValue);
					sumDivident += pow(actMean, 2);
					//second sum over all elements
					for (int depth2 = 0; depth2 < inputMatrix.shape()[2]; depth2++) {
						for (int row2 = 0; row2 < inputMatrix.shape()[0]; row2++) {
							for (int col2 = 0; col2 < inputMatrix.shape()[1]; col2++) {
								if ((depth2 != depth || row2 != row || col2 != col) && !(std::isnan(inputMatrix[row2][col2][depth2]))) {
									//calculate euclidean distance between voxels
									rowDiff = abs(row2 - row);
									colDiff = abs(col2 - col);
									depthDiff = abs(depth2 - depth);
									dist = 1 / calculateEuclideanDistance(rowDiff, colDiff, depthDiff);
									//sum up weighting factor (eucl distance) * grey level -mean value * actmean
									sum += dist*(inputMatrix[row2][col2][depth2] - meanValue) * actMean;
									sumGeary += dist*pow((actualElement - inputMatrix[row2][col2][depth2]), 2);
									sumDistances += dist;
								}
							}
						}
					}
				}
			}
		}
	}
	moransI = nrElements*sum / (sumDistances*sumDivident);
	gearysC = (nrElements - 1)*sumGeary / (2 * sumDistances*sumDivident);
}

/*!
Calculate the volume density of the aligned enclosing ellipsoid.
*/
template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateVolDensityAEE() {
	volDensityAEE = 3 * volume / (pi*majorAxisLength*minorAxisLength*leastAxisLength / 2);
}

/*!
In order to calculate the legendre polynom, the binomial coefficients are needed.
*/
template <class T, size_t R>
int MorphologicalFeatures<T, R>::binomialCoefficient(int n, int k) {
	int prod = 1;
	if (n == k) {
		prod = prod * 1;
	}
	else {
		for (int j = 1; j < k; j++) {
			prod = prod * (n + 1 - j) / j;
		}
	}
	return prod;
}

/*!
In order to calculate the area density of the minimal enclosing ellipsoid, the legendre polynom is needed.
*/
template <class T, size_t R>
float MorphologicalFeatures<T, R>::legendrePolynom(float x, int exponent) {
	float sum = 0;
	T Pn = x;
	T Pnmin1=1;
	if (exponent == 0) {
		sum = 1;
	}
	else if (exponent == 1) {
		sum = x;
	}
	else if(exponent>1) {
		for (int pol = 1; pol < exponent; pol++) {
			sum = float(2 * pol + 1) / float(pol + 1) * x * Pn - float(pol) / float(pol + 1) * Pnmin1;
			Pnmin1 = Pn;
			Pn = sum;
			
		}
	}
	return sum;
}

/*!
The area density of the aligned enclosing ellipsoid is calculated in this function. \n
Herefore, the minor, major and least axis length are needed which have already been calculated before.\n
Furthermore the legendre polynom is used in this function. As it is enough to sum the first 20 parts of this polynom,
we only calculate them.
*/
template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateAreaDensityAEE() {
	float alpha = sqrt(1 - pow((minorAxisLength / majorAxisLength), 2));
	float beta = sqrt(1 - pow((leastAxisLength / majorAxisLength), 2));
	float areaEllipse = pi*majorAxisLength*minorAxisLength;
	float x = (pow(alpha, 2) + pow(beta, 2)) / (2 * alpha*beta);
	areaDensityAEE = 0;
	for (int nu = 0; nu < 21; nu++) {
		T legendre = legendrePolynom(x, nu);
		//areaDensityAEE += pi * minorAxisLength * majorAxisLength*  boost::math::legendre_p(nu, x) * pow((alpha * beta), nu) / (1 - 4 * pow(nu, 2));
		areaDensityAEE += pi * minorAxisLength * majorAxisLength * legendrePolynom(x, nu) * pow((alpha * beta), nu) / (1 - 4 * pow(nu, 2));
	}
	areaDensityAEE = surface / areaDensityAEE;
}

/*!
The volume density of the minimal enclosing ellipsoid is calculated in this function. \n
Herefore, the minor, major and least axis length are needed which have already been calculated before.
*/
template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateVolDensityMEE() {
	//I divide by 24 and not by 3 because I have to divide majorAxisLength etc all by 2 in order to obtain a/b/c
	volDensityMEE = (pi * majorAxisLength*minorAxisLength*leastAxisLength) / 6;
	volDensityMEE = volume / volDensityMEE;
}



/*!
The area density of the minimal enclosing ellipsoid is calculated in this function. \n
Herefore, the minor, major and least axis length are needed which have already been calculated before.
*/
template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateAreaDensityMEE() {
	float alpha = sqrt(1 - pow((minorAxisLength / majorAxisLength), 2));
	float beta = sqrt(1 - pow((leastAxisLength / majorAxisLength), 2));
	float areaEllipse = pi*majorAxisLength*minorAxisLength;
	areaDensityMEE = 3 * volume / (4 * pi*majorAxisLength*minorAxisLength*leastAxisLength);
}


template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateAllMorphologicalFeatures(MorphologicalFeatures<T, R> &morphFeatures, Image<float, 3> imageAttr, ConfigFile config) {
	const typename ImageType::SpacingType& inputSpacing = imageAttr.image->GetSpacing();
	mask = imageAttr.mask;
	imageSpacingX = inputSpacing[0];
	imageSpacingY = inputSpacing[1];
	imageSpacingZ = inputSpacing[2];
	//morphFeatures.getSurface(imageAttr.mask);
	if (config.useReSegmentation == 1) {
		morphFeatures.calculateApproximateVolume(imageAttr.imageMatrixOriginal, imageAttr.vectorOfMatrixElements);
	}
	else {
		morphFeatures.calculateApproximateVolume(imageAttr.imageMatrix, imageAttr.vectorOfMatrixElements);
	}
	morphFeatures.getLabelObjectFeatures(imageAttr.mask);
	morphFeatures.getBoundingBoxValues(imageAttr.mask);
	morphFeatures.calculateMajorAxisLength();
	morphFeatures.calculateMinorAxisLength();
	morphFeatures.calculateLeastAxisLength();
	morphFeatures.calculateSurface2Volume();
	morphFeatures.calculateCompactness1();
	morphFeatures.calculateCompactness2();
	morphFeatures.calculateSphericalDisproportion();
	morphFeatures.calculateSphericity();
	morphFeatures.calculateAsphericity();
	morphFeatures.calculateElongation();
	morphFeatures.calculateFlatness();
	morphFeatures.calculateIntegratedIntensity(imageAttr.vectorOfMatrixElements);
	morphFeatures.calculateCentreOfMassShift(imageAttr.imageMatrix, imageAttr.vectorOfMatrixElements);
	if (boost::size(imageAttr.vectorOfMatrixElements) > 7000) {
		std::cout << "Warning: The VOI is very big. To calculate MoransI and GearysC, the image is downsampled in order to speed up the calculations" << std::endl;
		int factor = std::floor(boost::size(imageAttr.vectorOfMatrixElements) / 7000) +1;
		ImageType::Pointer subImage = subsampleImage(imageAttr.image, factor);
		const typename ImageType::RegionType& inputRegion = subImage->GetLargestPossibleRegion();
		ImageType::Pointer subMask = subsampleImage(imageAttr.mask, factor);
		const typename ImageType::SizeType& inputSize = inputRegion.GetSize();
		
		const typename ImageType::SpacingType& spacing = subImage->GetSpacing();
		Image<T, R> subsampledImage(spacing[0], spacing[1], spacing[2]);
		boost::multi_array<T, R> imageMatrix = subsampledImage.get3Dimage(subImage, subMask, config);
		if (config.imageType=="PET" &&  config.useSUV == 1) {
			subsampledImage.calculateSUV(imageMatrix, config);
		}
		if (config.imageType == "PET" &&  config.useSUL == 1) {
			subsampledImage.calculateSUL(imageMatrix, config);
		}
		morphFeatures.calculateMoransI(imageMatrix);
	}
	else {
		morphFeatures.calculateMoransI(imageAttr.imageMatrix);
	}
	morphFeatures.calculateVolDensityAEE();
	morphFeatures.calculateAreaDensityAEE();
	morphFeatures.calculateAreaDensityMEE();
	morphFeatures.calculateVolDensityMEE();
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::writeCSVFileMorphological(MorphologicalFeatures<T, R> morph, string outputFolder)
{
	string csvName = outputFolder + "_morphologicalFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream morphCSV;
	morphCSV.open(name);
	vector<string> features;
	defineMorphologicalFeatures(features);

	vector<T> morphData;
	extractMorphologicalData(morphData, morph);
	for (int i = 0; i< morphData.size(); i++) {
		morphCSV << "Morphology" << "," << features[i] << ",";
		morphCSV << morphData[i];
		morphCSV << "\n";
	}
	morphCSV.close();
}


template <class T, size_t R>
void MorphologicalFeatures<T, R>::writeOneFileMorphological(MorphologicalFeatures<T, R> morph, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream morphCSV;
	morphCSV.open(name);
	vector<string> features;
	defineMorphologicalFeatures(features);

	vector<T> morphData;
	extractMorphologicalData(morphData, morph);
	for (int i = 0; i< morphData.size(); i++) {
		morphCSV << "Morphology" << "," << features[i] << ",";
		morphCSV << morphData[i];
		morphCSV << "\n";
	}
	morphCSV.close();
}
template <class T, size_t R>
void MorphologicalFeatures<T, R>::defineMorphologicalFeatures(vector<string> &features) {
	features.push_back("Volume");
	features.push_back("approximate volume");
	features.push_back("Surface");
	features.push_back("Surface to volume ratio");
	features.push_back("Compactness1");
	features.push_back("Compactness2");
	features.push_back("Spherical disproportion");
	features.push_back("sphericity");
	features.push_back("asphericity");
	features.push_back("center of mass shift");
	features.push_back("maximum 3D diameter");
	features.push_back("major axis length");
	features.push_back("minor axis length");
	features.push_back("least axis length");
	features.push_back("elongation");
	features.push_back("flatness");
	features.push_back("vol density AABB");
	features.push_back("area density AABB");
	//features.push_back("vol density OMBB");
	//features.push_back("area density OMBB");
	features.push_back("vol density AEE");
	//features.push_back("area density AEE");
	//features.push_back("vol density MEE");
	//features.push_back("area density MEE");
	features.push_back("integrated intensity");
	features.push_back("Morans I");
	features.push_back("Gearys C");
	features.push_back("center of mass shift");
	features.push_back("maximum 3D diameter");
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::extractMorphologicalData(vector<T> &morphData, MorphologicalFeatures<T, R> morphFeatures) {
	morphData.push_back(morphFeatures.volume);
	morphData.push_back(morphFeatures.appVolume);
	morphData.push_back(morphFeatures.surface);
	morphData.push_back(morphFeatures.surface2volumeRatio);
	morphData.push_back(morphFeatures.compactness1);
	morphData.push_back(morphFeatures.compactness2);
	morphData.push_back(morphFeatures.sphericalDisproportion);
	morphData.push_back(morphFeatures.sphericity);
	morphData.push_back(morphFeatures.asphericity);
	morphData.push_back(morphFeatures.centerOfMassShift);
	morphData.push_back(morphFeatures.maximumDiameter);
	morphData.push_back(morphFeatures.majorAxisLength);
	morphData.push_back(morphFeatures.minorAxisLength);
	morphData.push_back(morphFeatures.leastAxisLength);
	morphData.push_back(morphFeatures.elongation);
	morphData.push_back(morphFeatures.flatness);
	morphData.push_back(morphFeatures.volDensityAABB);
	morphData.push_back(morphFeatures.areaDensityAABB);
	//morphData.push_back(morphFeatures.volDensityOMBB);
	//morphData.push_back(morphFeatures.areaDensityOMBB);
	morphData.push_back(morphFeatures.volDensityAEE);
	//morphData.push_back(morphFeatures.areaDensityAEE);
	//morphData.push_back(morphFeatures.volDensityMEE);
	//morphData.push_back(morphFeatures.areaDensityMEE);
	morphData.push_back(morphFeatures.integratredInt);
	morphData.push_back(morphFeatures.moransI);
	morphData.push_back(morphFeatures.gearysC);

}

#endif // MORPHOLOGICALFEATURES_H_INCLUDED
