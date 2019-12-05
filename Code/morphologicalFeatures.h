#ifndef MORPHOLOGICALFEATURES_H_INCLUDED
#define MORPHOLOGICALFEATURES_H_INCLUDED

/*! \file */


#include <cmath>
#include <boost/geometry.hpp>
#include "matrixFunctions.h"
#include "itkBinaryDilateImageFilter.h"
#include <algorithm>
#include <vector>
#include "itkTetrahedronCell.h"
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
#include "itkBinaryFillholeImageFilter.h"
#include "itkImageMaskSpatialObject.h"

#include "itkCastImageFilter.h"

#include "image.h"
#include "itkTypes.h"

#include "itkChangeInformationImageFilter.h"



//using namespace itkTypes;
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
	typedef boost::accumulators::features <tag::mean> Features;
	typedef accumulator_set <T, Features> Accumulator;
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
	itk::Vector< float, 3> principalMoments;
	typename ImageType::PointType origin;

	float imageSpacingX;
	float imageSpacingY;
	float imageSpacingZ;
	float calculateSurface(vector<vector<float> > vec);

	//two help functions that are needed in order to calculate some of the features
	//calculate the euclidean distance
	float calculateEuclideanDistance(int rowDiff, int colDiff, int depthDiff);
	//convert the rows, cols and depths to cartesian coordinates
	vector<float>  convertToCoordinates(int row, int col, int depth);
	//this function subsamples the image 
	//for big tumors this has to be done before MoransI and GearysC can be calculated
	ImageType::Pointer subsampleImage(ImageType::Pointer image, float factor);

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

	ImageType::Pointer changeMaskSpacingToImageSpacing(ImageType::Pointer image, ImageType::Pointer mask);
	ImageType::Pointer thresholdMask(ImageType::Pointer mask, ConfigFile config);
	void defineMorphologicalFeatures(vector<string> &features);
	void extractMorphologicalData(vector<T> &morphData, MorphologicalFeatures<T, R> morphFeatures);
	float getSurface(ImageType::Pointer mask);
public:
	//constructor
	MorphologicalFeatures() {
	}
	//destructor
	~MorphologicalFeatures() {
	}
	vector<T> cartCoordinateROI;
	vector<T> voxVolumeROI;
	const float pi = 3.141592653589793238463;
	boost::multi_array<vector<T>, R> coordinatesMatrix;
	void defineMorphologicalFeaturesOntology(vector<string> &features);
	void calculateAllMorphologicalFeatures(MorphologicalFeatures<T, R> &morphFeatures, Image<float, 3> imageAttr, ConfigFile config);
	void writeCSVFileMorphological(MorphologicalFeatures<T, R> morph, string outputFolder, ConfigFile config);
	void writeOneFileMorphological(MorphologicalFeatures<T, R> morph, ConfigFile config);
	

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
	xcoord = row*imageSpacingX + float(imageSpacingX) / 2;
	ycoord = col*imageSpacingY + float(imageSpacingY) / 2;
	zcoord = depth * imageSpacingZ + float(imageSpacingZ) / 2;
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
ImageType::Pointer MorphologicalFeatures<T, R>::subsampleImage(ImageType::Pointer imageR, float factor) {
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
	newWidth = floor((float)oldWidth  / (float)factor);
	newHeight = floor((float)oldHeight / (float)factor);
	newDepth = floor((float)oldDepth/ (float)factor);
	double outputSpacing[3];
	outputSpacing[0] = inputSpacing[0] * factor;
	outputSpacing[1] = inputSpacing[1] * factor;
	outputSpacing[2] = inputSpacing[2] * factor;
	imageSpacingX = outputSpacing[0];
	imageSpacingY = outputSpacing[1];
	imageSpacingZ = outputSpacing[2];
	itk::Size<3> outputSize = { { newWidth, newHeight, newDepth } };
	Image<float, 3> imageMask(outputSize[0], outputSize[1], outputSize[2]);
	image = imageMask.getResampledImage(imageR, outputSpacing, outputSize, "Linear", 1);
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
			volume += (0.998*labelObject->GetPhysicalSize()+ 0.999*labelObject->GetPhysicalSize()+0.9998*labelObject->GetPhysicalSize())/3;
			maximumDiameter += labelObject->GetFeretDiameter();
			if (labelNr == 1) {
				principalMoments = labelObject->GetPrincipalMoments();
			}
			else {
				principalMoments += labelObject->GetPrincipalMoments();
			}
			nrPixels += labelObject->GetNumberOfPixels();
			const itk::Vector<float, R> bbSize = labelObject->GetOrientedBoundingBoxSize();
			itk::Size<R> bbSizeITK = { { bbSize[0], bbSize[1], bbSize[2] } };
			calculateVADensity(volDensityOMBB, areaDensityOMBB, bbSizeITK);
#else
			volume += labelObject->GetPhysicalSize();

			principalMoments = labelObject->GetPrincipalMoments();
			nrPixels += labelObject->GetNumberOfPixels();
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
#ifdef _WIN32
	typename typedef itk::ImageMaskSpatialObject<R> ImageMaskSpatialObject;
	typename typedef ImageMaskSpatialObject::ImageType ImageTypeSpatial;
	typename typedef ImageTypeSpatial::RegionType RegionType;
	
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
#endif
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
	surface2volumeRatio = surface / volume;
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateCompactness1() {
	compactness1 = volume / (sqrt(pi)*pow(surface, 1.5));
}

template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateCompactness2() {
	compactness2 = 36 * pi*pow(volume, 2) / pow(surface, 3);
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
	float centerX = 0;
	float centerY = 0;
	float centerZ = 0;
	float x = 0;
	float y = 0;
	float z = 0;
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
//	std::cout << "matrix size" << inputMatrix.shape()[0] << " " << inputMatrix.shape()[1] << " " << inputMatrix.shape()[2] << std::endl;
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
ImageType::Pointer MorphologicalFeatures<T, R>::changeMaskSpacingToImageSpacing(ImageType::Pointer image, ImageType::Pointer mask) {
	using FilterType = itk::ChangeInformationImageFilter< ImageType >;
	FilterType::Pointer filter = FilterType::New();
	ImageType::PointType     origin = image->GetOrigin();
	ImageType::SpacingType   spacing = image->GetSpacing();
	ImageType::DirectionType direction = image->GetDirection();


	filter->SetOutputSpacing(spacing);
	filter->ChangeSpacingOn();

	filter->SetOutputOrigin(origin);
	filter->ChangeOriginOn();

	filter->SetOutputDirection(direction);
	filter->ChangeDirectionOn();

	filter->SetInput(mask);
	filter->Update();
	ImageType::Pointer newMask = filter->GetOutput();
	return newMask;
}

template <class T, size_t R>
ImageType::Pointer MorphologicalFeatures<T, R>::thresholdMask(ImageType::Pointer mask, ConfigFile config) {
	itk::ImageRegionConstIterator<ImageType> countVoxels(mask, mask->GetLargestPossibleRegion());
	
	countVoxels.GoToBegin();
	float maxValue = 0;
	while (!countVoxels.IsAtEnd())
	{   
		if (countVoxels.Get() > maxValue)
		{
			maxValue = countVoxels.Get();
			//std::cout << maxValue << std::endl;
		}
		++countVoxels;
	}
	using FilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
	
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(mask);
	//filter->SetLowerThreshold(0);
	filter->SetLowerThreshold((config.threshold)*maxValue);
	filter->SetOutsideValue(0);
	filter->SetInsideValue(1);
	filter->Update();
	ImageType::Pointer updatedMask =filter->GetOutput();
	return updatedMask;


}


template<class T, size_t R>
float MorphologicalFeatures<T, R>::getSurface(ImageType::Pointer mask) {
	const typename ImageType::SpacingType& inputSpacing = mask->GetSpacing();
	if (inputSpacing[0] > 3) {
		float factor = (2 / inputSpacing[0]);

		mask = subsampleImage(mask, factor);
	}
	using binaryFilterType = itk::BinaryFillholeImageFilter<ImageType>;
	binaryFilterType::Pointer filter = binaryFilterType::New();
	filter->SetInput(mask);
	filter->Update();
	//to calculate the surface we have to convert our mask to a mesh
	typename MeshSourceType::Pointer meshSource = MeshSourceType::New();
	const PixelType objectValue = static_cast<PixelType>(1);
	meshSource->SetObjectValue(objectValue);
	meshSource->SetInput(mask);

	try {
		meshSource->Update();
	}
	catch (itk::ExceptionObject &exp) {
		std::cerr << "Exception throwing during update()" << std::endl;
	}
	//using MeshType = itk::Mesh<double, 3>;
	
	//ensure that all cells are triangles (I do not think that is really necessary, maybe
	//it can be deleted in another version
	for (MeshType::CellsContainerIterator it = meshSource->GetOutput()->GetCells()->Begin();
		it != meshSource->GetOutput()->GetCells()->End(); ++it) {
		MeshType::CellAutoPointer cell;
		meshSource->GetOutput()->GetCell(it->Index(), cell);
		if (3 != cell->GetNumberOfPoints()) {
			std::cerr << "ERROR: All cells must be trianglar." << std::endl;
		}
	}
	
	//to calculate the surface, we have to iterate over all cells of our mesh
	Celliterator cellIterator = meshSource->GetOutput()->GetCells()->Begin();
	Celliterator cellEnd = meshSource->GetOutput()->GetCells()->End();

	typedef itk::SimplexMesh< float, 3 > TSimplex;
	typedef itk::TriangleMeshToSimplexMeshFilter< MeshType, TSimplex > TConvert;
	// Convert the triangle mesh to a simplex mesh.
	TConvert::Pointer convert = TConvert::New();
	//convert->SetInput(meshSource->GetOutput());
	//convert->Update();
	//TSimplex::Pointer testMesh = convert->GetOutput();
	//typedef itk::SimplexMeshVolumeCalculator<TSimplex> VolumeCalculatorType;

	//VolumeCalculatorType::Pointer volumeCalculator = VolumeCalculatorType::New();
	//volumeCalculator->SetSimplexMesh(convert->GetOutput());
	//volumeCalculator->Compute();
	
//	std::cout << "volume" << volumeCalculator->GetVolume() << std::endl;
	//convert->Update();
	//to really get the coordinates of the points of the cell, we need a polyDataReader
	//iterate over all cells and store the indices of the points in a vector
	//this vector is used later to calculate the surface of each cell
#ifdef _WIN32
	MeshType::PointType pointOnMesh;
	bool pointExists;
	float testtest = 0;
	while (cellIterator != cellEnd) {
		CellType *cell = cellIterator.Value();
		typedef CellType::PointIdIterator PointIdIterator;
		PointIdIterator pointIdIter = cell->PointIdsBegin();
		PointIdIterator pointIdend = cell->PointIdsEnd();
		//testtest = cell->IsExplicitBoundary();
		while (pointIdIter != pointIdend) {
			
			pointExists = meshSource->GetOutput()->GetPoint(*pointIdIter, &pointOnMesh);
			if (pointExists) {

				coordinates.push_back(pointOnMesh[0]);
				coordinates.push_back(pointOnMesh[1]);
				coordinates.push_back(pointOnMesh[2]);
				vec.push_back(coordinates);
				coordinates.clear();
			}
			++pointIdIter;
		}
		++cellIterator;
	}
	surface = calculateSurface(vec);

#else
	ReaderTypeVTK::Pointer polyDataReader = ReaderTypeVTK::New();
	PointType pointOnMesh;
	bool pointExists;

	//now iterate over all cells and store the coordinates of each cell in a vector
	while (cellIterator != cellEnd) {
		CellType *cell = cellIterator.Value();
		typedef CellType::PointIdIterator PointIdIterator;
		PointIdIterator pointIdIter = cell->PointIdsBegin();
		PointIdIterator pointIdend = cell->PointIdsEnd();
		while (pointIdIter != pointIdend) {
			pointExists = meshSource->GetOutput()->GetPoint(*pointIdIter, &pointOnMesh);
			if (pointExists) {
				coordinates.push_back(pointOnMesh[0]);
				coordinates.push_back(pointOnMesh[1]);
				coordinates.push_back(pointOnMesh[2]);
				vec.push_back(coordinates);
				coordinates.clear();
			}
			++pointIdIter;
		}
		++cellIterator;
	}
	surface = calculateSurface(vec);
#endif
	return surface;
}

//calculate the surface of a mesh
//as input the function takes a vector that consists of the points of each cells
//these points are used to calculate the surface of each cell
//the surfaces of the different cells are added
template<class T, size_t R>
float MorphologicalFeatures<T, R>::calculateSurface(vector<vector<float> > vec) {
	//vector<float> a;
	//vector<float> b;
	//vector<float> c;
	vector<float> ab;
	vector<float> ac;
	vector<float> bc;
	//vector<float> cross;
	double m_Volume, m_VolumeX, m_VolumeY, m_VolumeZ = 0.0;
	double m_Area, m_Kx, m_Ky, m_Kz = 0.0;
	double m_Wxyz, m_Wxy, m_Wxz, m_Wyz, m_Muncx, m_Muncy, m_Muncz = 0;
	double area;
	double a, b, c, s;
	int m_NumberOfTriangles = 0;
	surface = 0;
	double i[3] , j[3], k[3], u[3], absu[3], length;
	double ii[3], jj[3], kk[3];
	double xavg, yavg, zavg;
	vector<T> tmpPoint;
	double p1[3], p2[3], p3[3];
	double cross[3];
	vector<T> centerPoint = vec[0];
	float surface = 0;
	for (int nr = 0; nr < (vec.size() / 3) - 1; nr++) {
		tmpPoint = vec[3 * nr];
		p1[0] = tmpPoint[0] -centerPoint[0];
		p1[1] = tmpPoint[1] -centerPoint[1];
		p1[2] = tmpPoint[2] -centerPoint[2];
		//std::cout << "P1" <<" "<< p1[0] << " " << p1[1] << " " << p1[2] << " ";
		tmpPoint = vec[3 * nr + 1];
		p2[0] = tmpPoint[0] -centerPoint[0];
		p2[1] = tmpPoint[1] -centerPoint[1];
		p2[2] = tmpPoint[2] -centerPoint[2];
		//std::cout << "P2" <<" "<< p2[0] << " " << p2[1] << " " << p2[2] << " ";
		tmpPoint = vec[3 * nr + 2];
		p3[0] = tmpPoint[0] -centerPoint[0];
		p3[1] = tmpPoint[1] -centerPoint[1];
		p3[2] = tmpPoint[2] -centerPoint[2];
		//std::cout << "P3" <<" "<< p3[0] << " " << p3[1] << " " << p3[2] << std::endl;
		//cross[0] =p1[0]* (p2[1] * p3[2] - p3[2] * p2[1]);
		//cross[1] = p1[1] * (p3[2] * p2[0] - p3[0] * p2[1]);
		//cross[2] = p1[2] * (p2[0] * p3[1] - p3[0] * p2[1]);

		//std::cout << "vol" << (cross[0] + cross[1] + cross[2]) / 6 << std::endl;
		// Get i j k vectors ...
  ////
		i[0] = (p2[0] -p1[0]); j[0] = (p2[1] - p1[1]); k[0] =  (p2[2] - p1[2]);
		i[1] = (p3[0] - p1[0]); j[1] = (p3[1] - p1[1]); k[1] =  (p3[2] - p1[2]);
		i[2] = (p3[0] - p2[0]); j[2] =  (p3[1] - p2[1]); k[2] =  (p3[2] - p2[2]);

		// Cross product between two vectors, to determine normal vector
		//
		u[0] = (j[0] * k[1] - k[0] * j[1]);
		u[1] = (k[0] * i[1] - i[0] * k[1]);
		u[2] = (i[0] * j[1] - j[0] * i[1]);

		// Normalize normal
		//
		length = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
		if (length != 0.0)
		{
			u[0] /= length;
			u[1] /= length;
			u[2] /= length;
		}
		else
		{
			u[0] = u[1] = u[2] = 0.0;
		}
			//// Determine max unit normal component...
		////
		absu[0] = std::fabs(u[0]); absu[1] = std::fabs(u[1]); absu[2] = std::fabs(u[2]);
		if ((absu[0] > absu[1]) && (absu[0] > absu[2]))
		{
			m_Muncx++;
		}
		else if ((absu[1] > absu[0]) && (absu[1] > absu[2]))
		{
			m_Muncy++;
		}
		else if ((absu[2] > absu[0]) && (absu[2] > absu[1]))
		{
			m_Muncz++;
		}
		else if (itk::Math::AlmostEquals(absu[0], absu[1]) && itk::Math::AlmostEquals(absu[0], absu[2]))
		{
			m_Wxyz++;
		}
		else if (itk::Math::AlmostEquals(absu[0], absu[1]) && (absu[0] > absu[2]))
		{
			m_Wxy++;
		}
		else if (itk::Math::AlmostEquals(absu[0], absu[2]) && (absu[0] > absu[1]))
		{
			m_Wxz++;
		}
		else if (itk::Math::AlmostEquals(absu[1], absu[2]) && (absu[0] < absu[2]))
		{
			m_Wyz++;
		}
		
		// This is reduced to ...
		ii[0] = i[0] * i[0]; ii[1] = i[1] * i[1]; ii[2] = i[2] * i[2];
		jj[0] = j[0] * j[0]; jj[1] = j[1] * j[1]; jj[2] = j[2] * j[2];
		kk[0] = k[0] * k[0]; kk[1] = k[1] * k[1]; kk[2] = k[2] * k[2];

		
		// Area of a triangle using Heron's formula...
		a = std::sqrt(ii[1] + jj[1] + kk[1]);
		b = std::sqrt(ii[0] + jj[0] + kk[0]);
		c = std::sqrt(ii[2] + jj[2] + kk[2]);
		s = 0.5 * (a + b + c);
		area = std::sqrt(std::fabs(s * (s - a) * (s - b) * (s - c)));
		// Volume elements ...
		zavg = (p1[2] + p2[2] + p3[2]) / 3.0;
		yavg = (p1[1] + p2[1] + p3[1]) / 3.0;
		xavg = (p1[0] + p2[0] + p3[0]) / 3.0;
		m_VolumeX += (area * (double)u[2] * (double)zavg);
		m_VolumeY += (area * (double)u[1] * (double)yavg);
		m_VolumeZ += (area * (double)u[0] * (double)xavg);
		//S	rface += (cross[0] + cross[1] + cross[2]) / 6;
		surface += area;
		m_NumberOfTriangles++;
		
	}
	// Compute fraction of elements that primarily point along the x, y
  // and z directions
	m_Kx = (m_Muncx + (m_Wxyz / 3.0) + ((m_Wxy + m_Wxz) / 2.0)) / m_NumberOfTriangles;
	m_Ky = (m_Muncy + (m_Wxyz / 3.0) + ((m_Wxy + m_Wyz) / 2.0)) / m_NumberOfTriangles;
	m_Kz = (m_Muncz + (m_Wxyz / 3.0) + ((m_Wxz + m_Wyz) / 2.0)) / m_NumberOfTriangles;

	m_Volume = (m_Kx * m_VolumeX
		+ m_Ky * m_VolumeY
		+ m_Kz * m_VolumeZ);
	m_Volume = std::fabs(m_Volume);
	volume = m_Volume;
	return surface;
}


template <class T, size_t R>
void MorphologicalFeatures<T, R>::calculateAllMorphologicalFeatures(MorphologicalFeatures<T, R> &morphFeatures, Image<float, 3> imageAttr, ConfigFile config) {
	const typename ImageType::SpacingType& inputSpacing = imageAttr.image->GetSpacing();
	mask = imageAttr.mask;
	
	mask = changeMaskSpacingToImageSpacing(imageAttr.image, mask);
	ImageType::SpacingType   spacingMask = mask->GetSpacing();
	imageSpacingX = inputSpacing[0];
	imageSpacingY = inputSpacing[1];
	imageSpacingZ = inputSpacing[2];
	ImageType::Pointer updatedMask = thresholdMask(mask, config);
	//morphFeatures.getSurface(imageAttr.mask);
	if (config.useReSegmentation == 1) {
		morphFeatures.calculateApproximateVolume(imageAttr.imageMatrixOriginal, imageAttr.vectorOfMatrixElements);
	}
	else {
		morphFeatures.calculateApproximateVolume(imageAttr.imageMatrix, imageAttr.vectorOfMatrixElements);
	}
	itk::ImageRegionConstIterator<ImageType> countVoxels(mask, mask->GetLargestPossibleRegion());

	
	morphFeatures.getLabelObjectFeatures(updatedMask);
	surface = getSurface(updatedMask);
	morphFeatures.getBoundingBoxValues(updatedMask);
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
		int factor = std::floor(10/(imageSpacingX));
		
		vector<T> subsampleElements;
		
		ImageType::Pointer subImage = subsampleImage(imageAttr.image, factor);
		const typename ImageType::RegionType& inputRegion = subImage->GetLargestPossibleRegion();
		ImageType::Pointer subMask = subsampleImage(imageAttr.mask, factor);
		const typename ImageType::SizeType& inputSize = inputRegion.GetSize();
		const typename ImageType::SpacingType& spacing = subImage->GetSpacing();
		Image<T, R> subsampledImage(spacing[0], spacing[1], spacing[2]);
		boost::multi_array<T, R> imageMatrix = subsampledImage.get3Dimage(subImage, subMask, config);
		subsampleElements = subsampledImage.getVectorOfMatrixElementsNotNAN(imageMatrix);
		Accumulator acc;
		for_each(subsampleElements.begin(), subsampleElements.end(), boost::bind<void>(boost::ref(acc), _1));
		meanValue = mean(acc);
		if (config.imageType == "PET" &&  config.useSUV == 1) {
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
void MorphologicalFeatures<T, R>::writeCSVFileMorphological(MorphologicalFeatures<T, R> morph, string outputFolder, ConfigFile config)
{
	string csvName = outputFolder + "_morphologicalFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream morphCSV;
	if (config.imageType != "PET") {
		morphCSV.open(name);
	}
	else {
		morphCSV.open(name, std::ios_base::app);
	}
	
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
void MorphologicalFeatures<T, R>::writeOneFileMorphological(MorphologicalFeatures<T, R> morph, ConfigFile config) {
	string csvName;
	if (config.getOneCSVFile == 1) {
		csvName = config.outputFolder + ".csv";		
	}
	else if (config.ontologyOutput == 1) {
		csvName = config.outputFolder + "/feature_table.csv";
	}
	
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';
	ofstream morphCSV;
	if (config.imageType != "PET") {
		morphCSV.open(name);
	}
	else {
		morphCSV.open(name, std::ios_base::app);
	}
	vector<string> features;
	vector<T> morphData;
	extractMorphologicalData(morphData, morph);
	if (config.getOneCSVFile == 1) {
		defineMorphologicalFeatures(features);
		for (int i = 0; i < morphData.size(); i++) {
			morphCSV << "Morphology" << "," << features[i] << ",";
			morphCSV << morphData[i];
			morphCSV << "\n";
		}
		morphCSV.close();
	}
	else if (config.ontologyOutput == 1) {
		defineMorphologicalFeaturesOntology(features);
		for (int i = 0; i < morphData.size(); i++) {
			
			morphCSV << config.patientID<<","<<config.patientLabel<<","<< features[i] << ",";
			morphCSV << morphData[i]<<"," << config.featureParameterSpaceName<<","<<config.calculationSpaceName;
			morphCSV << "\n";
		}
		morphCSV.close();
	}
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
}

template <class T, size_t R>
void  MorphologicalFeatures<T, R>::defineMorphologicalFeaturesOntology(vector<string> &features) {
	features.push_back("Fmorph.vol");
	features.push_back("Fmorph.approx.vol");
	features.push_back("Fmorph.area");
	features.push_back("Fmorph.av");
	features.push_back("Fmorph.comp.1");
	features.push_back("Fmorph.comp.2");
	features.push_back("Fmorph.sph.dispr");
	features.push_back("Fmorph.sphericity");
	features.push_back("Fmorph.asphericity");	
	features.push_back("Fmorph.com");	
	features.push_back("Fmorph.diam");
	features.push_back("Fmorph.pca.major");
	features.push_back("Fmorph.pca.minor");
	features.push_back("Fmorph.pca.least");
	features.push_back("Fmorph.pca.elongation");
	features.push_back("Fmorph.pca.flatness");	
	features.push_back("Fmorph.v.dens.aabb");
	features.push_back("Fmorph.a.dens.aabb");
	//features.push_back("Fmorph.v.dens.ombb");
	//features.push_back("Fmorph.a.dens.ombb");
	features.push_back("Fmorph.v.dens.aee");
	//features.push_back("Fmorph.a.dens.aee");
	//features.push_back("Fmorph.v.dens.mee");
	//features.push_back("Fmorph.a.dens.mee");
	features.push_back("Fmorph.integ.int");		
	features.push_back("Fmorph.moransI");
	features.push_back("Fmorph.gearysC");
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
