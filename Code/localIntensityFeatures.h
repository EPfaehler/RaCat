#ifndef LOCALINTENSITYFEATURES_H_INCLUDED
#define LOCALINTENSITYFEATURES_H_INCLUDED

#include "itkConvolutionImageFilter.h"
#include "image.h"
#include "morphologicalFeatures.h"
#include "matrixFunctions.h"


/*! \file */
/*!
In the class LocalIntensityFeatures only the local and global peak are calculated. \n
In order to calculate global/local peak, a convolutional matrix is calculated. Multiplying by this matrix
gives the circle around the present voxel. \n
Herefore, dependent on voxel size, the size of the convolutional matrix is calculated.\n
The calculation of the feature values is done before discretizing the matrix values to a user specified bin number.
*/

template<class T, size_t R>
class LocalIntensityFeatures {
private:
	T maxElement;
	T localIntensityPeak;
	T globalIntensityPeak;
	float spacingX;
	float spacingY;
	float spacingZ;
	float volVoxel;
	MorphologicalFeatures<T, R> morpho;
	typedef boost::multi_array<T, R> locMatrix;
	float pi = 3.1415926535;
	float originalRadius = 6.2;
	float voxelSize[3];
	//get indices of all maximal elements
	vector<vector<int> > getIndexOfMax(boost::multi_array<T, R> inputMatrix);
	//vector where all values inside the circle are stored
	vector<T> intValuesInCircle;
	//get the size of the convolutional matrix
	void getConvMatrixSize(ImageType::Pointer mask, int &nrVoxelsDirection, float spacing, float radius);
	//fill the convolutional matrix
	void fillConvMatrix(boost::multi_array<T, R> &matrix, ImageType::Pointer mask);
	void fillVector(vector<double> &index, boost::multi_array<T, R> convMatrix);
	ImageType::Pointer calculatePeaks(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> localIntMatrix, boost::multi_array<T, R> convolutionalMatrix, ImageType::Pointer image);
	boost::multi_array<T, R> calculateConvolutionMatrix(ImageType::Pointer mask);
	float calculatePeakValues(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> convolutionalMatrix, vector<int> nrValues, ImageType::Pointer image);
	void calculateLocalIntensityPeak(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> peakMatrix);
	void calculateGlobalIntensityPeak(boost::multi_array<T, R> peakMatrix, boost::multi_array<T, R> resegmentedMask);

	void defineLocalIntenseFeatures(vector<string> &features);
	void extractLocalIntenseData(vector<T> &localIntData, LocalIntensityFeatures<T, R> localIntFeatures);

public:
	LocalIntensityFeatures() {
	}
	~LocalIntensityFeatures() {
	}
	void calculateAllLocalIntensityFeatures(LocalIntensityFeatures<T, R> &localInt, Image<T, R> imageAttr, ConfigFile config);
	void writeCSVFileLocalIntensity(LocalIntensityFeatures<T, R> localInt, string outputfolder);
	void writeOneFileLocalInt(LocalIntensityFeatures<T, R> localInt, string outputfolder);

};


/*!
In the function getIndexOfMax all indices of the voxels that yield the maximum intensity value are stored in a vector. \n
Herefore, a search is performed in the input matrix in order to look for the maximum element.\n
@parameter[in]: boost::multi_array inputMatrix
@parameter[out]: vector containing indices
*/
template<class T, size_t R>
vector<vector<int> > LocalIntensityFeatures<T, R>::getIndexOfMax(boost::multi_array<T, R> inputMatrix) {
	vector<int> actIndex;
	//vector where all indices of the maximal indices are stored (if there are more than one)
	vector<vector<int> > allMaxIndices;
	for (int depth = 0; depth < inputMatrix.shape()[2]; ++depth) {
		for (int row = 0; row < inputMatrix.shape()[0]; ++row) {
			for (int col = 0; col < inputMatrix.shape()[1]; ++col) {
				//if element is maximal element, put index in array
				if (inputMatrix[row][col][depth] == maxElement) {
					actIndex.push_back(row);
					actIndex.push_back(col);
					actIndex.push_back(depth);
					allMaxIndices.push_back(actIndex);
					actIndex.clear();
				}
			}
		}
	}
	return allMaxIndices;
}


/*!
In the function getConvMatrixSize, it is determined which size the convolutional matrix will have
@parameter[in]: ImageType::Pointer mask: image mask
@parameter[in]: int nrVoxelsDirection: as reference: how many vel we go in one direction
@parameter[in]: float spacing: image spacing
@parameter[in]: float radius: radius of circle
*/
template<class T, size_t R>
void LocalIntensityFeatures<T, R>::getConvMatrixSize(ImageType::Pointer mask, int &nrVoxelsDirection, float spacing, float radius) {
	//now get the number of voxels which are fitting 100% in the matrix for all directions
	nrVoxelsDirection = floor((radius-double(spacing)/2.0)  / spacing);
	//nrVoxelsDirection = floor(radius/ spacing);
	if (fmod(radius, nrVoxelsDirection) != 0) {
		//if we have some fractional voxels in x direction, the convolutional matrix has also to include the element of the fractional voxel
		nrVoxelsDirection += 1;
	}
	//*2 because we have to go in two directions (have radius, not diameter) +1 to also include center voxel
	nrVoxelsDirection = 2 * nrVoxelsDirection + 1;
}



/*!
In the function fillConvMatrix, we fill the convolutional matrix\n
For this, we start in the center of the matrix. In the center of the matrix, the radius is the actual radius of the sphere.\n
Then we go one slice in x-direction. The radius changes, as we want to have a sphere with center in the center of the matrix.\n
The matrix is filled slice by slice.\n
@parameter[in]: boost multi_array convMatrix: the convolutional matrix is given as reference
@parameter[in]: ImageType::Pointer mask: the image mask
*/
template<class T, size_t R>
void LocalIntensityFeatures<T, R>::fillConvMatrix(boost::multi_array<T, R> &convMatrix, ImageType::Pointer mask) {
	const typename ImageType::SpacingType& inputSpacing = mask->GetSpacing();
	//get input spacing
	vector<double> voxelSpacingX = { { inputSpacing[1], inputSpacing[2] } };
	vector<double> voxelSpacingY = { { inputSpacing[0], inputSpacing[2] } };
	vector<double> voxelSpacingZ = { { inputSpacing[0], inputSpacing[1] } };
	//float volVoxel = inputSpacing[0] * inputSpacing[1] * inputSpacing[2];
	float frac;
	//get the center index
	vector<double> indexOfCenter;
	
	fillVector(indexOfCenter, convMatrix);
	vector<double> actPos = { {0,0,0} };
	double distFromCenter;
	int nrVoxels = 0;
	int test = 0;
	for (int depth = 0; depth < convMatrix.shape()[2]; ++depth) {
		for (int row = 0; row < convMatrix.shape()[0]; ++row) {
			for (int col = 0; col < convMatrix.shape()[1]; ++col) {
				actPos[0] = double(row)*inputSpacing[0] + inputSpacing[0]/2;
				actPos[1] = double(col)*inputSpacing[1] + inputSpacing[1]/2 ;
				actPos[2] = double(depth)*inputSpacing[2] + inputSpacing[2]/2 ;
				distFromCenter = std::sqrt(std::pow((actPos[0] - (indexOfCenter[0]+0.5)*inputSpacing[0] ), 2) + std::pow((actPos[1]  - (indexOfCenter[1]+0.5) * inputSpacing[1] ), 2) + std::pow((actPos[2] - (indexOfCenter[2]+0.5) * inputSpacing[2] ), 2));
				if (distFromCenter <= 6.2) {
					convMatrix[row][col][depth] = 1;
				}
				
			}
		}
	}

}



//fill a vector with indices
template<class T, size_t R>
void LocalIntensityFeatures<T, R>::fillVector(vector<double> &index, boost::multi_array<T, R> convMatrix) {
	index.push_back(double(convMatrix.shape()[0] -1 )/ 2.0);
	index.push_back(double(convMatrix.shape()[1] -1)/ 2.0);
	index.push_back(double(convMatrix.shape()[2] -1)/ 2.0);
}


/*!
In the function calculateConvolutionMatrix, the convolutional matrix is filled using the function fillConvMatrix
@parameter[in]: ImageType mask: image mask
@parameter[out]: boost multi_array convolutional matrix
*/
template<class T, size_t R>
boost::multi_array<T, R> LocalIntensityFeatures<T, R>::calculateConvolutionMatrix(ImageType::Pointer mask) {
	const typename ImageType::SpacingType& inputSpacing = mask->GetSpacing();
	//original radius of 1 cm3 sphere in mm
	float radius = originalRadius;
	//now get the number of voxels which are fitting 100% in the matrix for all directions

	int nrVoxelsDirectionX;
	int nrVoxelsDirectionY;
	int nrVoxelsDirectionZ;
	//get the size of the convolution matrix
	getConvMatrixSize(mask, nrVoxelsDirectionX, inputSpacing[0], originalRadius);
	getConvMatrixSize(mask, nrVoxelsDirectionY, inputSpacing[1], originalRadius);
	getConvMatrixSize(mask, nrVoxelsDirectionZ, inputSpacing[2], originalRadius);

	locMatrix convolutionMatrix(boost::extents[nrVoxelsDirectionX][nrVoxelsDirectionY][nrVoxelsDirectionZ]);
	fillConvMatrix(convolutionMatrix, mask);

	return convolutionMatrix;

}



/*!
In the function calculatePeaks calculates a matrix. \n
Every matrix element is the corresponding peak value to the element in the input matrix. \n
@parameter[in]: boost multi_array input matrix
@parameter[in]: boost multi_array convolutional matrix
@parameter[out]: boost multi_array peak matrix
*/
template<class T, size_t R>
ImageType::Pointer LocalIntensityFeatures<T, R>::calculatePeaks(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> localIntMatrix, boost::multi_array<T, R> convolutionalMatrix, ImageType::Pointer image) {
//boost::multi_array<T, R> LocalIntensityFeatures<T, R>::calculatePeaks(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> localIntMatrix, boost::multi_array<T, R> convolutionalMatrix, ImageType::Pointer image) {
	int peakMatrixX = inputMatrix.shape()[0];
	int peakMatrixY = inputMatrix.shape()[1];
	int peakMatrixZ = inputMatrix.shape()[2];
	vector<T> convVector;
	for (int depth = 0; depth < convolutionalMatrix.shape()[2]; depth++) {
		for (int row = 0; row < convolutionalMatrix.shape()[0]; row++) {
			for (int col = 0; col < convolutionalMatrix.shape()[1]; col++) {
				convVector.push_back(convolutionalMatrix[row][col][depth]);
			}
		}
	}		
	float *kernelArray = &convVector[0];
	unsigned int dimKernel[] = { convolutionalMatrix.shape()[0], convolutionalMatrix.shape()[1], convolutionalMatrix.shape()[2] };
	ImageType::Pointer kernel = converArray2Image(kernelArray, dimKernel, voxelSize);
	using FilterType = itk::ConvolutionImageFilter<ImageType>;
	FilterType::Pointer convolutionFilter = FilterType::New();
	convolutionFilter->SetNormalize(TRUE);
	convolutionFilter->SetInput(image);
	convolutionFilter->SetKernelImage(kernel);
	convolutionFilter->Update();
	ImageType::Pointer peakImage = convolutionFilter->GetOutput();

	
	return peakImage;
}




/*!
In the function calculateLocalIntensityPeak the local intensity peak is calculated using the peak matrix. \n
@parameter[in]: boost multi_array input matrix
@parameter[in]: boost multi_array peak matrix
*/
template<class T, size_t R>
void LocalIntensityFeatures<T, R>::calculateLocalIntensityPeak(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> peakMatrix) {
	vector<vector<int> > allMaxIndices = getIndexOfMax(inputMatrix);
	vector<int> actualIndex;
	T mean = 0;
	//iterate over all indices that have the intensity value equal to the maximum element
	//put all elements that are lying inside this circle on a vector
	for (int i = 0; i < allMaxIndices.size(); i++) {
		actualIndex = allMaxIndices[i];
		mean += peakMatrix[actualIndex[0]][actualIndex[1]][actualIndex[2]];
	}
	localIntensityPeak = mean / allMaxIndices.size();
}

/*!
In the function calculateGlobalIntensityPeak the global intensity peak is calculated using the peak matrix. \n
@parameter[in]: boost multi_array peak matrix
*/
template<class T, size_t R>
void LocalIntensityFeatures<T, R>::calculateGlobalIntensityPeak(boost::multi_array<T, R> peakMatrix, boost::multi_array<T, R> resegmentedMask) {
	vector<T> peakValues;
	for (int depth = 0; depth < peakMatrix.shape()[2]; depth++) {
		for (int row = 0; row < peakMatrix.shape()[0]; row++) {
			for (int col = 0; col < peakMatrix.shape()[1]; col++) {
				if (!std::isnan(peakMatrix[row][col][depth]) && !std::isnan(resegmentedMask[row][col][depth])) {
					peakValues.push_back(peakMatrix[row][col][depth]);

				}
			}
		}
	}
	globalIntensityPeak = *max_element(peakValues.begin(), peakValues.end());

}

template<class T, size_t R>
void LocalIntensityFeatures<T, R>::calculateAllLocalIntensityFeatures(LocalIntensityFeatures<T, R> &localInt, Image<T,R> imageAttr, ConfigFile config) {
	const typename ImageType::SpacingType& inputSpacing = imageAttr.mask->GetSpacing();
	voxelSize[0] = inputSpacing[0];
	voxelSize[1] = inputSpacing[1];
	voxelSize[2] = inputSpacing[2];
	maxElement = *max_element(imageAttr.vectorOfMatrixElements.begin(), imageAttr.vectorOfMatrixElements.end());
	vector<vector<int> > allMaxIndices = getIndexOfMax(imageAttr.imageMatrix);
	boost::multi_array<T, R> convMatrix = localInt.calculateConvolutionMatrix(imageAttr.mask);
	//boost::multi_array<T, R> peakMatrix = calculatePeaks(imageAttr.imageMatrix, imageAttr.imageMatrixLocalInt, convMatrix, imageAttr.image);
	ImageType::Pointer peakImage = calculatePeaks(imageAttr.imageMatrix, imageAttr.imageMatrixLocalInt, convMatrix, imageAttr.image);
	const typename ImageType::RegionType region = peakImage->GetBufferedRegion();

	const typename ImageType::SizeType imageSize = region.GetSize();
	Image<T, R> imageElement(imageSize[0], imageSize[1], imageSize[2]);
	boost::multi_array<T, R> peakMatrix = imageElement.get3Dimage(peakImage, imageAttr.mask, config);
	//	
	localInt.calculateLocalIntensityPeak(imageAttr.imageMatrix, peakMatrix);
	localInt.calculateGlobalIntensityPeak(peakMatrix, imageAttr.imageMatrix);

}


template<class T, size_t R>
void LocalIntensityFeatures<T, R>::writeCSVFileLocalIntensity(LocalIntensityFeatures<T, R> localInt, string outputFolder) {
	string csvName = outputFolder + "_localIntensity.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream localIntCSV;
	localIntCSV.open(name);
	vector<string> features;
	defineLocalIntenseFeatures(features);

	vector<T> localIntData;
	extractLocalIntenseData(localIntData, localInt);
	for (int i = 0; i< localIntData.size(); i++) {
		localIntCSV << "Local intensity" << "," << features[i] << ",";
		localIntCSV << localIntData[i];
		localIntCSV << "\n";
	}
	localIntCSV.close();
}

template <class T, size_t R>
void LocalIntensityFeatures<T, R>::writeOneFileLocalInt(LocalIntensityFeatures<T, R> localInt, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream localIntCSV;
	localIntCSV.open(name, std::ios_base::app);
	vector<string> features;
	defineLocalIntenseFeatures(features);

	vector<T> localIntData;
	extractLocalIntenseData(localIntData, localInt);
	for (int i = 0; i< localIntData.size(); i++) {
		localIntCSV << "Local intensity" << "," << features[i] << ",";
		localIntCSV << localIntData[i];
		localIntCSV << "\n";
	}
	localIntCSV.close();
}


template <class T, size_t R>
void LocalIntensityFeatures<T, R>::defineLocalIntenseFeatures(vector<string> &features) {
	features.push_back("local intensity peak");
	features.push_back("global intensity peak");

}

template <class T, size_t R>
void LocalIntensityFeatures<T, R>::extractLocalIntenseData(vector<T> &localIntData, LocalIntensityFeatures<T, R> localIntFeatures) {
	localIntData.push_back(localIntFeatures.localIntensityPeak);
	localIntData.push_back(localIntFeatures.globalIntensityPeak);
}

#endif
