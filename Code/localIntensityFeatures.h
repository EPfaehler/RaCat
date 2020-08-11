#ifndef LOCALINTENSITYFEATURES_H_INCLUDED
#define LOCALINTENSITYFEATURES_H_INCLUDED

#include "itkConvolutionImageFilter.h"
#include "image.h"

#include "morphologicalFeatures.h"
#include "matrixFunctions.h"
#include <cmath>


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
	vector<itk::Index<R>  > getIndexOfMax(ImageType::Pointer image, ImageType::Pointer mask);
	//vector where all values inside the circle are stored
	vector<T> intValuesInCircle;
	//get the size of the convolutional matrix
	void getConvMatrixSize(ImageType::Pointer mask, int &nrVoxelsDirection, float spacing, float radius);
	//fill the convolutional matrix
	void fillConvMatrix(boost::multi_array<T, R> &matrix, ImageType::Pointer mask);
	void fillVector(vector<float> &index, boost::multi_array<T, R> convMatrix);
	ImageType::Pointer calculatePeaks(boost::multi_array<T, R> convolutionalMatrix, ImageType::Pointer image);
	boost::multi_array<T, R> calculateConvolutionMatrix(ImageType::Pointer mask);
	void calculateLocalIntensityPeak(ImageType::Pointer peakMatrix, ImageType::Pointer image, ImageType::Pointer mask );
	void calculateGlobalIntensityPeak(ImageType::Pointer peakMatrix, ImageType::Pointer mask);

	void defineLocalIntenseFeatures(vector<string> &features);
	void defineLocalIntenseFeaturesOntology(vector<string> &features);
	void extractLocalIntenseData(vector<T> &localIntData, LocalIntensityFeatures<T, R> localIntFeatures);

public:
	LocalIntensityFeatures() {
	}
	~LocalIntensityFeatures() {
	}
	T localIntensityPeak = NAN;
	T globalIntensityPeak = NAN;
	void calculateAllLocalIntensityFeatures(LocalIntensityFeatures<T, R> &localInt, ImageType::Pointer image, ImageType::Pointer mask, ConfigFile config);
	void writeCSVFileLocalIntensity(LocalIntensityFeatures<T, R> localInt, string outputfolder);
	void writeOneFileLocalInt(LocalIntensityFeatures<T, R> localInt, ConfigFile config);
	void writeCSVFileLocalIntensityPET(LocalIntensityFeatures<T, R> localInt, ConfigFile config);
	void writeOneFileLocalIntPET(LocalIntensityFeatures<T, R> localInt, ConfigFile config);

};


/*!
In the function getIndexOfMax all indices of the voxels that yield the maximum intensity value are stored in a vector. \n
Herefore, a search is performed in the input matrix in order to look for the maximum element.\n
@parameter[in]: boost::multi_array inputMatrix
@parameter[out]: vector containing indices
*/
template<class T, size_t R>
vector<itk::Index<R> > LocalIntensityFeatures<T, R>::getIndexOfMax(ImageType::Pointer image, ImageType::Pointer mask) {
	itk::ImageRegionConstIterator<ImageType> maskIt(mask, mask->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<ImageType> imageIt(image, image->GetLargestPossibleRegion());
	maskIt.GoToBegin();
	imageIt.GoToBegin();
	float maxValue = 0;
	itk::Index<R> actIndex;
	PointType actCoordinates;
	int count = 0;
	vector<itk::Index<R> > allMaxIndices;
	while (!maskIt.IsAtEnd())
	{

		if (maskIt.Get() > 0 && imageIt.Get() > maxValue)
		{
			maxValue = imageIt.Get();
			actIndex = maskIt.GetIndex();

		}
		++maskIt;
		++imageIt;
	}
	maskIt.GoToBegin();
	imageIt.GoToBegin();
	while (!maskIt.IsAtEnd())
	{

		if (maskIt.Get() > 0 && imageIt.Get() == maxValue)
		{
			actIndex = maskIt.GetIndex();
			allMaxIndices.push_back(actIndex);

		}
		++maskIt;
		++imageIt;
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
	nrVoxelsDirection = floor((radius-(float(spacing)/2.0))  / spacing);
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

	//get the center index
	vector<float> indexOfCenter;
	
	fillVector(indexOfCenter, convMatrix);
	vector<float> actPos = { {0,0,0} };
	float distFromCenter;
	int nrVoxels = 0;
	int test = 0;
	for (int depth = 0; depth < convMatrix.shape()[2]; ++depth) {
		for (int row = 0; row < convMatrix.shape()[0]; ++row) {
			for (int col = 0; col < convMatrix.shape()[1]; ++col) {
				actPos[0] = float(row)*inputSpacing[0] + inputSpacing[0]/2 - (indexOfCenter[0] + 0.5)*inputSpacing[0] ;
				actPos[1] = float(col)*inputSpacing[1] + inputSpacing[1]/2 - (indexOfCenter[1] + 0.5)*inputSpacing[1] ;
				actPos[2] = float(depth)*inputSpacing[2] + inputSpacing[2]/2 - (indexOfCenter[2] + 0.5)*inputSpacing[2];
				distFromCenter = std::sqrt(std::pow((actPos[0] - (indexOfCenter[0]+0.5)*inputSpacing[0] ), 2) + std::pow((actPos[1]  - (indexOfCenter[1]+0.5) * inputSpacing[1] ), 2) + std::pow((actPos[2] - (indexOfCenter[2]+0.5) * inputSpacing[2] ), 2));
				if(abs(actPos[0]) <= 6.2 && abs(actPos[1]) <= 6.2 && abs(actPos[2])<=6.2){
					convMatrix[row][col][depth] = 1;
				}	
			}
		}
	}

}



//fill a vector with indices
template<class T, size_t R>
void LocalIntensityFeatures<T, R>::fillVector(vector<float> &index, boost::multi_array<T, R> convMatrix) {
	index.push_back(float(convMatrix.shape()[0] -1 )/ 2.0);
	index.push_back(float(convMatrix.shape()[1] -1)/ 2.0);
	index.push_back(float(convMatrix.shape()[2] -1)/ 2.0);
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
ImageType::Pointer LocalIntensityFeatures<T, R>::calculatePeaks(boost::multi_array<T, R> convolutionalMatrix, ImageType::Pointer image) {
//boost::multi_array<T, R> LocalIntensityFeatures<T, R>::calculatePeaks(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> localIntMatrix, boost::multi_array<T, R> convolutionalMatrix, ImageType::Pointer image) {
	
	//assign the image matrix
	const typename ImageType::RegionType& imageRegion = image->GetLargestPossibleRegion();
	const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
	int peakMatrixX = imageRegionSize[0];
	int peakMatrixY = imageRegionSize[1];
	int peakMatrixZ = imageRegionSize[2];

	const typename ImageType::SpacingType& inputSpacing = image->GetSpacing();
	unsigned int dimImage[] = { convolutionalMatrix.shape()[0],convolutionalMatrix.shape()[1],convolutionalMatrix.shape()[2] };
	float voxelSizeImage[] = { inputSpacing[0],inputSpacing[1],inputSpacing[2] };

	ImportFilterType::SizeType size;
	//set image size
	size[0] = convolutionalMatrix.shape()[0];
	size[1] = convolutionalMatrix.shape()[1];
	size[2] = convolutionalMatrix.shape()[2];

	int nrVoxelsPET = size[0] * size[1] * size[2];
	//start with importing the image
	ImageType::RegionType region;

	ImageType::IndexType start;
	start.Fill(0);

	region.SetIndex(start);



	region.SetSize(size);

	ImageType::Pointer imageNew = ImageType::New();
	imageNew->SetRegions(region);
	imageNew->Allocate();
	imageNew->FillBuffer(itk::NumericTraits<T>::Zero);

	ImageType::IndexType pixelIndex;
	region.SetIndex(start);
	region.SetSize(size);

	for (int depth = 0; depth < convolutionalMatrix.shape()[2] ; depth++) {
		for (int row = 0; row < convolutionalMatrix.shape()[0] ; row++) {
			for (int col = 0; col < convolutionalMatrix.shape()[1] ; col++) {
				pixelIndex[0] = row;
				pixelIndex[1] = col ;
				pixelIndex[2] = depth;
				imageNew->SetPixel(pixelIndex, convolutionalMatrix[row][col][depth]);
			}
		}
	}
	using FilterType = itk::ConvolutionImageFilter<ImageType>;
	FilterType::Pointer convolutionFilter = FilterType::New();
	convolutionFilter->SetNormalize(true);
	convolutionFilter->SetInput(image);
	convolutionFilter->SetKernelImage(imageNew);
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
void LocalIntensityFeatures<T, R>::calculateLocalIntensityPeak(ImageType::Pointer peakMatrix, ImageType::Pointer image, ImageType::Pointer mask) {
	vector<itk::Index<R> > allMaxIndices = getIndexOfMax(image, mask);
	T mean = 0;
	itk::Index<R> actualIndex;
	//iterate over all indices that have the intensity value equal to the maximum element
	//put all elements that are lying inside this circle on a vector
	for (int i = 0; i < allMaxIndices.size(); i++) {
		actualIndex = allMaxIndices[i];
		mean += peakMatrix->GetPixel(actualIndex);
	}
	localIntensityPeak = mean / allMaxIndices.size();
}

/*!
In the function calculateGlobalIntensityPeak the global intensity peak is calculated using the peak matrix. \n
@parameter[in]: boost multi_array peak matrix
*/
template<class T, size_t R>
void LocalIntensityFeatures<T, R>::calculateGlobalIntensityPeak(ImageType::Pointer peakMatrix, ImageType::Pointer mask) {
	itk::ImageRegionConstIterator<ImageType> maskIt(mask, mask->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<ImageType> imageIt(peakMatrix, peakMatrix->GetLargestPossibleRegion());
	maskIt.GoToBegin();
	imageIt.GoToBegin();
	float maxValue = 0;
	
	while (!maskIt.IsAtEnd())
	{
		if (maskIt.Get() > 0 && imageIt.Get() > maxValue)
		{
			maxValue = imageIt.Get();
		}
		++maskIt;
		++imageIt;
	}
	globalIntensityPeak = maxValue;

}

template<class T, size_t R>
void LocalIntensityFeatures<T, R>::calculateAllLocalIntensityFeatures(LocalIntensityFeatures<T, R> &localInt, ImageType::Pointer image, ImageType::Pointer mask, ConfigFile config) {
	const typename ImageType::SpacingType& inputSpacing = mask->GetSpacing();
	voxelSize[0] = inputSpacing[0];
	voxelSize[1] = inputSpacing[1];
	voxelSize[2] = inputSpacing[2];
	const typename ImageType::RegionType& imageRegion = mask->GetLargestPossibleRegion();
	const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();

	vector<itk::Index<R> > allMaxIndices = getIndexOfMax(image, mask);
	boost::multi_array<T, R> convMatrix = localInt.calculateConvolutionMatrix(image);
	//boost::multi_array<T, R> peakMatrix = calculatePeaks(imageAttr.imageMatrix, imageAttr.imageMatrixLocalInt, convMatrix, imageAttr.image);
	ImageType::Pointer peakImage = calculatePeaks(convMatrix, image);
	
	localInt.calculateLocalIntensityPeak(peakImage, image, mask);
	localInt.calculateGlobalIntensityPeak(peakImage, mask);
	localInt.localIntensityPeak = localInt.localIntensityPeak;
	localInt.globalIntensityPeak = localInt.globalIntensityPeak;
	if (config.imageType == "PET" && config.useSUV == 1) {
		float correctionParam;
		if (config.correctionParam != 0) {
			correctionParam = config.correctionParam;
		}
		else {
			correctionParam = config.patientWeight / (config.initActivity * 1000);
		}
		localInt.localIntensityPeak = localInt.localIntensityPeak*correctionParam;
		localInt.globalIntensityPeak = localInt.globalIntensityPeak*correctionParam;
	}
	peakImage = nullptr;
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
void LocalIntensityFeatures<T, R>::writeOneFileLocalInt(LocalIntensityFeatures<T, R> localInt, ConfigFile config) {
	string csvName;
	if (config.csvOutput == 1) {
		csvName = config.outputFolder + ".csv";
	}
	else if (config.ontologyOutput == 1) {
		csvName = config.outputFolder + "/feature_table.csv";
	}
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream localIntCSV;
	localIntCSV.open(name, std::ios_base::app);
	vector<string> features;
	

	vector<T> localIntData;
	extractLocalIntenseData(localIntData, localInt);
	if(config.getOneCSVFile == 1) {
		defineLocalIntenseFeatures(features);
		for (int i = 0; i < localIntData.size(); i++) {
			localIntCSV << "Local intensity" << "," << features[i] << ",";
			localIntCSV << localIntData[i];
			localIntCSV << "\n";
		}
		localIntCSV.close();
	}
	else if (config.ontologyOutput == 1) {
		defineLocalIntenseFeaturesOntology(features);
		for (int i = 0; i < localIntData.size(); i++) {
			localIntCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			localIntCSV << localIntData[i] << "," << config.featureParameterSpaceName << "," << config.calculationSpaceName;
			localIntCSV << "\n";
		}
		localIntCSV.close();
	}
}


template<class T, size_t R>
void LocalIntensityFeatures<T, R>::writeCSVFileLocalIntensityPET(LocalIntensityFeatures<T, R> localInt, ConfigFile config) {
	string csvName = config.outputFolder + "_localIntensity.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream localIntCSV;
	localIntCSV.open(name);
	vector<string> features;
	defineLocalIntenseFeatures(features);

	vector<T> localIntData;
	extractLocalIntenseData(localIntData, localInt);
	for (int i = 0; i < localIntData.size(); i++) {
		if (config.imageType == "PET") {
			localIntCSV << "PET Uptake Metrics" << "," << features[i] << ",";
		}
		else {
			localIntCSV << "Exact Metrics" << "," << features[i] << ",";
		}
		localIntCSV << localIntData[i];
		localIntCSV << "\n";
	}
	localIntCSV.close();
}

template <class T, size_t R>
void LocalIntensityFeatures<T, R>::writeOneFileLocalIntPET(LocalIntensityFeatures<T, R> localInt, ConfigFile config) {
	string csvName;
	if (config.csvOutput == 1) {
		csvName = config.outputFolder + ".csv";
	}
	else if (config.ontologyOutput == 1) {
		csvName = config.outputFolder + "/feature_table.csv";
	}
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream localIntCSV;
	localIntCSV.open(name, std::ios_base::app);
	
	vector<string> features;


	vector<T> localIntData;
	extractLocalIntenseData(localIntData, localInt);
	if (config.getOneCSVFile == 1) {
		defineLocalIntenseFeatures(features);
		for (int i = 0; i < localIntData.size(); i++) {
			if (config.imageType == "PET") {
				localIntCSV << "PET Uptake Metrics" << "," << features[i] << ",";
			}
			else {
				localIntCSV << "Exact Metrics" << "," << features[i] << ",";
			}
			localIntCSV << localIntData[i];
			localIntCSV << "\n";
		}
		localIntCSV.close();
	}
	else if (config.ontologyOutput == 1) {
		defineLocalIntenseFeaturesOntology(features);
		for (int i = 0; i < localIntData.size(); i++) {
			localIntCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			localIntCSV << localIntData[i] << "," << config.featureParameterSpaceName << "," << config.calculationSpaceName;
			localIntCSV << "\n";
		}
		localIntCSV.close();
	}
}


template <class T, size_t R>
void LocalIntensityFeatures<T, R>::defineLocalIntenseFeatures(vector<string> &features) {
	features.push_back("local intensity peak");
	features.push_back("global intensity peak");

}

template <class T, size_t R>
void LocalIntensityFeatures<T, R>::defineLocalIntenseFeaturesOntology(vector<string> &features) {
	features.push_back("Floc.peak.local");
	features.push_back("Floc.peak.global");

}

template <class T, size_t R>
void LocalIntensityFeatures<T, R>::extractLocalIntenseData(vector<T> &localIntData, LocalIntensityFeatures<T, R> localIntFeatures) {
	localIntData.push_back(localIntFeatures.localIntensityPeak);
	localIntData.push_back(localIntFeatures.globalIntensityPeak);
}

#endif
