#ifndef DISPERSITYFEATURES_H_INCLUDED
#define DISPERSITYFEATURES_H_INCLUDED

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
class DispersityFeatures {
private:
	using PointType = itk::Point<T, R>;
	typedef itk::ConnectedComponentImageFilter <intImage, intImage> ConnectedComponentFilterType;

	typedef itk::LabelImageToShapeLabelMapFilter<intImage> LabelImageToShapeLabelMapFilterType;
	typedef int LabelType;
	typedef itk::ShapeLabelObject<LabelType, R> ShapeLabelObjectType;
	typedef itk::LabelMap<ShapeLabelObjectType> LabelMapType;
	float correctionParam = 0;
	float DmaxBulk = NAN;
	float DmaxPatient = NAN;
	int nrLesions = NAN;
	float spreadBulk = NAN;
	float spreadPatient = NAN;
	float volMaxPatient = NAN;
	float volSpreadPatient = NAN;
	float volSpreadBulk = NAN;
	float maxSpreadPatient = NAN;
	float maxSpreadBulk = NAN;
	float maxSpreadHot = NAN;
	float maxDistPatient = NAN;
	float maxDistBulk = NAN;
	float peakSpreadPatient = NAN;
	float peakSpreadBulk = NAN;
	float peakSpreadHot = NAN;
	float peakDistPatient = NAN;
	float peakDistBulk = NAN;
	float ratioPeakBulk = NAN;
	float ratioPeakPatient = NAN;
	float ratioVolPatient = NAN;
	float voxelSize[3];
	void getLabelObjectFeatures(ImageType::Pointer image, ImageType::Pointer mask, ConfigFile config);
	itk::Point<T, R> getMaxIndices(ImageType::Pointer image, ImageType::Pointer mask, int labelNr, vector<float> &maxValues);
	void getBulkFeatures(vector<int> volume, vector<PointType> indices);
	void getTotalFeatures(vector<int> volume, vector<PointType> indices);
	float getEuclideanDist(PointType p1, PointType p2);
	void getVolumeDispersity(vector<int> volume);
	void getMaxDispersity(vector<T> maxValues, vector<int> volume);
	float getPEAKvalues(ImageType::Pointer mask, ImageType::Pointer image, ConfigFile config);
	void defineDispersityFeatures(vector<string> &features);
	ImageType::Pointer setMaskOneLesion(ImageType::Pointer mask, int labelNr);
	void extractDispersityData(vector<T> &dispData, DispersityFeatures<T, R> dispFeatures);
	void getPeakDispersityFeatures(vector<int> volume, vector<float> peakVector);
public:
	void writeOneFileDispersity(DispersityFeatures<T, R> disp, ConfigFile config);
	void writeCSVFileDispersity(DispersityFeatures<T, R> disp, string outputFolder, ConfigFile config);
	void calculateAllDispersityFeatures(DispersityFeatures<T, R> &dispFeatures, ImageType::Pointer image, ImageType::Pointer mask, ConfigFile config);
};
/*!
In the function getLabelObjectFeatures a set of morphological features is calculated using
the label image to shape label map filter
*/
template<class T, size_t R>
itk::Point<T, R> DispersityFeatures<T, R>::getMaxIndices(ImageType::Pointer image, ImageType::Pointer mask, int labelNr, vector<float> &maxValues) {


	itk::ImageRegionConstIterator<ImageType> maskIt(mask, mask->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator<ImageType> imageIt(image, image->GetLargestPossibleRegion());
	maskIt.GoToBegin();
	imageIt.GoToBegin();
	float maxValue = 0;
	itk::Index<R> actIndex;
	PointType actCoordinates;
	
	while (!maskIt.IsAtEnd()){
	
		if (maskIt.Get() == labelNr && imageIt.Get()* correctionParam > maxValue)
		{
			maxValue = imageIt.Get() * correctionParam;
			actIndex = maskIt.GetIndex();

		}

		++maskIt;
		++imageIt;
	}
	if (maxValue == 0) {
		maskIt.GoToBegin();
		imageIt.GoToBegin();
		while (!maskIt.IsAtEnd())
		{
			if (maskIt.Get() > 0)
			{
				maxValue = maskIt.Get() * correctionParam;
				actIndex = maskIt.GetIndex();
			}
			++maskIt;

		}

	}
	const itk::Index<R> test = actIndex;
	image->TransformIndexToPhysicalPoint(test, actCoordinates);
	maxValues.push_back(maxValue);
	return actCoordinates;
}


//setMaskOneLesion: Get a mask that contains only the mask of one lesion (set all other mask points to 0)
//this method is only applied if there is more than one lesion present in the mask
template<class T, size_t R>
ImageType::Pointer DispersityFeatures<T, R>::setMaskOneLesion(ImageType::Pointer mask, int labelNr) {
	ImageType::Pointer copyMask = ImageType::New();
	copyMask->SetRegions(mask->GetLargestPossibleRegion());
	copyMask->Allocate();
	itk::ImageRegionConstIterator<ImageType> maskIt(mask, mask->GetLargestPossibleRegion());
	itk::ImageRegionIterator<ImageType> copyIt(copyMask, copyMask->GetLargestPossibleRegion());
	maskIt.GoToBegin();
	copyIt.GoToBegin();

	while (!maskIt.IsAtEnd())
	{
		if (maskIt.Get() == labelNr) {
			copyIt.Set(100);
		}
		else {
			copyIt.Set(0);
		}
		++maskIt;
		++copyIt;
	}
	return copyMask;
}

//calculate the euclidean distance of two ITK points
template<class T, size_t R>
float DispersityFeatures<T, R>::getEuclideanDist(PointType p1, PointType p2) {
	float dist = std::sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p2[1])*(p1[1] - p2[1]) + (p1[2] - p2[2])*(p1[2] - p2[2]));
	
	return dist;
}

//get global SUV peak of each lesion in order to calculate the 'Peak' dispersity features
template<class T, size_t R>
float DispersityFeatures<T, R>::getPEAKvalues(ImageType::Pointer mask, ImageType::Pointer image, ConfigFile config) {
	LocalIntensityFeatures<float, 3> localInt;
	RegionType boundingBoxRegion = getBoundingBoxMask(mask);
	mask = getImageMasked(mask, boundingBoxRegion);
	image = getImageMasked(image, boundingBoxRegion);
	localInt.calculateAllLocalIntensityFeatures(localInt, image, mask, config);
	return localInt.globalIntensityPeak;
}


//get distance of other lesions to largest lesion in the body
template<class T, size_t R>
void DispersityFeatures<T, R>::getBulkFeatures(vector<int> volume, vector<PointType> indices) {
	int maxVolume = 0;
	PointType indexMax;
	//get index of biggest tumor
	for (int i = 0; i < volume.size(); i++) {
		if (volume[i] > maxVolume) {
			maxVolume = volume[i];
			indexMax = indices[i];
		}
	}

	//calculate maxDistance and spread
	spreadBulk = 0;
	DmaxBulk = 0;
	float actDist = 0;
	for (int i = 0; i < volume.size(); i++) {
		if (volume[i] != maxVolume) {
			actDist = getEuclideanDist(indexMax, indices[i]);
			if (actDist > DmaxBulk) {
				DmaxBulk = actDist;
			}
			spreadBulk += actDist;
		}
	}


}

template<class T, size_t R>
void DispersityFeatures<T, R>::getPeakDispersityFeatures(vector<int> volume, vector<float> peakVector) {
	peakSpreadPatient = 0;
	peakSpreadBulk = 0;
	peakDistPatient = 0;
	peakDistBulk = 0;
	int maxVol = 0;
	int actVolume = 0;
	int bulkIndex = 0;
	for (int i = 0; i < boost::size(volume); i++) {
		actVolume = volume[i];
		if (volume[i] > maxVol) {
			maxVol = volume[i];
			bulkIndex = i;
		}
	}

	float maxBulk = peakVector[bulkIndex];
	float totalMax = maxBulk;
	float minValue = maxBulk;
	for (int i = 0; i < boost::size(peakVector); i++) {
		if (i != bulkIndex) {
			peakSpreadBulk += abs(peakVector[i] - maxBulk);
		}
		if (peakVector[i] < minValue && peakVector[i]>0) {
			minValue = peakVector[i];
		}
		if (peakVector[i] > maxBulk) {
			totalMax = peakVector[i];
		}
	}
	peakDistBulk = maxBulk - minValue;
	peakDistPatient = totalMax - minValue;
	ratioPeakBulk = maxBulk/minValue;
	ratioPeakPatient = totalMax/minValue;
	peakSpreadPatient = 0;
	peakSpreadHot = 0;
	for (int i = 0; i < boost::size(peakVector); i++) {
		for (int j = i; j < boost::size(peakVector); j++) {
			peakSpreadPatient += abs(peakVector[i] - peakVector[j]);
		}
		if (peakVector[i] != totalMax) {
			peakSpreadHot += totalMax - peakVector[i];
		}
	}
}

//get features describing the total tumor spread
template<class T, size_t R>
void DispersityFeatures<T, R>::getTotalFeatures(vector<int> volume, vector<PointType> indices) {

	//calculate maxDistance and spread
	spreadPatient = 0;
	DmaxPatient = 0;
	float actDist = 0;
	for (int i = 0; i < boost::size(volume); i++) {
		for (int j = i; j < boost::size(volume); j++) {
			actDist = getEuclideanDist(indices[i], indices[j]);
			if (actDist > DmaxPatient) {
				DmaxPatient = actDist;
			}
			spreadPatient += actDist;
		}
	}
}

//get the dispersity features desribing volume differences in patient
template<class T, size_t R>
void DispersityFeatures<T, R>::getVolumeDispersity(vector<int> volume) {
	int maxVolume = 0;
	//get index of biggest tumor
	for (int i = 0; i < volume.size(); i++) {
		if (volume[i] > maxVolume) {
			maxVolume = volume[i];
		}
	}
	int minVolume = maxVolume;
	for (int i = 0; i < boost::size(volume); i++) {
		if (volume[i] < minVolume) {
			minVolume = volume[i];
		}
	}
	volMaxPatient = maxVolume - minVolume;
	ratioVolPatient = maxVolume / minVolume;
	
	volSpreadBulk = 0;
	int actVolume = 0;
	for (int i = 0; i < boost::size(volume); i++) {
		actVolume = volume[i];
		if (actVolume < maxVolume) {
			volSpreadBulk += maxVolume - actVolume;
		}
	}
	volSpreadPatient = 0;
	for (int i = 0; i < boost::size(volume); i++) {
		for (int j = i; j < boost::size(volume); j++) {
			actVolume = volume[i];
			volSpreadPatient += abs(volume[i] - volume[j]);
		}
	}
}

//get the dispersity features desribing SUV maximum differences in patient
template<class T, size_t R>
void DispersityFeatures<T, R>::getMaxDispersity(vector<T> maxValues, vector<int> volume) {
	maxSpreadPatient = 0;
	maxSpreadBulk = 0;
	maxDistPatient = 0;
	maxDistBulk = 0;
	int maxVol = 0;
	int actVolume = 0;
	int bulkIndex = 0;
	for (int i = 0; i < boost::size(volume); i++) {
		actVolume = volume[i];
		if (volume[i] > maxVol) {
			maxVol = volume[i];
			bulkIndex = i;
		}
	}

	float maxBulk = maxValues[bulkIndex];
	float totalMax = maxBulk;
	float minValue = maxBulk;
	for (int i = 0; i < boost::size(maxValues); i++) {
		if (i != bulkIndex) {
			maxSpreadBulk += abs(maxValues[i] - maxBulk);
		}
		if (maxValues[i] < minValue && maxValues[i]>0) {
			minValue = maxValues[i];
		}
		if (maxValues[i] > maxBulk) {
			totalMax = maxValues[i];
		}
	}
	maxDistBulk = maxBulk - minValue;
	maxDistPatient = totalMax - minValue;

	maxSpreadPatient = 0;
	maxSpreadHot = 0;
	for (int i = 0; i < boost::size(maxValues); i++) {
		for (int j = i; j < boost::size(maxValues); j++) {
			maxSpreadPatient += abs(maxValues[i] - maxValues[j]);	
		}
		if (maxValues[i] != totalMax) {
			maxSpreadHot += totalMax - maxValues[i];
		}
	}
}

//calculate all dispersity features by first:
//labelling the mask such that each lesion gets a different mask value
//i.e. lesion 1 is labeled with 1 in the mask, lesion 2 with 2 etc. 
//this 'labeled' mask is then further used to calculate the different types of dispersity features
template<class T, size_t R>
void DispersityFeatures<T, R>::getLabelObjectFeatures(ImageType::Pointer image, ImageType::Pointer mask, ConfigFile config) {


	//With the cast filter type function, the mask is converted to an int-image
	//because otherwise the Label image to shape label map filter is not working
	typename CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(mask);
	//With the ConnectedComponentFilter from the itk library the connected components of the mask are calculated
	typename ConnectedComponentFilterType::Pointer connectedComponentImageFilter = ConnectedComponentFilterType::New();
	connectedComponentImageFilter->SetFullyConnected(true);
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
	using LabelMapToLabelImageFilterType = itk::LabelMapToLabelImageFilter<LabelMapType, ImageType>;
	LabelMapToLabelImageFilterType::Pointer labelImageConverter = LabelMapToLabelImageFilterType::New();
	labelImageConverter->SetInput(labelMap);
	labelImageConverter->Update();
	ImageType::Pointer labeledMask = labelImageConverter->GetOutput();
	//For every connected component a labelObject is created
	//Because we are only interested in the object with the label one, we check the label number
	DmaxBulk = 0;
	DmaxPatient = 0;
	nrLesions = 1;
	spreadBulk = 0;
	spreadPatient = 0;
	volMaxPatient = 0;
	volSpreadPatient = 0;
	volSpreadBulk = 0;
	maxSpreadPatient = 0;
	maxSpreadBulk = 0;
	maxSpreadHot = 0;
	maxDistPatient = 0;
	maxDistBulk = 0;
	peakSpreadPatient = 0;
	peakSpreadBulk = 0;
	peakSpreadHot = 0;
	peakDistPatient = 0;
	peakDistBulk = 0;
	ratioPeakBulk = 0;
	ratioPeakPatient = 0;
	ratioVolPatient = 0;
	vector<int> volumes;
	vector< PointType> actIndices;
	vector<float> peakVector;
	int nrPixels;
	int nrObjects = labelMap->GetNumberOfLabelObjects();
	if (nrObjects > 1) {
		//-1 because background is also an object
		nrLesions = nrObjects;
		vector<T> maxValues;
		for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); n++) {
			ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
			int labelNr = labelObject->GetLabel();

			//if (labelNr > 0) {
			nrPixels = labelObject->GetNumberOfPixels();
			volumes.push_back(nrPixels);

			actIndices.push_back(getMaxIndices(image, labeledMask, labelNr, maxValues));
			ImageType::Pointer actMask = setMaskOneLesion(labeledMask, labelNr);
			peakVector.push_back(getPEAKvalues(actMask, image, config));
			//}
		}
		getBulkFeatures(volumes, actIndices);
		getTotalFeatures(volumes, actIndices);
		getMaxDispersity(maxValues, volumes);
		getVolumeDispersity(volumes);
		getPeakDispersityFeatures(volumes, peakVector);


	}
	//delete object from memory
	labeledMask = nullptr;
	labelMap = nullptr;
}



template <class T, size_t R>
void DispersityFeatures<T, R>::calculateAllDispersityFeatures(DispersityFeatures<T, R> &dispFeatures, ImageType::Pointer image, ImageType::Pointer mask, ConfigFile config) {
	const typename ImageType::SpacingType& inputSpacing = image->GetSpacing();
	voxelSize[0] = inputSpacing[0];
	voxelSize[1] = inputSpacing[1];
	voxelSize[2] = inputSpacing[2];
	if (config.correctionParam != 0) {
		correctionParam = config.correctionParam;
	}
	else {
		correctionParam = config.patientWeight / (config.initActivity * 1000);
	}
	getLabelObjectFeatures(image, mask, config);
	volMaxPatient = volMaxPatient*voxelSize[0]* voxelSize[1]* voxelSize[2];
	volSpreadPatient = volSpreadPatient * voxelSize[0] * voxelSize[1] * voxelSize[2];
	volSpreadBulk = volSpreadBulk * voxelSize[0] * voxelSize[1] * voxelSize[2];
}

template <class T, size_t R>
void DispersityFeatures<T, R>::writeCSVFileDispersity(DispersityFeatures<T, R> disp, string outputFolder, ConfigFile config)
{
	string csvName = outputFolder + "_dispersityFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream dispCSV;

	dispCSV.open(name, std::ios_base::app);


	vector<string> features;
	defineDispersityFeatures(features);

	vector<T> dispData;
	extractDispersityData(dispData, disp);
	for (int i = 0; i < dispData.size(); i++) {
		dispCSV << "Dispersity" << "," << features[i] << ",";
		dispCSV << dispData[i];
		dispCSV << "\n";
	}
	dispCSV.close();
}


template <class T, size_t R>
void DispersityFeatures<T, R>::writeOneFileDispersity(DispersityFeatures<T, R> disp, ConfigFile config) {
	string csvName;
	if (config.getOneCSVFile == 1) {
		csvName = config.outputFolder + ".csv";
	}


	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';
	ofstream dispCSV;
	dispCSV.open(name, std::ios_base::app);

	vector<string> features;
	vector<T> dispData;
	extractDispersityData(dispData, disp);
	if (config.getOneCSVFile == 1) {
		defineDispersityFeatures(features);
		for (int i = 0; i < dispData.size(); i++) {
			dispCSV << "Dispersity" << "," << features[i] << ",";
			dispCSV << dispData[i];
			dispCSV << "\n";
		}
		dispCSV.close();
	}

}
template <class T, size_t R>
void DispersityFeatures<T, R>::defineDispersityFeatures(vector<string> &features) {
	features.push_back("NumberLesions");
	features.push_back("DmaxBulk");
	features.push_back("SpreadBulk");
	features.push_back("DmaxPatient");
	features.push_back("SpreadPatient");
	features.push_back("VolSpreadBulk");
	//features.push_back("ratioPeakBulk");
	//features.push_back("ratioPeakPatient");
	//features.push_back("ratioVolPatient");
	features.push_back("DvolPatient");
	features.push_back("VolSpreadPatient");
	features.push_back("DSUVmaxBulk");
	features.push_back("DSUVmaxSumBulk");	
	features.push_back("DSUVmaxPatient");
	features.push_back("DSUVmaxSumPatient");
	features.push_back("DSUVmaxSumHot");
	features.push_back("DSUVpeakBulk");
	features.push_back("DSUVpeakSumBulk");
	features.push_back("DSUVpeakPatient");
	features.push_back("DSUVpeakSumPatient");
	features.push_back("DSUVpeakSumHot");
	
}


template <class T, size_t R>
void DispersityFeatures<T, R>::extractDispersityData(vector<T> &dispData, DispersityFeatures<T, R> dispFeatures) {
	dispData.push_back(dispFeatures.nrLesions);
	dispData.push_back(dispFeatures.DmaxBulk);
	dispData.push_back(dispFeatures.spreadBulk);
	dispData.push_back(dispFeatures.DmaxPatient);
	dispData.push_back(dispFeatures.spreadPatient);
	dispData.push_back(dispFeatures.volSpreadBulk);
	//dispData.push_back(dispFeatures.ratioPeakBulk);
	//dispData.push_back(dispFeatures.ratioPeakPatient);
	//dispData.push_back(dispFeatures.ratioVolPatient);
	dispData.push_back(dispFeatures.volMaxPatient);
	dispData.push_back(dispFeatures.volSpreadPatient);
	dispData.push_back(dispFeatures.maxDistBulk);
	dispData.push_back(dispFeatures.maxSpreadBulk);
	dispData.push_back(dispFeatures.maxDistPatient);
	dispData.push_back(dispFeatures.maxSpreadPatient);
	dispData.push_back(dispFeatures.maxSpreadHot);
	dispData.push_back(dispFeatures.peakDistBulk);
	dispData.push_back(dispFeatures.peakSpreadBulk);
	dispData.push_back(dispFeatures.peakDistPatient);
	dispData.push_back(dispFeatures.peakSpreadPatient);
	dispData.push_back(dispFeatures.peakSpreadHot);
	
}
#endif // DISPERSITYFEATURES_H_INCLUDED
