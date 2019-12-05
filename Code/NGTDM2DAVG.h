#ifndef NGTDM2DAVG_H_INCLUDED
#define NGTDM2DAVG_H_INCLUDED

#include "NGTDM2DMRG.h"

/*! \file */
/*!
The class NGTDM inherites from the class NGTDM. The feature calculations are done according to the definitions in this class. \n
These matrices combine the sum of grey level differences of voxels with intensity value i and the average
discretised grey levels of a neighborhood with distance dist from the actual voxel. \n
The average grey level within a neighborhood is defined as:
\f$ A_{i} = \frac{1}{W} \sum_{k_{x} = -dist}^{dist} \sum_{k_{y} = -dist}^{dist} \sum_{k_{z} = -dist}^{dist} X_{dgl}(j_{x} + k_{x}, j_{y} + k_{y}, j_{z} + k_{z}) \f$ \n
where \f$ k_{x}, k_{y},k_{z} !=0 \f$ and \f$ W = (2dist+1)\f$\n
Now, let \f$ n_{i} \f$ be the number of voxels with grey level i that have a complete neighborhood. \n
The entry in the NGTDM matrix is then: \f$ s_{i} = \sum^{n_{i}} i-A_{i}, if n_{i} >0 \f$ and 0 otherwise.\n
*/


template <class T, size_t R>
class NGTDM2DAVG : NGTDMFeatures2DMRG<T, R> {

private:
	NGTDMFeatures2DMRG<T, R> ngtdm;
	//these are the values to determine if the user wants to calculate novel/uncommon features
	vector<float> actualSpacing;
	string normNGTDM;
	//distance of neighborhood defined by the user
	int dist;

	void extractNGTDMData2DAVG(vector<T> &ngtdmData, NGTDM2DAVG<T, R> ngtdmFeatures);


	boost::multi_array<float, 2> getNGTDMatrix2DAVG(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> sumNeighborHoods, int depth, ConfigFile config);
public:
	void getProbability(vector<T> elementsOfWholeNeighborhood, boost::multi_array<float, 2> &ngtdMatrix);
	void calculateAllNGTDMFeatures2DAVG(NGTDM2DAVG<T, R> &ngtdmFeatures, Image<T, R> imageAttr, boost::multi_array<T, R> sumNeighborHoods, vector<float> spacing, ConfigFile config);
	void writeCSVFileNGTDM2DAVG(NGTDM2DAVG<T, R> ngtdm, string outputFolder);
	void writeOneFileNGTDM2DAVG(NGTDM2DAVG<T, R> ngtdm, ConfigFile config, int &parameterSpaceNr);
};



/*!
\brief getNGTDMatrix2DAVG
In this function the NGTDM is filled.
@param[in] boost multi array inputMatrix: matrix filled with intensity values
@param[in] int dist: size of neighborhood
@param[in] int depth: actual number of slice
@param[out] boost multi array: filled NGTD matrix
The function fills the NGTDMatrix with the corresponding values
*/
template <class T, size_t R>
boost::multi_array<float, 2> NGTDM2DAVG<T, R>::getNGTDMatrix2DAVG(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> sumNeighborHoods, int depth, ConfigFile config) {
	typedef boost::multi_array<float, 2>  ngtdmat;
	int sizeMatrix = this->diffGreyLevels.size();
	ngtdmat NGTDMatrix(boost::extents[sizeMatrix][3]);
	int indexOfElement[3] = { 0,0,0 };
	std::vector<T> elementsOfWholeNeighborhoods;
	T sumOfActualNeighborhood;
	T s = 0;
	T actualElement;
	int posActualElement;
	for (int row = 0; row < inputMatrix.shape()[0]; row++) {
		for (int col = 0; col < inputMatrix.shape()[1]; col++) {
			indexOfElement[0] = row;
			indexOfElement[1] = col;
			indexOfElement[2] = depth;
			//get actual Element if it is the centre of a whole neighborhood
			actualElement = inputMatrix[row][col][depth];
			if (!std::isnan(actualElement)) {
				//get the sum of the actual neighborhood
				sumOfActualNeighborhood = sumNeighborHoods[row][col][depth];
				//get the s_i value
				s = abs(actualElement - sumOfActualNeighborhood);
				//get the position of the actual Element in the diffGreyLevel vector and so also in the matrix
				posActualElement = std::find(this->diffGreyLevels.begin(), this->diffGreyLevels.end(), actualElement) - this->diffGreyLevels.begin();
				//add the s_i value to the right element in the matrix
				NGTDMatrix[posActualElement][2] += s;
				NGTDMatrix[posActualElement][0] += 1;
				//save the actual Element in a vector, so we can calculate later the probabilities
				elementsOfWholeNeighborhoods.push_back(actualElement);
			}
		}
	}
	getProbability(elementsOfWholeNeighborhoods, NGTDMatrix);
	elementsOfWholeNeighborhoods.clear();
	return NGTDMatrix;
}



/*
\brief getProbability
@param[in] vector<T> elementsOfWholeNeighborhood: vector containing all elements of the neighborhood
@param[in] boost::multi_array<float, 2> ngtdMatrix: ngtd matrix

The function calculates the probability matrix of the NGTD matrix, as many features are calculated using the probabilities.
*/
template <class T, size_t R>
void NGTDM2DAVG<T, R>::getProbability(vector<T> elementsOfWholeNeighborhood, boost::multi_array<float, 2> &ngtdMatrix) {
	float numItem = 0;
	for (int actElementIndex = 0; actElementIndex<boost::size(this->diffGreyLevels); actElementIndex++) {
		numItem += ngtdMatrix[actElementIndex][0];
	}
	for (int actElementIndex = 0; actElementIndex<boost::size(this->diffGreyLevels); actElementIndex++) {
		ngtdMatrix[actElementIndex][1] = ngtdMatrix[actElementIndex][0] / numItem;
	}

}




template <class T, size_t R>
void NGTDM2DAVG<T, R>::calculateAllNGTDMFeatures2DAVG(NGTDM2DAVG<T, R> &ngtdmFeatures, Image<T, R> imageAttr, boost::multi_array<T,R> sumNeighborHoods, vector<float> spacing, ConfigFile config) {
	this->diffGreyLevels = imageAttr.diffGreyLevels;
	//fill these values with the values set by user
	actualSpacing = spacing;
	normNGTDM = config.normNGTDM;
	dist = config.dist;

	float sumCoarseness = 0;
	float sumContrast = 0;
	float sumBusyness = 0;
	float sumComplexity = 0;
	float sumStrength = 0;

	int totalDepth = imageAttr.imageMatrix.shape()[2];
	for (int depth = 0; depth < totalDepth; depth++) {
		boost::multi_array<float, 2> ngtdm = ngtdmFeatures.getNGTDMatrix2DAVG(imageAttr.imageMatrix, sumNeighborHoods, depth, config);
		ngtdmFeatures.calculateCoarseness(ngtdm);

		sumCoarseness += this->coarseness;
		ngtdmFeatures.calculateContrast(ngtdm);
		sumContrast += this->contrast;
		ngtdmFeatures.calculateBusyness(ngtdm);
		sumBusyness += this->busyness;
		ngtdmFeatures.calculateStrength(ngtdm);
		sumStrength += this->strength;
		ngtdmFeatures.calculateComplexity(ngtdm);
		sumComplexity += this->complexity;
	}

	this->coarseness = sumCoarseness / totalDepth;
	this->contrast = sumContrast / totalDepth;
	this->busyness = sumBusyness / totalDepth;
	this->strength = sumStrength / totalDepth;
	this->complexity = sumComplexity / totalDepth;

}

template <class T, size_t R>
void NGTDM2DAVG<T, R>::writeCSVFileNGTDM2DAVG(NGTDM2DAVG<T, R> ngtdmFeatures, string outputFolder)
{
	string csvName = outputFolder + "_ngtdmFeatures2avg.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream ngtdmCSV;
	ngtdmCSV.open(name);
	vector<string> features;
	ngtdm.defineNGTDMFeatures2DMRG(features);

	vector<T> ngtdmData;
	extractNGTDMData2DAVG(ngtdmData, ngtdmFeatures);
	for (int i = 0; i< ngtdmData.size(); i++) {
		ngtdmCSV << "ngtdmFeatures2avg" << "," << features[i] << ",";
		ngtdmCSV << ngtdmData[i];
		ngtdmCSV << "\n";
	}
	ngtdmCSV.close();
}

template <class T, size_t R>
void NGTDM2DAVG<T, R>::writeOneFileNGTDM2DAVG(NGTDM2DAVG<T, R> ngtdmFeatures, ConfigFile config, int &parameterSpaceNr) {
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

	ofstream ngtdmCSV;
	ngtdmCSV.open(name, std::ios_base::app);
	vector<string> features;

	vector<T> ngtdmData;
	extractNGTDMData2DAVG(ngtdmData, ngtdmFeatures);
	
	if (config.csvOutput == 1) {
		ngtdm.defineNGTDMFeatures2DMRG(features);
		for (int i = 0; i < ngtdmData.size(); i++) {
			ngtdmCSV << "ngtdmFeatures2avg" << "," << features[i] << ",";
			ngtdmCSV << ngtdmData[i];
			ngtdmCSV << "\n";
		}
	}
	else if (config.ontologyOutput == 1) {
		ngtdm.defineNGTDMFeatures2DMRGOntology(features);
		string featParamSpaceTable = config.outputFolder + "/FeatureParameterSpace_table.csv";
		char * featParamSpaceTableName = new char[featParamSpaceTable.size() + 1];
		std::copy(featParamSpaceTable.begin(), featParamSpaceTable.end(), featParamSpaceTableName);
		featParamSpaceTableName[featParamSpaceTable.size()] = '\0';

		ofstream featSpaceTable;
		parameterSpaceNr += 1;
		string parameterSpaceName = "FeatureParameterSpace_" + std::to_string(parameterSpaceNr);
		featSpaceTable.open(featParamSpaceTableName, std::ios_base::app);
		featSpaceTable << parameterSpaceName << "," << "2DVmrg" << "," << config.imageSpaceName << "," << config.interpolationMethod << "\n";
		featSpaceTable.close();

		for (int i = 0; i < ngtdmData.size(); i++) {
			ngtdmCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			ngtdmCSV << ngtdmData[i] << "," << parameterSpaceName << "," << config.calculationSpaceName;
			ngtdmCSV << "\n";
		}

	}
	ngtdmCSV.close();
}


template <class T, size_t R>

void NGTDM2DAVG<T, R>::extractNGTDMData2DAVG(vector<T> &ngtdmData, NGTDM2DAVG<T, R> ngtdmFeatures) {
	ngtdmData.push_back(ngtdmFeatures.coarseness);
	ngtdmData.push_back(ngtdmFeatures.contrast);
	ngtdmData.push_back(ngtdmFeatures.busyness);
	ngtdmData.push_back(ngtdmFeatures.complexity);
	ngtdmData.push_back(ngtdmFeatures.strength);

}

#endif // NGTDM2DAVG_H_INCLUDED