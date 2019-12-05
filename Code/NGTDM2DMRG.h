#ifndef NGTDM2DMRG_H_INCLUDED
#define NGTDM2DMRG_H_INCLUDED

#include <iostream>
#include "boost/multi_array.hpp"
#include "image.h"

using namespace std;

/*! \file */
/*!
The class NGTDM is the class of the Neighborhood Grey Tone Difference Matrices. \n
These matrices combine the sum of grey level differences of voxels with intensity value i and the average
discretised grey levels of a neighborhood with distance dist from the actual voxel. \n
The average grey level within a neighborhood is defined as:
\f$ A_{i} = \frac{1}{W} \sum_{k_{x} = -dist}^{dist} \sum_{k_{y} = -dist}^{dist} \sum_{k_{z} = -dist}^{dist} X_{dgl}(j_{x} + k_{x}, j_{y} + k_{y}, j_{z} + k_{z}) \f$ \n
where \f$ k_{x}, k_{y},k_{z} !=0 \f$ and \f$ W = (2dist+1)\f$\n
Now, let \f$ n_{i} \f$ be the number of voxels with grey level i that have a complete neighborhood. \n
The entry in the NGTDM matrix is then: \f$ s_{i} = \sum^{n_{i}} i-A_{i}, if n_{i} >0 \f$ and 0 otherwise.\n
*/

template <class T, size_t R>
class NGTDMFeatures2DMRG {

private:
	//values for weighting the entries (novel and uncommon features)
	vector<float> actualSpacing;
	string normNGTDM;
	int dist;
	boost::multi_array<float, 2> getNGTDMatrix(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> neighborHoodSum);

	void extractNGTDMData(vector<T> &ngtdmData, NGTDMFeatures2DMRG<T, R> NGTDMFeatures2DMRG);

public:
	vector<T> diffGreyLevels;
	double coarseness;
	double contrast;
	double busyness;
	double complexity;
	double strength;

	void getProbability(vector<T> elementsOfWholeNeighborhood, boost::multi_array<float, 2> &ngtdMatrix);
	double calculateSumSi(boost::multi_array<float, 2> ngtdm);
	double calculateSumSiPi(boost::multi_array<float, 2> ngtdm);
	int getNGP(boost::multi_array<float, 2> ngtdm);
	int getNV(boost::multi_array<float, 2> ngtdm);

	void calculateStrength(boost::multi_array<float, 2> ngtdm);
	void calculateComplexity(boost::multi_array<float, 2> ngtdm);
	void calculateCoarseness(boost::multi_array<float, 2> ngtdm);
	void calculateContrast(boost::multi_array<float, 2> ngtdm);
	void calculateBusyness(boost::multi_array<float, 2> ngtdm);
	void calculateAllNGTDMFeatures2DMRG(NGTDMFeatures2DMRG<T, R> &ngtdm, Image<T, R> imageAttr, boost::multi_array<T,R> neighborHoodSum, vector<float> spacing, ConfigFile config);
	void writeCSVFileNGTDM(NGTDMFeatures2DMRG<T, R> ngtdm, string outputFolder);
	void writeOneFileNGTDM(NGTDMFeatures2DMRG<T, R> ngtdm, ConfigFile config, int &parameterSpaceNr);
	void defineNGTDMFeatures2DMRG(vector<string> &features);
	void defineNGTDMFeatures2DMRGOntology(vector<string> &features);
};

/*!
\brief getNGTDMatrix
@param boost multi array inputMatrix: matrix filled with intensity values
@param int dist: size of neighborhood
@param[out] boost multi array: filled NGTD matrix

The function fills the NGTDMatrix with the corresponding values
*/
template <class T, size_t R>
boost::multi_array<float, 2> NGTDMFeatures2DMRG<T, R>::getNGTDMatrix(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> neighborHoodSum) {
	typedef boost::multi_array<float, 2>  ngtdmat;
	int sizeMatrix = diffGreyLevels.size();
	ngtdmat NGTDMatrix(boost::extents[sizeMatrix][3]);
	int indexOfElement[3] = { 0,0,0 };
	std::vector<T> elementsOfWholeNeighborhoods;
	T sumOfActualNeighborhood;
	T s;
	T actualElement;
	int posActualElement;
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row < inputMatrix.shape()[0]; row++) {
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				indexOfElement[0] = row;
				indexOfElement[1] = col;
				indexOfElement[2] = depth;
				//get actual Element if it is the centre of a whole neighborhood
				actualElement = inputMatrix[row][col][depth];
				if (!std::isnan(actualElement)) {
					//get the sum of the actual neighborhood
					sumOfActualNeighborhood = neighborHoodSum[row][col][depth];
					//get the s_i value
					s = abs(actualElement - sumOfActualNeighborhood);
					//get the position of the actual Element in the diffGreyLevel vector and so also in the matrix
					posActualElement = std::find(diffGreyLevels.begin(), diffGreyLevels.end(), actualElement) - diffGreyLevels.begin();
					//add the s_i value to the right element in the matrix
					NGTDMatrix[posActualElement][2] += s;
					NGTDMatrix[posActualElement][0] += 1;
					//save the actual Element in a vector, so we can calculate later the probabilities
					elementsOfWholeNeighborhoods.push_back(actualElement);
				}
			}
		}
	}
	getProbability(elementsOfWholeNeighborhoods, NGTDMatrix);
	return NGTDMatrix;
}


/*
\brief getProbability
@param[in] vector<T> elementsOfWholeNeighborhood: vector containing all elements of the neighborhood
@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix

The function calculates the probability matrix of the NGTD matrix, as many features are calculated using the probabilities.
*/
template <class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::getProbability(vector<T> elementsOfWholeNeighborhood, boost::multi_array<float, 2> &ngtdMatrix) {
	int numItem = 0;
	for (int actElementIndex = 0; actElementIndex<diffGreyLevels.size(); actElementIndex++) {
		numItem += ngtdMatrix[actElementIndex][0];
	}
	for (int actElementIndex = 0; actElementIndex < diffGreyLevels.size(); actElementIndex++) {
		ngtdMatrix[actElementIndex][1] = ngtdMatrix[actElementIndex][0] / numItem;
	}
}

/*
\brief getNGP

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
@param[out] int ngp : number of discretized values with probability >0
This function returns the number of discretized grey values that have a probability >0
*/
template<class T, size_t R>
int NGTDMFeatures2DMRG<T, R>::getNGP(boost::multi_array<float, 2> ngtdm) {
	int ngp = 0;
	for (int i = 0; i< ngtdm.shape()[0]; i++) {
		if (ngtdm[i][0] != 0) {
			ngp += 1;
		}
	}
	return ngp;
}

/*
\brief getNV

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
@param[out] int nv : sum over all numbers of voxels in a neighborhood
This function returns the sum of all numbers of voxels in each neighborhood
*/
template<class T, size_t R>
int NGTDMFeatures2DMRG<T, R>::getNV(boost::multi_array<float, 2> ngtdm) {
	int nv = 0;
	for (int i = 0; i< ngtdm.shape()[0]; i++) {
		if (!isnan(ngtdm[i][0])) {
			nv += ngtdm[i][0];
		}
	}
	return nv;
}

/*
\brief calculateSumSiPi

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
@param[out] int sumSiPi : sum over all products of column entries
*/
template<class T, size_t R>
double NGTDMFeatures2DMRG<T, R>::calculateSumSiPi(boost::multi_array<float, 2> ngtdm) {
	double sumSiPi = 0;
	for (int row = 0; row<ngtdm.shape()[0]; row++) {
		if (!isnan(ngtdm[row][1] * ngtdm[row][2])) {
			sumSiPi += ngtdm[row][1] * ngtdm[row][2];
		}
	}
	return sumSiPi;
}

/*
\brief calculateSumSi

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
@param[out] int sumSi : sum over all si
*/
template<class T, size_t R>
double NGTDMFeatures2DMRG<T, R>::calculateSumSi(boost::multi_array<float, 2> ngtdm) {
	double sumSi = 0;
	for (int row = 0; row<ngtdm.shape()[0]; row++) {
		if (!isnan(ngtdm[row][2])) {
			sumSi += ngtdm[row][2];
		}
	}
	return sumSi;
}


/*
\brief calculateCoarseness

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
This function calculate the coarseness of the VOI. The coarseness is an indicator of the change of intensity in the VOI. \n
\f$ F_{coarseness} = \frac{1}{\sum{i=1}^{N_{g}}p_{i}s{i}}\f$
*/
template<class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::calculateCoarseness(boost::multi_array<float, 2> ngtdm) {
	double sumSiPi = calculateSumSiPi(ngtdm);
	if (sumSiPi != 0) {
		coarseness = 1 / sumSiPi;
	}
	else {
		coarseness = 0;
	}
}

/*
\brief calculateContrast

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
This function calculate the contrast of the VOI. The value gives information about the range of grey levels, as well as the spatial frequency of intensity changes.\n
\f$ F_{contrast} = ( \frac{1}{N_{g, p}(N_{g, p}-1)}\sum{i=1}^{N_{g}}\sum{j=1}^{N_{g}}p_{i}p{j}(i-j)^{2} )\frac{1}\sum{i=1}^{N_{g}} s_{i}}\f$
*/
template<class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::calculateContrast(boost::multi_array<float, 2> ngtdm) {
	int ng = getNGP(ngtdm);
	int nv = getNV(ngtdm);
	contrast = 0;
	for (int row = 0; row<ngtdm.shape()[0]; row++) {
		for (int rowj = 0; rowj<ngtdm.shape()[0]; rowj++) {
			if (!std::isnan(ngtdm[row][1] * ngtdm[rowj][1] * pow((diffGreyLevels[row] - diffGreyLevels[rowj]), 2))) {
				contrast += ngtdm[row][1] * ngtdm[rowj][1] * pow((diffGreyLevels[row] - diffGreyLevels[rowj]), 2);
			}
		}
	}
	double sumSi = calculateSumSi(ngtdm);
	contrast = contrast*sumSi;
	if (ng > 1 && nv != 0) {
		contrast = contrast / (ng*(ng - 1)*nv);
	}

}
/*
\brief calculateBusyness

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
This function calculate the busyness of the VOI. It gives information about rapid changes in the intensity values from one voxel to the other.\n
\f$ F_{busyness} = ( \frac{\sum{i=1}^{N_{g}}p_{i}s_{i}}{\sum{j=1}^{N_{g}\sum{j=1}^{N_{g}|ip_{i}-jp_{j}|}\f$
*/
template<class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::calculateBusyness(boost::multi_array<float, 2> ngtdm) {
	double sumSiPi = calculateSumSiPi(ngtdm);
	double denominator = 0;
	int ng = getNGP(ngtdm);
	if (ng > 1) {
		for (int row = 0; row < ngtdm.shape()[0]; row++) {
			for (int rowj = 0; rowj < ngtdm.shape()[0]; rowj++) {
				if (ngtdm[row][1] != 0 && ngtdm[rowj][1] != 0) {
					denominator += abs(this->diffGreyLevels[row] * ngtdm[row][1] - ngtdm[rowj][1] * this->diffGreyLevels[rowj]);
				}
			}
		}
		if (denominator > 0) {
			busyness = sumSiPi / denominator;
		}
		else {
			busyness = 0;
		}
	}
	else {
		busyness = 0;
	}
}

/*
\brief calculateComplexity

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
This function calculates the complexity of the VOI. It gives also information about rapid changes in the intensity values from one voxel to the other.\n
\f$ F_{complexity} = ( \frac{1}{N_{v}} \sum{i=1}^{N_{g}}\sum{j=1}^{N_{g}}|i-j|\frac{p_{i}s_{i}+p_{j}s_{j}}{p_{i}+p_{j}}\f$
*/
template<class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::calculateComplexity(boost::multi_array<float, 2> ngtdm) {
	int ng = getNGP(ngtdm);
	int nv = getNV(ngtdm);
	complexity = 0;
	double nominator;
	double denominator;
	for (int row = 0; row < ngtdm.shape()[0]; row++) {
		for (int rowj = 0; rowj < ngtdm.shape()[0]; rowj++) {
			if (ngtdm[row][0] != 0 && ngtdm[rowj][0] != 0) {
				nominator = abs(diffGreyLevels[row] - diffGreyLevels[rowj])*(ngtdm[row][1] * ngtdm[row][2] + ngtdm[rowj][2] * ngtdm[rowj][1]);
				denominator = ngtdm[row][1] + ngtdm[rowj][1];
				if (denominator > 0 && !isnan(nominator) && !isnan(nominator)) {
					complexity += nominator / denominator;
				}
			}
		}
	}
	if (nv != 0) {
		complexity = complexity / nv;
	}
	else {
		complexity = 0;
	}
}


/*
\brief calculateStrength

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
This function calculates the strength of the VOI. It gives also information about rapid changes in the intensity values from one voxel to the other.\n
\f$ F_{strength} = ( \frac{\sum{i=1}^{N_{g}}\sum{j=1}^{N_{g}}(p_{i}+p_{j})j)^{2}}{\sum_{i=1}{N_{g}}s_{i}}\f$
*/
template<class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::calculateStrength(boost::multi_array<float, 2> ngtdm) {
	strength = 0;
	double sumSi = calculateSumSi(ngtdm);
	for (int row = 0; row < ngtdm.shape()[0]; row++) {
		for (int rowj = 0; rowj < ngtdm.shape()[0]; rowj++) {
			if (ngtdm[row][0] != 0 && ngtdm[rowj][0] != 0) {
				strength += (ngtdm[row][1] + ngtdm[rowj][1]) * pow(diffGreyLevels[row] - diffGreyLevels[rowj], 2);
			}
		}
	}
	if (sumSi != 0) {
		strength = strength / sumSi;
	}
	else {
		strength = 0;
	}
}


template <class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::calculateAllNGTDMFeatures2DMRG(NGTDMFeatures2DMRG<T, R> &ngtdm, Image<T, R> imageAttr, boost::multi_array<T, R> neighborHoodSum, vector<float> spacing, ConfigFile config) {
	dist = config.dist;
	this->diffGreyLevels = imageAttr.diffGreyLevels;
	actualSpacing = spacing;
	normNGTDM = config.normNGTDM;
	boost::multi_array<float, 2> ngtdmMatrix = getNGTDMatrix(imageAttr.imageMatrix, neighborHoodSum);
	calculateCoarseness(ngtdmMatrix);
	calculateContrast(ngtdmMatrix);
	calculateBusyness(ngtdmMatrix);
	calculateComplexity(ngtdmMatrix);
	calculateStrength(ngtdmMatrix);
}

template <class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::writeCSVFileNGTDM(NGTDMFeatures2DMRG<T, R> ngtdmFeatures, string outputFolder)
{
	string csvName = outputFolder + "_ngtdmFeatures2Dmrg.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';
	ofstream ngtdmCSV;
	ngtdmCSV.open(name);
	vector<string> features;
	defineNGTDMFeatures2DMRG(features);
	vector<T> ngtdmData;
	extractNGTDMData(ngtdmData, ngtdmFeatures);
	for (int i = 0; i< ngtdmData.size(); i++) {
		ngtdmCSV << "ngtdmFeatures2Dmrg" << "," << features[i] << ",";
		ngtdmCSV << ngtdmData[i];
		ngtdmCSV << "\n";
	}
	ngtdmCSV.close();
}

template <class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::writeOneFileNGTDM(NGTDMFeatures2DMRG<T, R> ngtdmFeatures, ConfigFile config, int &parameterSpaceNr) {
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
	extractNGTDMData(ngtdmData, ngtdmFeatures);
	
	if (config.csvOutput == 1) {
		defineNGTDMFeatures2DMRG(features);
		for (int i = 0; i < ngtdmData.size(); i++) {
			ngtdmCSV << "ngtdmFeatures2Dmrg" << "," << features[i] << ",";
			ngtdmCSV << ngtdmData[i];
			ngtdmCSV << "\n";
		}
	}
	else if (config.ontologyOutput == 1) {
		defineNGTDMFeatures2DMRGOntology(features);
		string featParamSpaceTable = config.outputFolder + "/FeatureParameterSpace_table.csv";
		char * featParamSpaceTableName = new char[featParamSpaceTable.size() + 1];
		std::copy(featParamSpaceTable.begin(), featParamSpaceTable.end(), featParamSpaceTableName);
		featParamSpaceTableName[featParamSpaceTable.size()] = '\0';

		ofstream featSpaceTable;
		featSpaceTable.open(featParamSpaceTableName, std::ios_base::app);
		parameterSpaceNr += 1;
		string parameterSpaceName = "FeatureParameterSpace_" + std::to_string(parameterSpaceNr);
		featSpaceTable << parameterSpaceName << "," << "2Dmrg" << "," << config.imageSpaceName << "," << config.interpolationMethod << "\n";
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
void NGTDMFeatures2DMRG<T, R>::defineNGTDMFeatures2DMRG(vector<string> &features) {
	features.push_back("coarseness");
	features.push_back("contrast");
	features.push_back("busyness");
	features.push_back("complexity");
	features.push_back("strength");
}

template <class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::defineNGTDMFeatures2DMRGOntology(vector<string> &features) {
	features.push_back("Fngt.coarseness");
	features.push_back("Fngt.contrast");
	features.push_back("Fngt.busyness");
	features.push_back("Fngt.complexity");
	features.push_back("Fngt.strength");
}



template <class T, size_t R>
void NGTDMFeatures2DMRG<T, R>::extractNGTDMData(vector<T> &ngtdmData, NGTDMFeatures2DMRG<T, R> ngtdmFeatures) {
	ngtdmData.push_back(ngtdmFeatures.coarseness);
	ngtdmData.push_back(ngtdmFeatures.contrast);
	ngtdmData.push_back(ngtdmFeatures.busyness);
	ngtdmData.push_back(ngtdmFeatures.complexity);
	ngtdmData.push_back(ngtdmFeatures.strength);
}

#endif // NGTDM2DMRG_H_INCLUDED
