#ifndef NGTDM3D_H_INCLUDED
#define NGTDM3D_H_INCLUDED

#include "NGTDM2DMRG.h"

/*! \file */
/*!
The class NGTDM3D inherites from the class NGTDM. Here the matrix is calculated looking at the 3D VOI and not slice by slice. \n
All other definitions are the same. \n
These matrices combine the sum of grey level differences of voxels with intensity value i and the average
discretised grey levels of a neighborhood with distance dist from the actual voxel. \n
The average grey level within a neighborhood is defined as:
\f$ A_{i} = \frac{1}{W} \sum_{k_{x} = -dist}^{dist} \sum_{k_{y} = -dist}^{dist} \sum_{k_{z} = -dist}^{dist} X_{dgl}(j_{x} + k_{x}, j_{y} + k_{y}, j_{z} + k_{z}) \f$ \n
where \f$ k_{x}, k_{y},k_{z} !=0 \f$ and \f$ W = (2dist+1)\f$\n
Now, let \f$ n_{i} \f$ be the number of voxels with grey level i that have a complete neighborhood. \n
The entry in the NGTDM matrix is then: \f$ s_{i} = \sum^{n_{i}} i-A_{i}, if n_{i} >0 \f$ and 0 otherwise.\n
*/


template <class T, size_t R = 3>
class NGTDMFeatures3D : NGTDMFeatures2DMRG<T, R> {

private:
	NGTDMFeatures2DMRG<T, R> ngtdm;
	//these are the values in case we want to calculate novel/uncommon features
	vector<double> actualSpacing;
	string normNGTDM;
	//distance defined by the user
	int dist;

	void extractNGTDMData3D(vector<T> &ngtdmData, NGTDMFeatures3D<T, R> ngtdmFeatures);
	boost::multi_array<double, 2> getNGTDMatrix3D(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> sumMatrix);

public:
	void getProbability(vector<T> elementsOfWholeNeighborhood, boost::multi_array<double, 2> &ngtdMatrix);
	void calculateAllNGTDMFeatures3D(NGTDMFeatures3D<T, R> &ngtdm, boost::multi_array<T, R> sumMatrix, Image<T, R> imageAttr, vector<double> spacing, ConfigFile config);
	void writeCSVFileNGTDM3D(NGTDMFeatures3D<T, R> ngtdmFeatures, string outputFolder);
	void writeOneFileNGTDM3D(NGTDMFeatures3D<T, R> ngtdmFeatures, string outputFolder);


};


/*!
\brief getNGTDMatrix3D
In this function the NGTDM is filled for the 3D case.
@param[in] boost multi array inputMatrix: matrix filled with intensity values
@param[in] boost multi array inputMatrix: sumMatrix calculated before, this matrix contains for every voxel the sum of the neighborhood
@param[out] boost multi array: filled NGTD matrix

*/
template <class T, size_t R>
boost::multi_array<double, 2> NGTDMFeatures3D<T, R>::getNGTDMatrix3D(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> sumMatrix) {
	typedef boost::multi_array<double, 2>  ngtdmat;
	int sizeMatrix = this->diffGreyLevels.size();
	ngtdmat NGTDMatrix(boost::extents[sizeMatrix][3]);
	std::vector<T> elementsOfWholeNeighborhoods;
	T sumOfActualNeighborhood;
	T s;
	T actualElement;
	int posActualElement;
	//calculate NGTDM for every element of the original matrix
	for (int depth = 0; depth<inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row<inputMatrix.shape()[0]; row++) {
			for (int col = 0; col<inputMatrix.shape()[1]; col++) {
				
				//get actual Element if it is the centre of a whole neighborhood
				actualElement = inputMatrix[row][col][depth];
				if (!std::isnan(actualElement)) {
					//get the sum of the actual neighborhood
					sumOfActualNeighborhood = sumMatrix[row][col][depth ];
					//get the s_i value
					s = abs(actualElement - sumOfActualNeighborhood);
					posActualElement = std::find(this->diffGreyLevels.begin(), this->diffGreyLevels.end(), actualElement) - this->diffGreyLevels.begin();
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



/*!
\brief getProbability
In this function the second row is filled with the probability.
@param[inelementsOfWholeNeighborhood elementsOfWholeNeighborhood: matrix filled with intensity values
@param[in] boost multi array inputMatrix: sumMatrix calculated before, this matrix contains for every voxel the sum of the neighborhood
@param[out] boost multi array: filled NGTD matrix

*/
template <class T, size_t R>
void NGTDMFeatures3D<T, R>::getProbability(vector<T> elementsOfWholeNeighborhood, boost::multi_array<double, 2> &ngtdMatrix) {
	int numItem = 0;
	for (int actElementIndex = 0; actElementIndex<this->diffGreyLevels.size(); actElementIndex++) {
		numItem += ngtdMatrix[actElementIndex][0];
	}
	for (int actElementIndex = 0; actElementIndex<this->diffGreyLevels.size(); actElementIndex++) {
		ngtdMatrix[actElementIndex][1] = ngtdMatrix[actElementIndex][0] / numItem;
	}
}


template <class T, size_t R>
void NGTDMFeatures3D<T, R>::calculateAllNGTDMFeatures3D(NGTDMFeatures3D<T, R> &ngtdm, boost::multi_array<T, R> sumMatrix, Image<T, R> imageAttr, vector<double> spacing, ConfigFile config) {
	this->diffGreyLevels = imageAttr.diffGreyLevels;
	actualSpacing = spacing;
	normNGTDM = config.normNGTDM;
	dist = config.dist;
	boost::multi_array<double, 2> ngtdMatrix = getNGTDMatrix3D(imageAttr.imageMatrix, sumMatrix);

	ngtdm.calculateCoarseness(ngtdMatrix);
	ngtdm.calculateContrast(ngtdMatrix);
	ngtdm.calculateBusyness(ngtdMatrix);
	ngtdm.calculateComplexity(ngtdMatrix);
	ngtdm.calculateStrength(ngtdMatrix);
}

template <class T, size_t R>
void NGTDMFeatures3D<T, R>::writeCSVFileNGTDM3D(NGTDMFeatures3D<T, R> ngtdmFeatures, string outputFolder)
{
	string csvName = outputFolder + "/ngtdmFeatures3D.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream ngtdmCSV;
	ngtdmCSV.open(name);
	vector<string> features;
	ngtdm.defineNGTDMFeatures2DMRG(features);

	vector<T> ngtdmData;
	extractNGTDMData3D(ngtdmData, ngtdmFeatures);
	for (int i = 0; i< ngtdmData.size(); i++) {
		ngtdmCSV << "ngtdmFeatures3D" << "," << features[i] << ",";
		ngtdmCSV << ngtdmData[i];
		ngtdmCSV << "\n";
	}
	ngtdmCSV.close();
}


template <class T, size_t R>
void NGTDMFeatures3D<T, R>::writeOneFileNGTDM3D(NGTDMFeatures3D<T, R> ngtdmFeatures, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream ngtdmCSV;
	ngtdmCSV.open(name, std::ios_base::app);
	vector<string> features;
	ngtdm.defineNGTDMFeatures2DMRG(features);

	vector<T> ngtdmData;
	extractNGTDMData3D(ngtdmData, ngtdmFeatures);
	for (int i = 0; i< ngtdmData.size(); i++) {
		ngtdmCSV << "ngtdmFeatures3D" << "," << features[i] << ",";
		ngtdmCSV << ngtdmData[i];
		ngtdmCSV << "\n";
	}
	ngtdmCSV.close();
}


template <class T, size_t R>

void NGTDMFeatures3D<T, R>::extractNGTDMData3D(vector<T> &ngtdmData, NGTDMFeatures3D<T, R> ngtdmFeatures) {
	ngtdmData.push_back(ngtdmFeatures.coarseness);
	ngtdmData.push_back(ngtdmFeatures.contrast);
	ngtdmData.push_back(ngtdmFeatures.busyness);
	ngtdmData.push_back(ngtdmFeatures.complexity);
	ngtdmData.push_back(ngtdmFeatures.strength);

}
#endif // NGTDM3D_H_INCLUDED
