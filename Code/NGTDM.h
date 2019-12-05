#ifndef NGTDM_H_INCLUDED
#define NGTDM_H_INCLUDED

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

template <class T,  size_t R>
class NGTDMFeatures{

    private:
		//values for weighting the entries (novel and uncommon features)
		vector<float> actualSpacing;
		string normNGTDM;
		int dist;
		boost::multi_array<float, 2> getNGTDMatrix(boost::multi_array<T, R> inputMatrix);
		T getNeighborhood(boost::multi_array<T, R> inputMatrix, int *indexOfElement);
        void extractNGTDMData(vector<T> &ngtdmData, NGTDMFeatures<T, R> ngtdmFeatures);
        
    public:
		vector<T> diffGreyLevels;
		float coarseness;
		float contrast;
		float busyness;
		float complexity;
		float strength;

        void getProbability(vector<T> elementsOfWholeNeighborhood, boost::multi_array<float, 2> &ngtdMatrix);
        float calculateSumSi(boost::multi_array<float, 2> ngtdm);
        float calculateSumSiPi(boost::multi_array<float, 2> ngtdm);
        int getNGP(boost::multi_array<float, 2> ngtdm);
        int getNV(boost::multi_array<float, 2> ngtdm);

        void calculateStrength(boost::multi_array<float, 2> ngtdm);
        void calculateComplexity(boost::multi_array<float, 2> ngtdm);
        void calculateCoarseness(boost::multi_array<float, 2> ngtdm);
        void calculateContrast(boost::multi_array<float, 2> ngtdm);
        void calculateBusyness(boost::multi_array<float, 2> ngtdm);
        void calculateAllNGTDMFeatures(NGTDMFeatures<T,R> &ngtdm, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<float> spacing, ConfigFile config);
        void writeCSVFileNGTDM(NGTDMFeatures<T, R> ngtdm, string outputFolder);
		void writeOneFileNGTDM(NGTDMFeatures<T, R> ngtdm, string outputFolder);
		void defineNGTDMFeatures(vector<string> &features);
};

/*!
\brief getNGTDMatrix
@param boost multi array inputMatrix: matrix filled with intensity values
@param int dist: size of neighborhood
@param[out] boost multi array: filled NGTD matrix

The function fills the NGTDMatrix with the corresponding values
*/
template <class T, size_t R>
boost::multi_array<float, 2> NGTDMFeatures<T, R>::getNGTDMatrix(boost::multi_array<T,R> inputMatrix){
    typedef boost::multi_array<float, 2>  ngtdmat;
    int sizeMatrix = diffGreyLevels.size();
    ngtdmat NGTDMatrix(boost::extents[sizeMatrix][3]);
    int indexOfElement[3] = {0,0,0};
    std::vector<T> elementsOfWholeNeighborhoods;
    T sumOfActualNeighborhood;
    T s;
    T actualElement;
    int posActualElement;
    for(int depth = 0; depth < inputMatrix.shape()[2]; depth++){
        for(int row = 0; row < inputMatrix.shape()[0]; row++){
            for(int col = 0; col < inputMatrix.shape()[1]; col++){
                indexOfElement[0] = row;
                indexOfElement[1] = col;
                indexOfElement[2] = depth;
                //get actual Element if it is the centre of a whole neighborhood
                actualElement = inputMatrix[row][col][depth];
                if(!std::isnan(actualElement)){
	               //get the sum of the actual neighborhood
                    sumOfActualNeighborhood = getNeighborhood(inputMatrix, indexOfElement);
                    //get the s_i value
                    s = abs(actualElement - sumOfActualNeighborhood);
                    //get the position of the actual Element in the diffGreyLevel vector and so also in the matrix
                    posActualElement = std::find(diffGreyLevels.begin(), diffGreyLevels.end(), actualElement) - diffGreyLevels.begin();
                   //add the s_i value to the right element in the matrix
                    NGTDMatrix[posActualElement][2] += s;
                    NGTDMatrix[posActualElement][0] +=1;
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
\brief getNeighborhood
@param[in] boost multi array inputMatrix: matrix filled with intensity values
@param[in] int indexOfElement: index of the actual element(for which the neighborhood is calculated)
@param[in] int dist: size of neighborhood
@param[out] T sum: average sume of all elements in the neighborhood except the center

The function get the values of all neighbors except the center of a certain element in a certain distance dist.
It calculates the average sum of all these elements 
*/
template <class T, size_t R>
T NGTDMFeatures<T, R>::getNeighborhood(boost::multi_array<T,R> inputMatrix, int *indexOfElement){
    vector<T> neighborhood;
    T actElement;
	float weight;
    for(int i = -dist; i < dist+1; i++){
        for(int j = -dist; j < dist+1; j++){
            if(i!=0 || j!=0){
                if(indexOfElement[0]+i > -1 && indexOfElement[0]+i < inputMatrix.shape()[0] && indexOfElement[1]+j > -1 && indexOfElement[1]+j < inputMatrix.shape()[1] ){
                    actElement=inputMatrix[indexOfElement[0]+i][indexOfElement[1]+j][indexOfElement[2]];
                    if(!std::isnan(actElement)){
						weight = calculateWeight2D(i, j, normNGTDM, actualSpacing);
                        neighborhood.push_back(weight * actElement);
                    }
                }
            }
        }
    }
    T sum = accumulate(neighborhood.begin(), neighborhood.end(), 0);
    sum=sum/neighborhood.size();
    return sum;
}
/*
\brief getProbability
@param[in] vector<T> elementsOfWholeNeighborhood: vector containing all elements of the neighborhood
@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix

The function calculates the probability matrix of the NGTD matrix, as many features are calculated using the probabilities. 
*/
template <class T, size_t R>
void NGTDMFeatures<T, R>::getProbability(vector<T> elementsOfWholeNeighborhood, boost::multi_array<float, 2> &ngtdMatrix){
    int numItem = 0;
    for(int actElementIndex =0; actElementIndex<diffGreyLevels.size(); actElementIndex++){
        numItem += ngtdMatrix[actElementIndex][0];
    }
    for(int actElementIndex = 0; actElementIndex < diffGreyLevels.size(); actElementIndex++){
        ngtdMatrix[actElementIndex][1] = ngtdMatrix[actElementIndex][0]/numItem;
    }
}

/*
\brief getNGP

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
@param[out] int ngp : number of discretized values with probability >0
This function returns the number of discretized grey values that have a probability >0
*/
template<class T, size_t R>
int NGTDMFeatures<T, R>::getNGP(boost::multi_array<float, 2> ngtdm){
    int ngp = 0;
    for(int i = 0; i< ngtdm.shape()[0]; i++){
        if(ngtdm[i][0]!=0){
            ngp +=1;
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
int NGTDMFeatures<T, R>::getNV(boost::multi_array<float, 2> ngtdm){
    int nv =0;
    for(int i = 0; i< ngtdm.shape()[0]; i++){
        nv+=ngtdm[i][0];
    }
    return nv;
}

/*
\brief calculateSumSiPi

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
@param[out] int sumSiPi : sum over all products of column entries
*/
template<class T, size_t R>
float NGTDMFeatures<T, R>::calculateSumSiPi(boost::multi_array<float, 2> ngtdm) {
	float sumSiPi = 0;
	for (int row = 0; row<ngtdm.shape()[0]; row++) {
		sumSiPi += ngtdm[row][1] * ngtdm[row][2];
	}
	return sumSiPi;
}

/*
\brief calculateSumSi

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
@param[out] int sumSi : sum over all si
*/
template<class T, size_t R>
float NGTDMFeatures<T, R>::calculateSumSi(boost::multi_array<float, 2> ngtdm) {
	float sumSi = 0;
	for (int row = 0; row<ngtdm.shape()[0]; row++) {
		sumSi += ngtdm[row][2];
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
void NGTDMFeatures<T, R>::calculateCoarseness(boost::multi_array<float, 2> ngtdm){
    float sumSiPi = calculateSumSiPi(ngtdm);
    coarseness = 1/sumSiPi;
}

/*
\brief calculateContrast

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
This function calculate the contrast of the VOI. The value gives information about the range of grey levels, as well as the spatial frequency of intensity changes.\n
\f$ F_{contrast} = ( \frac{1}{N_{g, p}(N_{g, p}-1)}\sum{i=1}^{N_{g}}\sum{j=1}^{N_{g}}p_{i}p{j}(i-j)^{2} )\frac{1}\sum{i=1}^{N_{g}} s_{i}}\f$ 
*/
template<class T, size_t R>
void NGTDMFeatures<T, R>::calculateContrast(boost::multi_array<float, 2> ngtdm){
    int ng = getNGP(ngtdm);
    int nv = getNV(ngtdm);
    contrast = 0;
    for(int row=0; row<ngtdm.shape()[0]; row++){
        for(int rowj = 0; rowj<ngtdm.shape()[0]; rowj++){
            contrast += ngtdm[row][1]*ngtdm[rowj][1]*pow((diffGreyLevels[row]-diffGreyLevels[rowj]), 2);
        }
    }
    float sumSi = calculateSumSi(ngtdm);
    contrast = contrast*sumSi;
    contrast = contrast/(ng*(ng-1)*nv);
}
/*
\brief calculateBusyness

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
This function calculate the busyness of the VOI. It gives information about rapid changes in the intensity values from one voxel to the other.\n
\f$ F_{busyness} = ( \frac{\sum{i=1}^{N_{g}}p_{i}s_{i}}{\sum{j=1}^{N_{g}\sum{j=1}^{N_{g}|ip_{i}-jp_{j}|}\f$ 
*/
template<class T, size_t R>
void NGTDMFeatures<T, R>::calculateBusyness(boost::multi_array<float, 2> ngtdm){
    float sumSiPi = calculateSumSiPi(ngtdm);
    float denominator = 0;
    int ng = getNGP(ngtdm);
    if(ng > 1){
        for(int row = 0; row < ngtdm.shape()[0]; row++){
            for(int rowj = 0; rowj < ngtdm.shape()[0]; rowj++){
                if (ngtdm[row][1] != 0 && ngtdm[rowj][1] != 0){
                    denominator += abs(this->diffGreyLevels[row]* ngtdm[row][1] - ngtdm[rowj][1]*this->diffGreyLevels[rowj]);
                }
            }
        }
        busyness = sumSiPi/denominator;
    }
    else{
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
void NGTDMFeatures<T, R>::calculateComplexity(boost::multi_array<float, 2> ngtdm){
    int ng = getNGP(ngtdm);
    int nv = getNV(ngtdm);
    complexity = 0;
    float nominator;
    float denominator;
    for(int row = 0; row < ngtdm.shape()[0]; row++){
        for(int rowj = 0; rowj < ngtdm.shape()[0]; rowj++){
            if (ngtdm[row][0] != 0 && ngtdm[rowj][0] != 0){
                nominator = abs(diffGreyLevels[row]-diffGreyLevels[rowj])*(ngtdm[row][1]* ngtdm[row][2] + ngtdm[rowj][2]*ngtdm[rowj][1]);
                denominator = ngtdm[row][1]+ngtdm[rowj][1];
                complexity += nominator/denominator;
            }
        }
    }
    complexity = complexity/nv;
}


/*
\brief calculateStrength

@param[in] boost::multi_array<double, 2> ngtdMatrix: ngtd matrix
This function calculates the strength of the VOI. It gives also information about rapid changes in the intensity values from one voxel to the other.\n
\f$ F_{strength} = ( \frac{\sum{i=1}^{N_{g}}\sum{j=1}^{N_{g}}(p_{i}+p_{j})j)^{2}}{\sum_{i=1}{N_{g}}s_{i}}\f$
*/
template<class T, size_t R>
void NGTDMFeatures<T, R>::calculateStrength(boost::multi_array<float, 2> ngtdm){
    strength = 0;
    float sumSi = calculateSumSi(ngtdm);
        for(int row = 0; row < ngtdm.shape()[0]; row++){
        for(int rowj = 0; rowj < ngtdm.shape()[0]; rowj++){
            if (ngtdm[row][0] != 0 && ngtdm[rowj][0] != 0){
                strength += (ngtdm[row][1] + ngtdm[rowj][1]) * pow(diffGreyLevels[row]-diffGreyLevels[rowj],2);
            }
        }
    }
    strength = strength/sumSi;
}


template <class T, size_t R>
void NGTDMFeatures<T, R>::calculateAllNGTDMFeatures(NGTDMFeatures<T,R> &ngtdmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<float> spacing, ConfigFile config){
	dist = config.dist;
	this->diffGreyLevels = diffGrey;
	actualSpacing = spacing;
	normNGTDM = config.normNGTDM;
    boost::multi_array<float, 2> ngtdm =getNGTDMatrix(inputMatrix);
    calculateCoarseness(ngtdm);
    calculateContrast(ngtdm);
    calculateBusyness(ngtdm);
    calculateComplexity(ngtdm);
    calculateStrength(ngtdm);
}

template <class T,  size_t R>
void NGTDMFeatures<T,R>::writeCSVFileNGTDM(NGTDMFeatures<T,R> ngtdmFeatures, string outputFolder)
{
    string csvName = outputFolder + "/ngtdmFeatures2Dmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';
    ofstream ngtdmCSV;
    ngtdmCSV.open (name);
    vector<string> features;
    defineNGTDMFeatures(features);
    vector<T> ngtdmData;
    extractNGTDMData(ngtdmData, ngtdmFeatures);
    for(int i = 0; i< ngtdmData.size(); i++){
        ngtdmCSV <<"ngtdmFeatures2Dmrg"<<","<< features[i] <<",";
        ngtdmCSV << ngtdmData[i];
        ngtdmCSV << "\n";
    }
    ngtdmCSV.close();
}

template <class T, size_t R>
void NGTDMFeatures<T, R>::writeOneFileNGTDM(NGTDMFeatures<T, R> ngtdmFeatures, string outputFolder) {
	string csvName = outputFolder + "/radiomicsFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';
	ofstream ngtdmCSV;
	ngtdmCSV.open(name, std::ios_base::app);
	vector<string> features;
	defineNGTDMFeatures(features);
	vector<T> ngtdmData;
	extractNGTDMData(ngtdmData, ngtdmFeatures);
	for (int i = 0; i< ngtdmData.size(); i++) {
		ngtdmCSV << "ngtdmFeatures2Dmrg" << "," << features[i] << ",";
		ngtdmCSV << ngtdmData[i];
		ngtdmCSV << "\n";
	}
	ngtdmCSV.close();
}


template <class T,  size_t R>
void NGTDMFeatures<T,R>::defineNGTDMFeatures(vector<string> &features){
    features.push_back("coarseness");
    features.push_back("contrast");
    features.push_back("busyness");
    features.push_back("complexity");
    features.push_back("strength");
}



template <class T,  size_t R>
void NGTDMFeatures<T,R>::extractNGTDMData(vector<T> &ngtdmData, NGTDMFeatures<T, R> ngtdmFeatures){
    ngtdmData.push_back(ngtdmFeatures.coarseness);
    ngtdmData.push_back(ngtdmFeatures.contrast);
    ngtdmData.push_back(ngtdmFeatures.busyness);
    ngtdmData.push_back(ngtdmFeatures.complexity);
    ngtdmData.push_back(ngtdmFeatures.strength);
}

#endif // NGTDM_H_INCLUDED
