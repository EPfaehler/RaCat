#ifndef NGLDMFEATURES2DMRG_H_INCLUDED
#define NGLDMFEATURES2DMRG_H_INCLUDED
#include "GLRLMFeatures.h"


/*! \file */
/*!
The class NGLDM is the class of the Neighborhood Grey Level Dependence Matrices. \n
A neighborhood are all voxels around one voxel within distance dist. \n
A voxel is dependent from the other, if \f$ |X_{c}-X_{m}|<a\f$, where \f$X_{c}\f$ is the center voxel and \f$X_{m}\f$ are the 
other voxels in the neighborhood. The number of dependent voxels in a neighborhood is counted. \n
\f$ a \f$ is called the coarseness parameter. \n
The neighborhoods are checked for every voxel. \n
\f$ s_ = s(i,j) \f$ is the number of neighborhoods with center voxel with grey level i and dependece k = j-1 \n

The most features are already defined in the GLRLM features.
*/
template <class T,  size_t R>
class NGLDMFeatures2DMRG : public GLRLMFeatures<T,R>{
    private:
		vector<T> sumProbRows;
		vector<T> sumProbCols;
		vector<float> rowSums;
		vector<float> colSums;

		int dist;
		int coarseParam;

        void extractNGLDMData(vector<T> &ngldmData, NGLDMFeatures2DMRG<T, R> ngldmFeatures);
        boost::multi_array<float, 2> getMatrix(boost::multi_array<T,R> inputMatrix);

		

    public:
        float dependenceCountEnergy;
		int findIndex(vector<T> array, int size, T target);
        void writeCSVFileNGLDM2DMRG(NGLDMFeatures2DMRG<T,R> ngldmFeat, string outputFolder);
		void writeOneFileNGLDM2DMRG(NGLDMFeatures2DMRG<T, R> ngldmFeat, ConfigFile config, int &parameterSpaceNr);
        void calculateAllNGLDMFeatures2DMRG(NGLDMFeatures2DMRG<T,R> &ngldmFeatures, Image<T, R> imageAttr, boost::multi_array<T, R> ngldmMatrix, ConfigFile config);
        void calculateDependenceCountEnergy(boost::multi_array<float,2> probMatrix);
        void defineNGLDMFeatures(vector<string> &features);
		void defineNGLDMFeaturesOntology(vector<string> &features);

};

/*!
\brief findIndex
@param[in] vector array: vector containing elements
@param[in] int size: size of vector
@param[in] T target: element which has to be searched in the vector
@param[out] position of target element in vector
This function finds a specific element in an array and returns the index of the target element.
In our case, it is used to find the position of the actual grey level in the vector containing all grey levels.
*/
template <class T, size_t R>
int NGLDMFeatures2DMRG<T, R>::findIndex(vector<T> array, int size, T target) {
	int i = 0;
	while ((i < size) && (array[i] != target)) i++;
	return (i < size) ? (i) : (-1);
}

/*!
\brief getMatrix
@param boost multi array ngldmMatrix: Matrix containing NGLDM elements for every slice, this matrix just has to be merged to a 2D matrix
@param[out] boost multi array: filled NGLD matrix

The function fills the NGLDMatrix with the corresponding values using the function getNeighborGreyLevels. It 
checks voxel by voxel the neighborhood in the distance that is set by the user.
*/
template <class T, size_t R>
boost::multi_array<float, 2> NGLDMFeatures2DMRG<T, R>::getMatrix(boost::multi_array<T,R> ngldmMatrix){
	boost::multi_array<float, 2> NGLDMatrix(boost::extents[ngldmMatrix.shape()[0]][ngldmMatrix.shape()[1]]);
	//check every element of the VOI
	for (int depth = 0; depth < ngldmMatrix.shape()[2]; depth++) {
		for (int row = 0; row < ngldmMatrix.shape()[0]; row++) {
			for (int col = 0; col < ngldmMatrix.shape()[1]; col++) {
				NGLDMatrix[row][col] += ngldmMatrix[row][col][depth];
			}

		}
	}
    return NGLDMatrix;
}



/*!
\brief calculateDependenceCountEnergy
@param[in] boost multi array probMatrix: probability matrix filled

The function calculate the dependence count energy: \f$ F_{countEnergy} = \sum_{i=1}^{N_{g} \sum_{j=1}^{N_{g} p_{ij}^{2} \f$.
*/
template <class T, size_t R>
void NGLDMFeatures2DMRG<T, R>::calculateDependenceCountEnergy(boost::multi_array<float,2> probMatrix){
    dependenceCountEnergy=0;
    for(int row=0; row<probMatrix.shape()[0]; row++){
        for(int col=0; col<probMatrix.shape()[1]; col++){
			if (!std::isnan(pow(probMatrix[row][col], 2))) {
				dependenceCountEnergy += pow(probMatrix[row][col], 2);
			}
        }

    }
}


template <class T, size_t R>
void NGLDMFeatures2DMRG<T, R>::calculateAllNGLDMFeatures2DMRG(NGLDMFeatures2DMRG<T,R> &ngldmFeatures, Image<T, R> imageAttr, boost::multi_array<T, R> ngldmMatrix, ConfigFile config){
    this->diffGreyLevels = imageAttr.diffGreyLevels;
	ngldmFeatures.getConfigValues(config);
	//get values from config file
	coarseParam = config.coarsenessParam;
	this->dist = config.distNGLDM;
	//get NGLDM 
    boost::multi_array<float,2> NGLDM=ngldmFeatures.getMatrix(ngldmMatrix);
    float totalSum = ngldmFeatures.calculateTotalSum(NGLDM);

    rowSums=ngldmFeatures.calculateRowSums(NGLDM);
    colSums = ngldmFeatures.calculateColSums(NGLDM);
    ngldmFeatures.calculateShortRunEmphasis(rowSums, totalSum);
    ngldmFeatures.calculateLongRunEmphasis(rowSums, totalSum);
    ngldmFeatures.calculateLowGreyEmph(colSums, totalSum);
    ngldmFeatures.calculateHighGreyEmph(colSums, totalSum);
    ngldmFeatures.calculateShortRunLow(NGLDM, totalSum);
    ngldmFeatures.calculateShortRunHigh(NGLDM, totalSum);
    ngldmFeatures.calculateLongRunLowEmph(NGLDM, totalSum);
    ngldmFeatures.calculateLongRunHighEmph(NGLDM, totalSum);
    ngldmFeatures.calculateGreyNonUniformity(colSums, totalSum);
    ngldmFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
    ngldmFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
    ngldmFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
    ngldmFeatures.calculateRunPercentage3D(imageAttr.vectorOfMatrixElements, totalSum, 1);
    boost::multi_array<float,2> probMatrix = ngldmFeatures.calculateProbMatrix(NGLDM, totalSum);
    float meanGrey = ngldmFeatures.calculateMeanProbGrey(probMatrix);
    ngldmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);

    float meanRun = ngldmFeatures.calculateMeanProbRun(probMatrix);
    ngldmFeatures.calculateRunLengthVar(probMatrix, meanRun);
    ngldmFeatures.calculateRunEntropy(probMatrix);
    ngldmFeatures.calculateDependenceCountEnergy(probMatrix);
}

template <class T, size_t R>
void NGLDMFeatures2DMRG<T, R>::extractNGLDMData(vector<T> &ngldmData, NGLDMFeatures2DMRG<T, R> ngldmFeatures){

    ngldmData.push_back(ngldmFeatures.shortRunEmphasis);
    ngldmData.push_back(ngldmFeatures.longRunEmphasis);
    ngldmData.push_back(ngldmFeatures.lowGreyEmph);
    ngldmData.push_back(ngldmFeatures.highGreyEmph);
    ngldmData.push_back(ngldmFeatures.shortRunLow);
    ngldmData.push_back(ngldmFeatures.shortRunHigh);
    ngldmData.push_back(ngldmFeatures.longRunLowEmph);
    ngldmData.push_back(ngldmFeatures.longRunHighEmph);
    ngldmData.push_back(ngldmFeatures.greyNonUniformity);
    ngldmData.push_back(ngldmFeatures.greyNonUniformityNorm);
    ngldmData.push_back(ngldmFeatures.runLengthNonUniformity);   //check if I have to give the matrix
    ngldmData.push_back(ngldmFeatures.runLengthNonUniformityNorm);
    ngldmData.push_back(ngldmFeatures.runPercentage);
    ngldmData.push_back(ngldmFeatures.greyLevelVar);
    ngldmData.push_back(ngldmFeatures.runLengthVar);
    ngldmData.push_back(ngldmFeatures.runEntropy);
    ngldmData.push_back(ngldmFeatures.dependenceCountEnergy);

}

template <class T, size_t R>
void NGLDMFeatures2DMRG<T, R>::writeCSVFileNGLDM2DMRG(NGLDMFeatures2DMRG<T,R> ngldmFeat, string outputFolder)
{
    string csvName = outputFolder + "_ngldmFeatures2Dmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream ngldmCSV;
    ngldmCSV.open (name);
    vector<string> features;
    defineNGLDMFeatures(features);

    vector<T> ngldmData;
    extractNGLDMData(ngldmData, ngldmFeat);
    for(int i = 0; i< ngldmData.size(); i++){
        ngldmCSV << "ngldmFeatures2Dmrg"<<","<<features[i] <<",";
        ngldmCSV << ngldmData[i];
        ngldmCSV << "\n";
    }
    ngldmCSV.close();
}

template <class T, size_t R>
void NGLDMFeatures2DMRG<T, R>::writeOneFileNGLDM2DMRG(NGLDMFeatures2DMRG<T, R> ngldmFeat, ConfigFile config, int &parameterSpaceNr) {
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

	ofstream ngldmCSV;
	ngldmCSV.open(name, std::ios_base::app);
	vector<string> features;

	vector<T> ngldmData;
	extractNGLDMData(ngldmData, ngldmFeat);
	
	if (config.csvOutput == 1) {
		defineNGLDMFeatures(features);
		for (int i = 0; i < ngldmData.size(); i++) {
			ngldmCSV << "ngldmFeatures2Dmrg" << "," << features[i] << ",";
			ngldmCSV << ngldmData[i];
			ngldmCSV << "\n";
		}
	}
	else if (config.ontologyOutput == 1) {
		features.clear();
		defineNGLDMFeaturesOntology(features);
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

		for (int i = 0; i < ngldmData.size(); i++) {
			ngldmCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			ngldmCSV << ngldmData[i] << "," << parameterSpaceName << "," << config.calculationSpaceName;
			ngldmCSV << "\n";
		}

	}
	ngldmCSV.close();
}

template <class T, size_t R>
void NGLDMFeatures2DMRG<T, R>::defineNGLDMFeatures(vector<string> &features){
    features.push_back("Low dependence emphasis");
    features.push_back("High dependence emphasis");
    features.push_back("Low grey level count emphasis");
    features.push_back("High grey level count emphasis");
    features.push_back("Low dependence low grey level emphasis");
    features.push_back("Low dependence high grey level emphasis");
    features.push_back("High dependence low grey level emphasis");
    features.push_back("High dependence high grey level emphasis");
    features.push_back("Grey level non uniformity");
    features.push_back("Grey level non uniformity normalized");
    features.push_back("Dependence count non uniformity");
    features.push_back("Dependence count non uniformity normalized");
    features.push_back("Dependence count percentage");
    features.push_back("Grey level variance");
    features.push_back("Dependence count variance");
    features.push_back("Dependence count entropy");
    features.push_back("dependence Count Energy");

}

template <class T, size_t R>
void NGLDMFeatures2DMRG<T, R>::defineNGLDMFeaturesOntology(vector<string> &features) {
	features.push_back("Fngl.lde");
	features.push_back("Fngl.hde");
	features.push_back("Fngl.lgce");
	features.push_back("Fngl.hgce");
	features.push_back("Fngl.ldlge");
	features.push_back("Fngl.ldhge");
	features.push_back("Fngl.hdlge");
	features.push_back("Fngl.hdhge");
	features.push_back("Fngl.glnu");
	features.push_back("Fngl.glnu.norm");
	features.push_back("Fngl.dcnu");
	features.push_back("Fngl.dcnu.norm");
	features.push_back("dependence count percentage");
	features.push_back("Fngl.gl.var");
	features.push_back("Fngl.dc.var");
	features.push_back("Fngl.dc.entr");
	features.push_back("Fngl.dc.energy");

}
#endif // NGLDMFEATURES2DMRG_H_INCLUDED
