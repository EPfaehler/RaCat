#ifndef GLRLMFEATURES3D_H_INCLUDED
#define GLRLMFEATURES3D_H_INCLUDED

#include "GLRLMFeatures3DAVG.h"

/*! \file */

/*!
The class GLRLMFeatures3D herites from the class GLRLMFeatures. \n
To calculate the run length it goes also in the depth. \n
It calculates the feature values for every angle and calculates then
the mean value of the features. \n
All feature calculations are defined in the class GLRLMFeatures. \n
*/

template <class T,  size_t R=3>
class GLRLMFeatures3D : GLRLMFeatures<T,R>{
private:
    GLRLMFeatures<T,R> glrlm;
    GLRLMFeatures3DAVG<T,R> glrlm3D;
	vector<float> actualSpacing;
	string normGLRLM;
	float totalSum;
    boost::multi_array<float, 2> createGLRLMatrix3D(boost::multi_array<T, R> inputMatrix);
    void extractGLRLMData3D(vector<T> &glrlmData, GLRLMFeatures3D<T, R> glrlmFeatures);
    void fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<float, 2> &glrlMatrix, int directionX, int directionY, int directionZ);
    int maxRunLength;
    int totalNrVoxels;
public:
	GLRLMFeatures3D(){}
	~GLRLMFeatures3D(){}
    void calculateAllGLRLMFeatures3D(GLRLMFeatures3D<T,R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElem, vector<float> spacing, ConfigFile config);
    void writeCSVFileGLRLM3D(GLRLMFeatures3D<T,R> glrlmFeat, string outputFolder);
	void writeOneFileGLRLM3D(GLRLMFeatures3D<T, R> glrlmFeat, ConfigFile config, int &parameterSpaceNr);

};



/*!
The method getNeighbors3D stores all neighbor pairs for the desired angle and the actual input matrix in a vector
@param[in] inputMatrix: the original matrix of the VOI
@param[in] directionX
@param[in] directionY
@param[in] directionZ
The direction-parameters determine in which direction the run length is calculated
*/
template <class T, size_t R>
void GLRLMFeatures3D<T, R>::fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<float, 2> &glrlMatrix, int directionX, int directionY, int directionZ){
	float actGreyLevel = 0;
	float actElement = 0;
    int runLength=0;
    int maxRowNr = inputMatrix.shape()[0];
    int maxColNr = inputMatrix.shape()[1];
    int maxDepth = inputMatrix.shape()[2];
    //have a look at the image-matrix slide by slide (2D)
    //look for every grey level separately in every image slide
    for(int actGreyIndex=0; actGreyIndex<this->diffGreyLevels.size(); actGreyIndex++){
        //get the grey level we are interested at the moment
        actGreyLevel = this->diffGreyLevels[actGreyIndex];
        //now look at the image matrix slide by slide
        for(int depth = 0; depth < maxDepth; depth++){
            for(int row = 0; row < maxRowNr; row++){
                for(int column = 0; column < maxColNr; column++){
//                  //at the beginning the run length =0
                    runLength = 0;
                    //get the actual matrix element
                    int actDepth;
                    if(directionZ > 0){
                        actDepth = depth;
                    }
                    else{
                        actDepth = maxDepth - depth -1;
                    }
                    int actRow;
                    if(directionY < 0){
                        actRow = maxRowNr-row-1;
                    }
                    else{
                        actRow = row;
                    }
                    actElement = inputMatrix[actRow][column][actDepth];
                    //if the actual matrix element is the same as the actual gre level
                    if(actElement == actGreyLevel){
                        //set the run length to 1
                        runLength = 1;
                        //to avoid to take an element more than once, set the element to NAN
                        inputMatrix[actRow][column][actDepth] = NAN;
                        //now look at the matrix element in the actual direction (depends on the
                        //angle we are interested at the moment
                        int colValue = column + directionX;
                        int rowValue = actRow + directionY;
                        int depthValue = actDepth + directionZ;
                        //now have a look at the following elements in the desired direction
                        //stop as soon as we look at an element diifferent from our actual element
                        while(colValue < maxColNr && rowValue > -1 && rowValue < maxRowNr && colValue > -1 && depthValue < maxDepth && depthValue > -1 && inputMatrix[rowValue][colValue][depthValue] == actGreyLevel ){
                            //for every element we find, count the runLength
                            runLength += 1;
                            inputMatrix[rowValue][colValue][depthValue]=NAN;
                            //go further in the desired direction
                            colValue += directionX;
                            rowValue += directionY;

                            depthValue += directionZ;
                      }
                   }
                   //as soon as we cannot find an element in the desired direction, count one up in the desired
                   //position of the glrl-matrix
                    if(runLength > 0 && runLength < maxRunLength+1){
                        glrlMatrix[actGreyIndex][runLength-1]+=1;
                    }
                }
            }
        }
    }
}


/*!
The method createGLRLMatrix3D sums up all matrices of the different directions
@param[in] inputMatrix: the original matrix of the VOI
*/
template <class T, size_t R>
boost::multi_array<float, 2> GLRLMFeatures3D<T, R>::createGLRLMatrix3D(boost::multi_array<T,R> inputMatrix){
    typedef boost::multi_array<float, 2> glrlmat;

    int directionX;
    int directionY;
    int directionZ;;

    maxRunLength = glrlm.getMaxRunLength(inputMatrix);
    int sizeMatrix = this->diffGreyLevels.size();
    glrlmat GLRLMatrix(boost::extents[sizeMatrix][maxRunLength]);
    glrlmat sum(boost::extents[sizeMatrix][maxRunLength]);
	float weight;
    for(int i = 0; i < 13; i++){
        glrlm3D.getXYdirections3D(directionX, directionY, directionZ, i);
        glrlmat GLRLMatrix(boost::extents[sizeMatrix][maxRunLength]);
        fill3DMatrices(inputMatrix, GLRLMatrix, directionX, directionY, directionZ);
		weight = calculateWeight3D(directionX, directionY, directionZ, normGLRLM, actualSpacing);
		multSkalarMatrix(GLRLMatrix, weight);
		matrixSum(sum, GLRLMatrix);
    }

    return sum;
}

template <class T, size_t R>
void GLRLMFeatures3D<T, R>::calculateAllGLRLMFeatures3D(GLRLMFeatures3D<T,R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElement, vector<float> spacing, ConfigFile config){
    this->diffGreyLevels = diffGrey;
	actualSpacing = spacing;
	normGLRLM = config.normGLRLM;
	boost::multi_array<float,2> glrlMatrix = glrlmFeatures.createGLRLMatrix3D(inputMatrix);
	glrlmFeatures.getConfigValues(config);
    totalSum=glrlmFeatures.calculateTotalSum(glrlMatrix);
    vector<float> rowSums =glrlmFeatures.calculateRowSums(glrlMatrix);
    vector<float> colSums=glrlmFeatures.calculateColSums(glrlMatrix);
    glrlmFeatures.calculateShortRunEmphasis(rowSums, totalSum);
    glrlmFeatures.calculateLongRunEmphasis(rowSums, totalSum);
    glrlmFeatures.calculateLowGreyEmph(colSums, totalSum);
    glrlmFeatures.calculateHighGreyEmph(colSums, totalSum);
    glrlmFeatures.calculateShortRunLow(glrlMatrix, totalSum);
    glrlmFeatures.calculateShortRunHigh(glrlMatrix, totalSum);
    glrlmFeatures.calculateLongRunLowEmph(glrlMatrix, totalSum);
    glrlmFeatures.calculateLongRunHighEmph(glrlMatrix, totalSum);
    glrlmFeatures.calculateGreyNonUniformity(colSums, totalSum);
    glrlmFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
    glrlmFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
    glrlmFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
    glrlmFeatures.calculateRunPercentage3D(vectorMatrElement, totalSum, 13);
    boost::multi_array<float,2> probMatrix = glrlmFeatures.calculateProbMatrix(glrlMatrix, totalSum);
	float meanGrey = glrlmFeatures.calculateMeanProbGrey(probMatrix);

    glrlmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);

	float meanRun = glrlmFeatures.calculateMeanProbRun(probMatrix);
    glrlmFeatures.calculateRunLengthVar(probMatrix, meanRun);

    glrlmFeatures.calculateRunEntropy(probMatrix);

}

template <class T, size_t R>
void GLRLMFeatures3D<T, R>::extractGLRLMData3D(vector<T> &glrlmData, GLRLMFeatures3D<T, R> glrlmFeatures){

    glrlmData.push_back(glrlmFeatures.shortRunEmphasis);
    glrlmData.push_back(glrlmFeatures.longRunEmphasis);
    glrlmData.push_back(glrlmFeatures.lowGreyEmph);
    glrlmData.push_back(glrlmFeatures.highGreyEmph);
    glrlmData.push_back(glrlmFeatures.shortRunLow);
    glrlmData.push_back(glrlmFeatures.shortRunHigh);
    glrlmData.push_back(glrlmFeatures.longRunLowEmph);
    glrlmData.push_back(glrlmFeatures.longRunHighEmph);
    glrlmData.push_back(glrlmFeatures.greyNonUniformity);
    glrlmData.push_back(glrlmFeatures.greyNonUniformityNorm);
    glrlmData.push_back(glrlmFeatures.runLengthNonUniformity);   //check if I have to give the matrix
    glrlmData.push_back(glrlmFeatures.runLengthNonUniformityNorm);
    glrlmData.push_back(glrlmFeatures.runPercentage);
    glrlmData.push_back(glrlmFeatures.greyLevelVar);
    glrlmData.push_back(glrlmFeatures.runLengthVar);
    glrlmData.push_back(glrlmFeatures.runEntropy);

}

template <class T, size_t R>
void GLRLMFeatures3D<T, R>::writeCSVFileGLRLM3D(GLRLMFeatures3D<T,R> glrlmFeat, string outputFolder)
{
    string csvName = outputFolder + "_GLRLMFeatures3Dmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glrlmCSV;
    glrlmCSV.open (name);
    vector<string> features;
    glrlmFeat.defineGLRLMFeatures(features);

    vector<T> glrlmData;
    extractGLRLMData3D(glrlmData, glrlmFeat);
    for(int i = 0; i< glrlmData.size(); i++){
        glrlmCSV <<"GLRLMFeatures3Dmrg"<<","<< features[i] <<",";
        glrlmCSV << glrlmData[i];
        glrlmCSV << "\n";
    }
    glrlmCSV.close();
}


template <class T, size_t R>
void GLRLMFeatures3D<T, R>::writeOneFileGLRLM3D(GLRLMFeatures3D<T, R> glrlmFeat, ConfigFile config, int &parameterSpaceNr) {
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

	ofstream glrlmCSV;
	glrlmCSV.open(name, std::ios_base::app);
	vector<string> features;

	vector<T> glrlmData;
	extractGLRLMData3D(glrlmData, glrlmFeat);
	if (config.csvOutput == 1) {
		glrlmFeat.defineGLRLMFeatures(features);
		for (int i = 0; i < glrlmData.size(); i++) {
			glrlmCSV << "GLRLMFeatures3Dmrg" << "," << features[i] << ",";
			glrlmCSV << glrlmData[i];
			glrlmCSV << "\n";
		}
	}
	else if (config.ontologyOutput == 1) {
		glrlmFeat.defineGLRLMFeaturesOntology(features);
		string featParamSpaceTable = config.outputFolder + "/FeatureParameterSpace_table.csv";
		char * featParamSpaceTableName = new char[featParamSpaceTable.size() + 1];
		std::copy(featParamSpaceTable.begin(), featParamSpaceTable.end(), featParamSpaceTableName);
		featParamSpaceTableName[featParamSpaceTable.size()] = '\0';

		ofstream featSpaceTable;
		parameterSpaceNr += 1;
		string parameterSpaceName = "FeatureParameterSpace_" + std::to_string(parameterSpaceNr);
		featSpaceTable.open(featParamSpaceTableName, std::ios_base::app);
		featSpaceTable << parameterSpaceName << "," << "3Dmrg" << "," << config.imageSpaceName << "," << config.interpolationMethod << "\n";
		featSpaceTable.close();

		for (int i = 0; i < glrlmData.size(); i++) {
			glrlmCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			glrlmCSV << glrlmData[i] << "," << parameterSpaceName << "," << config.calculationSpaceName;
			glrlmCSV << "\n";
		}

	}
	glrlmCSV.close();
}

#endif // GLRLMFEATURES3D_H_INCLUDED
