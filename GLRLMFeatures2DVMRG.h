#ifndef GLRLMFEATURES2DVMRG_H_INCLUDED
#define GLRLMFEATURES2DVMRG_H_INCLUDED

#include "GLRLMFeatures.h"

/*! \file */

/*!
The class GLRLFeatures2DFullMerge herites from the class GLRLMFeatures. \n
It merges the matrices of all slices and calculates afterwards the features from the merged matrix. \n
All feature calculations are defined in the class GLRLMFeatures. \n
This class only contains the calculations of the merged matrix.
*/


template <class T,  size_t R=3>
class GLRLMFeatures2DVMRG : GLRLMFeatures<T,R>  {
    private:

        GLRLMFeatures<T,R> glrlm;

        double totalSum;

        int directionX;
        int directionY;

        int maxRunLength;
        vector<double> actualSpacing;
        string normGLRLM;
		vector<double> emphasisValues;
        boost::multi_array<double,2> createGLRLMatrixVMRG(boost::multi_array<T, R> inputMatrix, int depth, int ang);
        void extractGLRLMDataVMRG(vector<T> &glrlmData, GLRLMFeatures2DVMRG<T, R> glrlmFeatures);


    public:
		GLRLMFeatures2DVMRG(){
		}
		~GLRLMFeatures2DVMRG(){
		}
        float powRow;
        float powCol;

        void calculateAllGLRLMFeatures2DVMRG(GLRLMFeatures2DVMRG<T,R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElem, vector<double> spacing, ConfigFile config);
        void writeCSVFileGLRLM2DVMRG(GLRLMFeatures2DVMRG<T,R> glrlmFeat, string outputFolder);
        void writeOneFileGLRLM2DVMRG(GLRLMFeatures2DVMRG<T, R> glrlmFeat, string outputFolder);
		void fill2DMatrices2DVMRG(boost::multi_array<T, R> inputMatrix, boost::multi_array<double,2> &glrlMatrix, int depth, int ang);

};


/*!
In the method createGLRLMatrixFullMerge the GLRLM-matrix for the given angle and a given slice is calculated \n
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[in] : int angle: angle that is used to calculate the run length
@param[in] : int depth: number of the actual slice
@param[out]: GLCM-matrix
*/
template <class T, size_t R>
boost::multi_array<double,2> GLRLMFeatures2DVMRG<T, R>::createGLRLMatrixVMRG(boost::multi_array<T,R> inputMatrix, int depth, int ang){
    int sizeMatrix = this->diffGreyLevels.size();
    glrlm.getXYDirections(directionX, directionY, 180);

    boost::multi_array<double,2> GLRLMatrix(boost::extents[sizeMatrix][maxRunLength]);

    fill2DMatrices2DVMRG(inputMatrix, GLRLMatrix, depth, ang);

    return GLRLMatrix;

}

/*!
In the method fill2DMatrices2DFullMerge the matrix is filled for the given image slice and angle
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: the GLCMatrix that will be filled with the corresponding values
@param[in]: number of actual slice
@param[in]: direction that determines how the matrix will be filled

The function works as follows: \n
It determines for every element the run length of the specific angle. Is an element included in the calculation of
a run length, it is not considered anymore \n
As next step, the GLRLMatrix is filled: \n
for every grey level, the run length is determined \n
he corresponding element of the GLRLMatrix is increased
*/
template <class T, size_t R>
void GLRLMFeatures2DVMRG<T, R>::fill2DMatrices2DVMRG(boost::multi_array<T, R> inputMatrix, boost::multi_array<double,2> &glrlMatrix, int depth, int ang){

    T actElement = 0;
	int actGreyIndex;
    int runLength=0;
    int maxRowNr = inputMatrix.shape()[0];
    int maxColNr = inputMatrix.shape()[1];
    glrlm.getXYDirections(directionX, directionY, ang);
    //have a look at the image-matrix slide by slide (2D)
    //now look at the image matrix slide by slide
	for (int row = 0; row < maxRowNr; row++) {
		for (int column = 0; column < maxColNr; column++) {
			//at the beginning the run length =0
			runLength = 0;
			//get the actual matrix element
			actElement = inputMatrix[maxRowNr - row - 1][column][depth];
			//if the actual matrix element is the same as the actual gre level
			if (!isnan(actElement)) {
				//set the run length to 1
				runLength = 1;
				//to avoid to take an element more than once, set the element to NAN
				inputMatrix[maxRowNr - row - 1][column][depth] = NAN;
				actGreyIndex = glrlm.findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actElement);
				//now look at the matrix element in the actual direction (depends on the
				//angle we are interested at the moment
				int colValue = column + directionX;
				int rowValue = maxRowNr - 1 - (row + directionY);
				//now have a look at the following elements in the desired direction
				//stop as soon as we look at an element diifferent from our actual element
				while (colValue<maxColNr && rowValue>-1 && colValue > -1 && inputMatrix[rowValue][colValue][depth] == actElement) {
					//for every element we find, count the runLength
					runLength += 1;
					inputMatrix[rowValue][colValue][depth] = NAN;
					//go further in the desired direction
					colValue += 1 * directionX;
					rowValue -= 1 * directionY;
				}
			}
			//as soon as we cannot find an element in the desired direction, count one up in the desired
			//position of the glrl-matrix
			if (runLength > 0 && runLength < glrlMatrix.shape()[1] + 1) {
				glrlMatrix[actGreyIndex][runLength - 1] += 1;
			}
		}
	}
}

template <class T, size_t R>
void GLRLMFeatures2DVMRG<T, R>::calculateAllGLRLMFeatures2DVMRG(GLRLMFeatures2DVMRG<T,R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElem, vector<double> spacing, ConfigFile config){
    this->diffGreyLevels = diffGrey;

    actualSpacing = spacing;
    normGLRLM = config.normGLRLM;

    T sumShortRunEmphasis = 0;
    T sumLongRunEmphasis = 0;
    T sumLowGreyEmph = 0;
    T sumHighGreyEmph = 0;
    T sumShortRunLow = 0;
    T sumShortRunHigh = 0;
    T sumLongRunLowEmph = 0;
    T sumLongRunHighEmph = 0;
    T sumGreyNonUniformity = 0;
    T sumGreyNonUniformityNorm = 0;
    T sumRunLengthNonUniformity = 0;
    T sumRunLengthNonUniformityNorm = 0;
    T sumRunPercentage = 0;
    T sumGreyLevelVar = 0;
    T sumRunLengthVar = 0;
    T sumRunEntropy = 0;

    vector<double> rowSums;
    vector<double> colSums;

    double meanGrey;
    double meanRun;

    int totalDepth = inputMatrix.shape()[2];

    maxRunLength = glrlm.getMaxRunLength(inputMatrix);

    boost::multi_array<double, 2> sum(boost::extents[this->diffGreyLevels.size()][maxRunLength]);
    float weight;
    int ang;
    for(int depth = 0; depth < totalDepth; depth++){
        for(int i = 0; i < 4; i++){
            ang = 180-i*45;
            boost::multi_array<double,2> glrlMatrix = glrlmFeatures.createGLRLMatrixVMRG(inputMatrix, depth, ang);
            weight = calculateWeight2D(directionX, directionY, normGLRLM, actualSpacing);
            multSkalarMatrix(glrlMatrix, weight);
            matrixSum(sum, glrlMatrix);
        }
    }
	glrlmFeatures.getConfigValues(config);
	

    totalSum = glrlmFeatures.calculateTotalSum(sum);
    rowSums = glrlmFeatures.calculateRowSums(sum);

    colSums = glrlmFeatures.calculateColSums(sum);

    boost::multi_array<double,2> probMatrix = glrlmFeatures.calculateProbMatrix(sum, totalSum);
    meanGrey = glrlmFeatures.calculateMeanProbGrey(probMatrix);
    meanRun = glrlmFeatures.calculateMeanProbRun(probMatrix);
    glrlmFeatures.calculateShortRunEmphasis(rowSums, totalSum);
    glrlmFeatures.calculateLongRunEmphasis(rowSums, totalSum);
    glrlmFeatures.calculateLowGreyEmph(colSums, totalSum);
    glrlmFeatures.calculateHighGreyEmph(colSums, totalSum);
    glrlmFeatures.calculateShortRunLow(sum, totalSum);
    glrlmFeatures.calculateShortRunHigh(sum, totalSum);
    glrlmFeatures.calculateLongRunLowEmph(sum, totalSum);
    glrlmFeatures.calculateLongRunHighEmph(sum, totalSum);
    glrlmFeatures.calculateGreyNonUniformity(colSums, totalSum);
    glrlmFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
    glrlmFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
    glrlmFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
    glrlmFeatures.calculateRunPercentage3D(vectorMatrElem, totalSum, 4);
    glrlmFeatures.calculateGreyLevelVar(sum, meanGrey);
    glrlmFeatures.calculateRunLengthVar(sum, meanRun);
    glrlmFeatures.calculateRunEntropy(sum);
}


template <class T, size_t R>
void GLRLMFeatures2DVMRG<T, R>::writeCSVFileGLRLM2DVMRG(GLRLMFeatures2DVMRG<T, R> glrlmFeat2D, string outputFolder)
{
    string csvName = outputFolder + "_GLRLMFeatures2Dvmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glrlmCSV;
    glrlmCSV.open (name);
    vector<string> features;
    glrlm.defineGLRLMFeatures(features);

    vector<T> glrlmData;
    extractGLRLMDataVMRG(glrlmData, glrlmFeat2D);
    for(int i = 0; i< glrlmData.size(); i++){
        glrlmCSV << "GLRLMFeatures2Dvmrg"<<","<<features[i] <<",";
        glrlmCSV << glrlmData[i];
        glrlmCSV << "\n";
    }
    glrlmCSV.close();
}

template <class T, size_t R>
void GLRLMFeatures2DVMRG<T, R>::writeOneFileGLRLM2DVMRG(GLRLMFeatures2DVMRG<T, R> glrlmFeat2D, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glrlmCSV;
	glrlmCSV.open(name, std::ios_base::app);
	vector<string> features;
	glrlm.defineGLRLMFeatures(features);

	vector<T> glrlmData;
	extractGLRLMDataVMRG(glrlmData, glrlmFeat2D);
	for (int i = 0; i< glrlmData.size(); i++) {
		glrlmCSV << "GLRLMFeatures2Dvmrg" << "," << features[i] << ",";
		glrlmCSV << glrlmData[i];
		glrlmCSV << "\n";
	}
	glrlmCSV.close();
}


template <class T, size_t R>
void GLRLMFeatures2DVMRG<T, R>::extractGLRLMDataVMRG(vector<T> &glrlmData, GLRLMFeatures2DVMRG<T, R> glrlmFeatures){

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


#endif // GLRLMFEATURES2DVMRG_H_INCLUDED
