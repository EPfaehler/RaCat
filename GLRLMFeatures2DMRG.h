#ifndef GLRLMFEATURES2DMRG_H_INCLUDED
#define GLRLMFEATURES2DMRG_H_INCLUDED

#include "GLRLMFeatures.h"

/*! \file */

/*!
The class GLRLMFeatures2DMRG herites from the class GLRLMFeatures. \n
It merges the matrices of every slice separately and calculates afterwards the features from the merged matrix. \n
Afterwards the mean value of the features for every slice is calculated. \n
The difference between this class and the other GLRLMFeature-classes is only the type of merging of the matrix. \n
All feature calculations are defined in the class GLRLMFeatures. \n
This class only contains the calculations of the merged matrix.
*/

template <class T,  size_t R=3>
class GLRLMFeatures2DMRG : GLRLMFeatures<T,R>  {
    private:

        GLRLMFeatures<T,R> glrlm;

        double totalSum;

        typedef boost::multi_array<double,2> glrlmMat;

        int directionX;
        int directionY;

        int maxRunLength;

		vector<double> actualSpacing;
		string normGLRLM;

        boost::multi_array<double,2> createGLRLMatrixMRG(boost::multi_array<T, R> inputMatrix, int depth);
        void extractGLRLMDataMRG(vector<T> &glrlmData, GLRLMFeatures2DMRG<T, R> glrlmFeatures);
        void fill2DMatrices2DMRG(boost::multi_array<T, R> inputMatrix, boost::multi_array<double,2> &glrlMatrix, int depth, int ang);


    public:
		GLRLMFeatures2DMRG() {
		}
		~GLRLMFeatures2DMRG() {
		}
        void calculateAllGLRLMFeatures2DMRG(GLRLMFeatures2DMRG<T,R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<double> spacing, ConfigFile config);
        void writeCSVFileGLRLM2DMRG(GLRLMFeatures2DMRG<T,R> glrlmFeat, string outputFolder);
		void writeOneFileGLRLM2DMRG(GLRLMFeatures2DMRG<T, R> glrlmFeat, string outputFolder);

};


/*!
In the method createGLRLMatrixMRG the GLRLM-matrix for given slice is calculated \n
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[in] : int depth: number of the actual slice
@param[out]: GLCM-matrix
*/
template <class T, size_t R>
boost::multi_array<double,2> GLRLMFeatures2DMRG<T, R>::createGLRLMatrixMRG(boost::multi_array<T,R> inputMatrix, int depth){

    int sizeMatrix = this->diffGreyLevels.size();
    glrlmMat sum(boost::extents[sizeMatrix][maxRunLength]);
    int ang;
	float weight;
    for(int i = 0; i < 4; i++){
        ang = 180-i*45;
        glrlmMat GLRLMatrix(boost::extents[sizeMatrix][maxRunLength]);
        fill2DMatrices2DMRG(inputMatrix, GLRLMatrix, depth, ang);
		weight = calculateWeight2D(directionX, directionY, normGLRLM, actualSpacing);
		multSkalarMatrix(GLRLMatrix, weight);
        matrixSum(sum, GLRLMatrix);
    }

    return sum;
}

/*!
In the method fill2DMatrices2DMRG the matrix is filled for the given image slice and angle
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: the GLCMatrix that will be filled with the corresponding values
@param[in]: number of actual slice
@param[in]: direction that determines how the matrix will be filled

The function works analog to the function in GLRLMFeatures2DFullMerge
*/
template <class T, size_t R>
void GLRLMFeatures2DMRG<T, R>::fill2DMatrices2DMRG(boost::multi_array<T, R> inputMatrix, boost::multi_array<double,2> &glrlMatrix, int depth, int ang){

    T actElement = 0;
	int actGreyIndex;
    int runLength=0;
    int maxRowNr = inputMatrix.shape()[0];
    int maxColNr = inputMatrix.shape()[1];
    glrlm.getXYDirections(directionX, directionY, ang);
    //have a look at the image-matrix slide by slide (2D)
    //look for every grey level separately in every image slide
    //now look at the image matrix slide by slide
	for(int row = 0; row<maxRowNr;row++){
        for(int column = 0; column<maxColNr; column++){
			//at the beginning the run length =0
            runLength =0;
           //get the actual matrix element
			actElement = inputMatrix[maxRowNr-row-1][column][depth];
            //if the actual matrix element is the same as the actual gre level
			if (!isnan(actElement) ){
				//set the run length to 1
				runLength = 1;
				actGreyIndex = glrlm.findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actElement);
				//to avoid to take an element more than once, set the element to NAN
				inputMatrix[maxRowNr - row - 1][column][depth] = NAN;
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
            if(runLength > 0 && runLength < glrlMatrix.shape()[1] + 1){
				glrlMatrix[actGreyIndex][runLength-1] += 1;
            }
       }  
    }
}


template <class T, size_t R>
void GLRLMFeatures2DMRG<T, R>::calculateAllGLRLMFeatures2DMRG(GLRLMFeatures2DMRG<T,R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<double> spacing, ConfigFile config){

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
	glrlmFeatures.getConfigValues(config);
    for(int depth = 0; depth < totalDepth; depth++){
        boost::multi_array<double,2> glrlMatrix = createGLRLMatrixMRG(inputMatrix, depth);

        totalSum = glrlmFeatures.calculateTotalSum(glrlMatrix);
        rowSums = glrlmFeatures.calculateRowSums(glrlMatrix);
        colSums = glrlmFeatures.calculateColSums(glrlMatrix);

        boost::multi_array<double,2> probMatrix = glrlmFeatures.calculateProbMatrix(glrlMatrix, totalSum);
        meanGrey = glrlmFeatures.calculateMeanProbGrey(probMatrix);
        meanRun = glrlmFeatures.calculateMeanProbRun(probMatrix);

		
		glrlmFeatures.calculateShortRunEmphasis(rowSums, totalSum);
		
		sumShortRunEmphasis += this->shortRunEmphasis;
		
        glrlmFeatures.calculateLongRunEmphasis(rowSums, totalSum);
		sumLongRunEmphasis += this->longRunEmphasis;
		
        glrlmFeatures.calculateLowGreyEmph(colSums, totalSum);
		sumLowGreyEmph += this->lowGreyEmph;
		glrlmFeatures.calculateHighGreyEmph(colSums, totalSum);
		
		sumHighGreyEmph += this->highGreyEmph;
		
        glrlmFeatures.calculateShortRunLow(glrlMatrix, totalSum);
		sumShortRunLow += this->shortRunLow;
		glrlmFeatures.calculateShortRunHigh(glrlMatrix, totalSum);
		
		sumShortRunHigh += this->shortRunHigh;
		
        glrlmFeatures.calculateLongRunLowEmph(glrlMatrix, totalSum);
		sumLongRunLowEmph += this->longRunLowEmph;
		
        glrlmFeatures.calculateLongRunHighEmph(glrlMatrix, totalSum);
		
		sumLongRunHighEmph += this->longRunHighEmph;
		
        glrlmFeatures.calculateGreyNonUniformity(colSums, totalSum);
		
		sumGreyNonUniformity += this->greyNonUniformity;
		glrlmFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
		sumGreyNonUniformityNorm += this->greyNonUniformityNorm;
		
        glrlmFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
		sumRunLengthNonUniformity += this->runLengthNonUniformity;
		
        glrlmFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
		sumRunLengthNonUniformityNorm += this->runLengthNonUniformityNorm;
		glrlmFeatures.calculateRunPercentage(inputMatrix, depth, totalSum, 4);
		sumRunPercentage += this->runPercentage;
		
        glrlmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);
		
		sumGreyLevelVar += this->greyLevelVar;
        glrlmFeatures.calculateRunLengthVar(probMatrix, meanRun);
		sumRunLengthVar += this->runLengthVar;
	    glrlmFeatures.calculateRunEntropy(probMatrix);
		sumRunEntropy += this->runEntropy;
		

    }

    this->shortRunEmphasis = sumShortRunEmphasis/ totalDepth;
    this->longRunEmphasis = sumLongRunEmphasis/ totalDepth;
    this->lowGreyEmph = sumLowGreyEmph/ totalDepth;
    this->highGreyEmph = sumHighGreyEmph/ totalDepth;
    this->shortRunLow = sumShortRunLow/ totalDepth;
    this->shortRunHigh = sumShortRunHigh/ totalDepth;
    this->longRunLowEmph = sumLongRunLowEmph/ totalDepth;
    this->longRunHighEmph = sumLongRunHighEmph/ totalDepth;
    this->greyNonUniformity = sumGreyNonUniformity/ totalDepth;
    this->greyNonUniformityNorm = sumGreyNonUniformityNorm/ totalDepth;
    this->runLengthNonUniformity = sumRunLengthNonUniformity/ totalDepth;
    this->runLengthNonUniformityNorm = sumRunLengthNonUniformityNorm/ totalDepth;
    this->runPercentage = sumRunPercentage/ totalDepth;

    this->greyLevelVar = sumGreyLevelVar/ totalDepth;
    this->runLengthVar = sumRunLengthVar/ totalDepth;
    this->runEntropy = sumRunEntropy/ totalDepth;

}



template <class T, size_t R>
void GLRLMFeatures2DMRG<T, R>::writeCSVFileGLRLM2DMRG(GLRLMFeatures2DMRG<T, R> glrlmFeat, string outputFolder)
{
    string csvName = outputFolder + "_GLRLMFeatures2DWmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glrlmCSV;
    glrlmCSV.open (name);
    vector<string> features;
    glrlm.defineGLRLMFeatures(features);

    vector<T> glrlmData;
    extractGLRLMDataMRG(glrlmData, glrlmFeat);
    for(int i = 0; i< glrlmData.size(); i++){
        glrlmCSV <<"GLRLMFeatures2DWmrg"<<","<< features[i] <<",";
        glrlmCSV << glrlmData[i];
        glrlmCSV << "\n";
    }
    glrlmCSV.close();
}

template <class T, size_t R>
void GLRLMFeatures2DMRG<T, R>::writeOneFileGLRLM2DMRG(GLRLMFeatures2DMRG<T, R> glrlmFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glrlmCSV;
	glrlmCSV.open(name, std::ios_base::app);
	vector<string> features;
	glrlm.defineGLRLMFeatures(features);

	vector<T> glrlmData;
	extractGLRLMDataMRG(glrlmData, glrlmFeat);
	for (int i = 0; i< glrlmData.size(); i++) {
		glrlmCSV << "GLRLMFeatures2DWmrg" << "," << features[i] << ",";
		glrlmCSV << glrlmData[i];
		glrlmCSV << "\n";
	}
	glrlmCSV.close();
}


template <class T, size_t R>
void GLRLMFeatures2DMRG<T, R>::extractGLRLMDataMRG(vector<T> &glrlmData, GLRLMFeatures2DMRG<T, R> glrlmFeatures){

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




#endif // GLRLMFEATURES2DMRG_H_INCLUDED
