#ifndef GLRLMFEATURES2DDMRG_H_INCLUDED
#define GLRLMFEATURES2DDMRG_H_INCLUDED


#include "GLRLMFeatures2DVMRG.h"
#include "GLRLMFeatures2DAVG.h"
/*! \file */
/*!
The class GLCMFeatures2DWOMerge inherits from the matrix GLCMFeatures. \n
It does not merge the matrices before feature calculation. \n
For every slice a GLCMatrix is calculated and from every of this matrices all features are extracted. \n
Then the average value of all features is calculated.
*/
template <class T, size_t R = 3>
class GLRLMFEATURES2DDMRG : GLRLMFeatures<T, R> {

private:

	GLRLMFeatures<T, R> glrlm;
	GLRLMFeatures2DVMRG<T, R> glrlm2DFullMerge;

	double totalSum;

	typedef boost::multi_array<double, 2> glrlmMat;

	int directionX;
	int directionY;

	int maxRunLength;
	string normGLRLM;
	vector<double> actualSpacing;
	boost::multi_array<double, 2> createGLRLMatrix2DDMRG(boost::multi_array<T, R> inputMatrix, int ang);
	void extractGLRLMData2DDMRG(vector<T> &glrlmData, GLRLMFEATURES2DDMRG<T, R> glrlmFeatures);
	void fill2DMatrices2DDMRG(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glrlMatrix, int depth, int ang);
	void calculateRunPercentage2DDMRG(boost::multi_array<T, R> inputMatrix, int depth, double totalSum, int nrNeighbor);
public:
	GLRLMFEATURES2DDMRG() {
	}
	~GLRLMFEATURES2DDMRG() {
	}
	void calculateAllGLRLMFeatures2DDMRG(GLRLMFEATURES2DDMRG<T, R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffFGrey, vector<double> spacing, ConfigFile config);
	void writeCSVFileGLRLM2DDMRG(GLRLMFEATURES2DDMRG<T, R> glrlmFeat, string outputFolder);
	void writeOneFileGLRLM2DDMRG(GLRLMFEATURES2DDMRG<T, R> glrlmFeat, string outputFolder);

};

/*!
In the method createGLRLMatrixW=Merge the GLRLM-matrix for given slice is calculated \n
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[in] : int depth: number of the actual slice
@param[in] : int angle: angle
@param[out]: GLCM-matrix
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLRLMFEATURES2DDMRG<T, R>::createGLRLMatrix2DDMRG(boost::multi_array<T, R> inputMatrix, int ang) {
	GLRLMFeatures2DAVG<T, R> glrlm2Davg;
	int sizeMatrix = this->diffGreyLevels.size();
	glrlmMat sum(boost::extents[sizeMatrix][this->maxRunLength]);
	float weight;
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		
		glrlmMat GLRLMatrix(boost::extents[sizeMatrix][maxRunLength]);
		glrlm2Davg.fill2DMatrices2DAVG(inputMatrix, GLRLMatrix, this->diffGreyLevels,depth, ang);
		

		weight = calculateWeight2D(directionX, directionY, normGLRLM, actualSpacing);
		multSkalarMatrix(GLRLMatrix, weight);
		matrixSum(sum, GLRLMatrix);
	}

	return sum;
}
	

template <class T, size_t R>
void GLRLMFEATURES2DDMRG<T, R>::fill2DMatrices2DDMRG(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glrlMatrix, int depth, int ang) {

	T actGreyLevel = 0;
	T actElement = 0;
	int runLength = 0;
	int maxRowNr = inputMatrix.shape()[0];
	int maxColNr = inputMatrix.shape()[1];
	glrlm.getXYDirections(directionX, directionY, ang);
	//have a look at the image-matrix slide by slide (2D)
	//look for every grey level separately in every image slide
	int actGreyIndex;
	//get the grey level we are interested at the moment

	for (int row = 0; row<maxRowNr; row++) {
		for (int column = 0; column<maxColNr; column++) {
			//         //at the beginning the run length =0
			runLength = 0;
			//get the actual matrix element
			actElement = inputMatrix[maxRowNr - row - 1][column][depth];
			actGreyIndex = glrlm.findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actElement);
			//if the actual matrix element is the same as the actual gre level
			if (!std::isnan(actElement)) {
				//set the run length to 1
				runLength = 1;
				//to avoid to take an element more than once, set the element to NAN
				inputMatrix[maxRowNr - row - 1][column][depth] = NAN;
				////          //now look at the matrix element in the actual direction (depends on the
				//angle we are interested at the moment
				int colValue = column + directionX;
				int rowValue = maxRowNr - 1 - (row + directionY);
				//now have a look at the following elements in the desired direction
				//stop as soon as we look at an element diifferent from our actual element
				while (colValue<maxColNr && rowValue>-1 && colValue>-1 && inputMatrix[rowValue][colValue][depth] == actElement) {
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
void GLRLMFEATURES2DDMRG<T, R>::calculateRunPercentage2DDMRG(boost::multi_array<T, R> inputMatrix, int depth, double totalSum, int nrNeighbor) {
	int totalNrVoxels = 0;
	for (int depth = 0; depth < inputMatrix.shape()[2];  depth++) {
		for (int row = 0; row < inputMatrix.shape()[0]; row++) {
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				if (!std::isnan(inputMatrix[row][col][depth])) {
					totalNrVoxels++;
				}
			}
		}
	}
	
	
	if ((totalNrVoxels)*nrNeighbor != 0) {
		runPercentage = totalSum / ((totalNrVoxels)*nrNeighbor);
	}
	else {
		runPercentage = 0;
	}
}

template <class T, size_t R>
void GLRLMFEATURES2DDMRG<T, R>::calculateAllGLRLMFeatures2DDMRG(GLRLMFEATURES2DDMRG<T, R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<double> spacing, ConfigFile config) {
	this->diffGreyLevels = diffGrey;
	normGLRLM = config.normGLRLM;
	actualSpacing = spacing;
	glrlmFeatures.getConfigValues(config);
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

	int ang;
	for (int i = 0; i < 4; i++) {
		ang = 180 - i * 45;
		boost::multi_array<double, 2> glrlMatrix = createGLRLMatrix2DDMRG(inputMatrix, ang);

		totalSum = glrlmFeatures.calculateTotalSum(glrlMatrix);
		rowSums = glrlmFeatures.calculateRowSums(glrlMatrix);
		colSums = glrlmFeatures.calculateColSums(glrlMatrix);

		boost::multi_array<double, 2> probMatrix = glrlmFeatures.calculateProbMatrix(glrlMatrix, totalSum);
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
		glrlmFeatures.calculateRunPercentage2DDMRG(inputMatrix, 0, totalSum, 1);
		sumRunPercentage += this->runPercentage;
		glrlmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);

		sumGreyLevelVar += this->greyLevelVar;

		glrlmFeatures.calculateRunLengthVar(probMatrix, meanRun);
		sumRunLengthVar += this->runLengthVar;

		glrlmFeatures.calculateRunEntropy(probMatrix);
		sumRunEntropy += this->runEntropy;


	}


	this->shortRunEmphasis = sumShortRunEmphasis / 4;
	this->longRunEmphasis = sumLongRunEmphasis /  4;
	this->lowGreyEmph = sumLowGreyEmph / 4;
	this->highGreyEmph = sumHighGreyEmph / 4;
	this->shortRunLow = sumShortRunLow / 4;
	this->shortRunHigh = sumShortRunHigh / 4;
	this->longRunLowEmph = sumLongRunLowEmph / 4;
	this->longRunHighEmph = sumLongRunHighEmph / 4;
	this->greyNonUniformity = sumGreyNonUniformity / 4;
	this->greyNonUniformityNorm = sumGreyNonUniformityNorm / 4;
	this->runLengthNonUniformity = sumRunLengthNonUniformity / 4;
	this->runLengthNonUniformityNorm = sumRunLengthNonUniformityNorm / 4;
	this->runPercentage = sumRunPercentage / 4;

	this->greyLevelVar = sumGreyLevelVar / 4;
	this->runLengthVar = sumRunLengthVar / 4;
	this->runEntropy = sumRunEntropy / 4;
}



template <class T, size_t R>
void GLRLMFEATURES2DDMRG<T, R>::writeCSVFileGLRLM2DDMRG(GLRLMFEATURES2DDMRG<T, R> glrlmFeat, string outputFolder)
{
	string csvName = outputFolder + "_GLRLMFeatures2DDmrg.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glrlmCSV;
	glrlmCSV.open(name);
	vector<string> features;
	glrlm.defineGLRLMFeatures(features);

	vector<T> glrlmData;
	extractGLRLMData2DDMRG(glrlmData, glrlmFeat);
	for (int i = 0; i< glrlmData.size(); i++) {
		glrlmCSV << "GLRLMFeatures2DDmrg" << "," << features[i] << ",";
		glrlmCSV << glrlmData[i];
		glrlmCSV << "\n";
	}
	glrlmCSV.close();
}

template <class T, size_t R>
void GLRLMFEATURES2DDMRG<T, R>::writeOneFileGLRLM2DDMRG(GLRLMFEATURES2DDMRG<T, R> glrlmFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glrlmCSV;
	glrlmCSV.open(name, std::ios_base::app);
	vector<string> features;
	glrlm.defineGLRLMFeatures(features);

	vector<T> glrlmData;
	extractGLRLMData2DDMRG(glrlmData, glrlmFeat);
	for (int i = 0; i< glrlmData.size(); i++) {
		glrlmCSV << "GLRLMFeatures2DDmrg" << "," << features[i] << ",";
		glrlmCSV << glrlmData[i];
		glrlmCSV << "\n";
	}
	glrlmCSV.close();
}


template <class T, size_t R>
void GLRLMFEATURES2DDMRG<T, R>::extractGLRLMData2DDMRG(vector<T> &glrlmData, GLRLMFEATURES2DDMRG<T, R> glrlmFeatures) {

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





#endif // GLRLMFEATURES2DDMRG_H_INCLUDED
