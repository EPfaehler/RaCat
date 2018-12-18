#ifndef GLCMFEATURES2DDMRG_H_INCLUDED
#define GLCMFEATURES2DDMRG_H_INCLUDED
/*! \file */


#include "GLCMFeatures2DAVG.h"

/*!
The class GLCMFeatures2DAVG inherits from the matrix GLCMFeatures. \n
It does not merge the matrices before feature calculation. \n
For every slice a GLCMatrix is calculated and from every of this matrices all features are extracted. \n
Then the average value of all features is calculated.
*/
template <class T, size_t R = 3>
class GLCMFeatures2DDMRG : GLCMFeatures2DAVG<T, R> {
private:
	GLCMFeatures<T, R> glcmComb;
	int sizeMatrix;

	void extractGLCMDataDMRG(vector<T> &glcmData, GLCMFeatures2DDMRG<T, R> glcmFeatures);
	//void fill2DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int depth, int angle);
	boost::multi_array<double, 2> calculateMatrix2DDMRG(boost::multi_array<T, R> inputMatrix, int angle);

	vector<T> diagonalProbabilities;
	vector<T> crossProbabilities;
	vector<T> sumProbRows;
	vector<T> sumProbCols;
	T HX;
	T HXY;
	T HXY1;
	T HXY2;

public:
	GLCMFeatures2DDMRG() {
	}
	~GLCMFeatures2DDMRG() {
	}
	void calculateAllGLCMFeatures2DDMRG(GLCMFeatures2DDMRG<T, R> &glcmFeat, boost::multi_array<T, R> inputMatrix, float maxIntensity);
	void writeCSVFileGLCM2DDMRG(GLCMFeatures2DDMRG<T, R> glcmFeat, string outputFolder);
	void writeOneFileGLCM2DDMRG(GLCMFeatures2DDMRG<T, R> glcmFeat, string outputFolder);
};


/*!
In the method calculateMatrix the GLCM-matrices for every direction are calculated, summed up and in the end the sum of this
matrices is divided by the sum of the elements (= nr. of neighbor pairs) to obtain a matrix which contains the probabilities
for the occurence of every neighbor pair.
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: GLCM-matrix
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLCMFeatures2DDMRG<T, R>::calculateMatrix2DDMRG(boost::multi_array<T, R> inputMatrix, int angle) {
	typedef boost::multi_array<double, 2> glcmat;
	GLCMFeatures2DAVG<T, R> glcmAVG;
	glcmat sum(boost::extents[sizeMatrix][sizeMatrix]);
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);
		glcmat inverseMatrix(boost::extents[sizeMatrix][sizeMatrix]);
		glcmAVG.fill2DMatrices(inputMatrix, GLCMatrix, depth, angle);
		inverse(GLCMatrix, inverseMatrix);
		matrixSum(sum, GLCMatrix);
		matrixSum(sum, inverseMatrix);
	}
	//calculate the sum of all matrix elements (= the number of neighbor-pairs in the matrix
	double sumMatrElement = accumulate(sum.origin(), sum.origin() + sum.num_elements(), 0);
	//divide the whole matrix by the sum to obtain matrix elements representing the probabilities
	//of the occurence of a neighbor pair
	if (sumMatrElement != 0) {
		transform(sum.origin(), sum.origin() + sum.num_elements(),
			sum.origin(), bind2nd(std::divides<double>(), int(sumMatrElement)));
	}
	return sum;
}



template <class T, size_t R>
void GLCMFeatures2DDMRG<T, R>::calculateAllGLCMFeatures2DDMRG(GLCMFeatures2DDMRG<T, R> &glcmFeatures, boost::multi_array<T, R> inputMatrix, float maxIntensity) {
	
	T sumJointMaximum = 0;
	T sumJointAverage = 0;
	T sumJointVariance = 0;
	T sumJointEntropy = 0;
	T sumDiffAverage = 0;
	T sumDiffVariance = 0;
	T sumDiffEntropy = 0;
	T sumSumAverage = 0;
	T sumSumVariance = 0;
	T sumSumEntropy = 0;
	T sumAngSecMoment = 0;
	T sumContrast = 0;
	T sumDissimilarity = 0;
	T sumInverseDiff = 0;
	T sumInverseDiffNorm = 0;
	T sumInverseDiffMom = 0;
	T sumInverseDiffMomNorm = 0;
	T sumInverseVar = 0;
	T sumAutoCorrelation = 0;
	T sumCorrelation = 0;
	T sumClusterTendency = 0;
	T sumClusterShade = 0;
	T sumClusterProminence = 0;
	T sumFirstMCorrelation = 0;
	T sumSecondMCorrelation = 0;

	int totalDepth = inputMatrix.shape()[2];

	double sumMatrElement;
	int ang;
	sizeMatrix = maxIntensity;
	for (int i = 0; i < 4; i++) {
		ang = 180 - i * 45;
		boost::multi_array<double, 2> GLCMatrix = glcmFeatures.calculateMatrix2DDMRG(inputMatrix, ang);

		glcmFeatures.calculateJointMaximum(GLCMatrix);
		sumJointMaximum += this->jointMaximum;
		glcmFeatures.calculateJointAverage(GLCMatrix);
		sumJointAverage += this->jointAverage;
		glcmFeatures.calculateJointVariance(GLCMatrix, this->jointAverage);
		sumJointVariance += this->jointVariance;
		glcmFeatures.calculateJointEntropy(GLCMatrix);
		sumJointEntropy += this->jointEntropy;
		glcmFeatures.getDiagonalProbabilities(GLCMatrix);
		glcmFeatures.getCrossProbabilities(GLCMatrix);
		glcmFeatures.calculateDiffAverage();
		sumDiffAverage += this->diffAverage;
		glcmFeatures.calculateDiffVariance(this->diffAverage);
		sumDiffVariance += this->diffVariance;
		glcmFeatures.calculateDiffEntropy();
		sumDiffEntropy += this->diffEntropy;
		glcmFeatures.calculateSumAverage();
		sumSumAverage += this->sumAverage;
		glcmFeatures.calculateSumVariance(this->sumAverage);
		sumSumVariance += this->sumVariance;
		glcmFeatures.calculateSumEntropy();
		sumSumEntropy += this->sumEntropy;
		glcmFeatures.calculateAngSecMoment(GLCMatrix);
		sumAngSecMoment += this->angSecMoment;
		glcmFeatures.calculateContrast(GLCMatrix);
		sumContrast += this->contrast;
		glcmFeatures.calculateDissimilarity(GLCMatrix);
		sumDissimilarity += this->dissimilarity;
		glcmFeatures.calculateInverseDiff(GLCMatrix);
		sumInverseDiff += this->inverseDiff;
		glcmFeatures.calculateInverseDiffNorm(GLCMatrix, this->inverseDiff);
		sumInverseDiffNorm += this->inverseDiffNorm;
		glcmFeatures.calculateInverseDiffMom(GLCMatrix);
		sumInverseDiffMom += this->inverseDiffMom;
		glcmFeatures.calculateInverseDiffMomNorm(GLCMatrix);
		sumInverseDiffMomNorm += this->inverseDiffMomNorm;
		glcmFeatures.calculateInverseVariance(GLCMatrix);
		sumInverseVar += this->inverseVar;
		glcmFeatures.calculateCorrelation(GLCMatrix);
		sumCorrelation += this->correlation;
		glcmFeatures.calculateAutoCorrelation(GLCMatrix);
		sumAutoCorrelation += this->autoCorrelation;
		glcmFeatures.calculateClusterProminence(GLCMatrix);
		sumClusterProminence += this->clusterProminence;
		glcmFeatures.calculateClusterShade(GLCMatrix);
		sumClusterShade += this->clusterShade;
		glcmFeatures.calculateClusterTendency(GLCMatrix);
		sumClusterTendency += this->clusterTendency;
		glcmFeatures.calculateFirstMCorrelation(GLCMatrix);
		sumFirstMCorrelation += this->firstMCorrelation;
		glcmFeatures.calculateSecondMCorrelation(GLCMatrix);
		sumSecondMCorrelation += this->secondMCorrelation;
	}

	this->jointMaximum = sumJointMaximum / 4;
	this->jointAverage = sumJointAverage / 4;
	this->jointVariance = sumJointVariance / 4;
	this->jointEntropy = sumJointEntropy / 4;
	this->diffAverage = sumDiffAverage / 4;
	this->diffVariance = sumDiffVariance / 4;
	this->diffEntropy = sumDiffEntropy / 4;
	this->sumAverage = sumSumAverage / 4;
	this->sumVariance = sumSumVariance / 4;
	this->sumEntropy = sumSumEntropy / 4;
	this->angSecMoment = sumAngSecMoment / 4;
	this->contrast = sumContrast / 4;
	this->dissimilarity = sumDissimilarity / 4;
	this->inverseDiff = sumInverseDiff / 4;
	this->inverseDiffNorm = sumInverseDiffNorm / 4;
	this->inverseDiffMom = sumInverseDiffMom / 4;
	this->inverseDiffMomNorm = sumInverseDiffMomNorm / 4;
	this->inverseVar = sumInverseVar / 4;
	this->correlation = sumCorrelation / 4;
	this->autoCorrelation = sumAutoCorrelation / 4;
	this->clusterProminence = sumClusterProminence / 4;
	this->clusterShade = sumClusterShade / 4;
	this->clusterTendency = sumClusterTendency / 4;
	this->firstMCorrelation = sumFirstMCorrelation / 4;
	this->secondMCorrelation = sumSecondMCorrelation / 4;

}

template <class T, size_t R>
void GLCMFeatures2DDMRG<T, R>::writeCSVFileGLCM2DDMRG(GLCMFeatures2DDMRG<T, R> glcmFeat, string outputFolder)
{
	string csvName = outputFolder + "_glcmFeatures2DDmrg.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glcmCSV;
	glcmCSV.open(name);
	vector<string> features;
	glcmComb.defineGLCMFeatures(features);

	vector<T> glcmData;
	extractGLCMDataDMRG(glcmData, glcmFeat);
	for (int i = 0; i< glcmData.size(); i++) {
		glcmCSV << "glcmFeatures2DDmrg" << "," << features[i] << ",";
		glcmCSV << glcmData[i];
		glcmCSV << "\n";
	}
	glcmCSV.close();
}

template <class T, size_t R>
void GLCMFeatures2DDMRG<T, R>::writeOneFileGLCM2DDMRG(GLCMFeatures2DDMRG<T, R> glcmFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glcmCSV;
	glcmCSV.open(name, std::ios_base::app);
	vector<string> features;
	glcmComb.defineGLCMFeatures(features);

	vector<T> glcmData;
	extractGLCMDataDMRG(glcmData, glcmFeat);
	for (int i = 0; i< glcmData.size(); i++) {
		glcmCSV << "glcmFeatures2DDmrg" << "," << features[i] << ",";
		glcmCSV << glcmData[i];
		glcmCSV << "\n";
	}
	glcmCSV.close();

}


template <class T, size_t R>
void GLCMFeatures2DDMRG<T, R>::extractGLCMDataDMRG(vector<T> &glcmData, GLCMFeatures2DDMRG<T, R> glcmFeatures) {

	glcmData.push_back(glcmFeatures.jointMaximum);
	glcmData.push_back(glcmFeatures.jointAverage);
	glcmData.push_back(glcmFeatures.jointVariance);
	glcmData.push_back(glcmFeatures.jointEntropy);
	glcmData.push_back(glcmFeatures.diffAverage);
	glcmData.push_back(glcmFeatures.diffVariance);
	glcmData.push_back(glcmFeatures.diffEntropy);
	glcmData.push_back(glcmFeatures.sumAverage);
	glcmData.push_back(glcmFeatures.sumVariance);
	glcmData.push_back(glcmFeatures.sumEntropy);
	glcmData.push_back(glcmFeatures.angSecMoment);
	glcmData.push_back(glcmFeatures.contrast);
	glcmData.push_back(glcmFeatures.dissimilarity);
	glcmData.push_back(glcmFeatures.inverseDiff);
	glcmData.push_back(glcmFeatures.inverseDiffNorm);
	glcmData.push_back(glcmFeatures.inverseDiffMom);
	glcmData.push_back(glcmFeatures.inverseDiffMomNorm);
	glcmData.push_back(glcmFeatures.inverseVar);
	glcmData.push_back(glcmFeatures.correlation);
	glcmData.push_back(glcmFeatures.autoCorrelation);
	glcmData.push_back(glcmFeatures.clusterTendency);
	glcmData.push_back(glcmFeatures.clusterShade);
	glcmData.push_back(glcmFeatures.clusterProminence);
	glcmData.push_back(glcmFeatures.firstMCorrelation);
	glcmData.push_back(glcmFeatures.secondMCorrelation);
}






#endif // GLCMFEATURES2DDMRG_H_INCLUDED