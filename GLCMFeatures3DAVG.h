#ifndef GLCMFEATURES3DAVG_H_INCLUDED
#define GLCMFEATURES3DAVG_H_INCLUDED

#include "GLCMFeatures3DMRG.h"

/*! \file */

/*!
The class GLCMFeatures3DAVG herites from the class GLCMFeatures. \n
It considers 13 neighbors to calculate the cooccurrence features. \n
It calculates a GLCM matrix for every angle and extracts the feature from every matrix.
Then the mean value over all these features is calculated.\n
All feature calculations are defined in the class GLCMFeatures. \n
This class only contains the calculations of the merged matrix.
*/

template <class T,  size_t R>
class GLCMFeatures3DAVG : GLCMFeatures<T, R>{
private:
	GLCMFeatures<T, R> glcm;
	int sizeMatrix;
	typedef boost::multi_array<double, 2>  glcmat;
	void defineGLCMFeatures3DAVG(vector<string> &features);
	void extractGLCMData3D(vector<T> &glcmData, GLCMFeatures3DAVG<T, R> glcmFeatures);
	void fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &sum, int angle, int directionZ);
	boost::multi_array<double, 2> getMatrixSum(boost::multi_array<T, R> inputMatrix, int ang, int directionZ);
	//store different grey levels in vector
	vector<T> diffGreyLevels;
	vector<T> diagonalProbabilities;
	vector<T> crossProbabilities;
	vector<T> sumProbRows;
	vector<T> sumProbCols;
	T HX;
	T HXY;
	T HXY1;
	T HXY2;

public:
	GLCMFeatures3DAVG() {
	}
	~GLCMFeatures3DAVG() {
	}
	void writeCSVFileGLCM3DAVG(GLCMFeatures3DAVG<T, R> glcmFeat, string outputFolder);
	void writeOneFileGLCM3DAVG(GLCMFeatures3DAVG<T, R> glcmFeat, string outputFolder);
	void calculateAllGLCMFeatures3DAVG(GLCMFeatures3DAVG<T, R> &glcmFeat, boost::multi_array<T, R> inputMatrix, float maxIntensity);
   
};





/*!
In the method fill3DMatrices the matrix is filled for all directions
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: the GLCMatrix that will be filled with the corresponding values
@param[in]: direction that determines how the matrix will be filled
@param[in] directionZ: determines the direction in the depth
The function works analog to the 2D method
*/
template <class T, size_t R>
void GLCMFeatures3DAVG<T, R>::fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int angle, int directionZ) {
	std::vector<std::pair<T, T> > neighbours;
	//get the vector of the nieghbour pairs
	neighbours = glcm.getNeighbours3D(inputMatrix, angle, directionZ);
	std::pair<T, T> actNeighbour;
	//iterate over the neighbor-vector
	for (int neighbourNr = 0; neighbourNr<neighbours.size(); neighbourNr++) {
		actNeighbour = neighbours[neighbourNr];
		if (!std::isnan(actNeighbour.first) && !std::isnan(actNeighbour.second)) {
			glcMatrix[actNeighbour.first - 1][actNeighbour.second - 1] += 1;

		}
	}
}


/*!
In the method getMatrixSum calculates the sum from one GLCM matrix and its inverse
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: the GLCMatrix that will be filled with the corresponding values
@param[in]: direction that determines how the matrix will be filled
@param[in] directionZ: determines the direction in the depth
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLCMFeatures3DAVG<T, R>::getMatrixSum(boost::multi_array<T, R> inputMatrix, int ang, int directionZ) {

	glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);

	fill3DMatrices(inputMatrix, GLCMatrix, ang, directionZ);

	return GLCMatrix;
}

template <class T, size_t R>
void GLCMFeatures3DAVG<T, R>::calculateAllGLCMFeatures3DAVG(GLCMFeatures3DAVG<T, R> &glcmFeatures, boost::multi_array<T, R> inputMatrix, float maxIntensity) {

	int ang;
	int directionZ = 0;

	sizeMatrix = maxIntensity;

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


	glcmat sum(boost::extents[sizeMatrix][sizeMatrix]);
	glcmat inverseMatrix(boost::extents[sizeMatrix][sizeMatrix]);
	for (int i = 0; i<5; i++) {

		ang = 180 - i * 45;
		if (ang > 0) {
			for (int j = 0; j<3; j++) {

				directionZ = -1 + j;
				glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);

				fill3DMatrices(inputMatrix, GLCMatrix, ang, directionZ);

				inverse(GLCMatrix, inverseMatrix);

				sum = GLCMatrix;
				matrixSum(sum, inverseMatrix);

				//calculate the sum of all matrix elements (= the number of neighbor-pairs in the matrix
				double sumMatrElement = accumulate(sum.origin(), sum.origin() + sum.num_elements(), 0);
				//divide the whole matrix by the sum to obtain matrix elements representing the probabilities
				//of the occurence of a neighbor pair
				if (sumMatrElement != 0) {
					transform(sum.origin(), sum.origin() + sum.num_elements(),
						sum.origin(), bind2nd(std::divides<double>(), int(sumMatrElement)));
				}

				glcmFeatures.calculateJointMaximum(sum);
				sumJointMaximum += this->jointMaximum;
				glcmFeatures.calculateJointAverage(sum);
				sumJointAverage += this->jointAverage;
				glcmFeatures.calculateJointVariance(sum, this->jointAverage);
				sumJointVariance += this->jointVariance;
				glcmFeatures.calculateJointEntropy(sum);
				sumJointEntropy += this->jointEntropy;
				glcmFeatures.getDiagonalProbabilities(sum);
				glcmFeatures.getCrossProbabilities(sum);
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
				glcmFeatures.calculateAngSecMoment(sum);
				sumAngSecMoment += this->angSecMoment;
				glcmFeatures.calculateContrast(sum);
				sumContrast += this->contrast;
				glcmFeatures.calculateDissimilarity(sum);
				sumDissimilarity += this->dissimilarity;
				glcmFeatures.calculateInverseDiff(sum);
				sumInverseDiff += this->inverseDiff;
				glcmFeatures.calculateInverseDiffNorm(sum, this->inverseDiff);
				sumInverseDiffNorm += this->inverseDiffNorm;
				glcmFeatures.calculateInverseDiffMom(sum);
				sumInverseDiffMom += this->inverseDiffMom;
				glcmFeatures.calculateInverseDiffMomNorm(sum);
				sumInverseDiffMomNorm += this->inverseDiffMomNorm;
				glcmFeatures.calculateInverseVariance(sum);
				sumInverseVar += this->inverseVar;
				glcmFeatures.calculateCorrelation(sum);
				sumCorrelation += this->correlation;
				glcmFeatures.calculateAutoCorrelation(sum);
				sumAutoCorrelation += this->autoCorrelation;
				glcmFeatures.calculateClusterProminence(sum);
				sumClusterProminence += this->clusterProminence;
				glcmFeatures.calculateClusterShade(sum);
				sumClusterShade += this->clusterShade;
				glcmFeatures.calculateClusterTendency(sum);
				sumClusterTendency += this->clusterTendency;
				glcmFeatures.calculateFirstMCorrelation(sum);
				sumFirstMCorrelation += this->firstMCorrelation;
				glcmFeatures.calculateSecondMCorrelation(sum);
				sumSecondMCorrelation += this->secondMCorrelation;
			}
		}

		else {
			directionZ = 1;
			boost::multi_array<double, 2> GLCMatrix = glcmFeatures.getMatrixSum(inputMatrix, ang, directionZ);
			inverse(GLCMatrix, inverseMatrix);
			sum = GLCMatrix;
			matrixSum(sum, inverseMatrix);
			//calculate the sum of all matrix elements (= the number of neighbor-pairs in the matrix
			double sumMatrElement = accumulate(sum.origin(), sum.origin() + sum.num_elements(), 0);
			//divide the whole matrix by the sum to obtain matrix elements representing the probabilities
			//of the occurence of a neighbor pair
			transform(sum.origin(), sum.origin() + sum.num_elements(),
				sum.origin(), bind2nd(std::divides<double>(), int(sumMatrElement)));
			glcmFeatures.calculateJointMaximum(sum);
			sumJointMaximum += this->jointMaximum;
			glcmFeatures.calculateJointAverage(sum);
			sumJointAverage += this->jointAverage;
			glcmFeatures.calculateJointVariance(sum, this->jointAverage);
			sumJointVariance += this->jointVariance;
			glcmFeatures.calculateJointEntropy(sum);
			sumJointEntropy += this->jointEntropy;
			glcmFeatures.getDiagonalProbabilities(sum);
			glcmFeatures.getCrossProbabilities(sum);
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
			glcmFeatures.calculateAngSecMoment(sum);
			sumAngSecMoment += this->angSecMoment;
			glcmFeatures.calculateContrast(sum);
			sumContrast += this->contrast;
			glcmFeatures.calculateDissimilarity(sum);
			sumDissimilarity += this->dissimilarity;
			glcmFeatures.calculateInverseDiff(sum);
			sumInverseDiff += this->inverseDiff;
			glcmFeatures.calculateInverseDiffNorm(sum, this->inverseDiff);
			sumInverseDiffNorm += this->inverseDiffNorm;
			glcmFeatures.calculateInverseDiffMom(sum);
			sumInverseDiffMom += this->inverseDiffMom;
			glcmFeatures.calculateInverseDiffMomNorm(sum);
			sumInverseDiffMomNorm += this->inverseDiffMomNorm;
			glcmFeatures.calculateInverseVariance(sum);
			sumInverseVar += this->inverseVar;
			glcmFeatures.calculateCorrelation(sum);
			sumCorrelation += this->correlation;
			glcmFeatures.calculateAutoCorrelation(sum);
			sumAutoCorrelation += this->autoCorrelation;
			glcmFeatures.calculateClusterProminence(sum);
			sumClusterProminence += this->clusterProminence;
			glcmFeatures.calculateClusterShade(sum);
			sumClusterShade += this->clusterShade;
			glcmFeatures.calculateClusterTendency(sum);
			sumClusterTendency += this->clusterTendency;
			glcmFeatures.calculateFirstMCorrelation(sum);
			sumFirstMCorrelation += this->firstMCorrelation;
			glcmFeatures.calculateSecondMCorrelation(sum);
			sumSecondMCorrelation += this->secondMCorrelation;
		}
	}
	this->jointMaximum = sumJointMaximum / 13;
	this->jointAverage = sumJointAverage / 13;
	this->jointVariance = sumJointVariance / 13;
	this->jointEntropy = sumJointEntropy / 13;
	this->diffAverage = sumDiffAverage / 13;
	this->diffVariance = sumDiffVariance / 13;
	this->diffEntropy = sumDiffEntropy / 13;
	this->sumAverage = sumSumAverage / 13;
	this->sumVariance = sumSumVariance / 13;
	this->sumEntropy = sumSumEntropy / 13;
	this->angSecMoment = sumAngSecMoment / 13;
	this->contrast = sumContrast / 13;
	this->dissimilarity = sumDissimilarity / 13;
	this->inverseDiff = sumInverseDiff / 13;
	this->inverseDiffNorm = sumInverseDiffNorm / 13;
	this->inverseDiffMom = sumInverseDiffMom / 13;
	this->inverseDiffMomNorm = sumInverseDiffMomNorm / 13;
	this->inverseVar = sumInverseVar / 13;
	this->correlation = sumCorrelation / 13;
	this->autoCorrelation = sumAutoCorrelation / 13;
	this->clusterProminence = sumClusterProminence / 13;
	this->clusterShade = sumClusterShade / 13;
	this->clusterTendency = sumClusterTendency / 13;
	this->firstMCorrelation = sumFirstMCorrelation / 13;
	this->secondMCorrelation = sumSecondMCorrelation / 13;
}



template <class T, size_t R>
void GLCMFeatures3DAVG<T, R>::writeCSVFileGLCM3DAVG(GLCMFeatures3DAVG<T,R> glcmFeat, string outputFolder)
{
    string csvName = outputFolder + "_glcmFeatures3Davg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glcmCSV;
    glcmCSV.open (name);
    vector<string> features;
    defineGLCMFeatures3DAVG(features);

    vector<T> glcmData;
    extractGLCMData3D(glcmData, glcmFeat);
    for(int i = 0; i< glcmData.size(); i++){
        glcmCSV <<"glcmFeatures3Davg"<<","<< features[i] <<",";
        glcmCSV << glcmData[i];
        glcmCSV << "\n";
    }
    glcmCSV.close();
}

template <class T, size_t R>
void GLCMFeatures3DAVG<T, R>::writeOneFileGLCM3DAVG(GLCMFeatures3DAVG<T, R> glcmFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glcmCSV;
	glcmCSV.open(name, std::ios_base::app);
	vector<string> features;
	defineGLCMFeatures3DAVG(features);

	vector<T> glcmData;
	extractGLCMData3D(glcmData, glcmFeat);
	for (int i = 0; i< glcmData.size(); i++) {
		glcmCSV << "glcmFeatures3Davg" << "," << features[i] << ",";
		glcmCSV << glcmData[i];
		glcmCSV << "\n";
	}
	glcmCSV.close();
}

template <class T, size_t R>
void GLCMFeatures3DAVG<T, R>::defineGLCMFeatures3DAVG(vector<string> &features){
    features.push_back("joint maximun");
    features.push_back("joint average");
    features.push_back("joint variance");
    features.push_back("joint entropy");
    features.push_back("difference average");
    features.push_back("difference variance");
    features.push_back("difference entropy");
    features.push_back("sum average");
    features.push_back("sum variance");
    features.push_back("sum entropy");
    features.push_back("angular second moment");
    features.push_back("contrast");
    features.push_back("dissimilarity");
    features.push_back("inverse difference");
    features.push_back("inverse difference normalised");
    features.push_back("inverse difference moment");
    features.push_back("inverse difference moment normalised");
    features.push_back("inverse variance");
    features.push_back("correlation");
    features.push_back("autocorrelation");
    features.push_back("cluster tendency");
    features.push_back("cluster shade");
    features.push_back("cluster prominence");
    features.push_back("first measure of information correlation");
    features.push_back("second measure of information correlation");
}

template <class T, size_t R>
void GLCMFeatures3DAVG<T, R>::extractGLCMData3D(vector<T> &glcmData, GLCMFeatures3DAVG<T, R> glcmFeatures){
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

#endif // GLCMFEATURES3DAVG_H_INCLUDED
