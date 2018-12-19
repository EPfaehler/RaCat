#ifndef GLCMFEATURES3DMRG_H_INCLUDED
#define GLCMFEATURES3DMRG_H_INCLUDED

#include "GLCMFeatures.h"

/*! \file */

/*!
The class GLCMFeatures3DMRG herites from the class GLCMFeatures. \n
It considers 13 neighbors to calculate the cooccurrence features. \n
It calculates the feature values for every angle and calculates then
the mean value of the features. \n
All feature calculations are defined in the class GLCMFeatures. \n
This class only contains the calculations of the merged matrix.
*/
template <class T,  size_t R>
class GLCMFeatures3DMRG : GLCMFeatures<T, R>{
private:
	typedef boost::multi_array<double, 2>  glcmat;

	string normGLCM;
	vector<double> actualSpacing;
	void defineGLCMFeatures3DMRG(vector<string> &features);
	void extractGLCMData3D(vector<T> &glcmData, GLCMFeatures3DMRG<T, R> glcmFeatures);
	void fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int angle, int directionZ);
	boost::multi_array<double, 2> getMatrixSum(boost::multi_array<T, R> inputMatrix, float maxIntensity);
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
	GLCMFeatures3DMRG() {
	}
	~GLCMFeatures3DMRG() {
	}
	void writeCSVFileGLCM3DMRG(GLCMFeatures3DMRG<T, R> glcmFeat, string outputFolder);
	void writeOneFileGLCM3DMRG(GLCMFeatures3DMRG<T, R> glcmFeat, string outputFolder);
	void calculateAllGLCMFeatures3DMRG(GLCMFeatures3DMRG<T, R> &glcmFeat, boost::multi_array<T, R> inputMatrix, float maxIntensity, vector<double> spacing, ConfigFile config);
};





/*!
In the method fill3DMatrices the GLCM matrix is filled with the according values
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: glcm matrix: as reference: matrix that has to be filled
@param[in]: int angle: angle for which the actual glcm matrix is calculated
@param[in]: in directionZ: in which z direction we are going
*/
template <class T, size_t R>
void GLCMFeatures3DMRG<T, R>::fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int angle, int directionZ) {

	float weight;
	int directionX;
	int directionY;
	//define in which direction you have to look for a neighbor
	GLCMFeatures<T, R> glcm;
	glcm.getXYDirections(directionX, directionY, angle);
	//get the vector of the nieghbour pairs
	std::vector<std::pair<T, T> > neighbours = glcm.getNeighbours3D(inputMatrix, angle, directionZ);
	std::pair<T, T> actNeighbour;
	//iterate over the neighbor-vector
	for (int neighbourNr = 0; neighbourNr<neighbours.size(); neighbourNr++) {
		actNeighbour = neighbours[neighbourNr];

		if (!std::isnan(actNeighbour.first) && !std::isnan(actNeighbour.second)) {
			glcMatrix[actNeighbour.first - 1][actNeighbour.second - 1] += 1;
		}
	}
	weight = calculateWeight3D(directionX, directionY, directionZ, normGLCM, actualSpacing);
	multSkalarMatrix(glcMatrix, weight);

}

/*!
In the method getMatrixSum calculates the sum of all calculated GLCM matrices
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: float maxValue: maximum intensity value of VOI,
@param[out]: boost multi_array: summed GLCM matrices
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLCMFeatures3DMRG<T, R>::getMatrixSum(boost::multi_array<T, R> inputMatrix, float maxIntensity) {

	int sizeMatrix = maxIntensity;

	int ang;
	int directionZ = 0;
	glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);
	glcmat sum(boost::extents[sizeMatrix][sizeMatrix]);
	glcmat inverseMatrix(boost::extents[sizeMatrix][sizeMatrix]);

	for (int i = 0; i<5; i++) {
		ang = 180 - i * 45;
		if (ang >0) {
			for (int j = 0; j<3; j++) {
				glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);
				directionZ = -1 + j;
				fill3DMatrices(inputMatrix, GLCMatrix, ang, directionZ);
				inverse(GLCMatrix, inverseMatrix);
				matrixSum(sum, GLCMatrix);
				matrixSum(sum, inverseMatrix);

			}
		}
		else {
			directionZ = 1;
			glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);

			fill3DMatrices(inputMatrix, GLCMatrix, ang, directionZ);
			inverse(GLCMatrix, inverseMatrix);

			matrixSum(sum, GLCMatrix);
			matrixSum(sum, inverseMatrix);
		}
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
void GLCMFeatures3DMRG<T, R>::calculateAllGLCMFeatures3DMRG(GLCMFeatures3DMRG<T, R> &GLCMFeatures3DMRG, boost::multi_array<T, R> inputMatrix, float maxIntensity, vector<double> spacing, ConfigFile config) {
	//get which norm should be used in the calculation of the GLCM features
	normGLCM = config.normGLCM;
	actualSpacing = spacing;

	boost::multi_array<double, 2> GLCM180 = GLCMFeatures3DMRG.getMatrixSum(inputMatrix, maxIntensity);

	GLCMFeatures3DMRG.calculateJointMaximum(GLCM180);
	GLCMFeatures3DMRG.calculateJointAverage(GLCM180);
	GLCMFeatures3DMRG.calculateJointVariance(GLCM180, this->jointAverage);
	GLCMFeatures3DMRG.calculateJointEntropy(GLCM180);
	GLCMFeatures3DMRG.getDiagonalProbabilities(GLCM180);
	GLCMFeatures3DMRG.getCrossProbabilities(GLCM180);
	GLCMFeatures3DMRG.calculateDiffAverage();
	GLCMFeatures3DMRG.calculateDiffVariance(this->diffAverage);
	GLCMFeatures3DMRG.calculateDiffEntropy();
	GLCMFeatures3DMRG.calculateSumAverage();
	GLCMFeatures3DMRG.calculateSumVariance(this->sumAverage);
	GLCMFeatures3DMRG.calculateSumEntropy();
	GLCMFeatures3DMRG.calculateAngSecMoment(GLCM180);
	GLCMFeatures3DMRG.calculateContrast(GLCM180);
	GLCMFeatures3DMRG.calculateDissimilarity(GLCM180);
	GLCMFeatures3DMRG.calculateInverseDiff(GLCM180);
	GLCMFeatures3DMRG.calculateInverseDiffNorm(GLCM180, this->inverseDiff);
	GLCMFeatures3DMRG.calculateInverseDiffMom(GLCM180);
	GLCMFeatures3DMRG.calculateInverseDiffMomNorm(GLCM180);
	GLCMFeatures3DMRG.calculateInverseVariance(GLCM180);
	GLCMFeatures3DMRG.calculateCorrelation(GLCM180);
	GLCMFeatures3DMRG.calculateClusterProminence(GLCM180);
	GLCMFeatures3DMRG.calculateClusterShade(GLCM180);
	GLCMFeatures3DMRG.calculateClusterTendency(GLCM180);
	GLCMFeatures3DMRG.calculateFirstMCorrelation(GLCM180);
	GLCMFeatures3DMRG.calculateSecondMCorrelation(GLCM180);
	GLCMFeatures3DMRG.calculateJointMaximum(GLCM180);
	GLCMFeatures3DMRG.calculateJointAverage(GLCM180);
	GLCMFeatures3DMRG.calculateJointVariance(GLCM180, this->jointAverage);
	GLCMFeatures3DMRG.calculateJointEntropy(GLCM180);
	GLCMFeatures3DMRG.getDiagonalProbabilities(GLCM180);
	GLCMFeatures3DMRG.getCrossProbabilities(GLCM180);
	GLCMFeatures3DMRG.calculateDiffAverage();
	GLCMFeatures3DMRG.calculateDiffVariance(this->diffAverage);
	GLCMFeatures3DMRG.calculateDiffEntropy();
	GLCMFeatures3DMRG.calculateSumAverage();
	GLCMFeatures3DMRG.calculateSumVariance(this->sumAverage);
	GLCMFeatures3DMRG.calculateSumEntropy();
	GLCMFeatures3DMRG.calculateAngSecMoment(GLCM180);
	GLCMFeatures3DMRG.calculateContrast(GLCM180);
	GLCMFeatures3DMRG.calculateDissimilarity(GLCM180);
	GLCMFeatures3DMRG.calculateInverseDiff(GLCM180);
	GLCMFeatures3DMRG.calculateInverseDiffNorm(GLCM180, this->inverseDiff);
	GLCMFeatures3DMRG.calculateInverseDiffMom(GLCM180);
	GLCMFeatures3DMRG.calculateInverseDiffMomNorm(GLCM180);
	GLCMFeatures3DMRG.calculateInverseVariance(GLCM180);
	GLCMFeatures3DMRG.calculateCorrelation(GLCM180);
	GLCMFeatures3DMRG.calculateAutoCorrelation(GLCM180);
	GLCMFeatures3DMRG.calculateClusterProminence(GLCM180);
	GLCMFeatures3DMRG.calculateClusterShade(GLCM180);
	GLCMFeatures3DMRG.calculateClusterTendency(GLCM180);
	GLCMFeatures3DMRG.calculateFirstMCorrelation(GLCM180);
	GLCMFeatures3DMRG.calculateSecondMCorrelation(GLCM180);
}


template <class T, size_t R>
void GLCMFeatures3DMRG<T, R>::writeCSVFileGLCM3DMRG(GLCMFeatures3DMRG<T,R> glcmFeat, string outputFolder)
{
    string csvName = outputFolder + "_glcmFeatures3DWmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glcmCSV;
    glcmCSV.open (name);
    vector<string> features;
    defineGLCMFeatures3DMRG(features);

    vector<T> glcmData;
    extractGLCMData3D(glcmData, glcmFeat);
    for(int i = 0; i< glcmData.size(); i++){
        glcmCSV << "glcmFeatures3DWmrg"<<","<<features[i] <<",";
        glcmCSV << glcmData[i];
        glcmCSV << "\n";
    }
    glcmCSV.close();
}

template <class T, size_t R>
void GLCMFeatures3DMRG<T, R>::writeOneFileGLCM3DMRG(GLCMFeatures3DMRG<T, R> glcmFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glcmCSV;
	glcmCSV.open(name, std::ios_base::app);
	vector<string> features;
	defineGLCMFeatures3DMRG(features);

	vector<T> glcmData;
	extractGLCMData3D(glcmData, glcmFeat);
	for (int i = 0; i< glcmData.size(); i++) {
		glcmCSV << "glcmFeatures3DWmrg" << "," << features[i] << ",";
		glcmCSV << glcmData[i];
		glcmCSV << "\n";
	}
	glcmCSV.close();
}

template <class T, size_t R>
void GLCMFeatures3DMRG<T, R>::defineGLCMFeatures3DMRG(vector<string> &features){
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
void GLCMFeatures3DMRG<T, R>::extractGLCMData3D(vector<T> &glcmData, GLCMFeatures3DMRG<T, R> glcmFeatures){
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


#endif // GLCMFEATURES3DMRG_H_INCLUDED
