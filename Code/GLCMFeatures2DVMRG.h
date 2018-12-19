#ifndef GLCMFEATURES2DVMRG_H_INCLUDED
#define GLCMFEATURES2DVMRG_H_INCLUDED

#include "GLCMFeatures.h"

/*! \file */

/*!
The class GLCMFeatures2DVMRG herites from the class GLCMFeatures. \n
It merges the matrices of all slices and calculates afterwards the features from the merged matrix. \n
All feature calculations are defined in the class GLCMFeatures. \n
This class only contains the calculations of the merged matrix.
*/
template <class T,  size_t R=3>
class GLCMFeatures2DVMRG : GLCMFeatures<T,R>  {
     private:
        GLCMFeatures<T, R> glcmComb;
        string normGLCM;
        vector<double> actualSpacing;
        void extractGLCMDataVMRG(vector<T> &glcmData, GLCMFeatures2DVMRG<T, R> glcmFeatures);
        void fill2DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int depth, int angle);
        boost::multi_array<double, 2> calculateMatrix2DVMRG( boost::multi_array<T, R> inputMatrix, float maxIntensity);

        int N_g;
        vector<T> diagonalProbabilities;
        vector<T> crossProbabilities;
        vector<T> sumProbRows;
        vector<T> sumProbCols;
        T HX;
        T HXY;
        T HXY1;
        T HXY2;

    public:
		GLCMFeatures2DVMRG() {
		}
		~GLCMFeatures2DVMRG() {
		}
        void calculateAllGLCMFeatures2DVMRG(GLCMFeatures2DVMRG<T,R> &glcmFeat, boost::multi_array<T,R> inputMatrix, float maxIntensity, vector<double> spacing, ConfigFile config);
        void writeCSVFileGLCM2DVMRG(GLCMFeatures2DVMRG<T,R> glcmFeat, string outputFolder);
		void writeOneFileGLCM2DVMRG(GLCMFeatures2DVMRG<T, R> glcmFeat, string outputFolder);
};


/*!
In the method calculateMatrix the GLCM-matrices for every direction are calculated, summed up and in the end the sum of this
matrices is divided by the sum of the elements (= nr. of neighbor pairs) to obtain a matrix which contains the probabilities
for the occurence of every neighbor pair.
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: GLCM-matrix
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLCMFeatures2DVMRG<T,R>::calculateMatrix2DVMRG( boost::multi_array<T, R> inputMatrix, float maxIntensity){
    typedef boost::multi_array<double, 2> glcmat;
    int sizeMatrix = maxIntensity;
    int ang;
	

    glcmat sum(boost::extents[sizeMatrix][sizeMatrix]);
    for(int depth = 0; depth < inputMatrix.shape()[2]; depth++){
        for(int i = 0; i<4; i++){
            glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);
            glcmat inverseMatrix(boost::extents[sizeMatrix][sizeMatrix]);
            ang = 180-i*45;
            fill2DMatrices(inputMatrix, GLCMatrix, depth,ang);
            inverse(GLCMatrix, inverseMatrix);
            matrixSum(sum, GLCMatrix);
            matrixSum(sum, inverseMatrix);
        }

    }
    //calculate the sum of all matrix elements (= the number of neighbor-pairs in the matrix
    double sumMatrElement = accumulate(sum.origin(), sum.origin() + sum.num_elements(), 0);
    //divide the whole matrix by the sum to obtain matrix elements representing the probabilities
    //of the occurence of a neighbor pair
    transform( sum.origin(), sum.origin() + sum.num_elements(),
                    sum.origin(),  bind2nd(std::divides<double>(),int(sumMatrElement)));
    return sum;

}


/*!
In the method fill2DMatrices the matrix is filled for all directions
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: the GLCMatrix that will be filled with the corresponding values
@param[in]: direction that determines how the matrix will be filled

The function works as follows: \n
it calculates for all matrix elements the corresponding neighbors for the specified angle and stores
these neighbors in a pair: the pair contains as first element the "first neighbor" and as second element the
"second" neighbor; this is done by the method getNeighbors2D. \n
As next step, the GLCMatrix is filled: for every neighborpair the position in the GLCMatrix is determined \n
the value on this position of the GLCMatrix is increased +1
*/
template <class T, size_t R>
void GLCMFeatures2DVMRG<T, R>::fill2DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int depth, int angle){
    //vector in which all neihbor pairs are stored
    std::vector<std::pair<T, T> > neighbours;
    float weight;
    //define directions in order to calculate the weights
    int directionX;
    int directionY;
    //define in which direction you have to look for a neighbor
    glcmComb.getXYDirections(directionX, directionY, angle);
     //fill this vector
     neighbours = glcmComb.getNeighbours2D(inputMatrix, depth,  directionX, directionY);
     std::pair<T, T> actNeighbour;
     //iterate over the neighbor-vector
     for(int neighbourNr =0; neighbourNr<neighbours.size(); neighbourNr++){
        actNeighbour = neighbours[neighbourNr];
        if(!std::isnan(actNeighbour.first) && !std::isnan(actNeighbour.second)){
            glcMatrix[actNeighbour.first-1][actNeighbour.second-1] += 1;
        }
     }
	 weight = calculateWeight2D(directionX, directionY, normGLCM, actualSpacing);
     multSkalarMatrix(glcMatrix, weight);
}




template <class T, size_t R>
void GLCMFeatures2DVMRG<T, R>::calculateAllGLCMFeatures2DVMRG(GLCMFeatures2DVMRG<T,R> &glcmFeatures, boost::multi_array<T, R> inputMatrix, float maxIntensity, vector<double> spacing, ConfigFile config){
    //get which norm should be used in the calculation of the GLCM features
    normGLCM = config.normGLCM;
    actualSpacing = spacing;
    boost::multi_array<double,2> glcmVMRG= glcmFeatures.calculateMatrix2DVMRG(inputMatrix, maxIntensity);
	glcmFeatures.calculateJointMaximum(glcmVMRG);
    glcmFeatures.calculateJointAverage(glcmVMRG);
    glcmFeatures.calculateJointVariance(glcmVMRG, this->jointAverage);
    glcmFeatures.calculateJointEntropy(glcmVMRG);
    glcmFeatures.getDiagonalProbabilities(glcmVMRG);
    glcmFeatures.getCrossProbabilities(glcmVMRG);
    glcmFeatures.calculateDiffAverage();
    glcmFeatures.calculateDiffVariance(this->diffAverage);
    glcmFeatures.calculateDiffEntropy();
    glcmFeatures.calculateSumAverage();
    glcmFeatures.calculateSumVariance(this->sumAverage);
    glcmFeatures.calculateSumEntropy();
    glcmFeatures.calculateAngSecMoment(glcmVMRG);
    glcmFeatures.calculateContrast(glcmVMRG);
    glcmFeatures.calculateDissimilarity(glcmVMRG);
    glcmFeatures.calculateInverseDiff(glcmVMRG);
    glcmFeatures.calculateInverseDiffNorm(glcmVMRG, this->inverseDiff);
    glcmFeatures.calculateInverseDiffMom(glcmVMRG);
    glcmFeatures.calculateInverseDiffMomNorm(glcmVMRG);
    glcmFeatures.calculateInverseVariance(glcmVMRG);
    glcmFeatures.calculateCorrelation(glcmVMRG);
    glcmFeatures.calculateAutoCorrelation(glcmVMRG);
    glcmFeatures.calculateClusterProminence(glcmVMRG);
    glcmFeatures.calculateClusterShade(glcmVMRG);
    glcmFeatures.calculateClusterTendency(glcmVMRG);
    glcmFeatures.calculateFirstMCorrelation(glcmVMRG);
    glcmFeatures.calculateSecondMCorrelation(glcmVMRG);
}

template <class T, size_t R>
void GLCMFeatures2DVMRG<T, R>::writeCSVFileGLCM2DVMRG(GLCMFeatures2DVMRG<T,R> glcmFeat, string outputFolder)
{
    string csvName = outputFolder + "_glcmFeatures2Dvmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glcmCSV;
    glcmCSV.open (name);
    vector<string> features;
    glcmComb.defineGLCMFeatures(features);

    vector<T> glcmData;
    extractGLCMDataVMRG(glcmData, glcmFeat);
    for(int i = 0; i< glcmData.size(); i++){
        glcmCSV << "glcmFeatures2Dvmrg" << "," << features[i] <<",";
        glcmCSV << glcmData[i];
        glcmCSV << "\n";
    }
    glcmCSV.close();
}

template <class T, size_t R>
void GLCMFeatures2DVMRG<T, R>::writeOneFileGLCM2DVMRG(GLCMFeatures2DVMRG<T, R> glcmFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glcmCSV;
	glcmCSV.open(name, std::ios_base::app);
	vector<string> features;
	glcmComb.defineGLCMFeatures(features);

	vector<T> glcmData;
	extractGLCMDataVMRG(glcmData, glcmFeat);
	for (int i = 0; i< glcmData.size(); i++) {
		glcmCSV << "glcmFeatures2Dvmrg"<<","<<features[i] << ",";
		glcmCSV << glcmData[i];
		glcmCSV << "\n";
	}
	glcmCSV.close();
}


template <class T, size_t R>
void GLCMFeatures2DVMRG<T, R>::extractGLCMDataVMRG(vector<T> &glcmData, GLCMFeatures2DVMRG<T, R> glcmFeatures){

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





#endif // GLCMFEATURES2DVMRG_H_INCLUDED
