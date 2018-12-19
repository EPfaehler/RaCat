#ifndef GLCMFEATURES2DWMERGE_H_INCLUDED
#define GLCMFEATURES2DWMERGE_H_INCLUDED

#include "GLCMFeatures.h"

/*! \file */

/*!
The class GLCMFeatures2DWMerge herites from the class GLCMFeatures. \n
It merges the matrices of every slice separately and calculates afterwards the features from the merged matrix. \n
Afterwards the mean value of the features for every slice is calculated. \n
The difference between this class and the other GLCMFeature-classes is only the type of merging of the matrix. \n
All feature calculations are defined in the class GLCMFeatures. \n
This class only contains the calculations of the merged matrix.
*/
template <class T,  size_t R=3>
class GLCMFeatures2DWMerge : GLCMFeatures<T,R>  {
     private:

        string normGLCM;
        vector<double> actualSpacing;
        GLCMFeatures<T, R> glcmComb;

        void extractGLCMDataWMerge(vector<T> &glcmData, GLCMFeatures2DWMerge<T, R> glcmFeatures);
        void fill2DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int depth, int angle);
        std::vector<std::pair<T, T> > getNeighbours2D(boost::multi_array<T, R> inputMatrix, int depth, int directionX, int directionY);
        boost::multi_array<double, 2> calculateMatrix2DWMerge( boost::multi_array<T, R> inputMatrix, int depth, float maxIntensity);

        vector<T> diagonalProbabilities;
        vector<T> crossProbabilities;
        vector<T> sumProbRows;
        vector<T> sumProbCols;
        T HX;
        T HXY;
        T HXY1;
        T HXY2;

    public:
		GLCMFeatures2DWMerge() {
		}
		~GLCMFeatures2DWMerge() {
		}
        void calculateAllGLCMFeatures2DWMerge(GLCMFeatures2DWMerge<T,R> &glcmFeat, boost::multi_array<T,R> inputMatrix, float maxIntensity, vector<double> spacing, ConfigFile config);
        void writeCSVFileGLCM2DWMerge(GLCMFeatures2DWMerge<T,R> glcmFeat, string outputFolder);
		void writeOneFileGLCM2DWMerge(GLCMFeatures2DWMerge<T, R> glcmFeat, string outputFolder);
};


/*!
In the method calculateMatrix the GLCM-matrices for every direction are calculated, summed up and in the end the sum of this
matrices is divided by the sum of the elements (= nr. of neighbor pairs) to obtain a matrix which contains the probabilities
for the occurence of every neighbor pair.
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: GLCM-matrix
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLCMFeatures2DWMerge<T,R>::calculateMatrix2DWMerge( boost::multi_array<T, R> inputMatrix, int depth, float maxIntensity){
    typedef boost::multi_array<double, 2> glcmat;
    int sizeMatrix = maxIntensity;
    int ang;
    glcmat sum(boost::extents[sizeMatrix][sizeMatrix]);
    for(int i = 0; i<4; i++){
            ang = 180-i*45;
            glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);
            glcmat inverseMatrix(boost::extents[sizeMatrix][sizeMatrix]);
            fill2DMatrices(inputMatrix, GLCMatrix, depth,ang);
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


/*!
In the method fill2DMatrices the matrix is filled for all directions
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: the GLCMatrix that will be filled with the corresponding values
@param[in]: direction that determines how the matrix will be filled

The function works as follows:\\
it calculates for all matrix elements the corresponding neighbors for the specified angle and stores
these neighbors in a pair: the pair contains as first element the "first neighbor" and as second element the
"second" neighbor; this is done by the method getNeighbors2D.\\
As next step, the GLCMatrix is filled: for every neighborpair the position in the GLCMatrix is determined\\
the value on this position of the GLCMatrix is increased +1
*/
template <class T, size_t R>
void GLCMFeatures2DWMerge<T, R>::fill2DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int depth, int angle){
    //vector in which all neihbor pairs are stored
    std::vector<std::pair<T, T> > neighbours;
    float weight;
    //define directions in order to calculate the weights
    int directionX;
    int directionY;
    //define in which direction you have to look for a neighbor
    glcmComb.getXYDirections(directionX, directionY, angle);
     //fill this vector
     neighbours = getNeighbours2D(inputMatrix, depth,  directionX, directionY);
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

/*!
The method getNeighbors2D stores all neighbor pairs for the desired angle and the actual input matrix in a vector
@param[in] inputMatrix: the original matrix of the VOI
@param[in] angle : the actual angle
@param[out] neighbors: vector containing all neighbor pairs of the actual direction

The method works as follows:\\
It looks for every matrix elements at the neighbors of the desired direction and makes a pair of the actual matrix element and
its neighbors. These pairs are stored in a vector and are returned.
*/

template <class T, size_t R>
std::vector<std::pair<T, T> > GLCMFeatures2DWMerge<T, R>::getNeighbours2D(boost::multi_array<T, R> inputMatrix, int depth, int directionX, int directionY){
    //store the neighbours as pairs in a vector
    std::vector<std::pair<T, T> > neighbours;

    T neighbour1;
    T neighbour2;
    int maxRowNumber = inputMatrix.shape()[0];
    int maxColNumber = inputMatrix.shape()[1];


    for(int row = 0; row < maxRowNumber; row++){
        //find now the neighbors, row by row
        for(int col = 0; col < maxColNumber-directionX*(directionX+1)/2; col++){
            if(row-directionY > -1  && col+directionX > -1){
                neighbour1 = inputMatrix[row][col][depth];
                neighbour2 = inputMatrix[row-directionY][col+directionX][depth];

                neighbours.push_back(std::make_pair(neighbour1, neighbour2));
            }
        }
    }

    return neighbours;
}


template <class T, size_t R>
void GLCMFeatures2DWMerge<T, R>::calculateAllGLCMFeatures2DWMerge(GLCMFeatures2DWMerge<T,R> &glcmFeatures, boost::multi_array<T, R> inputMatrix, float maxIntensity, vector<double> spacing, ConfigFile config){

    //get which norm should be used in the calculation of the GLCM features
    normGLCM = config.normGLCM;
    actualSpacing = spacing;
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


    for(int depth = 0; depth < totalDepth; depth++){
          boost::multi_array<double,2> GLCM180= glcmFeatures.calculateMatrix2DWMerge(inputMatrix, depth, maxIntensity);
          glcmFeatures.calculateJointMaximum(GLCM180);
          sumJointMaximum += this->jointMaximum;
          glcmFeatures.calculateJointAverage(GLCM180);
          sumJointAverage += this->jointAverage;
          glcmFeatures.calculateJointVariance(GLCM180, this->jointAverage);
          sumJointVariance += this->jointVariance;
          glcmFeatures.calculateJointEntropy(GLCM180);
          sumJointEntropy += this->jointEntropy;
          glcmFeatures.getDiagonalProbabilities(GLCM180);
          glcmFeatures.getCrossProbabilities(GLCM180);
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
          glcmFeatures.calculateAngSecMoment(GLCM180);
          sumAngSecMoment += this->angSecMoment;
          glcmFeatures.calculateContrast(GLCM180);
          sumContrast += this->contrast;
          glcmFeatures.calculateDissimilarity(GLCM180);
          sumDissimilarity += this->dissimilarity;
          glcmFeatures.calculateInverseDiff(GLCM180);
          sumInverseDiff += this->inverseDiff;
          glcmFeatures.calculateInverseDiffNorm(GLCM180, this->inverseDiff);
          sumInverseDiffNorm += this->inverseDiffNorm;
          glcmFeatures.calculateInverseDiffMom(GLCM180);
          sumInverseDiffMom += this->inverseDiffMom;
          glcmFeatures.calculateInverseDiffMomNorm(GLCM180);
          sumInverseDiffMomNorm += this->inverseDiffMomNorm;
          glcmFeatures.calculateInverseVariance(GLCM180);
          sumInverseVar += this->inverseVar;
          glcmFeatures.calculateCorrelation(GLCM180);
          sumCorrelation += this->correlation;
          glcmFeatures.calculateAutoCorrelation(GLCM180);
          sumAutoCorrelation += this->autoCorrelation;
          glcmFeatures.calculateClusterProminence(GLCM180);
          sumClusterProminence += this->clusterProminence;
          glcmFeatures.calculateClusterShade(GLCM180);
          sumClusterShade += this->clusterShade;
          glcmFeatures.calculateClusterTendency(GLCM180);
          sumClusterTendency += this->clusterTendency;
          glcmFeatures.calculateFirstMCorrelation(GLCM180);
          sumFirstMCorrelation += this->firstMCorrelation;
          glcmFeatures.calculateSecondMCorrelation(GLCM180);
          sumSecondMCorrelation += this->secondMCorrelation;
    }
    this->jointMaximum = sumJointMaximum/totalDepth;
    this->jointAverage = sumJointAverage/totalDepth;
    this->jointVariance = sumJointVariance/totalDepth;
    this->jointEntropy = sumJointEntropy/totalDepth;
    this->diffAverage = sumDiffAverage/totalDepth;
    this->diffVariance = sumDiffVariance/totalDepth;
    this->diffEntropy = sumDiffEntropy/totalDepth;
    this->sumAverage = sumSumAverage/totalDepth;
    this->sumVariance = sumSumVariance/totalDepth;
    this->sumEntropy = sumSumEntropy/totalDepth;
    this->angSecMoment = sumAngSecMoment/totalDepth;
    this->contrast = sumContrast/totalDepth;
    this->dissimilarity = sumDissimilarity/totalDepth;
    this->inverseDiff = sumInverseDiff/totalDepth;
    this->inverseDiffNorm = sumInverseDiffNorm/totalDepth;
    this->inverseDiffMom = sumInverseDiffMom/totalDepth;
    this->inverseDiffMomNorm = sumInverseDiffMomNorm/totalDepth;
    this->inverseVar = sumInverseVar/totalDepth;
    this->correlation = sumCorrelation/totalDepth;
    this->autoCorrelation = sumAutoCorrelation/totalDepth;
    this->clusterProminence = sumClusterProminence/totalDepth;
    this->clusterShade = sumClusterShade/totalDepth;
    this->clusterTendency = sumClusterTendency/totalDepth;
    this->firstMCorrelation = sumFirstMCorrelation/totalDepth;
    this->secondMCorrelation = sumSecondMCorrelation/totalDepth;

}

template <class T, size_t R>
void GLCMFeatures2DWMerge<T, R>::writeCSVFileGLCM2DWMerge(GLCMFeatures2DWMerge<T,R> glcmFeat, string outputFolder)
{
    string csvName = outputFolder + "/glcmFeatures2Dmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glcmCSV;
    glcmCSV.open (name);
    vector<string> features;
    glcmComb.defineGLCMFeatures(features);

    vector<T> glcmData;
    extractGLCMDataWMerge(glcmData, glcmFeat);
    for(int i = 0; i< glcmData.size(); i++){
        glcmCSV << "glcmFeatures2Dmrg" << "," << features[i] <<",";
        glcmCSV << glcmData[i];
        glcmCSV << "\n";
    }
    glcmCSV.close();
}

template <class T, size_t R>
void GLCMFeatures2DWMerge<T, R>::writeOneFileGLCM2DWMerge(GLCMFeatures2DWMerge<T, R> glcmFeat, string outputFolder) {
	string csvName = outputFolder + "/radiomicsFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glcmCSV;
	glcmCSV.open(name, std::ios_base::app);
	vector<string> features;
	glcmComb.defineGLCMFeatures(features);

	vector<T> glcmData;
	extractGLCMDataWMerge(glcmData, glcmFeat);
	for (int i = 0; i< glcmData.size(); i++) {
		glcmCSV << "glcmFeatures2Dmrg"<<","<< features[i] << ",";
		glcmCSV << glcmData[i];
		glcmCSV << "\n";
	}
	glcmCSV.close();
}


template <class T, size_t R>
void GLCMFeatures2DWMerge<T, R>::extractGLCMDataWMerge(vector<T> &glcmData, GLCMFeatures2DWMerge<T, R> glcmFeatures){

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



#endif // GLCMFEATURES2DWMERGE_H_INCLUDED
