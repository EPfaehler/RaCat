#ifndef GLCMFEATURES2DWOMERGE_H_INCLUDED
#define GLCMFEATURES2DWOMERGE_H_INCLUDED
/*! \file */


#include "GLCMFeatures.h"

/*!
The class GLCMFeatures2DWOMerge inherits from the matrix GLCMFeatures. \n
It does not merge the matrices before feature calculation. \n
For every slice a GLCMatrix is calculated and from every of this matrices all features are extracted. \n
Then the average value of all features is calculated.
*/
template <class T,  size_t R=3>
class GLCMFeatures2DWOMerge : GLCMFeatures<T,R>  {
     private:
        GLCMFeatures<T, R> glcmComb;
        int sizeMatrix;

        void extractGLCMDataWOMerge(vector<T> &glcmData, GLCMFeatures2DWOMerge<T, R> glcmFeatures);
        void fill2DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<float, 2> &glcMatrix, int depth, int angle);
        std::vector<std::pair<T, T> > getNeighbours2D(boost::multi_array<T, R> inputMatrix, int depth, int angle);
        boost::multi_array<float, 2> calculateMatrix2DWOMerge( boost::multi_array<T, R> inputMatrix, int depth, int angle);

        vector<T> diagonalProbabilities;
        vector<T> crossProbabilities;
        vector<T> sumProbRows;
        vector<T> sumProbCols;
        T HX;
        T HXY;
        T HXY1;
        T HXY2;

    public:
		GLCMFeatures2DWOMerge() {
		}
		~GLCMFeatures2DWOMerge() {
		}
        void calculateAllGLCMFeatures2DWOMerge(GLCMFeatures2DWOMerge<T,R> &glcmFeat, boost::multi_array<T,R> inputMatrix, float maxIntensity);
        void writeCSVFileGLCM2DWOMerge(GLCMFeatures2DWOMerge<T,R> glcmFeat, string outputFolder);
		void writeOneFileGLCM2DWOMerge(GLCMFeatures2DWOMerge<T, R> glcmFeat, string outputFolder);
};


/*!
In the method calculateMatrix the GLCM-matrices for every direction are calculated, summed up and in the end the sum of this
matrices is divided by the sum of the elements (= nr. of neighbor pairs) to obtain a matrix which contains the probabilities
for the occurence of every neighbor pair.
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: GLCM-matrix
*/
template <class T, size_t R>
boost::multi_array<float, 2> GLCMFeatures2DWOMerge<T,R>::calculateMatrix2DWOMerge( boost::multi_array<T, R> inputMatrix, int depth, int angle){
    typedef boost::multi_array<float, 2> glcmat;

    glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);

    fill2DMatrices(inputMatrix, GLCMatrix, depth,angle);

    return GLCMatrix;

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
void GLCMFeatures2DWOMerge<T, R>::fill2DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<float, 2> &glcMatrix, int depth, int angle){
    //vector in which all neihbor pairs are stored
     std::vector<std::pair<T, T> > neighbours;
     //fill this vector
     neighbours = getNeighbours2D(inputMatrix, depth,  angle);
     std::pair<T, T> actNeighbour;
     //iterate over the neighbor-vector
     for(int neighbourNr =0; neighbourNr<neighbours.size(); neighbourNr++){
        actNeighbour = neighbours[neighbourNr];

        if(!std::isnan(actNeighbour.first) && !std::isnan(actNeighbour.second)){
            glcMatrix[actNeighbour.first-1][actNeighbour.second-1] += 1;
        }
     }

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
std::vector<std::pair<T, T> > GLCMFeatures2DWOMerge<T, R>::getNeighbours2D(boost::multi_array<T, R> inputMatrix, int depth, int angle){
    //store the neighbours as pairs in a vector
    std::vector<std::pair<T, T> > neighbours;
    int directionX;
    int directionY;
    //define in which direction you have to look for a neighbor
    glcmComb.getXYDirections(directionX, directionY, angle);

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
void GLCMFeatures2DWOMerge<T, R>::calculateAllGLCMFeatures2DWOMerge(GLCMFeatures2DWOMerge<T,R> &glcmFeatures, boost::multi_array<T, R> inputMatrix, float maxIntensity){
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
    for(int depth = 0; depth < totalDepth; depth++){
        for(int i = 0; i < 4; i++){
          ang = 180-i*45;
          boost::multi_array<float,2> sum(boost::extents[sizeMatrix][sizeMatrix]) ;
          boost::multi_array<float,2> inverseMatrix(boost::extents[sizeMatrix][sizeMatrix]);
          boost::multi_array<float,2> GLCMatrix = glcmFeatures.calculateMatrix2DWOMerge(inputMatrix, depth, ang);

          sum = GLCMatrix;
          inverse(GLCMatrix, inverseMatrix);
          matrixSum(sum, inverseMatrix);

          //calculate the sum of all matrix elements (= the number of neighbor-pairs in the matrix
          sumMatrElement = accumulate(sum.origin(), sum.origin() + sum.num_elements(), 0);
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
    this->jointMaximum = sumJointMaximum/(4*totalDepth);
    this->jointAverage = sumJointAverage/(totalDepth*4);
    this->jointVariance = sumJointVariance/(totalDepth*4);
    this->jointEntropy = sumJointEntropy/(totalDepth*4);
    this->diffAverage = sumDiffAverage/(totalDepth*4);
    this->diffVariance = sumDiffVariance/(totalDepth*4);
    this->diffEntropy = sumDiffEntropy/(totalDepth*4);
    this->sumAverage = sumSumAverage/(totalDepth*4);
    this->sumVariance = sumSumVariance/(totalDepth*4);
    this->sumEntropy = sumSumEntropy/(totalDepth*4);
    this->angSecMoment = sumAngSecMoment/(totalDepth*4);
    this->contrast = sumContrast/(totalDepth*4);
    this->dissimilarity = sumDissimilarity/(totalDepth*4);
    this->inverseDiff = sumInverseDiff/(totalDepth*4);
    this->inverseDiffNorm = sumInverseDiffNorm/(totalDepth*4);
    this->inverseDiffMom = sumInverseDiffMom/(totalDepth*4);
    this->inverseDiffMomNorm = sumInverseDiffMomNorm/(totalDepth*4);
    this->inverseVar = sumInverseVar/(totalDepth*4);
    this->correlation = sumCorrelation/(totalDepth*4);
    this->autoCorrelation = sumAutoCorrelation/(totalDepth*4);
    this->clusterProminence = sumClusterProminence/(totalDepth*4);
    this->clusterShade = sumClusterShade/(totalDepth*4);
    this->clusterTendency = sumClusterTendency/(totalDepth*4);
    this->firstMCorrelation = sumFirstMCorrelation/(totalDepth*4);
    this->secondMCorrelation = sumSecondMCorrelation/(totalDepth*4);

}

template <class T, size_t R>
void GLCMFeatures2DWOMerge<T, R>::writeCSVFileGLCM2DWOMerge(GLCMFeatures2DWOMerge<T,R> glcmFeat, string outputFolder)
{
    string csvName = outputFolder + "/glcmFeatures2Davg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glcmCSV;
    glcmCSV.open (name);
    vector<string> features;
    glcmComb.defineGLCMFeatures(features);

    vector<T> glcmData;
    extractGLCMDataWOMerge(glcmData, glcmFeat);
    for(int i = 0; i< glcmData.size(); i++){
        glcmCSV << "glcmFeatures2Davg"<<","<<features[i] <<",";
        glcmCSV << glcmData[i];
        glcmCSV << "\n";
    }
    glcmCSV.close();
}

template <class T, size_t R>
void GLCMFeatures2DWOMerge<T, R>::writeOneFileGLCM2DWOMerge(GLCMFeatures2DWOMerge<T, R> glcmFeat, string outputFolder) {
	string csvName = outputFolder + "/radiomicsFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glcmCSV;
	glcmCSV.open(name, std::ios_base::app);
	vector<string> features;
	glcmComb.defineGLCMFeatures(features);

	vector<T> glcmData;
	extractGLCMDataWOMerge(glcmData, glcmFeat);
	for (int i = 0; i< glcmData.size(); i++) {
		glcmCSV << "glcmFeatures2Davg" << "," << features[i] << ",";
		glcmCSV << glcmData[i];
		glcmCSV << "\n";
	}
	glcmCSV.close();

}


template <class T, size_t R>
void GLCMFeatures2DWOMerge<T, R>::extractGLCMDataWOMerge(vector<T> &glcmData, GLCMFeatures2DWOMerge<T, R> glcmFeatures){

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






#endif // GLCMFEATURES2DWOMERGE_H_INCLUDED
