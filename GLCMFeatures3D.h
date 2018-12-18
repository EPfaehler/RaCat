#ifndef GLCMFEATURES3D_H_INCLUDED
#define GLCMFEATURES3D_H_INCLUDED

#include "GLCMFeatures.h"
/*! \file */
/*!
In the class GLCMFeatures3D the calculation of the GLC-matrix is modified to 3D mode. \n
So now every voxel has 13 direct neighbors. \n
All feature calculations stay the same as in the 2D approach.
*/
template <class T,  size_t R>
class GLCMFeatures3D : GLCMFeatures<T, R>{
    private:
        typedef boost::multi_array<double, 2>  glcmat;
        void createGLCMMatrix(boost::multi_array<T, R> inputMatrix);
        void defineGLCMFeatures3D(vector<string> &features);
        void extractGLCMData3D(vector<T> &glcmData, GLCMFeatures3D<T, R> glcmFeatures);
        void generate3DMatrices(boost::multi_array<T,R> inputMatrix, glcmat &GLCM3D, int angle, int directionZ);
        std::vector<std::pair<T, T> > getNeighbours3D(boost::multi_array<T, R> inputMatrix, int angle, int directionZ);
        void fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int angle, int directionZ);
        boost::multi_array<double, 2> getMatrixSum( boost::multi_array<T, R> inputMatrix);
        //store different grey levels in vector
        vector<T> diffGreyLevels;

        vector<T> diagonalProbabilities;
        vector<T> crossProbabilities;
        Image<T,R> image{Image<T,R>(4,5,4)};
        vector<T> sumProbRows;
        vector<T> sumProbCols;
        T HX;
        T HXY;
        T HXY1;
        T HXY2;
    public:
        void writeCSVFileGLCM3D(GLCMFeatures3D<T,R> glcmFeat);

        void calculateAllGLCMFeatures3D(GLCMFeatures3D<T,R> &glcmFeat, boost::multi_array<T,R> inputMatrix);
};






template <class T, size_t R>
std::vector<std::pair<T, T> > GLCMFeatures3D<T, R>::getNeighbours3D(boost::multi_array<T, R> inputMatrix, int angle, int directionZ){
//    store the neighbours as pairs in a vector
    std::vector<std::pair<T, T> > neighbours;
    int directionX;
    int directionY;

    //define in which direction you have to look for a neighbor
    GLCMFeatures<T,R> glcm;
    glcm.getXYDirections(directionX, directionY, angle);

    T neighbour1;
    T neighbour2;
    int maxRowNumber = inputMatrix.shape()[0];
    int maxColNumber = inputMatrix.shape()[1];
    int maxDepthNumber = inputMatrix.shape()[2];
    //if we have a 3D image go also in the depth

    for(int depth=0; depth<maxDepthNumber-directionZ*(directionZ+1)/2; depth++){
        for(int row = 0; row<maxRowNumber; row++){
            //find now the neighbors, row by row
            for(int col=0; col<maxColNumber-directionX*(directionX+1)/2; col++){
                neighbour1 = inputMatrix[row][col][depth];
                if(row-directionY>-1&&col+directionX>-1&&depth+directionZ>-1){
                    neighbour2 = inputMatrix[row-directionY][col+directionX][depth+directionZ];
                    neighbours.push_back(std::make_pair(neighbour1, neighbour2));
                }
            }

        }
    }
    return neighbours;
}


//fill now the GLCMatrices with the right values
template <class T, size_t R>
void GLCMFeatures3D<T, R>::fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix, int angle, int directionZ){
     std::vector<std::pair<T, T> > neighbours;
     //get the vector of the nieghbour pairs
     neighbours = getNeighbours3D(inputMatrix, angle, directionZ);
     std::pair<T, T> actNeighbour;
     //iterate over the neighbor-vector
     for(int neighbourNr =0; neighbourNr<neighbours.size(); neighbourNr++){
        actNeighbour = neighbours[neighbourNr];
        if(!std::isnan(actNeighbour.first) && !std::isnan(actNeighbour.second)){
            glcMatrix[actNeighbour.first -1][actNeighbour.second-1]+=1;
     }
}
}


template <class T, size_t R>
boost::multi_array<double, 2> GLCMFeatures3D<T,R>::getMatrixSum( boost::multi_array<T, R> inputMatrix){

    int sizeMatrix= 6;
    int ang;
    int directionZ=0;
    glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);
    glcmat sum(boost::extents[sizeMatrix][sizeMatrix]);
    glcmat inverseMatrix(boost::extents[sizeMatrix][sizeMatrix]);

    for(int i = 0; i<4; i++){
        ang = 180-i*45;
        for (int j=0; j<3; j++){
            directionZ=-1+j;
            fill3DMatrices(inputMatrix, GLCMatrix, ang, directionZ);
            inverse(GLCMatrix, inverseMatrix);
            matrixSum(sum, GLCMatrix);
            matrixSum(sum, inverseMatrix);
       }
    }
    //calculate the sum of all matrix elements (= the number of neighbor-pairs in the matrix
    double sumMatrElement =accumulate(sum.origin(), sum.origin()+sum.num_elements(), 0);
    //divide the whole matrix by the sum to obtain matrix elements representing the probabilities
    //of the occurence of a neighbor pair
    transform( sum.origin(), sum.origin() + sum.num_elements(),
                    sum.origin(),  bind2nd(std::divides<double>(),int(sumMatrElement)));
    return sum;
}

template <class T, size_t R>
void GLCMFeatures3D<T, R>::calculateAllGLCMFeatures3D(GLCMFeatures3D<T,R> &glcmFeatures3D, boost::multi_array<T, R> inputMatrix){


    boost::multi_array<double,2> glcm3D=glcmFeatures3D.getMatrixSum(inputMatrix);
    glcmFeatures3D.calculateJointMaximum(glcm3D);
    glcmFeatures3D.calculateJointAverage(glcm3D);
    glcmFeatures3D.calculateJointVariance(glcm3D, this->jointAverage);
    glcmFeatures3D.calculateJointEntropy(glcm3D);
    glcmFeatures3D.getDiagonalProbabilities(glcm3D);
    glcmFeatures3D.getCrossProbabilities(glcm3D);
    glcmFeatures3D.calculateDiffAverage();
    glcmFeatures3D.calculateDiffVariance(this->diffAverage);
    glcmFeatures3D.calculateDiffEntropy();
    glcmFeatures3D.calculateSumAverage();
    glcmFeatures3D.calculateSumVariance(this->sumAverage);
    glcmFeatures3D.calculateSumEntropy();
    glcmFeatures3D.calculateAngSecMoment(glcm3D);
    glcmFeatures3D.calculateContrast(glcm3D);
    glcmFeatures3D.calculateDissimilarity(glcm3D);
    glcmFeatures3D.calculateInverseDiff(glcm3D);
    glcmFeatures3D.calculateInverseDiffNorm(glcm3D, this->inverseDiff);
    glcmFeatures3D.calculateInverseDiffMom(glcm3D);
    glcmFeatures3D.calculateInverseDiffMomNorm(glcm3D);
    glcmFeatures3D.calculateInverseVariance(glcm3D);
    glcmFeatures3D.calculateCorrelation(glcm3D);
    glcmFeatures3D.calculateClusterProminence(glcm3D);
    glcmFeatures3D.calculateClusterShade(glcm3D);
    glcmFeatures3D.calculateClusterTendency(glcm3D);
    glcmFeatures3D.calculateFirstMCorrelation(glcm3D);
    glcmFeatures3D.calculateSecondMCorrelation(glcm3D);
    glcmFeatures3D.calculateJointMaximum(glcm3D);
    glcmFeatures3D.calculateJointAverage(glcm3D);
    glcmFeatures3D.calculateJointVariance(glcm3D, this->jointAverage);
    glcmFeatures3D.calculateJointEntropy(glcm3D);
    glcmFeatures3D.getDiagonalProbabilities(glcm3D);
    glcmFeatures3D.getCrossProbabilities(glcm3D);
    glcmFeatures3D.calculateDiffAverage();
    glcmFeatures3D.calculateDiffVariance(this->diffAverage);
    glcmFeatures3D.calculateDiffEntropy();
    glcmFeatures3D.calculateSumAverage();
    glcmFeatures3D.calculateSumVariance(this->sumAverage);
    glcmFeatures3D.calculateSumEntropy();
    glcmFeatures3D.calculateAngSecMoment(glcm3D);
    glcmFeatures3D.calculateContrast(glcm3D);
    glcmFeatures3D.calculateDissimilarity(glcm3D);
    glcmFeatures3D.calculateInverseDiff(glcm3D);
    glcmFeatures3D.calculateInverseDiffNorm(glcm3D, this->inverseDiff);
    glcmFeatures3D.calculateInverseDiffMom(glcm3D);
    glcmFeatures3D.calculateInverseDiffMomNorm(glcm3D);
    glcmFeatures3D.calculateInverseVariance(glcm3D);
    glcmFeatures3D.calculateCorrelation(glcm3D);
    glcmFeatures3D.calculateClusterProminence(glcm3D);
    glcmFeatures3D.calculateClusterShade(glcm3D);
    glcmFeatures3D.calculateClusterTendency(glcm3D);
    glcmFeatures3D.calculateFirstMCorrelation(glcm3D);
    glcmFeatures3D.calculateSecondMCorrelation(glcm3D);
  }

template <class T, size_t R>
void GLCMFeatures3D<T, R>::writeCSVFileGLCM3D(GLCMFeatures3D<T,R> glcmFeat)
{
    ofstream glcmCSV;
    glcmCSV.open ("glcmFeatures3D.csv");
    vector<string> features;
    defineGLCMFeatures3D(features);

    vector<T> glcmData;
    extractGLCMData3D(glcmData, glcmFeat);
    for(int i = 0; i< glcmData.size(); i++){
        glcmCSV << features[i] <<";";
        glcmCSV << glcmData[i];
        glcmCSV << "\n";
    }
    glcmCSV.close();
}

template <class T, size_t R>
void GLCMFeatures3D<T, R>::defineGLCMFeatures3D(vector<string> &features){
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
void GLCMFeatures3D<T, R>::extractGLCMData3D(vector<T> &glcmData, GLCMFeatures3D<T, R> glcmFeatures){
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
    glcmData.push_back(glcmFeatures.autoCorrelation);
    glcmData.push_back(glcmFeatures.correlation);
    glcmData.push_back(glcmFeatures.clusterTendency);
    glcmData.push_back(glcmFeatures.clusterShade);
    glcmData.push_back(glcmFeatures.clusterProminence);
    glcmData.push_back(glcmFeatures.firstMCorrelation);
    glcmData.push_back(glcmFeatures.secondMCorrelation);
}

#endif // GLCMFEATURES3D_H_INCLUDED
