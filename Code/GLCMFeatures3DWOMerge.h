#ifndef GLCMFEATURES3DWOMERGE_H_INCLUDED
#define GLCMFEATURES3DWOMERGE_H_INCLUDED

#include "GLCMFeatures3DWMerge.h"

/*! \file */

/*!
The class GLCMFeatures3DWOMerge herites from the class GLCMFeatures. \n
It considers 13 neighbors to calculate the cooccurrence features. \n
It calculates a GLCM matrix for every angle and extracts the feature from every matrix.
Then the mean value over all these features is calculated.\n
All feature calculations are defined in the class GLCMFeatures. \n
This class only contains the calculations of the merged matrix.
*/

template <class T,  size_t R>
class GLCMFeatures3DWOMerge : GLCMFeatures<T, R>{
    private:
        typedef boost::multi_array<float, 2>  glcmat;

        string normGLCM;
        vector<float> actualSpacing;
        void defineGLCMFeatures3DWOMerge(vector<string> &features);
        void extractGLCMData3D(vector<T> &glcmData, GLCMFeatures3DWOMerge<T, R> glcmFeatures);
        void fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<float, 2> &glcMatrix, int angle, int directionZ);
        boost::multi_array<float, 2> getMatrixSum( boost::multi_array<T, R> inputMatrix, float maxIntensity);
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
        void writeCSVFileGLCM3D(GLCMFeatures3DWOMerge<T,R> glcmFeat, string outputFolder);

        void calculateAllGLCMFeatures3DWOMerge(GLCMFeatures3DWOMerge<T,R> &glcmFeat, boost::multi_array<T,R> inputMatrix, float maxIntensity, vector<float> spacing, ConfigFile config);
};





//fill now the GLCMatrices with the right values
template <class T, size_t R>
void GLCMFeatures3DWOMerge<T, R>::fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<float, 2> &glcMatrix, int angle, int directionZ){
    
    float weight;
    int directionX;
    int directionY;
    //define in which direction you have to look for a neighbor
    GLCMFeatures<T,R> glcm;
    glcm.getXYDirections(directionX, directionY, angle);
     //get the vector of the nieghbour pairs
	GLCMFeatures3DWMerge<T, R> glcm3D;
	std::vector<std::pair<T, T> > neighbours = glcm3D.getNeighbours3D(inputMatrix, angle, directionZ);
     std::pair<T, T> actNeighbour;
     //iterate over the neighbor-vector
     for(int neighbourNr =0; neighbourNr<neighbours.size(); neighbourNr++){
        actNeighbour = neighbours[neighbourNr];

        if(!std::isnan(actNeighbour.first) && !std::isnan(actNeighbour.second)){
            glcMatrix[actNeighbour.first -1][actNeighbour.second-1]+=1;
        }
    }
    weight = calculateWeight3D(directionX, directionY, directionZ, normGLCM, actualSpacing);
	multSkalarMatrix(glcMatrix, weight);
	
}


template <class T, size_t R>
boost::multi_array<float, 2> GLCMFeatures3DWOMerge<T,R>::getMatrixSum( boost::multi_array<T, R> inputMatrix, float maxIntensity){

    int sizeMatrix= maxIntensity;

    int ang;
    int directionZ=0;
    glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);
    glcmat sum(boost::extents[sizeMatrix][sizeMatrix]);
    glcmat inverseMatrix(boost::extents[sizeMatrix][sizeMatrix]);

    for(int i = 0; i<5; i++){
        ang = 180-i*45;
        if(ang >0){
            for (int j=0; j<3; j++){
				glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);
                directionZ=-1+j;
                fill3DMatrices(inputMatrix, GLCMatrix, ang, directionZ);
                inverse(GLCMatrix, inverseMatrix);
                matrixSum(sum, GLCMatrix);
                matrixSum(sum, inverseMatrix);

           }
        }
        else{
            directionZ = 1;
			glcmat GLCMatrix(boost::extents[sizeMatrix][sizeMatrix]);

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
	if (sumMatrElement != 0) {
		transform(sum.origin(), sum.origin() + sum.num_elements(),
			sum.origin(), bind2nd(std::divides<double>(), int(sumMatrElement)));
	}
    return sum;
}

template <class T, size_t R>
void GLCMFeatures3DWOMerge<T, R>::calculateAllGLCMFeatures3DWOMerge(GLCMFeatures3DWOMerge<T,R> &GLCMFeatures3DWOMerge, boost::multi_array<T, R> inputMatrix, float maxIntensity, vector<float> spacing, ConfigFile config){
    //get which norm should be used in the calculation of the GLCM features
    normGLCM = config.normGLCM;
    actualSpacing = spacing;

    boost::multi_array<float,2> GLCM180=GLCMFeatures3DWOMerge.getMatrixSum(inputMatrix, maxIntensity);
    GLCMFeatures3DWOMerge.calculateJointMaximum(GLCM180);
    GLCMFeatures3DWOMerge.calculateJointAverage(GLCM180);
    GLCMFeatures3DWOMerge.calculateJointVariance(GLCM180, this->jointAverage);
    GLCMFeatures3DWOMerge.calculateJointEntropy(GLCM180);
    GLCMFeatures3DWOMerge.getDiagonalProbabilities(GLCM180);
    GLCMFeatures3DWOMerge.getCrossProbabilities(GLCM180);
    GLCMFeatures3DWOMerge.calculateDiffAverage();
    GLCMFeatures3DWOMerge.calculateDiffVariance(this->diffAverage);
    GLCMFeatures3DWOMerge.calculateDiffEntropy();
    GLCMFeatures3DWOMerge.calculateSumAverage();
    GLCMFeatures3DWOMerge.calculateSumVariance(this->sumAverage);
    GLCMFeatures3DWOMerge.calculateSumEntropy();
    GLCMFeatures3DWOMerge.calculateAngSecMoment(GLCM180);
    GLCMFeatures3DWOMerge.calculateContrast(GLCM180);
    GLCMFeatures3DWOMerge.calculateDissimilarity(GLCM180);
    GLCMFeatures3DWOMerge.calculateInverseDiff(GLCM180);
    GLCMFeatures3DWOMerge.calculateInverseDiffNorm(GLCM180, this->inverseDiff);
    GLCMFeatures3DWOMerge.calculateInverseDiffMom(GLCM180);
    GLCMFeatures3DWOMerge.calculateInverseDiffMomNorm(GLCM180);
    GLCMFeatures3DWOMerge.calculateInverseVariance(GLCM180);
    GLCMFeatures3DWOMerge.calculateCorrelation(GLCM180);
    GLCMFeatures3DWOMerge.calculateClusterProminence(GLCM180);
    GLCMFeatures3DWOMerge.calculateClusterShade(GLCM180);
    GLCMFeatures3DWOMerge.calculateClusterTendency(GLCM180);
    GLCMFeatures3DWOMerge.calculateFirstMCorrelation(GLCM180);
    GLCMFeatures3DWOMerge.calculateSecondMCorrelation(GLCM180);
    GLCMFeatures3DWOMerge.calculateJointMaximum(GLCM180);
    GLCMFeatures3DWOMerge.calculateJointAverage(GLCM180);
    GLCMFeatures3DWOMerge.calculateJointVariance(GLCM180, this->jointAverage);
    GLCMFeatures3DWOMerge.calculateJointEntropy(GLCM180);
    GLCMFeatures3DWOMerge.getDiagonalProbabilities(GLCM180);
    GLCMFeatures3DWOMerge.getCrossProbabilities(GLCM180);
    GLCMFeatures3DWOMerge.calculateDiffAverage();
    GLCMFeatures3DWOMerge.calculateDiffVariance(this->diffAverage);
    GLCMFeatures3DWOMerge.calculateDiffEntropy();
    GLCMFeatures3DWOMerge.calculateSumAverage();
    GLCMFeatures3DWOMerge.calculateSumVariance(this->sumAverage);
    GLCMFeatures3DWOMerge.calculateSumEntropy();
    GLCMFeatures3DWOMerge.calculateAngSecMoment(GLCM180);
    GLCMFeatures3DWOMerge.calculateContrast(GLCM180);
    GLCMFeatures3DWOMerge.calculateDissimilarity(GLCM180);
    GLCMFeatures3DWOMerge.calculateInverseDiff(GLCM180);
    GLCMFeatures3DWOMerge.calculateInverseDiffNorm(GLCM180, this->inverseDiff);
    GLCMFeatures3DWOMerge.calculateInverseDiffMom(GLCM180);
    GLCMFeatures3DWOMerge.calculateInverseDiffMomNorm(GLCM180);
    GLCMFeatures3DWOMerge.calculateInverseVariance(GLCM180);
    GLCMFeatures3DWOMerge.calculateCorrelation(GLCM180);
    GLCMFeatures3DWOMerge.calculateAutoCorrelation(GLCM180);
    GLCMFeatures3DWOMerge.calculateClusterProminence(GLCM180);
    GLCMFeatures3DWOMerge.calculateClusterShade(GLCM180);
    GLCMFeatures3DWOMerge.calculateClusterTendency(GLCM180);
    GLCMFeatures3DWOMerge.calculateFirstMCorrelation(GLCM180);
    GLCMFeatures3DWOMerge.calculateSecondMCorrelation(GLCM180);
  }

template <class T, size_t R>
void GLCMFeatures3DWOMerge<T, R>::writeCSVFileGLCM3D(GLCMFeatures3DWOMerge<T,R> glcmFeat, string outputFolder)
{
    string csvName = outputFolder + "/glcmFeatures3DWOMerge.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glcmCSV;
    glcmCSV.open (name);
    vector<string> features;
    defineGLCMFeatures3DWOMerge(features);

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
void GLCMFeatures3DWOMerge<T, R>::defineGLCMFeatures3DWOMerge(vector<string> &features){
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
void GLCMFeatures3DWOMerge<T, R>::extractGLCMData3D(vector<T> &glcmData, GLCMFeatures3DWOMerge<T, R> glcmFeatures){
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

#endif // GLCMFEATURES3DWOMERGE_H_INCLUDED
