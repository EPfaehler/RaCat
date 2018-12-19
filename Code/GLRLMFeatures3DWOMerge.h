#ifndef GLRLMFEATURES3DWOMERGE_H_INCLUDED
#define GLRLMFEATURES3DWOMERGE_H_INCLUDED

#include "GLRLMFeatures.h"

/*! \file */

/*!
The class GLRLMFeatures3DWOMerge herites from the class GLRLMFeatures. \n
A GLRLM matrix is calculated for every angle separately. From all these matrices the feature values are then calculated 
and the mean value is calculated.\n
All feature calculations are defined in the class GLRLMFeatures. \n
*/
template <class T,  size_t R=3>
class GLRLMFeatures3DWOMerge : GLRLMFeatures<T,R>{
private:
    GLRLMFeatures<T,R> glrlm;
    double totalSum;
    boost::multi_array<double, 2> createGLRLMatrix3D(boost::multi_array<T, R> inputMatrix, int ang);
    void extractGLRLMData3D(vector<T> &glrlmData, GLRLMFeatures3DWOMerge<T, R> glrlmFeatures);

    void fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glrlMatrix, int directionX, int directionY, int directionZ);
    void getMaxRunLength3D(boost::multi_array<T, R> inputMatrix);
    int maxRunLength;

    int totalNrVoxels;

public:
	GLRLMFeatures3DWOMerge(){}
	~GLRLMFeatures3DWOMerge(){}
    void calculateAllGLRLMFeatures3DWOMerge(GLRLMFeatures3DWOMerge<T,R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElem, ConfigFile config);
    void writeCSVFileGLRLM3DWOMerge(GLRLMFeatures3DWOMerge<T,R> glrlmFeat, string outputFolder);
	void writeOneFileGLRLM3DWOMerge(GLRLMFeatures3DWOMerge<T, R> glrlmFeat, string outputFolder);
    void getXYdirections3D(int &directionX, int &directionY, int &directionZ, int ang);

};

template <class T, size_t R>
void GLRLMFeatures3DWOMerge<T, R>::getMaxRunLength3D(boost::multi_array<T, R> inputMatrix){
    int tempMax = std::max(inputMatrix.shape()[0], inputMatrix.shape()[1]);
    maxRunLength = std::max(tempMax, int(inputMatrix.shape()[2]));
}

//fill now the GLRLMatrices with the right values
template <class T, size_t R>
void GLRLMFeatures3DWOMerge<T, R>::fill3DMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glrlMatrix, int directionX, int directionY, int directionZ){
     double actGreyLevel = 0;
    double actElement = 0;
    int runLength=0;
    int maxRowNr = inputMatrix.shape()[0];
    int maxColNr = inputMatrix.shape()[1];
    int maxDepth = inputMatrix.shape()[2];
    //have a look at the image-matrix slide by slide (2D)
    //look for every grey level separately in every image slide
    for(int actGreyIndex=0; actGreyIndex<this->diffGreyLevels.size(); actGreyIndex++){
        //get the grey level we are interested at the moment
        actGreyLevel = this->diffGreyLevels[actGreyIndex];
        //now look at the image matrix slide by slide
        for(int depth = 0; depth < maxDepth; depth++){
            for(int row = 0; row < maxRowNr;row++){
                for(int column = 0; column < maxColNr; column++){

//                  //at the beginning the run length =0
                    runLength =0;
                    //get the actual matrix element
                    int actDepth;
                    if(directionZ >0){
                        actDepth = depth;
                    }
                    else{
                        actDepth = maxDepth - depth -1;

                    }
                    int actRow;
                    if(directionY < 0){
                        actRow = maxRowNr-row-1;
                    }
                    else{
                        actRow = row;
                    }
                    actElement = inputMatrix[actRow][column][actDepth];

                    //if the actual matrix element is the same as the actual gre level
                    if(actElement == actGreyLevel){
                        //set the run length to 1
                        runLength = 1;

                        //to avoid to take an element more than once, set the element to NAN
                        inputMatrix[actRow][column][actDepth] = NAN;

////                    //now look at the matrix element in the actual direction (depends on the
                        //angle we are interested at the moment
                        int colValue = column + directionX;
                        int rowValue = actRow + directionY;


                        int depthValue = actDepth + directionZ;

                        //now have a look at the following elements in the desired direction
                        //stop as soon as we look at an element diifferent from our actual element
                        while(colValue < maxColNr && rowValue > -1 && rowValue < maxRowNr && colValue > -1 && depthValue < maxDepth && depthValue > -1 && inputMatrix[rowValue][colValue][depthValue] == actGreyLevel ){
                            //for every element we find, count the runLength
                            runLength += 1;
                            inputMatrix[rowValue][colValue][depthValue]=NAN;
                            //go further in the desired direction
                            colValue +=1*directionX;
                            rowValue +=1*directionY;

                            depthValue +=1*directionZ;
                      }
                   }
                   //as soon as we cannot find an element in the desired direction, count one up in the desired
                   //position of the glrl-matrix
                    if(runLength>0&&runLength<maxRunLength+1){
                        glrlMatrix[actGreyIndex][runLength-1]+=1;
                    }
                }
            }
        }
    }
}

template <class T, size_t R>
void GLRLMFeatures3DWOMerge<T, R>::getXYdirections3D(int &directionX, int &directionY, int &directionZ, int ang){
    if(ang == 0){
        directionX = 0;
        directionY = 0;
        directionZ = 1;
    }
    if(ang == 1){
        directionX = 0;
        directionY = 1;
        directionZ = -1;
    }
    if(ang == 2){
        directionX = 0;
        directionY = 1;
        directionZ = 0;
    }
    if(ang == 3){
        directionX = 0;
        directionY = 1;
        directionZ = 1;
    }
    if(ang == 4){
        directionX = 1;
        directionY = -1;
        directionZ = -1;
    }
    if(ang == 5){
        directionX = 1;
        directionY = -1;
        directionZ = 0;
    }
    if(ang == 6){
        directionX = 1;
        directionY = -1;
        directionZ = 1;
    }
    if(ang == 7){
        directionX = 1;
        directionY = 0;
        directionZ = -1;
    }
    if(ang == 8){
        directionX = 1;
        directionY = 0;
        directionZ = 0;
    }
    if(ang == 9){
        directionX = 1;
        directionY = 0;
        directionZ = 1;
    }
    if(ang == 10){
        directionX = 1;
        directionY = 1;
        directionZ = -1;
    }
    if(ang == 11){
        directionX = 1;
        directionY = 1;
        directionZ = 0;
    }
    if(ang == 12){
        directionX = 1;
        directionY = 1;
        directionZ = 1;
    }

}

//create the GLRL-matrices for all the angles
template <class T, size_t R>
boost::multi_array<double, 2> GLRLMFeatures3DWOMerge<T, R>::createGLRLMatrix3D(boost::multi_array<T,R> inputMatrix, int ang){
    typedef boost::multi_array<double, 2> glrlmat;

    int directionX;
    int directionY;
    int directionZ;


//    define the glcMatrices
//
//    fill the 2D Matrices
    getMaxRunLength3D(inputMatrix);
    int sizeMatrix = this->diffGreyLevels.size();

    glrlmat GLRLMatrix(boost::extents[sizeMatrix][maxRunLength]);
    glrlmat sum(boost::extents[sizeMatrix][maxRunLength]);

    getXYdirections3D(directionX, directionY, directionZ, ang);
    fill3DMatrices(inputMatrix, GLRLMatrix, directionX, directionY, directionZ);


    return GLRLMatrix;

}

template <class T, size_t R>
void GLRLMFeatures3DWOMerge<T, R>::calculateAllGLRLMFeatures3DWOMerge(GLRLMFeatures3DWOMerge<T,R> &glrlmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElem, ConfigFile config){
    this->diffGreyLevels = diffGrey;
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

    int directionZ = -1;


    for(int i = 0; i < 13; i++){


        boost::multi_array<double,2> glrlMatrix = createGLRLMatrix3D(inputMatrix, i);

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
        glrlmFeatures.calculateRunPercentage3D(vectorMatrElem, totalSum, 1);
        sumRunPercentage += this->runPercentage;

        glrlmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);
        sumGreyLevelVar += this->greyLevelVar;
        glrlmFeatures.calculateRunLengthVar(probMatrix, meanRun);
        sumRunLengthVar += this->runLengthVar;
        glrlmFeatures.calculateRunEntropy(probMatrix);
        sumRunEntropy += this->runEntropy;
    }

    this->shortRunEmphasis = sumShortRunEmphasis/13;
    this->longRunEmphasis = sumLongRunEmphasis/13;
    this->lowGreyEmph = sumLowGreyEmph/13;
    this->highGreyEmph = sumHighGreyEmph/13;
    this->shortRunLow = sumShortRunLow/13;
    this->shortRunHigh = sumShortRunHigh/13;
    this->longRunLowEmph = sumLongRunLowEmph/13;
    this->longRunHighEmph = sumLongRunHighEmph/13;
    this->greyNonUniformity = sumGreyNonUniformity/13;
    this->greyNonUniformityNorm = sumGreyNonUniformityNorm/13;
    this->runLengthNonUniformity = sumRunLengthNonUniformity/13;
    this->runLengthNonUniformityNorm = sumRunLengthNonUniformityNorm/13;
    this->runPercentage = sumRunPercentage/13;
    this->greyLevelVar = sumGreyLevelVar/13;
    this->runLengthVar = sumRunLengthVar/13;
    this->runEntropy = sumRunEntropy/13;

}

template <class T, size_t R>
void GLRLMFeatures3DWOMerge<T, R>::extractGLRLMData3D(vector<T> &glrlmData, GLRLMFeatures3DWOMerge<T, R> glrlmFeatures){

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

template <class T, size_t R>
void GLRLMFeatures3DWOMerge<T, R>::writeCSVFileGLRLM3DWOMerge(GLRLMFeatures3DWOMerge<T,R> glrlmFeat, string outputFolder)
{
    string csvName = outputFolder + "/GLRLMFeatures3Davg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream glrlmCSV;
    glrlmCSV.open (name);
    vector<string> features;
    glrlmFeat.defineGLRLMFeatures(features);

    vector<T> glrlmData;
    extractGLRLMData3D(glrlmData, glrlmFeat);
    for(int i = 0; i< glrlmData.size(); i++){
        glrlmCSV <<"GLRLMFeatures3Davg"<<","<< features[i] <<",";
        glrlmCSV << glrlmData[i];
        glrlmCSV << "\n";
    }
    glrlmCSV.close();
}


template <class T, size_t R>
void GLRLMFeatures3DWOMerge<T, R>::writeOneFileGLRLM3DWOMerge(GLRLMFeatures3DWOMerge<T, R> glrlmFeat, string outputFolder) {
	string csvName = outputFolder + "/radiomicsFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream glrlmCSV;
	glrlmCSV.open(name, std::ios_base::app);
	vector<string> features;
	glrlmFeat.defineGLRLMFeatures(features);

	vector<T> glrlmData;
	extractGLRLMData3D(glrlmData, glrlmFeat);
	for (int i = 0; i< glrlmData.size(); i++) {
		glrlmCSV << "GLRLMFeatures3Davg" << "," << features[i] << ",";
		glrlmCSV << glrlmData[i];
		glrlmCSV << "\n";
	}
	glrlmCSV.close();

}

#endif // GLRLMFEATURES3DWOMERGE_H_INCLUDED
