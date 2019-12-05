#ifndef GLSZMFEATURES2DAVG_H_INCLUDED
#define GLSZMFEATURES2DAVG_H_INCLUDED

#include "GLSZMFeatures2D.h"

/*! \file */

/*!
The class GLSZMFeatures2DAVG herites from the class GLSZMFeatures2D, because the feature calculations are the same.
Only the matrix calculation is different\n
All feature calculations are defined in the class GLRLMFeatures. \n
This class calculates a GLSZM matrix for every slice of the VOI and extracts the feature values from every slice. Then a mean value of all
these feature values is calculated.\n
\n
For grey level size zone matrices, groups of connected voxels with a specific grey value and size are grouped.
A voxel is connected with another voxel if they have the same grey level.
*/

template <class T,  size_t R>
class  GLSZMFeatures2DAVG : public GLSZMFeatures2DMRG<T, R>{
    private:
        void extractGLSZMData(vector<T> &GLSZMData, GLSZMFeatures2DAVG<T, R> GLSZMFeatures);
		GLSZMFeatures2DMRG<T, R> glszm2D;
		GLRLMFeatures<T, R> glrlm;
        void fill2DGLSZMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<float, 2> &glcMatrix, int depth);
        boost::multi_array<float, 2> getGLSZMMatrix( boost::multi_array<T, R> inputMatrix, int depth);
        int getBiggestZoneNr(boost::multi_array<T, R> inputMatrix);
        int maxZoneSize;
        vector<float> rowSums;
        vector<float> colSums;
    public:
        void getNeighbors(boost::multi_array<T, R> &inputMatrix, T actElement, vector<vector< int> > &matrixIndices);
        void getALLXYDirections(int &directionX, int &directionY, int angle);

        void calculateAllGLSZMFeatures2DAVG(GLSZMFeatures2DAVG<T,R> &GLSZMFeat, boost::multi_array<T,R> inputMatrix, vector<T> diffGrey, ConfigFile config);
        void writeCSVFileGLSZM2DAVG(GLSZMFeatures2DAVG<T,R> GLSZMFeat, string outputFolder);
		void writeOneFileGLSZM2DAVG(GLSZMFeatures2DAVG<T, R> GLSZMFeat, ConfigFile config, int &parameterSpaceNr);

};


/*!
In the method fill2DGLSZMatrices the GLSZM matrices are filled using the matrix filled with the intensity values of the VOI.
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: GLSZM is given as reference

It checks the neighborhood of every element of the input matrix. Is an element already considered as an element
of a neighborhood, it is set to NAN.\n
The size of the neighborhoods are stored in the matrix.
*/
template <class T, size_t R>
void GLSZMFeatures2DAVG<T, R>::fill2DGLSZMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<float, 2> &GLSZMatrix, int depth){
     T actualElement;
     vector<vector<int> > matrixIndices;
     vector<int> actualIndex;
	 int actualGreyIndex;
     //do this grey level by grey level
	 maxZoneSize = 0;
	 int tempZoneSize = 0;
     //go element by element through the image
     //look at the neighbors of every element and check if they have the same grey value
     //set every element already seen to NAN
	 for (int row = 0; row < inputMatrix.shape()[0]; row++) {
		 for (int col = 0; col < inputMatrix.shape()[1]; col++) {
			 actualElement = inputMatrix[row][col][depth];
			 if (!isnan(actualElement)) {
				 actualGreyIndex = glrlm.findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actualElement);
				 inputMatrix[row][col][depth] = NAN;
				 actualIndex.push_back(row);
				 actualIndex.push_back(col);
				 actualIndex.push_back(depth);
				 matrixIndices.push_back(actualIndex);
				 actualIndex.clear();
				 getNeighbors(inputMatrix, actualElement, matrixIndices);
			 }
			 tempZoneSize = matrixIndices.size();
			 if (tempZoneSize > maxZoneSize) {
				 maxZoneSize = tempZoneSize;
			 }
			 if (matrixIndices.size() > 0) {
				 GLSZMatrix[actualGreyIndex][matrixIndices.size() - 1] += 1;
				 matrixIndices.clear();
			 }
		 }
	 }
    
}



template <class T, size_t R>
int GLSZMFeatures2DAVG<T, R>::getBiggestZoneNr(boost::multi_array<T, R> inputMatrix){
     T actualElement;
     vector<vector<int> > matrixIndices;
     vector<int> actualIndex;
     maxZoneSize=0;
     int tempZoneSize=0;
	 int actualGreyIndex;
	 for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {

		 for (int row = 0; row < inputMatrix.shape()[0]; row++) {
			 for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				 actualElement = inputMatrix[row][col][depth];
				 if (!isnan(actualElement)) {
					 actualGreyIndex = glszm2D.findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actualElement);
					 inputMatrix[row][col][depth] = NAN;
					 actualIndex.push_back(row);
					 actualIndex.push_back(col);
					 actualIndex.push_back(depth);
					 matrixIndices.push_back(actualIndex);
					 actualIndex.clear();
					 getNeighbors(inputMatrix, actualElement, matrixIndices);

				 }
			 }
			 tempZoneSize = matrixIndices.size();
			 if (tempZoneSize > maxZoneSize) {
				 maxZoneSize = tempZoneSize;
			 }
			 matrixIndices.clear();
		 }

	 }
    return maxZoneSize;
}

/*!
In the method getNeighbors the number of voxels in the biggest zone is determined
@param[in] inputMatrix: the original matrix of the VOI
@param[out]: int biggestZoneNr: the nr of voxels in the biggest zone

The function works as follows: \n
It gets the zone sizes for every intensity value and compares them \n
The biggest zone size is stored in the int value maxZoneNumber
*/
template <class T, size_t R>
void GLSZMFeatures2DAVG<T, R>::getNeighbors(boost::multi_array<T, R> &inputMatrix, T actElement, vector<vector< int> > &matrixIndices){
    vector<int> actIndex;
    vector<int> newIndex;
    vector<int> tempIndex;
    int neighborDist =1;
    int directionX;
    int directionY;
    vector<vector< int> > tempMatrixIndices = matrixIndices;
    int ang;

    //look at all indices which are saved in the vector matrix indices
    while(tempMatrixIndices.size()>0){
        //get the last element of the vector and delete it
        actIndex = tempMatrixIndices.back();
        tempMatrixIndices.pop_back();
        tempIndex = actIndex;
        for(int i=0; i<8; i++){
            ang = 360-i*45;
            getALLXYDirections(directionX, directionY, ang);
            actIndex = tempIndex;
            if(actIndex[0]+directionX>-1 && actIndex[0]+directionX<inputMatrix.shape()[0] && actIndex[1]+directionY>-1 &&actIndex[1]+directionY<inputMatrix.shape()[1] && actElement == inputMatrix[actIndex[0]+directionX][actIndex[1]+directionY][actIndex[2]]){
                    inputMatrix[actIndex[0]+directionX][actIndex[1]+directionY][actIndex[2]] = NAN;
                    actIndex[0] += directionX;
                    actIndex[1] += directionY;
                    tempMatrixIndices.push_back(actIndex);
                    matrixIndices.push_back(actIndex);
            }
        }
    }
}


//

template <class T, size_t R>
void GLSZMFeatures2DAVG<T, R>::getALLXYDirections(int &directionX, int &directionY, int angle){
    //if we only go in the depth
    if(angle==360){
        directionX=0;
        directionY=1;
    }
        //if angle is 180 degrees, only look in x direction
    else if(angle==315){
        directionX=1;
        directionY=1;
    }
    //if angle is 90 degrees only look in y direction
    else if(angle==270){
        directionX=1;
        directionY=0;
    }
    //if angle is in 45 degrees direction look in x and y direction
    else if(angle==225){
        directionX=1;
        directionY=-1;
    }
    else if(angle==180){
        directionX=0;
        directionY=-1;
    }
    else if(angle==135){
        directionX=-1;
        directionY=-1;
    }
    else if(angle==90){
        directionX=-1;
        directionY=0;
    }
    else if(angle==45){
        directionX=-1;
        directionY=1;
    }
    else if(angle==0){
        directionX=0;
        directionY=0;
    }

    else{
        std::cout<<" incorrect angle"<<std::endl;
    }
}

/*!
In the method getGLSZMMatrix the GLSZM matrices with the right size are generated and filled using the fill2DGLSZM function.
@param[in] inputMatrix: the original matrix of the VOI
@param[in] depth: number of actual slice for which the GLSZM should be calculated
@param[out]: filled GLSZM

This function only initiates a GLSZM with the right size.
*/
template <class T, size_t R>
boost::multi_array<float, 2> GLSZMFeatures2DAVG<T,R>::getGLSZMMatrix( boost::multi_array<T, R> inputMatrix, int depth){
    typedef boost::multi_array<float, 2>  GLSZMat;
    //maxZoneSize=getBiggestZoneNr(inputMatrix);
    int sizeMatrix= (this->diffGreyLevels).size();
	maxZoneSize = inputMatrix.shape()[0] * inputMatrix.shape()[1];
    GLSZMat GLSZMatrix(boost::extents[sizeMatrix][maxZoneSize]);
    fill2DGLSZMatrices(inputMatrix, GLSZMatrix, depth);

	GLSZMat smallGLSZM(boost::extents[sizeMatrix][maxZoneSize]);
	for (int row = 0; row < sizeMatrix; row++) {
		for (int col = 0; col < maxZoneSize; col++) {
			smallGLSZM[row][col] = GLSZMatrix[row][col];
		}
	}
    return smallGLSZM;

}

template <class T, size_t R>
void GLSZMFeatures2DAVG<T, R>::calculateAllGLSZMFeatures2DAVG(GLSZMFeatures2DAVG<T,R> &GLSZMFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, ConfigFile config){
	GLSZMFeatures.getConfigValues(config);
	this->diffGreyLevels = diffGrey;
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

    int totalDepth = inputMatrix.shape()[2];


    vector<float> rowSums;
    vector<float> colSums;

    float meanGrey;
    float meanRun;

    for(int depth = 0; depth < totalDepth; depth++){
        boost::multi_array<float,2> GLSZM=GLSZMFeatures.getGLSZMMatrix(inputMatrix, depth);
		
		
        float totalSum = GLSZMFeatures.calculateTotalSum(GLSZM);

        totalSum = GLSZMFeatures.calculateTotalSum(GLSZM);
        rowSums = GLSZMFeatures.calculateRowSums(GLSZM);
        colSums = GLSZMFeatures.calculateColSums(GLSZM);

        boost::multi_array<float,2> probMatrix = GLSZMFeatures.calculateProbMatrix(GLSZM, totalSum);
        meanGrey = GLSZMFeatures.calculateMeanProbGrey(probMatrix);
        meanRun = GLSZMFeatures.calculateMeanProbRun(probMatrix);


        GLSZMFeatures.calculateShortRunEmphasis(rowSums, totalSum);
        sumShortRunEmphasis += this->shortRunEmphasis;

        GLSZMFeatures.calculateLongRunEmphasis(rowSums, totalSum);
        sumLongRunEmphasis += this->longRunEmphasis;
        GLSZMFeatures.calculateLowGreyEmph(colSums, totalSum);
        sumLowGreyEmph += this->lowGreyEmph;
        GLSZMFeatures.calculateHighGreyEmph(colSums, totalSum);
        sumHighGreyEmph += this->highGreyEmph;
        GLSZMFeatures.calculateShortRunLow(GLSZM, totalSum);
        sumShortRunLow += this->shortRunLow;
        GLSZMFeatures.calculateShortRunHigh(GLSZM, totalSum);
        sumShortRunHigh += this->shortRunHigh;
        GLSZMFeatures.calculateLongRunLowEmph(GLSZM, totalSum);
        sumLongRunLowEmph += this->longRunLowEmph;
        GLSZMFeatures.calculateLongRunHighEmph(GLSZM, totalSum);
        sumLongRunHighEmph += this->longRunHighEmph;
        GLSZMFeatures.calculateGreyNonUniformity(colSums, totalSum);
        sumGreyNonUniformity += this->greyNonUniformity;
        GLSZMFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
        sumGreyNonUniformityNorm += this->greyNonUniformityNorm;
        GLSZMFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
        sumRunLengthNonUniformity += this->runLengthNonUniformity;
        GLSZMFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
        sumRunLengthNonUniformityNorm += this->runLengthNonUniformityNorm;
        GLSZMFeatures.calculateRunPercentage(inputMatrix, depth, totalSum, 1);
        sumRunPercentage += this->runPercentage;

        GLSZMFeatures.calculateGreyLevelVar(probMatrix, meanGrey);
        sumGreyLevelVar += this->greyLevelVar;
        GLSZMFeatures.calculateRunLengthVar(probMatrix, meanRun);
        sumRunLengthVar += this->runLengthVar;
        GLSZMFeatures.calculateRunEntropy(probMatrix);
        sumRunEntropy += this->runEntropy;

        }

        this->shortRunEmphasis = sumShortRunEmphasis/(totalDepth);
        this->longRunEmphasis = sumLongRunEmphasis/(totalDepth);
        this->lowGreyEmph = sumLowGreyEmph/(totalDepth);
        this->highGreyEmph = sumHighGreyEmph/(totalDepth);
        this->shortRunLow = sumShortRunLow/(totalDepth);
        this->shortRunHigh = sumShortRunHigh/(totalDepth);
        this->longRunLowEmph = sumLongRunLowEmph/(totalDepth);
        this->longRunHighEmph = sumLongRunHighEmph/(totalDepth);
        this->greyNonUniformity = sumGreyNonUniformity/(totalDepth);
        this->greyNonUniformityNorm = sumGreyNonUniformityNorm/(totalDepth);
        this->runLengthNonUniformity = sumRunLengthNonUniformity/(totalDepth);
        this->runLengthNonUniformityNorm = sumRunLengthNonUniformityNorm/(totalDepth);
        this->runPercentage = sumRunPercentage/(totalDepth);

        this->greyLevelVar = sumGreyLevelVar/(totalDepth);
        this->runLengthVar = sumRunLengthVar/(totalDepth);
        this->runEntropy = sumRunEntropy/(totalDepth);

  }

  template <class T, size_t R>
void GLSZMFeatures2DAVG<T, R>::extractGLSZMData(vector<T> &GLSZMData, GLSZMFeatures2DAVG<T, R> GLSZMFeatures){

    GLSZMData.push_back(GLSZMFeatures.shortRunEmphasis);
    GLSZMData.push_back(GLSZMFeatures.longRunEmphasis);
    GLSZMData.push_back(GLSZMFeatures.lowGreyEmph);
    GLSZMData.push_back(GLSZMFeatures.highGreyEmph);
    GLSZMData.push_back(GLSZMFeatures.shortRunLow);
    GLSZMData.push_back(GLSZMFeatures.shortRunHigh);
    GLSZMData.push_back(GLSZMFeatures.longRunLowEmph);
    GLSZMData.push_back(GLSZMFeatures.longRunHighEmph);
    GLSZMData.push_back(GLSZMFeatures.greyNonUniformity);
    GLSZMData.push_back(GLSZMFeatures.greyNonUniformityNorm);
    GLSZMData.push_back(GLSZMFeatures.runLengthNonUniformity);   //check if I have to give the matrix
    GLSZMData.push_back(GLSZMFeatures.runLengthNonUniformityNorm);
    GLSZMData.push_back(GLSZMFeatures.runPercentage);
    GLSZMData.push_back(GLSZMFeatures.greyLevelVar);
    GLSZMData.push_back(GLSZMFeatures.runLengthVar);
    GLSZMData.push_back(GLSZMFeatures.runEntropy);

}

template <class T, size_t R>
void GLSZMFeatures2DAVG<T, R>::writeCSVFileGLSZM2DAVG(GLSZMFeatures2DAVG<T,R> GLSZMFeat, string outputFolder)
{
    string csvName = outputFolder + "_GLSZMFeatures2Davg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream GLSZMCSV;
    GLSZMCSV.open (name);
    vector<string> features;
    glszm2D.defineGLSZMFeatures(features);

    vector<T> GLSZMData;
    extractGLSZMData(GLSZMData, GLSZMFeat);
    for(int i = 0; i< GLSZMData.size(); i++){
        GLSZMCSV <<"GLSZMFeatures2Davg"<<","<< features[i] <<",";
        GLSZMCSV << GLSZMData[i];
        GLSZMCSV << "\n";
    }
    GLSZMCSV.close();
}

template <class T, size_t R>
void GLSZMFeatures2DAVG<T, R>::writeOneFileGLSZM2DAVG(GLSZMFeatures2DAVG<T, R> GLSZMFeat, ConfigFile config, int &parameterSpaceNr) {
	string csvName;
	if (config.csvOutput == 1) {
		csvName = config.outputFolder + ".csv";
	}
	else if (config.ontologyOutput == 1) {
		csvName = config.outputFolder + "/feature_table.csv";
	}
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream GLSZMCSV;
	GLSZMCSV.open(name, std::ios_base::app);
	vector<string> features;

	vector<T> GLSZMData;
	extractGLSZMData(GLSZMData, GLSZMFeat);
	
	if (config.csvOutput == 1) {
		glszm2D.defineGLSZMFeatures(features);
		for (int i = 0; i < GLSZMData.size(); i++) {
			GLSZMCSV << "GLSZMFeatures2Davg" << "," << features[i] << ",";
			GLSZMCSV << GLSZMData[i];
			GLSZMCSV << "\n";
		}
	}
	else if (config.ontologyOutput == 1) {
		glszm2D.defineGLSZMFeaturesOntology(features);
		string featParamSpaceTable = config.outputFolder + "/FeatureParameterSpace_table.csv";
		char * featParamSpaceTableName = new char[featParamSpaceTable.size() + 1];
		std::copy(featParamSpaceTable.begin(), featParamSpaceTable.end(), featParamSpaceTableName);
		featParamSpaceTableName[featParamSpaceTable.size()] = '\0';

		ofstream featSpaceTable;
		featSpaceTable.open(featParamSpaceTableName, std::ios_base::app);
		parameterSpaceNr += 1;
		string parameterSpaceName = "FeatureParameterSpace_" + std::to_string(parameterSpaceNr);
		featSpaceTable << parameterSpaceName << "," << "2DAVG" << "," << config.imageSpaceName << "," << config.interpolationMethod << "\n";
		featSpaceTable.close();

		for (int i = 0; i < GLSZMData.size(); i++) {
			GLSZMCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			GLSZMCSV << GLSZMData[i] << "," << parameterSpaceName << "," << config.calculationSpaceName;
			GLSZMCSV << "\n";
		}

	}
	GLSZMCSV.close();
}

#endif // GLSZMFEATURES2DAVG_H_INCLUDED
