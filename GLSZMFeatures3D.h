#ifndef GLSZMFEATURES3D_H_INCLUDED
#define GLSZMFEATURES3D_H_INCLUDED

#include "GLSZMFeatures2D.h"

/*! \file */

/*!
The class GLSZMFeatures3D herites from the class GLSZMFeatures2D, because the feature calculations are the same.
Only the matrix calculation is different\n
All feature calculations are defined in the class GLRLMFeatures. \n
This class calculates a GLSZM matrix considering 13 neighbors(3D approach)\n
\n
For grey level size zone matrices, groups of connected voxels with a specific grey value and size are grouped.
A voxel is connected with another voxel if they have the same grey level.
*/

template <class T,  size_t R>
class GLSZMFeatures3D : public GLSZMFeatures2DMRG<T, R>{
    private:
		GLSZMFeatures2DMRG<T, R> GLSZM2D;
		GLRLMFeatures<T, R> glrlm;
		int maxZoneSize;
		vector<double> rowSums;
		vector<double> colSums;

        void extractGLSZMData3D(vector<T> &GLSZMData, GLSZMFeatures3D<T, R> GLSZMFeatures);
        boost::multi_array<double, 2> getGLSZMMatrix3D(boost::multi_array<T,R> inputMatrix, vector<T> vectorMatrElem);
        void fill3DGLSZMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix);
        int getBiggestZoneNr3D(vector<T> vectorMatrElem);

    public:
        void getNeighbors3D(boost::multi_array<T, R> &inputMatrix, T actElement, vector<vector< int> > &matrixIndices);
        void calculateAllGLSZMFeatures3D(GLSZMFeatures3D<T,R> &GLSZMFeat, Image<T, R> imageAttr, ConfigFile config);
        void writeCSVFileGLSZM3D(GLSZMFeatures3D<T,R> GLSZMFeat, string outputFolder);
		void writeOneFileGLSZM3D(GLSZMFeatures3D<T, R> GLSZMFeat, string outputFolder);

};





/*!
In the method fill3DGLSZMatrices the GLSZM matrices are filled using the matrix filled with the intensity values of the VOI.
@param[in] inputMatrix: the original matrix of the VOI
@param[in]: GLSZM is given as reference

It checks the neighborhood of every element of the input matrix. Is an element already considered as an element
of a neighborhood, it is set to NAN.\n
The size of the neighborhoods are stored in the matrix.
*/
template <class T, size_t R>
void GLSZMFeatures3D<T, R>::fill3DGLSZMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &GLSZMatrix) {
	T actualElement;
	vector<vector<int> > matrixIndices;
	vector<int> actualIndex;
	T actGreyElement;
	int actGreyIndex = 0;
	//do this grey level by grey level
	maxZoneSize = 0;
	int tempZoneSize = 0;
	//go element by element through the image
	//look at the neighbors of every element and check if they have the same grey value
	//set every element already seen to NAN
	for (int depth = 0; depth<inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row<inputMatrix.shape()[0]; row++) {
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				actualElement = inputMatrix[row][col][depth];
				if (!isnan(actualElement)) {
					actGreyIndex = glrlm.findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actualElement);
					inputMatrix[row][col][depth] = NAN;
					actualIndex.push_back(row);
					actualIndex.push_back(col);
					actualIndex.push_back(depth);
					matrixIndices.push_back(actualIndex);
					actualIndex.clear();
					getNeighbors3D(inputMatrix, actualElement, matrixIndices);
				}
				tempZoneSize = matrixIndices.size();
				if (tempZoneSize>maxZoneSize) {
					maxZoneSize = tempZoneSize;
				}
				if (matrixIndices.size() > 0) {

					GLSZMatrix[actGreyIndex][matrixIndices.size() - 1] += 1;
					matrixIndices.clear();
				}
			}
		}
	}

}



template <class T, size_t R>
int GLSZMFeatures3D<T, R>::getBiggestZoneNr3D(vector<T> vectorMatrElem) {
	T actGreyElement;
	vector<vector<int> > matrixIndices;
	vector<int> actualIndex;
	maxZoneSize = 0;
	vectorMatrElement.sort();

	/*int tempZoneSize = 0;
	int actGreyIndex = 0;

	for (int depth = 0; depth<inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row<inputMatrix.shape()[0]; row++) {
			for (int col = 0; col<inputMatrix.shape()[1]; col++) {
				actGreyElement = inputMatrix[row][col][depth];
				if (!isnan(actGreyElement)) {
					inputMatrix[row][col][depth] = NAN;
					actualIndex.push_back(row);
					actualIndex.push_back(col);
					actualIndex.push_back(depth);
					matrixIndices.push_back(actualIndex);
					actualIndex.clear();
					getNeighbors3D(inputMatrix, actGreyElement, matrixIndices);

				}
				tempZoneSize = matrixIndices.size();
				if (tempZoneSize>maxZoneSize) {
					maxZoneSize = tempZoneSize;
				}
				matrixIndices.clear();
			}			
		}
	}*/
	return maxZoneSize;
}




/*!
In the method getNeighbors3D the number of voxels in the biggest zone is determined
@param[in] inputMatrix: the original matrix of the VOI
@param[out]: int biggestZoneNr: the nr of voxels in the biggest zone

The function works as follows: \n
It gets the zone sizes for every intensity value and compares them \n
The biggest zone size is stored in the int value maxZoneNumber
*/
template <class T, size_t R>
void GLSZMFeatures3D<T, R>::getNeighbors3D(boost::multi_array<T, R> &inputMatrix, T actElement, vector<vector< int> > &matrixIndices){
    vector<int> actIndex;
    vector<int> newIndex;
    vector<int> tempIndex;
    int neighborDist = 1;
    T actMatrElement;
    int borderZ = 0;
    int directionX;
    int directionY;
    int directionZ = -1;
    vector<vector< int> > tempMatrixIndices = matrixIndices;
    int ang;
    //look at all indices which are saved in the vector matrix indices
    while(tempMatrixIndices.size()>0){
        //get the last element of the vector and delete it
        actIndex.clear();
        actIndex = tempMatrixIndices.back();
		tempMatrixIndices.pop_back();

        for(int i = 0; i < 9; i++){
			tempIndex = actIndex;

			//go in 45 steps from angle 360 to 0
            ang = 360 - i*45;			
            GLSZM2D.getALLXYDirections(directionX, directionY, ang);
            for(int j = 0; j < 3; j++){
                directionZ = -1 + j;
				
				if (directionZ != 0 || directionX != 0 || directionY != 0){
					tempIndex[0] = actIndex[0]+ directionX;
					tempIndex[1] =  actIndex[1]+directionY;
					tempIndex[2] = actIndex[2]+directionZ;

					if(tempIndex[0] > -1 && tempIndex[0] < inputMatrix.shape()[0] && tempIndex[1] > -1 && tempIndex[1] < inputMatrix.shape()[1] && tempIndex[2] > -1 && tempIndex[2]<inputMatrix.shape()[2]){
						actMatrElement = inputMatrix[tempIndex[0]][tempIndex[1]][tempIndex[2]];
	
						if (actMatrElement == actElement&&!isnan(actMatrElement)) {
							
							inputMatrix[tempIndex[0]][tempIndex[1]][tempIndex[2]] = NAN;
							tempMatrixIndices.push_back(tempIndex);
							matrixIndices.push_back(tempIndex);
						}
                    }
                }
			}
        }
    }
}



/*!
In the method getGLSZMMatrix3D the GLSZM matrices with the right size are generated and filled using the fill2DGLSZM function.
@param[in] inputMatrix: the original matrix of the VOI
@param[in] depth: number of actual slice for which the GLSZM should be calculated
@param[out]: filled GLSZM

This function only initiates a GLSZM with the right size.
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLSZMFeatures3D<T,R>::getGLSZMMatrix3D( boost::multi_array<T, R> inputMatrix, vector<T> vectorMatrElem){
    typedef boost::multi_array<double, 2>  GLSZMat;
	//get the maximal zone size in order to determine the size of the matrix
    //maxZoneSize=getBiggestZoneNr3D(inputMatrix);
	maxZoneSize = inputMatrix.shape()[0] * inputMatrix.shape()[1] * inputMatrix.shape()[2];
    int sizeMatrix= (this->diffGreyLevels).size();
    GLSZMat GLSZMatrix(boost::extents[sizeMatrix][maxZoneSize]);
    fill3DGLSZMatrices(inputMatrix, GLSZMatrix);

	GLSZMat smallGLSZM(boost::extents[sizeMatrix][maxZoneSize]);
	for (int row = 0; row < sizeMatrix; row++) {
		for (int col = 0; col < maxZoneSize; col++) {
			smallGLSZM[row][col] = GLSZMatrix[row][col];
		}
	}
    return smallGLSZM;

}

template <class T, size_t R>
void GLSZMFeatures3D<T, R>::calculateAllGLSZMFeatures3D(GLSZMFeatures3D<T,R> &GLSZMFeatures, Image<T,R> imageAttr, ConfigFile config){
    this->diffGreyLevels = imageAttr.diffGreyLevels;
	GLSZMFeatures.getConfigValues(config);
    boost::multi_array<double,2> GLSZM=GLSZMFeatures.getGLSZMMatrix3D(imageAttr.imageMatrix, imageAttr.vectorOfMatrixElements);

    double totalSum = GLSZMFeatures.calculateTotalSum(GLSZM);
    rowSums=GLSZMFeatures.calculateRowSums(GLSZM);
    colSums = GLSZMFeatures.calculateColSums(GLSZM);


    GLSZMFeatures.calculateShortRunEmphasis(rowSums, totalSum);
    GLSZMFeatures.calculateLongRunEmphasis(rowSums, totalSum);
    GLSZMFeatures.calculateLowGreyEmph(colSums, totalSum);
    GLSZMFeatures.calculateHighGreyEmph(colSums, totalSum);
    GLSZMFeatures.calculateShortRunLow(GLSZM, totalSum);
    GLSZMFeatures.calculateShortRunHigh(GLSZM, totalSum);
    GLSZMFeatures.calculateLongRunLowEmph(GLSZM, totalSum);
    GLSZMFeatures.calculateLongRunHighEmph(GLSZM, totalSum);
    GLSZMFeatures.calculateGreyNonUniformity(colSums, totalSum);
    GLSZMFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
    GLSZMFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
    GLSZMFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
    GLSZMFeatures.calculateRunPercentage3D(imageAttr.vectorOfMatrixElements, totalSum, 1);
    boost::multi_array<double,2> probMatrix = GLSZMFeatures.calculateProbMatrix(GLSZM, totalSum);
    double meanGrey = GLSZMFeatures.calculateMeanProbGrey(probMatrix);

    GLSZMFeatures.calculateGreyLevelVar(probMatrix, meanGrey);

    double meanRun = GLSZMFeatures.calculateMeanProbRun(probMatrix);
    GLSZMFeatures.calculateRunLengthVar(probMatrix, meanRun);
    GLSZMFeatures.calculateRunEntropy(probMatrix);

  }

  template <class T, size_t R>
void GLSZMFeatures3D<T, R>::extractGLSZMData3D(vector<T> &GLSZMData, GLSZMFeatures3D<T, R> GLSZMFeatures){

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
void GLSZMFeatures3D<T, R>::writeCSVFileGLSZM3D(GLSZMFeatures3D<T,R> GLSZMFeat, string outputFolder)
{
    string csvName = outputFolder + "_GLSZMFeatures3D.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream GLSZMCSV;
    GLSZMCSV.open (name);
    vector<string> features;
    GLSZM2D.defineGLSZMFeatures(features);

    vector<T> GLSZMData;
    extractGLSZMData3D(GLSZMData, GLSZMFeat);
    for(int i = 0; i< GLSZMData.size(); i++){
        GLSZMCSV <<"GLSZMFeatures3D"<< ","<< features[i] <<",";
        GLSZMCSV << GLSZMData[i];
        GLSZMCSV << "\n";
    }
    GLSZMCSV.close();
}

template <class T, size_t R>
void GLSZMFeatures3D<T, R>::writeOneFileGLSZM3D(GLSZMFeatures3D<T, R> GLSZMFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream GLSZMCSV;
	GLSZMCSV.open(name, std::ios_base::app);
	vector<string> features;
	GLSZM2D.defineGLSZMFeatures(features);

	vector<T> GLSZMData;
	extractGLSZMData3D(GLSZMData, GLSZMFeat);
	for (int i = 0; i< GLSZMData.size(); i++) {
		GLSZMCSV << "GLSZMFeatures3D" << "," << features[i] << ",";
		GLSZMCSV << GLSZMData[i];
		GLSZMCSV << "\n";
	}
	GLSZMCSV.close();
}

#endif // GLSZMFEATURES3D_H_INCLUDED
