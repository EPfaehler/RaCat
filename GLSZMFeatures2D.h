#ifndef GLSZMFEATURES2DMRG_H_INCLUDED
#define GLSZMFEATURES2DMRG_H_INCLUDED

#include "GLRLMFeatures.h"


/*! \file */

/*!
The class GLSZMFeatures2DMRG herites from the class GLRLMFeatures, because the feature calculations are the same.
Only the matrix calculation is different\n
All feature calculations are defined in the class GLRLMFeatures. \n
This class only contains the calculations of the 2D matrix. \n
\n
For grey level size zone matrices, groups of connected voxels with a specific grey value and size are grouped. 
A voxel is connected with another voxel if they have the same grey level. 
*/

template <class T,  size_t R>
class GLSZMFeatures2DMRG : public GLRLMFeatures<T, R>{
    private:
		//the biggest number of  voxels one zone (it is needed to determine the matrix size)		
		int maxZoneSize;
		vector<double> rowSums;
		vector<double> colSums;
		GLRLMFeatures<T, R> glrlm;
		int getBiggestZoneNr(boost::multi_array<T, R> inputMatrix);
		boost::multi_array<double, 2> getGLSZMMatrix(boost::multi_array<T, R> inputMatrix);
        void fill2DGLSZMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &glcMatrix);
		void extractGLSZMData(vector<T> &GLSZMData, GLSZMFeatures2DMRG<T, R> GLSZMFeatures);

    public:

        void defineGLSZMFeatures(vector<string> &features);

        void getNeighbors(boost::multi_array<T, R> &inputMatrix, T actElement, vector<vector< int> > &matrixIndices);
        void getALLXYDirections(int &directionX, int &directionY, int angle);
        void calculateAllGLSZMFeatures2DMRG(GLSZMFeatures2DMRG<T,R> &GLSZMFeat, boost::multi_array<T,R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElem, ConfigFile config);
        void writeCSVFileGLSZM(GLSZMFeatures2DMRG<T,R> GLSZMFeat, string outputFolder);
		void writeOneFileGLSZM(GLSZMFeatures2DMRG<T, R> GLSZMFeat, string outputFolder);

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
void GLSZMFeatures2DMRG<T, R>::fill2DGLSZMatrices(boost::multi_array<T, R> inputMatrix, boost::multi_array<double, 2> &GLSZMatrix){
	T actElement;
	vector<vector<int> > matrixIndices;
	vector<int> actualIndex;
	int actGreyIndex;
	maxZoneSize = 0;
	int tempZoneSize = 0;
	//go element by element through the image
	//look at the neighbors of every element and check if they have the same grey value
	//set every element already seen to NAN
	for (int depth = 0; depth<inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row<inputMatrix.shape()[0]; row++) {
			for (int col = 0; col<inputMatrix.shape()[1]; col++) {
				actElement = inputMatrix[row][col][depth];
				if (!isnan(actElement)) {
					actGreyIndex = glrlm.findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actElement);
					inputMatrix[row][col][depth] = NAN;
					actualIndex.push_back(row);
					actualIndex.push_back(col);
					actualIndex.push_back(depth);
					matrixIndices.push_back(actualIndex);
					actualIndex.clear();
					getNeighbors(inputMatrix, actElement, matrixIndices);
				}
				//the number of the voxels of the actual zone 
				tempZoneSize = matrixIndices.size();
				//if this number is bigger than the actual max zone size
				if (tempZoneSize>maxZoneSize) {
					maxZoneSize = tempZoneSize;
				}
				if (matrixIndices.size()>0) {
					GLSZMatrix[actGreyIndex][matrixIndices.size() - 1] += 1;
					matrixIndices.clear();
				}
			}

		}
	}
}





/*!
In the method getBiggestZoneNr the number of voxels in the biggest zone is determined
@param[in] inputMatrix: the original matrix of the VOI
@param[out]: int biggestZoneNr: the nr of voxels in the biggest zone

The function works as follows: \n
It gets the zone sizes for every intensity value and compares them \n
The biggest zone size is stored in the int value maxZoneNumber
*/
template <class T, size_t R>
int GLSZMFeatures2DMRG<T, R>::getBiggestZoneNr(boost::multi_array<T, R> inputMatrix){
	T actualElement;
	vector<vector<int> > matrixIndices;
	vector<int> actualIndex;
	maxZoneSize = 0;
	int tempZoneSize = 0;
	int actualGreyIndex;
	//calculate biggest zone for every intensity value
		//get the neighbors for every element in the matrix
		for (int depth = 0; depth<inputMatrix.shape()[2]; depth++) {
			for (int row = 0; row<inputMatrix.shape()[0]; row++) {
				for (int col = 0; col<inputMatrix.shape()[1]; col++) {
					actualElement = inputMatrix[row][col][depth];
					if (!isnan(actualElement)) {
						actualGreyIndex = findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actualElement);
						inputMatrix[row][col][depth] = NAN;
						actualIndex.push_back(row);
						actualIndex.push_back(col);
						actualIndex.push_back(depth);
						matrixIndices.push_back(actualIndex);
						actualIndex.clear();
						getNeighbors(inputMatrix, actualElement, matrixIndices);
					}
				}
				//the number of the voxels of the actual zone 
				tempZoneSize = matrixIndices.size();
				//if this number is bigger than the actual max zone size
				if (tempZoneSize>maxZoneSize) {
					maxZoneSize = tempZoneSize;
				}

			}
			matrixIndices.clear();
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
void GLSZMFeatures2DMRG<T, R>::getNeighbors(boost::multi_array<T, R> &inputMatrix, T actElement, vector<vector< int> > &matrixIndices){
    vector<int> actIndex;
    vector<int> newIndex;
    vector<int> tempIndex;
    int neighborDist =1;
    T actMatrElement;

    int directionX;
    int directionY;
    vector<vector< int> > tempMatrixIndices = matrixIndices;
    int ang;
    //look at all indices which are saved in the vector matrix indices
    while(tempMatrixIndices.size()>0){
        //get the last element of the vector and delete it
        actIndex = tempMatrixIndices.back();
        tempMatrixIndices.pop_back();
        
        for(int i=0; i<8; i++){
			//temporal index
			tempIndex = actIndex;
            ang = 360-i*45;
            getALLXYDirections(directionX, directionY, ang);
			tempIndex[0] += directionX;
			tempIndex[1] += directionY;
            if(tempIndex[0]>-1 && tempIndex[0]<inputMatrix.shape()[0] && tempIndex[1]>-1 && tempIndex[1]<inputMatrix.shape()[1]){
				
				actMatrElement = inputMatrix[tempIndex[0]][tempIndex[1]][actIndex[2]];
                if(actMatrElement == actElement){
                    inputMatrix[tempIndex[0]][tempIndex[1]][actIndex[2]] = NAN;

                    tempMatrixIndices.push_back(tempIndex);
                    matrixIndices.push_back(tempIndex);
                }
            }
        }
    }
}


//

template <class T, size_t R>
void GLSZMFeatures2DMRG<T, R>::getALLXYDirections(int &directionX, int &directionY, int angle){
    if(angle==360){
        directionX=0;
        directionY=1;
    }
    else if(angle==315){
        directionX=1;
        directionY=1;
    }
    else if(angle==270){
        directionX=1;
        directionY=0;
    }
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
@param[out]: filled GLSZM

This function only initiates a GLSZM with the right size.
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLSZMFeatures2DMRG<T,R>::getGLSZMMatrix( boost::multi_array<T, R> inputMatrix){
    typedef boost::multi_array<double, 2>  GLSZMat;
    //maxZoneSize=getBiggestZoneNr(inputMatrix);
	maxZoneSize = inputMatrix.shape()[0] * inputMatrix.shape()[1];
    int sizeMatrix= (this->diffGreyLevels).size();
    GLSZMat GLSZMatrix(boost::extents[sizeMatrix][maxZoneSize]);
    fill2DGLSZMatrices(inputMatrix, GLSZMatrix);

	GLSZMat smallGLSZM(boost::extents[sizeMatrix][maxZoneSize]);
	for (int row = 0; row < sizeMatrix; row++) {
		for (int col = 0; col < maxZoneSize; col++) {
			smallGLSZM[row][col] = GLSZMatrix[row][col];
		}
	}
    return smallGLSZM;
}

template <class T, size_t R>
void GLSZMFeatures2DMRG<T, R>::calculateAllGLSZMFeatures2DMRG(GLSZMFeatures2DMRG<T,R> &GLSZMFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElem, ConfigFile config){
	this->diffGreyLevels = diffGrey;
	GLSZMFeatures.getConfigValues(config);
	boost::multi_array<double, 2> GLSZM = GLSZMFeatures.getGLSZMMatrix(inputMatrix);

	double totalSum = GLSZMFeatures.calculateTotalSum(GLSZM);

	rowSums = GLSZMFeatures.calculateRowSums(GLSZM);
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
	GLSZMFeatures.calculateRunPercentage3D(vectorMatrElem, totalSum, 1);
	boost::multi_array<double, 2> probMatrix = GLSZMFeatures.calculateProbMatrix(GLSZM, totalSum);
	double meanGrey = GLSZMFeatures.calculateMeanProbGrey(probMatrix);

	GLSZMFeatures.calculateGreyLevelVar(probMatrix, meanGrey);
	double meanRun = GLSZMFeatures.calculateMeanProbRun(probMatrix);
	GLSZMFeatures.calculateRunLengthVar(probMatrix, meanRun);
	GLSZMFeatures.calculateRunEntropy(probMatrix);
}

template <class T, size_t R>
void GLSZMFeatures2DMRG<T, R>::extractGLSZMData(vector<T> &GLSZMData, GLSZMFeatures2DMRG<T, R> GLSZMFeatures){

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
void GLSZMFeatures2DMRG<T, R>::writeCSVFileGLSZM(GLSZMFeatures2DMRG<T,R> GLSZMFeat, string outputFolder)
{
	string csvName = outputFolder + "_GLSZMFeatures2Dvmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';
    ofstream GLSZMCSV;
    GLSZMCSV.open (name);
    vector<string> features;
    GLSZMFeat.defineGLSZMFeatures(features);
    vector<T> GLSZMData;
    extractGLSZMData(GLSZMData, GLSZMFeat);
    for(int i = 0; i< GLSZMData.size(); i++){
        GLSZMCSV <<"GLSZMFeatures2Dvmrg"<<","<< features[i] <<",";
        GLSZMCSV << GLSZMData[i];
        GLSZMCSV << "\n";
    }
    GLSZMCSV.close();
}

template <class T, size_t R>
void GLSZMFeatures2DMRG<T, R>::writeOneFileGLSZM(GLSZMFeatures2DMRG<T, R> GLSZMFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';
	ofstream GLSZMCSV;
	GLSZMCSV.open(name, std::ios_base::app);
	vector<string> features;
	GLSZMFeat.defineGLSZMFeatures(features);
	vector<T> GLSZMData;
	extractGLSZMData(GLSZMData, GLSZMFeat);
	for (int i = 0; i< GLSZMData.size(); i++) {
		GLSZMCSV << "GLSZMFeatures2Dvmrg" << "," << features[i] << ",";
		GLSZMCSV << GLSZMData[i];
		GLSZMCSV << "\n";
	}
	GLSZMCSV.close();
}

template <class T, size_t R>
void GLSZMFeatures2DMRG<T, R>::defineGLSZMFeatures(vector<string> &features){
    features.push_back("small zone emphasis");
    features.push_back("Large zone emphasis");
    features.push_back("Low grey level zone emphasis");
    features.push_back("High grey level zone emphasis");
    features.push_back("Small zone low grey level emphasis");
    features.push_back("Small zone high grey level emphasis");
    features.push_back("Large zone low grey level emphasis");
    features.push_back("Large zone high grey level emphasis");
    features.push_back("Grey level non uniformity GLSZM");
    features.push_back("Grey level non uniformity normalized GLSZM");
    features.push_back("Zone size non uniformity");
    features.push_back("Zone size non uniformity normalized");
    features.push_back("Zone percentage GLSZM");
    features.push_back("Grey level variance GLSZM");
    features.push_back("Zone size variance");
    features.push_back("Zone size entropy");

}


#endif // GLSZMFEATURES2D_H_INCLUDED
