#ifndef GLDZMFEATURES2DAVG_H_INCLUDED
#define GLDZMFEATURES2DAVG_H_INCLUDED

#include "GLDZMFeatures3D.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodOperatorImageFunction.h"
#include "getNeighborhoodMatrices.h"
/*! \file */
/*!
The class GLDZMFeatures2DAVG is the class of the Grey Level Distance Zone Matrices, it inheritates from 
the class GLDZMFeatures2D. \n
For further explanation look at GLDZMFeatures2D file.\n
For every slice a GLDZM is calculated and from every of this matrices all features are extracted. \n
Then the average value of all features is calculated.
*/

template <class T,  size_t R=3>
class GLDZMFeatures2DAVG : public GLSZMFeatures2DMRG<T,R>{
    private:
		vector<T> diagonalProbabilities;
		vector<T> crossProbabilities;
		vector<T> sumProbRows;
		vector<T> sumProbCols;
		int totalNrZones;
		vector<float> rowSums;
		vector<float> colSums;
		GLSZMFeatures2DMRG<T, R> GLSZM2D;
        GLDZMFeatures2D<T,R> GLDZM2D;
		GLDZMFeatures3D<T, R> GLDZM3D;
		void fillConvolutionalVectorGLD(vector<T> &convolutionalVector);
		int checkNeighbors(boost::multi_array<T, R> &distanceMap, boost::multi_array<T, R> &inputMatrix, vector< int> actIndex, int actualDistance);
		void extractGLDZMData2DAVG(vector<T> &gldzmData, GLDZMFeatures2DAVG<T, R> gldzmFeatures);
        boost::multi_array<float, 2> fillMatrix(boost::multi_array<T,R> inputMatrix, boost::multi_array<T, R> distanceMap, boost::multi_array<float, 2>  &gldzmat, int depth);
        boost::multi_array<float, 2> getMatrix( boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap, int depth);
		int checkNeighborsNAN(boost::multi_array<T, R> &inputMatrix, boost::multi_array<T, R> &tempMatrix, vector< int> actIndex, int actDist);
     public:
		 boost::multi_array<T, R> generateDistanceMap(boost::multi_array<T, R> inputMatrix, Image<T, R> imageAttr, boost::multi_array<T, R> &distanceMap, ConfigFile config);

        void writeCSVFileGLDZM2DAVG(GLDZMFeatures2DAVG<T,R> gldzmFeat, string outputFolder);
		void writeOneFileGLDZM2DAVG(GLDZMFeatures2DAVG<T, R> gldzmFeat, ConfigFile config, int &parameterSpaceNr);
        void calculateAllGLDZMFeatures2DAVG(GLDZMFeatures2DAVG<T,R> &gldzmFeat, Image<T,R> imageAttr, boost::multi_array<T, R> distanceMap, ConfigFile config);
};

template <class T, size_t R>
void GLDZMFeatures2DAVG<T, R>::fillConvolutionalVectorGLD(vector<T> &convolutionalVector) {
	for (int row = 0; row < 3; row++) {
		for (int col = 0; col < 3; col++) {

					if (col == row ) {
						convolutionalVector.push_back(float(0.0));
					}
					else if (col == 0 && row == 2) {
						convolutionalVector.push_back(float(0.0));
					}
					else if (col == 2 && row == 0) {
						convolutionalVector.push_back(float(0.0));
					}
					
					else {
						convolutionalVector.push_back(1);
					}

				}
		}

}
 

/*!
In the method generateDistanceMap the distance map is generated, taking
the matrix of the VOI as input. \n
According to the position of a value in the matrix, the corresponding distance is saved in the distance map
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: distance map
*/
template <class T, size_t R>
boost::multi_array<T, R> GLDZMFeatures2DAVG<T, R>::generateDistanceMap(boost::multi_array<T,R> inputMatrix, Image<T, R> imageAttr, boost::multi_array<T,R> &distanceMap, ConfigFile config){
	boost::multi_array<float, 3> distanceMapTmp(boost::extents[distanceMap.shape()[0]][distanceMap.shape()[1]][distanceMap.shape()[2]]);
	
	//construct convolutional kernel
	vector<float> convolutionalVector;
	fillConvolutionalVectorGLD(convolutionalVector);
	float *kernelArray = &convolutionalVector[0];
	unsigned int dimKernel[] = { 3,3,1 };
	float voxelSize[] = { 1,1,1 };
	ImageType::Pointer kernel = converArray2Image(kernelArray, dimKernel, voxelSize);
	//assign the image matrix
	boost::multi_array<T, 3> actualMatrix = inputMatrix;
	//how many grey levels do I have
	int sizeGreyLevels = (imageAttr.diffGreyLevels).size();
	//the original (not discretize image!)
	ImageType::Pointer actImage = imageAttr.image;

	ImageType::Pointer imageNew = constructEmptyNewImage(imageAttr, actualMatrix, 0);
	ImageType::IndexType pixelIndex;
	//fill image with 0 and 1 so transform it to a matrix
	for (int depth = 0; depth < actualMatrix.shape()[2] ; depth++) {
		for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
			for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
				if (!std::isnan(actualMatrix[row - 1][col - 1][depth ])) {
					pixelIndex[0] = row;
					pixelIndex[1] = col;
					pixelIndex[2] = depth;
					imageNew->SetPixel(pixelIndex, 1);
				}
			}
		}
	}
	//make copy of mask
	ImageType::Pointer maskOriginal = imageNew.GetPointer();
	const typename ImageType::RegionType& imageRegion = imageNew->GetLargestPossibleRegion();
	const typename ImageType::SizeType& imageSize = imageRegion.GetSize();
	//boost::multi_array<T, 3> sumMatrix;
	//make a convolution of the image
	ImageType::Pointer sumImage = convolutionImage(imageNew, kernel);
	Image<T, 3> imageElement(imageSize[0], imageSize[1], imageSize[2]);
	boost::multi_array<T, 3> sumMatrix = imageElement.get3Dimage(sumImage, maskOriginal, config);
	int dist = 1;
	int stop = 0;
	int flag;
	//go slice by slice (2D)
	for (int depth = 0; depth < distanceMap.shape()[2]; depth++) {
		dist = 1;
		stop = 0;
		while (stop == 0) {
			flag = 0;
			//get values of image and mask
			sumImage = convolutionImage(imageNew, kernel);
			sumMatrix = imageElement.get3Dimage(sumImage, maskOriginal, config);
			sumImage = nullptr;
			//imageNew = constructEmptyNewImage(imageAttr, actualMatrix, 1);
		
			//for every matrix element we access the neighborhood
			for (int row = 0; row < distanceMap.shape()[0]; row++) {
				for (int col = 0; col < distanceMap.shape()[1]; col++) {

					pixelIndex[0] = row + 1;
					pixelIndex[1] = col + 1;
					pixelIndex[2] = depth;
					
					//get actual Element if it is the centre of a whole neighborhood
					if (sumMatrix[row + 1][col + 1][depth ] == 4) {
						flag = 1;
					}
					
					else if (sumMatrix[row + 1][col + 1][depth ] < 4) {//&&sumMatrix[row][col][depth]!=0 && !std::isnan(sumMatrix[row][col][depth])) {
						distanceMap[row][col][depth] = dist;
						imageNew->SetPixel(pixelIndex, 0);
					}
				}
			}
			if (flag == 0) {
				stop = 1;
			}
			dist++;
		}
	
	}
	return distanceMap;
}


template <class T, size_t R>
int GLDZMFeatures2DAVG<T, R>::checkNeighborsNAN(boost::multi_array<T, R> &inputMatrix, boost::multi_array<T, R> &tempMatrix, vector< int> actIndex, int actDist) {
	T actMatrElement;
	vector<int> tempIndex;
	int directionX;
	int directionY;
	int distance;
	int actMatrElementDist;
	//look at all indices which are saved in the vector matrix indices
	tempIndex = actIndex;
	int actPoint, distBorder2, ang;
	float index;
	vector<int> actDirection;
	int row, col;
	vector<vector<int>> directionVector = { { 1,0 },{ -1,0 },{ 0,1 },{ 0,-1 } };
	//get the points at border of bounding box
	actPoint = std::min(actIndex[0], actIndex[1]);
	distBorder2 = inputMatrix.shape()[1] - actIndex[1] - 1;
	actPoint = std::min(actPoint, distBorder2);
	distBorder2 = inputMatrix.shape()[0] - actIndex[0] - 1;
	actPoint = std::min(actPoint, distBorder2);
	if (actPoint == 0 && !isnan(tempMatrix[actIndex[0]][actIndex[1]][actIndex[2]])) {
		distance = 1;
		tempMatrix[actIndex[0]][actIndex[1]][actIndex[2]] = NAN;
	}
	
	for (int i = 0; i < directionVector.size(); i++) {
		actDirection = directionVector[i];
		row = actDirection[0];
		col = actDirection[1];
		if ( actIndex[0] + row < inputMatrix.shape()[0] && actIndex[0] + row >-1 && actIndex[1] + col<inputMatrix.shape()[1] && actIndex[1] + col>-1 && !std::isnan(inputMatrix[actIndex[0]][actIndex[1]][actIndex[2]])) {
			index = inputMatrix[actIndex[0] + row][actIndex[1] + col][actIndex[2]] ;
			if (std::isnan(index)) {
				distance = actDist;
				tempMatrix[actIndex[0]][actIndex[1]][actIndex[2]] = NAN;
				break;
			}
		}
			
	}

	return distance;
}

/*!
In the method fillMatrix the GLDZM matrix is filled, taking the original matrix of the VOI as input. \n
The GLDZM matrix is given as reference and filled in the function
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[in]: boost::multi_array<T, 3> GLDZM: GLDZM matrix
*/
template <class T, size_t R>
boost::multi_array<float, 2> GLDZMFeatures2DAVG<T, R>::fillMatrix(boost::multi_array<T,R> inputMatrix, boost::multi_array<T, R> distanceMap, boost::multi_array<float, 2>  &gldzmat, int depth){
	
	vector<vector<int> > matrixIndices;
    vector<int> actualIndex;
    T actualElement;
    int minDistance;
    totalNrZones=0;
	int actualGreyIndex;

    for(int row = 0; row<inputMatrix.shape()[0]; row++){
        for(int col=0; col<inputMatrix.shape()[1];col++){
			actualElement = inputMatrix[row][col][depth];
			if (!isnan(actualElement)) {
				actualGreyIndex = GLSZM2D.findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actualElement);
                inputMatrix[row][col][depth]=NAN;
                actualIndex.push_back(row);
                actualIndex.push_back(col);
                actualIndex.push_back(depth);
                matrixIndices.push_back(actualIndex);
                actualIndex.clear();
                GLSZM2D.getNeighbors(inputMatrix, actualElement, matrixIndices);
            }
            if(matrixIndices.size()>0){
				minDistance = GLDZM2D.getMinimalDistance(distanceMap, matrixIndices);
                matrixIndices.clear();
				if (minDistance > 0) {
					gldzmat[actualGreyIndex][minDistance - 1] += 1;
				}
            }
        }
    }
    return gldzmat;
}


template <class T, size_t R>
int GLDZMFeatures2DAVG<T, R>::checkNeighbors(boost::multi_array<T, R> &distanceMap, boost::multi_array<T, R> &tempMatrix, vector< int> actIndex, int actualDistance) {
	T actMatrElement;
	int distance = 0;
	int directionX;
	int directionY;
	int ang;
	int actMatrElementDist;
	//look at all indices which are saved in the vector matrix indices
	for (int i = 0; i < 8; i++) {
		ang = 360 - i * 45;
		GLSZM2D.getALLXYDirections(directionX, directionY, ang);
		if (actIndex[0] + directionX > -1 && actIndex[0] + directionX < tempMatrix.shape()[0] && actIndex[1] + directionY > -1 && actIndex[1] + directionY < tempMatrix.shape()[1]) {
			actMatrElementDist = distanceMap[actIndex[0] + directionX][actIndex[1] + directionY][actIndex[2]];
			if (actMatrElementDist == actualDistance) {
				distance = actualDistance + 1;
				tempMatrix[actIndex[0]][actIndex[1]][actIndex[2]] = NAN;
				break;
			}
		}
	}
	return distance;
}
/*!
In the method getMatrix the GLDZM matrix is generated and filled using the function fillMatrix.
The function is mainly used get the size of the GLDZM matrix. \n
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: boost::multi_array<T, 3> GLDZM: GLDZM matrix
*/
template <class T, size_t R>
boost::multi_array<float, 2> GLDZMFeatures2DAVG<T,R>::getMatrix( boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap, int depth){
    typedef boost::multi_array<float, 2>  gldzmat;

    int sizeGreyLevels = (this->diffGreyLevels).size();

	int nrCols = std::min(std::ceil(float(inputMatrix.shape()[0]) / 2), std::ceil(float(inputMatrix.shape()[1]) / 2))+1;
    gldzmat GLDZMatrix(boost::extents[sizeGreyLevels][nrCols]);
	
	fillMatrix(inputMatrix, distanceMap, GLDZMatrix, depth);

    return GLDZMatrix;

}




template <class T, size_t R>
void GLDZMFeatures2DAVG<T, R>::calculateAllGLDZMFeatures2DAVG(GLDZMFeatures2DAVG<T, R> &gldzmFeatures, Image<T, R> imageAttr, boost::multi_array<T, R> distanceMap, ConfigFile config) {
	this->diffGreyLevels = imageAttr.diffGreyLevels;
	int totalDepth = imageAttr.imageMatrix.shape()[2];
	gldzmFeatures.getConfigValues(config);
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
	T sumRunLengthNonUniformityNorm = 0;
	T sumRunLengthNonUniformity = 0;
	T sumRunPercentage = 0;
	T sumGreyLevelVar = 0;
	T sumRunLengthVar = 0;
	T sumRunEntropy = 0;
	
    for(int depth = 0; depth < totalDepth; depth++){
		boost::multi_array<float, 2> GLDZM = gldzmFeatures.getMatrix(imageAttr.imageMatrix, distanceMap, depth);
        float totalSum = gldzmFeatures.calculateTotalSum(GLDZM);
        rowSums=gldzmFeatures.calculateRowSums(GLDZM);
        colSums = gldzmFeatures.calculateColSums(GLDZM);

        gldzmFeatures.calculateShortRunEmphasis(rowSums, totalSum);
        sumShortRunEmphasis += this->shortRunEmphasis;
        gldzmFeatures.calculateLongRunEmphasis(rowSums, totalSum);
        sumLongRunEmphasis += this->longRunEmphasis;
        gldzmFeatures.calculateLowGreyEmph(colSums, totalSum);
        sumLowGreyEmph += this->lowGreyEmph;
        gldzmFeatures.calculateHighGreyEmph(colSums, totalSum);
        sumHighGreyEmph += this->highGreyEmph;
        gldzmFeatures.calculateShortRunLow(GLDZM, totalSum);
        sumShortRunLow += this->shortRunLow;
        gldzmFeatures.calculateShortRunHigh(GLDZM, totalSum);
        sumShortRunHigh += this->shortRunHigh;
        gldzmFeatures.calculateLongRunLowEmph(GLDZM, totalSum);
        sumLongRunLowEmph += this->longRunLowEmph;
        gldzmFeatures.calculateLongRunHighEmph(GLDZM, totalSum);
        sumLongRunHighEmph += this->longRunHighEmph;
        gldzmFeatures.calculateGreyNonUniformity(colSums, totalSum);
        sumGreyNonUniformity += this->greyNonUniformity;
        gldzmFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
        sumGreyNonUniformityNorm += this->greyNonUniformityNorm;
        gldzmFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
        sumRunLengthNonUniformityNorm += this->runLengthNonUniformityNorm;
        gldzmFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
        sumRunLengthNonUniformity += this->runLengthNonUniformity;
        gldzmFeatures.calculateRunPercentage(imageAttr.imageMatrix, depth, totalSum, 1);
        sumRunPercentage += this->runPercentage;
        boost::multi_array<float,2> probMatrix = gldzmFeatures.calculateProbMatrix(GLDZM, totalSum);
        float meanGrey = gldzmFeatures.calculateMeanProbGrey(probMatrix);

        gldzmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);
        sumGreyLevelVar += this->greyLevelVar;
        float meanRun = gldzmFeatures.calculateMeanProbRun(probMatrix);
        gldzmFeatures.calculateRunLengthVar(probMatrix, meanRun);
        sumRunLengthVar += this->runLengthVar;
        gldzmFeatures.calculateRunEntropy(probMatrix);
        sumRunEntropy += this->runEntropy;
    }
    this->shortRunEmphasis = sumShortRunEmphasis/totalDepth;
    this->longRunEmphasis = sumLongRunEmphasis/totalDepth;
    this->lowGreyEmph = sumLowGreyEmph/totalDepth;
    this->highGreyEmph = sumHighGreyEmph/totalDepth;
    this->shortRunLow = sumShortRunLow/totalDepth;
    this->shortRunHigh = sumShortRunHigh/totalDepth;
    this->longRunLowEmph = sumLongRunLowEmph/totalDepth;
    this->longRunHighEmph = sumLongRunHighEmph/totalDepth;
    this->greyNonUniformity = sumGreyNonUniformity/totalDepth;
    this->greyNonUniformityNorm = sumGreyNonUniformityNorm/totalDepth;
    this->runLengthNonUniformityNorm = sumRunLengthNonUniformityNorm/totalDepth;
    this->runLengthNonUniformity = sumRunLengthNonUniformity/totalDepth;
    this->runPercentage = sumRunPercentage/totalDepth;
    this->runLengthVar = sumRunLengthVar/totalDepth;
    this->runEntropy = sumRunEntropy/totalDepth;
}

template <class T, size_t R>
void GLDZMFeatures2DAVG<T, R>::extractGLDZMData2DAVG(vector<T> &gldzmData, GLDZMFeatures2DAVG<T, R> gldzmFeatures){

    gldzmData.push_back(gldzmFeatures.shortRunEmphasis);
    gldzmData.push_back(gldzmFeatures.longRunEmphasis);
    gldzmData.push_back(gldzmFeatures.lowGreyEmph);
    gldzmData.push_back(gldzmFeatures.highGreyEmph);
    gldzmData.push_back(gldzmFeatures.shortRunLow);
    gldzmData.push_back(gldzmFeatures.shortRunHigh);
    gldzmData.push_back(gldzmFeatures.longRunLowEmph);
    gldzmData.push_back(gldzmFeatures.longRunHighEmph);
    gldzmData.push_back(gldzmFeatures.greyNonUniformity);
    gldzmData.push_back(gldzmFeatures.greyNonUniformityNorm);
    gldzmData.push_back(gldzmFeatures.runLengthNonUniformity);   //check if I have to give the matrix
    gldzmData.push_back(gldzmFeatures.runLengthNonUniformityNorm);
    gldzmData.push_back(gldzmFeatures.runPercentage);
    gldzmData.push_back(gldzmFeatures.greyLevelVar);
    gldzmData.push_back(gldzmFeatures.runLengthVar);
    gldzmData.push_back(gldzmFeatures.runEntropy);

}

template <class T, size_t R>
void GLDZMFeatures2DAVG<T, R>::writeCSVFileGLDZM2DAVG(GLDZMFeatures2DAVG<T,R> gldzmFeat, string outputFolder)
{
    string csvName = outputFolder + "_gldzmFeatures2Davg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream gldzmCSV;
    gldzmCSV.open (name);
    vector<string> features;
	GLDZM2D.defineGLDZMFeatures(features);

    vector<T> gldzmData;
    extractGLDZMData2DAVG(gldzmData, gldzmFeat);
	for (int i = 0; i< gldzmData.size(); i++) {
		gldzmCSV << "gldzmFeatures2Davg" << "," << features[i] << ",";
		gldzmCSV << gldzmData[i];
		gldzmCSV << "\n";
	}
	
    gldzmCSV.close();
}

template <class T, size_t R>
void GLDZMFeatures2DAVG<T, R>::writeOneFileGLDZM2DAVG(GLDZMFeatures2DAVG<T, R> gldzmFeat, ConfigFile config, int &parameterSpaceNr) {
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

	ofstream gldzmCSV;
	gldzmCSV.open(name, std::ios_base::app);
	vector<string> features;
	

	vector<T> gldzmData;
	extractGLDZMData2DAVG(gldzmData, gldzmFeat);
	
	if (config.csvOutput == 1) {
		GLDZM2D.defineGLDZMFeatures(features);
		for (int i = 0; i < gldzmData.size(); i++) {
			gldzmCSV << "gldzmFeatures2Davg" << "," << features[i] << ",";
			gldzmCSV << gldzmData[i];
			gldzmCSV << "\n";
		}
	}
	else if (config.ontologyOutput == 1) {
		features.clear();
		GLDZM2D.defineGLDZMFeatures(features);
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

		for (int i = 0; i < gldzmData.size(); i++) {
			gldzmCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			gldzmCSV << gldzmData[i] << "," << parameterSpaceName << "," << config.calculationSpaceName;
			gldzmCSV << "\n";
		}

	}
	gldzmCSV.close();

}

#endif // GLDZMFEATURES2DAVG_H_INCLUDED
