#ifndef GLDZMFEATURES3D_H_INCLUDED
#define GLDZMFEATURES3D_H_INCLUDED

#include "GLDZMFeatures2DMRG.h"
#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include <itkIsoContourDistanceImageFilter.h>
#include "getNeighborhoodMatrices.h"
/*! \file */
/*!
The class GLDZM is the class of the Grey Level Distance Zone Matrices for the 3D approach. \n
It combines the grey level size zone matrices with a distance map. Voxels are considered as connected,
when they have the same grey value. \n
The distance to the edge is also defined according to 6-connectedness. The distance of a voxel to the
outer border is defined as the number of edges that have to be crossed to reach the edge of the VOI. \n
*/


template <class T,  size_t R=3>
class GLDZMFeatures3D : public GLDZMFeatures2D<T,R>{
    private:
        GLSZMFeatures2DMRG<T, R> GLSZM2D;
        GLDZMFeatures2D<T,R> GLDZM2D;
        GLSZMFeatures3D<T, R> GLSZM3D;
		GLRLMFeatures<T, R> glrlm;
		vector<T> diagonalProbabilities;
		vector<T> crossProbabilities;
		vector<T> sumProbRows;
		vector<T> sumProbCols;
		int totalNrZones;
		vector<float> rowSums;
		vector<float> colSums;

        void extractGLDZMData3D(vector<T> &gldzmData, GLDZMFeatures3D<T, R> gldzmFeatures);
		
		boost::multi_array<T, R> generateDistanceMap3D(Image<T,R> imageAttr, boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap2D, boost::multi_array<T, R> &distanceMap, ConfigFile config);
        void fillMatrix3D(boost::multi_array<T,R> inputMatrix, boost::multi_array<T, R> distanceMatrix, boost::multi_array<float, 2>  &gldzmat);
        boost::multi_array<float, 2> getMatrix3D( boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMatrix);
		int getMaxDistance(boost::multi_array<T, R> inputMatrix);
		void fillConvolutionalVectorGLD(vector<T> &convolutionalVector, int is3D);
     public:
        void writeCSVFileGLDZM3D(GLDZMFeatures3D<T,R> gldzmFeat, string outputFolder);
		void writeOneFileGLDZM3D(GLDZMFeatures3D<T, R> gldzmFeat, ConfigFile config, int &parameterSpaceNr);
        void calculateAllGLDZMFeatures3D(GLDZMFeatures3D<T,R> &gldzmFeat, boost::multi_array<T, R> distanceMap, Image<T,R> imageAttr, ConfigFile config);

};

template <class T, size_t R>
int GLDZMFeatures3D<T, R>::getMaxDistance(boost::multi_array<T, R> inputMatrix) {
	vector<float> distances;
	distances.push_back(ceil(inputMatrix.shape()[0] / 2));
	distances.push_back(ceil(inputMatrix.shape()[1] / 2));
	distances.push_back(ceil(inputMatrix.shape()[2] / 2));
	int maximalDistance = *max_element(distances.begin(), distances.end());
	return maximalDistance;
}

template <class T, size_t R>
void GLDZMFeatures3D<T, R>::fillConvolutionalVectorGLD(vector<T> &convolutionalVector, int is3D) {
	if (is3D == 1) {
		for (int depth = 0; depth < 3; depth++) {
			for (int row = 0; row < 3; row++) {
				for (int col = 0; col < 3; col++) {

					if (depth == 0 && col == 0 && row == 0) {
						convolutionalVector.push_back(float(0.0));
					}
					else if ((depth == 0 && col == 0)|| (depth == 0 && col == 2)||(depth ==0&&col ==1&&row!=1)) {
						convolutionalVector.push_back(float(0.0));
					}
					else if (depth == 2 && col == 2 && row == 2) {
						convolutionalVector.push_back(float(0.0));
					}
					else if ((depth == 2 && col == 0) || (depth == 2 && col == 2) || (depth == 2 && col == 1 && row != 1)) {
						convolutionalVector.push_back(float(0.0));
					}
					else if (depth == 1 && col == 2 && row == 0||(depth == 1 && col == 0 && row == 2)) {
						convolutionalVector.push_back(float(0.0));
					}
					else if ((depth == 1 && col == row)) {
						convolutionalVector.push_back(float(0.0));
					}
					else {
						convolutionalVector.push_back(1);
					}
					
				}
			}

		}
	}
	else {
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				if (col*row != 1) {
					convolutionalVector.push_back(1);
				}
				else {
					convolutionalVector.push_back(float(0.0));
				}
			}
		}
	}
}

template <class T, size_t R>
boost::multi_array<T, R> GLDZMFeatures3D<T, R>::generateDistanceMap3D(Image<T, R> imageAttr, boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap2D, boost::multi_array<T, R> &distanceMap3D, ConfigFile config) {
	int dist = 1;
	//construct convolutional kernel
	vector<float> convolutionalVector;
	fillConvolutionalVectorGLD(convolutionalVector, 1);
	float *kernelArray = &convolutionalVector[0];
	unsigned int dimKernel[] = { 3,3,3 };
	float voxelSize[] = { 1,1,1 };
	ImageType::Pointer kernel = converArray2Image(kernelArray, dimKernel, voxelSize);
	//assign the image matrix
	boost::multi_array<T, 3> actualMatrix = inputMatrix;
	//how many grey levels do I have
	int sizeGreyLevels = (imageAttr.diffGreyLevels).size();
	//the original (not discretize image!)
	ImageType::Pointer actImage = imageAttr.image;
	
	ImageType::Pointer imageNew = constructEmptyNewImage(imageAttr, actualMatrix, 1);
	ImageType::IndexType pixelIndex;
	//fill image with 0 and 1 so transform it to a matrix
	for (int depth = 1; depth < actualMatrix.shape()[2] + 1; depth++) {
		for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
			for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
				if (!std::isnan(actualMatrix[row-1][col-1][depth-1])) {
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
	boost::multi_array<T, 3> sumMatrix = imageElement.get3Dimage(sumImage, maskOriginal, config);;
	
	int stop = 0;
	int flag;
	while (stop ==0) {
		//get values of image and mask
		sumImage = convolutionImage(imageNew, kernel);
		sumMatrix = imageElement.get3Dimage(sumImage, maskOriginal, config);
		sumImage = nullptr;
		//imageNew = constructEmptyNewImage(imageAttr, actualMatrix, 1);
		flag = 0;
		//for every matrix element we access the neighborhood
		for (int depth = 0; depth < distanceMap3D.shape()[2]; depth++) {
			for (int row =0; row < distanceMap3D.shape()[0]; row++) {
				for (int col = 0; col < distanceMap3D.shape()[1]; col++) {


					pixelIndex[0] = row+1;
					pixelIndex[1] = col+1;
					pixelIndex[2] = depth+1;
					
						//get actual Element if it is the centre of a whole neighborhood
				//std:cout << sumMatrix[row+1][col+1][depth+1] << " ";
					if (sumMatrix[row+1][col+1][depth+1] == 6) {
						flag=1;

					}
					else if (sumMatrix[row + 1][col + 1][depth + 1] <6 ){//&&sumMatrix[row][col][depth]!=0 && !std::isnan(sumMatrix[row][col][depth])) {
						distanceMap3D[row][col][depth] = dist;
					//std:cout << sumMatrix[row+1][col+1][depth+1] << " " << distanceMap3D[row][col][depth] << " ";
						//->SetPixel(pixelIndex, dist);
						imageNew->SetPixel(pixelIndex, 0);


					}
					

				}//std::cout << "new" << std::endl;
			}
		}
		
		if (flag == 0) {
			stop = 1;
		}
		dist++;
	}
	
	
	return distanceMap3D;

}


/*!
In the method fillMatrix3D the GLDZM matrix is filled, taking the original matrix of the VOI as input. \n
The GLDZM matrix is given as reference and filled in the function
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[in]: boost::multi_array<T, 3> GLDZM: GLDZM matrix
*/
template <class T, size_t R>
void GLDZMFeatures3D<T, R>::fillMatrix3D(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap, boost::multi_array<float, 2>  &gldzmat) {
	vector<vector<int> > matrixIndices;
	vector<int> actualIndex;
	T actElement;
	int minDistance;
	int actGreyIndex;

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
					GLSZM3D.getNeighbors3D(inputMatrix, actElement, matrixIndices);
					
				}
				if (matrixIndices.size()>0) {
					minDistance = GLDZM2D.getMinimalDistance(distanceMap, matrixIndices);
					matrixIndices.clear();
					gldzmat[actGreyIndex][minDistance - 1] += 1;
					totalNrZones += 1;
				}
			}
		}
	}
}


/*!
In the method getMatrix3D the GLDZM matrix is generated and filled using the function fillMatrix.
The function is mainly used get the size of the GLDZM matrix. \n
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: boost::multi_array<T, 3> GLDZM: GLDZM matrix
*/
template <class T, size_t R>
boost::multi_array<float, 2> GLDZMFeatures3D<T,R>::getMatrix3D( boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap){
    typedef boost::multi_array<float, 2>  gldzmat;
	//get the number of different grey levels in the VOI
    int sizeGreyLevels = (this->diffGreyLevels).size();
	//in order to get the maximal distance to the  borders, calculate the distance from
	//the center voxel to the border
    int colSize = ceil(double(inputMatrix.shape()[0])/2);
    if(ceil(double(inputMatrix.shape()[1])/2)<colSize){
        colSize = ceil(double(inputMatrix.shape()[1])/2);
    }
    if(ceil(double(inputMatrix.shape()[2])/2)<colSize){
        colSize = ceil(double(inputMatrix.shape()[2])/2);
    }
    gldzmat GLDZMatrix(boost::extents[sizeGreyLevels][colSize]);
	fillMatrix3D(inputMatrix, distanceMap, GLDZMatrix);

    return GLDZMatrix;
}

template <class T, size_t R>
void GLDZMFeatures3D<T, R>::calculateAllGLDZMFeatures3D(GLDZMFeatures3D<T,R> &gldzmFeatures, boost::multi_array<T, R> distanceMap2D, Image<T,R> imageAttr, ConfigFile config){
	gldzmFeatures.getConfigValues(config);
	this->diffGreyLevels = imageAttr.diffGreyLevels;
	boost::multi_array<T, R> distanceMap(boost::extents[imageAttr.imageMatrix.shape()[0]][imageAttr.imageMatrix.shape()[1]][imageAttr.imageMatrix.shape()[2]]);
	if (config.useReSegmentation == 0 && config.excludeOutliers== 0) {
		generateDistanceMap3D(imageAttr, imageAttr.imageMatrix, distanceMap2D, distanceMap, config);
	}
	else {
		generateDistanceMap3D(imageAttr, imageAttr.imageMatrixOriginal, distanceMap2D, distanceMap, config);
	}
    boost::multi_array<float,2> GLDZM=gldzmFeatures.getMatrix3D(imageAttr.imageMatrix, distanceMap);
	
	float totalSum = gldzmFeatures.calculateTotalSum(GLDZM);
    rowSums=gldzmFeatures.calculateRowSums(GLDZM);
    colSums = gldzmFeatures.calculateColSums(GLDZM);

    gldzmFeatures.calculateShortRunEmphasis(rowSums, totalSum);
    gldzmFeatures.calculateLongRunEmphasis(rowSums, totalSum);
    gldzmFeatures.calculateLowGreyEmph(colSums, totalSum);
    gldzmFeatures.calculateHighGreyEmph(colSums, totalSum);
    gldzmFeatures.calculateShortRunLow(GLDZM, totalSum);
    gldzmFeatures.calculateShortRunHigh(GLDZM, totalSum);
    gldzmFeatures.calculateLongRunLowEmph(GLDZM, totalSum);
    gldzmFeatures.calculateLongRunHighEmph(GLDZM, totalSum);
    gldzmFeatures.calculateGreyNonUniformity(colSums, totalSum);
    gldzmFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
    gldzmFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
    gldzmFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
    gldzmFeatures.calculateRunPercentage3D(imageAttr.vectorOfMatrixElements, totalSum, 1);
    boost::multi_array<float,2> probMatrix = gldzmFeatures.calculateProbMatrix(GLDZM, totalSum);
	float meanGrey = gldzmFeatures.calculateMeanProbGrey(probMatrix);

    gldzmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);

	float meanRun = gldzmFeatures.calculateMeanProbRun(probMatrix);
    gldzmFeatures.calculateRunLengthVar(probMatrix, meanRun);
    gldzmFeatures.calculateRunEntropy(probMatrix);
}

template <class T, size_t R>
void GLDZMFeatures3D<T, R>::extractGLDZMData3D(vector<T> &gldzmData, GLDZMFeatures3D<T, R> gldzmFeatures){

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
void GLDZMFeatures3D<T, R>::writeCSVFileGLDZM3D(GLDZMFeatures3D<T,R> gldzmFeat, string outputFolder)
{
    string csvName = outputFolder + "_gldzmFeatures3D.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream gldzmCSV;
    gldzmCSV.open (name);
    vector<string> features;
    GLDZM2D.defineGLDZMFeatures(features);

    vector<T> gldzmData;
    extractGLDZMData3D(gldzmData, gldzmFeat);
    for(int i = 0; i< gldzmData.size(); i++){
        gldzmCSV <<"gldzmFeatures3D"<<","<< features[i] <<",";
        gldzmCSV << gldzmData[i];
        gldzmCSV << "\n";
    }
    gldzmCSV.close();
}

template <class T, size_t R>
void GLDZMFeatures3D<T, R>::writeOneFileGLDZM3D(GLDZMFeatures3D<T, R> gldzmFeat, ConfigFile config, int &parameterSpaceNr) {
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
	extractGLDZMData3D(gldzmData, gldzmFeat);
	
	if (config.csvOutput == 1) {
		GLDZM2D.defineGLDZMFeatures(features);
		for (int i = 0; i < gldzmData.size(); i++) {
			gldzmCSV << "gldzmFeatures3D" << "," << features[i] << ",";
			gldzmCSV << gldzmData[i];
			gldzmCSV << "\n";
		}
	}
	else if (config.ontologyOutput == 1) {
		features.clear();
		GLDZM2D.defineGLDZMFeaturesOntology(features);
		string featParamSpaceTable = config.outputFolder + "/FeatureParameterSpace_table.csv";
		char * featParamSpaceTableName = new char[featParamSpaceTable.size() + 1];
		std::copy(featParamSpaceTable.begin(), featParamSpaceTable.end(), featParamSpaceTableName);
		featParamSpaceTableName[featParamSpaceTable.size()] = '\0';

		ofstream featSpaceTable;
		featSpaceTable.open(featParamSpaceTableName, std::ios_base::app);
		parameterSpaceNr += 1;
		string parameterSpaceName = "FeatureParameterSpace_" + std::to_string(parameterSpaceNr);
		featSpaceTable << parameterSpaceName << "," << "3Dmrg" << "," << config.imageSpaceName << "," << config.interpolationMethod << "\n";
		featSpaceTable.close();

		for (int i = 0; i < gldzmData.size(); i++) {
			gldzmCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			gldzmCSV << gldzmData[i] << "," << parameterSpaceName << "," << config.calculationSpaceName;
			gldzmCSV << "\n";
		}

	}
	gldzmCSV.close();
}

#endif // GLDZMFEATURES3D_H_INCLUDED
