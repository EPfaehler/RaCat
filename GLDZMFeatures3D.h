#ifndef GLDZMFEATURES3D_H_INCLUDED
#define GLDZMFEATURES3D_H_INCLUDED

#include "GLDZMFeatures2DMRG.h"
#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include <itkIsoContourDistanceImageFilter.h>
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

		vector<T> diagonalProbabilities;
		vector<T> crossProbabilities;
		vector<T> sumProbRows;
		vector<T> sumProbCols;
		int totalNrZones;
		vector<double> rowSums;
		vector<double> colSums;

        void extractGLDZMData3D(vector<T> &gldzmData, GLDZMFeatures3D<T, R> gldzmFeatures);
		
		boost::multi_array<T, R> generateDistanceMap3D(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap2D, boost::multi_array<T, R> &distanceMap);
        void fillMatrix3D(boost::multi_array<T,R> inputMatrix, boost::multi_array<T, R> distanceMatrix, boost::multi_array<double, 2>  &gldzmat);
        boost::multi_array<double, 2> getMatrix3D( boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMatrix);
     public:
        void writeCSVFileGLDZM3D(GLDZMFeatures3D<T,R> gldzmFeat, string outputFolder);
		void writeOneFileGLDZM3D(GLDZMFeatures3D<T, R> gldzmFeat, string outputFolder);
        void calculateAllGLDZMFeatures3D(GLDZMFeatures3D<T,R> &gldzmFeat, boost::multi_array<T, R> distanceMap, Image<T,R> imageAttr, ConfigFile config);
		int checkNeighbors3DNAN(boost::multi_array<T, R> distanceMap, boost::multi_array<T, R> &tempMatrix, vector< int> actIndex, int actDistance);
		int checkNeighbors3D(boost::multi_array<T, R> &distanceMap, boost::multi_array<T, R> &tempMatrix, vector< int> actIndex, int actualDistance);

};

template <class T, size_t R>
boost::multi_array<T, R> GLDZMFeatures3D<T, R>::generateDistanceMap3D(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap2D, boost::multi_array<T, R> &distanceMap3D) {
	int dist1;
	
	
	for (int row = 0; row < inputMatrix.shape()[0]; row++) {
		for (int col = 0; col < inputMatrix.shape()[1]; col++) {
			dist1 = 1;
			for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
				if (!std::isnan(inputMatrix[row][col][depth])) {
					distanceMap3D[row][col][depth] = dist1;
					dist1++;
				}
				else {
					dist1 = 1;
				}
			}
		}
	}

	for (int row = 0; row < inputMatrix.shape()[0]; row++) {
		for (int col = 0; col < inputMatrix.shape()[1]; col++) {
			dist1 = 1;
			for (int depth = 0; depth <inputMatrix.shape()[2]; depth++) {
				if (!std::isnan(inputMatrix[row][col][inputMatrix.shape()[2] - depth - 1])) {
					distanceMap3D[row][col][inputMatrix.shape()[2] - depth - 1] = dist1 < distanceMap3D[row][col][inputMatrix.shape()[2] - depth - 1] ? dist1 : distanceMap3D[row][col][inputMatrix.shape()[2] - depth - 1];
					dist1++;
				}
				else {
					dist1 = 1;
				}
			}
		}
	}
	for (int depth = 0; depth <inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row < inputMatrix.shape()[0]; row++) {
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {

				if (!std::isnan(inputMatrix[row][col][depth])) {
					distanceMap3D[row][col][depth] = distanceMap3D[row][col][depth] < distanceMap2D[row][col][depth] ? distanceMap3D[row][col][depth] : distanceMap2D[row][col][depth];
				}
				else {
					distanceMap3D[row][col][depth] = NAN;
				}
			}
		}
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
void GLDZMFeatures3D<T, R>::fillMatrix3D(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap, boost::multi_array<double, 2>  &gldzmat) {
	vector<vector<int> > matrixIndices;
	vector<int> actualIndex;
	T actElement;
	int minDistance;
	
	for (int i = 0; i < (this->diffGreyLevels).size(); i++) {
		actElement = (this->diffGreyLevels)[i];

		for (int depth = 0; depth<inputMatrix.shape()[2]; depth++) {
			for (int row = 0; row<inputMatrix.shape()[0]; row++) {
				for (int col = 0; col<inputMatrix.shape()[1]; col++) {
					if (inputMatrix[row][col][depth] == actElement) {
						inputMatrix[row][col][depth] = NAN;
						actualIndex.push_back(row);
						actualIndex.push_back(col);
						actualIndex.push_back(depth);
						matrixIndices.push_back(actualIndex);
						actualIndex.clear();
						GLSZM3D.getNeighbors3D(inputMatrix, actElement, matrixIndices);
					}
				}
				if (matrixIndices.size()>0) {
					minDistance = GLDZM2D.getMinimalDistance(distanceMap, matrixIndices);
					matrixIndices.clear();
					gldzmat[i][minDistance - 1] += 1;
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
boost::multi_array<double, 2> GLDZMFeatures3D<T,R>::getMatrix3D( boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap){
    typedef boost::multi_array<double, 2>  gldzmat;
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
		generateDistanceMap3D(imageAttr.imageMatrix, distanceMap2D, distanceMap);
	}
	else {
		generateDistanceMap3D(imageAttr.imageMatrixOriginal, distanceMap2D, distanceMap);
	}
	
    boost::multi_array<double,2> GLDZM=gldzmFeatures.getMatrix3D(imageAttr.imageMatrix, distanceMap);

    double totalSum = gldzmFeatures.calculateTotalSum(GLDZM);
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
    boost::multi_array<double,2> probMatrix = gldzmFeatures.calculateProbMatrix(GLDZM, totalSum);
    double meanGrey = gldzmFeatures.calculateMeanProbGrey(probMatrix);

    gldzmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);

    double meanRun = gldzmFeatures.calculateMeanProbRun(probMatrix);
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
void GLDZMFeatures3D<T, R>::writeOneFileGLDZM3D(GLDZMFeatures3D<T, R> gldzmFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream gldzmCSV;
	gldzmCSV.open(name, std::ios_base::app);
	vector<string> features;
	GLDZM2D.defineGLDZMFeatures(features);

	vector<T> gldzmData;
	extractGLDZMData3D(gldzmData, gldzmFeat);
	for (int i = 0; i< gldzmData.size(); i++) {
		gldzmCSV << "gldzmFeatures3D" << "," << features[i] << ",";
		gldzmCSV << gldzmData[i];
		gldzmCSV << "\n";
	}
	gldzmCSV.close();
}

#endif // GLDZMFEATURES3D_H_INCLUDED
