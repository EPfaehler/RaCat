#ifndef GLDZMFEATURES2DAVG_H_INCLUDED
#define GLDZMFEATURES2DAVG_H_INCLUDED

#include "GLDZMFeatures3D.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodOperatorImageFunction.h"
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
		vector<double> rowSums;
		vector<double> colSums;
		GLSZMFeatures2DMRG<T, R> GLSZM2D;
        GLDZMFeatures2D<T,R> GLDZM2D;
		GLDZMFeatures3D<T, R> GLDZM3D;
		int checkNeighbors(boost::multi_array<T, R> &distanceMap, boost::multi_array<T, R> &inputMatrix, vector< int> actIndex, int actualDistance);
		void extractGLDZMData2DAVG(vector<T> &gldzmData, GLDZMFeatures2DAVG<T, R> gldzmFeatures);
        boost::multi_array<double, 2> fillMatrix(boost::multi_array<T,R> inputMatrix, boost::multi_array<T, R> distanceMap, boost::multi_array<double, 2>  &gldzmat, int depth);
        boost::multi_array<double, 2> getMatrix( boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap, int depth);
		int checkNeighborsNAN(boost::multi_array<T, R> &inputMatrix, boost::multi_array<T, R> &tempMatrix, vector< int> actIndex, int actDist);
     public:
		 boost::multi_array<T, R> generateDistanceMap(boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> &distanceMap);

        void writeCSVFileGLDZM2DAVG(GLDZMFeatures2DAVG<T,R> gldzmFeat, string outputFolder);
		void writeOneFileGLDZM2DAVG(GLDZMFeatures2DAVG<T, R> gldzmFeat, string outputFolder);
        void calculateAllGLDZMFeatures2DAVG(GLDZMFeatures2DAVG<T,R> &gldzmFeat, Image<T,R> imageAttr, boost::multi_array<T, R> distanceMap, ConfigFile config);
};


 

/*!
In the method generateDistanceMap the distance map is generated, taking
the matrix of the VOI as input. \n
According to the position of a value in the matrix, the corresponding distance is saved in the distance map
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: distance map
*/
template <class T, size_t R>
boost::multi_array<T, R> GLDZMFeatures2DAVG<T, R>::generateDistanceMap(boost::multi_array<T,R> inputMatrix, boost::multi_array<T,R> &distanceMap){
	boost::multi_array<float, 3> distanceMapTmp(boost::extents[distanceMap.shape()[0]][distanceMap.shape()[1]][distanceMap.shape()[2]]);
	int dist1 = 1;
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		
		for (int col = 0; col < inputMatrix.shape()[1]; col++) {
			dist1 = 1;
			for (int row = 0; row < inputMatrix.shape()[0]; row++) {
				if (!std::isnan(inputMatrix[row][col][depth])) {
					distanceMap[row][col][depth] = dist1;
					dist1++;
				}
				else {
					dist1 = 1;
				}
			}
		}
	}
	for (int depth = 0; depth <inputMatrix.shape()[2]; depth++) {
		for (int col = 0; col < inputMatrix.shape()[1]; col++) {
			dist1 = 1;
			for (int row = 0; row < inputMatrix.shape()[0]; row++) {
				if (!std::isnan(inputMatrix[inputMatrix.shape()[0] - row - 1][col][depth])) {
					distanceMap[inputMatrix.shape()[0] - row - 1][col][depth] = dist1 < distanceMap[inputMatrix.shape()[0] - row - 1][col][depth] ? dist1 : distanceMap[inputMatrix.shape()[0] - row - 1][col][depth];
					dist1++;
				}
				else {
					dist1 = 1;
				}
			}
		}
	}
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row < inputMatrix.shape()[0]; row++) {
			dist1 = 1;
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				if (!std::isnan(inputMatrix[row][col][depth])) {
					distanceMapTmp[row][col][depth] = dist1;
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
			dist1 = 1;
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				if (!std::isnan(inputMatrix[row][inputMatrix.shape()[1] - col - 1][depth])) {
					distanceMapTmp[row][inputMatrix.shape()[1] - col - 1][depth] = dist1 < distanceMapTmp[row][inputMatrix.shape()[1] - col - 1][depth] ? dist1 : distanceMapTmp[row][inputMatrix.shape()[1] - col - 1][depth];
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
					distanceMap[row][col][depth] = distanceMap[row][col][depth] < distanceMapTmp[row][col][depth] ? distanceMap[row][col][depth] : distanceMapTmp[row][col][depth];
				}
				else {
					distanceMap[row][col][depth] = NAN;
				}
			}
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
boost::multi_array<double, 2> GLDZMFeatures2DAVG<T, R>::fillMatrix(boost::multi_array<T,R> inputMatrix, boost::multi_array<T, R> distanceMap, boost::multi_array<double, 2>  &gldzmat, int depth){
	
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
boost::multi_array<double, 2> GLDZMFeatures2DAVG<T,R>::getMatrix( boost::multi_array<T, R> inputMatrix, boost::multi_array<T, R> distanceMap, int depth){
    typedef boost::multi_array<double, 2>  gldzmat;

    int sizeGreyLevels = (this->diffGreyLevels).size();

	int nrCols = std::min(std::ceil(double(inputMatrix.shape()[0]) / 2), std::ceil(double(inputMatrix.shape()[1]) / 2))+1;
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
		boost::multi_array<double, 2> GLDZM = gldzmFeatures.getMatrix(imageAttr.imageMatrix, distanceMap, depth);
        double totalSum = gldzmFeatures.calculateTotalSum(GLDZM);
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
        boost::multi_array<double,2> probMatrix = gldzmFeatures.calculateProbMatrix(GLDZM, totalSum);
        double meanGrey = gldzmFeatures.calculateMeanProbGrey(probMatrix);

        gldzmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);
        sumGreyLevelVar += this->greyLevelVar;
        double meanRun = gldzmFeatures.calculateMeanProbRun(probMatrix);
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
    for(int i = 0; i< gldzmData.size(); i++){
        gldzmCSV <<"gldzmFeatures2Davg"<<","<< features[i] <<",";
        gldzmCSV << gldzmData[i];
        gldzmCSV << "\n";
    }
    gldzmCSV.close();
}

template <class T, size_t R>
void GLDZMFeatures2DAVG<T, R>::writeOneFileGLDZM2DAVG(GLDZMFeatures2DAVG<T, R> gldzmFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream gldzmCSV;
	gldzmCSV.open(name, std::ios_base::app);
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

#endif // GLDZMFEATURES2DAVG_H_INCLUDED
