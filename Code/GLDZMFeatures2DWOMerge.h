#ifndef GLDZMFEATURES2DWOMERGE_H_INCLUDED
#define GLDZMFEATURES2DWOMERGE_H_INCLUDED

#include "GLDZMFeatures2D.h"

/*! \file */
/*!
The class GLDZMFeatures2DWOMerge is the class of the Grey Level Distance Zone Matrices, it inheritates from 
the class GLDZMFeatures2D. \n
For further explanation look at GLDZMFeatures2D file.\n
For every slice a GLDZM is calculated and from every of this matrices all features are extracted. \n
Then the average value of all features is calculated.
*/

template <class T,  size_t R=3>
class GLDZMFeatures2DWOMerge : public GLSZMFeatures2D<T,R>{
    private:
		vector<T> diagonalProbabilities;
		vector<T> crossProbabilities;
		vector<T> sumProbRows;
		vector<T> sumProbCols;
		int totalNrZones;
		vector<double> rowSums;
		vector<double> colSums;
        GLSZMFeatures2D<T, R> GLSZM2D;
        GLDZMFeatures2D<T,R> GLDZM2D;
        void extractGLDZMData2DWOMerge(vector<T> &gldzmData, GLDZMFeatures2DWOMerge<T, R> gldzmFeatures);
        boost::multi_array<T, R> generateDistanceMap(boost::multi_array<T,R> inputMatrix);
        boost::multi_array<double, 2> fillMatrix(boost::multi_array<T,R> inputMatrix, boost::multi_array<double, 2>  &gldzmat, int depth);
        boost::multi_array<double, 2> getMatrix( boost::multi_array<T, R> inputMatrix, int depth);

     public:
        void writeCSVFileGLDZM2DWOMerge(GLDZMFeatures2DWOMerge<T,R> gldzmFeat, string outputFolder);
		void writeOneFileGLDZM2DWOMerge(GLDZMFeatures2DWOMerge<T, R> gldzmFeat, string outputFolder);
        void calculateAllGLDZMFeatures2DWOMerge(GLDZMFeatures2DWOMerge<T,R> &gldzmFeat, boost::multi_array<T,R> inputMatrix, vector<T> diffGrey, ConfigFile config);
};


 

/*!
In the method generateDistanceMap the distance map is generated, taking
the matrix of the VOI as input. \n
According to the position of a value in the matrix, the corresponding distance is saved in the distance map
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: distance map
*/
template <class T, size_t R>
boost::multi_array<T, R> GLDZMFeatures2DWOMerge<T, R>::generateDistanceMap(boost::multi_array<T,R> inputMatrix){
	typedef boost::multi_array<T, R> multi_array;
	multi_array distanceMap(boost::extents[inputMatrix.shape()[0]][inputMatrix.shape()[1]][inputMatrix.shape()[2]]);
	vector<int> actualIndex;
	for (int depth = 0; depth<inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row<inputMatrix.shape()[0]; row++) {
			for (int col = 0; col<inputMatrix.shape()[1]; col++) {
				actualIndex.push_back(row);
				actualIndex.push_back(col);
				actualIndex.push_back(depth);

				distanceMap[row][col][depth] = GLDZM2D.checkNeighbors(distanceMap, inputMatrix, actualIndex);
				actualIndex.clear();
			}
		}
	}
	return distanceMap;
}

/*!
In the method fillMatrix the GLDZM matrix is filled, taking the original matrix of the VOI as input. \n
The GLDZM matrix is given as reference and filled in the function
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[in]: boost::multi_array<T, 3> GLDZM: GLDZM matrix
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLDZMFeatures2DWOMerge<T, R>::fillMatrix(boost::multi_array<T,R> inputMatrix, boost::multi_array<double, 2>  &gldzmat, int depth){
    
	boost::multi_array<T,R> distanceMap = generateDistanceMap(inputMatrix);
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
                GLSZM2D.getNeighbors(inputMatrix, actualElement, matrixIndices, false);
            }
            if(matrixIndices.size()>0){
				minDistance = GLDZM2D.getMinimalDistance(distanceMap, matrixIndices);
                matrixIndices.clear();
                gldzmat[actualGreyIndex][minDistance-1]+=1;
            }
        }
    }
    
    return gldzmat;
}

/*!
In the method getMatrix the GLDZM matrix is generated and filled using the function fillMatrix.
The function is mainly used get the size of the GLDZM matrix. \n
@param[in]: boost::multi_array<T, 3> inputMatrix: original matrix of the VOI
@param[out]: boost::multi_array<T, 3> GLDZM: GLDZM matrix
*/
template <class T, size_t R>
boost::multi_array<double, 2> GLDZMFeatures2DWOMerge<T,R>::getMatrix( boost::multi_array<T, R> inputMatrix, int depth){
    typedef boost::multi_array<double, 2>  gldzmat;

    int sizeGreyLevels = (this->diffGreyLevels).size();
	int nrCols = std::min(std::ceil(double(inputMatrix.shape()[0]) / 2), std::ceil(double(inputMatrix.shape()[1]) / 2))+1;
    gldzmat GLDZMatrix(boost::extents[sizeGreyLevels][nrCols]);
    fillMatrix(inputMatrix, GLDZMatrix, depth);
    return GLDZMatrix;

}




template <class T, size_t R>
void GLDZMFeatures2DWOMerge<T, R>::calculateAllGLDZMFeatures2DWOMerge(GLDZMFeatures2DWOMerge<T,R> &gldzmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, ConfigFile config){
    this->diffGreyLevels = diffGrey;
    int totalDepth = inputMatrix.shape()[2];
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
        boost::multi_array<double,2> GLDZM=gldzmFeatures.getMatrix(inputMatrix, depth);
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
        gldzmFeatures.calculateRunPercentage(inputMatrix, depth, totalSum, 1);
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
void GLDZMFeatures2DWOMerge<T, R>::extractGLDZMData2DWOMerge(vector<T> &gldzmData, GLDZMFeatures2DWOMerge<T, R> gldzmFeatures){

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
void GLDZMFeatures2DWOMerge<T, R>::writeCSVFileGLDZM2DWOMerge(GLDZMFeatures2DWOMerge<T,R> gldzmFeat, string outputFolder)
{
    string csvName = outputFolder + "/gldzmFeatures2Davg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream gldzmCSV;
    gldzmCSV.open (name);
    vector<string> features;
    GLDZM2D.defineGLDZMFeatures(features);

    vector<T> gldzmData;
    extractGLDZMData2DWOMerge(gldzmData, gldzmFeat);
    for(int i = 0; i< gldzmData.size(); i++){
        gldzmCSV <<"gldzmFeatures2Davg"<<","<< features[i] <<",";
        gldzmCSV << gldzmData[i];
        gldzmCSV << "\n";
    }
    gldzmCSV.close();
}

template <class T, size_t R>
void GLDZMFeatures2DWOMerge<T, R>::writeOneFileGLDZM2DWOMerge(GLDZMFeatures2DWOMerge<T, R> gldzmFeat, string outputFolder) {
	string csvName = outputFolder + "/radiomicsFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream gldzmCSV;
	gldzmCSV.open(name, std::ios_base::app);
	vector<string> features;
	GLDZM2D.defineGLDZMFeatures(features);

	vector<T> gldzmData;
	extractGLDZMData2DWOMerge(gldzmData, gldzmFeat);
	for (int i = 0; i< gldzmData.size(); i++) {
		gldzmCSV << "gldzmFeatures2Davg" << "," << features[i] << ",";
		gldzmCSV << gldzmData[i];
		gldzmCSV << "\n";
	}
	gldzmCSV.close();

}

#endif // GLDZMFEATURES2DWOMERGE_H_INCLUDED
