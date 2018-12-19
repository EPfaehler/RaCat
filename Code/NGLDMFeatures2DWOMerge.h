#ifndef NGLDMFEATURES2DWOMERGE_H_INCLUDED
#define NGLDMFEATURES2DWOMERGE_H_INCLUDED

 #include "NGLDMFeatures.h"

/*! \file */
/*!
The class NGLDM2DWOMerge inherites from the class NGLDMFeatures. Here the 
matrices are calculated slice by slice. For every slice the feature values are calculated.\n
*/

template <class T,  size_t R>
class NGLDMFeatures2DWOMerge : public NGLDMFeatures<T,R>{
    private:
		int dist;
		int coarseParam;

		vector<T> sumProbRows;
		vector<T> sumProbCols;
		vector<double> rowSums;
		vector<double> colSums;

        NGLDMFeatures<T, R> ngldm;
        void extractNGLDMData2DWOMerge(vector<T> &ngldmData, NGLDMFeatures2DWOMerge<T, R> ngldmFeatures);
        boost::multi_array<double, 2> getMatrix(boost::multi_array<T,R> inputMatrix, int depth);
        int getNeighborGreyLevels(boost::multi_array<T, R> inputMatrix, vector<int> actualIndex);

    public:
        //double dependenceCountEnergy;
        void writeCSVFileNGLDM(NGLDMFeatures2DWOMerge<T,R> ngldmFeat, string outputFolder);
		void writeOneFileNGLDM(NGLDMFeatures2DWOMerge<T, R> ngldmFeat, string outputFolder);
        void calculateAllNGLDMFeatures2DWOMerge2D(NGLDMFeatures2DWOMerge<T,R> &ngldmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, ConfigFile config);

};

/*!
\brief getMatrix
@param boost multi array inputMatrix: matrix filled with intensity values
@param[out] boost multi array: filled NGLD matrix

The function fills the NGLDMatrix with the corresponding values using the function getNeighborGreyLevels. It
checks voxel by voxel the neighborhood in the distance that is set by the user.
*/
template <class T, size_t R>
boost::multi_array<double, 2> NGLDMFeatures2DWOMerge<T, R>::getMatrix(boost::multi_array<T,R> inputMatrix, int depth){
    typedef boost::multi_array<double, 2>  ngldmat;
    vector<int> actualIndex;
    T actualElement;
    int sizeGreyLevels = (this->diffGreyLevels).size();
    T actualGreyLevel;
	int actualGreyIndex;
    int ngldmnr;
    //define the NGLDMarices; col.size=9 because we have 8 neighbors
    ngldmat NGLDMatrix(boost::extents[sizeGreyLevels][9]);

    for(int row =0; row<inputMatrix.shape()[0]; row++){
		for(int col =0; col<inputMatrix.shape()[1]; col++){
			actualElement = inputMatrix[row][col][depth];
            if(!std::isnan(actualElement)){
				actualIndex.push_back(row);
                actualIndex.push_back(col);
                actualIndex.push_back(depth);
                ngldmnr = getNeighborGreyLevels(inputMatrix, actualIndex);
				actualIndex.clear();
				actualGreyIndex = ngldm.findIndex(this->diffGreyLevels, sizeGreyLevels, actualElement);
				NGLDMatrix[actualGreyIndex][ngldmnr] = NGLDMatrix[actualGreyIndex][ngldmnr] + 1;
            }  
        }
    }
    return NGLDMatrix;
}

/*!
\brief getNeighborGreyLevels
@param[in] boost multi array inputMatrix: matrix filled with intensity values
@param[in] vector actualIndex : vector of the actual index
@param[out] ngldmnr : number of elements in neighborhood

The function checks the grey levels in the neighborhood and counts how many of them have the same intensity value.
*/
template <class T, size_t R>
int NGLDMFeatures2DWOMerge<T, R>::getNeighborGreyLevels(boost::multi_array<T, R> inputMatrix, vector<int> actualIndex){
    int actualX = actualIndex[0];
    int actualY = actualIndex[1];
    int actualZ = actualIndex[2];
    int ngldmnr = 0;
    T tempElement;
    T actualElement = inputMatrix[actualX][actualY][actualZ];
    int leftBorderY = actualY - dist;
    if(leftBorderY<0){
        leftBorderY =0;
    }
    int leftBorderX;
    while(leftBorderY < 1+actualY + dist && leftBorderY < inputMatrix.shape()[1]){
        leftBorderX = actualX -dist;
        if(leftBorderX<0){
            leftBorderX=0;
        }
        while(leftBorderX < 1+actualX + dist && leftBorderX < inputMatrix.shape()[0]){
            tempElement = inputMatrix[leftBorderX][leftBorderY][actualZ];
            if(actualElement == tempElement ){
                ngldmnr +=1;
            }
            leftBorderX += 1;
        }
        leftBorderY += 1;
    }
    //because with this method I also counted the actual Element, I have to subtract 1
    ngldmnr -= 1;
    return ngldmnr;
}


template <class T, size_t R>
void NGLDMFeatures2DWOMerge<T, R>::calculateAllNGLDMFeatures2DWOMerge2D(NGLDMFeatures2DWOMerge<T,R> &ngldmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, ConfigFile config){
	ngldmFeatures.getConfigValues(config);
	//get config values
	coarseParam = config.coarsenessParam;
	dist = config.distNGLDM;
	this->diffGreyLevels = diffGrey;

    int totalDepth = inputMatrix.shape()[2];

    double meanRun;
	double meanGrey;
    double totalSum;

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
    T sumDependenceCountEnergy = 0;

    for(int depth = 0; depth < totalDepth; depth++){

        boost::multi_array<double,2> NGLDM=ngldmFeatures.getMatrix(inputMatrix, depth);
        totalSum = ngldmFeatures.calculateTotalSum(NGLDM);

        boost::multi_array<double,2> probMatrix = ngldmFeatures.calculateProbMatrix(NGLDM, totalSum);

        meanGrey = ngldmFeatures.calculateMeanProbGrey(probMatrix);
        meanRun = ngldmFeatures.calculateMeanProbRun(probMatrix);

        rowSums = ngldmFeatures.calculateRowSums(NGLDM);
        colSums = ngldmFeatures.calculateColSums(NGLDM);

        ngldmFeatures.calculateShortRunEmphasis(rowSums, totalSum);
        sumShortRunEmphasis += this->shortRunEmphasis;
        ngldmFeatures.calculateLongRunEmphasis(rowSums, totalSum);
        sumLongRunEmphasis += this->longRunEmphasis;
        ngldmFeatures.calculateLowGreyEmph(colSums, totalSum);
        sumLowGreyEmph += this->lowGreyEmph;
        ngldmFeatures.calculateHighGreyEmph(colSums, totalSum);
        sumHighGreyEmph += this->highGreyEmph;
        ngldmFeatures.calculateShortRunLow(NGLDM, totalSum);
        sumShortRunLow += this->shortRunLow;
        ngldmFeatures.calculateShortRunHigh(NGLDM, totalSum);
        sumShortRunHigh += this->shortRunHigh;
        ngldmFeatures.calculateLongRunLowEmph(NGLDM, totalSum);
        sumLongRunLowEmph += this->longRunLowEmph;
        ngldmFeatures.calculateLongRunHighEmph(NGLDM, totalSum);
        sumLongRunHighEmph += this->longRunHighEmph;
        ngldmFeatures.calculateGreyNonUniformity(colSums, totalSum);
        sumGreyNonUniformity += this->greyNonUniformity;
        ngldmFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
        sumGreyNonUniformityNorm += this->greyNonUniformityNorm;
        ngldmFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
        sumRunLengthNonUniformityNorm += this->runLengthNonUniformityNorm;
        ngldmFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
        sumRunLengthNonUniformity += this->runLengthNonUniformity;
        ngldmFeatures.calculateRunPercentage(inputMatrix, depth, totalSum, 1);
        sumRunPercentage += this->runPercentage;
        ngldmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);
        sumGreyLevelVar += this->greyLevelVar;
        ngldmFeatures.calculateRunLengthVar(probMatrix, meanRun);
        sumRunLengthVar += this->runLengthVar;
        ngldmFeatures.calculateRunEntropy(probMatrix);
        sumRunEntropy += this->runEntropy;
        ngldmFeatures.calculateDependenceCountEnergy(probMatrix);
        sumDependenceCountEnergy += this->dependenceCountEnergy;
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
    this->runLengthNonUniformity = sumRunLengthNonUniformity/totalDepth;
    this->runLengthNonUniformityNorm = sumRunLengthNonUniformityNorm/totalDepth;
    this->runPercentage = sumRunPercentage/totalDepth;
    this->greyLevelVar = sumGreyLevelVar/totalDepth;
    this->runLengthVar = sumRunLengthVar/totalDepth;
    this->runEntropy = sumRunEntropy/totalDepth;
    this->dependenceCountEnergy = sumDependenceCountEnergy/totalDepth;
}

template <class T, size_t R>
void NGLDMFeatures2DWOMerge<T, R>::extractNGLDMData2DWOMerge(vector<T> &ngldmData, NGLDMFeatures2DWOMerge<T, R> ngldmFeatures){

    ngldmData.push_back(ngldmFeatures.shortRunEmphasis);
    ngldmData.push_back(ngldmFeatures.longRunEmphasis);
    ngldmData.push_back(ngldmFeatures.lowGreyEmph);
    ngldmData.push_back(ngldmFeatures.highGreyEmph);
    ngldmData.push_back(ngldmFeatures.shortRunLow);
    ngldmData.push_back(ngldmFeatures.shortRunHigh);
    ngldmData.push_back(ngldmFeatures.longRunLowEmph);
    ngldmData.push_back(ngldmFeatures.longRunHighEmph);
    ngldmData.push_back(ngldmFeatures.greyNonUniformity);
    ngldmData.push_back(ngldmFeatures.greyNonUniformityNorm);
    ngldmData.push_back(ngldmFeatures.runLengthNonUniformity);   
    ngldmData.push_back(ngldmFeatures.runLengthNonUniformityNorm);
    ngldmData.push_back(ngldmFeatures.runPercentage);
    ngldmData.push_back(ngldmFeatures.greyLevelVar);
    ngldmData.push_back(ngldmFeatures.runLengthVar);
    ngldmData.push_back(ngldmFeatures.runEntropy);
    ngldmData.push_back(ngldmFeatures.dependenceCountEnergy);
}

template <class T, size_t R>
void NGLDMFeatures2DWOMerge<T, R>::writeCSVFileNGLDM(NGLDMFeatures2DWOMerge<T,R> ngldmFeat, string outputFolder)
{
    string csvName = outputFolder + "/ngldmFeatures2Davg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream ngldmCSV;
    ngldmCSV.open (name);
    vector<string> features;
    ngldm.defineNGLDMFeatures(features);

    vector<T> ngldmData;
    extractNGLDMData2DWOMerge(ngldmData, ngldmFeat);
    for(int i = 0; i< ngldmData.size(); i++){
        ngldmCSV << "ngldmFeatures2Davg"<< ","<<features[i] <<",";
        ngldmCSV << ngldmData[i];
        ngldmCSV << "\n";
    }
    ngldmCSV.close();
}

template <class T, size_t R>
void NGLDMFeatures2DWOMerge<T, R>::writeOneFileNGLDM(NGLDMFeatures2DWOMerge<T, R> ngldmFeat, string outputFolder) {
	string csvName = outputFolder + "/radiomicsFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream ngldmCSV;
	ngldmCSV.open(name, std::ios_base::app);
	vector<string> features;
	ngldm.defineNGLDMFeatures(features);

	vector<T> ngldmData;
	extractNGLDMData2DWOMerge(ngldmData, ngldmFeat);
	for (int i = 0; i< ngldmData.size(); i++) {
		ngldmCSV << "ngldmFeatures2Davg" << "," << features[i] << ",";
		ngldmCSV << ngldmData[i];
		ngldmCSV << "\n";
	}
	ngldmCSV.close();

}

#endif // NGLDMFEATURES2DWOMERGE_H_INCLUDED
