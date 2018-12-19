#ifndef NGLDMFEATURES_H_INCLUDED
#define NGLDMFEATURES_H_INCLUDED
#include "GLRLMFeatures.h"


/*! \file */
/*!
The class NGLDM is the class of the Neighborhood Grey Level Dependence Matrices. \n
A neighborhood are all voxels around one voxel within distance dist. \n
A voxel is dependent from the other, if \f$ |X_{c}-X_{m}|<a\f$, where \f$X_{c}\f$ is the center voxel and \f$X_{m}\f$ are the 
other voxels in the neighborhood. The number of dependent voxels in a neighborhood is counted. \n
\f$ a \f$ is called the coarseness parameter. \n
The neighborhoods are checked for every voxel. \n
\f$ s_ = s(i,j) \f$ is the number of neighborhoods with center voxel with grey level i and dependece k = j-1 \n

The most features are already defined in the GLRLM features.
*/
template <class T,  size_t R>
class NGLDMFeatures : public GLRLMFeatures<T,R>{
    private:
		vector<T> sumProbRows;
		vector<T> sumProbCols;
		vector<double> rowSums;
		vector<double> colSums;

		int dist;
		int coarseParam;

        void extractNGLDMData(vector<T> &ngldmData, NGLDMFeatures<T, R> ngldmFeatures);
        boost::multi_array<double, 2> getMatrix(boost::multi_array<T,R> inputMatrix);
        int getNeighborGreyLevels(boost::multi_array<T, R> inputMatrix, vector<int> actualIndex);
		

    public:
        double dependenceCountEnergy;
		int findIndex(vector<T> array, int size, T target);
        void writeCSVFileNGLDM(NGLDMFeatures<T,R> ngldmFeat, string outputFolder);
		void writeOneFileNGLDM(NGLDMFeatures<T, R> ngldmFeat, string outputFolder);
        void calculateAllNGLDMFeatures2D(NGLDMFeatures<T,R> &ngldmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElem, ConfigFile config);
        void calculateDependenceCountEnergy(boost::multi_array<double,2> probMatrix);
        void defineNGLDMFeatures(vector<string> &features);

};


template <class T, size_t R>
int NGLDMFeatures<T, R>::findIndex(vector<T> array, int size, T target) {
	int i = 0;
	while ((i < size) && (array[i] != target)) i++;
	return (i < size) ? (i) : (-1);
}

/*!
\brief getMatrix
@param boost multi array inputMatrix: matrix filled with intensity values
@param[out] boost multi array: filled NGLD matrix

The function fills the NGLDMatrix with the corresponding values using the function getNeighborGreyLevels. It 
checks voxel by voxel the neighborhood in the distance that is set by the user.
*/
template <class T, size_t R>
boost::multi_array<double, 2> NGLDMFeatures<T, R>::getMatrix(boost::multi_array<T,R> inputMatrix){
    typedef boost::multi_array<double, 2>  ngldmat;
    vector<int> actualIndex;
    T actualElement;
	//get number of different grey levels in the VOI
    int sizeGreyLevels = (this->diffGreyLevels).size();
    T actualGreyLevel;
	//here the position of the grey level in the different grey level-array will be stored
	int actualGreyIndex;
    int ngldmnr;
    //define the NGLDMarices; col.size=9 because we have 8 neighbors
    ngldmat NGLDMatrix(boost::extents[sizeGreyLevels][9]);
	//check every element of the VOI
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		for (int row = 0; row < inputMatrix.shape()[0]; row++) {
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				actualElement = inputMatrix[row][col][depth];
				if (!std::isnan(actualElement)) {
					actualGreyIndex = findIndex(this->diffGreyLevels, boost::size(this->diffGreyLevels), actualElement);
					actualIndex.push_back(row);
					actualIndex.push_back(col);
					actualIndex.push_back(depth);
					ngldmnr = getNeighborGreyLevels(inputMatrix, actualIndex);
					actualIndex.clear();
					//update NGLDM
					NGLDMatrix[actualGreyIndex][ngldmnr] = NGLDMatrix[actualGreyIndex][ngldmnr] + 1;
				}
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
int NGLDMFeatures<T, R>::getNeighborGreyLevels(boost::multi_array<T, R> inputMatrix, vector<int> actualIndex){
    int actualX = actualIndex[0];
    int actualY = actualIndex[1];
    int actualZ = actualIndex[2];
    int ngldmnr = 0;
    T tempElement;
	//get the actual grey value
    T actualElement = inputMatrix[actualX][actualY][actualZ];
	//leftBorder
    int leftBorderY = actualY - dist;
    if(leftBorderY<0){
        leftBorderY =0;
    }
	//go from the left to the right border
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


/*!
\brief calculateDependenceCountEnergy
@param[in] boost multi array probMatrix: probability matrix filled

The function calculate the dependence count energy: \f$ F_{countEnergy} = \sum_{i=1}^{N_{g} \sum_{j=1}^{N_{g} p_{ij}^{2} \f$.
*/
template <class T, size_t R>
void NGLDMFeatures<T, R>::calculateDependenceCountEnergy(boost::multi_array<double,2> probMatrix){
    dependenceCountEnergy=0;
    for(int row=0; row<probMatrix.shape()[0]; row++){
        for(int col=0; col<probMatrix.shape()[1]; col++){
            dependenceCountEnergy += pow(probMatrix[row][col], 2);
        }

    }
}


template <class T, size_t R>
void NGLDMFeatures<T, R>::calculateAllNGLDMFeatures2D(NGLDMFeatures<T,R> &ngldmFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGrey, vector<T> vectorMatrElem, ConfigFile config){
    this->diffGreyLevels = diffGrey;
	ngldmFeatures.getConfigValues(config);
	//get values from config file
	coarseParam = config.coarsenessParam;
	this->dist = config.distNGLDM;
	//get NGLDM 
    boost::multi_array<double,2> NGLDM=ngldmFeatures.getMatrix(inputMatrix);
    double totalSum = ngldmFeatures.calculateTotalSum(NGLDM);

    rowSums=ngldmFeatures.calculateRowSums(NGLDM);
    colSums = ngldmFeatures.calculateColSums(NGLDM);
    ngldmFeatures.calculateShortRunEmphasis(rowSums, totalSum);
    ngldmFeatures.calculateLongRunEmphasis(rowSums, totalSum);
    ngldmFeatures.calculateLowGreyEmph(colSums, totalSum);
    ngldmFeatures.calculateHighGreyEmph(colSums, totalSum);
    ngldmFeatures.calculateShortRunLow(NGLDM, totalSum);
    ngldmFeatures.calculateShortRunHigh(NGLDM, totalSum);
    ngldmFeatures.calculateLongRunLowEmph(NGLDM, totalSum);
    ngldmFeatures.calculateLongRunHighEmph(NGLDM, totalSum);
    ngldmFeatures.calculateGreyNonUniformity(colSums, totalSum);
    ngldmFeatures.calculateGreyNonUniformityNorm(colSums, totalSum);
    ngldmFeatures.calculateRunLengthNonUniformityNorm(rowSums, totalSum);
    ngldmFeatures.calculateRunLengthNonUniformity(rowSums, totalSum);
    ngldmFeatures.calculateRunPercentage3D(vectorMatrElem, totalSum, 1);
    boost::multi_array<double,2> probMatrix = ngldmFeatures.calculateProbMatrix(NGLDM, totalSum);
    double meanGrey = ngldmFeatures.calculateMeanProbGrey(probMatrix);
    ngldmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);

    double meanRun = ngldmFeatures.calculateMeanProbRun(probMatrix);
    ngldmFeatures.calculateRunLengthVar(probMatrix, meanRun);
    ngldmFeatures.calculateRunEntropy(probMatrix);
    ngldmFeatures.calculateDependenceCountEnergy(probMatrix);
}

template <class T, size_t R>
void NGLDMFeatures<T, R>::extractNGLDMData(vector<T> &ngldmData, NGLDMFeatures<T, R> ngldmFeatures){

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
    ngldmData.push_back(ngldmFeatures.runLengthNonUniformity);   //check if I have to give the matrix
    ngldmData.push_back(ngldmFeatures.runLengthNonUniformityNorm);
    ngldmData.push_back(ngldmFeatures.runPercentage);
    ngldmData.push_back(ngldmFeatures.greyLevelVar);
    ngldmData.push_back(ngldmFeatures.runLengthVar);
    ngldmData.push_back(ngldmFeatures.runEntropy);
    ngldmData.push_back(ngldmFeatures.dependenceCountEnergy);

}

template <class T, size_t R>
void NGLDMFeatures<T, R>::writeCSVFileNGLDM(NGLDMFeatures<T,R> ngldmFeat, string outputFolder)
{
    string csvName = outputFolder + "/ngldmFeatures2Dmrg.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream ngldmCSV;
    ngldmCSV.open (name);
    vector<string> features;
    defineNGLDMFeatures(features);

    vector<T> ngldmData;
    extractNGLDMData(ngldmData, ngldmFeat);
    for(int i = 0; i< ngldmData.size(); i++){
        ngldmCSV << "ngldmFeatures2Dmrg"<<","<<features[i] <<",";
        ngldmCSV << ngldmData[i];
        ngldmCSV << "\n";
    }
    ngldmCSV.close();
}

template <class T, size_t R>
void NGLDMFeatures<T, R>::writeOneFileNGLDM(NGLDMFeatures<T, R> ngldmFeat, string outputFolder) {
	string csvName = outputFolder + "/radiomicsFeatures.csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream ngldmCSV;
	ngldmCSV.open(name, std::ios_base::app);
	vector<string> features;
	defineNGLDMFeatures(features);

	vector<T> ngldmData;
	extractNGLDMData(ngldmData, ngldmFeat);
	for (int i = 0; i< ngldmData.size(); i++) {
		ngldmCSV << "ngldmFeatures2Dmrg" << "," << features[i] << ",";
		ngldmCSV << ngldmData[i];
		ngldmCSV << "\n";
	}
	ngldmCSV.close();
}

template <class T, size_t R>
void NGLDMFeatures<T, R>::defineNGLDMFeatures(vector<string> &features){
    features.push_back("Low dependence emphasis");
    features.push_back("High dependence emphasis");
    features.push_back("Low grey level count emphasis");
    features.push_back("High grey level count emphasis");
    features.push_back("Low dependence low grey level emphasis");
    features.push_back("Low dependence high grey level emphasis");
    features.push_back("High dependence low grey level emphasis");
    features.push_back("High dependence high grey level emphasis");
    features.push_back("Grey level non uniformity");
    features.push_back("Grey level non uniformity normalized");
    features.push_back("Dependence count non uniformity");
    features.push_back("Dependence count non uniformity normalized");
    features.push_back("Dependence count percentage");
    features.push_back("Grey level variance");
    features.push_back("Dependence count variance");
    features.push_back("Dependence count entropy");
    features.push_back("dependence Count Energy");

}
#endif // NGLDMFEATURES_H_INCLUDED
