#ifndef NGLDMFEATURES3D_H_INCLUDED
#define NGLDMFEATURES3D_H_INCLUDED

#include "NGLDMFeatures2DMRG.h"

/*! \file */
/*!
The class NGLDMFeatures3D inherites from the class NGLDMFeatures. The feature calculation are the same, only this matrix checks 3D neighborhoods. \n
A neighborhood are all voxels around one voxel within distance dist. \n
A voxel is dependent from the other, if \f$ |X_{c}-X_{m}|<a\f$, where \f$X_{c}\f$ is the center voxel and \f$X_{m}\f$ are the
other voxels in the neighborhood. The number of dependent voxels in a neighborhood is counted. \n
\f$ a \f$ is called the coarseness parameter. \n
The neighborhoods are checked for every voxel. \n
\f$ s_ = s(i,j) \f$ is the number of neighborhoods with center voxel with grey level i and dependece k = j-1 \n

The most features are already defined in the GLRLM features.
*/

template <class T,  size_t R>
class NGLDMFeatures3D : public NGLDMFeatures2DMRG<T,R>{
    private:

		int dist;
		int coarseParam;

		vector<T> sumProbRows;
		vector<T> sumProbCols;
		vector<double> rowSums;
		vector<double> colSums;
		NGLDMFeatures2DMRG<T, R> ngldm;
		void defineNGLDMFeatures3D(vector<string> &features);
        void extractNGLDMData3D(vector<T> &ngldmData, NGLDMFeatures3D<T, R> ngldmFeatures);

    public:
        void writeCSVFileNGLDM3D(NGLDMFeatures3D<T,R> ngldmFeat, string outputFolder);
		void writeOneFileNGLDM3D(NGLDMFeatures3D<T, R> ngldmFeat, string outputFolder);
        void calculateAllNGLDMFeatures3D(NGLDMFeatures3D<T,R> &ngldmFeatures, boost::multi_array<double, 2> ngldm3DMatrix, Image<T, R> imageAttr, ConfigFile config);
};



/*!
\brief calculateAllNGLDMFeatures3D
@param[in] NGLDMFeatures3D ngldmFeatures: element of class ngldmFeatures
@param[in] multi_array NGLDM: NGLDM matrix calculated before
@param[out] Image imageAttr : element image containing all information about the current image

This function calculates all NGLDM features.
*/
template <class T, size_t R>
void NGLDMFeatures3D<T, R>::calculateAllNGLDMFeatures3D(NGLDMFeatures3D<T,R> &ngldmFeatures, boost::multi_array<double, 2> NGLDM, Image<T,R> imageAttr, ConfigFile config){
	ngldmFeatures.getConfigValues(config);

	coarseParam = config.coarsenessParam;
	dist = config.distNGLDM;

	this->diffGreyLevels = imageAttr.diffGreyLevels;
	if (boost::size(diffGreyLevels) > 1) {
		double totalSum = ngldmFeatures.calculateTotalSum(NGLDM);
		rowSums = ngldmFeatures.calculateRowSums(NGLDM);
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
		ngldmFeatures.calculateRunPercentage3D(imageAttr.vectorOfMatrixElements, totalSum, 1);
		boost::multi_array<double, 2> probMatrix = ngldmFeatures.calculateProbMatrix(NGLDM, totalSum);
		double meanGrey = ngldmFeatures.calculateMeanProbGrey(probMatrix);

		ngldmFeatures.calculateGreyLevelVar(probMatrix, meanGrey);

		double meanRun = ngldmFeatures.calculateMeanProbRun(probMatrix);
		ngldmFeatures.calculateRunLengthVar(probMatrix, meanRun);
		ngldmFeatures.calculateRunEntropy(probMatrix);
		ngldmFeatures.calculateDependenceCountEnergy(probMatrix);
	}
	else {
		std::cout << "There is only one grey level in the image, the NGLDM 3D features cannot be calculated" << std::endl;
	}
}



template <class T, size_t R>
void NGLDMFeatures3D<T, R>::writeCSVFileNGLDM3D(NGLDMFeatures3D<T,R> ngldmFeat, string outputFolder)
{
    string csvName = outputFolder + "_ngldmFeatures3D.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream ngldmCSV;
    ngldmCSV.open (name);
    vector<string> features;
    defineNGLDMFeatures3D(features);

    vector<T> ngldmData;
    extractNGLDMData3D(ngldmData, ngldmFeat);
    for(int i = 0; i< ngldmData.size(); i++){
        ngldmCSV <<"ngldmFeatures3D"<<","<< features[i] <<",";
        ngldmCSV << ngldmData[i];
        ngldmCSV << "\n";
    }
    ngldmCSV.close();
}

template <class T, size_t R>
void NGLDMFeatures3D<T, R>::writeOneFileNGLDM3D(NGLDMFeatures3D<T, R> ngldmFeat, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream ngldmCSV;
	ngldmCSV.open(name, std::ios_base::app);
	vector<string> features;
	defineNGLDMFeatures3D(features);

	vector<T> ngldmData;
	extractNGLDMData3D(ngldmData, ngldmFeat);
	for (int i = 0; i< ngldmData.size(); i++) {
		ngldmCSV << "ngldmFeatures3D" <<","<< features[i] << ",";
		ngldmCSV << ngldmData[i];
		ngldmCSV << "\n";
	}
	ngldmCSV.close();
}

template <class T, size_t R>
void NGLDMFeatures3D<T, R>::extractNGLDMData3D(vector<T> &ngldmData, NGLDMFeatures3D<T, R> ngldmFeatures){

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
void NGLDMFeatures3D<T, R>::defineNGLDMFeatures3D(vector<string> &features){
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
#endif // NGLDMFEATURES3D_H_INCLUDED
