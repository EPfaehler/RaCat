#ifndef GLRLMFEATURES_H_INCLUDED
#define GLRLMFEATURES_H_INCLUDED

#include "matrixFunctions.h"
#include "image.h"
#include "helpFunctions.h"
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

/*! \file */
/*!
The Grey Level Run Length Matrices also define a set of textural features. \n
The GLRL-matrices also investigate the distribution of the grey levels in the image. In this matrices the
run length of a grey level in a direction \f$ \Delta \f$ is counted. \n
Let \f$ N_g \f$ be the number of discretized grey levels present in the image. \n
Let \f$ N_r \f$ be the maximal run length present in the image. \n
Let \f$ M_{\Delta} \f$ be the \f$ N_g \times N_r \f$ GLRM matrix of direction \f$ \Delta \f$ \n
In this matrix every element represents how often an intensity element i occurred j-times consecutively. \n
The directions in which the run length are counted are the same as for the GLCM matrices.
The row number of the GLRLM matrix is representing the intensity value of the voxels and the
column number is representing the run length of these voxels. \n
E.g. the matrix element \f$ r_{ij} = r(i, j) \f$ is the number of occurences where discretized grey levels
occur j-times consecutively. \n
Let \f$ r_{i.} = \sum_{j = i} ^{N_r} r_{ij} \f$ be the marginal sum of the runs of the run lengths j for grey level i \n
Let \f$ r_{.j} = \sum_{i = 1} ^{N_g} r_{ij} \f$ be the marginal sum of the runs over the grey levels i for run length j \n

The feature values are calculated after calculating the GLRM matrices. \n
Also here are several methods possible to merge the matrices while calculating the features. The methods are the same as with the
GLCM matrices. \n

*/


template <class T,  size_t R=3>
class GLRLMFeatures  {
    private:
        typedef boost::multi_array<float,2> mat;
        float totalSum;
        void extractGLRLMData(vector<T> &glrlmData, GLRLMFeatures<T, R> glrlmFeatures);
        int totalNrVoxels;
        //store different grey levels in vector

    public:


        int maxRunLength;
        float powRow;
        float powCol;
		int calculateExtEmph;
        void defineGLRLMFeatures(vector<string> &features);
		void defineGLRLMFeaturesOntology(vector<string> &features);
        void getXYDirections(int &directionX, int &directionY, int angle);

        vector<T> diffGreyLevels;
        GLRLMFeatures(){

        }

        float shortRunEmphasis=NAN;
        float longRunEmphasis = NAN;
        float lowGreyEmph = NAN;
        float highGreyEmph = NAN;
        float shortRunLow = NAN;
        float shortRunHigh = NAN;
        float longRunLowEmph = NAN;
        float longRunHighEmph = NAN;
        float greyNonUniformity = NAN;
        float greyNonUniformityNorm = NAN;
        float runLengthNonUniformity = NAN;
        float runLengthNonUniformityNorm = NAN;
        float runPercentage = NAN;
        float greyLevelVar = NAN;
        float runLengthVar = NAN;
        float runEntropy = NAN;

        vector<float> calculateRowSums(boost::multi_array<float,2> glrlmatrix);
        vector<float> calculateColSums(boost::multi_array<float,2> glrlmatrix);
		int findIndex(vector<T> array, int size, T target);

		void getConfigValues(ConfigFile config);
		void setEmphasisValues(int extEmph, float powRow, float powCol);
        float calculateTotalSum(boost::multi_array<float,2> glrlMatrix);
        int getMaxRunLength(boost::multi_array<T, R> inputMatrix);
		
        void calculateShortRunEmphasis(vector<float> colSums, float totalSum);
        void calculateLongRunEmphasis(vector<float> colSums, float totalSum);
        void calculateLowGreyEmph(vector<float> colSums, float totalSum);
        void calculateHighGreyEmph(vector<float> colSums, float totalSum);
        void calculateShortRunLow(boost::multi_array<float,2> glrlmatrix, float totalSum);
        void calculateShortRunHigh(boost::multi_array<float,2> glrlmatrix, float totalSum);
        void calculateLongRunLowEmph(boost::multi_array<float,2> glrlmatrix, float totalSum);
        void calculateLongRunHighEmph(boost::multi_array<float,2> glrlmatrix, float totalSum);
        void calculateGreyNonUniformity(vector<float> colSums, float totalSum);
        void calculateGreyNonUniformityNorm(vector<float> colSums, float totalSum);
        void calculateRunLengthNonUniformityNorm(vector<float> rowSums, float totalSum);
        void calculateRunLengthNonUniformity(vector<float> rowSums, float totalSum);

        int calculateTotalNrVoxels(boost::multi_array<T,R> inputMatrix, int depth);
        void calculateTotalNrVoxels3D(vector<T> vectorMatrElement);
        void calculateRunPercentage(boost::multi_array<T,R> inputMatrix, int depth, float totalSum, int nrNeighbor);
        void calculateRunPercentage3D(vector<T> vectorMatrElement, float totalSum, int nrNeighbor);
        boost::multi_array<float,2> calculateProbMatrix(boost::multi_array<float,2> glrlmatrix, float totalSum);
        float calculateMeanProbGrey(boost::multi_array<float,2> probMatrix);
        void calculateGreyLevelVar(boost::multi_array<float,2> probMatrix, float mean);
        float calculateMeanProbRun(boost::multi_array<float,2> probMatrix);
        void calculateRunLengthVar(boost::multi_array<float,2> probMatrix, float meanRun);
        void calculateRunEntropy(boost::multi_array<float,2> probMatrix);
		
};
template <class T, size_t R>
int GLRLMFeatures<T, R>::findIndex(vector<T> array, int size, T target) {
	int i = 0;
	while ((i < size) && (array[i] != target)) i++;
	return (i < size) ? (i) : (-1);
}


template<class T, size_t R>
void GLRLMFeatures<T, R>::getConfigValues(ConfigFile config) {
	if (config.extendedEmphasis == 1) {
		setEmphasisValues(config.extendedEmphasis, config.powerRow, config.powerCol);
	}
	else if (config.extendedEmphasis == 0 ) {
		setEmphasisValues(0, float(1), float(1));
	}
	else {
		std::cout << "The value ExtendedEmphasisFeatures.CalculateExtendedEmph is not correct. Please us 1 or 0. It is automatically set to 0" << std::endl;
		setEmphasisValues(0, float(1), float(1));
	}
}
//set the exponential values to the user specified values
//this is part of the novel and uncommon feature section
template<class T, size_t R>
void GLRLMFeatures<T, R>::setEmphasisValues(int extEmph, float row, float col) {
	calculateExtEmph = extEmph;
	powRow = row;
	powCol = col;
}

//from the GLRL-Matrix calculate the probability matrx
//do this by dividing every matrix elemnt with the total nr. of voxels
template <class T, size_t R>
boost::multi_array<float,2> GLRLMFeatures<T, R>::calculateProbMatrix(boost::multi_array<float,2> glrlmatrix, float totalSum){
	boost::multi_array<float,2> probMatrix=glrlmatrix;
    transform( probMatrix.origin(), probMatrix.origin() + probMatrix.num_elements(),
                    probMatrix.origin(),  bind2nd(std::divides<float>(),int(totalSum)));
    return probMatrix;
}

/*!
\brief calculateMeanProbGrey
@param boost::multi_array<float,2> probMatrix matrix filled with the probabilities

calculates the mean probability of the appearance of every grey level
TODO change bordwers in for loop (probMatrix.shape())
*/
template <class T, size_t R>
float GLRLMFeatures<T, R>::calculateMeanProbGrey(boost::multi_array<float,2> probMatrix){
    float mean=0;
    for(int i=0; i<probMatrix.shape()[0]; i++){
        for(int j=0; j<probMatrix.shape()[1]; j++){
			if (!std::isnan(diffGreyLevels[i] * probMatrix[i][j])) {
				mean += diffGreyLevels[i] * probMatrix[i][j];
			}
			else {
				mean += 0;
			}
        }
    }
    return mean;
}

//calcuöate the mean probability of the runlength
template <class T, size_t R>
float GLRLMFeatures<T, R>::calculateMeanProbRun(boost::multi_array<float,2> probMatrix){
    float mean = 0;
    for(int i = 0; i < probMatrix.shape()[0]; i++){
        for(int j = 0; j < probMatrix.shape()[1]; j++){
			if (!std::isnan(j*probMatrix[i][j])) {
				mean += j*probMatrix[i][j];
			}
			else {
				mean += 0;
			}
        }
    }
    return mean;
}

/*!
\brief getMaxRunLength
@param boost::multi_array<T, R> inputMatrix
get the maximal run length \n
The maximal run length is the maximal size of one dimension
*/
template <class T, size_t R>
int GLRLMFeatures<T, R>::getMaxRunLength(boost::multi_array<T, R> inputMatrix){
    maxRunLength = std::max(inputMatrix.shape()[0], inputMatrix.shape()[1]);
    maxRunLength = std::max(maxRunLength, int(inputMatrix.shape()[2]));
    return maxRunLength;
}



/*!
calculate the sum of all matrix elements
*/
template<class T, size_t R>
float GLRLMFeatures<T, R>::calculateTotalSum(boost::multi_array<float,2> glrlmatrix){
    T sum = 0;
    sum = accumulate(glrlmatrix.origin(), glrlmatrix.origin() + glrlmatrix.num_elements(), 0 );
    return sum;
}

/*!
\brief calculateRowSums
@param boost::multi_array<float,2> glrlmatrix : GLRM matrix

calculates the sum of rows and stores them in the vector rowSums
*/
template<class T, size_t R>
vector<float> GLRLMFeatures<T,R>::calculateRowSums(boost::multi_array<float,2> glrlmatrix){
    vector<float> rowSums;
    rowSums.clear();
    int sum =0;
    for(int col = 0; col < glrlmatrix.shape()[1]; col++){
        sum = 0;
        for(int row = 0; row < glrlmatrix.shape()[0]; row++){
			sum += glrlmatrix[row][col];

        }
        rowSums.push_back(sum);

    }
    return rowSums;
}


/*!

\brief calculateColSums
@param boost::multi_array<float,2> glrlmatrix : GLRM matrix

calculates the sum of columns and stores them in the vector colSums
*/
template<class T, size_t R>
vector<float> GLRLMFeatures<T,R>::calculateColSums(boost::multi_array<float,2> glrlmatrix){
    int sum = 0;

    vector<float> colSums;
    colSums.clear();
    for(int row=0; row<glrlmatrix.shape()[0]; row++){
        sum =0;
        for(int col=0; col<glrlmatrix.shape()[1]; col++){
            sum += glrlmatrix[row][col];
            }
        colSums.push_back(sum);
    }
    return colSums;
}


/*!
\brief getXYDirections
@param int directionX
@param int directionY
@param int angle

The function gets directionX and directionY as reference. Depending on the angle value,
the parameter are set: \n
angle == 180 :  go one pixel/voxel in x-direction; no move in y-direction \n
angle == 90 :  no move in x-direction; go one pixel/voxel in y-direction \n
angle == 45 : go one pixel/voxel in x-direction; go one pixel/voxel in y direction \n
angle == 135 : go minus one pixel/voxel in x-direction; one pixel/voxel in y direction
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::getXYDirections(int &directionX, int &directionY, int angle){
    if(angle==0){
        directionX=0;
        directionY=0;

    }
        //if angle is 180 degrees, only look in x direction
    if(angle==180){
        directionX=1;
        directionY=0;
    }
    //if angle is 90 degrees only look in y direction
    else if(angle==90){
        directionX=0;
        directionY=1;
    }
    //if angle is in 45 degrees direction look in x and y direction
    else if(angle==45){
        directionX=1;
        directionY=1;
    }
    else if(angle==135){
        directionX=-1;
        directionY=1;
    }
    else{
        std::cout<<"Incorrect angle!"<<std::endl;
    }
}

/*!
\brief calculateTotalNrVoxels
@param boost::multi_array<float,2> glrlmatrix : GLRLM matrix

calculates the total number of voxels of one slice of the matrix

*/
template <class T, size_t R>
int GLRLMFeatures<T, R>::calculateTotalNrVoxels(boost::multi_array<T,R> inputMatrix, int depth){
    vector<T> vectorSliceElements;
    for(int row = 0; row < inputMatrix.shape()[0]; row++){
        for(int col = 0; col < inputMatrix.shape()[1]; col++){
            if(!std::isnan(inputMatrix[row][col][depth])){
                vectorSliceElements.push_back(inputMatrix[row][col][depth]);
            }
        }
    }
    totalNrVoxels=boost::size(vectorSliceElements);
	return totalNrVoxels;
}

template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateTotalNrVoxels3D(vector<T> vectorMatrElement){
    totalNrVoxels=boost::size(vectorMatrElement);
}

/*!
\brief calculateShortRunEmphasis
@param vector<float> rowSums : vector of the sums of the rows
@param float totalSum : sum of all matrix elements

This feature emphasizes the short runs. The higher the value, the more short runs are in the matrix.
*/

template<class T, size_t R>
void GLRLMFeatures<T, R>::calculateShortRunEmphasis(vector<float> rowSums, float totalSum){
    shortRunEmphasis = 0;
	if (totalSum != 0) {
		for(int j=0; j<rowSums.size(); j++){
			if (calculateExtEmph == 0) {
				if (!std::isnan(rowSums[j] / pow(j + 1, 2)) ){

					shortRunEmphasis += rowSums[j] / pow(j + 1, 2);
				}
				else {
					shortRunEmphasis += 0;
				}
			}
			else {
				if (!std::isnan(rowSums[j] / pow(j + 1, powRow))) {
					shortRunEmphasis += rowSums[j] / pow(j + 1, powRow);
				}
				else {
					shortRunEmphasis += 0;
				}
			}
		}

		shortRunEmphasis = shortRunEmphasis / totalSum;
	}
}


/*!
\brief calculateLongRunEmphasis
@param vector<float> rowSums : vector of the sums of the rows
@param float totalSum : sum of all matrix elements

This feature emphasizes the long runs. The higher the value, the more long runs are in the matrix.
*/

template<class T, size_t R>
void GLRLMFeatures<T, R>::calculateLongRunEmphasis(vector<float> rowSums, float totalSum){
    longRunEmphasis=0;
	if (totalSum != 0) {
		for(int j=0; j<rowSums.size(); j++){
			if (calculateExtEmph == 0) {
				if (!std::isnan(rowSums[j] * pow(j + 1, 2))) {
					longRunEmphasis += rowSums[j] * pow(j + 1, 2);
				}
				else {
					longRunEmphasis += 0;
				}
			}
			else {
				if (!std::isnan(rowSums[j] * pow(j + 1, powRow))) {
					longRunEmphasis += rowSums[j] * pow(j + 1, powRow);
				}
			}
		}
	
		longRunEmphasis = longRunEmphasis / totalSum;
	}
}


/*!
\brief calculateLowGreyEmph
@param vector<float> colSums : vector of the sums of the columns
@param float totalSum : sum of all matrix elements

This feature emphasizes the low grey levels. The higher the value, the more low grey levels are in the matrix.
*/
template<class T, size_t R>
void GLRLMFeatures<T, R>::calculateLowGreyEmph(vector<float> colSums, float totalSum){
    lowGreyEmph=0;
	if (totalSum != 0) {
		for(int i=0; i<colSums.size(); i++){
			if (calculateExtEmph == 0) {
				if (diffGreyLevels[i] != 0) {
					lowGreyEmph += colSums[i] / pow(diffGreyLevels[i], 2);
				}
			}
			else {
				if (diffGreyLevels[i] != 0) {
					lowGreyEmph += colSums[i] / pow(diffGreyLevels[i], powCol);
				}
			}
		}
		lowGreyEmph = lowGreyEmph / totalSum;
	}
}


/*!
\brief calculateHighGreyEmph
@param vector<float> colSums : vector of the sums of the columns
@param float totalSum : sum of all matrix elements

This feature emphasizes the high grey levels. The higher the value, the more high grey levels are in the matrix.
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateHighGreyEmph(vector<float> colSums, float totalSum){
    highGreyEmph=0;
	if (totalSum != 0) {
		for(int i=0; i<colSums.size(); i++){
			if (calculateExtEmph == 0) {
				if (!std::isnan(colSums[i] * pow(diffGreyLevels[i], 2))) {
					highGreyEmph += colSums[i] * pow(diffGreyLevels[i], 2);
				}

			}
			else {
				if (!std::isnan(colSums[i] * pow(diffGreyLevels[i], powCol))) {
					highGreyEmph += colSums[i] * pow(diffGreyLevels[i], powCol);
				}

			}
		}

		highGreyEmph = highGreyEmph / totalSum;
	}
}


/*!
\brief calculateShortRunLow
@param boost::multi_array<float,2> glrlmatrix : GLCM matrix
@param float totalSum : sum of all matrix elements

This feature emphasizes the low grey levels which habe a short run. The higher the value, the more low grey levels with short runs are in the matrix.
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateShortRunLow(boost::multi_array<float,2> glrlmatrix, float totalSum){
    shortRunLow = 0;
	if (totalSum != 0) {
		for(int row = 0; row < glrlmatrix.shape()[0]; row++){
			for(int col = 1; col < glrlmatrix.shape()[1]+1; col++){
				if (calculateExtEmph == 0) {
					if (!std::isnan(glrlmatrix[row][col - 1] / (pow(diffGreyLevels[row], 2)*pow(col, 2)))) {
						shortRunLow += glrlmatrix[row][col - 1] / (pow(diffGreyLevels[row], 2)*pow(col, 2));
					}
					else {
						shortRunLow += 0;
					}
				}
				else {
					if (!std::isnan(glrlmatrix[row][col - 1] / (pow(diffGreyLevels[row], 2)*pow(col, powCol)))) {
						shortRunLow += glrlmatrix[row][col - 1] / (pow(diffGreyLevels[row], powRow)*pow(col, powCol));
					}
				}
			}
		}

		shortRunLow = shortRunLow / totalSum;
	}
}


/*!
\brief calculateShortRunHigh
@param boost::multi_array<float,2> glrlmatrix : GLCM matrix
@param float totalSum : sum of all matrix elements

This feature emphasizes the high grey levels which habe a short run. The higher the value, the more high grey levels with short runs are in the matrix.
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateShortRunHigh(boost::multi_array<float,2> glrlmatrix, float totalSum){
    shortRunHigh = 0;
	if (totalSum != 0) {
		for(int row = 0; row < glrlmatrix.shape()[0]; row++){
			for(int col = 1; col < glrlmatrix.shape()[1]+1; col++){
				if (calculateExtEmph == 0) {
					if (!std::isnan(pow(diffGreyLevels[row], 2)*glrlmatrix[row][col - 1] / pow(col, 2))) {
						shortRunHigh += pow(diffGreyLevels[row], 2)*glrlmatrix[row][col - 1] / pow(col, 2);
					}
				}
				else {
					if (!std::isnan(pow(diffGreyLevels[row], 2)*glrlmatrix[row][col - 1] / pow(col, powCol))) {
						shortRunHigh += pow(diffGreyLevels[row], powRow)*glrlmatrix[row][col - 1] / pow(col, powCol);
					}
				
				}
			}
		}
	
		shortRunHigh = shortRunHigh / totalSum;
	}
}

/*!
\brief calculateLongRunLowEmph
@param boost::multi_array<float,2> glrlmatrix : GLCM matrix
@param float totalSum : sum of all matrix elements

This feature emphasizes the low grey levels which habe a long run. The higher the value, the more low grey levels with long runs are in the matrix.
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateLongRunLowEmph(boost::multi_array<float,2> glrlmatrix, float totalSum){
    longRunLowEmph = 0;
	if (totalSum != 0) {
		for(int row = 0; row < glrlmatrix.shape()[0]; row++){
			for(int col = 1; col < glrlmatrix.shape()[1]+1; col++){
				if(diffGreyLevels[row]!=0){
					if (calculateExtEmph == 0) {
						if (!std::isnan(pow(col, 2)*glrlmatrix[row][col - 1] / pow(diffGreyLevels[row], 2))) {
							longRunLowEmph += pow(col, 2)*glrlmatrix[row][col - 1] / pow(diffGreyLevels[row], 2);
						}
					}
					else {
						if (!std::isnan(pow(col, 2)*glrlmatrix[row][col - 1] / pow(diffGreyLevels[row], powRow))) {
							longRunLowEmph += pow(col, powCol)*glrlmatrix[row][col - 1] / pow(diffGreyLevels[row], powRow);
						}
					}
				}
			}
		}
	
		longRunLowEmph = longRunLowEmph / totalSum;
	}
}


/*!
\brief calculateLongRunHighEmph
@param boost::multi_array<float,2> glrlmatrix : GLCM matrix
@param float totalSum : sum of all matrix elements

This feature emphasizes the high grey levels which habe a long run. The higher the value, the more high grey levels with long runs are in the matrix.
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateLongRunHighEmph(boost::multi_array<float,2> glrlmatrix, float totalSum){
    longRunHighEmph=0;
	if (totalSum != 0) {
		for (int row = 0; row < glrlmatrix.shape()[0]; row++) {
			for (int col = 1; col < glrlmatrix.shape()[1] + 1; col++) {
				if (calculateExtEmph == 0) {
					if (!std::isnan(pow(col, 2)*pow(diffGreyLevels[row], 2)*glrlmatrix[row][col - 1])) {
						longRunHighEmph += pow(col, 2)*pow(diffGreyLevels[row], 2)*glrlmatrix[row][col - 1];
					}
				}
				else {
					if (!std::isnan(pow(col, 2)*pow(diffGreyLevels[row], powRow)*glrlmatrix[row][col - 1])) {
						longRunHighEmph += pow(col, powCol)*pow(diffGreyLevels[row], powRow)*glrlmatrix[row][col - 1];
					}
				}
			}
		}
		longRunHighEmph = longRunHighEmph / totalSum;
	}
}

/*!
\brief calculateGreyNonUniformity
@param vector<float> colSums : vector of the column sums
@param float totalSum : sum of all matrix elements

This features is a measure for the distribution of the grey levels in the image matrix. \n
The more equally distrbuted the runs of the grey levels are, the lower is the value.
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateGreyNonUniformity(vector<float> colSums, float totalSum){
    greyNonUniformity = 0;
    greyNonUniformity = for_each(colSums.begin(), colSums.end(), square_accumulate<float>()).result();
	if (totalSum != 0) {
		greyNonUniformity = greyNonUniformity / totalSum;
	}
	else {
		greyNonUniformity = 0;
	}
}


/*!
\brief calculateGreyNonUniformityNorm
@param vector<float> colSums : vector of the column sums
@param float totalSum : sum of all matrix elements

This features is a normalized version of the grey-non-uniformity feature.
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateGreyNonUniformityNorm(vector<float> colSums, float totalSum){
    greyNonUniformityNorm = for_each(colSums.begin(), colSums.end(), square_accumulate<float>()).result();
	if (totalSum != 0) {
		greyNonUniformityNorm = greyNonUniformityNorm / pow(totalSum, 2);
	}
	else {
		greyNonUniformityNorm = 0;
	}
		
}


/*!
\brief calculateRunLengthNonUniformity
@param vector<float> colSums : vector of the column sums
@param float totalSum : sum of all matrix elements

This feature is a measurement for the distribution of the run length. \n
The lower this value is, the more equally the run length are distributed.
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateRunLengthNonUniformity(vector<float> rowSums, float totalSum){
    runLengthNonUniformity=for_each(rowSums.begin(), rowSums.end(), square_accumulate<float>()).result();
	if (totalSum != 0) {
		runLengthNonUniformity = runLengthNonUniformity / totalSum;
	}
	else {
		runLengthNonUniformity = 0;
	}
}


/*!
\brief calculateRunLengthNonUniformityNorm
@param vector<float> colSums : vector of the column sums
@param float totalSum : sum of all matrix elements

This is a normalised version of the run-length non uniformity feature.
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateRunLengthNonUniformityNorm(vector<float> rowSums, float totalSum){
    runLengthNonUniformityNorm=for_each(rowSums.begin(), rowSums.end(), square_accumulate<float>()).result();
	if (totalSum != 0) {
		runLengthNonUniformityNorm = runLengthNonUniformityNorm / pow(totalSum, 2);
	}
	else {
		runLengthNonUniformityNorm = 0;
	}
}


/*!
\brief calculateRunPercentage
@param boost::multi_array<float,2> glrlmatrix : GLRLM matrix
@param float totalSum : sum of all matrix elements

calculates the fraction of runs appearing in the matrix and potential runs
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateRunPercentage(boost::multi_array<T,R> inputMatrix, int depth, float totalSum, int nrNeighbor){
	//if (depth == 0) {
	//	totalNrVoxels = inputMatrix.shape()[0] * inputMatrix.shape()[1] * inputMatrix.shape()[2];
	//}
	//else {
		totalNrVoxels = calculateTotalNrVoxels(inputMatrix, depth);
	//}
	if ((totalNrVoxels)*nrNeighbor != 0) {
		runPercentage = totalSum / ((totalNrVoxels)*nrNeighbor);
	}
	else {
		runPercentage = 0;
	}
}

template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateRunPercentage3D(vector<T> vectorMatrElement, float totalSum, int nrNeighbor){
    calculateTotalNrVoxels3D(vectorMatrElement);
	if ((totalNrVoxels)*nrNeighbor != 0) {
		runPercentage = totalSum / ((totalNrVoxels)*nrNeighbor);
	}
	else {
		runPercentage = 0;
	}
}

/*!
\brief calculateGreyLevelVar
@param boost::multi_array<float,2> probMatrix : probability matrix
@param float meanGrey : mean value of the grey levels

calculates the variance of grey levels \n
the lower the value, the more homogeneous is the region
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateGreyLevelVar(boost::multi_array<float,2> probMatrix, float meanGrey){
    greyLevelVar=0;
    for(int i=0; i<probMatrix.shape()[0]; i++){
        for(int j= 0; j<probMatrix.shape()[1]; j++){
			if (!std::isnan(pow((diffGreyLevels[i] - meanGrey), 2)*probMatrix[i][j])) {
				greyLevelVar += pow((diffGreyLevels[i] - meanGrey), 2)*probMatrix[i][j];
			}
        }
    }
}


/*!
\brief calculateRunLengthVar
@param boost::multi_array<float,2> probMatrix : probability matrix
@param float meanRun : mean value of the run length

calculates the variance of run length \n
the lower the value, the more homogeneous is the region
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateRunLengthVar(boost::multi_array<float,2> probMatrix, float meanRun){
    runLengthVar = 0;
    for(int i=0; i<probMatrix.shape()[0]; i++){
        for(int j= 0; j<probMatrix.shape()[1]; j++){
			if (!std::isnan(pow((j - meanRun), 2)*probMatrix[i][j])) {
				runLengthVar += pow((j - meanRun), 2)*probMatrix[i][j];
			}
        }
    }
}

/*!
\brief calculateRunEntropy
@param boost::multi_array<float,2> probMatrix : probability matrix

calculates the entropy of the probability matrix
*/
template <class T, size_t R>
void GLRLMFeatures<T, R>::calculateRunEntropy(boost::multi_array<float,2> probMatrix){
    runEntropy=0;
    for(int i=0; i<probMatrix.shape()[0]; i++){
        for(int j= 0; j<probMatrix.shape()[1]; j++){
            if(probMatrix[i][j]>0){
				if (!std::isnan(probMatrix[i][j] * log2(probMatrix[i][j]))) {
					runEntropy -= probMatrix[i][j] * log2(probMatrix[i][j]);
				}

            }
        }
    }
}


template <class T, size_t R>
void GLRLMFeatures<T, R>::extractGLRLMData(vector<T> &glrlmData, GLRLMFeatures<T, R> glrlmFeatures){

    glrlmData.push_back(glrlmFeatures.shortRunEmphasis);
    glrlmData.push_back(glrlmFeatures.longRunEmphasis);
    glrlmData.push_back(glrlmFeatures.lowGreyEmph);
    glrlmData.push_back(glrlmFeatures.highGreyEmph);
    glrlmData.push_back(glrlmFeatures.shortRunLow);
    glrlmData.push_back(glrlmFeatures.shortRunHigh);
    glrlmData.push_back(glrlmFeatures.longRunLowEmph);
    glrlmData.push_back(glrlmFeatures.longRunHighEmph);
    glrlmData.push_back(glrlmFeatures.greyNonUniformity);
    glrlmData.push_back(glrlmFeatures.greyNonUniformityNorm);
    glrlmData.push_back(glrlmFeatures.runLengthNonUniformity);   //check if I have to give the matrix
    glrlmData.push_back(glrlmFeatures.runLengthNonUniformityNorm);
    glrlmData.push_back(glrlmFeatures.runPercentage);
    glrlmData.push_back(glrlmFeatures.greyLevelVar);
    glrlmData.push_back(glrlmFeatures.runLengthVar);
    glrlmData.push_back(glrlmFeatures.runEntropy);

}



template <class T, size_t R>
void GLRLMFeatures<T, R>::defineGLRLMFeatures(vector<string> &features){
    features.push_back("short run emphasis");
    features.push_back("long runs emphasis");
    features.push_back("Low grey level run emphasis");
    features.push_back("High grey level run emphasis");
    features.push_back("Short run low grey level emphasis");
    features.push_back("Short run high grey level emphasis");
    features.push_back("Long run low grey level emphasis");
    features.push_back("Long run high grey level emphasis");
    features.push_back("Grey level non uniformity");
    features.push_back("Grey level non uniformity normalized");
    features.push_back("Run length non uniformity");
    features.push_back("Run length non uniformity normalized");
    features.push_back("Run percentage");
    features.push_back("Grey level variance");
    features.push_back("Run length variance");
    features.push_back("Run entropy");

}

template <class T, size_t R>
void GLRLMFeatures<T, R>::defineGLRLMFeaturesOntology(vector<string> &features) {
	features.push_back("Frlm.sre");
	features.push_back("Frlm.lre");
	features.push_back("Frlm.lgre");
	features.push_back("Frlm.hgre");
	features.push_back("Frlm.srlge");
	features.push_back("Frlm.srhge");
	features.push_back("Frlm.lrlge");
	features.push_back("Frlm.lrhge");
	features.push_back("Frlm.glnu");
	features.push_back("Frlm.glnu.norm");
	features.push_back("Frlm.rlnu");
	features.push_back("Frlm.rlnu.norm");
	features.push_back("Frlm.r.perc");
	features.push_back("Frlm.gl.var");
	features.push_back("Frlm.rl.var");
	features.push_back("Frlm.rl.entr");
	
}


#endif //GLRLMFEATURES_H_INCLUDED
