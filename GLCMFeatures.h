#ifndef GLCMFEATURES_H_INCLUDED
#define GLCMFEATURES_H_INCLUDED

#include <iostream>
#include "boost/multi_array.hpp"
#include "math.h"
#include <boost/range/algorithm.hpp>

#include "matrixFunctions.h"
#include "helpFunctions.h"
#include "image.h"

/*! \file */

using namespace std;

/*!


A Grey-Level-Co-occurence-Matrix (GLCM) is used to calculate the spacial dependence of grey levels in an image.
Before calculating the GLCM-matrix, the grey-values in the image have to be discretized.
The GLC-matrix expresses how combinations of the discretized grey-levels of neighbor pixels/voxels are distributed along one of the image directions.

Two different approaches are available: The 2D- and the 3D-approach \n

In 3D for every voxel we have 26 direct neighbors and 13 unique directions. \n

In 2D we only look at the neighbors slice by slice, ignoring the connectivity between slices. In the 2D approach, for every pixel, there are 8 neighbors
available.

In a GLC-matrix \f$ M_{\Delta} \f$ the number of rows and columns is equal to the number of gray levels \f$ N_g \f$ in the region of interest. \f$ \Delta \f$ is
the direction of the neighbor: Four different unique directions are possible: 0, 45, 90 and 135 degrees. For all of these directions a GLC-matrix is calculated. Every GLC-matrix is calculated as:
\f$ M_{\Delta} = M_{\delta} + M_{-\delta} \f$.

Every matrix element \f$ m_{ij} \f$ of \f$ M_{\Delta} \f$ is the relative frequency with which two pixels i and j occur within a given neighborhood ( i is the intensity of one element
and j is the intensity of its neighbor element).

After calculating the GLC-matrices, feature values can be calculated. Five different methods can be used to
generate the feature values. For all of these methods, a separate class is defined:

- GLCMFeaturesWOMerge \n
- GLCMFeaturesWMerge \n
- GLCMFeaturesFullMerge \n
- GLCMFeatures3DWMerge \n
- GLCMFeatures3DWOMerge

All this classes inherit from the class GLCMFeatures. In this class all feature calculations are defined.
*/

template <class T,  size_t R=3>
class GLCMFeatures{

    private:


        vector<T> sumProbRows;
        vector<T> sumProbCols;
        T HX;
        T HXY;
        T HXY1;
        T HXY2;
        int N_g;

    public:


		//!vectors where the diagonal and cross probabilities are stored
        vector<T> diagonalProbabilities;
        vector<T> crossProbabilities;
        //constructor
        GLCMFeatures(){
            //clear the vectors
            sumProbCols.clear();
            sumProbRows.clear();
            diagonalProbabilities.clear();
            crossProbabilities.clear();
        }
		~GLCMFeatures() {
		}

        //!define the different feature values of the GLCM matrix
        T jointMaximum;
        T jointAverage;
        T jointVariance;
        T jointEntropy;
        T diffAverage;
        T diffVariance;
        T diffEntropy;
        T sumAverage;
        T sumVariance;
        T sumEntropy;
        T angSecMoment;
        T contrast;
        T dissimilarity;
        T inverseDiff;
        T inverseDiffNorm;
        T inverseDiffMom;
        T inverseDiffMomNorm;
        T inverseVar;
        T meanRowProb;
        T meanColProb;
        T stdRowProb;
        T stdColProb;
        T autoCorrelation;
        T correlation;
        T clusterTendency;
        T clusterShade;
        T clusterProminence;
        T firstMCorrelation;
        T secondMCorrelation;


        void getXYDirections(int &directionX, int &directionY, int angle);
		std::vector<std::pair<T, T> > getNeighbours2D(boost::multi_array<T, R> inputMatrix, int depth, int directionX, int directionY);
		std::vector<std::pair<T, T> > getNeighbours3D(boost::multi_array<T, R> inputMatrix, int angle, int directionZ);
        void getDiagonalProbabilities(boost::multi_array<double, 2> &glcMatrix);
        void getCrossProbabilities(boost::multi_array<double, 2> &glcMatrix);
		void calculateMeanRowProb(boost::multi_array<double, 2> glcMatrix);
        void calculateMeanColProb(boost::multi_array<double, 2> glcMatrix);
        void calculateRowProb(boost::multi_array<double, 2> glcMatrix);
        void calculateColProb(boost::multi_array<double, 2> glcMatrix);


        void calculateJointMaximum(boost::multi_array<double, 2> glcMatrix);
        void calculateJointAverage(boost::multi_array<double, 2> glcMatrix);
        void calculateJointVariance(boost::multi_array<double, 2> glcMatrix, T jointAvg);
        void calculateJointEntropy(boost::multi_array<double, 2> glcMatrix);
        void calculateDiffAverage();
        void calculateDiffVariance(T diffAverage);
        void calculateDiffEntropy();
        void calculateSumAverage();
        void calculateSumVariance(T sumAverage);
        void calculateSumEntropy();
        void calculateAngSecMoment(boost::multi_array<double, 2> glcMatrix);
        void calculateContrast(boost::multi_array<double, 2> glcMatrix);
        void calculateDissimilarity(boost::multi_array<double, 2> glcMatrix);
        void calculateInverseDiff(boost::multi_array<double, 2> glcMatrix);
        void calculateInverseDiffNorm(boost::multi_array<double, 2> glcMatrix, T inverseDiff);
        void calculateInverseDiffMom(boost::multi_array<double, 2> glcMatrix);
        void calculateInverseDiffMomNorm(boost::multi_array<double, 2> glcMatrix);
        void calculateInverseVariance(boost::multi_array<double, 2> glcMatrix);
        void calculateCorrelation(boost::multi_array<double, 2> glcMatrix);


        void calculateAutoCorrelation(boost::multi_array<double, 2> glcMatrix);
        void calculateClusterTendency(boost::multi_array<double, 2> glcMatrix);
        void calculateClusterShade(boost::multi_array<double, 2> glcMatrix);
        void calculateClusterProminence(boost::multi_array<double, 2> glcMatrix);
        void calculateFirstMCorrelation(boost::multi_array<double, 2> glcMatrix);
        void calculateSecondMCorrelation(boost::multi_array<double, 2> glcMatrix);

        //store the feature values in a csv-file
        void defineGLCMFeatures(vector<string> &features);

};


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
void GLCMFeatures<T, R>::getXYDirections(int &directionX, int &directionY, int angle){
    //if we only go in the depth
    if(angle == 0){
        directionX = 0;
        directionY = 0;

    }
    if(angle == 180){
        directionX = 1;
        directionY = 0;
    }
    else if(angle == 90){
        directionX = 0;
        directionY = 1;
    }
    else if(angle == 45){
        directionX = 1;
        directionY = 1;
    }
    else if(angle == 135){
        directionX = -1;
        directionY = 1;
    }
}

/*!
The method getNeighbors2D stores all neighbor pairs for the desired angle and the actual input matrix in a vector
@param[in] inputMatrix: the original matrix of the VOI
@param[in] angle : the actual angle
@param[out] neighbors: vector containing all neighbor pairs of the actual direction

The method works as follows: \n
It looks for every matrix elements at the neighbors of the desired direction and makes a pair of the actual matrix element and
its neighbors. These pairs are stored in a vector and are returned.
*/
template <class T, size_t R>
std::vector<std::pair<T, T> > GLCMFeatures<T, R>::getNeighbours2D(boost::multi_array<T, R> inputMatrix, int depth, int directionX, int directionY) {
	//store the neighbours as pairs in a vector
	std::vector<std::pair<T, T> > neighbours;


	T neighbour1;
	T neighbour2;
	int maxRowNumber = inputMatrix.shape()[0];
	int maxColNumber = inputMatrix.shape()[1];

	for (int row = 0; row < maxRowNumber; row++) {
		//find now the neighbors, row by row
		for (int col = 0; col < maxColNumber - directionX*(directionX + 1) / 2; col++) {
			if (row - directionY > -1 && col + directionX > -1) {
				neighbour1 = inputMatrix[row][col][depth];
				neighbour2 = inputMatrix[row - directionY][col + directionX][depth];

				neighbours.push_back(std::make_pair(neighbour1, neighbour2));
			}
		}
	}

	return neighbours;
}

/*!
The method getNeighbors3D stores all neighbor pairs for the desired angle and the actual input matrix in a vector
@param[in] inputMatrix: the original matrix of the VOI
@param[in] angle : the actual angle
@param[in] directionZ: goes in the z-direction, adds the 3D calculation
@param[out] neighbors: vector containing all neighbor pairs of the actual direction
The method works as follows: \n
It looks for every matrix elements at the neighbors of the desired direction and makes a pair of the actual matrix element and
its neighbors. These pairs are stored in a vector and are returned.
*/
template <class T, size_t R>
std::vector<std::pair<T, T> > GLCMFeatures<T, R>::getNeighbours3D(boost::multi_array<T, R> inputMatrix, int angle, int directionZ) {
	//    store the neighbours as pairs in a vector
	std::vector<std::pair<T, T> > neighbours;
	int directionX;
	int directionY;
	//define in which direction you have to look for a neighbor
	GLCMFeatures<T, R> glcm;
	glcm.getXYDirections(directionX, directionY, angle);

	T neighbour1;
	T neighbour2;
	int maxRowNumber = inputMatrix.shape()[0];
	int maxColNumber = inputMatrix.shape()[1];
	int maxDepthNumber = inputMatrix.shape()[2];
	//if we have a 3D image go also in the depth
	for (int depth = 0; depth < maxDepthNumber - directionZ*(directionZ + 1) / 2; depth++) {
		for (int row = 0; row < maxRowNumber; row++) {
			//find now the neighbors, row by row
			for (int col = 0; col < maxColNumber - directionX*(directionX + 1) / 2; col++) {
				neighbour1 = inputMatrix[row][col][depth];
				if (row - directionY > -1 && col + directionX > -1 && depth + directionZ > -1) {
					neighbour2 = inputMatrix[row - directionY][col + directionX][depth + directionZ];
					neighbours.push_back(std::make_pair(neighbour1, neighbour2));
				}
			}
		}
	}
	return neighbours;
}


/*!
\brief calculateMeanRowProb
@param matrix glcMatrix: GLC-matrix

Calculate the mean row probability
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateMeanRowProb(boost::multi_array<double, 2> glcMatrix){
    calculateRowProb(glcMatrix);
    meanRowProb = 0;
	stdRowProb = 0;
    for(int k = 0; k < N_g; k++){
        meanRowProb += (k+1) * sumProbRows[k];
    }

    for(int k = 0; k < N_g; k++){
        stdRowProb += pow(((k+1)-meanRowProb),2)*sumProbRows[k];
    }
    stdRowProb = std::sqrt(stdRowProb);
}

/*!
\brief calculateMeanColProb
@param matrix glcMatrix: GLC-matrix

Calculate the mean column probability
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateMeanColProb(boost::multi_array<double, 2> glcMatrix){
    calculateColProb(glcMatrix);
    meanColProb = 0;
    for(int k = 0; k < N_g; k++){
        meanColProb += (k+1) * sumProbCols[k];

    }
    for(int k = 0; k < N_g; k++){
        stdColProb += pow(((k+1)-meanColProb),2)*sumProbCols[k];

    }
    stdColProb = std::sqrt(stdColProb);
}

/*!
\brief getDiagonalProbabilities
@param matrix glcMatrix: GLC-matrix

Calculate the diagonal probabilities and store them in vector diagonalProbabilities

formula see above
*/
template <class T, size_t R>
void GLCMFeatures<T, R>::getDiagonalProbabilities(boost::multi_array<double, 2> &glcMatrix){

    T diagProbability;
    diagonalProbabilities.clear();
    for(int k = 0; k < N_g; k++){
        diagProbability = 0;
        //get the diagonal elements
        for(int row = 0; row < N_g; row++){
            for(int col = 0; col < N_g; col++){
                //Kronecker delta
                if(k == abs(row-col)){
                    diagProbability += glcMatrix[row][col];
                }
            }
        }
        //store the diagonal probabilities in a vector
        diagonalProbabilities.push_back(diagProbability);
    }
}

/*!
\brief calculateRowProb
@param matrix glcMatrix: GLC-matrix

Calculate the sum of the row probabilities and store them in vector sumProbRows

formula see above
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateRowProb(boost::multi_array<double, 2> glcMatrix){
    T rowProb;
    sumProbRows.clear();
    for(int i = 0; i < N_g; i++ ){
        rowProb = 0;
        for(int j = 0; j < N_g; j++){
            rowProb += glcMatrix[i][j];
        }
        sumProbRows.push_back(rowProb);
    }
}

/*!
\brief calculateColProb
@param matrix glcMatrix: GLC-matrix

Calculate the sum of the col probabilities and store them in vector sumProbCOls

formula see above
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateColProb(boost::multi_array<double, 2> glcMatrix){
    T colProb;
    sumProbCols.clear();
    for(int j = 0; j < N_g; j++ ){
        colProb = 0;
        for(int i = 0; i < N_g; i++){
            colProb += glcMatrix[i][j];
        }
        sumProbCols.push_back(colProb);
    }
}

/*!
\brief getCrossProbabilities

@param matrix: GLC-matrix

Calculate the cross probabilities and store them in vector crossProbabilities
*/
template <class T, size_t R>
void GLCMFeatures<T, R>::getCrossProbabilities(boost::multi_array<double, 2> &glcMatrix){
    T crossProbability;
    crossProbabilities.clear();
    for(int k =0; k < 2*N_g; k++){
        crossProbability=0;
        for(int row = 0; row < N_g; row++){
            for(int col = 0; col < N_g; col++){
				//Kronecker delta
                if(k == abs(row + col)){
                    crossProbability += glcMatrix[row][col];
                }
            }
        }
        crossProbabilities.push_back(crossProbability);
    }
}

/*!
\brief calculateJointMaximum
@param matrix GLC-matrix

The joint maximum is the probability belonging to the neighbor pair which occurs the most in the VOI
*/
template <class T, size_t R>
void GLCMFeatures<T, R>::calculateJointMaximum(boost::multi_array<double, 2> glcMatrix){
    N_g = glcMatrix.shape()[0];
    jointMaximum = *max_element( glcMatrix.origin(), glcMatrix.origin() + glcMatrix.num_elements());
}

/*!
\brief calculateJointAverage
The joint average is the sum of joint probabilities. This sum is weighted by the grey level belonging to the
probabilities.
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateJointAverage(boost::multi_array<double, 2> glcMatrix){
    jointAverage = 0;
    for(int row = 0; row < N_g; row++){
        for(int col = 0; col < N_g; col++){
            jointAverage += (row + 1) * glcMatrix[row][col];
        }
    }
}

/*!
\brief calculateJointVariance
@param matrix GLC-matrix
@param jointAvg joint average

Calculates the joint variance using the joint average defined before \n
The joint variance is the variance of the numbers of neighbor pairs
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateJointVariance(boost::multi_array<double, 2> glcMatrix, T jointAvg){
    jointVariance = 0;
    for(int row = 0; row < N_g; row++){
        for(int col = 0; col < N_g; col++){
            jointVariance += pow((row+1-jointAvg),2) * glcMatrix[row][col];
        }
    }
}

/*!
\brief calculateJointEntropy
@param matrix GLC-matrix

Calculates the joint entropy \n
The joint entropy is a measurement for the uncertainity of numbers of neighbor pairs
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateJointEntropy(boost::multi_array<double, 2> glcMatrix){
    jointEntropy = 0;
    double actElement;
    for(int i = 0; i < N_g; i++){
        for(int j = 0; j < N_g; j++){
            actElement = glcMatrix[i][j];
            if(actElement != 0){
                jointEntropy -= actElement * log2(actElement);
            }
        }
    }
}

/*!
\brief calculateDiffAverage

The difference average is defined as the weighted average of the diagonal probabilities. \n
In this way the relationship between occurences of pairs with similar intensity values and differing intensity values is measured.

*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateDiffAverage(){
    diffAverage = 0;
    for(int diagElement = 0; diagElement < diagonalProbabilities.size(); diagElement++){
        diffAverage += diagonalProbabilities[diagElement] * diagElement;
    }
}


/*!
\brief calculateDiffVariance
@param diffAvg difference average
The difference variance is defined as the weighted variance of the diagonal probabilities. \n
It is a measure of heterogeneity. The more the intensity values are deviating from the mean, the higher is weight set to the intensity.
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateDiffVariance(T diffAvg){
    diffVariance = 0;
    for(int diagElement = 0; diagElement<diagonalProbabilities.size(); diagElement++){
        diffVariance += diagonalProbabilities[diagElement] * pow((diagElement - diffAvg),2);
    }
}

/*!
\brief calculateDiffEntropy

The difference entropy is defined as the entropy of the diagonal probabilities. \n
The different entropy feature gives information about the variablity in the differences of neighborhood intensities.
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateDiffEntropy(){
    diffEntropy = 0;
    double actElement = 0;
    for(int k = 0; k < diagonalProbabilities.size(); k++){
        actElement = diagonalProbabilities[k];
        if(actElement != 0){
            diffEntropy -= actElement*log2(actElement);
        }
    }
}


/*!
\brief calculateSumAverage

This method measures the relationship between occurences of low intensity pairs and high intensity pairs.
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateSumAverage(){
    sumAverage = 0;
    for(int crossElement = 1; crossElement < crossProbabilities.size()+1; crossElement++){
        sumAverage += crossProbabilities[crossElement-1]*(crossElement+1);
    }
}


/*!
\brief calculateSumVariance
@param sumAvg sum
The sum variance is defined as the weighted variance of the cross probabilities.  \n
Also this is a measurement of heterogeneity. The more the neighboring level pairs deviate from the mean, the higher they are weighted.
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateSumVariance(T sumAvg){
    sumVariance = 0;
    for(int crossElement = 1; crossElement < crossProbabilities.size()+1; crossElement++){
        sumVariance += crossProbabilities[crossElement-1] * pow((crossElement+1-sumAvg),2);
    }
}


/*!
\brief calculateSumEntropy
The sum entropy is a measurement of the differences in the intensity value pairs
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateSumEntropy(){
    sumEntropy = 0;
    double actElement = 0;
    for(int k = 0; k < crossProbabilities.size(); k++){
        actElement = crossProbabilities[k];
        if(actElement != 0 && !std::isnan(actElement)){
            sumEntropy -= actElement*log2(actElement);
        }
    }
}

/*!
\brief calculateAngSecMoment
@param matrix GLC-matrix

The angular second moment is the same as the energy of the probability distribution
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateAngSecMoment(boost::multi_array<double, 2> glcMatrix){
    angSecMoment = for_each(glcMatrix.origin(), glcMatrix.origin() + glcMatrix.num_elements(), square_accumulate<float>()).result();
}



/*!
\brief calculateContrast
@param matrix GLC-matrix

The contrast is a weighted sum of the elements of the GLC-matrix. The bigger the difference in the intensities of the neighbors
the higher the weight
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateContrast(boost::multi_array<double, 2> glcMatrix){
    contrast = 0;
    for(int i = 0; i < N_g; i++){
        for(int j = 0; j < N_g; j++){
            contrast += pow((i-j),2)*glcMatrix[i][j];
        }
    }
}


/*!
\brief calculateDissimilarity
@param matrix GLC-matrix

The calculation of the dissimilarity is similar to the calculation of the contrast \n
It represents the mean difference of the intensity values between neighboring pixels/voxels
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateDissimilarity(boost::multi_array<double, 2> glcMatrix){
    dissimilarity = 0;
    for(int i = 0; i < N_g; i++){
        for(int j = 0; j < N_g; j++){
            dissimilarity += abs(i-j)*glcMatrix[i][j];
        }
    }
}


/*!
\brief calculateInverseDiff
@param matrix GLC-matrix

The calculation of the inverse difference is the weighted sum of the GLC-matrix elements \n
The bigger the difference between the intensities of the neighboring pixels/voxels, the smaller the weight \n
The more equal the neighbor pairs are, the lower is the denominator and the higher is the value

*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateInverseDiff(boost::multi_array<double, 2> glcMatrix){
    inverseDiff = 0;
    for(int i = 0; i < N_g; i++){
        for(int j = 0; j < N_g; j++){
            inverseDiff +=glcMatrix[i][j]/(1+abs(i-j));
        }
    }
}


/*!
\brief calculateInverseDiffNorm
@param matrix GLC-matrix

The calculation of the inverse difference norm is similar to the inverse difference \n
The difference of the intensity differences is here normalised by the number of different grey levels \n
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateInverseDiffNorm(boost::multi_array<double, 2> glcMatrix, T inverseDiff){
    inverseDiffNorm = 0;
    for(int i = 0; i < N_g; i++){
        for(int j = 0; j < N_g; j++){
            inverseDiffNorm += glcMatrix[i][j]/(1+double(abs(i-j))/double(N_g));

        }
    }
}

/*!
\brief calculateInverseDiffMom
@param matrix GLC-matrix

The calculation of the inverse difference moment is similar to the inverse difference \n
It is also a weighted sum of the elements of the GLC-matrix. The higher the intensity values of a neighbor
pair differ, the smaller is the weight.
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateInverseDiffMom(boost::multi_array<double, 2> glcMatrix){
    inverseDiffMom = 0;
    for(int i = 0; i < N_g; i++){
        for(int j = 0; j < N_g; j++){
            inverseDiffMom += glcMatrix[i][j]/(1+pow((i-j),2));
        }
    }
}



/*!
\brief calculateInverseDiffMomNorm
@param matrix GLC-matrix

The calculation of the inverse difference moment norm is similar to the one of inverse difference moment \n
Here the difference of the intensities of neighbor pairs is normalised by the number of grey levels
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateInverseDiffMomNorm(boost::multi_array<double, 2> glcMatrix){
    inverseDiffMomNorm = 0;
    for(int i = 0; i < N_g; i++){
        for(int j = 0; j < N_g;j++){
            inverseDiffMomNorm += glcMatrix[i][j]/(1+pow((i-j),2)/pow(N_g,2));
        }
    }
}


/*!
\brief calculateInverseVariance
@param matrix GLC-matrix

The inverse variance is another measure if the image is locally homogen or not.
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateInverseVariance(boost::multi_array<double, 2> glcMatrix){
    inverseVar = 0;
    for(int i = 0; i < N_g;  i++){
        for(int j = i + 1; j < N_g; j++){
            inverseVar += 2*glcMatrix[i][j]/pow((i-j),2);
        }
    }
}

/*!
\brief calculateCorrelation
@param matrix GLC-matrix

The correlation shows the correlation between different grey values
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateCorrelation(boost::multi_array<double, 2> glcMatrix){
    calculateMeanColProb(glcMatrix);
    calculateMeanRowProb(glcMatrix);
    correlation = 0;
    for(int row = 1; row < N_g + 1; row++){
        for(int col = 1; col < N_g + 1; col++){
            correlation += (row - meanRowProb)*(col - meanRowProb)*glcMatrix[row-1][col-1];
        }
    }
	if (!isnan(correlation) && stdRowProb != 0 &&!isnan(stdRowProb)) {
		correlation = correlation / pow(stdRowProb, 2);
	}
}

/*!
\brief calculateAutoCorrelation
@param matrix GLC-matrix

AUto correlation measures the fineness or coarseness of the VOI.
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateAutoCorrelation(boost::multi_array<double, 2> glcMatrix){
    autoCorrelation = 0;
    for(int row = 1; row < N_g + 1; row++){
        for(int col = 1; col < N_g + 1; col++){
            autoCorrelation += glcMatrix[row-1][col-1]*(row)*(col);
        }
    }

}


/*!
\brief calculateClusterTendency
@param matrix GLC-matrix

The cluster tendency gives information about the formation of voxels with similar grey values in groups
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateClusterTendency(boost::multi_array<double, 2> glcMatrix){
    clusterTendency = 0;
    for(int i = 0; i < N_g; i++ ){
        for(int j = 0; j < N_g; j++){
            clusterTendency += pow((i+j+2-2*meanRowProb),2)*glcMatrix[i][j];
        }
    }
}

/*!
\brief calculateClusterShade
@param matrix GLC-matrix

The cluster shade measures the skewness and asymmetry of the VOI
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateClusterShade(boost::multi_array<double, 2> glcMatrix){
    clusterShade = 0;
    for(int i = 0; i < N_g; i++ ){
        for(int j = 0; j < N_g; j++){
            clusterShade += pow((i+j+2-2*meanRowProb),3)*glcMatrix[i][j];
        }
    }
}


/*!
\brief calculateClusterProminence
@param matrix GLC-matrix

The cluster prominence gives also information about the skewness and asymmetry of the VOI (higher value: more assymetry)
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateClusterProminence(boost::multi_array<double, 2> glcMatrix){
    clusterProminence = 0;
    for(int i = 0; i < N_g; i++ ){
        for(int j = 0; j < N_g; j++){
            clusterProminence += pow((i+j+2-2*meanRowProb),4)*glcMatrix[i][j];
        }
    }
}

/*!
\brief calculateFirstMCorrelation
@param matrix GLC-matrix

The first moment of correlation measures as well the homogeneity
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateFirstMCorrelation(boost::multi_array<double, 2> glcMatrix){
    HXY = 0;
    HX = 0;
    HXY1 = 0;
    double actualProbRows = 0;
    double actualProbCols = 0;
    for(int i = 0; i < N_g; i++ ){
        actualProbRows = sumProbRows[i];
        if(actualProbRows != 0){
            HX -= actualProbRows*log2(actualProbRows);
            for(int j = 0; j < N_g; j++){
                actualProbCols = sumProbCols[j];
                if(actualProbCols!= 0 && glcMatrix[i][j] != 0 &&!isnan(glcMatrix[i][j])&&actualProbRows!=0){
                    HXY -= glcMatrix[i][j]*log2(glcMatrix[i][j]);
                    HXY1 -= glcMatrix[i][j]*log2(actualProbRows*actualProbCols);
                }
            }
        }
    }
	if (HX > 0) {
		firstMCorrelation = (HXY - HXY1) / HX;
	}
}

/*!
\brief calculateSecondMCorrelation
@param matrix GLC-matrix

The second moment of correlation measures the similarity in intensity values for neighbor pairs
*/
template <class T, size_t R>
void GLCMFeatures<T,R>::calculateSecondMCorrelation(boost::multi_array<double, 2> glcMatrix){
    HXY2 = 0;
    for(int i = 0; i < N_g; i++ ){
      for(int j = 0; j < N_g; j++){
            if(sumProbRows[i]*sumProbCols[j] != 0){
                HXY2 -= sumProbRows[i]*sumProbCols[j]*log2(sumProbRows[i]*sumProbCols[j]);
            }
        }
    }
    secondMCorrelation = pow((1-exp(-2*(HXY2-HXY))),0.5);
}


template <class T, size_t R>
void GLCMFeatures<T, R>::defineGLCMFeatures(vector<string> &features){
    features.push_back("joint maximum");
    features.push_back("joint average");
    features.push_back("joint variance");
    features.push_back("joint entropy");
    features.push_back("difference average");
    features.push_back("difference variance");
    features.push_back("difference entropy");
    features.push_back("sum average");
    features.push_back("sum variance");
    features.push_back("sum entropy");
    features.push_back("angular second moment");
    features.push_back("contrast");
    features.push_back("dissimilarity");
    features.push_back("inverse difference");
    features.push_back("inverse difference normalised");
    features.push_back("inverse difference moment");
    features.push_back("inverse difference moment normalised");
    features.push_back("inverse variance");
    features.push_back("correlation");
    features.push_back("autocorrelation");
    features.push_back("cluster tendency");
    features.push_back("cluster shade");
    features.push_back("cluster prominence");
    features.push_back("first measure of information correlation");
    features.push_back("second measure of information correlation");
}



#endif // GLCMFEATURES_H_INCLUDED
