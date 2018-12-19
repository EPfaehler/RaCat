#include "matrixFunctions.h"

void matrixSum(boost::multi_array<double, 2> &matrix1, boost::multi_array<double, 2> matrix2){
    for(int i = 0; i<matrix1.shape()[0]; i++){
        for(int j = 0; j<matrix1.shape()[1]; j++){
                matrix1[i][j]+=matrix2[i][j];
        }
    }
}

void inverse(boost::multi_array<double, 2> matrix, boost::multi_array<double, 2> &inverseMatrix){
    for(int i = 0; i<matrix.shape()[0]; i++){
        for(int j = 0; j<matrix.shape()[1]; j++){
            inverseMatrix[j][i]=matrix[i][j];
        }
    }
}

void multSkalarMatrix(boost::multi_array<double, 2> &matrix, double weight){
    for(int i = 0; i<matrix.shape()[0]; i++){
        for(int j = 0; j<matrix.shape()[1]; j++){
            matrix[i][j]=weight*matrix[i][j];
        }
    }
}

float calculateWeight2D(int directionX, int directionY, string norm, vector<double> spacing){
	float weight;
    if(!norm.compare("Manhattan") || !norm.compare("manhattan")){
        weight = calculateManhattanNorm2D(directionX, directionY, spacing);
     }
     else if(!norm.compare("Chebyshev") || !norm.compare("chebyshev")){
        weight = 1;
     }
     else if(!norm.compare("Euclidean") || !norm.compare("euclidean")){
        weight = calculateEuclidianNorm2D(directionX, directionY, spacing);
     }
     else{
        std::cout<<"The chosen norm is not implemented, please check the spelling. The Chebyshev norm is chosen."<<std::endl;
        weight = 1;
     }
     return weight;
}

float calculateWeight3D(int directionX, int directionY, int directionZ, string norm, vector<double> spacing) {
	float weight;
	if (!norm.compare("Manhattan") || !norm.compare("manhattan")) {
		weight = calculateManhattanNorm3D(directionX, directionY, directionZ, spacing);
	}
	else if (!norm.compare("Chebyshev") || !norm.compare("chebyshev")) {
		weight = 1;
	}
	else if (!norm.compare("Euclidean") || !norm.compare("euclidean")) {
		weight = calculateEuclidianNorm3D(directionX, directionY, directionZ, spacing);
	}
	else {
		std::cout << "The chosen norm is not implemented, please check the spelling. The Chebyshev norm is chosen." << std::endl;
		weight = 1;
	}
	return weight;
}