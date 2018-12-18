/*!
\brief getNeighborhoodMatrix3D
In this function matrices necessary to fill the NGTDM and NGLDM matrices are generated for the 3D case. \n
I do this in one function to save computational time.  \n

The matrix for the NGTDM matrix contains the sum of all neighborhoods, while the matrix for the NGLDM matrix
contains the number of elements smaller than the coarseness factor. \n
NGTDM and NGLDM matrix are given as reference to this function. \n

@param[in] Image imageAttr: image attribute of actual image
@param[in] boost multi array ngtdm3D: matrix to store results for ngtdm calculations, given as reference
@param[in] boost multi array ngldm3D: matrix to store results for ngldm calculations, given as reference
@param[in] vector spacing: spacing of actual image, necessary to calculate weights
@param[in] ConfigFile config: ConfigFile to read the distances the user set for the NGTDM and NGLDM case
*/
template<typename T>
void getNeighborhoodMatrix3D(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm3D, boost::multi_array<double, 2> &ngldm3D, vector<double> spacing, ConfigFile config) {
	int indexOfElement[3] = { 0,0,0 };
	NGLDMFeatures3D<T, 3> ngldm;
	//get distance set for NGTDM case
	int dist = config.dist;
	T actualElement;
	vector<T> matrixValues;
	int actualGreyIndex;
	int sizeGreyLevels = (imageAttr.diffGreyLevels).size();
	//for every matrix element we access the neighborhood
	for (int depth = 0; depth < imageAttr.imageMatrix.shape()[2]; depth++) {
		for (int row = 0; row < imageAttr.imageMatrix.shape()[0]; row++) {
			for (int col = 0; col < imageAttr.imageMatrix.shape()[1]; col++) {
				//store actual index in a vector
				indexOfElement[0] = row;
				indexOfElement[1] = col;
				indexOfElement[2] = depth;
				//get actual Element if it is the centre of a whole neighborhood
				actualElement = imageAttr.imageMatrix[row][col][depth];
				if (!std::isnan(actualElement)) {
					//get the sum of the actual neighborhood
					matrixValues = getNeighborhood3D(imageAttr.imageMatrix, indexOfElement, spacing, config);
					ngtdm3D[row][col][depth] = matrixValues[0];
					actualGreyIndex = ngldm.findIndex(imageAttr.diffGreyLevels, sizeGreyLevels, actualElement);
					//only if distance of NGLDM and NGTDM matrices are the same, we can combine the calculations
					//if not we have to calculate the NGLDM matrices separately what means more computational time
					if (dist == config.distNGLDM) {
						ngldm3D[actualGreyIndex][matrixValues[1]] = ngldm3D[actualGreyIndex][matrixValues[1]] + 1;
					
					}
					else {
						//matrixValues = getNeighborhoodNGLDM(imageAttr.imageMatrix, indexOfElement, spacing, config);
						//actualGreyIndex = ngldm.findIndex(imageAttr.diffGreyLevels, sizeGreyLevels, actualElement);
						//ngldm2D[actualGreyIndex][matrixValues[1]][depth] = ngldm2D[actualGreyIndex][matrixValues[1]][depth] + 1;

					}

				}
				else {
					ngtdm3D[row][col][depth] = 0;
				}
			}
		}
	}
}
/*!
\brief getNeighborhood3D
In this function we access the neighborhood of the actual voxel (given by the actual index).  \n

The sum of the actual neighborhood is calculated, as well as the number of elements independent from actual element are counted.

@param[in] boost multi array inputMatrix: actual matrix
@param[in] int indexOfElement: index of actual element
@param[in] vector spacing: spacing of actual image, necessary to calculate weights
@param[in] ConfigFile config: ConfigFile to read the distances the user set for the NGTDM and NGLDM case
*/
template<typename T>
vector<T> getNeighborhood3D(boost::multi_array<T, 3> inputMatrix, int *indexOfElement,  vector<double> spacing, ConfigFile config) {
	std::vector<T> neighborhood;
	T actElement;
	float weight;
	int count = 0;
	int dist = config.distNGLDM;
	//access whole neighborhood in all 3 dimensions
	for (int k = -dist; k<dist + 1; k++) {
		for (int i = -dist; i<dist + 1; i++) {
			for (int j = -dist; j<dist + 1; j++) {

				if (i != 0 || j != 0 || k != 0) {
					//the element has to lie inside the matrix, check border cases
					if (indexOfElement[0] + i>-1 && indexOfElement[0] + i<inputMatrix.shape()[0] && indexOfElement[1] + j>-1 && indexOfElement[1] + j<inputMatrix.shape()[1] && indexOfElement[2] + k>-1 && indexOfElement[2] + k<inputMatrix.shape()[2]) {
						actElement = inputMatrix[indexOfElement[0] + i][indexOfElement[1] + j][indexOfElement[2] + k];
						if (!std::isnan(actElement)) {
							weight = calculateWeight3D(i, j, k, config.normNGTDM, spacing);
							neighborhood.push_back(actElement);
						}
					}
				}
			}
		}
	}
	T nrNeighbor = std::count(neighborhood.begin(), neighborhood.end(), inputMatrix[indexOfElement[0]][indexOfElement[1]][indexOfElement[2]]);

	T sum = accumulate(neighborhood.begin(), neighborhood.end(), 0);
	if (neighborhood.size() > 0) {
		sum = sum / neighborhood.size();
	}
	//return a vector containing first the sum and the number of elements 
	vector<T> returnValue = { sum, nrNeighbor };
	return returnValue;

}
/*!
\brief getNeighborhoodMatrix2D
In this function matrices necessary to fill the NGTDM and NGLDM matrices are generated for the 2D case. \n
I do this in one function to save computational time.  \n

The matrix for the NGTDM matrix contains the sum of all neighborhoods, while the matrix for the NGLDM matrix
contains the number of elements smaller than the coarseness factor. \n
NGTDM and NGLDM matrix are given as reference to this function. \n

@param[in] Image imageAttr: image attribute of actual image
@param[in] boost multi array ngtdm2D: matrix to store results for ngtdm calculations, given as reference
@param[in] boost multi array ngldm2D: matrix to store results for ngldm calculations, given as reference
@param[in] vector spacing: spacing of actual image, necessary to calculate weights
@param[in] ConfigFile config: ConfigFile to read the distances the user set for the NGTDM and NGLDM case
*/
template<typename T>
void getNeighborhoodMatrix2D(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm2D, boost::multi_array<T, 3> &ngldm2D, vector<double> spacing, ConfigFile config) {
	int indexOfElement[3] = { 0,0,0 };
	NGLDMFeatures2DMRG<T, 3> ngldm;
	int dist = config.dist;
	T actualElement;
	vector<T> matrixValues;
	int actualGreyIndex;
	int sizeGreyLevels = (imageAttr.diffGreyLevels).size();
	//we are checking the neighborhoods for every voxel
	for (int depth = 0; depth < imageAttr.imageMatrix.shape()[2]; depth++) {
		for (int row = 0; row < imageAttr.imageMatrix.shape()[0]; row++) {
			for (int col = 0; col < imageAttr.imageMatrix.shape()[1]; col++) {
				indexOfElement[0] = row;
				indexOfElement[1] = col;
				indexOfElement[2] = depth;
				//get actual Element if it is the centre of a whole neighborhood
				actualElement = imageAttr.imageMatrix[row][col][depth];
				if (!std::isnan(actualElement)) {
					//get the sum of the actual neighborhood
					matrixValues = getNeighborhood(imageAttr.imageMatrix, indexOfElement, spacing, config);
					ngtdm2D[row][col][depth] = matrixValues[0];
					actualGreyIndex = ngldm.findIndex(imageAttr.diffGreyLevels, sizeGreyLevels, actualElement);
					//only if the distances are the same for NGLDM and NGTDM, we can combine the calculations
					//if not, we have to calcute the matrices separately
					if (dist == config.distNGLDM) {
						ngldm2D[actualGreyIndex][matrixValues[1]][depth] = ngldm2D[actualGreyIndex][matrixValues[1]][depth] + 1;
					}
					else {
						matrixValues = getNeighborhoodNGLDM(imageAttr.imageMatrix, indexOfElement, spacing, config);
						actualGreyIndex = ngldm.findIndex(imageAttr.diffGreyLevels, sizeGreyLevels, actualElement);
						ngldm2D[actualGreyIndex][matrixValues[1]][depth] = ngldm2D[actualGreyIndex][matrixValues[1]][depth] + 1;

					}

				}
				else {
					ngtdm2D[row][col][depth] = 0;
				}
			}
		}
	}

}


/*!
\brief getNeighborhood
In this function we access the neighborhood of the actual voxel (given by the actual index). It is basically the same function as
getNeighborhood3D but for the 2D case. \n

The sum of the actual neighborhood is calculated, as well as the number of elements independent from actual element are counted.

@param[in] boost multi array inputMatrix: actual matrix
@param[in] int indexOfElement: index of actual element
@param[in] vector spacing: spacing of actual image, necessary to calculate weights
@param[in] ConfigFile config: ConfigFile to read the distances the user set for the NGTDM and NGLDM case
*/
template<typename T>
vector<T> getNeighborhood(boost::multi_array<T, 3> inputMatrix, int *indexOfElement, vector<double> spacing, ConfigFile config) {
	
	vector<T> neighborhood;
	T actElement;
	float weight;
	int dist = config.dist;
	//for every element access the whole neighborhood given by dist
	for (int i = -dist; i<dist + 1; i++) {
		for (int j = -dist; j<dist + 1; j++) {
			if (i != 0 || j != 0) {
				//get sure that we do not cross the borders of the matrix
				if (indexOfElement[0] + i > -1 && indexOfElement[0] + i < inputMatrix.shape()[0] && indexOfElement[1] + j > -1 && indexOfElement[1] + j < inputMatrix.shape()[1]) {
					actElement = inputMatrix[indexOfElement[0] + i][indexOfElement[1] + j][indexOfElement[2]];
					if (!std::isnan(actElement)) {
						
						weight = calculateWeight2D(i, j, config.normNGTDM, spacing);
						neighborhood.push_back(weight * actElement);
					}
				}
			}
		}
	}
	T nrNeighbor = std::count(neighborhood.begin(), neighborhood.end(), inputMatrix[indexOfElement[0]][indexOfElement[1]][indexOfElement[2]]);
	T sum = accumulate(neighborhood.begin(), neighborhood.end(), 0);
	if (neighborhood.size() > 0) {
		sum = sum / neighborhood.size();
	}
	//return a vector containing first the sum and the number of elements 
	vector<T> returnValue = { sum, nrNeighbor };
	return returnValue;
}

/*!
\brief getNeighborhoodNGLDM
If the distances given for NGLDM and NGTDM are different, the values have to be calculated separately. 
This comes with a higher computational cost.   \n

The number of elements independent from actual element are counted.

@param[in] boost multi array inputMatrix: actual matrix
@param[in] int indexOfElement: index of actual element
@param[in] vector spacing: spacing of actual image, necessary to calculate weights
@param[in] ConfigFile config: ConfigFile to read the distances the user set for the NGTDM and NGLDM case
*/
template<typename T>
vector<T> getNeighborhoodNGLDM(boost::multi_array<T, 3> inputMatrix, int *indexOfElement, vector<double> spacing, ConfigFile config) {

	vector<T> neighborhood;
	T actElement;
	float weight;
	int dist = config.distNGLDM;
	for (int i = -dist; i<dist + 1; i++) {
		for (int j = -dist; j<dist + 1; j++) {
			if (i != 0 || j != 0) {
				if (indexOfElement[0] + i > -1 && indexOfElement[0] + i < inputMatrix.shape()[0] && indexOfElement[1] + j > -1 && indexOfElement[1] + j < inputMatrix.shape()[1]) {
					actElement = inputMatrix[indexOfElement[0] + i][indexOfElement[1] + j][indexOfElement[2]];
					if (!std::isnan(actElement)) {
						weight = calculateWeight2D(i, j, config.normNGTDM, spacing);
						neighborhood.push_back(weight * actElement);
					}
				}
			}
		}
	}

	T nrNeighbor = std::count(neighborhood.begin(), neighborhood.end(), inputMatrix[indexOfElement[0]][indexOfElement[1]][indexOfElement[2]]);
	T sum = accumulate(neighborhood.begin(), neighborhood.end(), 0);
	if (neighborhood.size() > 0) {
		sum = sum / neighborhood.size();
	}
	//return a vector containing first the sum and the number of elements 
	vector<T> returnValue = { sum, nrNeighbor };
	return returnValue;
}

/*!
\brief getNeighborhoodNGLDM3D
If the distances given for NGLDM and NGTDM are different, the values have to be calculated separately.
This comes with a higher computational cost. This function is the equivalent to the function getNeighborhoodNGLDM, but for 
the 3D case.\n

The number of elements independent from actual element are counted.

@param[in] boost multi array inputMatrix: actual matrix
@param[in] int indexOfElement: index of actual element
@param[in] vector spacing: spacing of actual image, necessary to calculate weights
@param[in] ConfigFile config: ConfigFile to read the distances the user set for the NGTDM and NGLDM case
*/
template<typename T>
vector<T> getNeighborhoodNGLDM3D(boost::multi_array<T, 3> inputMatrix, int *indexOfElement, vector<double> spacing, ConfigFile config) {

	vector<T> neighborhood;
	T actElement;
	float weight;
	int dist = config.distNGLDM;
	for (int k = -dist, k < dist; k++) {
		for (int i = -dist; i < dist + 1; i++) {
			for (int j = -dist; j < dist + 1; j++) {
				if (i != 0 || j != 0 ||k!=0) {
					if (indexOfElement[0] + i > -1 && indexOfElement[0] + i < inputMatrix.shape()[0] && indexOfElement[1] + j > -1 && indexOfElement[1] + j < inputMatrix.shape()[1] && indexOfElement[2] + k> -1 && indexOfElement[2] + k < inputMatrix.shape()[2]) {
						actElement = inputMatrix[indexOfElement[0] + i][indexOfElement[1] + j][indexOfElement[2]+k];
						if (!std::isnan(actElement)) {
							weight = calculateWeight2D(i, j, config.normNGLDM, spacing);
							neighborhood.push_back(weight * actElement);
						}
					}
				}
			}
		}
	}

	T nrNeighbor = std::count(neighborhood.begin(), neighborhood.end(), inputMatrix[indexOfElement[0]][indexOfElement[1]][indexOfElement[2]]);
	//T nrNeighbor = neighborhood.size();
	T sum = accumulate(neighborhood.begin(), neighborhood.end(), 0);
	if (neighborhood.size() > 0) {
		sum = sum / neighborhood.size();
	}
	//return a vector containing first the sum and the number of elements 
	vector<T> returnValue = { sum, nrNeighbor };
	return returnValue;
}
