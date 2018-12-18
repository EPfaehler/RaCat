//void subtractScalarFromMatrix(boost::multi_array<float, 3> &inputMatrix, float scalar){
//    const long unsigned int* matrixShape = inputMatrix.shape();
//    for(int i=0; i<matrixShape[0]; i++)
//        for(int j=0; j<matrixShape[1]; j++)
//            for(int k=0; k<matrixShape[2]; k++){
//                inputMatrix[i][j][k]=inputMatrix[i][j][k] - scalar;
//            }
//}
//
////calculate the variance
//float calculateVariance3D(boost::multi_array<float, 3> inputMatrix, float mean){
//    //subtract the mean from every element of the matrix
//    subtractScalarFromMatrix(inputMatrix, mean);
//    float variance = calculateSum(inputMatrix, 2);
//    return variance;
//}
//
//float calculateSum(boost::multi_array<float, 3> inputMatrix, double exponent){
//    const long unsigned int* matrixShape = inputMatrix.shape();
//    int matrixNrElements = inputMatrix.num_elements();
//    float sum = 0;
//
//    for(int i=0; i<matrixShape[0]; i++)
//        for(int j=0; j<matrixShape[1]; j++)
//            for(int k=0; k<matrixShape[2]; k++){
//                sum += pow(inputMatrix[i][j][k],exponent);
//            }
//    sum = sum/matrixNrElements;
//    return sum;
//}
//
//float calculateSkewness(boost::multi_array<float, 3> inputMatrix, float mean, float variance){
//    subtractScalarFromMatrix(inputMatrix, mean);
//    float numerator = calculateSum(inputMatrix, 3);
//    float denumerator = pow(variance, 1.5);
//    float skewness = numerator/denumerator;
//    return skewness;
//}
//
//float calculateKurtosis(boost::multi_array<float, 3> inputMatrix, float mean, float variance){
//    subtractScalarFromMatrix(inputMatrix, mean);
//    float numerator = calculateSum(inputMatrix, 4);
//    float denumerator = pow(variance, 2);
//    float kurtosis = numerator/denumerator - 3;
//    return kurtosis;
//}
//
//float getMedian(boost::multi_array<float, 3> inputMatrix){
//  typedef boost::multi_array<float, 3> array_type;
//  typedef array_type::index index;
//  array_type helpMatrix(boost::extents[3][4][2], boost::fortran_storage_order());
//  helpMatrix = inputMatrix;
////  call_fortran_function(helpMatrix.data());
//
//}
//
//float getMinimum(boost::multi_array<float, 3> inputMatrix){
//    float minimum = *std::min_element( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements());
//    return minimum;
//
//}
//
//float getMaximum(boost::multi_array<float, 3> inputMatrix){
//    float maximum = *std::max_element( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements());
//    return maximum;
//}
//
//float getRange(float minimum, float maximum){
//  typedef boost::multi_array_types::index_range range;
//
//array_type::array_view<3>::type myview = A[ boost::indices[range()][range()][range(1,2)] ];
//
////get column:std::cout<< A[boost::indices[0][range()][0]][0];
////get row: std::cout<< A[boost::indices[range()][0][0]][2];
