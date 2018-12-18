#include "statisticalFeatures.h"

//in this part of the program the statistical radiomics features are defined

//calculate the mean of a matrix
//the matrix has to have the boost type of a multi_array,
//and it has to have 3 dimensions
// the size is not important
//template <class T>
//
//T& StatisticalFeatures<T>::calculateMean(boost::multi_array<T, 3> inputMatrix){
//    accumulator_set<double, stats<tag::mean > > acc;
//    acc = std::for_each( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), acc );
//    T meanMatrix = mean(acc);
//    return meanMatrix;
//}
//
//template <class T>
////calculate the variance
//T& StatisticalFeatures<T>::calculateVariance(boost::multi_array<T, 3> inputMatrix, T& mean){
//    accumulator_set<double, stats<tag::variance > > acc;
//    acc = std::for_each( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), acc );
//    T varianceMatrix = variance(acc);
//    return varianceMatrix;
//}
//
//template <class T>
//T& StatisticalFeatures<T>::calculateSkewness(boost::multi_array<T, 3> inputMatrix, T mean, T variance){
//    accumulator_set<double, stats<tag::skewness > > acc;
//    acc = std::for_each( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), acc );
//    T skewnessMatrix = skewness(acc);
//    return skewnessMatrix;
//}
//
//template <class T>
//T& StatisticalFeatures<T>::calculateKurtosis(boost::multi_array<T, 3> inputMatrix, T mean, T variance){
//    accumulator_set<double, stats<tag::kurtosis > > acc;
//    acc = std::for_each( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), acc );
//    T kurtosisMatrix = kurtosis(acc);
//    return kurtosisMatrix;
//}
//
//template <class T>
//T& StatisticalFeatures<T>::getMedian(boost::multi_array<T, 3> inputMatrix){
//    accumulator_set<double, stats<tag::median > > acc;
//    acc = std::for_each( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), acc );
//    T medianMatrix =median(acc);
//    return medianMatrix;
//}
//
//template <class T>
//T& StatisticalFeatures<T>::getMinimum(boost::multi_array<T, 3> inputMatrix){
//    T minimum = *std::min_element( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements());
//    return minimum;
//
//}
//
//template <class T>
//T& StatisticalFeatures<T>::getMaximum(boost::multi_array<T, 3> inputMatrix){
//    T maximum = *std::max_element( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements());
//    return maximum;
//}
//
//template <class T>
//T& StatisticalFeatures<T>::getRange(T minimum, T maximum){
//    return (maximum-minimum);
//}
//
//template <class T>
//T& StatisticalFeatures<T>::getPercentile(boost::multi_array<T, 3> inputMatrix, T probability){
//    accumulator_set<double, stats<tag::median > > acc;
//    acc = std::for_each( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), acc );
//    typedef accumulator_set<double, stats<tag::p_square_quantile> > accumulator_t;
//    accumulator_t acc0(quantile_probability = probability);
//    const long unsigned int* matrixShape = inputMatrix.shape();
//    for(int i = 0; i<matrixShape[0]; i++)
//        for(int j = 0; j<matrixShape[1]; j++)
//            for(int k = 0; k<matrixShape[2]; k++){
//            acc0(inputMatrix[i][j][k]);
//    }
//    return p_square_quantile(acc0);
//}
//
//template <class T>
//T& StatisticalFeatures<T>::get10percentile(boost::multi_array<T, 3> inputMatrix){
//    return getPercentile(inputMatrix, 0.1);
//}
//
//template <class T>
//T& StatisticalFeatures<T>::get90percentile(boost::multi_array<T, 3> inputMatrix){
//    return getPercentile(inputMatrix, 0.9);
//}
//
//template <class T>
//T& StatisticalFeatures<T>::getInterquartileRange(boost::multi_array<T, 3> inputMatrix){
//    T perc25 = getPercentile(inputMatrix, 0.25);
//    T perc75 = getPercentile(inputMatrix, 0.75);
//    return perc75 - perc25;
//}
//
//template <class T>
//T& StatisticalFeatures<T>::getQuartileCoeff(boost::multi_array<T, 3> inputMatrix){
//    T perc25 = getPercentile(inputMatrix, 0.25);
//    T perc75 = getPercentile(inputMatrix, 0.75);
//    return (perc75 - perc25)/(perc75 + perc25);
//
//}
//
//template <class T>
//T& StatisticalFeatures<T>::getCoeffOfVar(T variance, T mean){
//    return pow(variance, 0.5)/mean;
//}
//
//template <typename T>
//class square_accumulate
//{
//public:
//    square_accumulate(void) :
//      _sum(0)
//      {
//      }
//
//      const T& result(void) const
//      {
//          return _sum;
//      }
//
//      void operator()(const T& val)
//      {
//          _sum += val * val;
//      }
//
//private:
//    T _sum;
//};
//
//template <class T>
//T& StatisticalFeatures<T>::energy(boost::multi_array<T, 3> inputMatrix){
//    T sum = std::for_each(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), square_accumulate<float>()).result();
//    return sum;
//}
//
//template <class T>
//T& StatisticalFeatures<T>::rootMeanSquare(boost::multi_array<T, 3> inputMatrix, T energy){
//    T rootMean = pow(energy/inputMatrix.num_elements(),0.5);
//    return rootMean;
//}
//
//template <class T>
//T& StatisticalFeatures<T>::accumulate_func(T& accumulator, T& scalar){
//   return abs(accumulator - scalar);
// };
//
//template <class T>
//T& StatisticalFeatures<T>::meanAbsoluteDev(boost::multi_array<T, 3> inputMatrix, T mean){
//    T meanAbsDev = std::accumulate(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), 0, accumulate_func);
//    meanAbsDev = meanAbsDev/inputMatrix.num_elements();
//    return meanAbsDev;
//}
//
//
