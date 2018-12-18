//#include "greyLevelDiscretization.h"
//
//#include <iostream>
//using namespace std;
//
//template <typename T>
//class rounding
//{
//public:
//    rounding(void) :
//      _roundRes(0)
//      {
//      }
//
//      const T& result(void) const
//      {
//          return _roundRes;
//      }
//
//      void operator()(const T& val)
//      {
//      _roundRes = ceil(val);
//    }
//
//
//private:
//    T _roundRes;
//};
//
//
//
//double myfunction (double i) {
//  return ceil(i);
//}
//template <class T,  size_t R>
//double getInfo(T i){
//    std::cout<<i;
//}
//template <class T,  size_t R>
//void equalIntervalWidth(boost::multi_array<T, R> inputMatrix, double intervalWidth){
//
//    float minimumValue = *min_element( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements());
//    transform( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(),
//                    inputMatrix.origin(), bind2nd( minus<T>(), minimumValue) );
//    transform( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(),
//                    inputMatrix.origin(), bind2nd( std::divides<T>(), intervalWidth) );
//    replace(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), 0, 1);
//    std::transform(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), inputMatrix.origin(), myfunction);
//    //std::for_each(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), getInfo);
//}

//template <class T,  size_t R>
//void equalBinSize(boost::multi_array<T, R> inputMatrix){
//    float minimumValue = *min_element( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements());
//    float maximumValue = *max_element( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements());
//    transform( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(),
//                    inputMatrix.origin(), bind2nd( minus<T>(), minimumValue) );
//    float range = maximumValue-minimumValue;
//    transform( inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(),
//                    inputMatrix.origin(), bind2nd( std::divides<T>(), range) );
//    transform(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(),
//                        inputMatrix.origin(), bind2nd( std::multiplies<T>(), inputMatrix.num_elements()));
//    replace(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), 0, 1);
//    std::transform(inputMatrix.origin(), inputMatrix.origin() + inputMatrix.num_elements(), inputMatrix.origin(), myfunction);
//
//}
