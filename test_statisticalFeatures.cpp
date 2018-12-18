//
//
#include "statisticalFeatures.h"


#define BOOST_TEST_MODULE Stat_Features
#include <boost/test/unit_test.hpp>

using namespace std;
#include <iostream>

void getMatrix(boost::multi_array<float, 3> &imageMatrix){
    imageMatrix[0][0][0]=1;
    imageMatrix[1][0][0]=1;
    imageMatrix[2][0][0]=4;
    imageMatrix[3][0][0]=4;

    imageMatrix[0][1][0]=4;
    imageMatrix[1][1][0]=4;
    imageMatrix[2][1][0]=1;
    imageMatrix[3][1][0]=4;

    imageMatrix[0][2][0]=4;
    imageMatrix[1][2][0]=6;
    imageMatrix[2][2][0]=6;
    imageMatrix[3][2][0]=6;

    imageMatrix[0][3][0]=1;
    imageMatrix[1][3][0]=1;
    imageMatrix[2][3][0]=4;
    imageMatrix[3][3][0]=4;

    imageMatrix[0][4][0]=1;
    imageMatrix[1][4][0]=1;
    imageMatrix[2][4][0]=1;
    imageMatrix[3][4][0]=1;

    imageMatrix[0][0][1]=1;
    imageMatrix[1][0][1]=1;
    imageMatrix[2][0][1]=NAN;
    imageMatrix[3][0][1]=4;

    imageMatrix[0][1][1]=4;
    imageMatrix[1][1][1]=1;
    imageMatrix[2][1][1]=1;
    imageMatrix[3][1][1]=4;

    imageMatrix[0][2][1]=4;
    imageMatrix[1][2][1]=6;
    imageMatrix[2][2][1]=3;
    imageMatrix[3][2][1]=6;

    imageMatrix[0][3][1]=1;
    imageMatrix[1][3][1]=1;
    imageMatrix[2][3][1]=1;
    imageMatrix[3][3][1]=1;

    imageMatrix[0][4][1]=1;
    imageMatrix[1][4][1]=1;
    imageMatrix[2][4][1]=1;
    imageMatrix[3][4][1]=1;

    imageMatrix[0][0][2]=1;
    imageMatrix[1][0][2]=1;
    imageMatrix[2][0][2]=1;
    imageMatrix[3][0][2]=1;

    imageMatrix[0][1][2]=4;
    imageMatrix[1][1][2]=1;
    imageMatrix[2][1][2]=1;
    imageMatrix[3][1][2]=1;

    imageMatrix[0][2][2]=4;
    imageMatrix[1][2][2]=1;
    imageMatrix[2][2][2]=NAN;
    imageMatrix[3][2][2]=6;

    imageMatrix[0][3][2]=NAN;
    imageMatrix[1][3][2]=1;
    imageMatrix[2][3][2]=1;
    imageMatrix[3][3][2]=1;

    imageMatrix[0][4][2]=NAN;
    imageMatrix[1][4][2]=1;
    imageMatrix[2][4][2]=1;
    imageMatrix[3][4][2]=1;

    imageMatrix[0][0][3]=1;
    imageMatrix[1][0][3]=1;
    imageMatrix[2][0][3]=1;
    imageMatrix[3][0][3]=1;

    imageMatrix[0][1][3]=4;
    imageMatrix[1][1][3]=1;
    imageMatrix[2][1][3]=1;
    imageMatrix[3][1][3]=1;

    imageMatrix[0][2][3]=4;
    imageMatrix[1][2][3]=1;
    imageMatrix[2][2][3]=1;
    imageMatrix[3][2][3]=6;

    imageMatrix[0][3][3]=NAN;
    imageMatrix[1][3][3]=1;
    imageMatrix[2][3][3]=1;
    imageMatrix[3][3][3]=1;

    imageMatrix[0][4][3]=NAN;
    imageMatrix[1][4][3]=1;
    imageMatrix[2][4][3]=1;
    imageMatrix[3][4][3]=1;


}

void getVectorOfMatrixElementsNotNAN(vector<float> &vectorOfMatrixElements, boost::multi_array<float, 3> imageMatrix){

    for(int depth = 0; depth < imageMatrix.shape()[2]; depth++){
        for(int row = 0; row < imageMatrix.shape()[0]; row++){
            for(int col = 0; col < imageMatrix.shape()[1]; col++){
                if(!std::isnan(imageMatrix[row][col][depth])){
                    vectorOfMatrixElements.push_back(imageMatrix[row][col][depth]);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE (stat_features){
  StatisticalFeatures<float, 3> statFeatures;
  typedef boost::multi_array<float, 3> array_type;
  typedef array_type::index index;
  array_type A(boost::extents[4][5][4]);

   getMatrix(A);
   vector<float> vectorMatrElement;
   getVectorOfMatrixElementsNotNAN(vectorMatrElement, A);
  statFeatures.calculateAllStatFeatures(statFeatures, vectorMatrElement);
  BOOST_CHECK_CLOSE(double(statFeatures.meanValue), double(2.14865), 0.1);
  BOOST_CHECK_CLOSE(statFeatures.varianceValue, 3.04547, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.skewnessValue, 1.08382, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.kurtosisValue, -0.354619, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.medianValue, 1.04077, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.minimumValue, 1, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.maximumValue, 6, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.rangeValue, 5, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.percentile10, 1, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.percentile90, 5, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.interquartileRange, 2, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.quartileCoeff, 0.6, 1);
  BOOST_CHECK_CLOSE(statFeatures.coeffOfVar, 0.812198, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.energyValue, 567, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.rootMean, 2.76806, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.meanAbsDev ,1.55223, 0.1);
  BOOST_CHECK_CLOSE(statFeatures.robustMeanAbsDev, 1.11383, 0.1);
}

/*BOOST_AUTO_TEST_CASE(stat_min){
  StatisticalFeatures<float, 3> statFeatures;
  typedef boost::multi_array<float, 3> array_type;
  typedef array_type::index index;
  array_type A(boost::extents[5][4][4]);
  void getMatrix(boost::multi_array<float, 3> &A);

  statFeatures.getMinimum(A);
  BOOST_CHECK(statFeatures.minimumValue==2);
}

int main(){}*/



