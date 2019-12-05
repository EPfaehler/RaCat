

#ifndef TEST_STATISTICALFEATURES_H_INCLUDED
#define TEST_STATISTICALFEATURES_H_INCLUDED

//#include <boost/test/unit_test.hpp>
//#include "intensityHistogram.h"
//#define BOOST_TEST_MODULE MyTest
//
//
//
//int add( int i, int j ) { return i+j; }
//
//BOOST_AUTO_TEST_CASE( my_test )
//{
//    // seven ways to detect and report the same error:
//    BOOST_CHECK( add( 2,2 ) == 4 );        // #1 continues on error
//
//    BOOST_REQUIRE( add( 2,2 ) == 4 );      // #2 throws on error
//
//    if( add( 2,2 ) != 4 )
//      BOOST_ERROR( "Ouch..." );            // #3 continues on error
//
//    if( add( 2,2 ) != 4 )
//      BOOST_FAIL( "Ouch..." );             // #4 throws on error
//
//    if( add( 2,2 ) != 4 ) throw "Ouch..."; // #5 throws on error
//
//    BOOST_CHECK_MESSAGE( add( 2,2 ) == 4,  // #6 continues on error
//                         "add(..) result: " << add( 2,2 ) );
//
//    BOOST_CHECK_EQUAL( add( 2,2 ), 4 );	  // #7 continues on error
//}
//
//BOOST_AUTO_TEST_CASE (test_statisticalFeatures){
//
//typedef boost::multi_array<float, 3> array_type;
//  typedef array_type::index index;
//  array_type A(boost::extents[3][4][2]);
//
//  // Assign values to the elements
//  int values = 0;
//  for(index i = 0; i != 3; ++i)
//    for(index j = 0; j != 4; ++j)
//      for(index k = 0; k != 2; ++k)
//        A[i][j][k] = values++;
//
//  // Verify values
//  int verify = 0;
//  for(index i = 0; i != 3; ++i)
//    for(index j = 0; j != 4; ++j)
//      for(index k = 0; k != 2; ++k)
//        assert(A[i][j][k] == verify++);
//
//  StatisticalFeatures<float, 3> statFeatures;
//  statFeatures.calculateMean(A);
//  BOOST_CHECK(statFeatures.meanValue = 11);
//}
#include <boost/test/unit_test.hpp>

void free_test_function();
//test_suite* init_unit_test_suite( int, char* [] );


#endif // TEST_STATISTICALFEATURES_H_INCLUDED
