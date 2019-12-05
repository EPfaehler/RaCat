#ifndef STATISTICALFEATURES_INCLUDED
#define STATISTICALFEATURES_INCLUDED

#include "functional"
#include "image.h"
#include <iostream>
#include "boost/multi_array.hpp"
#include "math.h"
#include "matrixFunctions.h"
#include "helpFunctions.h"
#include "vectorFunctions.h"
#include <typeinfo>
#include <fstream>
#include <string>



#include <boost/bind.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/cstdlib.hpp>

using namespace boost::accumulators;


using namespace std;

/*! \file */

/*!
Statistical Features describe the distribution of grey levels within the Volume of interest (VOI).
This features are calculated without previous discretization. \n
The class Statistical Features calculates the statistical features of the whole 3D VOI. \n
Every features is an attribute of the class Statistical Features. \n
The feature values are calculated using the accumulator function of the boost library.
*/


template <class T,  size_t R>
class StatisticalFeatures{

    private:
        vector<T> vectorOfMatrixElem;
        typedef boost::accumulators::features <
        tag::mean, tag::variance, tag::median, tag::skewness, tag::kurtosis> Features;
        typedef accumulator_set <T, Features> Accumulator;
        Accumulator acc;


        void fillAccumulator(Accumulator &acc, vector<T> vectorMatrElem);
        void defineStatFeatures(vector<string> &features);
		void defineStatFeaturesOntology(vector<string> &features);
        void extractStatData(vector<T> &statData, StatisticalFeatures<T, R> statFeatures);


    public:
        //constructor
        StatisticalFeatures(){
        }
		~StatisticalFeatures() {
		}

        T meanValue;
        T varianceValue;
        T skewnessValue;
        T kurtosisValue;
        T medianValue;
        T minimumValue;
        T maximumValue;
        T rangeValue;
        T percentile10;
        T percentile90;
        T percentile25;
        T percentile75;
        T interquartileRange;
        T quartileCoeff;
        T coeffOfVar;
        T energyValue;
        T rootMean;
        T meanAbsDev;
        T medianAbsDev;
        T robustMeanAbsDev;


        void calculateMean(vector<T> vectorMatrElem);
        void calculateVariance();
        void calculateSkewness();
        void calculateKurtosis();
        void getMedian(std::vector<T> vectorMatrElement);
        void getMinimum(vector<T> matrixVector);
        void getMaximum(vector<T> matrixVector);
        void getRange();

        double getPercentile(vector<T> matrixVector, T probability);
        void get10percentile(vector<T> matrixVector);
        void get90percentile(vector<T> matrixVector);
        void getInterquartileRange(vector<T> matrixVector);
        void getQuartileCoeff();
        void getCoeffOfVar();
        void rootMeanSquare(vector<T> vectorMatrElem);
        void energy(vector<T> vectorMatrElem);
        void meanAbsoluteDev(vector<T> vectorMatrElem);
        void medianAbsoluteDev(vector<T> vectorMatrElem);
        void getRobustMeanAbsDev(vector<T> vectorMatrElem);

        void getGreaterElements(vector<T> &vectorOfMatrixElem, T &value);
        void getSmallerElements(vector<T> &vectorOfMatrixElem, T &value);

        void calculateAllStatFeatures(StatisticalFeatures<T,R> &stat, vector<T> vectorMatrElement);
        void writeCSVFileStatistic(StatisticalFeatures<T, R> stat, string outputFolder);
		void writeOneFileStatistic(StatisticalFeatures<T, R> stat, ConfigFile config);


};

/*!
\brief fillAccumulator
@param input: inputMatrix: matrix containing intensity values \n
@param input: Accumulator with different tags \n

The intensity values of the image are copied in the accumulator.
*/
template <class T,  size_t R>
void StatisticalFeatures<T, R>::fillAccumulator(Accumulator &acc, std::vector<T> vectorMatrElem){

    for_each( vectorMatrElem.begin(), vectorMatrElem.end(), boost::bind<void>( boost::ref(acc), _1 ) );
}

/*!
\brief calculateMean
\param[in] vectorOfMatrElement: array containing all intensitiy values within the VOI \n
The accumulator is filled calling the fillAccumulator function and the mean intensity value of the VOI is calculated \n
\f$
F_{mean} = \frac{1}{N_{v}}\sum{X_{gl}}
\f$
*/
template <class T,  size_t R>
void StatisticalFeatures<T, R>::calculateMean(std::vector<T> vectorMatrElement){
    fillAccumulator(acc, vectorMatrElement);
    meanValue = mean(acc);

}

/*!
\brief calculateVariance
The variance value is calculated \n
\f$
F_{var} = \frac{1}{N_{V}}\sum{(X_{gl}-\mu)^{2}}
\f$
*/
template <class T,  size_t R>
void StatisticalFeatures<T, R>::calculateVariance(){
    varianceValue = variance(acc);

}

/*!
\brief calculateSkewness
The skewness of the intensity distribution is calculated
\f$
F_{skew} = \frac{\frac{1}{N_{V}}\sum{(X_{gl}-\mu)^{3}}}{(\frac{1}{N_{V}}\sum{(X_{gl}-\mu)^{2}})^{\frac{3}{2}}}
\f$
*/
template <class T,  size_t R>
void StatisticalFeatures<T, R>::calculateSkewness(){
    skewnessValue = skewness(acc);
}

/*!
\brief calculateKurtosis
The kurtosis of the intensity distribution is calculated
\f$
F_{kur} = \frac{\frac{1}{N_{V}}\sum{(X_{gl}-\mu)^{4}}}{(\frac{1}{N_{V}}\sum{(X_{gl}-\mu)^{2}})^{2}}-3
\f$
*/
template <class T,  size_t R>
void StatisticalFeatures<T, R>::calculateKurtosis(){
    kurtosisValue = kurtosis(acc);
}


/*!
\brief getMedian
The median of the intensity distribution is calculated
*/
template <class T,  size_t R>
void StatisticalFeatures<T, R>::getMedian(std::vector<T> vectorMatrElement){
	vector<T> test = vectorMatrElement;
	int size = int(boost::size(test) / 2);
	std::sort(test.begin(), test.end());
    medianValue = test[size];
}

/*!
\brief getMinimum
The smallest element of all intensity values in the VOI is extracted as \f$
F_{min}
\f$
*/
template <class T,  size_t R>
void StatisticalFeatures<T, R>::getMinimum(vector<T> vectorOfMatrixElem){
    minimumValue = *min_element( vectorOfMatrixElem.begin(), vectorOfMatrixElem.end());
}

/*!
\brief getMaximum
The highest element of all intensity values in the VOI is extracted as \f$
F_{max}
\f$
*/
template <class T,  size_t R>
void StatisticalFeatures<T, R>::getMaximum(vector<T> vectorOfMatrixElem){
    maximumValue = *max_element( vectorOfMatrixElem.begin(), vectorOfMatrixElem.end());
}

/*!
\brief getRange
The range of the distribution is defined as:
\f$
F_{range} = F_{max} - F_{min}
\f$
*/
template <class T,  size_t R>
void StatisticalFeatures<T, R>::getRange(){
    rangeValue = (maximumValue-minimumValue);
}

/*!
\brief getPercentile
@param[in] matrixVector: vector containing all intensity values of the VOI
@param[in] probability: probability of the percentile, that should be calculated
The function is a help function to calculate percentiles of a certain probability \n
It uses the accumulator function of the boost library together with the p_square_quantile-tag
*/
template <class T,  size_t R>
double StatisticalFeatures<T, R>::getPercentile(vector<T> matrixVector, T probability){

	vector<T> vectorToSort = matrixVector;
	int size = int(probability*boost::size(vectorToSort));
	std::sort(vectorToSort.begin(), vectorToSort.end());
	return vectorToSort[size];
}

/*!
\brief get10Percentile
@param[in] matrixVector: vector containing all intensity values of the VOI
The function calculates the \f$ 10^{th} \f$ percentile \f$ P_{10} \f$.
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::get10percentile(vector<T> matrixVector){
    percentile10 = getPercentile(matrixVector, 0.1);
}

/*!
\brief get90Percentile
@param[in] matrixVector: vector containing all intensity values of the VOI
The function calculates the \f$ 90^{th} \f$ percentile \f$ P_{90} \f$.
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::get90percentile(vector<T> matrixVector){
    percentile90 =  getPercentile(matrixVector, 0.89);
    percentile90 = percentile90;
}

/*!
\brief getInterquartileRange
The interquartile range is defined as follows:
\f$
F_{intquarRange} = P_{75}-P_{25}
\f$
With \f$ P_{75} \f$ and \f$ P_{25} \f$ being the \f$ 75^{th} \f$ and \f$ 25^{th} \f$ percentile.
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::getInterquartileRange(vector<T> matrixVector){
    percentile25 = getPercentile(matrixVector, 0.25);
    percentile75 = getPercentile(matrixVector, 0.75);
    interquartileRange = (percentile75 - percentile25);
}

/*!
\brief getQuartileCoeff
The quartile coefficient of dispersion is defined as follows:
\f$
F_{quartCoeff} = \frac{P_{75}-P_{25}}{P_{75}+P_{25}}
\f$
With \f$ P_{75} \f$ and \f$ P_{25} \f$ being the \f$ 75^{th} \f$ and \f$ 25^{th} \f$ percentile. \n
It is another measurement for the dispersion of the distribution
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::getQuartileCoeff(){
    quartileCoeff = (percentile75 - percentile25)/(percentile75 + percentile25);
}


/*!
\brief getCoeffOfVar
The coefficient of variation is defined as follows:
\f$
F_{coeffOfVar} = \frac{ \sigma }{ \mu }
\f$
It measures for the dispersion of the distribution
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::getCoeffOfVar(){
    coeffOfVar = pow(this->varianceValue, 0.5)/this->meanValue;
}


/*!
\brief energy
The energy of the distribution is defined as follows:
\f$
F_{energy} = \sum{ X_{gl}^{2}}
\f$
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::energy(vector<T> vectorMatrElem){
    energyValue = for_each(vectorMatrElem.begin(), vectorMatrElem.end(), square_accumulate<float>()).result();
}

/*!
\brief rootMeanSquare
The root mean square is defined as follows:
\f$
F_{rootMeanSquare} = \sqrt{\frac{\sum{ X_{gl}^{2}}}{N_{V}}}
\f$
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::rootMeanSquare(vector<T> vectorMatrElem){
    rootMean = pow(this->energyValue/vectorMatrElem.size(),0.5);
}

/*!
\brief meanAbsoulteDev
The mean absolute deviation is calculated as:
\f$
F_{meanAbsoulteDev} =  \frac{1}{N_{V}}\sum{ X_{gl,j}^{2} - \mu }
\f$
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::meanAbsoluteDev(vector<T> vectorMatrElem){
    vector<T> tempVector = vectorMatrElem;
    transform( tempVector.begin(), tempVector.end(), tempVector.begin(),
                     std::bind2nd( minus<T>(), this->meanValue) );
    meanAbsDev = for_each(tempVector.begin(), tempVector.end(), sum_absol_value<float>()).result();
    meanAbsDev = meanAbsDev/tempVector.size();

}

/*!
\brief medianAbsoluteDev
The median absolute deviation is defined as:
\f$
F_{medianAbsoulteDev} =  \frac{1}{N_{V}}\sum{ X_{gl,j}^{2} - F_{median} }
\f$
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::medianAbsoluteDev(vector<T> vectorMatrElem){
    vector<T> tempVector = vectorMatrElem;
    transform( tempVector.begin(), tempVector.end(), tempVector.begin(),
                     std::bind2nd( minus<T>(), this->medianValue) );
    medianAbsDev = for_each(tempVector.begin(), tempVector.end(), sum_absol_value<float>()).result();
    medianAbsDev = medianAbsDev/tempVector.size();

}


/*!
\brief getRobustMeanAbsDev
Because outliers can have a big influence on the mean absolute deviation, the set of intensity values
included in the calculation of the robust mean can be limited to:
\f$
X_{10-90} =  \{ x \in X_{gl} | P_{10} (X_{gl}) \leq x \leq P_{90} (X_{gl}) \}
\f$
So, the set of grey levels is limited to the grey levels closer to the median of the distribution and the influence of
outliers is minimized. \n
\f$ X_{10-90} \f$ is the set of \f$ N_{10-90} \f$ grey level elements which lie in between (or are equal to) the \f$ 10^{th} \f$
and \f$ 90^{th} \f$ percentile. \n
The robust mean absolute deviation is calculated as follows:
\f$
F_{robmeanAbsDev} =  \frac{1}{N_{10-90}}\sum{ X_{gl10-90,j} -X_{gl10-90,j}}
\f$
*/
template <class T, size_t R>
void StatisticalFeatures<T, R>::getRobustMeanAbsDev(vector<T> vectorMatrElem){
    vector<T> tempVector = vectorMatrElem;
    getSmallerElements(tempVector, percentile10);
    getGreaterElements(tempVector, percentile90);

    accumulator_set<int, features<tag::mean> > accRobMean;
    accRobMean = for_each( tempVector.begin(), tempVector.end(), accRobMean );
    T meanValueRob = mean(accRobMean);
    transform(tempVector.begin(), tempVector.end(), tempVector.begin(),std::bind2nd( minus<T>(), meanValueRob) );
    robustMeanAbsDev = for_each(tempVector.begin(), tempVector.end(),sum_absol_value<float>()).result();
    robustMeanAbsDev = robustMeanAbsDev/tempVector.size();

}

/*!
\brief getSmallerElements
@param[in] vectorOfElem: vector of all elements in the VOI
@param[in] valueLimit: limit, in the end all values smaller than this value are stored in the vector
*/
template<class T, size_t R>
void StatisticalFeatures<T, R>::getSmallerElements(vector<T> &vectorOfElem, T &valueLimit){
    typename vector<T>::iterator lessThan;
    lessThan = std::remove_if(vectorOfElem.begin(), vectorOfElem.end(), std::bind2nd(less<T>(), valueLimit-1));
    vectorOfElem.erase(lessThan, vectorOfElem.end());
}

/*!
\brief getBiggerElements
@param[in] vectorOfElem: vector of all elements in the VOI
@param[in] valueLimit: limit, in the end all values higher than this value are stored in the vector
*/
template<class T, size_t R>
void StatisticalFeatures<T, R>::getGreaterElements(vector<T> &vectorOfElem, T &value){
    typename vector<T>::iterator greaterThan;
    greaterThan = std::remove_if(vectorOfElem.begin(), vectorOfElem.end(), std::bind2nd(greater<T>(), value));
    vectorOfElem.erase(greaterThan, vectorOfElem.end());
}




template <class T, size_t R>
void StatisticalFeatures<T, R>::calculateAllStatFeatures(StatisticalFeatures<T,R> &statFeatures, vector<T> vectorMatrElement){
   statFeatures.calculateMean(vectorMatrElement);
   statFeatures.getMedian(vectorMatrElement);
   statFeatures.calculateVariance();
   statFeatures.calculateSkewness();
   statFeatures.calculateKurtosis();

   statFeatures.getMinimum(vectorMatrElement);
   statFeatures.getMaximum(vectorMatrElement);
   statFeatures.getRange();
   statFeatures.get10percentile(vectorMatrElement);
   statFeatures.get90percentile(vectorMatrElement);
   statFeatures.getInterquartileRange(vectorMatrElement);
   statFeatures.getQuartileCoeff();
   statFeatures.getCoeffOfVar();
   statFeatures.energy(vectorMatrElement);
   statFeatures.rootMeanSquare(vectorMatrElement);
   statFeatures.meanAbsoluteDev(vectorMatrElement);
   statFeatures.medianAbsoluteDev(vectorMatrElement);
   statFeatures.getRobustMeanAbsDev(vectorMatrElement);


   
}

template <class T, size_t R>
void StatisticalFeatures<T, R>::writeCSVFileStatistic(StatisticalFeatures<T,R> stat, string outputFolder)
{

    string csvName = outputFolder + "_statisticalFeatures.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream statisticCSV;
    statisticCSV.open (name,std::ios_base::app);
    vector<string> features;
    defineStatFeatures(features);

    vector<T> statData;
    extractStatData(statData, stat);
    for(int i = 0; i< statData.size(); i++){
        statisticCSV << "Statistics" << "," <<features[i] <<",";
        statisticCSV << statData[i];
        statisticCSV << "\n";
    }
    statisticCSV.close();
}


template <class T, size_t R>
void StatisticalFeatures<T, R>::writeOneFileStatistic(StatisticalFeatures<T, R> stat, ConfigFile config)
{
	string csvName;
	if (config.csvOutput == 1) {
		csvName = config.outputFolder + ".csv";
	}
	else if (config.ontologyOutput == 1) {
		csvName = config.outputFolder + "/feature_table.csv";
	}
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream statisticCSV;
	statisticCSV.open(name, std::ios_base::app);
	vector<string> features;
	

	vector<T> statData;
	extractStatData(statData, stat);
	if (config.csvOutput == 1) {
		defineStatFeatures(features);
		for (int i = 0; i < statData.size(); i++) {
			statisticCSV << "Statistics" << "," << features[i] << ",";
			statisticCSV << statData[i];
			statisticCSV << "\n";
		}
		statisticCSV.close();
	}
	else if (config.ontologyOutput == 1) {
		defineStatFeaturesOntology(features);
		for (int i = 0; i < statData.size(); i++) {
			statisticCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			statisticCSV << statData[i] << "," << config.featureParameterSpaceName << "," << config.calculationSpaceName;
			statisticCSV << "\n";
		}
		statisticCSV.close();
	}
}

template <class T, size_t R>
void StatisticalFeatures<T, R>::defineStatFeaturesOntology(vector<string> &features) {
	features.push_back("Fstat.mean");
	features.push_back("Fstat.var");
	features.push_back("Fstat.skew");
	features.push_back("Fstat.kurt");
	features.push_back("Fstat.median");
	features.push_back("Fstat.min");		
	features.push_back("Fstat.P10");
	features.push_back("Fstat.P90");
	features.push_back("Fstat.max");
	features.push_back("Fstat.iqr");
	features.push_back("Fstat.range");
	features.push_back("Fstat.mad");
	features.push_back("Fstat.rmad");
	features.push_back("Fstat.medad");
	features.push_back("Fstat.cov");
	features.push_back("Fstat.qcod");
	features.push_back("Fstat.energy");
	features.push_back("Fstat.rms");

}

template <class T, size_t R>
void StatisticalFeatures<T, R>::defineStatFeatures(vector<string> &features){
    features.push_back("mean");
    features.push_back("variance");
    features.push_back("skewness");
    features.push_back("kurtosis");
    features.push_back("median");
    features.push_back("minimum");
	features.push_back("10th percentile");
	features.push_back("90th percentile");
    features.push_back("maximum");
	features.push_back("Interquartile range");
    features.push_back("range");
	features.push_back("Mean absolut deviation");
	features.push_back("Robust mean absolute deviation");
	features.push_back("Median absolute deviation");
	features.push_back("Coefficient of variation");
    features.push_back("Quartile coefficient");   
    features.push_back("Energy");
    features.push_back("Root mean");
    
}

template <class T, size_t R>
void StatisticalFeatures<T, R>::extractStatData(vector<T> &statData, StatisticalFeatures<T, R> statFeatures){
    statData.push_back(statFeatures.meanValue);
    statData.push_back(statFeatures.varianceValue);
    statData.push_back(statFeatures.skewnessValue);
    statData.push_back(statFeatures.kurtosisValue);
    statData.push_back(statFeatures.medianValue);
    statData.push_back(statFeatures.minimumValue);
	statData.push_back(statFeatures.percentile10);
	statData.push_back(statFeatures.percentile90);
    statData.push_back(statFeatures.maximumValue);
	statData.push_back(statFeatures.interquartileRange);
    statData.push_back(statFeatures.rangeValue);
	statData.push_back(statFeatures.meanAbsDev);
	statData.push_back(statFeatures.robustMeanAbsDev);
	statData.push_back(statFeatures.medianAbsDev);
	statData.push_back(statFeatures.coeffOfVar);
    statData.push_back(statFeatures.quartileCoeff);
    statData.push_back(statFeatures.energyValue);
    statData.push_back(statFeatures.rootMean);

}


#endif // STATISTICALFEATURES_INCLUDED
