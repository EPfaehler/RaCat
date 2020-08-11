#ifndef INTENSITYHISTOGRAM_INCLUDED
#define INTENSITYHISTOGRAM_INCLUDED



#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include "statisticalFeatures.h"
using namespace boost;
using namespace boost::accumulators;

/*! \file */
/*!
In the class IntensityHistogram the intensity histogram features are calculated. \n
The class inherits from the class StatisticalFeatures, as the majority of the features are the same. \n
The calculation of the feature values is done after discretizing the matrix values to a user specified bin number.
*/
template <class T,  size_t R>
class IntensityHistogram : StatisticalFeatures<T,R>{
    private:

        typedef accumulator_set<T, features<tag::density> > accIntensity;
		//vector containing the probabilities
        vector<T> probabilities;
		void getProbabilities(boost::multi_array<T, R> inputMatrix);
		void getSkewnessKurtosisIntHist();
		//feature values
        T entropy = NAN;
        T mode = NAN;
        T histUniformity = NAN;
        T maxHistGradient = NAN;
        T maxHistGradGreyValue = NAN;
        T minHistGradient = NAN;
        T minHistGradGreyValue = NAN;
		T kurtosisInt = NAN;
		T skewnessInt = NAN;
        vector<T> vectorOfMatrixElem;
        vector<double> nrElementsH;
        vector<T> diffGreyLevels;
        vector<double> probElements;

		vector<T> maxHistVecGradient;
		vector<T> minHistVecGradient;

		void getHistGradient();
		void extractIntenseData(vector<T> &intenseData, IntensityHistogram<T, R> intenseFeatures);
		
		void getNrElements(vector<double> &nrElements);
		void getHistUniformity();
		void getEntropy();
		void getMode();
		void getMaxHistGradient();
		void getMinHistGradient();
		void defineIntenseFeatures(vector<string> &features);
		void defineIntenseFeaturesOntology(vector<string> &features);

    public:
		IntensityHistogram() {
		}
		~IntensityHistogram() {
		}

		void calculateAllIntFeatures(IntensityHistogram<T, R> &intense, boost::multi_array<T, R> inputMatrix, vector<T> vectorOfMatrElements, vector<T> diffGrey);
		void writeCSVFileIntensity(IntensityHistogram<T, R> intensHist, string outputFolder);
		void writeOneFileIntensity(IntensityHistogram<T, R> intense, ConfigFile config);
};


/*!
The method getProbabilities calculates the probabilities for every element and stores the probabilities 
in the vector probElements \n
@param[in] inputMatrix: the original matrix of the VOI
*/
template <class T,  size_t R>
void IntensityHistogram<T,R>::getProbabilities(boost::multi_array<T, R> inputMatrix){
    //save matrix elements in a vector
    getNrElements(nrElementsH);
    probElements = nrElementsH;
    transform( probElements.begin(), probElements.end(),
                    probElements.begin(), bind2nd( std::divides<float>(), boost::size(vectorOfMatrixElem)) );
}

/*!
The method getMode calculates the mode of the distribution. \n
@param[in] inputMatrix: the original matrix of the VOI
*/
template <class T,  size_t R>
void IntensityHistogram<T,R>::getMode(){
	for (int i = 0; i < boost::size(probElements); i++) {
	}
    int indexMaxProb = std::distance(probElements.begin(), max_element(probElements.begin(), probElements.end()));
	sort(vectorOfMatrixElem.begin(), vectorOfMatrixElem.end());
    mode = vectorOfMatrixElem[indexMaxProb];
}

/*!
The method getNrElements gets for every grey level the amount of voxel with the specific grey level.\n
@param[in] vector nrElements: reference to a vector, where the number of different elements are stored
*/
template <class T,  size_t R>
void IntensityHistogram<T,R>::getNrElements(vector<double> &nrElements){
    //save matrix elements in a vector
    vector<T> tempMatrElements = vectorOfMatrixElem;
	T maxElement = *max_element(tempMatrElements.begin(), tempMatrElements.end());
    int nrActElement;
	int greyLevelIndex;
    std::sort(tempMatrElements.begin(), tempMatrElements.end());
    for(int nrGreyLevel = 0; nrGreyLevel < this->diffGreyLevels.size(); nrGreyLevel++){
		greyLevelIndex = this->diffGreyLevels[nrGreyLevel];
        nrActElement = 0;
        for(int matrIndex = 0; matrIndex < tempMatrElements.size(); matrIndex++ ){
			//std::cout << tempMatrElements[matrIndex] << " ";
            if(tempMatrElements[matrIndex] == greyLevelIndex){
                nrActElement+=1;
            }
        }
        nrElements.push_back(nrActElement);
    }
}

/*!
The method getEntropy calculates the entropy of the probabilities
*/
template <class T, size_t R>
void IntensityHistogram<T,R>::getEntropy(){
    entropy = 0;
    for(int i = 0; i<probElements.size(); i++){
		if (probElements[i] > 0) {
			entropy -= probElements[i] * std::log2(probElements[i]);
		}
    }

}

/*!
The method getHistUniformity calculates the uniformity of the histogram of the discretized grey levels.
*/
template<class T, size_t R>
void IntensityHistogram<T,R>::getHistUniformity(){
    histUniformity = for_each(probElements.begin(), probElements.end(), square_accumulate<T>()).result();
}

template<class T, size_t R>
void IntensityHistogram<T,R>::getHistGradient(){
	
	T actualGradient;
	vector<T> histGradient;
	for (int i = 1; i < boost::size(probElements); i++) {
		if (i == 1 ) {
			actualGradient = (nrElementsH[i ] - nrElementsH[i -1]);
		}
		else if (i == boost::size(probElements) - 1) {
			actualGradient = (nrElementsH[i] - nrElementsH[i - 1]);
		}
		else {
			actualGradient = (nrElementsH[i + 1] - nrElementsH[i-1])/2;
				
		}
		histGradient.push_back(actualGradient);
		
		
		
	}
	histGradient.push_back((nrElementsH[boost::size(probElements) - 1]- nrElementsH[boost::size(probElements) - 2]));
	maxHistVecGradient =histGradient;
	minHistVecGradient=histGradient;
	if (boost::size(maxHistVecGradient) == 0) {
		maxHistVecGradient.push_back(0);
		minHistVecGradient.push_back(0);
	}

}


template<class T, size_t R>
void IntensityHistogram<T,R>::getMaxHistGradient(){
	if (boost::size(maxHistVecGradient) == 1 && maxHistVecGradient[0] == 0) {
		maxHistGradient = 0;
		maxHistGradGreyValue = 0;
		std::cout << "The histogram gradients could not be calculated, check intensity values in VOI" << std::endl;
	}
	else {
		maxHistGradient = *max_element(maxHistVecGradient.begin(), maxHistVecGradient.end());
		maxHistGradGreyValue = max_element(maxHistVecGradient.begin(), maxHistVecGradient.end()) - maxHistVecGradient.begin();
		maxHistGradGreyValue = diffGreyLevels[maxHistGradGreyValue +1] ;
	}
}

template<class T, size_t R>
void IntensityHistogram<T,R>::getMinHistGradient(){
	if (boost::size(minHistVecGradient) == 1 && minHistVecGradient[0] == 0) {
		minHistGradient = 0;
		minHistGradGreyValue = 0;
		std::cout << "The histogram gradients could not be calculated, check intensity values in VOI" << std::endl;
	}
	else {
		minHistGradient = *min_element(minHistVecGradient.begin(), minHistVecGradient.end());
		minHistGradGreyValue = min_element(minHistVecGradient.begin(), minHistVecGradient.end()) - minHistVecGradient.begin();
		minHistGradGreyValue = diffGreyLevels[minHistGradGreyValue + 1] ;

	}
}

template <class T, size_t R>
void IntensityHistogram<T, R>::getSkewnessKurtosisIntHist() {
	
	vector<T> minusMean;
	minusMean = vectorOfMatrixElem;
	transform(minusMean.begin(), minusMean.end(), minusMean.begin(),
		std::bind2nd(minus<T>(), this->meanValue));
	skewnessInt = 0;
	kurtosisInt = 0;
	float tmpValue = 0;
	for (int i = 0; i < minusMean.size(); i++) {
		skewnessInt += pow(minusMean[i], 3);
		kurtosisInt += pow(minusMean[i], 4);
		tmpValue += pow(minusMean[i], 2);
	}
	float tmpValue2 = tmpValue / minusMean.size();
	skewnessInt = (skewnessInt / minusMean.size()) / pow(tmpValue2, 1.5);
	kurtosisInt = (kurtosisInt/minusMean.size())/ pow(tmpValue2, 2) -3;
	


}

template <class T,  size_t R>



void IntensityHistogram<T,R>::calculateAllIntFeatures(IntensityHistogram<T,R> &intense, boost::multi_array<T,R> inputMatrix, vector<T> vectorOfMatrElements, vector<T> diffGrey ){
  this->diffGreyLevels = diffGrey;
 
  vectorOfMatrixElem = vectorOfMatrElements;
  intense.getProbabilities(inputMatrix);
  intense.calculateMean(vectorOfMatrixElem);
  intense.calculateVariance();
  intense.calculateSkewness();
  getSkewnessKurtosisIntHist();
  intense.calculateKurtosis();
  intense.getMedian(vectorOfMatrixElem);
  intense.getMinimum(vectorOfMatrixElem);
  intense.getMaximum(vectorOfMatrixElem);
  intense.getMode();
  intense.getRange();
  intense.get10percentile(vectorOfMatrixElem);
  intense.get90percentile(vectorOfMatrixElem);
  intense.getInterquartileRange(vectorOfMatrixElem);
  intense.getQuartileCoeff();
  intense.getCoeffOfVar();
  intense.energy(vectorOfMatrixElem);
  intense.meanAbsoluteDev(vectorOfMatrixElem);
  intense.medianAbsoluteDev(vectorOfMatrixElem);
  intense.getRobustMeanAbsDev(vectorOfMatrixElem);
  
  intense.getEntropy();
  getHistGradient();
  intense.getHistUniformity();
  intense.getMaxHistGradient();
  intense.getMinHistGradient();
  vector<float>().swap(vectorOfMatrixElem);
  
}

template <class T,  size_t R>
void IntensityHistogram<T,R>::writeCSVFileIntensity(IntensityHistogram<T,R> intense, string outputFolder)
{
    string csvName = outputFolder + "_intensityHistogram.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream intenseCSV;
    intenseCSV.open (name, std::ios_base::app);
    vector<string> features;
    defineIntenseFeatures(features);

    vector<T> intenseData;
    extractIntenseData(intenseData, intense);
    for(int i = 0; i< intenseData.size(); i++){
        intenseCSV << "Intensity histogram" << "," << features[i] <<",";
        intenseCSV << intenseData[i];
        intenseCSV << "\n";
    }
    intenseCSV.close();
}

template <class T, size_t R>
void IntensityHistogram<T, R>::writeOneFileIntensity(IntensityHistogram<T, R> intense, ConfigFile config)
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

	ofstream intenseCSV;
	intenseCSV.open(name, std::ios_base::app);
	vector<string> features;
	

	vector<T> intenseData;
	extractIntenseData(intenseData, intense);
	if (config.csvOutput == 1) {
		defineIntenseFeatures(features);
		for (int i = 0; i < intenseData.size(); i++) {
			intenseCSV << "Intensity histogram" << "," << features[i] << ",";
			intenseCSV << intenseData[i];
			intenseCSV << "\n";
		}
	}
	else if (config.ontologyOutput == 1) {
		defineIntenseFeaturesOntology(features);
		for (int i = 0; i < intenseData.size(); i++) {
			intenseCSV << config.patientID << "," << config.patientLabel << "," << features[i] << ",";
			intenseCSV << intenseData[i] << "," << config.featureParameterSpaceName << "," << config.calculationSpaceName;
			intenseCSV << "\n";
		}
	}
	intenseCSV.close();
}

template <class T,  size_t R>

void IntensityHistogram<T,R>::defineIntenseFeatures(vector<string> &features){
    features.push_back("mean");
    features.push_back("variance");
    features.push_back("skewness");
    features.push_back("kurtosis");
    features.push_back("median");
    features.push_back("minimum");
	features.push_back("10th percentile");
	features.push_back("90th percentile");
    features.push_back("maximum");
    features.push_back("mode");
	features.push_back("Interquartile range");
    features.push_back("range");
	features.push_back("Mean absolut deviation");
	features.push_back("Robust mean absolute deviation");
	features.push_back("Median absolut deviation");
	features.push_back("Coefficient of variation");
    features.push_back("Quartile coefficient");
	features.push_back("Entropy");
	features.push_back("Uniformity");
    features.push_back("Energy");
    features.push_back("Maximum histogram gradient");
    features.push_back("Maximum histogram gradient grey level");
    features.push_back("Minimum histogram gradient");
    features.push_back("Minimum histogram gradient grey level");

}

template <class T, size_t R>

void IntensityHistogram<T, R>::defineIntenseFeaturesOntology(vector<string> &features) {

	features.push_back("Fih.mean");
	features.push_back("Fih.var");
	features.push_back("Fih.skew");
	features.push_back("Fih.kurt");
	features.push_back("Fih.median");
	features.push_back("Fih.min");
	features.push_back("Fih.P10");
	features.push_back("Fih.P90");
	features.push_back("Fih.max");
	features.push_back("Fih.mode");
	features.push_back("Fih.iqr");
	features.push_back("Fih.range");
	features.push_back("Fih.mad");
	features.push_back("Fih.rmad");
	features.push_back("Fih.medad");
	features.push_back("Fih.cov");
	features.push_back("Fih.qcod");
	features.push_back("Fih.entropy");
	features.push_back("Fih.uniformity");
	features.push_back("Fih.energy");
	features.push_back("Fih.max.grad");
	features.push_back("Fih.max.grad.gl");
	features.push_back("Fih.min.grad");
	features.push_back("Fih.min.grad.gl");

}

template <class T,  size_t R>

void IntensityHistogram<T,R>::extractIntenseData(vector<T> &intenseData, IntensityHistogram<T, R> intenseFeatures){
    intenseData.push_back(intenseFeatures.meanValue);
    intenseData.push_back(intenseFeatures.varianceValue);
    intenseData.push_back(intenseFeatures.skewnessInt);
    intenseData.push_back(intenseFeatures.kurtosisInt);
    intenseData.push_back(intenseFeatures.medianValue);
    intenseData.push_back(intenseFeatures.minimumValue);
	intenseData.push_back(intenseFeatures.percentile10);
	intenseData.push_back(intenseFeatures.percentile90);
    intenseData.push_back(intenseFeatures.maximumValue);
    intenseData.push_back(intenseFeatures.mode);
	intenseData.push_back(intenseFeatures.interquartileRange);
    intenseData.push_back(intenseFeatures.rangeValue);
	intenseData.push_back(intenseFeatures.meanAbsDev);
	intenseData.push_back(intenseFeatures.robustMeanAbsDev);
	intenseData.push_back(intenseFeatures.medianAbsDev);
	intenseData.push_back(intenseFeatures.coeffOfVar);
    intenseData.push_back(intenseFeatures.quartileCoeff);
	intenseData.push_back(intenseFeatures.entropy);
	intenseData.push_back(intenseFeatures.histUniformity);
    intenseData.push_back(intenseFeatures.energyValue);
    intenseData.push_back(intenseFeatures.maxHistGradient);
    intenseData.push_back(intenseFeatures.maxHistGradGreyValue);
    intenseData.push_back(intenseFeatures.minHistGradient);
    intenseData.push_back(intenseFeatures.minHistGradGreyValue);

}

#endif // INTENSITYHISTOGRAM_INCLUDED
