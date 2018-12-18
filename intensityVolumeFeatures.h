#ifndef INTENSITYVOLUMEFEATURES_H_INCLUDED
#define INTENSITYVOLUMEFEATURES_H_INCLUDED

/*! \file */
/*!
Intensity volume histogram features describe the relationship between the grey level i and the volume fraction which contains
at least grey level i or higher.\n
For this, the intensity histogram and the corresponding volume fractions have to be calculated. \n
If the VOI contains only one grey level, the intensity volume features are not calculated. 
*/

#include <iostream>
#include <vector>
#include "boost/multi_array.hpp"
#include "boost/range/combine.hpp"
#include "boost/foreach.hpp"
#include "image.h"
using namespace std;

template <class T,  size_t R>
class IntensityVolumeFeatures{
    private:
        vector<T> diffGreyLevels;
        T maxGreyLevel;
        T minGreyLevel;
        vector<T> greyLevelFraction;
        vector<T> fracVolume;

        T volAtIntFrac10;
        T volAtIntFrac90;
        T intAtVolFrac10;
        T intAtVolFrac90;

        T diffVolAtIntFrac;
        T diffIntAtVolFrac;

        T getVolumeAtIntFraction(double percent);
        T getIntAtVolFraction(double percent, vector<T> diffGreyLevels);
        void defineIntVolFeatures(vector<string> &features);
        void extractIntVolData(vector<T> &intVolData, IntensityVolumeFeatures<T, R> intVolFeatures);

    public:
        IntensityVolumeFeatures(){

        }
		~IntensityVolumeFeatures() {
			
		}
        void getFractionalVolume(boost::multi_array<T,R> inputMatrix, vector<T> vectorMatrElem);
        void getGreyLevelFraction(boost::multi_array<T,R> inputMatrix);
        void calculateAllIntensVolFeatures(IntensityVolumeFeatures<T,R> &intVolFeatures, boost::multi_array<T, R> inputMatrix, vector<T> vectorMatrElem);
        void writeCSVFileIntVol(IntensityVolumeFeatures<T,R> intVol, string outputFolder);
		void writeOneFileIntVol(IntensityVolumeFeatures<T, R> intVol, string outputFolder);
};


/*!
In the function getFractionalVolume the fractional volume of each grey level is calculated. \n
The vector fractional volume vector is filled in this function
@parameter[in]: boost multi_array input matrix: matrix containing intensity values of VOI
@parameter[in] vectorMatrElemen: vector containing all grey levels of VOI
*/
template <class T, size_t R>
void IntensityVolumeFeatures<T, R>::getFractionalVolume(boost::multi_array<T,R> inputMatrix, vector<T> vectorMatrElem){
    T actFracVolume;
    double nrElementsSmaller;
	double nrElementsNotNAN;

    for(int greyLevel = minGreyLevel; greyLevel < maxGreyLevel + 1; greyLevel ++){
        nrElementsSmaller = 0;
		nrElementsNotNAN = 0;
        for(int depth = 0; depth < inputMatrix.shape()[2]; depth++){
            for(int rows = 0; rows < inputMatrix.shape()[0]; rows++){
                for(int col = 0; col < inputMatrix.shape()[1]; col++){
					if(inputMatrix[rows][col][depth] < greyLevel && !isnan(inputMatrix[rows][col][depth])){
						//std::cout << inputMatrix[rows][col][depth] << std::endl;
                        nrElementsSmaller += 1;
                    }
					if (!isnan(inputMatrix[rows][col][depth])) {
						nrElementsNotNAN += 1;
					}
                }
            }
        }
        actFracVolume = 1 - nrElementsSmaller/nrElementsNotNAN;
        fracVolume.push_back(actFracVolume);
    }

}

/*!
In the function getGreyLevelFraction the grey level fraction is calculated and appended to the vector greyLevelFraction. \n
@parameter[in]: boost multi_array input matrix: matrix containing intensity values of VOI
@parameter[in] vectorMatrElemen: vector containing all grey levels of VOI
*/
template <class T, size_t R>
void IntensityVolumeFeatures<T, R>::getGreyLevelFraction(boost::multi_array<T,R> inputMatrix){
    T actGreyLevelFraction;

    for(int actGreyLevel = minGreyLevel; actGreyLevel < maxGreyLevel + 1; actGreyLevel++){
        actGreyLevelFraction = (actGreyLevel - minGreyLevel)/(maxGreyLevel - minGreyLevel);
        greyLevelFraction.push_back(actGreyLevelFraction);
    }

}

/*!
In the function getVolumeAtIntFraction calculates the volume at a certain intensity fraction for a certain percentage value.
@parameter[in] double percent: percentage value for which the volume fraction is calculated
*/
template <class T, size_t R>
T IntensityVolumeFeatures<T, R>::getVolumeAtIntFraction(double percent){
    vector<T> tempVector = greyLevelFraction;
    typename vector<T>::iterator it;
    typename vector<T>::iterator greaterThan;
    greaterThan = remove_if(tempVector.begin(), tempVector.end(), bind2nd(less<T>(), percent));
    tempVector.erase(greaterThan, tempVector.end());
    it=find(greyLevelFraction.begin(),greyLevelFraction.end(),tempVector[0]);
    int pos = distance(greyLevelFraction.begin(), it);
    return fracVolume[pos];
}


/*!
In the function getIntAtVolFraction calculates the intensity at a certain volume fraction for a certain percentage value.
@parameter[in] double percent: percentage value for which the volume fraction is calculated
*/
template <class T, size_t R>
T IntensityVolumeFeatures<T, R>::getIntAtVolFraction(double percent, vector<T> diffGreyLevels){
    vector<T> tempVector = fracVolume;
    typename vector<T>::iterator it;
    typename vector<T>::iterator greaterThan;

	int pos;
	if (tempVector[boost::size(tempVector) - 1] < percent) {
		greaterThan = remove_if(tempVector.begin(), tempVector.end(), bind2nd(greater<T>(), percent));
		tempVector.erase(greaterThan, tempVector.end());
		it = find(fracVolume.begin(), fracVolume.end(), tempVector[0]);
		pos = distance(fracVolume.begin(), it);

	}
	else {
		std::cout << "The frac volume is never smaller than 90 percent, error in intensity at volume fraction calculation" << std::endl;
		pos = 0;
	}
	return minGreyLevel + pos;

}


template <class T, size_t R>
void IntensityVolumeFeatures<T, R>::calculateAllIntensVolFeatures(IntensityVolumeFeatures<T,R> &intVolFeatures, boost::multi_array<T, R> inputMatrix, vector<T> diffGreyLevels){
	maxGreyLevel = 0;
	minGreyLevel = 10000000;
	for (int depth = 0; depth < inputMatrix.shape()[2]; depth++) {
		for (int rows = 0; rows < inputMatrix.shape()[0]; rows++) {
			for (int col = 0; col < inputMatrix.shape()[1]; col++) {
				if (inputMatrix[rows][col][depth] >maxGreyLevel) {
					maxGreyLevel = inputMatrix[rows][col][depth];
				}
				if (inputMatrix[rows][col][depth]<minGreyLevel){
					minGreyLevel = inputMatrix[rows][col][depth];
				}
			}
		}
	}
	if (maxGreyLevel != minGreyLevel) {
		intVolFeatures.getFractionalVolume(inputMatrix, diffGreyLevels);
		intVolFeatures.getGreyLevelFraction(inputMatrix);
		volAtIntFrac10 = intVolFeatures.getVolumeAtIntFraction(0.1);
		volAtIntFrac90 = intVolFeatures.getVolumeAtIntFraction(0.9);
		intAtVolFrac10 = intVolFeatures.getIntAtVolFraction(0.1, diffGreyLevels);
		intAtVolFrac90 = intVolFeatures.getIntAtVolFraction(0.9, diffGreyLevels);

		diffIntAtVolFrac = abs(intAtVolFrac90 - intAtVolFrac10);
		diffVolAtIntFrac = abs(volAtIntFrac90 - volAtIntFrac10);
	}
	else {
		std::cout << "The max and min value of the VOI are the same, the volume intensity features cannot be calculated." << std::endl;
	}
}

template <class T, size_t R>
void IntensityVolumeFeatures<T, R>::writeCSVFileIntVol(IntensityVolumeFeatures<T,R> intVol, string outputFolder)
{
    string csvName = outputFolder + "_intensityVolFeat.csv";
    char * name = new char[csvName.size() + 1];
    std::copy(csvName.begin(), csvName.end(), name);
    name[csvName.size()] = '\0';

    ofstream intVolCSV;
    intVolCSV.open (name);
    vector<string> features;
    defineIntVolFeatures(features);

    vector<T> intVolData;
    extractIntVolData(intVolData, intVol);
    for(int i = 0; i< intVolData.size(); i++){
        intVolCSV <<"intensity volume"<<","<< features[i] <<",";
        intVolCSV << intVolData[i];
        intVolCSV << "\n";
    }
    intVolCSV.close();
}

template <class T, size_t R>
void IntensityVolumeFeatures<T, R>::writeOneFileIntVol(IntensityVolumeFeatures<T, R> intVol, string outputFolder) {
	string csvName = outputFolder + ".csv";
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';

	ofstream intVolCSV;
	intVolCSV.open(name, std::ios_base::app);
	vector<string> features;
	defineIntVolFeatures(features);

	vector<T> intVolData;
	extractIntVolData(intVolData, intVol);
	for (int i = 0; i< intVolData.size(); i++) {
		intVolCSV << "intensity volume" << "," << features[i] << ",";
		intVolCSV << intVolData[i];
		intVolCSV << "\n";
	}
	intVolCSV.close();
}

template <class T, size_t R>
void IntensityVolumeFeatures<T, R>::defineIntVolFeatures(vector<string> &features){
    features.push_back("volume at int fraction 10");
    features.push_back("volume at int fraction 90");
    features.push_back("int at vol fraction 10");
    features.push_back("int at vol fraction 90");
	features.push_back("difference vol at int fraction");
    features.push_back("difference int at volume fraction");
   

}

template <class T, size_t R>
void IntensityVolumeFeatures<T, R>::extractIntVolData(vector<T> &intVolData, IntensityVolumeFeatures<T, R> intVolFeatures){
    intVolData.push_back(intVolFeatures.volAtIntFrac10);
    intVolData.push_back(intVolFeatures.volAtIntFrac90);
    intVolData.push_back(intVolFeatures.intAtVolFrac10);
    intVolData.push_back(intVolFeatures.intAtVolFrac90);
	intVolData.push_back(intVolFeatures.diffVolAtIntFrac);
    intVolData.push_back(intVolFeatures.diffIntAtVolFrac);
    
}
#endif // INTENSITYVOLUMEFEATURES_H_INCLUDED
