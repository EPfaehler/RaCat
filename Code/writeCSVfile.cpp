#include "writeCSVfile.h"



//void writeCSVFileStatistic(StatisticalFeatures<float,3> stat)
//{
//    ofstream statisticCSV;
//    statisticCSV.open ("example.csv");
//    vector<string> features;
//    defineStatFeatures(features);
//
//    vector<float> statData;
//    extractStatData(statData, stat);
//    std::cout<<statData.size();
//    for(int i = 0; i< statData.size(); i++){
//        statisticCSV << features[i] <<";";
//        statisticCSV << statData[i];
//        statisticCSV << "\n";
//    }
//    statisticCSV.close();
//}
//
//
//void defineStatFeatures(vector<string> &features){
//    features.push_back("mean");
//    features.push_back("variance");
//    features.push_back("skewness");
//    features.push_back("kurtosis");
//    features.push_back("median");
//    features.push_back("minimum");
//    features.push_back("maximum");
//    features.push_back("range");
//    features.push_back("10th percentile");
//    features.push_back("90th percentile");
//    features.push_back("25th percentile");
//    features.push_back("75th percentile");
//    features.push_back("Interquartile range");
//    features.push_back("Quartile coefficient");
//    features.push_back("Coefficient of variance");
//    features.push_back("Energy");
//    features.push_back("Root mean");
//    features.push_back("Mean absolut deviation");
//    features.push_back("Robust mean absolute deviation");
//}
//
//void extractStatData(vector<float> &statData, StatisticalFeatures<float, 3> statFeatures){
//    statData.push_back(statFeatures.meanValue);
//    statData.push_back(statFeatures.varianceValue);
//    statData.push_back(statFeatures.skewnessValue);
//    statData.push_back(statFeatures.kurtosisValue);
//    statData.push_back(statFeatures.medianValue);
//    statData.push_back(statFeatures.minimumValue);
//    statData.push_back(statFeatures.maximumValue);
//    statData.push_back(statFeatures.rangeValue);
//    statData.push_back(statFeatures.percentile10);
//    statData.push_back(statFeatures.percentile90);
//    statData.push_back(statFeatures.percentile25);
//    statData.push_back(statFeatures.percentile75);
//    statData.push_back(statFeatures.interquartileRange);
//    statData.push_back(statFeatures.quartileCoeff);
//    statData.push_back(statFeatures.coeffOfVar);
//    statData.push_back(statFeatures.energyValue);
//    statData.push_back(statFeatures.rootMean);
//    statData.push_back(statFeatures.meanAbsDev);
//    statData.push_back(statFeatures.robustMeanAbsDev);
//
//}

//void writeCSVFileStatistic(IntensityHistogram<float,3> intense)
//{
//    ofstream intenseCSV;
//    intenseCSV.open ("intensityHistogram.csv");
//    vector<string> features;
//    defineIntenseFeatures(features);
//
//    vector<float> intenseData;
//    extractIntenseData(intenseData, intense);
//    for(int i = 0; i< intenseData.size(); i++){
//        intenseCSV << features[i] <<";";
//        intenseCSV << intenseData[i];
//        intenseCSV << "\n";
//    }
//    intenseCSV.close();
//}


//void defineIntenseFeatures(vector<string> &features){
//    features.push_back("mean");
//    features.push_back("variance");
//    features.push_back("skewness");
//    features.push_back("kurtosis");
//    features.push_back("median");
//    features.push_back("minimum");
//    features.push_back("maximum");
//    features.push_back("range");
//    features.push_back("10th percentile");
//    features.push_back("90th percentile");
//    features.push_back("25th percentile");
//    features.push_back("75th percentile");
//    features.push_back("Interquartile range");
//    features.push_back("Quartile coefficient");
//    features.push_back("Coefficient of variance");
//    features.push_back("Energy");
//    features.push_back("Root mean");
//    features.push_back("Mean absolut deviation");
//    features.push_back("Robust mean absolute deviation");
//}

//void extractIntenseData(vector<float> &intenseData, IntensityHistogram<float, 3> intenseFeatures){
//    intenseData.push_back(intenseFeatures.meanValue);
//    intenseData.push_back(intenseFeatures.varianceValue);
//    intenseData.push_back(intenseFeatures.skewnessValue);
//    intenseData.push_back(intenseFeatures.kurtosisValue);
//    intenseData.push_back(intenseFeatures.medianValue);
//    intenseData.push_back(intenseFeatures.minimumValue);
//    intenseData.push_back(intenseFeatures.maximumValue);
//    intenseData.push_back(intenseFeatures.rangeValue);
//    intenseData.push_back(intenseFeatures.percentile10);
//    intenseData.push_back(intenseFeatures.percentile90);
//    intenseData.push_back(intenseFeatures.percentile25);
//    intenseData.push_back(intenseFeatures.percentile75);
//    intenseData.push_back(intenseFeatures.interquartileRange);
//    intenseData.push_back(intenseFeatures.quartileCoeff);
//    intenseData.push_back(intenseFeatures.coeffOfVar);
//    intenseData.push_back(intenseFeatures.energyValue);
//    intenseData.push_back(intenseFeatures.rootMean);
//    intenseData.push_back(intenseFeatures.meanAbsDev);
//    intenseData.push_back(intenseFeatures.robustMeanAbsDev);
//
//}
