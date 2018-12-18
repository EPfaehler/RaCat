#ifndef READCONFIGFILE_H_INCLUDED
#define READCONFIGFILE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <ctime>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#ifdef _WIN32
#include<windows.h>

#endif
#include <sys/types.h>
#include <sys/stat.h>

#include "string"
#include <fstream>

#include "featureCalculation.h"
using namespace std;

typedef boost::property_tree::ptree config;
/*! \file */


/*! @page config Setting up the config.ini file
In the config.ini file you can set the preprocessing steps you want to use.

\arg Smoothing: \n
If you want to use additional smoothing to the image, you can set a Gaussian kernel to the desired full-width-at-half-maximum. Smoothing is applied as soon as the Smoothing kernel is not equal to 0. \n

\arg Threshold (in %) for including voxels in the VOI:
Masks can contain different values. A mask can contain only 1s, values from 1-100, or other ranges. You can determine which voxels will be included in the final mask by setting a threshold. The program determines the
maximum value inside the mask and includes all values in the final mask which have a value higher than this threshold from the maximum values. The recommended value is 0.5.\n
The value is especially important when you apply resampling to the image and the mask. 

\arg Discretization method: \n
Before textural features are calculated, the image is normally discretized. You can choose between two different methods for discretization:
- use fixed bin width and/or \n
- use fixed nummer of bins \n
You can set the bin width and the number of bins you want to use.  \n

\arg Interpolation method: \n
If you want to up- or downsample the image, you can say it in the part "interpolation method". You can choose between up- or downsampling or if you want to sample the image to cubic voxels of 2 mm voxel size. \n
For interpolation, trilinear interpolation is used. \n

\arg Resegmentation: \n
If you want to resegment the image region, you can state that here. I.e. if you want to exclude values above/below a maximum/minimum value or if you want to exclude outliers from the segmentation.\n
Herefore, you have to set ReSegmentImage to 1 and set the minimum and maximum value you want to include. \n
For PET images, the minimum and maximum values have to be in kBq/ml. \n


\arg Image properties: \n
Here you can set the image type you are using. In general, the program is only interested if you have a PET image or not:  \n
If the image is a PET image, a patientInfo.ini file is required. In this file, you can set values typical for PET images. \n
If another image type is set, the patientInfo.ini file is not required and will not be read in. \n

\arg NGLDMParameters: \n
The neighbourhood grey level dependence matrix (NGLDM) captures the texture of the matrix. It compares the intensity values of neighboring voxels
with the intensity value of a center voxel. If the difference between these intensity values is smaller than a threshold (coarseness parameter), the
voxels are defined to be independent. The default coarseness parameter is 0, but you can change it to other values in
this section.\n
Furthermore, you can define the size of a neighborhood: dist is the distance of a voxel to the center voxel. All voxels with a distance smaller or equal
to this distance are regarded as part of the neighborhood. The default value is 1.\n

\arg NGTDM parameters: \n
You can set the size (distance to center voxel) of a neighborhood. The default value is also 1. \n

\arg Distance weight properties: \n
The contents of GLCM, GLRLM and NGTDM matrices can be updated by the distance between voxels. Herefore, different definitions of 
the distance can be used. Possible distances are: \n
- Chebyshev (default) \n
- Manhattan \n
- Euclidean norm \n
If another distance is set, the program gives an error. \n

\arg Extended Emphasis Features: \n
In the GLRLM, GLSZM, GLDZM and NGLDM matrices, you can set a particular emphasis on part of these matrices. In the formula of
the feature calculations, the emphasis can be set by different powers. The desired powers can be set here. \n
If the extended emphasis features should be calculated, the CalculateExtendedEmph value should be set to 1. \n

\arg Output information: \n
Here, you can set the output format. You can choose between csv or ontology output. If you want to have a csv file as output, you can decide if you want to have it as one csv file or 
one csv file per feature group.\n
 
The class ConfigFile has as attributes exactly the attributes which can be set in the config-file. \n
@page PatientInfo Setting up the patientInfo.ini file
The patientInfo.ini file is only required if a PET image is used as input. You can change the following parameters: \n
\arg  Patient weight\n
The patient weight is needed in order to calculate the SUV values.\n
The SUV values are only calculated if the parameter UseSUV is set to 1.\n
If the scaling parameter is set to a value \b not equal to 0, the weight is not needed.
\arg  Patient height \n
The patient height is needed in order to calculate the SUL values.\n
The SUL values are only calculated if the parameter UseSUL is set to 1. \n
If the scaling parameter is set to a value \b not equal to 0, the height is not needed.
\arg  ActivityMBq \n
The injected dose (in MBq) at scan \b START. \n
Please note that no decay correction is done by the program. Therefore the dose has to be the one at scan start. \n
If the scaling parameter is set to a value \b not equal to 0, the activity is ignored.
\arg Gender \n
he gender is needed in order to calculate the SUL values.\n
Possible Values are \b M or \b F. All other values are ignored. \n
The SUL values are only calculated if the parameter UseSUL is set to 1. \n
If the scaling parameter is set to a value \b not equal to 0, the height is not needed.
\arg UseSUV \n
The value has to  be set to 1, if the image should be converted in SUV. 
\arg UseSUL \n
The value has to  be set to 1, if the image should be converted to SUL.
\arg ScalingFactor \n
If the scaling factor is set to a value \b not equal to 0, all other values will be ignored and the image will be scaled by
this factor.\n
UseSUV has to be set to 1. 
*/


class ConfigFile{
    public:
        string fileName;
        inline config readIni(string iniFile);
        config pt;
		//!float for smoothing kernel if !=0, image is smoothed
		float smoothingKernel;
		//!threshold in %; all voxels which contain values higher than this percentage of the maximum mask value will be included in the mask
		float threshold;
		//!integer which states if exe is called by accurate tool
		int useAccurate;
		//!integer which states if .voi is dicom, nii or else
		int voiFile;
        //!integer which states if we use fixed bin width
        int useFixedBinWidth;
        float binWidth;
        //!integer which states if we use fixed number of bins
        int useFixedNrBins;
        int nrBins;
		//integer to state if IVH should be discretized separately
		int discretizeIVH;
		int discretizeIVHSeparated;
		int useFixedBinWidthIVH;
		int binWidthIVH;
		int useFixedNrBinsIVH;
		int nrBinsIVH;
		//! value to state if 2D or 3D interpolation, if set to 0 3D interpolation is applied
		int interpolation2D;
		string interpolationMethod;
        // !integer which states if we use down sampling
        int useDownSampling;
		int useUpSampling;
		int useSampling2mm;
		//! integer which states if resegmentation is used
		int useReSegmentation;
		int excludeOutliers;
		float minValueReSeg;
		float maxValueReSeg;
		// !norms set by the user to calculate the matrices
		string normGLCM;
		string normGLRLM;
		string normNGTDM;
		//! shall the extended emphasis features be calculated?
		int extendedEmphasis;
		int powerRow;
		int powerCol;
		//!parameters for NGLDM matrices
		int distNGLDM;
		int coarsenessParam;
		//!distance defined by user for NGTDM matrices
		int dist;
		string featureSelectionLocation;
		int calculateAllFeatures = 0;
		string patientInfoLocation;
        //!names of images and folders
        string imageName;
        string voiName;
        string imageType;
		//!name of folder, where values are stored
        string outputFolder;
		//!integer which states which kind of output
		int getOneCSVFile;
		int csvOutput;
		int ontologyOutput;
		//!parameters to calculate the SUV value in case we have a PET image; this values can be set in the patientInfo.ini file
        int useSUV;
        int useSUL;
        float patientWeight;
        float patientHeight;
		int malePatient;
		std::string gender;
        float initActivity;
        int minAfterInjection;
		//!the correction parameter overwrites the SUV values set before (if bigger than 1)
		//!all image values will be divided by it, instead of dividing by calculated SUV
		double correctionParam;
		//!get smoothing kernel
		void getSmoothingKernel();
		//!get threshold
		void getThreshold();
		//!get information if we are working with prj file
		void getAccurateState(string accState);
		//!get information with what kind of .voi we are working
		inline void getVoiState(string accState);
		//! get information about resegmentation
		inline void getResegmentationState();
		//!get the location of the featureSelection.ini
		void getFeatureSelectionLocation(string featurePath);
        //! read the discretization information
        void getDiscretizationInformation();
		//! read discretization information IVH
		void getDiscretizationInformationIVH();
		//! read distance weight properties
		void getDistanceWeightProperties();
		//! get extended emphasis information
		void getExtendedEmphasisInformation();
		//! get the NGLDM parameters
		void getNGLDMParameters();
		//! get the NGTDM distance value
		void getNGTDMdistanceValue();
        //! read the information about the folders where the images are saved
        void getImageFolder(string imageName, string voiName);
		//! read information about outputFolder location
        void getOutputInformation(string output);
		//! read information necessary to convert image from kBq/ml to SUV
        void getPETimageInformation(string imagePath, string patientInfoPath, ConfigFile config);
        //! read information about interpolation
        void getInterpolation();
		//! copies the config file to the output folder
		void copyConfigFile(string outputFolder);
		void createOutputFolder(ConfigFile &config);
		//! adjusts all values to the attributes in the ConfigFile
		void createConfigInfo(ConfigFile &config, string arguments[7]);
};

/*!
The method readIni reads the desired ini-file
@param[in]: string iniName: path + name of the iniFile (given user)
@param[out]: boost::property_tree pt: the property tree of the ini-file
*/
inline config ConfigFile::readIni(string iniName) {
	config pt;
	fileName = iniName;
	boost::property_tree::ini_parser::read_ini(iniName, pt);
	return pt;
}

/*!
The method getSmoothingKernel reads the value of the smoothing kernel. 
*/
inline void ConfigFile::getSmoothingKernel() {
	config pt = readIni(fileName);
	smoothingKernel = pt.get<float>("Smoothing.SmoothingKernel", 0.0);
}


/*!
The method getSmoothingKernel reads the value of the smoothing kernel.
*/
inline void ConfigFile::getThreshold() {
	config pt = readIni(fileName);
	threshold = pt.get<float>("ThresholdForVOI.threshold", 0.5);
}

/*!
The method getAccurateState sets the useAccurate value. \n
If we have a prj file, the value is set to 1. \n
If it is a dicom image, the value is set to 2. \n
Otherwise the value is 0.
*/
inline void ConfigFile::getAccurateState(string accState) {
	config pt = readIni(fileName);
	if (accState == "acc") {
		
		std::cout<< "The input image file is a prj file" << std::endl;
		useAccurate = 1;
	}
	else if (accState == "dic") {
		useAccurate = 2;
	}
	else {
		useAccurate = 0;
	}
}

/*!
The method getVoiState sets the voiFile value. It is 0 if we have a nifti or .voi file, if its a dicom image, the value is set to 2 
and if it is a rt struct, the value is set to 3.
*/
inline void ConfigFile::getVoiState(string voiState) {
	if (voiState == "dic") {
		voiFile = 2;
	}
	else if (voiState == "rts") {
		voiFile = 3;
	}
}

/*!
The method getResegmentationState reads the provided resampling information. 
*/
inline void ConfigFile::getResegmentationState() {
	config pt = readIni(fileName);
	useReSegmentation = pt.get<int>("ReSegmentation.ReSegmentImage");
	excludeOutliers = pt.get<int>("ReSegmentation.ExcludeOutliers");
	minValueReSeg = pt.get<float>("ReSegmentation.MinValueInReSegmentation");
	maxValueReSeg = pt.get<float>("ReSegmentation.MaxValueInReSegmentation");
}

/*!
The method getFeatureSelectionLocation reads the location of the feature selection file
*/
inline void ConfigFile::getFeatureSelectionLocation(string featPath) {
	config pt = readIni(fileName);
	
	if (featPath != "0") {
		featureSelectionLocation = featPath;
	}
	else {
		calculateAllFeatures = 1;
	}
}

/*!
The method getDiscretizationInformation reads the discretization information of the ini-file and
sets the attributes of the class Config to the equivalent values
*/
inline void ConfigFile::getDiscretizationInformation() {
	config pt = readIni(fileName);
	useFixedBinWidth = pt.get("Discretization.UseFixedBinWidth", 1);
	if (useFixedBinWidth != 0 && useFixedBinWidth != 1) {
		std::cout << "You inserted a value for useFixedBinWidth that is not 0 or 1, it will be set to 1" << std::endl;
		useFixedBinWidth = 1;
	}
	if (useFixedBinWidth == 1) {
		binWidth = pt.get<float>("Discretization.BinWidth", 0.25);
		if (binWidth < 0) {
			std::cout << "please fill in a valid number for the bin width" << std::endl;
		}
	}

	useFixedNrBins = pt.get("Discretization.UseFixedNrBins", 1);
	if (useFixedNrBins != 0 && useFixedNrBins != 1) {
		std::cout << "You inserted a value for UseFixedNrBins that is not 0 or 1, it will be set to 1" << std::endl;
		useFixedNrBins = 1;
	}
	if (useFixedNrBins == 1) {
		nrBins = pt.get("Discretization.NrBins", 64);
		if (nrBins < 0) {
			std::cout << "please fill in a valid number for the number of bins" << std::endl;
		}
	}
}

/*!
The method getDiscretizationInformationIVH reads the discretization information for the intensity volume histogram features and
sets the attributes of the class Config to the equivalent values
*/
inline void ConfigFile::getDiscretizationInformationIVH() {
	config pt = readIni(fileName);
	discretizeIVH = pt.get("DiscretizationIVH.DiscretizeIVH", 0);
	discretizeIVHSeparated = pt.get("DiscretizationIVH.DiscretizeIVHSeparated", 0);
	useFixedBinWidthIVH = pt.get("DiscretizationIVH.UseFixedBinWidthIVH", 1);
	if (useFixedBinWidthIVH != 0 && useFixedBinWidthIVH != 1 && discretizeIVHSeparated == 1) {
		std::cout << "You inserted a value for useFixedBinWidthIVH that is not 0 or 1, it will be set to 1" << std::endl;
		useFixedBinWidthIVH = 1;
	}
	if (useFixedBinWidthIVH == 1) {
		binWidthIVH = pt.get<float>("DiscretizationIVH.BinWidthIVH", 0.25);
		if (binWidthIVH < 0) {
			std::cout << "please fill in a valid number for the bin width" << std::endl;
		}
	}

	useFixedNrBinsIVH = pt.get("DiscretizationIVH.UseFixedNrBinsIVH", 1);
	if (useFixedNrBinsIVH != 0 && useFixedNrBinsIVH != 1 && discretizeIVHSeparated ==1) {
		std::cout << "You inserted a value for UseFixedNrBinsIVH that is not 0 or 1, it will be set to 1" << std::endl;
		useFixedNrBinsIVH = 1;
	}
	if (useFixedNrBinsIVH == 1) {
		nrBinsIVH = pt.get("DiscretizationIVH.NrBinsIVH", 64);
		if (nrBinsIVH < 0) {
			std::cout << "please fill in a valid number for the number of bins" << std::endl;
		}
	}
}


/*!
The method getInterpolation reads the interpolation information of the ini-file and
sets the attributes of the class Config to the equivalent values
*/
inline void ConfigFile::getInterpolation() {
	config pt = readIni(fileName);
	interpolation2D = pt.get("Interpolation.2DInterpolation", 0);
	interpolationMethod = pt.get<std::string>("Interpolation.InterpolationMethod");
	useDownSampling = pt.get("Interpolation.UseDownSampling2Cubic", 1);
	if (useDownSampling != 0 && useDownSampling != 1) {
		std::cout << "You inserted a value for UseDownSampling2Cubic that is not 0 or 1, it will be set to 1" << std::endl;
		useDownSampling = 1;
	}
	useUpSampling = pt.get("Interpolation.UseUpSampling2Cubic", 1);
	if (useUpSampling != 0 && useUpSampling != 1) {
		std::cout << "You inserted a value for UseUpSampling2Cubic that is not 0 or 1, it will be set to 1" << std::endl;
		useUpSampling = 1;
	}
	useSampling2mm = pt.get("Interpolation.UseSamplingTo2mmVoxel", 1);
	if (useSampling2mm != 0 && useSampling2mm != 1) {
		std::cout << "You inserted a value for useSampling2mm that is not 0 or 1, it will be set to 1" << std::endl;
		useSampling2mm = 1;
	}
	if (useSampling2mm + useDownSampling + useUpSampling > 1) {
		std::cout << "You chose more than one interpolation method" << std::endl;
		if ( useUpSampling ==1) {
			std::cout << "Only up sampling is performed. Please fill out another config file for the other interpolation methods" << std::endl;
			useSampling2mm = 0;
			useDownSampling = 0;
		}
		else if (useUpSampling != 1 && useDownSampling ==1) {
			std::cout << "Only down sampling is performed. Please fill out another config file for the other interpolation methods" << std::endl;
			useSampling2mm = 0;
			useUpSampling = 0;
		}
	}
}

/*!
The method getDistanceWeightProperties reads the information concerning the distance weights. \n
*/
inline void ConfigFile::getDistanceWeightProperties() {
	config pt = readIni(fileName);
	normGLCM = pt.get<std::string>("DistanceWeightProperties.NormGLCM");
	normGLRLM = pt.get<std::string>("DistanceWeightProperties.NormGLRLM");
	normNGTDM = pt.get<std::string>("DistanceWeightProperties.NormNGTDM");
}


/*!
The method getExtendedEmphasisInformation reads the information concerning the extended emphasis. \n
*/
inline void ConfigFile::getExtendedEmphasisInformation() {
	config pt = readIni(fileName);
	extendedEmphasis = pt.get("ExtendedEmphasisFeatures.CalculateExtendedEmph", 1);
	powerRow = pt.get("ExtendedEmphasisFeatures.PowerRow", 1);
	powerCol = pt.get("ExtendedEmphasisFeatures.PowerCol", 1);
}

inline void ConfigFile::getNGTDMdistanceValue() {
	config pt = readIni(fileName);
	dist = pt.get("NGTDMDistance.dist", 1);
}


inline void ConfigFile::getNGLDMParameters() {
	config pt = readIni(fileName);
	distNGLDM = pt.get("NGLDMParameters.dist", 1);
	coarsenessParam = pt.get("NGLDMParameters.coarseness", 0);
}

/*!
The method getImageFolder gets the image and voi path information the user provided in the command line. \n
Furthermore it reads the image type.
*/
inline void ConfigFile::getImageFolder(string image, string voi) {
	config pt = readIni(fileName);

	imageName = image;
	voiName = voi;
	imageType = pt.get<std::string>("ImageProperties.ImageType");
}

/*!
The method getOutputFolder gets the output folder path provided in the command line and reads the output information
*/
inline void ConfigFile::getOutputInformation(string output) {
	config pt = readIni(fileName);
	outputFolder = output;
	csvOutput = pt.get("OutputInformation.csvOutput", 1);
	ontologyOutput = pt.get("OutputInformation.ontologyOutput", 1);
	getOneCSVFile = pt.get("OutputInformation.GetOneCSVFile", 1);
	if (csvOutput == 0 && getOneCSVFile == 1) {
		std::cout << "The CSVOutput is set to 0, but getOneCSVFile is set to 1. The CSVOutput is still generated" << std::endl;
		csvOutput = 1;
	}
	else if (csvOutput == 0 && ontologyOutput == 0) {
		std::cout << "CSVOutput and ontologyOutput is set to 0. One CSVOutput will be generated" << std::endl;
		csvOutput = 1;
		getOneCSVFile = 1;
	}
}




/*!
The method getPETimageInformation reads the pet image information of the patientInfo.ini-file and
sets the attributes of the class Config to the equivalent values
*/
inline void ConfigFile::getPETimageInformation(string imagePath, string patientInfoPath, ConfigFile con) {
	config pt = readIni(fileName);
	if (pt.get<std::string>("ImageProperties.ImageType") == "PET" && patientInfoPath =="0") {
		std::cout << "The ImageType is set to PET but no patient info was given. Program stops." << std::endl;		
		exit(EXIT_FAILURE);
	}

	if (pt.get<std::string>("ImageProperties.ImageType") == "PET") {
		//read info from patient info ini file
		boost::property_tree::ptree patientInfo;
		boost::property_tree::ini_parser::read_ini(patientInfoLocation, patientInfo);
		useSUV = patientInfo.get("PatientInfo.UseSUV", 1);
		useSUL = patientInfo.get("PatientInfo.UseSUL", 0);
		correctionParam = patientInfo.get<double>("PatientInfo.ScalingFactor");

		if (correctionParam == double(0)) {
			patientWeight = patientInfo.get<float>("PatientInfo.PatientWeight");
			patientHeight = patientInfo.get<float>("PatientInfo.PatientHeight");

			initActivity = patientInfo.get<float>("PatientInfo.ActivityMBq");
			
			gender = patientInfo.get<std::string>("PatientInfo.Gender");
			
			if (gender.compare("M") == 0) {
				malePatient = 1;
			}
			else if (gender.compare("F") == 0) {
				malePatient = 0;
			}
			else {
				std::cout << "The gender could not be identified. The calculation of SUL is not possible" << std::endl;
			}
		}
	}
}

/*!
The method copyConfigFile copies the config file in the output folder
*/
inline void ConfigFile::copyConfigFile(string outputFolder) {
#ifdef _WIN32
	ifstream fin;
	fin.open(fileName, ios::in);
	string iniName = fileName;

	// Remove directory if present.
	// Do this before extension removal incase directory has a period character.
	const size_t last_slash_idx = iniName.find_last_of("\\/");
	if (std::string::npos != last_slash_idx)
	{
		iniName.erase(0, last_slash_idx + 1);
	}
	ofstream fout;
	string outputIniName = outputFolder + "_" + iniName;
	fout.open(outputIniName, ios::out);
	string content = "";
	int i;

	for (i = 0; fin.eof() != true; i++) // get content of infile
		content += fin.get();

	i--;
	content.erase(content.end() - 1);     // erase last character

	fin.close();

	fout << content;                 // output
	fout.close();
#else
	ifstream fin;
	fin.open(fileName.c_str(), ios::in);
	string iniName = fileName;

#endif
}

/*!
The method createConfigInfo fills the class ConfigFile with all information provided in the 
config.ini file. \n
The method just executes all ini-methods.
*/
inline void ConfigFile::createConfigInfo(ConfigFile &config, string arguments[8]) {
	
	config.readIni(arguments[0]);
	config.getThreshold();
	config.patientInfoLocation = arguments[6];
	config.getAccurateState(arguments[4]);
	config.getVoiState(arguments[5]);
	config.getImageFolder(arguments[1], arguments[2]);
	config.getResegmentationState();
	config.getOutputInformation(arguments[3]);
	config.getSmoothingKernel();
	config.getFeatureSelectionLocation(arguments[7]);
	config.getDiscretizationInformation();
	config.getDiscretizationInformationIVH();
	config.getPETimageInformation(arguments[1], arguments[6], config);
	config.getInterpolation();
	config.getDistanceWeightProperties();
	config.getExtendedEmphasisInformation();
	config.getNGLDMParameters();
	config.getNGTDMdistanceValue();
}

/*!
The method createOutputFolder creates the outputFolder, if it does not already exists.\n
If it exists, a warning message is printed on the screen, that the data will be overwritten
*/
inline void ConfigFile::createOutputFolder(ConfigFile &config) {
	string outputFolderTmp = config.outputFolder + string(".csv");
	char * writable = new char[outputFolderTmp.size() + 1];
	std::copy(outputFolderTmp.begin(), outputFolderTmp.end(), writable);
	writable[outputFolderTmp.size()] = '\0';
	char * writable2 = new char[config.outputFolder.size() + 1];
	std::copy(config.outputFolder.begin(), config.outputFolder.end(), writable2);
	writable2[config.outputFolder.size()] = '\0';
	ifstream f(writable);
	if (f.good() == 1) {
		std::time_t now = std::time(NULL);
		std::tm * ptm = std::localtime(&now);
								
		char buffer2[32];
				
		std::strftime(buffer2, 32, "%a%d%m%Y%H%M%S", ptm);
		std::string path = writable2+string(buffer2);
		config.outputFolder = path;
		std::cout << ("The chosen outputfile already exists. The data will be saved in the file %s .csv", path) << std::endl;
	}
}



#endif // READCONFIGFILE_H_INCLUDED
