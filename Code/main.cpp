#include <iostream>
#include <stdio.h>
using namespace std;

#include "softwareParameters.h"
#include "featureCalculation.h"

#include <filesystem>

int main (int argc, char* argv[]) {

/*! \mainpage RaCaT
Welcome to the documentation of RaCaT! \n
\n
RaCaT can be used with any kind of images (like CT, PET and MRI) in dicom, ecat, nrrd or nifti format.
Masks marking the region of interest can be imported as binary masks (also in dicom, ecat, nrrd or nifti format) or as rt-struct.\n
RaCat comes as an executable and does not require any further installations. As input it needs a configuration file where you can set the preprocessing steps (like discretization, resegmentation etc) you want to use.
For more explanation see the next sections.
\section steps_sec Getting started
In order to use the Radiomics Calculator, the Radiomics.exe has to be downloaded. In the folder ExampleFiles, you can find examples of the configuration files and other files you can give to the tool (optional). Every file
is explained in detail in the following section.
\arg Download and adapt config.ini \n
In the config.ini file you can set several preprocessing parameters. You can set e.g. if you want to use a discretization step (or not), if your image is a PET or CT image, if you want to apply resegmentation etc. More information about the config file can be found in the section "Setting up the config.ini file".  \n
\arg Optional: Download and adapt featureSelection.ini \n
RaCaT calculates a large number of feature values. The features are ordered in several groups, e.g. Statistical features, Grey-level-cooccurrence features etc. If you want to calculate only certain features, you can say this in
the featureSelection.ini file. This file is not required to make RaCaT work. So you only have to change and adapt it if you want to calculate only certain features. Then you also have to give the path where you saved your featureSelection file as parameter
to the executable. In the ExampleFiles folder, there are some examples for feature selection files and how you can call the tool if you want to include a feature selection file. If no featureSelection file is given, all features are calculated. \n
\arg Only for PET images: Download and adapt patientInfo.ini \n
If your image is a PET image, RaCaT needs the patient information so it can convert the intensity values in the image from Bq/ml to SUV. For this, you can download the patientInfo.ini file and fill in your patient parameters like weight etc.  \n

\arg Call the executable
After adapting these files, you can call the executable. As parameters it needs:
- path to the ini-file \n
- path to the image \n
- path to the mask \n
- path to folder where you want to save the images (if the folder does not exist, the folder is created) \n

A lot of examples on how to call the executable for PET/CT images, with or without featureSelection.ini can be found in the example folder.
IF the mask is a binary image in dicom/nifti/nrrd format, the executable can be called in the following way: \n

\code{.unparsed}
/path/to/.exe --ini /path/to/iniFile --img /path/to/prj/image file  --voi /path/to/voi(mask) file  --out /path/to/outputfolder
\endcode
\n

The output folder will be created automatically. If the folder already exists, a new folder will be created which has the name of the required output folder including the time stamp of the calculation. 
If the mask is a RT-struct, --voi has to be replaced by --rts:
\code{.unparsed}
/path/to/.exe --ini /path/to/iniFile --img /path/to/prj/image file  --rts /path/to/rt struct  --out /path/to/outputfolder
\endcode
\arg Help
With --h, a help is called, that gives an example of how the executable should be called.\n
\arg Example data
Example files for config, feature selection, and patient information are stored in the folder 'ExampleFiles'. Also batch-files with examples how to call the executable can be found in this folder. 
*/
    //read in the config file
    ConfigFile config;
	//the arguments given by the user are saved in a vector in a specific order
	string* arguments = new string[8] {"0", "0","0","0", "0", "0","0", "0"};
	string* nameArgument = new string[8]{ "ini", "img", "voi", "out" ,"acc", "voi", "pat", "fts"};
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		//get ini file path and store it in vector
		if (arg == "--ini") {
			i++;
			arguments[0]=std::string(argv[i]);
		}
		//get path to image and store it
		else if (arg == "--img") {
			i++;
			arguments[1] = std::string(argv[i]);
			if (arguments[1].substr(arguments[1].length() - 3) == "prj") {
				arguments[4] = "acc";
			}
			//if directory, its a dicomFile
			struct stat s;
			if (stat(arguments[1].c_str(), &s) == 0)
			{
				if (s.st_mode & S_IFDIR)
				{
					arguments[4] = "dic";
				}
			}
		}
		//get path to voi 
		else if (arg == "--voi") {
			i++;
			arguments[2]=std::string(argv[i]);
			//if it is a directory, also the voi is in dicom format
			struct stat s;
			if (stat(arguments[2].c_str(), &s) == 0)
			{
				if (s.st_mode & S_IFDIR)
				{
					arguments[5] = "dic";
				}
			}
		}
		//I let the user define if he has a rt struct as voi
		else if (arg == "--rts") {
			i++;
			arguments[2] = std::string(argv[i]);

			arguments[5] = "rts";
			
		}
		//get the path of the desired output folder
		else if (arg == "--out") {
			i++;
			arguments[3]=string(argv[i]);
		}
		//get path to patient info
		else if (arg == "--pat") {
			i++;
			arguments[6] = string(argv[i]);
		}
		//get path to feature selection
		else if (arg == "--fod") {
			i++;
			std::cout << "feat" << std::endl;
			arguments[7] = string(argv[i]);
			
		}
		//call the help
		else if (arg == "--h") {
			arguments[0] = "--h";
			i++;
		}

		
	}
	
	for (int i = 0; i < 4; i++) {
		if (arguments[i] == "0" && arguments[0]!="--h") {
			std::cout << "You forgot to set the input parameter " + nameArgument[i] << std::endl;
			std::cout << "Please check if you used two -: e.g. --ini/--out to call the function" << std::endl;
			return 0;
		}
		else if (arguments[0] == "--h") {
				std::cout << "help" << std::endl;
				std::cout << "The following parameters have to be set (IF mask is a binary image): \n --ini: location of ini file  \n --pat: location of patient info file (obligatory for PET images, for other images not required)"  << std::endl;
				std::cout << "--fod location of featureOutputDefinition file (optional, if not set, all feature values are calculated." << std::endl;
				std::cout << "--img: location of image or project file \n --voi: location of mask/voi file" << std::endl;
				std::cout << "--out: location of desired output folder" << std::endl;
				std::cout << "IF mask is RT struct, --voi has to be replaced by --rts: \n --ini: location of ini file  \n --pat: location of patient info file" << std::endl;
				std::cout << "--img: location of image or project file \n --rts: location of rt struct" << std::endl;

				return 0;
			
		}

	}

	if (arguments[6] == "0") {
		std::cout << "No patientInfo location has been set." << std::endl;
	}
	if (arguments[7] != "0") {
		std::cout << "A feature selection file has been set, only these features will be calculated." << std::endl;
	}
    config.createConfigInfo(config, arguments);
	string originalOutputname = config.outputFolder;
	//check the parameters set by the user
	//dependent on them, the name of the outpufolder is modified
	if (config.useFixedBinWidth == 1 && config.useFixedNrBins == 1 ) {
		config.useFixedBinWidth = 0;
		//change output folder name
		config.outputFolder = originalOutputname + "FXDBin";
		if (config.csvOutput == 1) {
			config.createOutputFile(config);
		}
		prepareDataForFeatureCalculation(config);
		config.copyConfigFile(config.outputFolder);
		
		config.useFixedBinWidth = 1;
		config.useFixedNrBins = 0;
		config.outputFolder = originalOutputname + "FXDWidth";
		if (config.csvOutput == 1) {
			config.createOutputFile(config);
		}
		prepareDataForFeatureCalculation(config);
		config.copyConfigFile(config.outputFolder);
		
	}
	else{
		if (config.csvOutput == 1) {
			config.createOutputFile(config);
		}
		prepareDataForFeatureCalculation(config);
		config.copyConfigFile(config.outputFolder);
	}
	return 0;
}


