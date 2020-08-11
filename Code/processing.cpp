//if tumor is too small, fill all values with NANs
//this means that I only write the .csv without calculating the features beforehand
void fillCSVwithNANs(ConfigFile config) {

	
	MorphologicalFeatures<float, 3> morphFeat;

	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		morphFeat.writeCSVFileMorphological(morphFeat, config.outputFolder, config);
	}
	else if (config.csvOutput == 1 && config.getOneCSVFile == 1) {
		morphFeat.writeOneFileMorphological(morphFeat, config);
	}


	LocalIntensityFeatures<float, 3> localIntFeat;

	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		localIntFeat.writeCSVFileLocalIntensity(localIntFeat, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		localIntFeat.writeOneFileLocalInt(localIntFeat, config);
	}

	StatisticalFeatures<float, 3> statFeatures;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		statFeatures.writeCSVFileStatistic(statFeatures, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		statFeatures.writeOneFileStatistic(statFeatures, config);
	}

	IntensityVolumeFeatures<float, 3> intVol;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		intVol.writeCSVFileIntVol(intVol, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		intVol.writeOneFileIntVol(intVol, config);
	}
	IntensityHistogram<float, 3> intensityHist;

	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		intensityHist.writeCSVFileIntensity(intensityHist, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		intensityHist.writeOneFileIntensity(intensityHist, config);
	}

	GLCMFeatures2DAVG<float, 3> glcm2DAVG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glcm2DAVG.writeCSVFileGLCM2DAVG(glcm2DAVG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glcm2DAVG.writeOneFileGLCM2DAVG(glcm2DAVG, config, config.featureParameterSpaceNr);
	}

	GLCMFeatures2DDMRG<float, 3> glcm2DDMRG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glcm2DDMRG.writeCSVFileGLCM2DDMRG(glcm2DDMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glcm2DDMRG.writeOneFileGLCM2DDMRG(glcm2DDMRG, config, config.featureParameterSpaceNr);
	}


	GLCMFeatures2DMRG<float, 3> glcm2DMRG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glcm2DMRG.writeCSVFileGLCM2DMRG(glcm2DMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glcm2DMRG.writeOneFileGLCM2DMRG(glcm2DMRG, config, config.featureParameterSpaceNr);
	}
	GLCMFeatures2DVMRG<float, 3> glcm2DVMRG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glcm2DVMRG.writeCSVFileGLCM2DVMRG(glcm2DVMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glcm2DVMRG.writeOneFileGLCM2DVMRG(glcm2DVMRG, config, config.featureParameterSpaceNr);
	}
	GLCMFeatures3DAVG<float, 3> glcmFeat3DAVGFeat;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glcmFeat3DAVGFeat.writeCSVFileGLCM3DAVG(glcmFeat3DAVGFeat, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glcmFeat3DAVGFeat.writeOneFileGLCM3DAVG(glcmFeat3DAVGFeat, config, config.featureParameterSpaceNr);
	}
	GLCMFeatures3DMRG<float, 3> glcm3DMRG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glcm3DMRG.writeCSVFileGLCM3DMRG(glcm3DMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glcm3DMRG.writeOneFileGLCM3DMRG(glcm3DMRG, config, config.featureParameterSpaceNr);
	}
	GLRLMFeatures2DAVG<float, 3> glrlm2DAVG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glrlm2DAVG.writeCSVFileGLRLM2DAVG(glrlm2DAVG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glrlm2DAVG.writeOneFileGLRLM2DAVG(glrlm2DAVG, config, config.featureParameterSpaceNr);
	}
	GLRLMFEATURES2DDMRG<float, 3> glrlm2DDMRG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glrlm2DDMRG.writeCSVFileGLRLM2DDMRG(glrlm2DDMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glrlm2DDMRG.writeOneFileGLRLM2DDMRG(glrlm2DDMRG, config, config.featureParameterSpaceNr);
	}
	GLRLMFeatures2DMRG<float, 3> glrlm2DMRG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glrlm2DMRG.writeCSVFileGLRLM2DMRG(glrlm2DMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glrlm2DMRG.writeOneFileGLRLM2DMRG(glrlm2DMRG, config, config.featureParameterSpaceNr);
	}
	GLRLMFeatures2DVMRG<float, 3> glrlm2DVMRG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glrlm2DVMRG.writeCSVFileGLRLM2DVMRG(glrlm2DVMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glrlm2DVMRG.writeOneFileGLRLM2DVMRG(glrlm2DVMRG, config, config.featureParameterSpaceNr);
	}
	GLRLMFeatures3DAVG<float, 3> glrlm3DAVG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glrlm3DAVG.writeCSVFileGLRLM3DAVG(glrlm3DAVG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glrlm3DAVG.writeOneFileGLRLM3DAVG(glrlm3DAVG, config, config.featureParameterSpaceNr);
	}
	GLRLMFeatures3D<float, 3> glrlm3DMRG;

	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glrlm3DMRG.writeCSVFileGLRLM3D(glrlm3DMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glrlm3DMRG.writeOneFileGLRLM3D(glrlm3DMRG, config, config.featureParameterSpaceNr);
	}
	GLSZMFeatures2DAVG<float, 3> glszm2DAVG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glszm2DAVG.writeCSVFileGLSZM2DAVG(glszm2DAVG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glszm2DAVG.writeOneFileGLSZM2DAVG(glszm2DAVG, config, config.featureParameterSpaceNr);
	}
	GLSZMFeatures2DMRG<float, 3> glszm2D;

	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glszm2D.writeCSVFileGLSZM(glszm2D, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glszm2D.writeOneFileGLSZM(glszm2D, config, config.featureParameterSpaceNr);
	}
	GLSZMFeatures3D<float, 3> glszm3D;

	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		glszm3D.writeCSVFileGLSZM3D(glszm3D, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		glszm3D.writeOneFileGLSZM3D(glszm3D, config, config.featureParameterSpaceNr);
	}

	NGTDM2DAVG<float, 3> NGTDM2DAVG;

	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		NGTDM2DAVG.writeCSVFileNGTDM2DAVG(NGTDM2DAVG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		NGTDM2DAVG.writeOneFileNGTDM2DAVG(NGTDM2DAVG, config, config.featureParameterSpaceNr);
	}

	NGTDMFeatures2DMRG<float, 3> ngtdm2DMRG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		ngtdm2DMRG.writeCSVFileNGTDM(ngtdm2DMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		ngtdm2DMRG.writeOneFileNGTDM(ngtdm2DMRG, config, config.featureParameterSpaceNr);
	}

	NGTDMFeatures3D<float, 3> ngtdm3D;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		ngtdm3D.writeCSVFileNGTDM3D(ngtdm3D, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		ngtdm3D.writeOneFileNGTDM3D(ngtdm3D, config, config.featureParameterSpaceNr);
	}
	GLDZMFeatures2DAVG<float, 3> gldzm2DAVG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		gldzm2DAVG.writeCSVFileGLDZM2DAVG(gldzm2DAVG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		gldzm2DAVG.writeOneFileGLDZM2DAVG(gldzm2DAVG, config, config.featureParameterSpaceNr);
	}

	GLDZMFeatures2D<float, 3> gldzm2D;

	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		gldzm2D.writeCSVFileGLDZM(gldzm2D, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		gldzm2D.writeOneFileGLDZM(gldzm2D, config, config.featureParameterSpaceNr);
	}

	GLDZMFeatures3D<float, 3> gldzm3D;

	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		gldzm3D.writeCSVFileGLDZM3D(gldzm3D, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		gldzm3D.writeOneFileGLDZM3D(gldzm3D, config, config.featureParameterSpaceNr);
	}
	NGLDMFeatures2DAVG<float, 3> ngldm2DAVG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		ngldm2DAVG.writeCSVFileNGLDM2DAVG(ngldm2DAVG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		ngldm2DAVG.writeOneFileNGLDM2DAVG(ngldm2DAVG, config, config.featureParameterSpaceNr);
	}

	NGLDMFeatures2DMRG<float, 3> ngldm2DMRG;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		ngldm2DMRG.writeCSVFileNGLDM2DMRG(ngldm2DMRG, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		ngldm2DMRG.writeOneFileNGLDM2DMRG(ngldm2DMRG, config, config.featureParameterSpaceNr);
	}

	NGLDMFeatures3D<float, 3> ngldm3D;
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		ngldm3D.writeCSVFileNGLDM3D(ngldm3D, config.outputFolder);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		ngldm3D.writeOneFileNGLDM3D(ngldm3D, config, config.featureParameterSpaceNr);
	}

}

ImageType::Pointer thresholdMask(ImageType *mask, float threshold) {
	itk::ImageRegionIterator<ImageType> countVoxels(mask, mask->GetLargestPossibleRegion());
	countVoxels.GoToBegin();

	float maxValue = 0;
	while (!countVoxels.IsAtEnd()){
		if (countVoxels.Get() > maxValue){
			maxValue = countVoxels.Get();

		}
		++countVoxels;
	}
	countVoxels.GoToBegin();
	while (!countVoxels.IsAtEnd()){
		if (countVoxels.Get() >= threshold * maxValue){
			countVoxels.Set(1);
		}
		else {
			countVoxels.Set(0);
		}
		++countVoxels;
	}
	return mask;
}

//write PET metrics to .csv
void writePETmetrics(float value, string nameVariable, ConfigFile config) {
	string csvName;

	if (config.getOneCSVFile == 1) {
		csvName = config.outputFolder + ".csv";
	}
	else if (config.ontologyOutput == 1) {
		csvName = config.outputFolder + "/feature_table.csv";
	}
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';
	ofstream morphCSV;


	morphCSV.open(name, std::ios_base::app);

	if (config.getOneCSVFile == 1) {
		if (config.imageType == "PET") {
			morphCSV << "PET Uptake Metrics" << "," << nameVariable << ", ";
		}
		else {
			morphCSV << "Exact Metrics" << "," << nameVariable << ", ";
		}

		morphCSV << value;
		morphCSV << "\n";


		morphCSV.close();
	}

}

//calculate the original volume before interpolation
float getOriginalVolume(ImageType::Pointer mask) {
	itk::ImageRegionConstIterator<ImageType> countVoxels(mask, mask->GetLargestPossibleRegion());
	countVoxels.GoToBegin();
	float maxValue = 0;
	while (!countVoxels.IsAtEnd())
	{
		if (countVoxels.Get() > maxValue)
		{
			maxValue = countVoxels.Get();
		}
		++countVoxels;
	}
	float volume = 0;
	countVoxels.GoToBegin();
	while (!countVoxels.IsAtEnd())
	{
		if (countVoxels.Get() >= 0.5*maxValue)
		{
			volume += 1;
		}
		++countVoxels;
	}
	return volume;
}


void calculatePETmetrics(ImageType::Pointer image, ImageType::Pointer mask, int volume, ConfigFile config) {
	//upsampled or downsampled include PET metrics
	RegionType boundingBoxRegionTmp = getBoundingBoxMask(mask);
	ImageType::Pointer tmpmaskFiltered = getImageMasked(mask, boundingBoxRegionTmp);
	ImageType::Pointer tmpimageFiltered = getImageMasked(image, boundingBoxRegionTmp);
	tmpmaskFiltered = thresholdMask(tmpmaskFiltered, config.threshold);
	const typename ImageType::RegionType& imageRegion = tmpimageFiltered->GetLargestPossibleRegion();
	const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
	Image<float, 3> imageTestTmp(imageRegionSize[0], imageRegionSize[1], imageRegionSize[2]);
	imageTestTmp.getImageAttributes(tmpimageFiltered, tmpmaskFiltered, config);
	LocalIntensityFeatures<float, 3> localInt;
	localInt.calculateAllLocalIntensityFeatures(localInt, tmpimageFiltered, tmpmaskFiltered, config);
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		localInt.writeCSVFileLocalIntensityPET(localInt, config);
	}
	else if ((config.csvOutput == 1 && config.getOneCSVFile == 1) || config.ontologyOutput == 1) {
		localInt.writeOneFileLocalIntPET(localInt, config);
	}
	std::string forLog = "Local intensity features without rebinning were calculated.";
	typedef boost::accumulators::features <
		tag::mean> Features;
	typedef accumulator_set <float, Features> Accumulator;
	Accumulator acc;
	for_each(imageTestTmp.vectorOfMatrixElements.begin(), imageTestTmp.vectorOfMatrixElements.end(), boost::bind<void>(boost::ref(acc), _1));
	float  maximumValue = *max_element(imageTestTmp.vectorOfMatrixElements.begin(), imageTestTmp.vectorOfMatrixElements.end());
	int length;
	float meanValue = mean(acc);
	float TLG = volume * meanValue;

	writePETmetrics(maximumValue, "Original max", config);
	writePETmetrics(meanValue, "Original mean", config);
	writePETmetrics(TLG, "Original TLG", config);
	tmpmaskFiltered = nullptr;
	tmpimageFiltered = nullptr;
	Image<float, 3> imageTestTmp2(0, 0, 0);
	imageTestTmp = imageTestTmp2;
}


void storePreInterpolationFeatures(ImageType::Pointer image, ImageType::Pointer maskImage, ConfigFile config) {
	//do the PET uptake metrics have to be calculated?
	const typename ImageType::SpacingType& inputSpacing = image->GetSpacing();
	std::string a = "1";
	int rebinning;
	if (config.calculateAllFeatures == 0) {
		boost::property_tree::ptree pt;
		boost::property_tree::ini_parser::read_ini(config.featureSelectionLocation, pt);
		rebinning = pt.get("ExactMetrics.CalculateExactMetrics", 0);
	}

	//calculate PET Uptake metrics if required

	if (config.calculateAllFeatures == 1 || rebinning == 1) {
		float volume = getOriginalVolume(maskImage);
		volume = volume * inputSpacing[0] * inputSpacing[1] * inputSpacing[2];
		calculatePETmetrics(image, maskImage, volume, config);
		writeExactVolume(volume, config);
	}
	std::cout << "Exact metrics calculated" << std::endl;
	//calculate dispersity features
	DispersityFeatures<float, 3> disp;

	int dispersity;
	//do the dispersity features have to be calculated?
	if (config.calculateAllFeatures == 0) {
		boost::property_tree::ptree ptDisp;

		boost::property_tree::ini_parser::read_ini(config.featureSelectionLocation, ptDisp);
		dispersity = ptDisp.get("DispersityFeatures.CalculateDispersityFeat", 1);
		std::cout << "DISP" << dispersity << std::endl;
	}

	if (dispersity == 1 || config.calculateAllFeatures == 1) {
		disp.calculateAllDispersityFeatures(disp, image, maskImage, config);
		string forLog = "Dispersity features were calculated.";
		writeLogFile(config.outputFolder, forLog);
		std::cout << "Dispersity features are calculated" << std::endl;
	}
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		disp.writeCSVFileDispersity(disp, config.outputFolder, config);
	}
	else if (config.csvOutput == 1 && config.getOneCSVFile == 1) {
		disp.writeOneFileDispersity(disp, config);
	}
}
//wirte exact volume (part of exact features) to csv
void writeExactVolume(float volume, ConfigFile config) {
	string csvName;
	if (config.getOneCSVFile == 1) {
		csvName = config.outputFolder + ".csv";
	}
	else if (config.ontologyOutput == 1) {
		csvName = config.outputFolder + "/feature_table.csv";
	}
	char * name = new char[csvName.size() + 1];
	std::copy(csvName.begin(), csvName.end(), name);
	name[csvName.size()] = '\0';
	ofstream morphCSV;

	morphCSV.open(name, std::ios_base::app);

	if (config.getOneCSVFile == 1) {
		if (config.imageType == "PET") {
			morphCSV << "PET Uptake Metrics" << "," << "ExactVolume,";
		}
		else {
			morphCSV << "Exact Metrics" << "," << "ExactVolume,";
		}
		morphCSV << volume;
		morphCSV << "\n";

		morphCSV.close();
	}


}