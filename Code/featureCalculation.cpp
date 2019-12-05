ImageType::Pointer thresholdMask(ImageType *mask, float threshold) {
	itk::ImageRegionConstIterator<ImageType> countVoxels(mask, mask->GetLargestPossibleRegion());
	countVoxels.GoToBegin();
	float maxValue = 0;
	while (!countVoxels.IsAtEnd())
	{
		if (countVoxels.Get() > maxValue)
		{
			maxValue = countVoxels.Get();
			//std::cout << maxValue << std::endl;
		}
		++countVoxels;
		
	}
	using PixelType = float;
	const auto LowerThreshold = static_cast<PixelType>((threshold)*maxValue);
	using FilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(mask);
	filter->SetLowerThreshold(LowerThreshold);
	filter->SetOutsideValue(0);
	filter->SetInsideValue(100);
	filter->Update();
	ImageType::Pointer updatedMask = filter->GetOutput();
	return updatedMask;


}



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
	if (config.imageType == "PET") {
		morphCSV.open(name, std::ios_base::app);
	}
	else {
		morphCSV.open(name);
	}

	if (config.getOneCSVFile == 1) {

		morphCSV << "PET Uptake Metrics" << "," << "ExactVolume,";
		morphCSV << volume;
		morphCSV << "\n";

		morphCSV.close();
	}


}

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

		morphCSV << "PET Uptake metrics" << "," << nameVariable<<", ";
		morphCSV << value;
		morphCSV << "\n";
		

		morphCSV.close();
	}

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
	localInt.calculateAllLocalIntensityFeatures(localInt, imageTestTmp, config);
	if (config.csvOutput == 1 && config.getOneCSVFile == 0) {
		localInt.writeCSVFileLocalIntensityPET(localInt, config.outputFolder);
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
	float meanValue = std::accumulate(imageTestTmp.vectorOfMatrixElements.begin(), imageTestTmp.vectorOfMatrixElements.end(), 0);
	int length;

	meanValue = meanValue / (imageTestTmp.vectorOfMatrixElements.size());
	//float meanValue = mean(acc);
	float TLG = volume * meanValue;

	writePETmetrics(maximumValue, "Original max", config);
	writePETmetrics(meanValue, "Original mean", config);
	writePETmetrics(TLG, "Original TLG", config);

}
//in this function the .voi file of the accurate project is read in
//First, the image and voxel size is read in from the .prj file
//then the .voi file is read in and transformed to an ITK image

ImageType::Pointer readVoiFilePET(string prjPath, string voiPath, ImageType *image, ConfigFile config, unsigned int(&dimPET)[3], float voxelSize[3]) {
	ImageType::Pointer maskImage;
	//parameters to store the dimension of the PET and CT image
	//the voxel size, halfLife of the tracer, the volume scale and the frame scale
	int nrVoxelsPET = int(dimPET[0] * dimPET[1] * dimPET[2]);
	//prjfile.close();
	//open the voi file
	ifstream voiFile;
	voiFile.open(voiPath, ios::in | ios::binary);
	if (voiFile.is_open() == false) {
		std::cout << "Cannot open voi file\n";
		exit(0);
	}

	else {
		//read the voi file
		streampos begin, end;
		begin = voiFile.tellg();
		voiFile.seekg(0, ios::end);
		end = voiFile.tellg();
		int fileSize = (end - begin);
		//vector to store all the elements
		vector<char> wholeFile(fileSize);

		voiFile.seekg(0, std::ios::beg);
		voiFile.read(&wholeFile[0], fileSize);

		vector< float> mask(fileSize);

		transform(wholeFile.begin(), wholeFile.end(), mask.begin(), [](auto& elem) {return ((float)(elem)); });

		vector< float> voi(wholeFile.begin() + (fileSize / 2), wholeFile.end());
		

		float *voiArray = &voi[0];

		maskImage = converArray2Image(voiArray, dimPET, voxelSize);
		using FilterType = itk::ChangeInformationImageFilter<ImageType>;
		FilterType::Pointer filter = FilterType::New();
		maskImage->Update();
		voiFile.close();
		
		//filter mask and image 
		ImageType::Pointer maskFiltered= maskImage;
		ImageType::Pointer imageFiltered= image;
		
		//do the PET uptake metrics have to be calculated?
		const typename ImageType::SpacingType& inputSpacing = image->GetSpacing();
		std::string a = "1";
		string rebinning;
		if (config.calculateAllFeatures == 0) {
			boost::property_tree::ptree pt;
			boost::property_tree::ini_parser::read_ini(config.featureSelectionLocation, pt);
			rebinning = pt.get<std::string>("PETUptakeMetrics.CalculatePETUptakeMetrics");
		}
		//calculate PET Uptake metrics if required
		
		if (config.imageType == "PET" && (config.calculateAllFeatures == 1 || rebinning == a)) {
			float volume = getOriginalVolume(maskImage);
			volume = volume * inputSpacing[0] * inputSpacing[1] * inputSpacing[2];
			calculatePETmetrics(image, maskImage, volume, config);
			writeExactVolume(volume, config);
		}
		//up or downsample the image
		if (config.useDownSampling == 1 || config.useUpSampling == 1 && (inputSpacing[0] != inputSpacing[1] || inputSpacing[0] != inputSpacing[2])) {

			const typename ImageType::RegionType& imageRegion = image->GetLargestPossibleRegion();
			const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
			double outputSpacing[3];
			vector<int> newImageSize = getImageSizeInterpolated(image, imageRegionSize, outputSpacing, config);

			itk::Size<3> outputSize = { { newImageSize[0], newImageSize[1], newImageSize[2] } };
			Image<float, 3> imageMask(newImageSize[0], newImageSize[1], newImageSize[2]);
			Image<float, 3> imageTest(newImageSize[0], newImageSize[1], newImageSize[2]);

			imageFiltered = imageTest.getResampledImage(imageFiltered, outputSpacing, outputSize, config.interpolationMethod, config.rebinning_centering);
			maskFiltered = imageTest.getResampledImage(maskFiltered, outputSpacing, outputSize, "linear", config.rebinning_centering);
			imageTest.createOntologyVoxelDimensionTable(config, voxelSize);

		}

		else if (config.useSamplingCubic == 1) {

			const typename ImageType::RegionType& imageRegion = image->GetLargestPossibleRegion();
			const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
			double outputSpacing[3];
			vector<int> newImageSize = getImageSizeInterpolated(image, imageRegionSize, outputSpacing, config);
			itk::Size<3> outputSize = { { newImageSize[0], newImageSize[1], newImageSize[2] } };
			Image<float, 3> imageMask(newImageSize[0], newImageSize[1], newImageSize[2]);
			Image<float, 3> imageTest(newImageSize[0], newImageSize[1], newImageSize[2]);
			imageFiltered = imageTest.getResampledImage(imageFiltered, outputSpacing, outputSize, config.interpolationMethod, config.rebinning_centering);
			maskFiltered = imageTest.getResampledImage(maskFiltered, outputSpacing, outputSize, "Linear", config.rebinning_centering);
			imageTest.getValueInMask(maskFiltered);
			imageTest.createOntologyVoxelDimensionTable(config, voxelSize);


		}
		
		maskFiltered = thresholdMask(maskFiltered, config.threshold);
		
		//get the region of the mask
		RegionType boundingBoxRegion = getBoundingBoxMask(maskFiltered);
		itk::Size<3> regionSize = boundingBoxRegion.GetSize();
		//shrink image and mask to the mask region
		imageFiltered = getImageMasked(imageFiltered, boundingBoxRegion);
		maskFiltered = getImageMasked(maskFiltered, boundingBoxRegion);
		
		//change spacing of mask to spacing of image (in order to calculate the morphological features)
		ImageType::Pointer maskNewSpacing = getMaskNewSpacing(imageFiltered, maskFiltered);

		std::cout << "bbsize" << regionSize[0] << " " << regionSize[1] << " " << regionSize[2] << std::endl;
		//shrink image and mask to the mask region
		//ImageType::Pointer imageFiltered = getImageMasked(image, boundingBoxRegion);
		//ImageType::Pointer maskFiltered = getImageMasked(mask, boundingBoxRegion);
		const typename ImageType::SpacingType& spacingVoxelDim = imageFiltered->GetSpacing();
		float voxelSize[3];
		for (int i = 0; i < 3; i++) {
			voxelSize[i] = (float)spacingVoxelDim[i];
		}
		Image<float, 3> imageVoxelDim(10, 10, 10);
		imageVoxelDim.createOntologyVoxelDimensionTable(config, voxelSize);
		//change spacing of mask to spacing of image (in order to calculate the morphological features)

		calculateFeaturesForConfig(imageFiltered, maskFiltered, config);

	}

	return maskImage;
}


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

void prepareDataForFeatureCalculation(ConfigFile config) {
	/*!
	In the function prepareDataForFeatureCalculation, first the image and the mask are read. For this, the ITK-library is used. \n
	After reading the mask, a bounding box from the region of interest is created. \n
	The region of this bounding box is extracted from the image and the mask, which leads to smaller subimages.
	From these subimages, image attributes are extracted.
	*/
	ImageType::Pointer image;
	ImageType::Pointer mask;
	writeImageData2Log(config);
	//read image and mask
	if (config.useAccurate == 0 || config.useAccurate == 2) {
		if (config.useAccurate == 0) {
			image = readImage(config.imageName);
			std::cout << "The input image is a nifti image" << std::endl;
			typename ImageType::DirectionType direction = image->GetDirection();
			using FlipImageFilterType = itk::FlipImageFilter<ImageType>;
			FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
			flipFilter->SetInput(image);
			FlipImageFilterType::FlipAxesArrayType flipAxes;
			flipAxes[0] = false;
			flipAxes[1] = false;
			flipAxes[2] = false;
			if (direction[0][0] == -1) {
				flipAxes[0] = true;
				
			}
			if (direction[1][1] == -1) {
				flipAxes[1] = true;
			}
			if (direction[2][2] == -1) {
				flipAxes[2] = true;
			}
			flipFilter->SetFlipAxes(flipAxes);
			flipFilter->Update();
			image = flipFilter->GetOutput();
			using FilterTypeInfo = itk::ChangeInformationImageFilter<ImageType>;
			FilterTypeInfo::Pointer filterInfo = FilterTypeInfo::New();
			filterInfo->SetInput(image);
			filterInfo->ChangeDirectionOn();
			//direction matrix in compliance with standard nifti files from UMCG
			using MatrixType = itk::Matrix<double, 3, 3>;
			MatrixType matrix;
			matrix[0][0] = 1;
			matrix[0][1] = 0;
			matrix[0][2] = 0.;

			matrix[1][0] = 0;
			matrix[1][1] = 1;
			matrix[1][2] = 0.;

			matrix[2][0] = 0;
			matrix[2][1] = 0;
			matrix[2][2] = 1;
			filterInfo->SetOutputDirection(matrix);
			filterInfo->Update();
			image = filterInfo->GetOutput();
		}
		else {
			image = readDicom(config.imageName);
		}
		if (config.voiFile == 2) {
			mask = readDicom(config.voiName);
		}
		else if (config.voiFile == 3) {
			mask = readRTstruct(config.voiName, image, mask);
			//I overwrite the image.
			//so I read it in again... not the best solution - check again
			if (config.useAccurate == 0) {
				image = readImage(config.imageName);
			}
			else {
				image = readDicom(config.imageName);
			}

		}
		else {
			mask = readImage(config.voiName);
			typename ImageType::DirectionType direction = mask->GetDirection();
			using FlipImageFilterType = itk::FlipImageFilter<ImageType>;
			FlipImageFilterType::Pointer flipFilterMask = FlipImageFilterType::New();
			flipFilterMask->SetInput(mask);
			FlipImageFilterType::FlipAxesArrayType flipAxesMask;
			flipAxesMask[0] = false;
			flipAxesMask[1] = false;
			flipAxesMask[2] = false;
			if (direction[0][0] == -1) {
				flipAxesMask[0] = true;

			}
			if (direction[1][1] == -1) {
				flipAxesMask[1] = true;
			}
			if (direction[2][2] == -1) {
				flipAxesMask[2] = true;
			}
			flipFilterMask->SetFlipAxes(flipAxesMask);
			flipFilterMask->Update();
			mask = flipFilterMask->GetOutput();
			mask->Update();
			using FilterTypeInfo = itk::ChangeInformationImageFilter<ImageType>;
			FilterTypeInfo::Pointer filterInfo = FilterTypeInfo::New();
			filterInfo->SetInput(mask);
			filterInfo->ChangeDirectionOn();
			//direction matrix in compliance with standard nifti files from UMCG
			using MatrixType = itk::Matrix<double, 3, 3>;
			MatrixType matrix;
			matrix[0][0] = 1;
			matrix[0][1] = 0;
			matrix[0][2] = 0.;

			matrix[1][0] = 0;
			matrix[1][1] = 1;
			matrix[1][2] = 0.;

			matrix[2][0] = 0;
			matrix[2][1] = 0;
			matrix[2][2] = 1;
			filterInfo->SetOutputDirection(matrix);
			filterInfo->Update();
			mask = filterInfo->GetOutput();

		}
		if (config.smoothingKernel > 0) {
			image = smoothImage(image, config.smoothingKernel);
		}
		
		
		//resample image if required
		//get the bounding box of the image voi
		//RegionType boundingBoxRegion = getBoundingBoxMask(mask);
		//filter mask and image 
		mask = thresholdMask(mask, config.threshold);

		ImageType::Pointer maskFiltered = mask;// = getImageMasked(mask, boundingBoxRegion);// = maskImage; //// = maskImage;//getImageMasked(maskImage, boundingBoxRegion);
		ImageType::Pointer imageFiltered = image;// = getImageMasked(image, boundingBoxRegion);// = image;// = image;//getImageMasked(image, boundingBoxRegion);
		const ImageType::DirectionType direction = mask->GetDirection();
		ImageType::ConstPointer output = mask;
		const typename ImageType::SpacingType& inputSpacing = image->GetSpacing();
		std::cout << "original voxel size" << inputSpacing[0] << " " << inputSpacing[1] << " " << inputSpacing[2] << std::endl;
		std::string a = "1";
		string rebinning;
		//check if PET uptake metrics should be calculated
		if (config.calculateAllFeatures == 0) {
			boost::property_tree::ptree pt;


			boost::property_tree::ini_parser::read_ini(config.featureSelectionLocation, pt);
			rebinning = pt.get<std::string>("PETUptakeMetrics.CalculatePETUptakeMetrics");
		}
		//calculate PET uptake metrics if required
		if (config.imageType == "PET" && (config.calculateAllFeatures == 1 || rebinning == a)) {
			float volume = getOriginalVolume(mask);
			volume = volume * inputSpacing[0] * inputSpacing[1] * inputSpacing[2];
			calculatePETmetrics(image, mask, volume, config);
			writeExactVolume(volume, config);
		}

		//now down or upsample the image
		if ((config.useDownSampling != 0 || config.useUpSampling != 0) && (inputSpacing[0] != inputSpacing[1] || inputSpacing[0] != inputSpacing[2])) {
			const typename ImageType::RegionType& imageRegion = image->GetLargestPossibleRegion();
			const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();

			double outputSpacing[3];
			vector<int> newImageSize = getImageSizeInterpolated(image, imageRegionSize, outputSpacing, config);
			itk::Size<3> outputSize = { { newImageSize[0], newImageSize[1], newImageSize[2] } };
			Image<float, 3> imageMask(newImageSize[0], newImageSize[1], newImageSize[2]);
			Image<float, 3> imageTest(newImageSize[0], newImageSize[1], newImageSize[2]);
			imageFiltered = imageTest.getResampledImage(imageFiltered, outputSpacing, outputSize, config.interpolationMethod, config.rebinning_centering);
			maskFiltered = imageTest.getResampledImage(maskFiltered, outputSpacing, outputSize, "Linear", config.rebinning_centering);
		}

		else if (config.useSamplingCubic == 1 ) {
			const typename ImageType::RegionType& imageRegion = image->GetLargestPossibleRegion();
			const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
			std::cout << "original image size" << imageRegionSize[0] << " " << imageRegionSize[1] << " " << imageRegionSize[2] << std::endl;
			double outputSpacing[3];
			vector<int> newImageSize = getImageSizeInterpolated(image, imageRegionSize, outputSpacing, config);
			itk::Size<3> outputSize = { { double(newImageSize[0]), double(newImageSize[1]), double(newImageSize[2]) } };
			Image<float, 3> imageMask(newImageSize[0], newImageSize[1], newImageSize[2]);
			std::cout << "after interpolation" << newImageSize[0] << " " << newImageSize[1] << " " << newImageSize[2] << std::endl;
			std::cout << outputSpacing[0] << " " << outputSpacing[1] << " " << outputSpacing[2] << std::endl;
			//Image<float, 3> imageTest(newImageSize[0], newImageSize[1], newImageSize[2]);
			imageFiltered = imageMask.getResampledImage(imageFiltered, outputSpacing, outputSize, config.interpolationMethod, config.rebinning_centering);
			maskFiltered = imageMask.getResampledImage(maskFiltered, outputSpacing, outputSize, "Linear", config.rebinning_centering);
			imageMask.getValueInMask(maskFiltered);
			std::cout << "afterInterpo" << std::endl;
		}


		//convert mask values to 1 (necessary after interpolation)
		maskFiltered = thresholdMask(maskFiltered, config.threshold);
		
		//get the region of the mask
		RegionType boundingBoxRegion = getBoundingBoxMask(maskFiltered);
		itk::Size<3> regionSize = boundingBoxRegion.GetSize();
		//shrink image and mask to the mask region
		imageFiltered = getImageMasked(imageFiltered, boundingBoxRegion);
		maskFiltered = getImageMasked(maskFiltered, boundingBoxRegion);
		//change spacing of mask to spacing of image (in order to calculate the morphological features)
		ImageType::Pointer maskNewSpacing = getMaskNewSpacing(imageFiltered, maskFiltered);
		
		std::cout << "bbsize" << regionSize[0] << " " << regionSize[1] <<" "<<  regionSize[2] << std::endl;
		//shrink image and mask to the mask region
		//ImageType::Pointer imageFiltered = getImageMasked(image, boundingBoxRegion);
		//ImageType::Pointer maskFiltered = getImageMasked(mask, boundingBoxRegion);
		const typename ImageType::SpacingType& spacingVoxelDim = imageFiltered->GetSpacing();
		float voxelSize[3];
		for (int i = 0; i < 3; i++) {
			voxelSize[i] = (float)spacingVoxelDim[i];
		}
		Image<float, 3> imageVoxelDim(10, 10, 10);
		imageVoxelDim.createOntologyVoxelDimensionTable(config, voxelSize);
		//change spacing of mask to spacing of image (in order to calculate the morphological features)
		
		calculateFeaturesForConfig(imageFiltered, maskFiltered, config);
	}
	else if (config.useAccurate == 1) {
		unsigned int dimPET[3];
		float voxelSize[3];

		image = readPrjFilePET(config.imageName, config.imageType, config.smoothingKernel, dimPET, voxelSize);
		mask = readVoiFilePET(config.imageName, config.voiName, image, config, dimPET, voxelSize);
	}
}


void calculateFeaturesForConfig(ImageType *imageFiltered, ImageType *maskNewSpacing, ConfigFile config) {
	//get size of shrinked image (in order to produce an image objectwith right size)
	const typename ImageType::RegionType regionFilter = imageFiltered->GetLargestPossibleRegion();
	const typename ImageType::SizeType imageSizeFilter = regionFilter.GetSize();
	vector<unsigned int> imageSize;
	imageSize.push_back(imageSizeFilter[0]);
	imageSize.push_back(imageSizeFilter[1]);
	imageSize.push_back(imageSizeFilter[2]);
	const typename ImageType::SpacingType& inputSpacing = imageFiltered->GetSpacing();
	vector<float> spacing = { float(inputSpacing[0]), float(inputSpacing[1]), float(inputSpacing[2]) };
	//now store the intensity values of thes selected region in an image object
	Image<float, 3> imageAttr(imageSize[0], imageSize[1], imageSize[2]);
	imageAttr.getImageAttributes(imageFiltered, maskNewSpacing, config);
	CalculateRelFeatures(imageAttr, imageFiltered, maskNewSpacing, config);
	if (config.useFixedBinWidth == 1 || config.useFixedNrBins == 1) {

		//we do not interpolate the image, so the image size is the same, however we still need to discretize it
		Image<float, 3> imageAttrDis(imageSize[0], imageSize[1], imageSize[2]);
		imageAttrDis.getImageAttributesDiscretized(imageFiltered, maskNewSpacing, config);
		calculateRelFeaturesDiscretized(imageAttrDis, spacing, config);
	}
	else {

		calculateRelFeaturesDiscretized(imageAttr, spacing, config);
	}
	std::cout << "The data is stored in the file " << config.outputFolder << std::endl;
}

void writeImageData2Log(ConfigFile config) {
	//write parameters to logFile
	string introduction = "The following software and images were used:";
	string date = "Software version " + version_date_nr;
	string imageNameForLog = "Image  " + config.imageName;
	string maskNameForLog = "Voi " + config.voiName;
	writeLogFile(config.outputFolder, introduction);
	writeLogFile(config.outputFolder, date);
	writeLogFile(config.outputFolder, imageNameForLog);
	writeLogFile(config.outputFolder, maskNameForLog);
	writeLogFile(config.outputFolder, string("#####################"));
	string configName = "The following config file was used:" + config.fileName;
	writeLogFile(config.outputFolder, configName);
	writeLogFile(config.outputFolder, string("#####################"));
	if (config.imageType == "PET") {

		string forLogFilePatInfo = "The following patient info file was used: " + config.patientInfoLocation;
		writeLogFile(config.outputFolder, forLogFilePatInfo);
		writeLogFile(config.outputFolder, string("Containing the following information:"));
		string logUseSUV = "UseSUV = " + to_string(config.useSUV);
		writeLogFile(config.outputFolder, logUseSUV);
		string logUseSUL = "UseSUL = " + to_string(config.useSUL);
		writeLogFile(config.outputFolder, logUseSUL);
		string logCorrection = "ScalingFactor = " + to_string(config.correctionParam);
		writeLogFile(config.outputFolder, logCorrection);
		string logPatientWeight = "PatientWeight = " + to_string(config.patientWeight);
		writeLogFile(config.outputFolder, logPatientWeight);
		string logPatientHeight = "PatientHeight = " + to_string(config.patientHeight);
		writeLogFile(config.outputFolder, logPatientHeight);
		string logGender = string("Gender = " + config.gender);
		writeLogFile(config.outputFolder, logGender);
		string logActivity = "ActivityMBq = " + to_string(config.initActivity);
		string unitsRescaling = "units_rescaling_factor = " + to_string(config.units_rescaling_factor);
		writeLogFile(config.outputFolder, unitsRescaling);
		writeLogFile(config.outputFolder, string("#####################"));
	}
}


