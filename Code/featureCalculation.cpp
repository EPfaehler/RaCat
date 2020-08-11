int getNrVoxels(ImageType *mask, float threshold) {
	int nrInMask = 0;
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
	countVoxels.GoToBegin();

	while (!countVoxels.IsAtEnd())
	{
		if (countVoxels.Get() > threshold*maxValue)
		{
			nrInMask += 1;
			//std::cout << maxValue << std::endl;
		}
		++countVoxels;
	}
	return nrInMask;
}



ImageType::Pointer flipNII(ImageType::Pointer mask) {
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
	return mask;
}
//in this function the .voi file of the accurate project is read in
//First, the image and voxel size is read in from the .prj file
//then the .voi file is read in and transformed to an ITK image


void readImageAndMask(ConfigFile config) {
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
	if (config.useAccurate == 0) {
		image = readImage(config.imageName);
		std::cout << "The input image is a nifti image" << std::endl;
		//to get same values for .voi/.prj and .nii, we need to flip the images
		image = flipNII(image);
		mask = readImage(config.voiName);
		mask = flipNII(mask);
	}
	if (config.useAccurate == 1) {
		unsigned int dimPET[3];
		float voxelSize[3];
		image = readPrjFilePET(config.imageName, config.imageType, config.smoothingKernel, dimPET, voxelSize);
		mask = readVoiFilePET(config.imageName, config.voiName, image, config, dimPET, voxelSize);
	}
	if (config.useAccurate == 2) {
		image = readDicom(config.imageName);

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
	}
	
		
	if (config.smoothingKernel > 0) {
		image = smoothImage(image, config.smoothingKernel);
	}
	
	//resample image if required
	//get the bounding box of the image voi
	//RegionType boundingBoxRegion = getBoundingBoxMask(mask);
	ImageType::Pointer maskFiltered = mask;
	ImageType::Pointer imageFiltered = image;
	int nrVoxelsInMask = getNrVoxels(maskFiltered, config.threshold);
	
	
	const typename ImageType::SpacingType& inputSpacing = image->GetSpacing();
	storePreInterpolationFeatures(image, mask, config);
	if (nrVoxelsInMask < 5) {
		fillCSVwithNANs(config);
	}
	else{
		//now down or upsample the image
		if (config.useSamplingCubic == 1 || config.useDownSampling != 0 || config.useUpSampling != 0) {
			const typename ImageType::RegionType& imageRegion = image->GetLargestPossibleRegion();
			const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
			double outputSpacing[3];
			vector<int> newImageSize = getImageSizeInterpolated(image, imageRegionSize, outputSpacing, config);
			itk::Size<3> outputSize = { { double(newImageSize[0]), double(newImageSize[1]), double(newImageSize[2]) } };
			Image<float, 3> imageMask(1, 1, 1);
			imageFiltered = imageMask.getResampledImage(imageFiltered, outputSpacing, outputSize, config.interpolationMethod, config.rebinning_centering);
			maskFiltered = imageMask.getResampledImage(maskFiltered, outputSpacing, outputSize, "Linear", config.rebinning_centering);
		}
		//convert mask values to 1 (necessary after interpolation)
		maskFiltered = thresholdMask(maskFiltered, config.threshold);
		//get the region of the mask
		RegionType boundingBoxRegion = getBoundingBoxMask(maskFiltered);
		itk::Size<3> regionSize = boundingBoxRegion.GetSize();
		//shrink image and mask to the mask region
		imageFiltered = getImageMasked(imageFiltered, boundingBoxRegion);
		maskFiltered = getImageMasked(maskFiltered, boundingBoxRegion);

		//for ontology table
		const typename ImageType::SpacingType& spacingVoxelDim = imageFiltered->GetSpacing();
		float voxelSize[3];
		for (int i = 0; i < 3; i++) {
			voxelSize[i] = (float)spacingVoxelDim[i];
		}
		Image<float, 3> imageVoxelDim(10, 10, 10);
		imageVoxelDim.createOntologyVoxelDimensionTable(config, voxelSize);
		mask = nullptr;
		image = nullptr;
		calculateFeatures(imageFiltered, maskFiltered, config);
		
	}	
}


void calculateFeatures(ImageType *imageFiltered, ImageType *maskNewSpacing, ConfigFile config) {
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
	CalculateRelFeatures(imageAttr, config);
	Image<float, 3> imageAttr2(0, 0, 0);
	imageAttr = imageAttr2;
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


