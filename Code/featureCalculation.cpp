//in this function the .voi file of the accurate project is read in
//First, the image and voxel size is read in from the .prj file
//then the .voi file is read in and transformed to an ITK image
ImageType::Pointer readVoiFilePET(string prjPath, string voiPath, ImageType *image, ConfigFile config, unsigned int(&dimPET)[3], float voxelSize[3]){
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
			maskImage->Update();
			voiFile.close();
			

			//get the bounding box of the image voi
			RegionType boundingBoxRegion = getBoundingBoxMask(maskImage);
			//filter mask and image 
			ImageType::Pointer maskFiltered = getImageMasked(maskImage, boundingBoxRegion);
			ImageType::Pointer imageFiltered = getImageMasked(image, boundingBoxRegion);
			const typename ImageType::SpacingType& inputSpacing = imageFiltered->GetSpacing();
			//std::cout << "spacingHERE" << inputSpacing[0] << " " << inputSpacing[1] << " " << inputSpacing[2] << std::endl;
			if (config.useDownSampling == 1 || config.useUpSampling == 1 && (inputSpacing[0] != inputSpacing[1] || inputSpacing[0] != inputSpacing[2] )) {
				const typename ImageType::RegionType& imageRegion = imageFiltered->GetLargestPossibleRegion();
				const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
				double outputSpacing[3];
				vector<int> newImageSize = getImageSizeInterpolated(imageFiltered, imageRegionSize, outputSpacing, config);
				itk::Size<3> outputSize = { { newImageSize[0], newImageSize[1], newImageSize[2] } };
				Image<float, 3> imageMask(newImageSize[0], newImageSize[1], newImageSize[2]);
				Image<float, 3> imageTest(newImageSize[0], newImageSize[1], newImageSize[2]);
				imageFiltered = imageTest.getResampledImage(imageFiltered, outputSpacing, outputSize, config.interpolationMethod);
				maskFiltered = imageTest.getResampledImage(maskFiltered, outputSpacing, outputSize, "Linear");
				//convert mask values to 1 (necessary after interpolation)
				//maskNew = maskValues2One(maskNew);
			}
			else if (config.useSampling2mm == 1 && (inputSpacing[0] != 2 && inputSpacing[1] != 2 && inputSpacing[2] != 2)) {
				const typename ImageType::RegionType& imageRegion = imageFiltered->GetLargestPossibleRegion();
				const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
				double outputSpacing[3];
				vector<int> newImageSize = getImageSizeInterpolated(imageFiltered, imageRegionSize, outputSpacing, config);
				itk::Size<3> outputSize = { { newImageSize[0], newImageSize[1], newImageSize[2] } };
				Image<float, 3> imageMask(newImageSize[0], newImageSize[1], newImageSize[2]);
				Image<float, 3> imageTest(newImageSize[0], newImageSize[1], newImageSize[2]);
				imageFiltered = imageTest.getResampledImage(imageFiltered, outputSpacing, outputSize, config.interpolationMethod);
				maskFiltered = imageTest.getResampledImage(maskFiltered, outputSpacing, outputSize, "Linear");
				//convert mask values to 1 (necessary after interpolation)
				//maskNew = maskValues2One(maskNew);
			}
			
			calculateFeaturesForConfig(imageFiltered, maskFiltered, config);
		}
	
	return maskImage;
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
		}
		if (config.smoothingKernel > 0) {
			image = smoothImage(image, config.smoothingKernel);
		}

		//resample image if required
		const typename ImageType::SpacingType& inputSpacing = image->GetSpacing();
		if ((config.useDownSampling != 0 || config.useUpSampling != 0) && (inputSpacing[0] != inputSpacing[1] || inputSpacing[0] != inputSpacing[2])) {
			const typename ImageType::RegionType& imageRegion = image->GetLargestPossibleRegion();
			const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
			double outputSpacing[3];
			vector<int> newImageSize = getImageSizeInterpolated(image, imageRegionSize, outputSpacing, config);
			itk::Size<3> outputSize = { { newImageSize[0], newImageSize[1], newImageSize[2] } };
			Image<float, 3> imageMask(newImageSize[0], newImageSize[1], newImageSize[2]);
			Image<float, 3> imageTest(newImageSize[0], newImageSize[1], newImageSize[2]);
			image = imageTest.getResampledImage(image, outputSpacing, outputSize, config.interpolationMethod);
			mask = imageTest.getResampledImage(mask, outputSpacing, outputSize, "Linear");
		}

		else if (config.useSampling2mm == 1 && (inputSpacing[0] != 2 && inputSpacing[1] != 2 && inputSpacing[2] != 2)) {
			const typename ImageType::RegionType& imageRegion = image->GetLargestPossibleRegion();			
			const typename ImageType::SizeType& imageRegionSize = imageRegion.GetSize();
			double outputSpacing[3];
			vector<int> newImageSize = getImageSizeInterpolated(image, imageRegionSize, outputSpacing, config);
			itk::Size<3> outputSize = { { newImageSize[0], newImageSize[1], newImageSize[2] } };
			Image<float, 3> imageMask(newImageSize[0], newImageSize[1], newImageSize[2]);
			Image<float, 3> imageTest(newImageSize[0], newImageSize[1], newImageSize[2]);
			image = imageTest.getResampledImage(image, outputSpacing, outputSize, config.interpolationMethod);
			mask = imageTest.getResampledImage(mask, outputSpacing, outputSize, "Linear");
		}
		//convert mask values to 1 (necessary after interpolation)
		mask = maskValues2One(mask);
		//get the region of the mask
		RegionType boundingBoxRegion = getBoundingBoxMask(mask);
		itk::Size<3> regionSize = boundingBoxRegion.GetSize();
		//shrink image and mask to the mask region
		ImageType::Pointer imageFiltered = getImageMasked(image, boundingBoxRegion);
		ImageType::Pointer maskFiltered = getImageMasked(mask, boundingBoxRegion);
		//change spacing of mask to spacing of image (in order to calculate the morphological features)
		ImageType::Pointer maskNewSpacing = getMaskNewSpacing(imageFiltered, maskFiltered);
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
	vector<double> spacing = { inputSpacing[0], inputSpacing[1], inputSpacing[2] };
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
		writeLogFile(config.outputFolder, string("#####################"));
	}
}


