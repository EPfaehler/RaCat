/*!
The method readImage reads the itk-image
@param[in]:string imageName: the name of the image that should be read
@param[out]: ITK-image
*/
ImageType::Pointer readImage(string imageName) {
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(imageName);
	try {
		reader->Update();
	}
	catch (itk::ExceptionObject &excp) {
		std::cerr << excp << std::endl;
	}
	ImageType::Pointer image = reader->GetOutput();

	return  reader->GetOutput();
}

/*!
The method getBoundingBoxMask generates the bounding box of the actual mask
@param[in]: string maskName: the name of the mask
@param[out]: bounding box region
*/
RegionType getBoundingBoxMask(ImageType *mask) {

	ImageMaskSpatialObject::Pointer maskSO = ImageMaskSpatialObject::New();
	CastFilterTypeChar::Pointer castFilterChar = CastFilterTypeChar::New();
	castFilterChar->SetInput(mask);
	castFilterChar->Update();
	charImage *castImage = castFilterChar->GetOutput();
	castImage->Update();
	maskSO->SetImage(castImage);
	maskSO->Update();
	//get the bounding box of the mask
	RegionType boundingBoxRegion = maskSO->GetAxisAlignedBoundingBoxRegion();
	return boundingBoxRegion;
}

/*!
The method getImageMasked cuts the image, so that only the part of the bounding box remains
@param[in]: ImageType image: ITK image
@param[out]: RegionType boundingBoxRegion: region of the mask bounding box
*/
ImageType::Pointer getImageMasked(ImageType *image, RegionType boundingBoxRegion) {
	//set the bounding box filter
	FilterType::Pointer filter = FilterType::New();
	filter->SetRegionOfInterest(boundingBoxRegion);
	filter->SetInput(image);
	filter->Update();
	ImageType::Pointer imageMask = filter->GetOutput();
	const typename ImageType::RegionType regionFilter = imageMask->GetLargestPossibleRegion();
	const typename ImageType::SizeType imageSizeFilter = regionFilter.GetSize();
	return imageMask;
}

/*!
The method getMaskNewSpacing changes the spacing of the mask, so that it fits with the spacing of the image
@param[in]: ImageType imageFiltered: ITK image
@param[in]: ImageType maskFiltered: ITK image
@param[out]: ImageType maskNewSpacing: ITK mask with new spacing
*/
ImageType::Pointer getMaskNewSpacing(ImageType *imageFiltered, ImageType *maskFiltered) {
	ChangeInfoFilter::Pointer changeInfo = ChangeInfoFilter::New();
	const typename ImageType::RegionType regionFilter = imageFiltered->GetLargestPossibleRegion();
	const typename ImageType::SizeType imageSizeFilter = regionFilter.GetSize();
	unsigned int oldWidth = imageSizeFilter[0];
	unsigned int oldHeight = imageSizeFilter[1];
	unsigned int oldDepth = imageSizeFilter[2];
	//change spacing of the mask
	const typename ImageType::SpacingType& inputSpacing = imageFiltered->GetSpacing();
	changeInfo->SetOutputSpacing(inputSpacing);
	changeInfo->ChangeSpacingOn();
	changeInfo->SetInput(maskFiltered);
	changeInfo->Update();
	ImageType::Pointer maskNewSpacing = changeInfo->GetOutput();
	return maskNewSpacing;
}

/*!
The method smoothImage smoothes the input image with the defined kernel using Gaussian smoothing. The kernel can be set in the .ini file.
If the kernel is not equal to 0, smoothing is applied.
@param[in]: ImageType image: ITK image
@param[out]:  ImageType image: ITK image: filtered ITK image
*/
ImageType::Pointer smoothImage(ImageType *image, float kernel) {
	using FilterType = itk::RecursiveGaussianImageFilter<ImageType, ImageType >;
	//set filter for each direction
	FilterType::Pointer filterX = FilterType::New();
	FilterType::Pointer filterY = FilterType::New();
	FilterType::Pointer filterZ = FilterType::New();
	filterX->SetDirection(0);   // 0 --> X direction
	filterY->SetDirection(1);   // 1 --> Y direction
	filterZ->SetDirection(2);   // 2 --> Z direction
	//Gaussian filter can smooth the image with using the Gaussian function or one of its derivatives
	//with setting order to zero order we use the Gaussian function (first order would be first derivative etc)
	filterX->SetOrder(FilterType::ZeroOrder);
	filterY->SetOrder(FilterType::ZeroOrder);
	filterZ->SetOrder(FilterType::ZeroOrder);
	//disable normalization flag
	filterX->SetNormalizeAcrossScale(false);
	filterY->SetNormalizeAcrossScale(false);
	filterZ->SetNormalizeAcrossScale(false);
	//set input
	filterX->SetInput(image);
	filterY->SetInput(filterX->GetOutput());
	filterZ->SetInput(filterY->GetOutput());
	//set smoothing kernel
	filterX->SetSigma(kernel);
	filterY->SetSigma(kernel);
	filterZ->SetSigma(kernel);
	filterZ->Update();
	return filterZ->GetOutput();
}

/*!
The method getImageSizeInterpolated determines the new image size after interpolation to cubic voxels
That is necessary, as I need the image size to determine the size of the boost::multi_array
@param[in]: ImageType imageFiltered: ITK image
@param[in]: vector<unsigned int> imageSize: original image size
@param[in]: ConfigFile config: config file with all information of the config.ini file
@param[out]: vector<unsigned int> newImageSize: image size after interpolation
*/
vector<int> getImageSizeInterpolated(ImageType *imageFiltered, ImageType::SizeType imageSize, double(&outputSpacing)[3], ConfigFile config) {

	const typename ImageType::SpacingType& inputSpacing = imageFiltered->GetSpacing();
#ifdef _WIN32
	vector<int> newImageSize = { 0, 0, 0 };
#else
	vector<int> newImageSize(3,0);
#endif // _WIN32
	//if the image should be downsampled, we use the minimum spacing value as new spacing
	//for this, we have first to get the minimum spacing value
	if (config.useDownSampling == 1) {
		int minimum = inputSpacing[0];
		if (inputSpacing[1]<minimum) {
			minimum = inputSpacing[1];
		}
		if (inputSpacing[2] < minimum) {
			minimum = inputSpacing[2];
		}
		outputSpacing[0] = minimum;
		outputSpacing[1] = minimum;
		if (config.interpolation2D == 0) {
			outputSpacing[2] = minimum;
			newImageSize[2] = (float)imageSize[2] * inputSpacing[2] / minimum;
		}
		else {
			outputSpacing[2] = inputSpacing[2];
			newImageSize[2] = (float)imageSize[2];
		}
		//calculate now the new image size
		newImageSize[0] = (float)imageSize[0] * inputSpacing[0] / minimum;
		newImageSize[1] = (float)imageSize[1] * inputSpacing[1] / minimum;
		
	}
	//if the image should  be up sampled, we use the maximum spacing value
	//for this, the maximum value has to be determined
	else if(config.useUpSampling == 1) {
		int maximum = inputSpacing[0];
		if (inputSpacing[1]>maximum) {
			maximum = inputSpacing[1];
		}
		if (inputSpacing[2] > maximum) {
			maximum = inputSpacing[2];
		}
		outputSpacing[0] = maximum;
		outputSpacing[1] = maximum;
		//get the new image size
		newImageSize[0] = (float)imageSize[0] * inputSpacing[0] / maximum;
		newImageSize[1] = (float)imageSize[1] * inputSpacing[1] / maximum;
		if (config.interpolation2D == 0) {
			outputSpacing[2] = maximum;
			newImageSize[2] = (float)imageSize[2] * inputSpacing[2] / maximum;
		}
		else {
			outputSpacing[2] = inputSpacing[2];
			newImageSize[2] = (float)imageSize[2] ;
		}
		
		
	}
	else if (config.useSamplingCubic == 1) {
		outputSpacing[0] = float(config.cubicVoxelSize);
		outputSpacing[1] = float(config.cubicVoxelSize);
		newImageSize[0] = ceil((float)imageSize[0] * inputSpacing[0] / float(config.cubicVoxelSize));
		newImageSize[1] = ceil((float)imageSize[1] * inputSpacing[1] / float(config.cubicVoxelSize));
		if (config.interpolation2D == 0) {
			outputSpacing[2] = float(config.cubicVoxelSize);
			newImageSize[2] = ceil((float)imageSize[2] * inputSpacing[2] / float(config.cubicVoxelSize));
		}
		else {
			outputSpacing[2] = inputSpacing[2];
			newImageSize[2] = (float)imageSize[2];
		}
		
		
	}
	return newImageSize;
}

/*!
\brief maskValues2One
After resampling, the mask can contain intensity values from 0-1
All values bigger (or equal) to 0.5 should be included in the mask
As the function getBoundingBox region only regards 1 as inside the mask, we have to transform the mask
@param[in] originalMask

*/
ImageType::Pointer maskValues2One(ImageType *originalMask) {
	Image<float, 3> tmpImage(1, 1, 1);
	float valueInMask = tmpImage.getValueInMask(originalMask);
	float thresholdValue = 0.5*valueInMask;
	std::cout <<"!!!"<< valueInMask<<0.5*valueInMask << std::endl;
	const typename ImageType::SpacingType& inputSpacing = originalMask->GetSpacing();
	const typename ImageType::RegionType& inputRegion = originalMask->GetLargestPossibleRegion();
	const typename ImageType::SizeType& inputSize = inputRegion.GetSize();
	int actPosition = 0;
	int row = 0;
	int col = 0;
	int depth = 0;
	ImageType::IndexType pixelIndex;

	while (actPosition < inputSize[0] * inputSize[1] * inputSize[2]) {
		if (row < inputSize[0] && depth < inputSize[2]) {
			pixelIndex[0] = row;
			pixelIndex[1] = col;
			pixelIndex[2] = depth;
			if (originalMask->GetBufferPointer()[actPosition] >= 0.5*valueInMask) {
				originalMask->SetPixel(pixelIndex, 100);
			}
			
			row = row + 1;
			if (row == inputSize[0] && col < inputSize[1]) {
				col = col + 1;
				row = 0;
			}
			if (col == inputSize[1] && depth < inputSize[2]) {
				col = 0;
				depth = depth + 1;
			}
			actPosition += 1;
		}
	}
	return originalMask;
}

