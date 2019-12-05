template<typename T>
void fillConvolutionalVector(vector<T> &convolutionalVector, int is3D) {
	if (is3D == 1) {
		for (int depth = 0; depth < 3; depth++) {
			for (int row = 0; row < 3; row++) {
				for (int col = 0; col < 3; col++) {
					if (col*depth*row != 1) {
						convolutionalVector.push_back(1);
					}
					else {
						convolutionalVector.push_back(float(0.0));
					}
				}
			}

		}
	}
	else {
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				if (col*row != 1) {
					convolutionalVector.push_back(1);
				}
				else {
					convolutionalVector.push_back(float(0.0));
				}
			}
		}
	}
}


template<typename T>
ImageType::Pointer constructEmptyNewImage(Image<T, 3> imageAttr, boost::multi_array<T, 3> actualMatrix, int is3D) {
	//get the actual image and matrix
	ImageType::Pointer actImage = imageAttr.image;
	//get image size and spacing
	const typename ImageType::SpacingType& inputSpacing = actImage->GetSpacing();
	const typename ImageType::RegionType regionFilter = actImage->GetLargestPossibleRegion();
	const typename ImageType::SizeType imageSizeFilter = regionFilter.GetSize();
	unsigned int dimImage[] = { actualMatrix.shape()[0],actualMatrix.shape()[1],actualMatrix.shape()[2] };
	float voxelSizeImage[] = { inputSpacing[0],inputSpacing[1],inputSpacing[2] };

	//construct a new image
	ImportFilterType::SizeType size;
	//set image size
	size[0] = actualMatrix.shape()[0] + 2;
	size[1] = actualMatrix.shape()[1] + 2;
	size[2] = actualMatrix.shape()[2] + 2*is3D;
	
	int nrVoxelsPET = size[0] * size[1] * size[2];
	//start with importing the image
	ImageType::RegionType region;
	ImageType::IndexType start;
	start.Fill(0);
	region.SetIndex(start);
	region.SetSize(size);
	//set image regions
	ImageType::Pointer imageNew = ImageType::New();
	imageNew->SetRegions(region);
	imageNew->Allocate();
	ImageType::IndexType pixelIndex;
	for (int depth = 0; depth < actualMatrix.shape()[2] + 2*is3D; depth++) {
		for (int row = 0; row < actualMatrix.shape()[0] + 2; row++) {
			for (int col = 0; col < actualMatrix.shape()[1] + 2; col++) {
				pixelIndex[0] = row;
				pixelIndex[1] = col;
				pixelIndex[2] = depth;
				imageNew->SetPixel(pixelIndex, 0);
			}
		}
	}

	region.SetIndex(start);
	region.SetSize(size);
	return imageNew;
}

ImageType::Pointer convolutionImage(ImageType::Pointer image, ImageType::Pointer kernel) {
	using FilterType = itk::ConvolutionImageFilter<ImageType>;
	FilterType::Pointer convolutionFilter = FilterType::New();
	convolutionFilter->SetInput(image);
	convolutionFilter->SetKernelImage(kernel);
	convolutionFilter->Update();
	return convolutionFilter->GetOutput();

}

template<typename T>
ImageType::Pointer constructNewImage(ImageType::Pointer actImage, boost::multi_array<T, 3> actualMatrix) {

	//get image size
	const typename ImageType::RegionType regionFilter = actImage->GetLargestPossibleRegion();
	const typename ImageType::SizeType imageSizeFilter = regionFilter.GetSize();

	const typename ImageType::SpacingType& inputSpacing = actImage->GetSpacing();
	unsigned int dimImage[] = { actualMatrix.shape()[0],actualMatrix.shape()[1],actualMatrix.shape()[2] };
	float voxelSizeImage[] = { inputSpacing[0],inputSpacing[1],inputSpacing[2] };
	vector<float> actImageVector;

	ImportFilterType::SizeType size;
	//set image size
	size[0] = actualMatrix.shape()[0] + 2;
	size[1] = actualMatrix.shape()[1] + 2;
	size[2] = actualMatrix.shape()[2] + 2;


	int nrVoxelsPET = size[0] * size[1] * size[2];
	//start with importing the image
	ImageType::RegionType region;

	ImageType::IndexType start;
	start.Fill(0);

	region.SetIndex(start);



	region.SetSize(size);

	ImageType::Pointer imageNew = ImageType::New();
	imageNew->SetRegions(region);
	imageNew->Allocate();
	imageNew->FillBuffer(itk::NumericTraits<T>::Zero);

	ImageType::IndexType pixelIndex;
	region.SetIndex(start);
	region.SetSize(size);

	for (int depth = 1; depth < actualMatrix.shape()[2] + 1; depth++) {
		for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
			for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
				pixelIndex[0] = row;
				pixelIndex[1] = col;
				pixelIndex[2] = depth;
				if (row*col*depth == 0) {
					imageNew->SetPixel(pixelIndex, 0);
				}
				else {
					if (actualMatrix[row - 1][col - 1][depth - 1] == actGreyLevel) {

						imageNew->SetPixel(pixelIndex, 1);
					}
					else {
						imageNew->SetPixel(pixelIndex, 0);
					}
				}

			}
		}
	}
	return imageNew;

}

template<typename T>
void getNGLDMatrix3D_convolution(Image<T, 3> imageAttr, boost::multi_array<float, 2> &ngldm3D, vector<float> spacing, ConfigFile config) {
	int sizeGreyLevels = (imageAttr.diffGreyLevels).size();
	vector<float> convolutionalVector;
	fillConvolutionalVector(convolutionalVector, 1);
	float *kernelArray = &convolutionalVector[0];
	unsigned int dimKernel[] = { 3,3,3 };
	float voxelSize[] = { 1,1,1 };
	NGLDMFeatures3D<T, 3> ngldm;
	int indexOfElement[3] = { 0,0,0 };
	//get distance set for NGTDM case
	int dist = config.distNGLDM;
	T actualElement;
	float matrixValue;
	int actualGreyIndex;

	ImageType::Pointer kernel = converArray2Image(kernelArray, dimKernel, voxelSize);

	//the original (not discretize image!)
	ImageType::Pointer actImage = imageAttr.image;
	//assign the image matrix
	boost::multi_array<T, 3> actualMatrix = imageAttr.imageMatrix;
	for (int greyLevelNr = 0; greyLevelNr < sizeGreyLevels; greyLevelNr++) {
		float actGreyLevel = (imageAttr.diffGreyLevels)[greyLevelNr];
		//how many grey levels do I have

		

		ImageType::Pointer imageNew = constructEmptyNewImage(imageAttr, actualMatrix, 1);
		ImageType::IndexType pixelIndex;
		for (int depth = 1; depth < actualMatrix.shape()[2] + 1; depth++) {
			for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
				for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
					pixelIndex[0] = row;
					pixelIndex[1] = col;
					pixelIndex[2] = depth;
					if (row*col*depth == 0) {
						imageNew->SetPixel(pixelIndex, 0);
					}
					else {
						if (actualMatrix[row - 1][col - 1][depth - 1] == actGreyLevel) {

							imageNew->SetPixel(pixelIndex, 1);
						}
						else {
							imageNew->SetPixel(pixelIndex, 0);
						}
					}

				}
			}
		}


		ImageType::Pointer maskNew = constructEmptyNewImage(imageAttr, actualMatrix, 1);
		
		Image<T, 3> mask1(actualMatrix.shape()[0], actualMatrix.shape()[1], actualMatrix.shape()[2]);

		boost::multi_array<T, 3> maskTest = mask1.get3Dimage(imageAttr.mask, imageAttr.mask, config);
		int count = 0;
		for (int depth = 1; depth < actualMatrix.shape()[2] + 1; depth++) {
			for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
				for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
					pixelIndex[0] = row;
					pixelIndex[1] = col;
					pixelIndex[2] = depth;
					if (row*col*depth == 0) {
						maskNew->SetPixel(pixelIndex, 0);
					}
					else {
						if (!std::isnan(actualMatrix[row - 1][col - 1][depth - 1])) {
							maskNew->SetPixel(pixelIndex, 1);
							count++;
						}

					}

				}
			}
		}
		const typename ImageType::RegionType& imageRegion = imageNew->GetLargestPossibleRegion();
		const typename ImageType::SizeType& imageSize = imageRegion.GetSize();


		ImageType::Pointer sumImage = convolutionImage(imageNew, kernel);
		ImageType::Pointer sumMask = convolutionImage(maskNew, kernel);


		//get values of image and mask
		Image<T, 3> imageElement(imageSize[0], imageSize[1], imageSize[2]);
		Image<T, 3> maskSum(imageSize[0], imageSize[1], imageSize[2]);
		boost::multi_array<T, 3> sumMatrix = imageElement.get3Dimage(sumImage, maskNew, config);
		boost::multi_array<T, 3> maskSummation = maskSum.get3Dimage(sumMask, maskNew, config);
		imageNew = nullptr;
		sumImage = nullptr;
		sumMask = nullptr;
		//for every matrix element we access the neighborhood
		for (int depth = 0; depth < actualMatrix.shape()[2]; depth++) {
			for (int row = 0; row < actualMatrix.shape()[0]; row++) {

				for (int col = 0; col < actualMatrix.shape()[1]; col++) {
					//store actual index in a vector
					indexOfElement[0] = row;
					indexOfElement[1] = col;
					indexOfElement[2] = depth;
					//get actual Element if it is the centre of a whole neighborhood
					actualElement = actualMatrix[row][col][depth];
					if (actualElement == actGreyLevel) {
						//get the sum of the actual neighborhood
						actualGreyIndex = ngldm.findIndex(imageAttr.diffGreyLevels, sizeGreyLevels, actualElement);
						//only if distance of NGLDM and NGTDM matrices are the same, we can combine the calculations
						//if not we have to calculate the NGLDM matrices separately what means more computational time
						ngldm3D[actualGreyIndex][sumMatrix[row + 1][col + 1][depth + 1]] = ngldm3D[actualGreyIndex][sumMatrix[row + 1][col + 1][depth + 1]] + 1;
					}
					
				}
			}
		}


	}
}



//generate the image that is later used for the convolution in order to get the NGTDM matrix
template<typename T>
ImageType::Pointer generateImageForConvolutionNGTDM(Image<T,3> imageAttr) {
	//assign the image matrix
	boost::multi_array<T, 3> actualMatrix = imageAttr.imageMatrix;
	//how many grey levels do I have
	int sizeGreyLevels = (imageAttr.diffGreyLevels).size();
	//the original (not discretize image!)
	ImageType::Pointer actImage = imageAttr.image;
	//generate new image with one padding layer
	
	ImageType::Pointer imageNew = constructEmptyNewImage(imageAttr, actualMatrix, 1);
	ImageType::IndexType pixelIndex;
	for (int depth = 1; depth < actualMatrix.shape()[2] + 1; depth++) {
		for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
			for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
				pixelIndex[0] = row;
				pixelIndex[1] = col;
				pixelIndex[2] = depth;
				if (row*col*depth == 0) {
					imageNew->SetPixel(pixelIndex, 0);
				}
				else {
					if (!std::isnan(actualMatrix[row - 1][col - 1][depth - 1])) {

						imageNew->SetPixel(pixelIndex, actualMatrix[row - 1][col - 1][depth - 1]);
					}
					else {
						imageNew->SetPixel(pixelIndex, 0);
					}
				}
			}
		}
	}
	return imageNew;
}


//generate the mask that is later used for the convolution in order to get the NGTDM matrix
template<typename T>
ImageType::Pointer generateMaskForConvolutionNGTDM(Image<T, 3> imageAttr, ConfigFile config) {
	//assign the image matrix
	boost::multi_array<T, 3> actualMatrix = imageAttr.imageMatrix;
	//the original (not discretize image!)
	ImageType::Pointer actImage = imageAttr.image;
		//generate mask image
	ImageType::Pointer maskNew = constructEmptyNewImage(imageAttr, actualMatrix, 1);
	
	Image<T, 3> mask1(actualMatrix.shape()[0], actualMatrix.shape()[1], actualMatrix.shape()[2]);

	boost::multi_array<T, 3> maskTest = mask1.get3Dimage(imageAttr.mask, imageAttr.mask, config);
	int count = 0;
	ImageType::IndexType pixelIndex;
	for (int depth = 1; depth < actualMatrix.shape()[2] + 1; depth++) {
		for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
			for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
				pixelIndex[0] = row;
				pixelIndex[1] = col;
				pixelIndex[2] = depth;
				if (row*col*depth == 0) {
					maskNew->SetPixel(pixelIndex, 0);
				}
				else {
					if (!std::isnan(actualMatrix[row - 1][col - 1][depth - 1])) {
						maskNew->SetPixel(pixelIndex, 1);
						count++;
					}

				}

			}
		}
	}
	return maskNew;
}

template<typename T>
void getNeighborhoodMatrix3D_convolution(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm3D, vector<float> spacing, ConfigFile config) {
	//construct convolutional kernel
	vector<float> convolutionalVector;
	fillConvolutionalVector(convolutionalVector, 1);
	
	float *kernelArray = &convolutionalVector[0];
	unsigned int dimKernel[] = { 3,3,3 };
	float voxelSize[] = { 1,1,1 };
	ImageType::Pointer kernel = converArray2Image(kernelArray, dimKernel, voxelSize);
	
	
	NGLDMFeatures3D<T, 3> ngldm;
	int indexOfElement[3] = { 0,0,0 };
	//get distance set for NGTDM case
	int dist = config.dist;
	T actualElement;
	float matrixValue;
	int actualGreyIndex;

	//assign the image matrix
	boost::multi_array<T, 3> actualMatrix = imageAttr.imageMatrix;
	//how many grey levels do I have
	int sizeGreyLevels = (imageAttr.diffGreyLevels).size();
	//the original (not discretize image!)
	ImageType::Pointer actImage = imageAttr.image;
	//generate new image with one padding layer
	

	ImageType::Pointer imageNew = generateImageForConvolutionNGTDM(imageAttr);
	ImageType::Pointer maskNew = generateMaskForConvolutionNGTDM(imageAttr, config);

	const typename ImageType::RegionType& imageRegion = imageNew->GetLargestPossibleRegion();
	const typename ImageType::SizeType& imageSize = imageRegion.GetSize();
	

	ImageType::Pointer sumImage = convolutionImage(imageNew, kernel);
	ImageType::Pointer sumMask = convolutionImage(maskNew, kernel);
	

	//get values of image and mask
	Image<T, 3> imageElement(imageSize[0], imageSize[1], imageSize[2]);
	Image<T, 3> maskSum(imageSize[0], imageSize[1], imageSize[2]);
	boost::multi_array<T, 3> sumMatrix = imageElement.get3Dimage(sumImage, maskNew, config);
	boost::multi_array<T, 3> maskSummation = maskSum.get3Dimage(sumMask, maskNew, config);

	imageNew = nullptr;
	sumImage = nullptr;
	sumMask = nullptr;
	//for every matrix element we access the neighborhood
	for (int depth = 0; depth < actualMatrix.shape()[2]; depth++) {
		for (int row = 0; row < actualMatrix.shape()[0]; row++) {

			for (int col = 0; col < actualMatrix.shape()[1]; col++) {
				//store actual index in a vector
				indexOfElement[0] = row;
				indexOfElement[1] = col;
				indexOfElement[2] = depth;
				//get actual Element if it is the centre of a whole neighborhood
				actualElement = actualMatrix[row][col][depth];
				if (!std::isnan(actualElement) ) {
					//get the sum of the actual neighborhood
					matrixValue = sumMatrix[row+1][col+1][depth+1]/ maskSummation[row+1][col+1][depth+1];// / maskSummation[row][col][depth];//getNeighborhood3D_convolution(actualMatrix, indexOfElement, spacing, config);
					//std::cout << sumMatrix[row + 1][col + 1][depth + 1] << " 	";														 //std::cout << sumMatrix[row][col][depth + 2] << std::endl;
					ngtdm3D[row][col][depth] = matrixValue;
					actualGreyIndex = ngldm.findIndex(imageAttr.diffGreyLevels, sizeGreyLevels, actualElement);
				}
				else {
					ngtdm3D[row][col][depth] = 0;
				}
			}
		}
	}
}



/*!
\brief getNeighborhoodMatrix2D
In this function matrices necessary to fill the NGTDM and NGLDM matrices are generated for the 2D case. \n
I do this in one function to save computational time.  \n

The matrix for the NGTDM matrix contains the sum of all neighborhoods, while the matrix for the NGLDM matrix
contains the number of elements smaller than the coarseness factor. \n
NGTDM and NGLDM matrix are given as reference to this function. \n

@param[in] Image imageAttr: image attribute of actual image
@param[in] boost multi array ngtdm2D: matrix to store results for ngtdm calculations, given as reference
@param[in] boost multi array ngldm2D: matrix to store results for ngldm calculations, given as reference
@param[in] vector spacing: spacing of actual image, necessary to calculate weights
@param[in] ConfigFile config: ConfigFile to read the distances the user set for the NGTDM and NGLDM case
*/
template<typename T>
void getNeighborhoodMatrix2D(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngtdm2D, vector<float> spacing, ConfigFile config) {
	//construct convolutional kernel
	vector<float> convolutionalVector;
	fillConvolutionalVector(convolutionalVector, 0);


	
	float *kernelArray = &convolutionalVector[0];
	unsigned int dimKernel[] = { 3,3,1 };
	float voxelSize[] = { 1,1,1 };
	ImageType::Pointer kernel = converArray2Image(kernelArray, dimKernel, voxelSize);
	int indexOfElement[3] = { 0,0,0 };
	NGLDMFeatures2DMRG<T, 3> ngldm;
	int dist = config.dist;
	T actualElement;

	//get distance set for NGTDM case
	float matrixValue;
	int actualGreyIndex;

	//assign the image matrix
	boost::multi_array<T, 3> actualMatrix = imageAttr.imageMatrix;
	//how many grey levels do I have
	int sizeGreyLevels = (imageAttr.diffGreyLevels).size();
	//the original (not discretize image!)
	ImageType::Pointer actImage = imageAttr.image;

	ImageType::Pointer imageNew = constructEmptyNewImage(imageAttr, actualMatrix, 0);
	
	ImageType::IndexType pixelIndex;

	for (int depth =0; depth < actualMatrix.shape()[2] ; depth++) {
		for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
			for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
				pixelIndex[0] = row;
				pixelIndex[1] = col;
				pixelIndex[2] = depth;
				if (row*col == 0) {
					imageNew->SetPixel(pixelIndex, 0);
				}
				else {
					if (!std::isnan(actualMatrix[row - 1][col - 1][depth])) {

						imageNew->SetPixel(pixelIndex, actualMatrix[row - 1][col - 1][depth]);
					}
					else {
						imageNew->SetPixel(pixelIndex, 0);
					}
				}
			}
		}
	}

	//generate mask image
	ImageType::Pointer maskNew = constructEmptyNewImage(imageAttr, actualMatrix, 0);
	

	Image<T, 3> mask1(actualMatrix.shape()[0], actualMatrix.shape()[1], actualMatrix.shape()[2]);

	boost::multi_array<T, 3> maskTest = mask1.get3Dimage(imageAttr.mask, imageAttr.mask, config);
	int count = 0;
	for (int depth = 0; depth < actualMatrix.shape()[2] ; depth++) {
		for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
			for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
				pixelIndex[0] = row;
				pixelIndex[1] = col;
				pixelIndex[2] = depth;
				if (row*col == 0) {
					maskNew->SetPixel(pixelIndex, 0);
				}
				else {
					if (!std::isnan(actualMatrix[row - 1][col - 1][depth])) {
						maskNew->SetPixel(pixelIndex, 1);
						count++;
					}

				}

			}
		}
	}
	const typename ImageType::RegionType& imageRegion = imageNew->GetLargestPossibleRegion();
	const typename ImageType::SizeType& imageSize = imageRegion.GetSize();


	ImageType::Pointer sumImage = convolutionImage(imageNew, kernel);
	ImageType::Pointer sumMask = convolutionImage(maskNew, kernel);

	imageNew = nullptr;

	//get values of image and mask
	Image<T, 3> imageElement(imageSize[0], imageSize[1], imageSize[2]);
	Image<T, 3> maskSum(imageSize[0], imageSize[1], imageSize[2]);
	boost::multi_array<T, 3> sumMatrix = imageElement.get3Dimage(sumImage, maskNew, config);
	boost::multi_array<T, 3> maskSummation = maskSum.get3Dimage(sumMask, maskNew, config);
	sumImage = nullptr;
	sumMask = nullptr;
	//for every matrix element we access the neighborhood
	for (int depth = 0; depth < actualMatrix.shape()[2]; depth++) {
		for (int row = 0; row < actualMatrix.shape()[0]; row++) {

			for (int col = 0; col < actualMatrix.shape()[1]; col++) {
				//store actual index in a vector
				indexOfElement[0] = row;
				indexOfElement[1] = col;
				indexOfElement[2] = depth;
				//get actual Element if it is the centre of a whole neighborhood
				actualElement = actualMatrix[row][col][depth];
				if (!std::isnan(actualElement)) {
					//get the sum of the actual neighborhood
					matrixValue = sumMatrix[row + 1][col + 1][depth ] / maskSummation[row + 1][col + 1][depth ];// / maskSummation[row][col][depth];//getNeighborhood3D_convolution(actualMatrix, indexOfElement, spacing, config);
					ngtdm2D[row][col][depth] = matrixValue;
					actualGreyIndex = ngldm.findIndex(imageAttr.diffGreyLevels, sizeGreyLevels, actualElement);
					//only if distance of NGLDM and NGTDM matrices are the same, we can combine the calculations
					//if not we have to calculate the NGLDM matrices separately what means more computational time

				}
				else {
					ngtdm2D[row][col][depth] = 0;

				}
			}

		}
	}
	maskNew = nullptr;
	
}





template<typename T>
void getNeighborhoodMatrix2DNGLDM(Image<T, 3> imageAttr, boost::multi_array<T, 3> &ngldm2D, vector<float> spacing, ConfigFile config) {
	int sizeGreyLevels = (imageAttr.diffGreyLevels).size();
	vector<float> convolutionalVector;
	fillConvolutionalVector(convolutionalVector, 0);
	float *kernelArray = &convolutionalVector[0];
	unsigned int dimKernel[] = { 3,3,1 };
	float voxelSize[] = { 1,1,1 };

	NGLDMFeatures3D<T, 3> ngldm;
	int indexOfElement[3] = { 0,0,0 };
	//get distance set for NGTDM case
	int dist = config.distNGLDM;
	T actualElement;
	float matrixValue;
	int actualGreyIndex;

	ImageType::Pointer actImage = imageAttr.image;
	boost::multi_array<T, 3> actualMatrix = imageAttr.imageMatrix;

	ImageType::Pointer kernel = converArray2Image(kernelArray, dimKernel, voxelSize);
	for (int greyLevelNr = 0; greyLevelNr < sizeGreyLevels; greyLevelNr++) {
		float actGreyLevel = (imageAttr.diffGreyLevels)[greyLevelNr];
		
		ImageType::Pointer imageNew = constructEmptyNewImage(imageAttr, actualMatrix, 0);
	
		ImageType::IndexType pixelIndex;
	
		for (int depth = 0; depth < actualMatrix.shape()[2]; depth++) {
			for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
				for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
					pixelIndex[0] = row;
					pixelIndex[1] = col;
					pixelIndex[2] = depth;
					if (row*col == 0) {
						imageNew->SetPixel(pixelIndex, 0);
					}
					else {
						if (actualMatrix[row - 1][col - 1][depth] == actGreyLevel) {

							imageNew->SetPixel(pixelIndex, 1);
						}
						else {
							imageNew->SetPixel(pixelIndex, 0);
						}
					}

				}
			}
		}


		ImageType::Pointer maskNew = constructEmptyNewImage(imageAttr, actualMatrix, 0);
		
		ImageType::IndexType pixelIndexMask;

		Image<T, 3> mask1(actualMatrix.shape()[0], actualMatrix.shape()[1], actualMatrix.shape()[2]);

		boost::multi_array<T, 3> maskTest = mask1.get3Dimage(imageAttr.mask, imageAttr.mask, config);
		int count = 0;
		for (int depth = 0; depth < actualMatrix.shape()[2] ; depth++) {
			for (int row = 1; row < actualMatrix.shape()[0] + 1; row++) {
				for (int col = 1; col < actualMatrix.shape()[1] + 1; col++) {
					pixelIndex[0] = row;
					pixelIndex[1] = col;
					pixelIndex[2] = depth;
					if (row*col == 0) {
						maskNew->SetPixel(pixelIndex, 0);
					}
					else {
						if (!std::isnan(actualMatrix[row - 1][col - 1][depth])) {
							maskNew->SetPixel(pixelIndex, 1);
							count++;
						}

					}

				}
			}
		}
		
		ImageType::Pointer sumImage = convolutionImage(imageNew, kernel);
		ImageType::Pointer sumMask = convolutionImage(maskNew, kernel);
		const typename ImageType::RegionType& imageRegion = imageNew->GetLargestPossibleRegion();
		const typename ImageType::SizeType& imageSize = imageRegion.GetSize();

		//get values of image and mask
		Image<T, 3> imageElement(imageSize[0], imageSize[1], imageSize[2]);
		Image<T, 3> maskSum(imageSize[0], imageSize[1], imageSize[2]);
		boost::multi_array<T, 3> sumMatrix = imageElement.get3Dimage(sumImage, maskNew, config);
		boost::multi_array<T, 3> maskSummation = maskSum.get3Dimage(sumMask, maskNew, config);
		imageNew = nullptr;
		sumImage = nullptr;
		sumMask = nullptr;
		//for every matrix element we access the neighborhood
		for (int depth = 0; depth < actualMatrix.shape()[2]; depth++) {
			for (int row = 0; row < actualMatrix.shape()[0]; row++) {

				for (int col = 0; col < actualMatrix.shape()[1]; col++) {
					//store actual index in a vector
					indexOfElement[0] = row;
					indexOfElement[1] = col;
					indexOfElement[2] = depth;
					//get actual Element if it is the centre of a whole neighborhood
					actualElement = actualMatrix[row][col][depth];
					if (actualElement == actGreyLevel) {
						//get the sum of the actual neighborhood
						actualGreyIndex = ngldm.findIndex(imageAttr.diffGreyLevels, sizeGreyLevels, actualElement);
						//only if distance of NGLDM and NGTDM matrices are the same, we can combine the calculations
						//if not we have to calculate the NGLDM matrices separately what means more computational time
						ngldm2D[actualGreyIndex][sumMatrix[row + 1][col + 1][depth]][depth] = ngldm2D[actualGreyIndex][sumMatrix[row + 1][col + 1][depth]][depth] + 1;
					}

				}
			}
		}
	}

}
