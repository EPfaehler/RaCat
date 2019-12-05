ImageType::Pointer reorderFile(ImageType::Pointer imageVoiFile) {
	const typename ImageType::RegionType regionFilter = imageVoiFile->GetLargestPossibleRegion();
	const typename ImageType::SizeType imageSizeFilter = regionFilter.GetSize();
	unsigned int imageSize[3] = { imageSizeFilter[0], imageSizeFilter[1],imageSizeFilter[2] };
	const typename ImageType::SpacingType& inputSpacing = imageVoiFile->GetSpacing();
	float voxelSize[3] = { inputSpacing [0], inputSpacing [1], inputSpacing[2]};
	boost::multi_array<float, 3> A(boost::extents[imageSize[0]][imageSize[1]][imageSize[2]]);
	int row = 0;
	int depth = 0;
	int col = 0;
	int actPosition = 0;
	//std::cout << "max mask" << maxValueInMask << std::endl;
	while (actPosition < imageSize[0] * imageSize[1] * imageSize[2]) {
		if (row < imageSize[0] && depth < imageSize[2]) {

			A[row][col][depth] = imageVoiFile->GetBufferPointer()[actPosition];
				
			row = row + 1;
		}

		if (row == imageSize[0] && col < imageSize[1]) {
			col = col + 1;
			row = 0;
		}
		if (col == imageSize[1] && depth < imageSize[2]) {
			col = 0;
			depth = depth + 1;
		}
		actPosition += 1;
	}
	vector<float> reorder;
	for (int row = 0; row < imageSize[0]; row++) {
		for (int col = 0; col < imageSize[1]; col++) {
			for (int depth = 0; depth < imageSize[2]; depth++) {
				reorder.push_back(A[imageSize[0] - row-1][col][imageSize[2] - depth-1]);
			}
		}
	}
	std::cout << "TESTEST" << boost::size(reorder) << std::endl;
	float *voiArray = &reorder[0];

	ImageType::Pointer reorderedImage = converArray2Image(voiArray, imageSize, voxelSize);
	return reorderedImage;
}

//convert an array to an ITK image
ImageType::Pointer converArray2Image(float *imageArray, unsigned int* dim, float *voxelSize) {
	
	ImageType::Pointer finalImage;
	//write the array in image
	ImportFilterType::Pointer importFilter = ImportFilterType::New();
	ImportFilterType::SizeType size;
	//set image size
	size[0] = dim[0];
	size[1] = dim[1];
	size[2] = dim[2];


	int nrVoxelsPET = dim[0] * dim[1] * dim[2];
	//start with importing the image
	ImportFilterType::IndexType start;
	start.Fill(0);
	//set image region
	ImportFilterType::RegionType region;
	region.SetIndex(start);
	region.SetSize(size);
	importFilter->SetRegion(region);

	float centerX = float((voxelSize[0] * dim[0]) / 2.0);
	float centerY = float((voxelSize[1] * dim[1]) / 2.0);

	float centerZ = float((voxelSize[2] * dim[2]) / 2.0);
	const itk::SpacePrecisionType origin[3] = { centerX, centerY, -centerZ };
	const itk::SpacePrecisionType spacing[3] = { voxelSize[0], voxelSize[1], voxelSize[2] };
	//const itk::SpacePrecisionType direction[3] = { voxelSize[0], voxelSize[1], voxelSize[2] };
	//direction matrix in compliance with standard nifti files from UMCG
	/*using MatrixType = itk::Matrix<double, 3, 3>;
	MatrixType matrix;
	matrix[0][0] = 1;
	matrix[0][1] = 0;
	matrix[0][2] = 0.;

	matrix[1][0] = 0;
	matrix[1][1] = -1;
	matrix[1][2] = 0.;

	matrix[2][0] = 0;
	matrix[2][1] = 0;
	matrix[2][2] = 1;*/
	
	importFilter->SetOrigin(origin);
	importFilter->SetSpacing(spacing);
	//importFilter->SetDirection(matrix);
	const bool importImageFilterWillOwnTheBuffer = true;
	importFilter->SetImportPointer(imageArray, nrVoxelsPET, importImageFilterWillOwnTheBuffer);
	importFilter->Update();

	finalImage = importFilter->GetOutput();
	finalImage->Register();
	
	using FlipImageFilterType = itk::FlipImageFilter<ImageType>;
	FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
	flipFilter->SetInput(finalImage);
	FlipImageFilterType::FlipAxesArrayType flipAxes;
	flipAxes[0] = false;
	flipAxes[1] = true;
	flipAxes[2] = false;
	flipFilter->SetFlipAxes(flipAxes);
	flipFilter->Update();

	return flipFilter->GetOutput();
}




//the project file of the accurate tool is read in
ImageType::Pointer readPrjFilePET(string prjPath, string imageType, float smoothingKernel, unsigned int(&dimPET)[3], float(&voxelSize)[3]) {
	//help array to read in the integer byte
	unsigned int dimCT[3];
	float voxelSizeCT[3];
	float halfLife, volumeScale, frameScale;
	ImageType::Pointer PETimage;
	//read in prj file
	ifstream prjfile;
	prjfile.open(prjPath, ios::binary);
	if (prjfile.is_open() == false) {
		std::cout << "Cannot open project file\n";
	}
	else {
		//get the image dimensions, voxel sizes and half life of the tracer
		getImageDimension(prjfile, dimPET);
		getVoxelSize(prjfile, voxelSize);
		getHalfLifeScale(prjfile, halfLife, volumeScale, frameScale);
		getImageDimension(prjfile, dimCT);
		getVoxelSize(prjfile, voxelSizeCT);
		int nrVoxelsPET = int(dimPET[0] * dimPET[1] * dimPET[2]);
		vector<float> petVector(nrVoxelsPET);
		vector<float> ctVector(nrVoxelsPET);
		getImageValues(prjfile, petVector);
		//if image is PET image get array with PET values
		if (imageType == "PET") {
			float *arr = &petVector[0];
			PETimage = converArray2Image(arr, dimPET, voxelSize);
			PETimage->Update();
		}

		//else get arrat with CT values
		else if (imageType == "CT") {
			getImageValues(prjfile, ctVector);
			float *arr = &ctVector[0];
			PETimage = converArray2Image(arr, dimCT, voxelSizeCT);
			std::cout << "image PET" << std::endl;
			PETimage->Update();
		}
		else {
			std::cout << "Unknown image type. Please check the image type in the config file" << std::endl;
		}

		itk::FixedArray<bool, 3> flipAxes;
		flipAxes[0] = false;
		flipAxes[1] = false;
		flipAxes[2] = false;
		typedef itk::FlipImageFilter <ImageType> FlipImageFilterType;
		FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
		flipFilter->SetInput(PETimage);
		flipFilter->SetFlipAxes(flipAxes);
		flipFilter->Update();
		PETimage = flipFilter->GetOutput();
		PETimage->Update();
		if (smoothingKernel > 0) {

			PETimage = smoothImage(PETimage, smoothingKernel);
		}
		prjfile.close();
	}
	return PETimage;
}

//to do: read in as array
void getHalfLifeScale(ifstream &inFile, float &halfLife, float &volumeScale, float &frameScale) {
	inFile.read((char*)&halfLife, sizeof(float));
	inFile.read((char*)&volumeScale, sizeof(float));
	inFile.read((char*)&frameScale, sizeof(float));
}

void getImageValues(ifstream &inFile, vector<float> &imageArray) {
	for (int i = 0; i < boost::size(imageArray); i++) {
		inFile.read((char*)&imageArray[i], sizeof(float));
	}
}
//get the dimension of the image
void getImageDimension(ifstream &inFile, unsigned int(&dim)[3]) {


	short test[3];
	for (int i = 0; i < 3; i++) {
	
		
		inFile.read((char*)&test[i], sizeof(short));

	}
	dim[0] = int(test[0]);
	dim[1] = int(test[1]);
	dim[2] = int(test[2]);
}
//get voxel size of PET/CT image
void getVoxelSize(ifstream &inFile, float(&voxelSize)[3]) {
	for (int i = 0; i < 3; i++) {
		inFile.read((char*)&voxelSize[i], sizeof(float));
	}

}
