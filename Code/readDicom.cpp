//read Dicom images from a folder
ImageType::Pointer readDicom(string path) {
	ImageType::Pointer finalImage;
	typedef itk::ImageSeriesReader< ImageType >        ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	typedef itk::GDCMImageIO       ImageIOType;
	ImageIOType::Pointer dicomIO = ImageIOType::New();
	reader->SetImageIO(dicomIO);
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	nameGenerator->SetUseSeriesDetails(true);
	nameGenerator->AddSeriesRestriction("0008|0021");
	nameGenerator->SetDirectory(path);
	using SeriesIdContainer = std::vector< std::string >;
	const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
	SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
	SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
	while (seriesItr != seriesEnd)
	{

		std::cout << seriesItr->c_str() << std::endl;
		++seriesItr;
	}
	std::string seriesIdentifier;
	seriesIdentifier = seriesUID.begin()->c_str();

	using FileNamesContainer = std::vector< std::string >;
	FileNamesContainer fileNames;
	fileNames = nameGenerator->GetFileNames(seriesIdentifier);
	reader->SetFileNames(fileNames);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
	}
	finalImage = reader->GetOutput();
	finalImage->Update();
	return finalImage;
}

//PolygonType::Pointer polygon;
//PolygonPointType p;

//GroupType::Pointer group = GroupType::New();

//PointType point;
//ITK spatial object to image filter segmentation definitions
//PolygonType::PointListType pointList;
//ImageSliceType::Pointer temp2Dimage = ImageSliceType::New();
ImageType::Pointer resetImage(ImageType *mask)
{
	ImageType::Pointer helpImg = mask;
	itk::ImageRegionIterator<ImageType> imageIt(helpImg, helpImg->GetLargestPossibleRegion());
	imageIt.GoToBegin();
	while (!imageIt.IsAtEnd())
	{
		imageIt.Set(0);
		++imageIt;
	}
	return helpImg;
}



//remove empty spaces from contour names
void trim(std::string &str)
{
	std::string temp;
	for (unsigned int i = 0; i < str.length(); i++)
		if (str[i] != ' ') temp += str[i];
	str = temp;
}
void reset2DImage(ImageSliceType::Pointer imageSlice)
{
	itk::ImageRegionIterator<ImageSliceType> imageIt(imageSlice, imageSlice->GetLargestPossibleRegion());
	imageIt.GoToBegin();
	while (!imageIt.IsAtEnd()) {
		imageIt.Set(0);
		++imageIt;
	}
}


void mergeImages(ImageSliceType::Pointer tempSlice, ImageType::Pointer finalImage, int iRequiredSlice) {
	float pixelValue = 0;
	ImageType::IndexType pixelIndex;
	ImageSliceType::IndexType sliceIndex;
	int iX = finalImage->GetLargestPossibleRegion().GetSize()[0];
	int iY = finalImage->GetLargestPossibleRegion().GetSize()[1];
	if (iRequiredSlice>0) {
		for (int i = 0; i<iX; i++)
			for (int j = 0; j<iY; j++)
			{
				pixelIndex[0] = i;
				pixelIndex[1] = j;
				sliceIndex[0] = i;
				sliceIndex[1] = j;

				pixelValue = tempSlice->GetPixel(sliceIndex);
				pixelIndex[2] = iRequiredSlice;

				//Disable hole filling (if required please uncomment the next line (and comment the following line)).  
				//if (pixelValue != 0)  finalImage->SetPixel(pixelIndex, pixelValue  );
				finalImage->SetPixel(pixelIndex, pow(finalImage->GetPixel(pixelIndex), pixelValue));
			}
	}
}



ImageType::Pointer readRTstruct(string voiPath, ImageType *imageDicom, ImageType *mask) {
	mask = imageDicom;
	mask->Update();
	const typename ImageType::RegionType region = imageDicom->GetBufferedRegion();
	resetImage(mask);
	//std::string templateFilename;
	//float dcmOrigin = 0;
	////take single slices, to set the image region to 1 slice by slice
	//ImageSliceType::Pointer temp2Dimage = ImageSliceType::New();
	////get a polygon
	//polygon = PolygonType::New();
	//ImageType::IndexType pixelIndex;
	//int iCurrentSlice = 0;
	//int iPointsOutsideBoundary = 0;
	//gdcm::Reader RTreader;
	//RTreader.SetFileName(voiPath.c_str());
	//if (!RTreader.Read())
	//{
	//	std::cout << "Problem reading file: " << voiPath << std::endl;
	//	return 0;
	//}
	////we need to create a temporary 2D slice as well...

	//ImageType::RegionType inputRegion = mask->GetLargestPossibleRegion();
	//FilterTypeRTS::Pointer filter = FilterTypeRTS::New();
	//ImageType::SizeType size = inputRegion.GetSize();
	////start at slice 0
	//size[2] = 0;
	//ImageType::IndexType start = inputRegion.GetIndex();
	//start[2] = 0;
	//ImageType::RegionType desiredRegion;
	//desiredRegion.SetSize(size);
	//desiredRegion.SetIndex(start);
	////get only one slice
	//filter->SetDirectionCollapseToIdentity();
	//filter->SetExtractionRegion(desiredRegion);
	//filter->SetInput(mask);
	//filter->Update();
	//temp2Dimage = filter->GetOutput();
	////get spacing of 2D image
	//const typename ImageType2D::SpacingType& inputSpacing = temp2Dimage->GetSpacing();
	////get the datasets of RT struct
	//const gdcm::DataSet& ds = RTreader.GetFile().GetDataSet();
	//std::cout << "Parsing: " << voiPath << std::endl;
	//gdcm::MediaStorage ms;
	//ms.SetFromFile(RTreader.GetFile());
	//std::cout << "media storage: " << ms << std::endl;
	//// (3006,0020) SQ (Sequence with explicit length #=4)      # 370, 1 StructureSetROISequence  
	//gdcm::Tag tssroisq(0x3006, 0x0020);
	//if (!ds.FindDataElement(tssroisq)) {
	//	std::cout << "Problem locating 0x3006,0x0020 - Is this a valid RT Struct file?" << std::endl;
	//	return 0;
	//}
	//gdcm::Tag troicsq(0x3006, 0x0039);
	//if (!ds.FindDataElement(troicsq)) {
	//	std::cout << "Problem locating 0x3006,0x0039 - Is this a valid RT Struct file?" << std::endl;
	//	return 0;
	//}
	//const gdcm::DataElement &roicsq = ds.GetDataElement(troicsq);

	//gdcm::SmartPointer<gdcm::SequenceOfItems> sqi = roicsq.GetValueAsSQ();
	//if (!sqi || !sqi->GetNumberOfItems()) {
	//	return 0;
	//}
	//const gdcm::DataElement &ssroisq = ds.GetDataElement(tssroisq);
	//gdcm::SmartPointer<gdcm::SequenceOfItems> ssqi = ssroisq.GetValueAsSQ();
	////get the number of structures in file
	//if (!ssqi || !ssqi->GetNumberOfItems()) {
	//	return 0;
	//}
	//std::cout << "Number of structures found:" << sqi->GetNumberOfItems() << std::endl;
	//if (sqi->GetNumberOfItems() > 1) {
	//	std::cout << "More than one structure was found in the RT struct file" << std::endl;
	//	std::cout << "Please change the RT struct so that there is only one tumor included" << std::endl;
	//	exit(EXIT_FAILURE);
	//}
	////loop through structures
	//for (unsigned int pd = 0; pd < sqi->GetNumberOfItems(); ++pd) {
	//	const gdcm::Item & item = sqi->GetItem(pd + 1); // Item start at #1
	//	gdcm::Attribute<0x3006, 0x0084> roinumber;
	//	const gdcm::DataSet& nestedds = item.GetNestedDataSet();
	//	roinumber.SetFromDataElement(nestedds.GetDataElement(roinumber.GetTag()));

	//	// find structure_set_roi_sequence corresponding to roi_contour_sequence (by comparing id numbers)
	//	unsigned int spd = 0;
	//	gdcm::Item & sitem = ssqi->GetItem(spd + 1);
	//	gdcm::DataSet& snestedds = sitem.GetNestedDataSet();
	//	gdcm::Attribute<0x3006, 0x0022> sroinumber;
	//	do {
	//		sitem = ssqi->GetItem(spd + 1);
	//		snestedds = sitem.GetNestedDataSet();
	//		sroinumber.SetFromDataElement(snestedds.GetDataElement(sroinumber.GetTag()));
	//		spd++;
	//	} while (sroinumber.GetValue() != roinumber.GetValue());

	//	gdcm::Tag stcsq(0x3006, 0x0026);
	//	if (!snestedds.FindDataElement(stcsq)) {
	//		std::cout << "Did not find sttsq data el " << stcsq << "   continuing..." << std::endl;
	//		continue; //return 0;
	//	}
	//	const gdcm::DataElement &sde = snestedds.GetDataElement(stcsq);

	//	//(3006,002a) IS [255\192\96]                              # 10,3 ROI Display Color
	//	gdcm::Tag troidc(0x3006, 0x002a);
	//	gdcm::Attribute<0x3006, 0x002a> color = {};
	//	if (nestedds.FindDataElement(troidc)) {
	//		const gdcm::DataElement &decolor = nestedds.GetDataElement(troidc);
	//		color.SetFromDataElement(decolor);
	//	}
	//	//(3006,0040) SQ (Sequence with explicit length #=8)      # 4326, 1 ContourSequence
	//	gdcm::Tag tcsq(0x3006, 0x0040);
	//	if (!nestedds.FindDataElement(tcsq)) {
	//		continue;
	//	}
	//	const gdcm::DataElement& csq = nestedds.GetDataElement(tcsq);

	//	gdcm::SmartPointer<gdcm::SequenceOfItems> sqi2 = csq.GetValueAsSQ();
	//	if (!sqi2 || !sqi2->GetNumberOfItems()) {
	//		std::cout << "csq: " << csq << std::endl;
	//		std::cout << "sqi2: " << *sqi2 << std::endl;
	//		std::cout << "Did not find sqi2 or no. items == 0   " << sqi2->GetNumberOfItems() << "   continuing..." << std::endl;
	//		continue;
	//	}
	//	int test = 0;
	//	//get number of regions in structure
	//	unsigned int nitems = sqi2->GetNumberOfItems();
	//	std::cout << "Structure " << pd << ". Number of regions: " << nitems << std::endl;
	//	std::string str_currentOrgan(sde.GetByteValue()->GetPointer(), sde.GetByteValue()->GetLength());
	//	//now loop through each item for this structure (eg one prostate region on a single slice is an item)
	//	for (unsigned int i = 0; i < nitems; ++i) {
	//		const gdcm::Item & item2 = sqi2->GetItem(i + 1); // Item start at #1
	//		const gdcm::DataSet& nestedds2 = item2.GetNestedDataSet();
	//		// (3006,0050) DS [43.57636\65.52504\-10.0\46.043102\62.564945\-10.0\49.126537\60.714... # 398,48 ContourData
	//		gdcm::Tag tcontourdata(0x3006, 0x0050);
	//		const gdcm::DataElement & contourdata = nestedds2.GetDataElement(tcontourdata);

	//		gdcm::Attribute<0x3006, 0x0050> at;
	//		at.SetFromDataElement(contourdata);
	//		const double* pts = at.GetValues();
	//		unsigned int npts = at.GetNumberOfValues() / 3;
	//		for (unsigned int j = 0; j < npts * 3; j += 3) {
	//			point[0] = pts[j + 0];
	//			point[1] = pts[j + 1];
	//			point[2] = pts[j + 2];
	//			//transform points to image co-ordinates
	//			if (!(mask->TransformPhysicalPointToIndex(point, pixelIndex))) {
	//				//Are there points outside the image boundary.  This may occur with automatically segmented objects such as benches or external body outlines?  
	//				iPointsOutsideBoundary++;
	//			}
	//			p.SetPosition(pixelIndex[0], pixelIndex[1], pixelIndex[2]);
	//			itk::Index<2> test;
	//			test[0] = pixelIndex[0];
	//			test[1] = pixelIndex[1];
	//			temp2Dimage->SetPixel(test, 1);
	//			p.SetRed(1);
	//			p.SetBlue(1);
	//			p.SetGreen(1);
	//			pointList.push_back(p);
	//		}

	//		// we have the points for a contour in a single slice.  We need to join these up and insert into the slice as polygon.
	//		iCurrentSlice = pixelIndex[2];
	//		SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
	//		//inserting points into slices
	//		reset2DImage(temp2Dimage);
	//		//need to create a 2D slice here, put the polygon on it, and insert it back into the 3D volume...  
	//		group->AddSpatialObject(polygon); //add a new polygon group

	//		try
	//		{
	//			itk::Point<float, 2> newP;


	//			polygon->SetPoints(pointList);  //so copy them to a polygon object


	//											/*imageFilter->SetInput(group);

	//											imageFilter->SetInsideValue(0);
	//											imageFilter->SetOutsideValue(1);
	//											imageFilter->SetSize(temp2Dimage->GetLargestPossibleRegion().GetSize());
	//											imageFilter->Update();
	//											temp2Dimage = imageFilter->GetOutput();*/
	//			itk::ImageRegionIterator<ImageType2D> it(temp2Dimage, temp2Dimage->GetLargestPossibleRegion());
	//			it.GoToBegin();
	//			while (!it.IsAtEnd())
	//			{
	//				newP[0] = float(it.GetIndex()[0]) + inputSpacing[0] * 0.5;
	//				newP[1] = float(it.GetIndex()[1]) + inputSpacing[1] * 0.5;
	//				if (polygon->IsInside(newP) == 0) {
	//					it.Set(1);
	//				}

	//				++it;
	//			}

	//		}
	//		catch (itk::ExceptionObject & err)
	//		{
	//			std::cerr << "Problem setting polygon->SetPoints for this region (non-planar)" << std::endl;
	//			std::cerr << err << std::endl;

	//		}

	//		//merge new polygon from temp image into the contour image
	//		mergeImages(temp2Dimage, mask, iCurrentSlice);

	//		//remove the polygon and clean up pointlist
	//		group->RemoveSpatialObject(polygon);
	//		pointList.clear();
	//	}

	//	if (iPointsOutsideBoundary > 0) {
	//		std::cout << " --" << iPointsOutsideBoundary << " contour points detected outside image boundary. Please check the output volume. ";
	//		iPointsOutsideBoundary = 0;
	//	}

	//}

	return mask;
}