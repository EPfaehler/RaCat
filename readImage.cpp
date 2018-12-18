//#include "readImage.h"
//
//void readImage(){
//    //determine what type the pixels have and which the dimension the image has
//    //typedef itk::Image< unsigned short, 3 > ImageType;
//    int dimension =3;
//
//    //create image
//    if(dimension == 2){
//
//        //instantiate image reader class
//        typedef itk::ImageFileReader< ImageType2D > ReaderType;
//        //create reader object
//        ReaderType::Pointer reader = ReaderType::New();
//        //set filename
//        reader->SetFileName("OSEM60M256H0_1_10_0.dcm");
//        ImageType2D::Pointer image = reader->GetOutput();
//        reader->Update();
//typename ImageType2D::RegionType wholeRegion = image->GetLargestPossibleRegion();
//    typename ImageType2D::SizeType imageSize = wholeRegion.GetSize();
//    int dede = ima
//    int dede = wholeRegion.Ge
//    //get the size of the image
//        get2Dimage(image);
//    }
//
//    else if(dimension == 3){
//        //instantiate image reader class
//        typedef itk::ImageFileReader< ImageType3D > ReaderType;
//        //create reader object
//        ReaderType::Pointer reader = ReaderType::New();
//        //set filename
//        reader->SetFileName("OSEM60M256H0_1_10_0.dcm");
//        ImageType3D::Pointer image = reader->GetOutput();
//        reader->Update();
//        get3Dimage(image);
//    }
////
////Const versions of iterators may read, but may not write pixel values
////typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
////typedef itk::ImageRegionIterator< ImageType > IteratorType;
//////
////ConstIteratorType constIterator( image, image->GetLargestPossibleRegion() );
////IteratorType out( image, image->GetRequestedRegion() );
////
////typedef boost::multi_array<double, 3> array_type;
////typedef array_type::index index;
////array_type A(boost::extents[imageSize[0]][imageSize[1]][imageSize[2]]);
////std::cout<<imageSize.GetSizeDimension();
//////std::cout<<imageSize[1];
//////std::cout<<imageSize[2];
////int k = 0;
////int i = 0;
////for ( constIterator.GoToBegin(); !constIterator.IsAtEnd(); ++constIterator)
////{
////
////   if(i<imageSize[0]){
////        A[i][k]=constIterator.Get();
////        i=i+1;
////    }
////    if(i==imageSize[0]){
////        k=k+1;
////        i=0;
////    }
////
////}
//
//}
//
//void get3Dimage(const ImageType3D* image){
//    ImageType3D::RegionType wholeRegion = image->GetLargestPossibleRegion();
//    ImageType3D::SizeType imageSize = wholeRegion.GetSize();
//    ConstIteratorType3D constIterator3D( image, image->GetLargestPossibleRegion()) ;
//    typedef boost::multi_array<double, 3> array_type;
//    typedef array_type::index index;
//    array_type A(boost::extents[imageSize[0]][imageSize[1]][imageSize[2]]);
////    //std::cout<<imageSize[1];
////    //std::cout<<imageSize[2];
//    int k = 0;
//    int i = 0;
//    int j = 0;
//    std::cout<<imageSize[1];
//    for ( constIterator3D.GoToBegin(); !constIterator3D.IsAtEnd(); ++constIterator3D)
//    {
//
//       if(i<imageSize[0]){
//            A[i][j][k]=constIterator3D.Get();
//            i=i+1;
//        }
//
//        if(i==imageSize[0] && j<imageSize[1]-1){
//            j=j+1;
//            i=0;
//        }
//        if(j==imageSize[1]){
//            j=0;
//            k=k+1;
//        }
//
//    }
//
//}
//
//void get2Dimage(const ImageType2D* image){
//    ImageType2D::RegionType wholeRegion = image->GetLargestPossibleRegion();
//    ImageType2D::SizeType imageSize = wholeRegion.GetSize();
//    ConstIteratorType2D constIterator2D( image, image->GetLargestPossibleRegion()) ;
//    typedef boost::multi_array<double, 2> array_type;
//    typedef array_type::index index;
//    array_type A(boost::extents[imageSize[0]][imageSize[1]]);
//    //std::cout<<imageSize[1];
//    //std::cout<<imageSize[2];
//    int k = 0;
//    int i = 0;
//    for ( constIterator2D.GoToBegin(); !constIterator2D.IsAtEnd(); ++constIterator2D)
//    {
//
//       if(i<imageSize[0]){
//            A[i][k]=constIterator2D.Get();
//            i=i+1;
//        }
//        if(i==imageSize[0]){
//            k=k+1;
//            i=0;
//        }
//
//    }
//
//}
