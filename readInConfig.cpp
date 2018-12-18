#include "readInFeatureSelection.h"
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "string"

using namespace std;





void readInFeatureSelection(EFoobar::Flags &featureFlags)
{
    std::string a = "1";

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("featureSelection.ini", pt);
    string statFeat = pt.get<std::string>("StatisticalFeatures.CalculateStatFeat");
    if(statFeat == a){
        featureFlags|= EFoobar::statisticalFeat;
    }
    string inthist= pt.get<std::string>("IntensityHistogramFeatures.CalculateIntensityHistogramFeat");
    if(inthist == a){
        featureFlags|= EFoobar::intensityHistFeat;
    }
    string glcmFeatures2DFullMerge = pt.get<std::string>("GLCMFeatures2DFullMerge.CalculateGLCMFeatures2DFullMerge");
    if(glcmFeatures2DFullMerge == a){
        featureFlags |= EFoobar::glcmFeat2DFullMerge;
    }
    string glcmFeatures2DWMerge = pt.get<std::string>("GLCMFeatures2DWMerge.CalculateGLCMFeatures2DWMerge");
    if(glcmFeatures2DWMerge == a){
        featureFlags |= EFoobar::glcmFeat2DWMerge;
    }
    string glcmFeatures2DWOMerge = pt.get<std::string>("GLCMFeatures2DWOMerge.CalculateGLCMFeatures2DWOMerge");
    if(glcmFeatures2DWOMerge == a){
        featureFlags |= EFoobar::glcmFeat2DWOMerge;
    }
    string glcmFeatures3DWMerge = pt.get<std::string>("GLCMFeatures3DWMerge.CalculateGLCMFeatures3DWMerge");
    if(glcmFeatures3DWMerge == a){
        featureFlags |= EFoobar::glcmFeat3DWMerge;
    }
    string glcmFeatures3DWOMerge = pt.get<std::string>("GLCMFeatures3DWOMerge.CalculateGLCMFeatures3DWOMerge");
    if(glcmFeatures3DWMerge == a){
        featureFlags |= EFoobar::glcmFeat3DWOMerge;
    }
    string glrlmFeatures2DFullMerge=pt.get<std::string>("GLRLMFeatures2DFullMerge.CalculateGLRLMFeatures2DFullMerge");
    if(glrlmFeatures2DFullMerge == a){
        featureFlags |= EFoobar::glrlmFeat2DFullMerge;
    }
    string glrlmFeatures2DWMerge=pt.get<std::string>("GLRLMFeatures2DWMerge.CalculateGLRLMFeatures2DWMerge");
    if(glrlmFeatures2DWMerge == a){
        featureFlags |= EFoobar::glrlmFeat2DWMerge;
    }
    string glrlmFeatures2DWOMerge=pt.get<std::string>("GLRLMFeatures2DWOMerge.CalculateGLRLMFeatures2DWOMerge");
    if(glrlmFeatures2DWOMerge == a){
        featureFlags |= EFoobar::glrlmFeat2DWOMerge;
    }
    string glrlmFeatures3DWMerge=pt.get<std::string>("GLRLMFeatures3DWMerge.CalculateGLRLMFeatures3DWMerge");
    if(glrlmFeatures3DWMerge == a){
        featureFlags |= EFoobar::glrlmFeat3DWMerge;
    }
    string glrlmFeatures3DWOMerge=pt.get<std::string>("GLRLMFeatures3DWOMerge.CalculateGLRLMFeatures3DWOMerge");
    if(glrlmFeatures3DWOMerge == a){
        featureFlags |= EFoobar::glrlmFeat3DWOMerge;
    }
    string gldzmFeatures2D=pt.get<std::string>("GLDZMFeatures2D.CalculateGLDZMFeatures2D");
    if(gldzmFeatures2D == a){
        featureFlags |= EFoobar::gldzmFeat2D;
    }
    string gldzmFeatures2DWOMerge=pt.get<std::string>("GLDZMFeatures2DWOMerge.CalculateGLDZMFeatures2DWOMerge");
    if(gldzmFeatures2DWOMerge == a){
        featureFlags |= EFoobar::gldzmFeat2DWOMerge;
    }
    string gldzmFeatures3D=pt.get<std::string>("GLDZMFeatures3D.CalculateGLDZMFeatures3D");
    if(gldzmFeatures3D == a){
        featureFlags |= EFoobar::gldzmFeat3D;
    }

    string glszmFeatures2DWMerge=pt.get<std::string>("GLSZMFeatures2DWMerge.CalculateGLSZMFeatures2DWMerge");
    if(glszmFeatures2DWMerge == a){
        featureFlags |= EFoobar::glszmFeat2DWMerge;
    }
    string glszmFeatures2DWOMerge=pt.get<std::string>("GLSZMFeatures2DWOMerge.CalculateGLSZMFeatures2DWOMerge");
    if(glszmFeatures2DWOMerge == a){
        featureFlags |= EFoobar::glszmFeat2DWOMerge;
    }
    string glszmFeatures3D=pt.get<std::string>("GLSZMFeatures3D.CalculateGLSZMFeatures3D");
    if(glszmFeatures3D== a){
        featureFlags |= EFoobar::glszmFeat3D;
    }

    string ngtdmFeatures2DWMerge=pt.get<std::string>("NGTDMFeatures2DWMerge.CalculateNGTDMFeatures2DWMerge");
    if(ngtdmFeatures2DWMerge== a){
        featureFlags |= EFoobar::ngtdmFeat2DWMerge;
    }
    string ngtdmFeatures2DWOMerge=pt.get<std::string>("NGTDMFeatures2DWOMerge.CalculateNGTDMFeatures2DWOMerge");
    if(ngtdmFeatures2DWOMerge== a){
        featureFlags |= EFoobar::ngtdmFeat2DWOMerge;
    }
    string ngtdmFeatures3D=pt.get<std::string>("NGTDMFeatures3D.CalculateNGTDMFeatures3D");
    if(ngtdmFeatures3D == a){
        featureFlags |= EFoobar::ngtdmFeat3D;
    }
    string intVolFeatures = pt.get<std::string>("IntensityVolume.CalculateIntensityVolume");
    if(intVolFeatures == a){
        featureFlags |= EFoobar::intVolFeat;
    }

}

void CalculateRelFeatures(EFoobar::Flags flags, boost::multi_array<float, 3> A, vector<float> vectorMatrElement, vector<float> diffGreyLevels)
{

    if (flags & EFoobar::statisticalFeat){
        StatisticalFeatures<float, 3> statFeatures;
        statFeatures.calculateAllStatFeatures(statFeatures, vectorMatrElement);
        statFeatures.writeCSVFileStatistic(statFeatures);
    }
//    if (flags & EFoobar::intensityHistFeat){
//        IntensityHistogram<float, 3> intensityHist;
//        intensityHist.calculateAllIntFeatures(intensityHist, A, vectorMatrElement, diffGreyLevels);
//        intensityHist.writeCSVFileIntensity(intensityHist);
//    }


}

void calculateRelFeaturesDiscretized(EFoobar::Flags flags, boost::multi_array<float, 3> A, vector<float> vectorMatrElement, vector<float> diffGreyLevels)
{

    if (flags & EFoobar::intensityHistFeat){
        IntensityHistogram<float, 3> intensityHist;
        intensityHist.calculateAllIntFeatures(intensityHist, A, vectorMatrElement, diffGreyLevels);
        intensityHist.writeCSVFileIntensity(intensityHist);
    }
    float maxIntensity = *max_element(vectorMatrElement.begin(), vectorMatrElement.end());
    if (flags & EFoobar::glcmFeat2DFullMerge){
        GLCMFeatures2DFullMerge<float, 3> glcm2DFullMerge;
        glcm2DFullMerge.calculateAllGLCMFeatures2DFullMerge(glcm2DFullMerge, A, maxIntensity);
        glcm2DFullMerge.writeCSVFileGLCM2DFullMerge(glcm2DFullMerge);
    }
    if (flags & EFoobar::glcmFeat2DWMerge){
        GLCMFeatures2DWMerge<float, 3> glcm2DWMerge;
        glcm2DWMerge.calculateAllGLCMFeatures2DWMerge(glcm2DWMerge, A, maxIntensity);
        glcm2DWMerge.writeCSVFileGLCM2DWMerge(glcm2DWMerge);
    }
    if (flags & EFoobar::glcmFeat2DWOMerge){
        GLCMFeatures2DWOMerge<float, 3> glcm2DWOMerge;
        glcm2DWOMerge.calculateAllGLCMFeatures2DWOMerge(glcm2DWOMerge, A, maxIntensity);
        glcm2DWOMerge.writeCSVFileGLCM2DWOMerge(glcm2DWOMerge);
    }
    if (flags & EFoobar::glcmFeat3DWMerge){
        GLCMFeatures3DWMerge<float, 3> glcm3DWMerge;
        glcm3DWMerge.calculateAllGLCMFeatures3DWMerge(glcm3DWMerge, A, maxIntensity);
        glcm3DWMerge.writeCSVFileGLCM3DWMerge(glcm3DWMerge);
    }
    if (flags & EFoobar::glcmFeat3DWOMerge){
        GLCMFeatures3DWOMerge<float, 3> glcmFeat3DWOMerge;
        glcmFeat3DWOMerge.calculateAllGLCMFeatures3DWOMerge(glcmFeat3DWOMerge, A, maxIntensity);
        glcmFeat3DWOMerge.writeCSVFileGLCM3D(glcmFeat3DWOMerge);
    }

    if (flags & EFoobar::glrlmFeat2DFullMerge){
        GLRLMFeatures2DFullMerge<float, 3> glrlm2DFullMerge;
        glrlm2DFullMerge.calculateAllGLRLMFeatures2DFullMerge(glrlm2DFullMerge, A, diffGreyLevels, vectorMatrElement);
        glrlm2DFullMerge.writeCSVFileGLRLM2DFullMerge(glrlm2DFullMerge);
    }
    if (flags & EFoobar::glrlmFeat2DWMerge){
        GLRLMFeatures2DWMerge<float, 3> glrlm2DWMerge;
        glrlm2DWMerge.calculateAllGLRLMFeatures2DWMerge(glrlm2DWMerge, A, diffGreyLevels);
        glrlm2DWMerge.writeCSVFileGLRLM2DWMerge(glrlm2DWMerge);
    }
    if (flags & EFoobar::glrlmFeat2DWOMerge){
        GLRLMFeatures2DWOMerge<float, 3> glrlm2DWOMerge;
        glrlm2DWOMerge.calculateAllGLRLMFeatures2DWOMerge(glrlm2DWOMerge, A, diffGreyLevels);
        glrlm2DWOMerge.writeCSVFileGLRLM2DWOMerge(glrlm2DWOMerge);
    }
    if (flags & EFoobar::glrlmFeat3DWMerge){
        GLRLMFeatures3D<float, 3> glrlm3DWMerge;
        glrlm3DWMerge.calculateAllGLRLMFeatures3D(glrlm3DWMerge, A, diffGreyLevels, vectorMatrElement);
        glrlm3DWMerge.writeCSVFileGLRLM3D(glrlm3DWMerge);
    }
    if (flags & EFoobar::glrlmFeat3DWOMerge){
        GLRLMFeatures3DWOMerge<float, 3> glrlm3DWOMerge;
        glrlm3DWOMerge.calculateAllGLRLMFeatures3DWOMerge(glrlm3DWOMerge, A, diffGreyLevels, vectorMatrElement);
        glrlm3DWOMerge.writeCSVFileGLRLM3DWOMerge(glrlm3DWOMerge);
    }
    if (flags & EFoobar::gldzmFeat2D){
        GLDZMFeatures2D<float, 3> gldzm2D;
        gldzm2D.calculateAllGLDZMFeatures2D(gldzm2D, A, diffGreyLevels, vectorMatrElement);
        gldzm2D.writeCSVFileGLDZM(gldzm2D);
    }
//
    if (flags & EFoobar::gldzmFeat2DWOMerge){
        GLDZMFeatures2DWOMerge<float, 3> gldzm2DWOMerge;
        gldzm2DWOMerge.calculateAllGLDZMFeatures2DWOMerge(gldzm2DWOMerge, A, diffGreyLevels);
        gldzm2DWOMerge.writeCSVFileGLDZM2DWOMerge(gldzm2DWOMerge);
    }
    if (flags & EFoobar::gldzmFeat3D){
        GLDZMFeatures3D<float, 3> gldzm3D;
        gldzm3D.calculateAllGLDZMFeatures3D(gldzm3D, A, diffGreyLevels, vectorMatrElement);
        gldzm3D.writeCSVFileGLDZM3D(gldzm3D);
    }
    if (flags & EFoobar::glszmFeat2DWMerge){
        GLSZMFeatures2D<float, 3> glszm2D;
        glszm2D.calculateAllGLSZMFeatures2D(glszm2D, A, diffGreyLevels, vectorMatrElement);
        glszm2D.writeCSVFileGLSZM(glszm2D);
    }
//
    if (flags & EFoobar::glszmFeat2DWOMerge){
        GLSZMFeatures2DWOMerge<float, 3> glszm2DWOMerge;
        glszm2DWOMerge.calculateAllGLSZMFeatures2DWOMerge(glszm2DWOMerge, A, diffGreyLevels);
        glszm2DWOMerge.writeCSVFileGLSZM(glszm2DWOMerge);
    }
    if (flags & EFoobar::glszmFeat3D){
        GLSZMFeatures3D<float, 3> glszm3D;
        glszm3D.calculateAllGLSZMFeatures3D(glszm3D, A, diffGreyLevels, vectorMatrElement);
        glszm3D.writeCSVFileGLSZM3D(glszm3D);
    }
////
    if (flags & EFoobar::ngldmFeat2DWMerge){
        NGLDMFeatures<float, 3> ngldm2DWMerge;
        ngldm2DWMerge.calculateAllNGLDMFeatures2D(ngldm2DWMerge, A, diffGreyLevels, vectorMatrElement, 0, 1);
        ngldm2DWMerge.writeCSVFileNGLDM(ngldm2DWMerge);
    }
    if (flags & EFoobar::ngldmFeat2DWOMerge){
        NGLDMFeatures2DWOMerge<float, 3> ngldm2DWOMerge;
        ngldm2DWOMerge.calculateAllNGLDMFeatures2DWOMerge2D(ngldm2DWOMerge, A, diffGreyLevels, 0, 1);
        ngldm2DWOMerge.writeCSVFileNGLDM(ngldm2DWOMerge);
    }
    if (flags & EFoobar::ngldmFeat3D){
        NGLDMFeatures3D<float, 3> ngldm3D;
        ngldm3D.calculateAllNGLDMFeatures3D(ngldm3D, A, diffGreyLevels, vectorMatrElement, 0, 1);
        ngldm3D.writeCSVFileNGLDM3D(ngldm3D);
    }
    if (flags & EFoobar::ngtdmFeat2DWMerge){
        NGTDMFeatures<float, 3> ngtdm2DWMerge;
        ngtdm2DWMerge.calculateAllNGTDMFeatures(ngtdm2DWMerge, A, diffGreyLevels);
        ngtdm2DWMerge.writeCSVFileNGTDM(ngtdm2DWMerge);
    }
    if (flags & EFoobar::ngtdmFeat2DWOMerge){
        NGTDM2DWOMerge<float, 3> NGTDM2DWOMerge;
        NGTDM2DWOMerge.calculateAllNGTDMFeatures2DWOMerge(NGTDM2DWOMerge, A, diffGreyLevels);
        NGTDM2DWOMerge.writeCSVFileNGTDM2DWOMerge(NGTDM2DWOMerge);
    }
    if (flags & EFoobar::ngtdmFeat3D){
        NGTDMFeatures3D<float, 3> ngtdm3D;
        ngtdm3D.calculateAllNGTDMFeatures3D(ngtdm3D, A, diffGreyLevels);
        ngtdm3D.writeCSVFileNGTDM3D(ngtdm3D);
    }
    if (flags & EFoobar::intVolFeat){
        IntensityVolumeFeatures<float, 3> intVol;
        intVol.calculateAllIntensVolFeatures(intVol, A, vectorMatrElement);
        intVol.writeCSVFileIntVol(intVol);
    }
}
