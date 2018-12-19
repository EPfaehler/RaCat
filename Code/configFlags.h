
#include "boost/multi_array.hpp"
#include "intensityHistogram.h"
#include "morphologicalFeatures.h"
#include "localIntensityFeatures.h"

#include "GLCMFeatures2DVMRG.h"
#include "GLCMFeatures2DMRG.h"
#include "GLCMFeatures2DAVG.h"
#include "GLCMFeatures2DDMRG.h"
#include "GLCMFeatures3DMRG.h"
#include "GLCMFeatures3DAVG.h"

#include "GLRLMFeatures2DVMRG.h"
#include "GLRLMFeatures2DDMRG.h"
#include "GLRLMFeatures2DMRG.h"
#include "GLRLMFeatures2DAVG.h"
#include "GLRLMFeatures3D.h"
#include "GLRLMFeatures3DAVG.h"

#include "GLSZMFeatures2D.h"
#include "GLSZMFeatures2DAVG.h"
#include "GLSZMFeatures3D.h"

#include "NGLDMFeatures2DMRG.h"
#include "NGLDMFeatures2DAVG.h"
#include "NGLDMFeatures3D.h"

#include "NGTDM2DMRG.h"
#include "NGTDM2DAVG.h"
#include "NGTDM3D.h"

#include "GLDZMFeatures2DMRG.h"
#include "GLDZMFeatures3D.h"
#include "GLDZMFeatures2DAVG.h"

#include "intensityVolumeFeatures.h"
#include "localIntensityFeatures.h"


//namespace EFoobar
//{
    enum Flags
    {
        statisticalFeat,
        intensityHistFeat    = 0x2,

        glcmFeat2DFullMerge = 0x3,
        glcmFeat2DWMerge = 0x4,
        glcmFeat2DWOMerge = 0x5,
		glcmFeat2DDMerge = 0x6,
        glcmFeat3DWMerge = 0x7,
        glcmFeat3DWOMerge = 0x8,

        glrlmFeat2DFullMerge = 0x9,
        glrlmFeat2DWMerge = 0x10,
        glrlmFeat2DWOMerge = 0x11,
        glrlmFeat3DWMerge = 0x12,
        glrlmFeat3DWOMerge = 0x13,

        gldzmFeat2D = 0x14,
        gldzmFeat2DWOMerge = 0x15,
        gldzmFeat3D = 0x16,

        glszmFeat2DWMerge = 0x17,
        glszmFeat2DWOMerge = 0x18,
        glszmFeat3D = 0x19,

        ngldmFeat2DWMerge = 0x20,
        ngldmFeat2DWOMerge = 0x21,
        ngldmFeat3D = 0x22,

        ngtdmFeat2DWMerge = 0x23,
        ngtdmFeat2DWOMerge = 0x24,
        ngtdmFeat3D = 0x25,

        intVolFeat = 0x26,
		morphologicalFeat = 0x27,
		localIntensityFeat = 0x28,

        FB_C    = 0x11,
    };
   
//}






