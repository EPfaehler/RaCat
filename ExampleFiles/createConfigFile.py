# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 15:35:07 2018

@author: pfaehlereag
"""
import os

#create a config file 
#the user can set the location where the gile will be saved, name of the config file
#and required number of bns
def createConfigFile(fileLocation, configName, nrBins):
    #generate configFile
    iniFile = open(fileLocation + '/' + configName + '.ini', 'w')
    iniFile.write('[Smoothing]')
    iniFile.write('\n')
    iniFile.write('SmoothingKernel = 0\n')
    iniFile.write('\n')
    iniFile.write('[ThresholdForVOI]')
    iniFile.write('\n')
    iniFile.write('threshold = 0.5\n')
    iniFile.write('\n')
    iniFile.write('[Discretization]\n')
    iniFile.write('UseFixedBinWidth = 0\n')
    iniFile.write('BinWidth = 0.25\n')
    iniFile.write('UseFixedNrBins = 1\n')
    #here the required number of bins is set
    iniFile.write('NrBins = %s \n\n' %nrBins)
    iniFile.write('[DiscretizationIVH]\n')
    iniFile.write('DiscretizeIVH = 0\n')
    iniFile.write('DiscretizeIVHSeparated = 0\n')
    iniFile.write('UseFixedBinWidthIVH = 0\n')
    iniFile.write('BinWidthIVH = 0.25\n')
    iniFile.write('UseFixedNrBinsIVH = 1\n')
    iniFile.write('NrBinsIVH = 64\n\n')
    iniFile.write('[Interpolation]\n')
    iniFile.write('Rebinning_centering = 0\n')
    iniFile.write('2DInterpolation = 0\n')
    iniFile.write('InterpolationMethod = linear\n')
    iniFile.write('UseDownSampling2Cubic = 0\n')
    iniFile.write('UseUpSampling2Cubic = 0\n')
    iniFile.write('UseSamplingTo2mmVoxel = 1\n\n')
    iniFile.write('[ReSegmentation]\n')
    iniFile.write('ReSegmentImage = 0\n')
    iniFile.write('ExcludeOutliers = 0\n')
    iniFile.write('MinValueInReSegmentation = 0\n')
    iniFile.write('MaxValueInReSegmentation = 0\n\n')
    iniFile.write('[ImageProperties]\n')
    iniFile.write('ImageType = PET\n\n')
    iniFile.write('[NGLDMParameters]\n')
    iniFile.write('dist = 1\n')
    iniFile.write('coarseness = 0\n\n')
    iniFile.write('[NGTDMDistance]\n')
    iniFile.write('dist = 1\n\n')
    iniFile.write('[DistanceWeightProperties]\n')
    iniFile.write('NormGLCM = Chebyshev\n')
    iniFile.write('NormGLRLM = Chebyshev\n')
    iniFile.write('NormNGTDM = Chebyshev\n\n')
    iniFile.write('[ExtendedEmphasisFeatures]\n')
    iniFile.write('CalculateExtendedEmph = 0\n')
    iniFile.write('PowerRow = 1\n')
    iniFile.write('PowerCol = 1\n\n')
    iniFile.write('[OutputInformation]\n')
    iniFile.write('GetOneCSVFile = 1\n')
    iniFile.write('csvOutput = 1\n')
    iniFile.write('OntologyOutput = 0\n\n')
    iniFile.close()
    
createConfigFile('C:/RadiomicsTool', 'config128', 128)
createConfigFile('C:/RadiomicsTool', 'config64', 64)
createConfigFile('C:/RadiomicsTool', 'config32', 32)


