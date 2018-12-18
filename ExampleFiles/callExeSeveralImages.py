# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 12:36:26 2018

@author: PfaehlerEAG
"""

import shlex, subprocess
import os
#in this file the executable is called several times with different images
#in this example, every image and corresponding mask and corresponding configuration file are stored in one folder
# all image folders from which features should be extracted are stored in one folder
# the mask image is named 'VOI...' and the names of the config files start with 'config'
#the user calls the function giving the folder name



def buildExecuteCommand(folderName):
    #go to main folder
    os.chdir(folderName)
    #go to every subfolder
    for subfolder in os.listdir('.'):
        
        if os.path.isdir(subfolder):
            #go to subfolder
            os.chdir(subfolder)
            for file in os.listdir('.'):
                #if file is mask
                if file[0:3]=='VOI':
                    voiPath = folderName + '/' + file
                #if file is image
                if file[0:3]!='VOI':
            
                    imagePath = folderName + '/' + file
                    outputFolderName = folderName + '/' + features
                #set config path
                iniPath = folderName + '/' + 'config.ini'
                patientInfoPath = folderName + '/' + 'patientInfo.ini'
                #please change here the path of the executable
                command = "C:/RadiomicsTool/radiomics/Release/RaCaT_v1.3.exe" +" --ini "+ iniPath + " --out " + outputFolderName+ " --img " + imagePath + " --voi " + voiPath + " --pat " + patientInfoPath
                args = shlex.split(command)
                subprocess.call(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
            os.chdir(folderName)