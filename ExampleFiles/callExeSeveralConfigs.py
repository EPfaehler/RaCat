# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 12:36:26 2018

@author: PfaehlerEAG
"""

import shlex, subprocess
import os
#in this file the executable is called several times with different config files
#in this example, the mask image is named 'VOI...' and the names of the config files start with 'config'
#the user calls the function giving the folder name
def buildExecuteCommandSeveralConfig(folderName):
    configNames = []
    outputFolderNames = []
    #go to folder
    os.chdir(folderName)
    #go over all diles in folder
    for file in os.listdir('.'):
        #if file is mask
        if file[0:3]=='VOI':
            voiPath = folderName + '/' + file
        #if file is image    
        if file[0:3]!='VOI' and file[0:6]!='config':
            
            imagePath = folderName + '/' + file
        #if file is a configuration file, file array with config names
        if file[0:6] == 'config':
            configPath = folderName + '/' + file
            configNames.append(configPath)
            
            outputFolder = folderName + '/output' + file[6:]
            outputFolderNames.append(outputFolder)
    #for every configuration, call executable
    #it is necessary to change the path of the executable in this part!!
    for nrConfig in range(len(configNames)):
        configPath = configNames[nrConfig]
        patientInfoPath = folderName + '/' + 'patientInfo.ini'
        outputFolderName = outputFolderNames[nrConfig]
        command = "C:/RadiomicsTool/radiomics/Release/RaCaT_v1.3.exe" +" --ini "+ configPath + " --out " + outputFolderName+ " --img " + imagePath + " --voi " + voiPath + " --pat " + patientInfoPath
        args = shlex.split(command)
            
        subprocess.call(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    os.chdir(folderName)




def buildExecuteCommand(folderName):
    os.chdir(folderName)
    for subfolder in os.listdir('.'):
        if os.path.isdir(subfolder):
            os.chdir(subfolder)
            for file in os.listdir('.'):
                if file[0:3]=='VOI':
                    voiPath = folderName + '/' + file
                    segPath = voiPath + '/' +seg
            
                if file[0:3]!='VOI':
            
                    imagePath = folderName + '/' + file
                    outputFolderName = folderName + '/' + features
                iniPath = folderName + '/' + 'config.ini'
                patientInfoPath = folderName + '/' + 'patientInfo.ini'
    
                command = "C:/RadiomicsTool/radiomics/Release/Radiomics_v1.1.exe" +" --ini "+ iniPath + " --out " + outputFolderName+ " --img " + imagePath + " --voi " + voiPath + " --pat " + patientInfoPath
                args = shlex.split(command)
                subprocess.call(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
            os.chdir(folderName)

    


