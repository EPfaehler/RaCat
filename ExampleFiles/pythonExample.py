# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 17:18:52 2018

@author: PfaehlerEAG
"""
import shlex, subprocess
import os


folderName = outputFolderNames[i] + ('tumor%d'%tumor)
command = "C:/RadiomicsTool/radiomics/Release/Radiomics_v1.1.exe" +" --ini "+ iniPath + " --out " + folderName+ " --img " + imageFiles[i] + " --voi " + voiNames[tumor] + " --pat " + patientInfoPath
args = shlex.split(command)
subprocess.call(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
