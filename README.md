
   

A detailed documentation about how to use RaCaT and how to run it from source can be found in the documentation folder. However, here a short overview:

In order to make the RaCaT run from source, follow the following steps:

1. Clone Radiomics repository

Download the files of this repository.

2. Get CMake
   CMake can be downloaded here: https://cmake.org/install/

3. Get ITK
     Windows
    Here we only describe the installation using Microsoft Visual Studio. For all other compilers, check: https://itk.org/Wiki/ITK/Getting_Started/Build/Windows
    - Launch CMake GUI 
    - Go to field "Where is the source code:", click "Browse Source..." and navigate to where you cloned the repository with Git.
    - Go to "Where to build the binaries:", select "Browse Build..." and select a place to build the ITK library. It should NOT be the same directory as the one where you cloned the repository.
    - Click "Configure", and then specify "Visual Studio x" as the generator for this project.
    - Choose your build options
    - Click "Generate".
    - Open Visual Studio x and open the ALL_BUILD project that is in the folder where you built the ITK project
    - Click F5 (run solution)
    3.2 Linux
    - Download the ITK source code: https://itk.org/ITK/resources/software.html
    - unpack tarball: sudo tar xvzf InsightToolkit-3.14.0.tar.gz
    - Create a directory, where ITK should be built, e.g: sudo mkdir /usr/local/itk/InsightToolkit-3.14.0/ITK-build
    - Go to this directory and run ccmake: sudo ccmake -DITK_USE_REVIEW=ON 
    - Press "c" to configure and "g" to generate
    - Then you can run ITK: sudo make
    - Then you can install ITK: sudo make install

4. Get boost version 1_70_0 (thats important as the newer version have some problems with ITK!)
   The Radiomics Toolbox needs also the boost library. Download the source code from here: http://www.boost.org/users/history/version_1_61_0.html and unpack the folder.
   We have to install some boost libraries:
   - Go to the folder, where you unpacked the boost libraries
   - Go to the directory tools\build\.
   - Run bootstrap.bat
   - Run b2 install --prefix=PREFIX where PREFIX is the directory where you want Boost.Build to be installed
   - Add PREFIX\bin to your PATH environment variable.

   Boost is installed.

5. Use CMakeList to configure, generate and compile Radiomics.exe
   In order to run the radiomics code with the libraries, the following CMakeLists.txt file is required:
   The /PATH/TO elements have to be replaced by the corresponding paths on your computer.
   
    cmake_minimum_required(VERSION 2.8)
    
    project(Radiomics)
    
    SET (BOOST_ROOT /PATH/TO/BOOST)
    
    SET (BOOST_LIBRARYDIR "/PATH/TO/BOOST/stage/lib")
    
    SET (BOOST_MIN_VERSION "1.55.0")
    
    set (Boost_NO_BOOST_CMAKE ON)
    
    FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} REQUIRED)
    
    if (NOT Boost_FOUND)
    
      message(FATAL_ERROR "Fatal error: Boost (version >= 1.55) required.")
    
    else()
     
     message(STATUS "Setting up BOOST")
     
     message(STATUS " Library  - ${Boost_LIBRARY_DIRS}")
     
     include_directories(${Boost_INCLUDE_DIRS})
     
     link_directories(${Boost_LIBRARY_DIRS})
    
    endif (NOT Boost_FOUND)
    
    find_package(ITK REQUIRED)
    
    include(${ITK_USE_FILE})
    
    add_executable(Radiomics MACOSX_BUNDLE /PATH/TO/radiomics-master/main.cpp)
    
    target_link_libraries(Radiomics
    
    ${Boost_LIBRARIES} ${Glue}  ${VTK_LIBRARIES} ${ITK_LIBRARIES})


   Open the CMAKE GUI.

6. Run Radiomics.exe
