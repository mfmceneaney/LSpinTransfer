cmake_minimum_required(VERSION 3.4)
project(LSpinTransfer)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#----------------------------------------------------------------------#
# CERN ROOT DEPENDENCY
find_package(ROOT CONFIG COMPONENTS RIO Net RooStats RooFit ROOTDataFrame)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------#
# YAML DEPENDENCY
add_subdirectory("yaml-cpp")

#----------------------------------------------------------------------#
# ARGPARSE DEPENDENCY
set(ARGPARSE_TEST_ENABLE "OFF") #NOTE: THIS IS NEEDED AT LEAST ON MACOS
add_subdirectory("argparse")

#----------------------------------------------------------------------#
# PROJECT DEPENDENCY
#message(${CMAKE_CURRENT_SOURCE_DIR}/)
#file(GLOB PROJECT_SRC_FILES ${PROJECT_SOURCE_DIR}/*.hpp)
#ADD_LIBRARY(AnalysisModule ${PROJECT_SRC_FILES})
include_directories(${PROJECT_SOURCE_DIR}/include)

#message(${CMAKE_CURRENT_SOURCE_DIR}/)
#file(GLOB PROJECT_SRC_FILES ${PROJECT_SOURCE_DIR}/*.hpp)
#ADD_LIBRARY(AnalysisModule ${PROJECT_SRC_FILES})

#----------------------------------------------------------------------#
# ADD EXECUTABLES
#add_executable(massfit massfit.cpp)
#target_link_libraries(massfit PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp argparse)
#add_executable(massfitmc massfitmc.cpp)
#target_link_libraries(massfitmc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp argparse)
add_executable(analysis analysis.cpp)
target_link_libraries(analysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp argparse)
