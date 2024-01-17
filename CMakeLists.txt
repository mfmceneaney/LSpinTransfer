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
#set(ARGPARSE_TEST_ENABLE "OFF") #NOTE: THIS IS NEEDED AT LEAST ON MACOS
#add_subdirectory("argparse")

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
add_executable(massfit massfit.cpp)
target_link_libraries(massfit PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitmc massfitmc.cpp)
target_link_libraries(massfitmc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitgauss massfitgauss.cpp)
target_link_libraries(massfitgauss PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitmcgauss massfitmcgauss.cpp)
target_link_libraries(massfitmcgauss PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitmcfixparams massfitmcfixparams.cpp)
target_link_libraries(massfitmcfixparams PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysis analysis.cpp)
target_link_libraries(analysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismc analysismc.cpp)
target_link_libraries(analysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismcfullinjection analysismcfullinjection.cpp)
target_link_libraries(analysismcfullinjection PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysisv2 analysisv2.cpp)
target_link_libraries(analysisv2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismcv2 analysismcv2.cpp)
target_link_libraries(analysismcv2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismcfullinjectionv2 analysismcfullinjectionv2.cpp)
target_link_libraries(analysismcfullinjectionv2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysisgausscbdiff analysisgausscbdiff.cpp)
target_link_libraries(analysisgausscbdiff PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(transverseanalysis transverseanalysis.cpp)
target_link_libraries(transverseanalysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
