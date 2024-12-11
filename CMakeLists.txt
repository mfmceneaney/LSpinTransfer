cmake_minimum_required(VERSION 3.4)
project(LSpinTransfer)

set (CMAKE_CXX_STANDARD 17)
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
#set executables directory
set(EXECUTABLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(test_splot ${EXECUTABLES_DIR}/test_splot.cpp)
target_link_libraries(test_splot PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(test_splot_lk ${EXECUTABLES_DIR}/test_splot_lk.cpp)
target_link_libraries(test_splot_lk PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(test_splot2d ${EXECUTABLES_DIR}/test_splot2d.cpp)
target_link_libraries(test_splot2d PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(test_splot_inject2d ${EXECUTABLES_DIR}/test_splot_inject2d.cpp)
target_link_libraries(test_splot_inject2d PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(massdistr ${EXECUTABLES_DIR}/massdistr.cpp)
target_link_libraries(massdistr PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massdistrmc ${EXECUTABLES_DIR}/massdistrmc.cpp)
target_link_libraries(massdistrmc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfit ${EXECUTABLES_DIR}/massfit.cpp)
target_link_libraries(massfit PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitmc ${EXECUTABLES_DIR}/massfitmc.cpp)
target_link_libraries(massfitmc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitgauss ${EXECUTABLES_DIR}/massfitgauss.cpp)
target_link_libraries(massfitgauss PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitmcgauss ${EXECUTABLES_DIR}/massfitmcgauss.cpp)
target_link_libraries(massfitmcgauss PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitmcfixparams ${EXECUTABLES_DIR}/massfitmcfixparams.cpp)
target_link_libraries(massfitmcfixparams PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitppimkp ${EXECUTABLES_DIR}/massfitppimkp.cpp)
target_link_libraries(massfitppimkp PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitppimkpmc ${EXECUTABLES_DIR}/massfitppimkpmc.cpp)
target_link_libraries(massfitppimkpmc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysisbincounts ${EXECUTABLES_DIR}/analysisbincounts.cpp)
target_link_libraries(analysisbincounts PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysis ${EXECUTABLES_DIR}/analysis.cpp)
target_link_libraries(analysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismc ${EXECUTABLES_DIR}/analysismc.cpp)
target_link_libraries(analysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismccosphidependence ${EXECUTABLES_DIR}/analysismccosphidependence.cpp)
target_link_libraries(analysismccosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismcfullinjection ${EXECUTABLES_DIR}/analysismcfullinjection.cpp)
target_link_libraries(analysismcfullinjection PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysisv2 ${EXECUTABLES_DIR}/analysisv2.cpp)
target_link_libraries(analysisv2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismcv2 ${EXECUTABLES_DIR}/analysismcv2.cpp)
target_link_libraries(analysismcv2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismcfullinjectionv2 ${EXECUTABLES_DIR}/analysismcfullinjectionv2.cpp)
target_link_libraries(analysismcfullinjectionv2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysisgausscbdiff ${EXECUTABLES_DIR}/analysisgausscbdiff.cpp)
target_link_libraries(analysisgausscbdiff PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysisgenericdiff ${EXECUTABLES_DIR}/analysisgenericdiff.cpp)
target_link_libraries(analysisgenericdiff PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)


add_executable(massfitlandau ${EXECUTABLES_DIR}/massfitlandau.cpp)
target_link_libraries(massfitlandau PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitmclandau ${EXECUTABLES_DIR}/massfitmclandau.cpp)
target_link_libraries(massfitmclandau PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysislandaucbdiff ${EXECUTABLES_DIR}/analysislandaucbdiff.cpp)
target_link_libraries(analysislandaucbdiff PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(transverseanalysis ${EXECUTABLES_DIR}/transverseanalysis.cpp)
target_link_libraries(transverseanalysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(bsaanalysis ${EXECUTABLES_DIR}/bsaanalysis.cpp)
target_link_libraries(bsaanalysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(bsaanalysismc ${EXECUTABLES_DIR}/bsaanalysismc.cpp)
target_link_libraries(bsaanalysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(alutbsaanalysis ${EXECUTABLES_DIR}/alutbsaanalysis.cpp)
target_link_libraries(alutbsaanalysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(alutbsaanalysismc ${EXECUTABLES_DIR}/alutbsaanalysismc.cpp)
target_link_libraries(alutbsaanalysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(lkbsaanalysis ${EXECUTABLES_DIR}/lkbsaanalysis.cpp)
target_link_libraries(lkbsaanalysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(lkbsaanalysismc ${EXECUTABLES_DIR}/lkbsaanalysismc.cpp)
target_link_libraries(lkbsaanalysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(dihadronbsaanalysismc ${EXECUTABLES_DIR}/dihadronbsaanalysismc.cpp)
target_link_libraries(dihadronbsaanalysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(dihadronbsaanalysismcbincounts ${EXECUTABLES_DIR}/dihadronbsaanalysismcbincounts.cpp)
target_link_libraries(dihadronbsaanalysismcbincounts PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(dihadronbsaanalysismc1 ${EXECUTABLES_DIR}/dihadronbsaanalysismc1.cpp)
target_link_libraries(dihadronbsaanalysismc1 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(dihadronbsaanalysismc1bincounts ${EXECUTABLES_DIR}/dihadronbsaanalysismc1bincounts.cpp)
target_link_libraries(dihadronbsaanalysismc1bincounts PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysis2 ${EXECUTABLES_DIR}/alutbsaanalysis2.cpp)
target_link_libraries(alutbsaanalysis2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(alutbsaanalysismc2 ${EXECUTABLES_DIR}/alutbsaanalysismc2.cpp)
target_link_libraries(alutbsaanalysismc2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysismc2cosphidependence ${EXECUTABLES_DIR}/alutbsaanalysismc2cosphidependence.cpp)
target_link_libraries(alutbsaanalysismc2cosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysismc2sinphidependence ${EXECUTABLES_DIR}/alutbsaanalysismc2sinphidependence.cpp)
target_link_libraries(alutbsaanalysismc2sinphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysis2Dmc2cosphidependence ${EXECUTABLES_DIR}/alutbsaanalysis2Dmc2cosphidependence.cpp)
target_link_libraries(alutbsaanalysis2Dmc2cosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysis2D2cosphidependence ${EXECUTABLES_DIR}/alutbsaanalysis2D2cosphidependence.cpp)
target_link_libraries(alutbsaanalysis2D2cosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysis3Dmc2cosphidependence ${EXECUTABLES_DIR}/alutbsaanalysis3Dmc2cosphidependence.cpp)
target_link_libraries(alutbsaanalysis3Dmc2cosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(bsaanalysisMC1DRooFitML ${EXECUTABLES_DIR}/bsaanalysisMC1DRooFitML.cpp)
target_link_libraries(bsaanalysisMC1DRooFitML PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(bsaanalysisMC2DRooFitML ${EXECUTABLES_DIR}/bsaanalysisMC2DRooFitML.cpp)
target_link_libraries(bsaanalysisMC2DRooFitML PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
