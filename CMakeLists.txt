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
add_executable(test_splot test_splot.cpp)
target_link_libraries(test_splot PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(test_splot_lk test_splot_lk.cpp)
target_link_libraries(test_splot_lk PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(test_splot2d test_splot2d.cpp)
target_link_libraries(test_splot2d PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(test_splot_inject2d test_splot_inject2d.cpp)
target_link_libraries(test_splot_inject2d PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(massdistr massdistr.cpp)
target_link_libraries(massdistr PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massdistrmc massdistrmc.cpp)
target_link_libraries(massdistrmc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
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
add_executable(massfitppimkp massfitppimkp.cpp)
target_link_libraries(massfitppimkp PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitppimkpmc massfitppimkpmc.cpp)
target_link_libraries(massfitppimkpmc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysisbincounts analysisbincounts.cpp)
target_link_libraries(analysisbincounts PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysis analysis.cpp)
target_link_libraries(analysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismc analysismc.cpp)
target_link_libraries(analysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysismccosphidependence analysismccosphidependence.cpp)
target_link_libraries(analysismccosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
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
add_executable(analysisgenericdiff analysisgenericdiff.cpp)
target_link_libraries(analysisgenericdiff PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)


add_executable(massfitlandau massfitlandau.cpp)
target_link_libraries(massfitlandau PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(massfitmclandau massfitmclandau.cpp)
target_link_libraries(massfitmclandau PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(analysislandaucbdiff analysislandaucbdiff.cpp)
target_link_libraries(analysislandaucbdiff PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(transverseanalysis transverseanalysis.cpp)
target_link_libraries(transverseanalysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(bsaanalysis bsaanalysis.cpp)
target_link_libraries(bsaanalysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(bsaanalysismc bsaanalysismc.cpp)
target_link_libraries(bsaanalysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(alutbsaanalysis alutbsaanalysis.cpp)
target_link_libraries(alutbsaanalysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(alutbsaanalysismc alutbsaanalysismc.cpp)
target_link_libraries(alutbsaanalysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(lkbsaanalysis lkbsaanalysis.cpp)
target_link_libraries(lkbsaanalysis PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(lkbsaanalysismc lkbsaanalysismc.cpp)
target_link_libraries(lkbsaanalysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(dihadronbsaanalysismc dihadronbsaanalysismc.cpp)
target_link_libraries(dihadronbsaanalysismc PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(dihadronbsaanalysismcbincounts dihadronbsaanalysismcbincounts.cpp)
target_link_libraries(dihadronbsaanalysismcbincounts PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(dihadronbsaanalysismc1 dihadronbsaanalysismc1.cpp)
target_link_libraries(dihadronbsaanalysismc1 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(dihadronbsaanalysismc1bincounts dihadronbsaanalysismc1bincounts.cpp)
target_link_libraries(dihadronbsaanalysismc1bincounts PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysis2 alutbsaanalysis2.cpp)
target_link_libraries(alutbsaanalysis2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
add_executable(alutbsaanalysismc2 alutbsaanalysismc2.cpp)
target_link_libraries(alutbsaanalysismc2 PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysismc2cosphidependence alutbsaanalysismc2cosphidependence.cpp)
target_link_libraries(alutbsaanalysismc2cosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysismc2sinphidependence alutbsaanalysismc2sinphidependence.cpp)
target_link_libraries(alutbsaanalysismc2sinphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysis2Dmc2cosphidependence alutbsaanalysis2Dmc2cosphidependence.cpp)
target_link_libraries(alutbsaanalysis2Dmc2cosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysis2D2cosphidependence alutbsaanalysis2D2cosphidependence.cpp)
target_link_libraries(alutbsaanalysis2D2cosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(alutbsaanalysis3Dmc2cosphidependence alutbsaanalysis3Dmc2cosphidependence.cpp)
target_link_libraries(alutbsaanalysis3Dmc2cosphidependence PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(bsaanalysisMC1DRooFitML bsaanalysisMC1DRooFitML.cpp)
target_link_libraries(bsaanalysisMC1DRooFitML PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)

add_executable(bsaanalysisMC2DRooFitML bsaanalysisMC2DRooFitML.cpp)
target_link_libraries(bsaanalysisMC2DRooFitML PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
