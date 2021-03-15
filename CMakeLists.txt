# Usage:
#   cmake -Ddeal.II_DIR=femocs/dealii/lib/cmake/deal.II -DParaView_DIR=paraview_build -DFEMOCS_DIR=femocs -DCMAKE_PREFIX_PATH="femocs/lib;femocs/GETELEC/lib" -S . -B build
#   cmake --build build

cmake_minimum_required(VERSION 3.16)

# General Settings

set(TARGET Demo)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

project(CatalyzingFemocs LANGUAGES CXX C)

# Dependencies

# femocs and its dependencies
find_package(deal.II 9.2 REQUIRED)
DEAL_II_INITIALIZE_CACHED_VARIABLES()
find_library(FEMOCS_LIB femocs)
find_library(TET_LIB tet)
find_library(GETELEC_LIB getelec)
find_library(SLATEC_LIB slatec)
find_library(PETSC_LIB petsc /usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/lib)

# ParaView Catalyst
find_package(MPI COMPONENTS C REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
find_package(ParaView 5.8.1 REQUIRED)

# Info messages
message("deal.II version: ${deal.II_VERSION}")
message("ParaView version: ${ParaView_VERSION}")
message("VTK version: ${VTK_VERSION}")

# CatalystAdaptor library
add_library(CatalystAdaptor src/CatalystAdaptor.cpp)
target_include_directories(CatalystAdaptor PRIVATE
        ${FEMOCS_DIR}/include
        ${FEMOCS_DIR}/lib
        ${FEMOCS_DIR}/dealii/include
        ${FEMOCS_DIR}/GELTELEC/modules
        /usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/include)
target_link_libraries(CatalystAdaptor
        INTERFACE
        VTK::PythonUsed
        PRIVATE
        ParaView::PythonCatalyst
        VTK::CommonDataModel
        ${MPI_C_LIBRARIES}
        ${FEMOCS_LIB}
        ${TET_LIB}
        ${GETELEC_LIB}
        ${SLATEC_LIB}
        ${PETSC_LIB})
add_definitions("-DUSE_CATALYST")

# Compilation flags

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
# remove all flags that control warnings
string(REGEX REPLACE " -W[A-Za-z0-9-]+" "" DEAL_II_CXX_FLAGS ${DEAL_II_CXX_FLAGS})
string(REPLACE "-pedantic" "" DEAL_II_CXX_FLAGS ${DEAL_II_CXX_FLAGS})
# remove optimization level flags
string(REGEX REPLACE "-O[0-9]" "" DEAL_II_CXX_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})
# remove warning about flag not present for f90
string(REPLACE "-Wno-unused-local-typedefs" "" DEAL_II_CXX_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})
set(DEAL_II_CXX_FLAGS "${DEAL_II_CXX_FLAGS} -O3 -w")

# CatalystAdaptor Demo

add_executable(${TARGET} src/main.cpp)
target_include_directories(${TARGET} PRIVATE
        ${FEMOCS_DIR}/include
        ${FEMOCS_DIR}/lib
        ${FEMOCS_DIR}/dealii/include
        ${FEMOCS_DIR}/GELTELEC/modules)
target_link_libraries(${TARGET}
        ${FEMOCS_LIB}
        ${TET_LIB}
        ${GETELEC_LIB}
        ${SLATEC_LIB}
        ${VTK_LIBRARIES}
        CatalystAdaptor)
DEAL_II_SETUP_TARGET(${TARGET})
vtk_module_autoinit(TARGETS ${TARGET} MODULES ${VTK_LIBRARIES})
