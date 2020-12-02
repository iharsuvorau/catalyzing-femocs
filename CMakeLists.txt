cmake_minimum_required(VERSION 3.16)

# General Settings

set(TARGET CatalystAdaptor)
set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_BUILD_TYPE Release)
set(FEMOCS_DIR ${CMAKE_SOURCE_DIR}/femocs)

project(CatalystAdaptor CXX)

# Dependencies

find_package(deal.II 9.2 REQUIRED HINTS ${FEMOCS_DIR}/dealii)
if (NOT deal.II_FOUND)
    message(FATAL_ERROR "Deal.II not found from " ${DEAL_II_DIR})
endif()
DEAL_II_INITIALIZE_CACHED_VARIABLES()

set(CMAKE_LIBRARY_PATH ${FEMOCS_DIR}/lib ${FEMOCS_DIR}/GETELEC/lib ${FEMOCS_DIR}/dealii/lib)
find_library(FEMOCS_LIB femocs)
find_library(TET_LIB tet)
find_library(GETELEC_LIB getelec)
find_library(SLATEC_LIB slatec)

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

# Building

add_executable(${TARGET} src/main.cpp)
target_include_directories(${TARGET} PRIVATE
  ${FEMOCS_DIR}/include
  ${FEMOCS_DIR}/lib
  ${FEMOCS_DIR}/dealii/include
  ${FEMOCS_DIR}/GELTELEC/modules
)
target_link_libraries(${TARGET}
  ${FEMOCS_LIB}
  ${TET_LIB}
  ${GETELEC_LIB}
  ${SLATEC_LIB}
)
DEAL_II_SETUP_TARGET(${TARGET})
