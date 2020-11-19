cmake_minimum_required (VERSION 2.8.8)

SET(TARGET "CatalystAdaptor")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# SET(CMAKE_Fortran_COMPILER gfortran)
# SET(CMAKE_C_COMPILER gcc)
# SET(CMAKE_CXX_COMPILER c++)
# SET(CMAKE_VERBOSE_MAKEFILE off)
# SET(CMAKE_BUILD_TYPE Release)

PROJECT(${TARGET})

# set explicitly number of threads used by Deal.II
# if not set, determination will be done automatically
# IF (N_DEALII GREATER 0)
#     MESSAGE(STATUS "Setting \# Deal.II threads to " ${N_DEALII})
#     SET(DEAL_II_NUM_THREADS N_DEALII)
# ENDIF()

# ==============================================
# ========== Path for input & output ===========
# ==============================================

# Specify folders where femocs, libfemocs and its .mod-file will be written
# Those paths are (not like others) specified in relation to the location of Makefile
#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${FEMOCS_REPO_PATH})

SET(FEMOCS_REPO_PATH ${CMAKE_SOURCE_DIR}/femocs)

SET(DEAL_II_DIR ${FEMOCS_REPO_PATH}/dealii)
SET(CGAL_DIR ${FEMOCS_REPO_PATH}/cgal)
SET(TET_DIR ${FEMOCS_REPO_PATH}/lib)
SET(GETELEC_DIR ${FEMOCS_REPO_PATH}/GETELEC/lib)
SET(SLATEC_DIR ${FEMOCS_REPO_PATH}/GETELEC/lib)
SET(FEMOCS_DIR ${FEMOCS_REPO_PATH}/lib)

SET(TET_PATH ${TET_DIR}/libtet.a)
SET(GETELEC_PATH ${GETELEC_DIR}/libgetelec.a)
SET(SLATEC_PATH ${SLATEC_DIR}/libslatec.a)
SET(FEMOCS_PATH ${FEMOCS_DIR}/libfemocs.a)
#SET(FEMOCS_FORTRAN_PATH ${FEMOCS_DIR}/libfemocs.mod)
# TODO: how to include libfemocs.mod?

# path to the directories with headers
IF (LIBCGAL)
    INCLUDE_DIRECTORIES(${FEMOCS_REPO_PATH}/include ${FEMOCS_REPO_PATH}/GETELEC/modules ${FEMOCS_REPO_PATH}/cgal/include)
ELSE()
    INCLUDE_DIRECTORIES(${FEMOCS_REPO_PATH}/include ${FEMOCS_REPO_PATH}/GETELEC/modules)
ENDIF()

# ==============================================
# =========== Locate dependencies ==============
# ==============================================

FIND_PACKAGE(deal.II 8.4 REQUIRED HINTS ${DEAL_II_DIR})
IF ( NOT deal.II_FOUND )
    MESSAGE(FATAL_ERROR "Deal.II not found from " ${DEAL_II_DIR})
ENDIF()

IF (LIBCGAL)
    FIND_PACKAGE(CGAL HINTS ${CGAL_DIR})
    IF ( CGAL_FOUND )
        MESSAGE(STATUS "Using CGAL from " ${CGAL_DIR})
        INCLUDE( ${CGAL_USE_FILE} )
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_CGAL=true")
    ELSE()
        MESSAGE(FATAL_ERROR "CGAL not found from " ${CGAL_DIR}) 
    ENDIF()
ELSE()
    MESSAGE(STATUS "Disabling usage of CGAL")
ENDIF()


# ==============================================
# ============= Compilation flags ==============
# ==============================================

# # in debug mode enable asserts
# MESSAGE(STATUS "Disabling FEMOCS asserts")

# # add OpenMP flag
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# # Specify warnings that are checked and that are omitted
# #SET(WARNINGS "-pedantic -Wall -Wextra -Wno-unused-local-typedefs -Wno-unused-parameter")
# SET(WARNINGS "")

# # flags should not be set in CMAKE_CXX_FLAGS, because deal.II overrides many of them.
# # The deal.II flags can be seen from dealii/share/deal.II/Make.global_options

# # remove all flags that control warnings
# STRING(REGEX REPLACE " -W[A-Za-z0-9-]+" "" DEAL_II_CXX_FLAGS ${DEAL_II_CXX_FLAGS})
# STRING(REPLACE "-pedantic" "" DEAL_II_CXX_FLAGS ${DEAL_II_CXX_FLAGS})

# # remove optimization level flags
# STRING(REGEX REPLACE "-O[0-9]" "" DEAL_II_CXX_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})

# # remove warning about flag not present for f90
# STRING(REPLACE "-Wno-unused-local-typedefs" "" DEAL_II_CXX_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})

# # set warning and optimization level
# SET(DEAL_II_CXX_FLAGS "${DEAL_II_CXX_FLAGS} -O3 -w")

# # print out resulting flags
# MESSAGE(STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
# MESSAGE(STATUS "DEAL_II_CXX_FLAGS = " ${DEAL_II_CXX_FLAGS})
# MESSAGE(STATUS "DEAL_II_CXX_FLAGS_RELEASE = " ${DEAL_II_CXX_FLAGS_RELEASE})
# MESSAGE(STATUS "DEAL_II_LINKER_FLAGS = " ${DEAL_II_LINKER_FLAGS})
# MESSAGE(STATUS "DEAL_II_LINKER_FLAGS_RELEASE = " ${DEAL_II_LINKER_FLAGS_RELEASE})
# MESSAGE(STATUS "DEAL_II_USER_DEFINITIONS = " ${DEAL_II_USER_DEFINITIONS})
# MESSAGE(STATUS "DEAL_II_USER_DEFINITIONS_RELEASE = " ${DEAL_II_USER_DEFINITIONS_RELEASE})

# ==============================================
# =========== Start building ===========
# ==============================================

DEAL_II_INITIALIZE_CACHED_VARIABLES()

ADD_LIBRARY(tet ${TET_DIR})
SET_TARGET_PROPERTIES(tet PROPERTIES LINKER_LANGUAGE CXX)

ADD_LIBRARY(getelec ${GETELEC_DIR})
SET_TARGET_PROPERTIES(getelec PROPERTIES LINKER_LANGUAGE Fortran)

ADD_LIBRARY(slatec ${SLATEC_DIR})
SET_TARGET_PROPERTIES(slatec PROPERTIES LINKER_LANGUAGE C)

ADD_LIBRARY(femocs ${FEMOCS_DIR})
SET_TARGET_PROPERTIES(femocs PROPERTIES LINKER_LANGUAGE CXX)

# ADD_LIBRARY(femocs_fortran ${FEMOCS_DIR})
# SET_TARGET_PROPERTIES(femocs_fortran PROPERTIES LINKER_LANGUAGE Fortran)

message("Library directories:")
message("tet " ${TET_PATH})
message("getelec " ${GETELEC_PATH})
message("slatec " ${SLATEC_PATH})
message("femocs " ${FEMOCS_PATH})

# Specify the path to main file separately
ADD_EXECUTABLE(${TARGET} src/main.cpp)
DEAL_II_SETUP_TARGET(${TARGET})
TARGET_LINK_LIBRARIES(${TARGET}
  ${TET_PATH}
  ${GETELEC_PATH}
  ${SLATEC_PATH}
  ${FEMOCS_PATH}
)

