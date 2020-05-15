# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(CCDDownloadExternal)

################################################################################
# Required libraries
################################################################################

# libigl
if(NOT TARGET igl)
  ccd_download_libigl()
  add_subdirectory(${CCD_EXTERNAL}/libigl EXCLUDE_FROM_ALL)
  # Set Eigen directory to the one in libigl (needed for EVCTCD)
  set(ENV{EIGEN3_INCLUDE_DIR} "${CCD_EXTERNAL}/libigl/external/eigen/")
endif()


ccd_download_geogram()
include(geogram)

#Indirect_Predicates
if(NOT TARGET indirect_predicates)
  ccd_download_indirect_predicates()

  set(INDIRECTPREDICATES_SOURCES
    ${CCD_EXTERNAL}/indirect_predicates/implicit_point.cpp
    ${CCD_EXTERNAL}/indirect_predicates/implicit_point.h
    ${CCD_EXTERNAL}/indirect_predicates/numerics.cpp
    ${CCD_EXTERNAL}/indirect_predicates/numerics.h
    ${CCD_EXTERNAL}/indirect_predicates/predicates/indirect_predicates.cpp
    ${CCD_EXTERNAL}/indirect_predicates/predicates/indirect_predicates.h
  )
  add_library(indirect_predicates  ${INDIRECTPREDICATES_SOURCES})
  target_include_directories(indirect_predicates PUBLIC ${CCD_EXTERNAL}/indirect_predicates ${CCD_EXTERNAL}/indirect_predicates/predicates)

  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_compile_options(indirect_predicates PRIVATE "/fp:strict")
  elseif(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(indirect_predicates PRIVATE "-frounding-math")
  else()
    message(WARNING "On CLang, there is no compiler flag which is required for ensuring the correctness of the algorithm.")
  endif()
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
  ccd_download_spdlog()
  add_library(spdlog INTERFACE)
  add_library(spdlog::spdlog ALIAS spdlog)
  target_include_directories(spdlog SYSTEM INTERFACE ${CCD_EXTERNAL}/spdlog/include)
endif()

# HDF5 Reader
if(NOT TARGET HighFive::HighFive)
  set(USE_EIGEN TRUE CACHE BOOL "Enable Eigen testing" FORCE)
  ccd_download_high_five()
  add_subdirectory(${CCD_EXTERNAL}/HighFive EXCLUDE_FROM_ALL)
  add_library(HighFive::HighFive ALIAS HighFive)
endif()
