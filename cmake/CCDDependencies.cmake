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

#Indirect_Predicates
if(NOT TARGET indirect_predicates)
  if(WIN32)
  #target_compile_definitions(${PROJECT_NAME}_bin PUBLIC NOMINMAX)
  else()
    ccd_download_indirect_predicates()
  endif()
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
