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
  download_libigl()
  add_subdirectory(${${PROJECT_NAME}_EXTERNAL}/libigl EXCLUDE_FROM_ALL)
  # Set Eigen directory to the one in libigl (needed for EVCTCD)
  set(ENV{EIGEN3_INCLUDE_DIR} "${${PROJECT_NAME}_EXTERNAL}/libigl/external/eigen/")
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
    download_spdlog()
    add_library(spdlog INTERFACE)
    add_library(spdlog::spdlog ALIAS spdlog)
    target_include_directories(spdlog SYSTEM INTERFACE ${${PROJECT_NAME}_EXTERNAL}/spdlog/include)
endif()

if(${PROJECT_NAME}_WITH_COMPARISONS)
  # Etienne Vouga's CTCD Library
  download_evctcd()
  add_subdirectory(${${PROJECT_NAME}_EXTERNAL}/EVCTCD)
  # These includes are PRIVATE for some reason
  target_include_directories(collisiondetection PUBLIC "${${PROJECT_NAME}_EXTERNAL}/EVCTCD/include")
  # Turn of floating point contraction for CCD robustness
  target_compile_options(collisiondetection PUBLIC "-ffp-contract=off")
  # Rename for convenience
  add_library(EVCTCD ALIAS collisiondetection)

  # Brochu et al. [2012] and Tang et al. [2014]
  download_exact_ccd()
  add_subdirectory(${${PROJECT_NAME}_EXTERNAL}/exact-ccd EXCLUDE_FROM_ALL)
  add_library(exact-ccd::exact-ccd ALIAS exact-ccd)

  # Rational implmentation of Brochu et al. [2012]
  download_rational_ccd()
  add_subdirectory(${${PROJECT_NAME}_EXTERNAL}/rational_ccd)
endif()
