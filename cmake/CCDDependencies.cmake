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
  
  add_library(indirect_predicates ${CCD_EXTERNAL}/indirect_predicates/predicates/indirect_predicates.cpp )
  target_include_directories(indirect_predicates PUBLIC ${CCD_EXTERNAL}/indirect_predicates)

endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
    ccd_download_spdlog()
    add_library(spdlog INTERFACE)
    add_library(spdlog::spdlog ALIAS spdlog)
    target_include_directories(spdlog SYSTEM INTERFACE ${CCD_EXTERNAL}/spdlog/include)
endif()

if(CCD_WITH_COMPARISONS AND NOT TARGET CCDWrapper)
  ccd_download_comparisons()
  # Force CCD Wrapper to not download Catch2
  set(CCD_WRAPPER_WITH_UNIT_TESTS OFF CACHE BOOL "" FORCE)
  add_subdirectory(${CCD_EXTERNAL}/comparisons EXCLUDE_FROM_ALL)
endif()
