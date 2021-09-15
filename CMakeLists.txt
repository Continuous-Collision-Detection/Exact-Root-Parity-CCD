# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(CCD_TOPLEVEL_PROJECT OFF)
else()
    set(CCD_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.14.0")
if(CCD_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build CCD is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()

# Include user-provided default options if available. We do that before the main
# `project()` so that we can define the C/C++ compilers from the option file.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/CCDOptions.cmake)
    message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/CCDOptions.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/CCDOptions.cmake)
endif()

################################################################################

project(CCD
        DESCRIPTION "A CCD library."
        LANGUAGES CXX)

# Project options
option(CCD_WITH_UNIT_TESTS "Build unit tests using Catch2"  ${CCD_TOPLEVEL_PROJECT})
option(CCD_BUILD_RATIONAL  "Build Rational Version"                             OFF)
option(CCD_ROUND_INPUTS    "Round the inputs before collision detection"        OFF)

# Set default minimum C++ standard
if(CCD_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# Configuration
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/ccd/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find/")

# Color output
include(ccd_use_colors)

# CCD utils
include(ccd_utils)

################################################################################
# Rational Library
################################################################################

if(CCD_BUILD_RATIONAL)
  add_subdirectory(src)
  prepend_current_path(SOURCES)
  ccd_copy_headers(${SOURCES})
  ccd_set_source_group(${SOURCES})

  add_library(ccd_rational ${SOURCES})
  add_library(ccd::rational ALIAS ccd_rational)

  target_include_directories(ccd_rational PUBLIC "${PROJECT_BINARY_DIR}/include")

  # ------------
  # Dependencies
  # ------------

  # Extra warnings
  include(ccd_warnings)
  target_link_libraries(ccd_rational PRIVATE ccd::warnings)

  # libigl
  include(eigen)
  include(libigl)
  target_link_libraries(ccd_rational PUBLIC
    Eigen3::Eigen
    igl::core
  )

  # GMP
  find_package(GMP)
  target_link_libraries(ccd_rational PUBLIC gmp::gmp)
endif()

################################################################################
# Floating-point Library
################################################################################
add_subdirectory(doubleccd)
prepend_current_path(DOUBLECCD_SOURCES)
doubleccd_copy_headers(${DOUBLECCD_SOURCES})
ccd_set_source_group(${DOUBLECCD_SOURCES})

add_library(ccd_double ${DOUBLECCD_SOURCES})
add_library(ccd::double ALIAS ccd_double)

target_include_directories(ccd_double PUBLIC "${PROJECT_BINARY_DIR}/include")
if(CCD_ROUND_INPUTS)
  target_compile_definitions(ccd_double PUBLIC "CCD_ROUND_INPUTS")
endif()

# Extra warnings
include(ccd_warnings)
target_link_libraries(ccd_double PRIVATE ccd::warnings)

# libigl
include(eigen)
include(libigl)
target_link_libraries(ccd_double PUBLIC
  Eigen3::Eigen
  igl::core
  igl::predicates
)

# GMP (TODO: remove this dependency for ccd::double)
find_package(GMP)
target_link_libraries(ccd_double PUBLIC gmp::gmp)

# For perturbation vertices
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # grant IEEE 754 compliance
  target_compile_options(ccd_double PUBLIC  "-frounding-math")
  message(STATUS "In Fixed Root Parity CCD, the -frounding-math flag is setted")
  # use intrinsic functions (CHECK WHAT TO DO FOR GCC !!!!!!!!)
  #	target_compile_options(ccd_double PUBLIC "/Oi")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # grant IEEE 754 compliance
  target_compile_options(ccd_double PUBLIC "/fp:strict")
endif()

################################################################################
# Executable
################################################################################

if(CCD_TOPLEVEL_PROJECT AND CCD_BUILD_RATIONAL)
  add_executable(ccd_bin
    # Add executable source here
    app/main.cpp
    # src/read_collision_data.cpp
  )

  # target_link_libraries(ccd_bin PUBLIC HighFive::HighFive)

  target_link_libraries(ccd_bin PUBLIC ccd_rational)
  target_link_libraries(ccd_bin PUBLIC ccd_double)
  # target_link_libraries(ccd_bin PUBLIC ccd_interval)
  # WARNING: This requires a check  as it will break compatability with non-AVX
  #          processors.
  # target_compile_options(CCD_interval PUBLIC "-mavx2")
  set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/")
  target_compile_definitions(ccd_bin PUBLIC CCD_DATA_PATH=\"${DATA_PATH}\")
endif()

################################################################################
# Tests
################################################################################

if(CCD_WITH_UNIT_TESTS)
    include(CTest)
    enable_testing()

    # Include Catch2 and provide function `catch_discover_tests` to register tests.
    include(catch2)
    FetchContent_GetProperties(catch2)
    include("${catch2_SOURCE_DIR}/contrib/Catch.cmake")

    add_subdirectory(tests)
endif()
