cmake_minimum_required(VERSION 3.8)
project(CCD)

################################################################################

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(CCD_TOPLEVEL_PROJECT ON)
else()
    set(CCD_TOPLEVEL_PROJECT OFF)
endif()

# project-options
option(CCD_WITH_UNIT_TESTS   "Build unit tests using Catch2"  ON)
option(CCD_BUILD_RATIONAL    "Buil Rational Version"          ON)


# libigl
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK              "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM               "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            OFF)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)
################################################################################




### Configuration
set(CCD_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Dependencies)
include(CCDUtils)

################################################################################
# Rational Library
################################################################################
if(CCD_BUILD_RATIONAL)
  add_subdirectory(src)
  prepend_current_path(SOURCES)
  CCD_copy_headers(${SOURCES})
  CCD_set_source_group(${SOURCES})

  add_library(${PROJECT_NAME}_rational ${SOURCES})
  target_include_directories(${PROJECT_NAME}_rational PUBLIC ${PROJECT_BINARY_DIR}/include)

  # libigl
  target_link_libraries(${PROJECT_NAME}_rational PUBLIC igl::core)

  # Logger
  target_link_libraries(${PROJECT_NAME}_rational PUBLIC spdlog::spdlog)

  # Extra warnings
  target_link_libraries(${PROJECT_NAME}_rational PRIVATE warnings::all)

  #GMP
  target_include_directories(${PROJECT_NAME}_rational PUBLIC ${GMP_INCLUDE_DIRS} src)
  target_link_libraries(${PROJECT_NAME}_rational PUBLIC ${GMP_LIBRARIES})
endif()

################################################################################
# Floating-point Library
################################################################################
# doubleCCD
add_subdirectory(doubleccd)
prepend_current_path(DOUBLECCD_SOURCES)
doubleCCD_copy_headers(${DOUBLECCD_SOURCES})
CCD_set_source_group(${DOUBLECCD_SOURCES})

add_library(${PROJECT_NAME}_double  ${DOUBLECCD_SOURCES})
target_include_directories(${PROJECT_NAME}_double PUBLIC ${PROJECT_BINARY_DIR}/include)

# libigl
target_link_libraries(${PROJECT_NAME}_double PUBLIC igl::core igl::predicates)

# Logger
target_link_libraries(${PROJECT_NAME}_double PUBLIC spdlog::spdlog)

# Extra warnings
target_link_libraries(${PROJECT_NAME}_double PRIVATE warnings::all)

#GMP
find_package(GMPECCD)
IF(NOT ${GMP_FOUND})
        MESSAGE(FATAL_ERROR "Cannot find GMP")
ENDIF()
target_include_directories(${PROJECT_NAME}_double PUBLIC doubleccd)
target_include_directories(${PROJECT_NAME}_double PUBLIC ${GMP_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_double PUBLIC ${GMP_LIBRARIES})


################################################################################
# Executable
################################################################################
if(CCD_TOPLEVEL_PROJECT AND CCD_BUILD_RATIONAL)
  add_executable(${PROJECT_NAME}_bin
    # Add executable source here
    src/main.cpp
    src/read_collision_data.cpp
  )

  target_link_libraries(${PROJECT_NAME}_bin PUBLIC HighFive::HighFive)

  # Boost
  set(Boost_USE_MULTITHREADED TRUE)
  find_package(Boost COMPONENTS filesystem REQUIRED)
  target_link_libraries(${PROJECT_NAME}_bin PUBLIC Boost::filesystem)

  target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_rational)
  target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_double)
  set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/")
  target_compile_definitions(${PROJECT_NAME}_bin PUBLIC -DCCD_DATA_PATH=\"${DATA_PATH}\")
endif()

################################################################################
# Tests
################################################################################

if(CCD_WITH_UNIT_TESTS)
  # Catch2
  if(NOT TARGET Catch2::Catch2)
      ccd_download_catch2()
      add_subdirectory(${CCD_EXTERNAL}/Catch2 catch2)
      list(APPEND CMAKE_MODULE_PATH ${CCD_EXTERNAL}/Catch2/contrib)
  endif()

  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()
