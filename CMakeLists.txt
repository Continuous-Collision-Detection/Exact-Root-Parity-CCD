cmake_minimum_required(VERSION 3.8)
project(CCD)

################################################################################

# project-options
option(BUILD_UNIT_TESTS  "Build unit-tests"  ON)

################################################################################

### Configuration
set(${PROJECT_NAME}_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Dependencies)

################################################################################
# Library
################################################################################

# add_library(${PROJECT_NAME}
#   # Add source files here
# )
#
# target_include_directories(${PROJECT_NAME} PUBLIC src)

################################################################################
# Executable
################################################################################

add_executable(${PROJECT_NAME}_bin
  # Add executable source here
  src/main.cpp
)

target_include_directories(${PROJECT_NAME}_bin PRIVATE src)
# target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME})

################################################################################
# Required Libraries
################################################################################

################################################################################
# Tests
################################################################################

if(BUILD_UNIT_TESTS)
    # Catch2
    if(NOT TARGET Catch2::Catch2)
        download_catch2()
        add_subdirectory(${${PROJECT_NAME}_EXTERNAL}/Catch2 catch2)
        list(APPEND CMAKE_MODULE_PATH ${${PROJECT_NAME}_EXTERNAL}/Catch2/contrib)
    endif()

    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()
