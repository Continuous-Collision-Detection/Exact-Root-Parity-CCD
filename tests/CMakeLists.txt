################################################################################
# Tests
################################################################################

add_executable(ccd_tests
    main.cpp
    test_ccd.cpp
)

################################################################################
# Required Libraries
################################################################################

if(CCD_BUILD_RATIONAL)
    target_link_libraries(ccd_tests PUBLIC ccd::rational)
    target_compile_definitions(ccd_tests PUBLIC CCD_TEST_RATIONAL)
endif()
target_link_libraries(ccd_tests PUBLIC ccd::double)

include(ccd_warnings)
target_link_libraries(ccd_tests PRIVATE ccd::warnings)

include(catch2)
target_link_libraries(ccd_tests PUBLIC Catch2::Catch2)

################################################################################
# Add Tests
################################################################################

foreach(source IN ITEMS ${test_sources})
   source_group("tests" FILES "${source}")
endforeach()

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
# include(Catch)
catch_discover_tests(ccd_tests)
