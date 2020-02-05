include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(CCD_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(CCD_EXTRA_OPTIONS "")
endif()

function(custom_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${CCD_EXTERNAL}/${name}
        DOWNLOAD_DIR ${CCD_EXTERNAL}/.cache/${name}
        QUIET
        ${CCD_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

# libigl
function(download_libigl)
  custom_download_project(libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG        56f129e4403d3b7b04ddd786745fd3d574e95e04
  )
endfunction()

# Logger
function(download_spdlog)
    custom_download_project(spdlog
       GIT_REPOSITORY https://github.com/gabime/spdlog.git
       GIT_TAG        v1.5.0
    )
endfunction()

# Catch2 for testing
function(download_catch2)
    custom_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.11.1
    )
endfunction()

## Comparisons

# Etienne Vouga's CTCD Library
function(download_evctcd)
  custom_download_project(EVCTCD
    GIT_REPOSITORY https://github.com/evouga/collisiondetection.git
    GIT_TAG        e5fe5c9767207df5047e375fb20180a665ae186f
  )
endfunction()

# exact-ccd (clone of Brochu et al. [2012] and Tang et al. [2014])
function(download_exact_ccd)
  custom_download_project(exact-ccd
    GIT_REPOSITORY https://github.com/jiangzhongshi/exact-ccd.git
    GIT_TAG        305bb6f0e57d399b283161dc3669c260f90fb7f5
  )
endfunction()

# Rational CCD (rational version of Brochu et al. [2012])
function(download_rational_ccd)
  custom_download_project(rational_ccd
    GIT_REPOSITORY https://github.com/teseoch/Exact-CDD.git
    GIT_TAG        94b05129efa2c33348622b2890bd464cc95f65a7
  )
endfunction()
