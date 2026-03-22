cmake_minimum_required(VERSION 3.15)
project(astro VERSION 0.1.0 LANGUAGES CXX C)

# Prevent in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not permitted.\n"
        "  mkdir build && cd build && cmake ..")
endif()

# Default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Defaulting to RelWithDebInfo build type.")
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# CSPICE — bundled static library
add_library(cspice STATIC IMPORTED GLOBAL)
set_target_properties(cspice PROPERTIES
    IMPORTED_LOCATION
        "${CMAKE_CURRENT_SOURCE_DIR}/libraries/cspice/libcspice.a"
    INTERFACE_INCLUDE_DIRECTORIES
        "${CMAKE_CURRENT_SOURCE_DIR}/libraries"
)

# GLM — prefer system install, fall back to FetchContent.
# Some system installs create 'glm' instead of 'glm::glm'; normalise to the
# namespaced form so the rest of the build can always use glm::glm.
find_package(glm QUIET)
if(NOT TARGET glm::glm)
    if(TARGET glm)
        # Older install without namespace — create an alias
        add_library(glm::glm ALIAS glm)
        message(STATUS "GLM found (aliased from 'glm' to 'glm::glm')")
    else()
        message(STATUS "GLM not found, fetching via FetchContent")
        include(FetchContent)
        FetchContent_Declare(
            glm
            GIT_REPOSITORY https://github.com/g-truc/glm.git
            GIT_TAG        1.0.1
        )
        FetchContent_MakeAvailable(glm)
    endif()
endif()

add_subdirectory(astro)

# Tests and examples are only built when astro is the top-level project,
# not when consumed as a submodule.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    enable_testing()
    add_subdirectory(test)
    add_subdirectory(example)
endif()
