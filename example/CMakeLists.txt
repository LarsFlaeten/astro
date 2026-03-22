# cxxopts — prefer system install, fall back to downloading the single header.
# We intentionally avoid FetchContent_MakeAvailable here: cxxopts' own
# CMakeLists adds an executable also named "example" which would conflict.
find_package(cxxopts QUIET)
if(NOT TARGET cxxopts::cxxopts)
    include(FetchContent)
    FetchContent_Declare(
        cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG        v2.2.0
    )
    FetchContent_GetProperties(cxxopts)
    if(NOT cxxopts_POPULATED)
        FetchContent_Populate(cxxopts)   # download only — do NOT call add_subdirectory
    endif()
    add_library(cxxopts::cxxopts INTERFACE IMPORTED GLOBAL)
    target_include_directories(cxxopts::cxxopts INTERFACE "${cxxopts_SOURCE_DIR}/include")
endif()

add_executable(example example.cpp)
target_link_libraries(example PRIVATE astro cxxopts::cxxopts)

add_executable(propagate propagate.cpp)
target_link_libraries(propagate PRIVATE astro cxxopts::cxxopts)
