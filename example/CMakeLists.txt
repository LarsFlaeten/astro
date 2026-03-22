# cxxopts — prefer system install, fall back to FetchContent
find_package(cxxopts QUIET)
if(NOT TARGET cxxopts::cxxopts)
    include(FetchContent)
    FetchContent_Declare(
        cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG        v2.2.1
    )
    FetchContent_MakeAvailable(cxxopts)
endif()

add_executable(example example.cpp)
target_link_libraries(example PRIVATE astro cxxopts::cxxopts)

add_executable(propagate propagate.cpp)
target_link_libraries(propagate PRIVATE astro cxxopts::cxxopts)
