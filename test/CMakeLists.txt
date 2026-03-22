# GTest — prefer system install, fall back to FetchContent
find_package(GTest CONFIG QUIET)
if(NOT GTest_FOUND)
    message(STATUS "GTest not found via find_package, fetching via FetchContent")
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.14.0
    )
    # Prevent GoogleTest from overriding our compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

add_executable(runTests
    tests_main.cpp
    testTime.cpp
    testState.cpp
    testUtil.cpp
    testSpiceCore.cpp
    testReferenceFrame.cpp
    testObserver.cpp
    testOrbit.cpp
    testOrbitElements.cpp
    testODE.cpp
    testInterpolate.cpp
    testNumInt.cpp
    testPCDM.cpp
)

target_link_libraries(runTests PRIVATE
    astro
    GTest::gtest
    GTest::gmock
)

include(GoogleTest)
gtest_discover_tests(runTests WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
