add_library(astro SHARED
    State.cpp
    Time.cpp
    Util.cpp
    SpiceCore.cpp
    ReferenceFrame.cpp
    Observer.cpp
    Orbit.cpp
    OrbitElements.cpp
    ODE.cpp
    Interpolate.cpp
    PCDM.cpp
    RKF45.cpp
    RKF78.cpp
)

target_compile_features(astro PUBLIC cxx_std_17)
target_compile_options(astro PRIVATE -Wall)

# Public include path: the parent of this directory, so consumers use
# #include "astro/State.h". The PRIVATE entry lets our own .cpp files
# use #include "State.h" without qualification.
target_include_directories(astro
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(astro
    PUBLIC  glm::glm
    PRIVATE -Wl,--whole-archive cspice -Wl,--no-whole-archive
)

# Install
include(GNUInstallDirs)
install(TARGETS astro
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES
    Exceptions.h
    Math.h
    State.h
    Time.h
    Util.h
    SpiceCore.h
    ReferenceFrame.h
    Observer.h
    Orbit.h
    OrbitElements.h
    ODE.h
    Interpolate.h
    PCDM.h
    Propagator.h
    RK1_4.h
    RKF45.h
    RKF78.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/astro
)
