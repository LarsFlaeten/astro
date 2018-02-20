set(LIBNAME "astro")

# Sources
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/libraries" )
file(GLOB SOURCE_FILES *.cpp */*.cpp)

# Libraries
set(SPICE_LIBS cspice)
#if(UNIX)
#	set(LIBS ${LIBS} rt)
#endif(UNIX)

message(STATUS "Using libs:  " ${LIBS})
message(STATUS "Include dirs:" ${PROJECT_SOURCE_DIR} ", " ${PROJECT_SOURCE_DIR}/libraries)

link_directories(${PROJECT_SOURCE_DIR}/libraries ${PROJECT_SOURCE_DIR}/libraries/cspice)

# Static or shared?
set(LIBTYPE SHARED)


# Create a library
add_library(${LIBNAME} ${LIBTYPE} ${SOURCE_FILES})
# We link cspice as PRIVATE, so the dependency is not propagated to the executable using this library
target_link_libraries(${LIBNAME} PRIVATE -Wl,--whole-archive ${SPICE_LIBS} -Wl,--no-whole-archive)

# Adds SO Version and subversion. To be added to ensure ABI/API compatibility.
#SET_TARGET_PROPERTIES(${LIBNAME} PROPERTIES SOVERSION ${ORK_VERSION_MAJOR} VERSION ${ORK_VERSION})

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

# Install headers
#FOREACH(subdir core math render resource scenegraph taskgraph ui util)
FILE(GLOB include-files *.h)
INSTALL(FILES ${include-files} DESTINATION include/${LIBNAME})
#INSTALL(FILES ${include-files} DESTINATION include/ork/${subdir})
#ENDFOREACH(subdir)

# Install library
INSTALL(TARGETS ${LIBNAME} LIBRARY DESTINATION ${LIB_INSTALL_DIR} ARCHIVE DESTINATION ${LIB_INSTALL_DIR})

