set(EXENAME example)
message(STATUS " ***** EXAMPLE ***** ")
# Sources
include_directories("${PROJECT_SOURCE_DIR}")
file(GLOB SOURCE_FILES *.cpp)

# Add defines for exe build
# add_definitions("-DXXX=")

message(STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR})

message(STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})


link_directories(${PROJECT_SOURCE_DIR}/lib)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17")
message(STATUS "Using CXX flags: " ${CMAKE_CXX_FLAGS})



add_executable(${EXENAME} ${SOURCE_FILES})

message(STATUS "INTERFACE_LINK_LIBRARIES: " ${INTERFACE_LINK_LIBRARIES} )
message(STATUS "LINK_INTERFACE_LIBRARIES: " ${LINK_INTERFACE_LIBRARIES} )
target_link_libraries(${EXENAME} astro)

