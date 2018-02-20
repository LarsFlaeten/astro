set(EXENAME example)

# Sources
include_directories("${PROJECT_SOURCE_DIR}")
file(GLOB SOURCE_FILES *.cpp)

# Add defines for exe build
# add_definitions("-DXXX=")

message(STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR})

message(STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})


link_directories(${PROJECT_SOURCE_DIR}/lib)


add_executable(${EXENAME} ${SOURCE_FILES})
target_compile_options(${EXENAME} PUBLIC -Wall)
message(STATUS "INTERFACE_LINK_LIBRARIES: " ${INTERFACE_LINK_LIBRARIES} )
message(STATUS "LINK_INTERFACE_LIBRARIES: " ${LINK_INTERFACE_LIBRARIES} )
target_link_libraries(${EXENAME} astro)

