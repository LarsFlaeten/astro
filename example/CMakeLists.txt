# Sources
include_directories("${PROJECT_SOURCE_DIR}")


message(STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR})

message(STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})


link_directories(${PROJECT_SOURCE_DIR}/lib)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17")
message(STATUS "Using CXX flags: " ${CMAKE_CXX_FLAGS})



add_executable(example example.cpp)
target_link_libraries(example astro)



add_executable(propagate propagate.cpp)
target_link_libraries(propagate astro)
