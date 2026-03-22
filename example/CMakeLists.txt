add_executable(example example.cpp)
target_link_libraries(example PRIVATE astro)

add_executable(propagate propagate.cpp)
target_link_libraries(propagate PRIVATE astro)
