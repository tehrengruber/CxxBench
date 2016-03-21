cmake_minimum_required (VERSION 2.8)
project (CxxBench)

find_package(Boost REQUIRED)

include_directories("${Boost_INCLUDE_DIRS}")

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-march=native")

# build examples
add_subdirectory("examples")

INSTALL(DIRECTORY include/
	DESTINATION "include"
)