add_executable(ex_nbody nbody.cpp)
#target_link_libraries(ex_nbody cxx_bench)
target_include_directories (ex_nbody PUBLIC ${PROJECT_SOURCE_DIR}/include)
