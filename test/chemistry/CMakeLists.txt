cmake_minimum_required(VERSION 3.6)

project(
  cam_chemistry_tests
  VERSION 0.1.0
  LANGUAGES Fortran
)

option(CAMCHEM_ENABLE_MEMCHECK "Enable memory checking with valgrind" ON)
option(CAMCHEM_ENABLE_MPI "Build with MPI support" OFF)
option(CAMCHEM_ENABLE_OPENMP "Build with OpenMP support" OFF)

include(cmake/dependencies.cmake)

# Copy test data
add_custom_target(copy_test_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/test/chemistry/data)

enable_testing()

add_subdirectory(cloud_aqueous_chemistry)