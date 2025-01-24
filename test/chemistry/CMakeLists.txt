cmake_minimum_required(VERSION 3.5)

project(
  cam_chemistry_tests
  VERSION 0.1.0
  LANGUAGES Fortran
)

include(FetchContent)

FetchContent_Declare(
  pfunit
  GIT_REPOSITORY https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
  GIT_TAG main
)

FetchContent_MakeAvailable(pfunit)

enable_testing()

add_subdirectory(cloud_aqueous_chemistry)