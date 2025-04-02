################################################################################
# Utility functions for creating tests

if(CAMCHEM_ENABLE_MEMCHECK)
  find_program(MEMORYCHECK_COMMAND "valgrind")
endif()

################################################################################
# build and add a standard test

function(create_standard_test)
  set(prefix TEST)
  set(singleValues NAME WORKING_DIRECTORY)
  set(multiValues SOURCES CPP_FLAGS)
  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix} " " "${singleValues}" "${multiValues}" ${ARGN})
  add_executable(test_${TEST_NAME} ${TEST_SOURCES})
  if(TEST_CPP_FLAGS)
    target_compile_options(test_${TEST_NAME} PRIVATE ${TEST_CPP_FLAGS})
  endif()
  set_target_properties(test_${TEST_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(test_${TEST_NAME} PRIVATE PkgConfig::netcdff PkgConfig::netcdfc)
  if(CAMCHEM_ENABLE_OPENMP)
    target_link_libraries(test_${TEST_NAME} PUBLIC OpenMP::OpenMP_Fortran)
  endif()
  if(NOT DEFINED TEST_WORKING_DIRECTORY)
    set(TEST_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  endif()
  add_camchem_test(${TEST_NAME} test_${TEST_NAME} "" ${TEST_WORKING_DIRECTORY})
endfunction(create_standard_test)

################################################################################
# Add a test

function(add_camchem_test test_name test_binary test_args working_dir)
  if(CAMCHEM_ENABLE_MPI)
    add_test(NAME ${test_name}
      COMMAND mpirun -v -np 2 ${CMAKE_BINARY_DIR}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  else()
    add_test(NAME ${test_name}
             COMMAND ${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  endif()
  set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --trace-children=yes --leak-check=full -s --gen-suppressions=all ${MEMCHECK_SUPPRESS}")
  set(memcheck "${MEMORYCHECK_COMMAND} ${MEMORYCHECK_COMMAND_OPTIONS}")
  separate_arguments(memcheck)
  if(CAMCHEM_ENABLE_MPI AND MEMORYCHECK_COMMAND AND CAMCHEM_ENABLE_MEMCHECK)
    add_test(NAME memcheck_${test_name}
      COMMAND mpirun -v -np 2 ${memcheck} ${CMAKE_BINARY_DIR}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  elseif(MEMORYCHECK_COMMAND AND CAMCHEM_ENABLE_MEMCHECK)
    add_test(NAME memcheck_${test_name}
             COMMAND ${memcheck} ${working_dir}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  endif()
endfunction(add_camchem_test)