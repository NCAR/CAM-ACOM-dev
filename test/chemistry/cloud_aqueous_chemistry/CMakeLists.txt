################################################################################
# Test utilities

include(${CMAKE_SOURCE_DIR}/cmake/test_util.cmake)

################################################################################
# New cloud aqueous chemistry code tests

create_standard_test(
  NAME cloud_aqueous_chemistry
  SOURCES
    test_cloud_aqueous_chemistry.F90
    cloud_chemistry_dependencies.F90
    ../file_io.F90
    ../../../src/chemistry/aerosol/cloud_aqueous_chemistry.F90
    ../../../src/chemistry/aerosol/cloud_utilities.F90
    ../../../src/chemistry/modal_aero/mam_clouds.F90
    ../../../src/chemistry/aerosol/mo_setsox.F90
    ../../../src/chemistry/aerosol/cldaero_mod.F90
    ../../../src/chemistry/modal_aero/sox_cldaero_mod.F90
  WORKING_DIRECTORY
    ${CMAKE_BINARY_DIR}/cloud_aqueous_chemistry)
