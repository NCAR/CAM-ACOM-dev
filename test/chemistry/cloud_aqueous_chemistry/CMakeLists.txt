################################################################################
# Test utilities

include(${CMAKE_SOURCE_DIR}/cmake/test_util.cmake)

################################################################################
# New cloud aqueous chemistry code tests

create_standard_test(
  NAME cloud_aqueous_chemistry_MAM
  SOURCES
    test_cloud_aqueous_chemistry_MAM.F90
    cloud_chemistry_dependencies_MAM.F90
    ../file_io.F90
    ../../../src/chemistry/aerosol/cloud_aqueous_chemistry.F90
    ../../../src/chemistry/aerosol/cloud_utilities.F90
    ../../../src/chemistry/bulk_aero/bam_clouds.F90
    ../../../src/chemistry/modal_aero/mam_clouds.F90
    ../../../src/chemistry/carma_aero/carma_clouds.F90
    ../../../src/chemistry/aerosol/mo_setsox.F90
    ../../../src/chemistry/aerosol/cldaero_mod.F90
    ../../../src/chemistry/modal_aero/sox_cldaero_mod.F90
  CPP_FLAGS
    -DUSE_MAM
  WORKING_DIRECTORY
    ${CMAKE_BINARY_DIR}/cloud_aqueous_chemistry)

create_standard_test(
  NAME cloud_aqueous_chemistry_BAM
  SOURCES
    test_cloud_aqueous_chemistry_BAM.F90
    cloud_chemistry_dependencies_BAM.F90
    ../file_io.F90
    ../../../src/chemistry/aerosol/cloud_aqueous_chemistry.F90
    ../../../src/chemistry/aerosol/cloud_utilities.F90
    ../../../src/chemistry/bulk_aero/bam_clouds.F90
    ../../../src/chemistry/modal_aero/mam_clouds.F90
    ../../../src/chemistry/carma_aero/carma_clouds.F90
    ../../../src/chemistry/aerosol/mo_setsox.F90
    ../../../src/chemistry/aerosol/cldaero_mod.F90
    ../../../src/chemistry/bulk_aero/sox_cldaero_mod.F90
  CPP_FLAGS
    -DUSE_BAM
  WORKING_DIRECTORY
    ${CMAKE_BINARY_DIR}/cloud_aqueous_chemistry)

create_standard_test(
  NAME cloud_aqueous_chemistry_CARMA
  SOURCES
    test_cloud_aqueous_chemistry_CARMA.F90
    cloud_chemistry_dependencies_CARMA.F90
    ../file_io.F90
    ../../../src/chemistry/aerosol/cloud_aqueous_chemistry.F90
    ../../../src/chemistry/aerosol/cloud_utilities.F90
    ../../../src/chemistry/bulk_aero/bam_clouds.F90
    ../../../src/chemistry/modal_aero/mam_clouds.F90
    ../../../src/chemistry/carma_aero/carma_clouds.F90
    ../../../src/chemistry/aerosol/mo_setsox.F90
    ../../../src/chemistry/aerosol/cldaero_mod.F90
    ../../../src/chemistry/carma_aero/sox_cldaero_mod.F90
  CPP_FLAGS
    -DUSE_CARMA
  WORKING_DIRECTORY
    ${CMAKE_BINARY_DIR}/cloud_aqueous_chemistry)

create_standard_test(
  NAME cloud_aqueous_chemistry_full
  SOURCES
    test_cloud_aqueous_chemistry_full.F90
    cloud_chemistry_dependencies_full.F90
    ../file_io.F90
    ../../../src/chemistry/aerosol/cloud_aqueous_chemistry.F90
    ../../../src/chemistry/aerosol/cloud_utilities.F90
    ../../../src/chemistry/bulk_aero/bam_clouds.F90
    ../../../src/chemistry/modal_aero/mam_clouds.F90
    ../../../src/chemistry/carma_aero/carma_clouds.F90
    ../../../src/chemistry/aerosol/mo_setsox.F90
    ../../../src/chemistry/aerosol/cldaero_mod.F90
    ../../../src/chemistry/modal_aero/sox_cldaero_mod.F90
  CPP_FLAGS
    -DUSE_MAM
  WORKING_DIRECTORY
    ${CMAKE_BINARY_DIR}/cloud_aqueous_chemistry)
