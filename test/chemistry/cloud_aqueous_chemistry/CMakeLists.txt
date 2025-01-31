# Cloud aqueous chemistry tests
add_pfunit_ctest(test_cloud_aqueous_chemistry
                 TEST_SOURCES
                 test_cloud_aqueous_chemistry.pf
                 cloud_chemistry_dependencies.F90
                 ../../../src/chemistry/aerosol/cloud_aqueous_chemistry.F90
                 ../../../src/chemistry/aerosol/cloud_utilities.F90
                 ../../../src/chemistry/modal_aero/mam_clouds.F90)
