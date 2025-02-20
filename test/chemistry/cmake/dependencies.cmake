find_package(PkgConfig REQUIRED)
include(FetchContent)

# ##############################################################################
# Memory check

if(CAMCHEM_ENABLE_MEMCHECK)
  find_file(
    MEMCHECK_SUPPRESS_FILE
    DOC "Suppression file for memory checking"
    NAMES openmpi-valgrind.supp
    PATHS /usr/share/openmpi /usr/lib64/openmpi/share
          /usr/lib64/openmpi/share/openmpi /usr/share)
  if(MEMCHECK_SUPPRESS_FILE)
    set(MEMCHECK_SUPPRESS
        "--suppressions=${PROJECT_SOURCE_DIR}/valgrind.supp --suppressions=${MEMCHECK_SUPPRESS_FILE}"
    )
  else()
    set(MEMCHECK_SUPPRESS
        "--suppressions=${PROJECT_SOURCE_DIR}/valgrind.supp")
  endif()
endif()

# ##############################################################################
# NetCDF library

find_package(PkgConfig REQUIRED)

pkg_check_modules(netcdff IMPORTED_TARGET REQUIRED netcdf-fortran)
pkg_check_modules(netcdfc IMPORTED_TARGET REQUIRED netcdf)
