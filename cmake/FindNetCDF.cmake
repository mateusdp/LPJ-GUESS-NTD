# - Try to find NetCDF
# Once done this will define
#  NETCDF_FOUND - System has NetCDF
#  NETCDF_INCLUDE_DIRS - The NetCDF include directories
#  NETCDF_LIBRARIES - The libraries needed to use NetCDF

find_path(NETCDF_INCLUDE_DIR netcdf.h)

find_library(NETCDF_C_LIBRARY netcdf)
find_library(NETCDF_CXX_LIBRARY netcdf_c++)

set(NETCDF_LIBRARIES ${NETCDF_C_LIBRARY})

# The C++ API currently doesn't come with the prebuilt binaries
# on Windows, so only include the C++ lib if found.
if(NETCDF_CXX_LIBRARY)
  set(NETCDF_LIBRARIES ${NETCDF_CXX_LIBRARY} ${NETCDF_LIBRARIES})
endif()

set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NetCDF  DEFAULT_MSG
  NETCDF_C_LIBRARY NETCDF_INCLUDE_DIR)

mark_as_advanced(NETCDF_INCLUDE_DIR NETCDF_C_LIBRARY NETCDF_CXX_LIBRARY)
