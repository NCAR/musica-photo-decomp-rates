cmake_minimum_required(VERSION 2.8)
project(MUSICA)
enable_language(Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

################################################################################
# NetCDF library

find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod"
  PATHS
    $ENV{NETCDF_HOME}/include
    /usr/lib/gfortran/modules
    /usr/lib64/gfortran/modules
    /opt/local/include)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS
    $ENV{NETCDF_HOME}/lib
    $ENV{NETCDF_HOME}/lib64
    /opt/local/lib)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS
    $ENV{NETCDF_HOME}/lib
    $ENV{NETCDF_HOME}/lib64
    /opt/local/lib)
  set(NETCDF_LIBS ${NETCDF_C_LIB} ${NETCDF_FORTRAN_LIB})
include_directories(${NETCDF_INCLUDE_DIR})

################################################################################
# copy photolysis input data to build directory

add_custom_target(copy_input_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/inputs ${CMAKE_BINARY_DIR}/inputs)

################################################################################
# musica-photo-decomp-rates library

add_subdirectory(src)
include_directories(${CMAKE_BINARY_DIR}/src)

################################################################################
# tests

enable_testing()
add_subdirectory(test)

################################################################################
