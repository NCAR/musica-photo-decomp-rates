# musica-photo-decomp-rates
Formerly known as TUV in MICM


To build the library and tests:

```
 mkdir build
 cd build
 ccmake ..
 make
```

To run the tests, from the `build/` folder:

```
 make test
```

Options are in test_nml


Notes for modeling2:

For GNU:

```
 export PATH="/opt/local/bin:$PATH"
 export LD_LIBRARY_PATH="/opt/local/lib64:/opt/local/lib"

 cmake3 -DCMAKE_CXX_COMPILER="/opt/local/bin/g++" \
        -DCMAKE_C_COMPILER="/opt/local/bin/gcc" \
        -DCMAKE_Fortran_COMPILER="/opt/local/bin/gfortran" \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_Fortran_FLAGS_DEBUG="-g -ggdb -ffpe-trap='invalid,zero,overflow' -finit-real=snan -fcheck=bounds" \
        ..

 cmake3 -DCMAKE_Fortran_COMPILER="/opt/local/bin/gfortran" ..

 make VERBOSE=1

 ctest -V

 ctest3 -V -R integration

```

For PGI:

```
 cmake3 -DCMAKE_C_COMPILER_ID="PGI" -DCMAKE_C_COMPILER="pgcc" \
        -DCMAKE_CXX_COMPILER_ID="PGI" -DCMAKE_CXX_COMPILER="pgc++" \
        -DCMAKE_Fortran_COMPILER_ID="PGI" -DCMAKE_Fortran_COMPILER="pgf90" \
        -DNETCDF_INCLUDE_DIR="/usr/local/netcdf-4.7.0/include" \
        -DNETCDF_C_LIB="/usr/local/netcdf-4.7.0/lib/libnetcdf.so" \
        -DNETCDF_FORTRAN_LIB="/usr/local/netcdf-4.7.0/lib/libnetcdff.so" \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_Fortran_FLAGS_DEBUG="-g -Ktrap=fp -Mbounds -Kieee -traceback" \
       ..


 env NETCDF_HOME="/usr/local/netcdf-4.7.0" cmake3 -DCMAKE_Fortran_COMPILER="pgf90" -DCMAKE_BUILD_TYPE=Debug ..

```
