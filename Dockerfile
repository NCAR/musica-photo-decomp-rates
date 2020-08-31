FROM fedora:33

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        gcc-c++ \
        netcdf-fortran-devel \
        cmake \
        make \
        diffutils \
    && dnf clean all

# copy the musica-photo-decomp-rates code
COPY . /musica-photo-decomp-rates/

# build the library and tests
RUN mkdir /build \
    && cd /build \
    && cmake -DCMAKE_BUILD_TYPE=Debug \
             -DCMAKE_Fortran_FLAGS_DEBUG="-g -ggdb -ffpe-trap='invalid,zero,overflow' -finit-real=snan -fcheck=bounds" \
             ../musica-photo-decomp-rates \
    && make
