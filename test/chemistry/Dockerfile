FROM fedora:41

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc \
        gcc-c++ \
        gcc-fortran \
        gdb \
        git \
        lapack-devel \
        lcov \
        m4 \
        make \
        netcdf-fortran-devel \
        python3 \
        python3-pip \
        valgrind \
    && dnf clean all

COPY . /cam

RUN cd /cam/test/chemistry \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make test

WORKDIR /cam/test/chemistry/build
