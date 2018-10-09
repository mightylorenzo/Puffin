# start from base
FROM ubuntu:16.04
MAINTAINER Lawrence Campbell <mightylorenzo@gmail.com>

# install system-wide deps for Puffin
RUN apt-get -yqq update && \
             apt -y install gfortran \
             libhdf5-mpi-dev \
             libhdf5-openmpi-dev \
             libfftw3-dev \
             libfftw3-mpi-dev \
             cmake \
             libopenmpi-dev \
             g++

# NOTE - 'WORKDIR' is essentially just 'cd'

# Create user for running Puffin


RUN useradd --create-home puffin_user
RUN mkdir /home/puffin_user/tmp/
RUN mkdir /home/puffin_user/tmp/puffin-src
RUN mkdir /home/puffin_user/tmp/puffin-build
RUN mkdir /home/puffin_user/tmp/puffin-test
RUN mkdir /home/puffin_user/tmp/pFUnit
RUN mkdir /home/puffin_user/built/
RUN mkdir /home/puffin_user/built/puffin
RUN mkdir /home/puffin_user/built/pfunit-parallel
RUN mkdir /home/puffin_user/project # for running in Puffin
COPY . /home/puffin_user/tmp/puffin-src
RUN chown -R puffin_user /home/puffin_user
# Grant sudo access without password
RUN echo 'puffin_user ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers


# Switch to Puffin user to compile Puffin
USER puffin_user

WORKDIR /home/puffin_user/tmp/puffin-build

ENV PATH="/usr/include/hdf5/openmpi:${PATH}"

# Run CMake:
RUN cmake -DCMAKE_INSTALL_PREFIX:PATH=/home/puffin_user/built/puffin \
          -DENABLE_PARALLEL:BOOL=TRUE \
          -DHdf5_MODULE_DIRS='/usr/include/hdf5/openmpi' \
          -DHdf5_LIBRARY_DIRS='/usr/lib/x86_64-linux-gnu/hdf5/openmpi;/usr/lib/x86_64-linux-gnu/' \
          -DHdf5_INCLUDE_DIRS='/usr/include/hdf5/openmpi' \
          -DHdf5_LIBRARY_NAMES='hdf5_openmpi_fortran;hdf5_openmpi' \
          -DHdf5_LIBRARIES='/usr/lib/x86_64-linux-gnu/libhdf5_openmpi_fortran.so;/usr/lib/x86_64-linux-gnu/libhdf5_openmpi.so' \
          -DHdf5_STLIBS='/usr/lib/x86_64-linux-gnu/libhdf5_openmpi_fortran.a;/usr/lib/x86_64-linux-gnu/libhdf5_openmpi.a' \
          /home/puffin_user/tmp/puffin-src

RUN make && make install

ENV OMP_NUM_THREADS=1
RUN mpirun -np 4 /home/puffin_user/tmp/puffin-build/test/testexe
WORKDIR /home/puffin_user/project

WORKDIR /home/puffin_user/tmp/puffin-test
RUN cp /home/puffin_user/built/puffin/examples/simple/1D/OptCommV165pp65-70/fig1/*.in .
RUN chmod u+x /home/puffin_user/tmp/puffin-src/docker-entrypoint.sh

ENTRYPOINT ["/home/puffin_user/tmp/puffin-src/docker-entrypoint.sh"]
CMD ["2", "f1main.in"]

