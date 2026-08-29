###
# Copyright (c) 2012-2018, University of Strathclyde
# License: BSD-3-Clause
# 25.6.2015 Jonathan Smith - CMake lists for PUFFIN Project
###

cmake_minimum_required (VERSION 3.20)
project (puffin)

enable_language (Fortran)

# The version number.
set (Puffin_VERSION_MAJOR 2)
set (Puffin_VERSION_MINOR 1)
set (Puffin_VERSION_PATCH 0a)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

option(ENABLE_PARALLEL "whether to build parallel PUFFIN" ON)
option(USE_MKL "Use Intel MKL" OFF)
option(ENABLE_TESTING "Enable testing" OFF)
option(ENABLE_BENCHMARKS "Enable performance benchmarks (see benchmark/README.md)" OFF)
option(ENABLE_LINT "Enable extra compiler warnings (-Wall -Wextra)" OFF)
option(USE_LLVM_OPENMP "macOS + gfortran only: link LLVM's libomp instead of GCC's libgomp. Opt-in; see the USE_LLVM_OPENMP block below before enabling" OFF)

######################################################################
#
# Find packages
#
######################################################################

if(ENABLE_PARALLEL)
  find_package(MPI REQUIRED)
  set(USE_MPI 1 CACHE BOOL "setting USE_MPI to 1" FORCE)
endif()

if (USE_MKL)
  # User is expected to provide MKLROOT or MKL_LIBRARY_DIRS
  if(NOT DEFINED MKL_LIBRARY_DIRS AND DEFINED ENV{MKLROOT})
    set(MKL_LIBRARY_DIRS "$ENV{MKLROOT}/lib/intel64")
    set(MKL_INCLUDE_DIRS "$ENV{MKLROOT}/include")
  endif()
  message(STATUS "Using MKL. MKL_LIBRARY_DIRS: ${MKL_LIBRARY_DIRS}")
else()
  find_package(FFTW3 REQUIRED)
endif()

find_package(HDF5 COMPONENTS Fortran HL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Doxygen)

if (DOXYGEN_FOUND)
    message(STATUS "Doxygen FOUND!!")
    message(STATUS "DOXYGEN IS : ${DOXYGEN_EXECUTABLE}")
endif()

######################################################################
#
# OpenMP runtime selection
#
# Default (USE_LLVM_OPENMP=OFF): whatever the compiler's own -fopenmp pulls
# in -- libgomp for GCC. This is the portable path and the one CI uses.
#
# USE_LLVM_OPENMP=ON links LLVM's libomp instead. This exists because
# Homebrew's libgomp on macOS has no futex: it puts idle threads to sleep
# between barriers and pays a kernel wake-up on every one, and it ignores
# both GOMP_SPINCOUNT and OMP_WAIT_POLICY, so there is no way to make it
# spin. Measured barrier cost on an M-series Mac, 200k barriers:
#
#   runtime                            2 threads    4 threads
#   GCC libgomp                           5613 ns     14407 ns
#   LLVM libomp (default blocktime)       3883 ns      8336 ns
#   LLVM libomp + KMP_BLOCKTIME=infinite   197 ns       283 ns
#
# gfortran-compiled objects link against libomp unchanged, because libomp
# implements the GOMP_* ABI -- no recompilation, only a different runtime.
#
# IMPORTANT: libomp is only worth linking if you also set
#
#   export KMP_BLOCKTIME=infinite
#
# at run time. Without it you get the middle row above, which is not enough
# to make threading pay. With it, and OMP_NUM_THREADS=2, the 1D benchmarks
# go from 10.92 s to 7.24 s (1d_serial) and 8.12 s to 6.75 s (1d_mpi). The
# 3D cases are still slower threaded than serial -- they carry 16 !$OMP
# ATOMIC updates per macroparticle in getSource_3D plus roughly twice the
# barriers -- so this is not yet a reason to thread 3D runs.
#
# Also note threaded runs are NOT bit-reproducible: the atomic scatter in
# getSource makes summation order thread-dependent. Divergence starts at
# machine epsilon and grows with step count (worst L2rel 1.3e-12 over the
# 6000-step 1D deck, but 5.5e-10 max relative on signal-carrying nodes by
# the final dump). The E2E suite asserts to 1e-10, so it must keep running
# at OMP_NUM_THREADS=1.
#
# This is a macOS-specific workaround. On Linux libgomp has futex-backed
# barriers and none of the above applies; do not enable it there.
#
######################################################################

if (USE_LLVM_OPENMP)
    if (NOT APPLE)
        message(FATAL_ERROR
            "USE_LLVM_OPENMP is a macOS-only workaround for Homebrew libgomp's "
            "sleeping barriers. On this platform the compiler's own OpenMP "
            "runtime is fine -- leave it OFF.")
    endif()
    if (NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR
            "USE_LLVM_OPENMP assumes gfortran (it relies on libomp providing the "
            "GOMP_* ABI that gfortran emits). Fortran compiler here is "
            "'${CMAKE_Fortran_COMPILER_ID}' -- leave it OFF.")
    endif()

    find_library(LIBOMP_LIBRARY
        NAMES omp
        HINTS ${LIBOMP_ROOT}/lib
              /opt/homebrew/opt/libomp/lib
              /usr/local/opt/libomp/lib
        DOC "LLVM OpenMP runtime (brew install libomp)")

    if (NOT LIBOMP_LIBRARY)
        message(FATAL_ERROR
            "USE_LLVM_OPENMP=ON but libomp was not found. Install it with "
            "'brew install libomp', or point LIBOMP_ROOT at its prefix.")
    endif()

    # -fopenmp must reach the COMPILER but not the LINKER: on the link line it
    # would drag in libgomp alongside libomp, and two OpenMP runtimes in one
    # process is undefined behaviour. add_compile_options is compile-only,
    # unlike CMAKE_Fortran_FLAGS which is used for both.
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>)

    # Repoint the imported target too, not just our own link lines. Dependencies
    # carry OpenMP::OpenMP_Fortran in their INTERFACE_LINK_LIBRARIES -- pFUnit
    # does -- and that resolves to libgomp, so without this the test binaries
    # end up linked against BOTH runtimes even though puffin itself is clean.
    #
    # This has to be a macro rather than a one-shot edit: find_package's
    # imported targets are directory-scoped, and PFUNITConfig.cmake calls
    # find_dependency(OpenMP) itself, which recreates the target pointing back
    # at libgomp in test/'s scope. So test/CMakeLists.txt calls this again
    # after its find_package(PFUNIT).
    macro(puffin_force_llvm_openmp)
        if (USE_LLVM_OPENMP AND TARGET OpenMP::OpenMP_Fortran)
            set_target_properties(OpenMP::OpenMP_Fortran PROPERTIES
                INTERFACE_LINK_LIBRARIES "${LIBOMP_LIBRARY}")
        endif()
    endmacro()

    puffin_force_llvm_openmp()

    message(STATUS "OpenMP runtime : LLVM libomp (${LIBOMP_LIBRARY})")
    message(STATUS "                 remember to set KMP_BLOCKTIME=infinite at run time")
elseif (OpenMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    message(STATUS "OpenMP runtime : compiler default (${CMAKE_Fortran_COMPILER_ID})")
endif()

if (ENABLE_LINT)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # -Wno-unused-dummy-argument: many routines share a common call signature
    # (lattice element interfaces, callback-style routines) where a given
    # implementation legitimately ignores some arguments.
    # -Wno-compare-reals: exact REAL(8) ==/ /= comparisons are used
    # deliberately (sentinel values, zero guards) in several places; see
    # doc/LINT_NOTES.md for the full list. Rewriting them risks changing
    # numerical output in a codebase with a bit-exact E2E test.
    # -Wno-aliasing: transforms.f90 passes the same buffer as FFTW's `in`
    # and `out` arguments deliberately, for an in-place transform; see
    # doc/LINT_NOTES.md.
    # -Wno-intrinsic-shadow: error_fn.f90 deliberately defines its own
    # erf/erfc, shadowing the Fortran 2008 intrinsics of the same name;
    # see doc/LINT_NOTES.md.
    # -Wno-conversion: legacy numerical/HDF5 code relies on implicit
    # narrowing (REAL->INTEGER truncation, INTEGER(8)->INTEGER(4)) as a
    # deliberate idiom throughout; see doc/LINT_NOTES.md.
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wno-unused-dummy-argument -Wno-compare-reals -Wno-aliasing -Wno-intrinsic-shadow -Wno-conversion")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all")
  else()
    message(WARNING "ENABLE_LINT requested but no warning flags known for Fortran compiler '${CMAKE_Fortran_COMPILER_ID}'")
  endif()
endif()

execute_process(COMMAND git --work-tree=${CMAKE_SOURCE_DIR} rev-parse --abbrev-ref HEAD
  RESULT_VARIABLE GIT_RESULT
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

if (NOT GIT_BRANCH)
  set(GIT_BRANCH "unknown branch, not in git")
endif()

execute_process(COMMAND git --work-tree=${CMAKE_SOURCE_DIR} rev-parse HEAD
  RESULT_VARIABLE GIT_RESULT
  OUTPUT_VARIABLE GIT_REVISION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

if (NOT GIT_REVISION)
  set(GIT_REVISION "unknown revision, not in git")
endif()

message(STATUS "git branch : ${GIT_BRANCH} ${GIT_REVISION}")
string(TIMESTAMP CONFIGURE_TIME_STAMP "%Y-%m-%d %H:%M")

######################################################################
#
# Variables for inclusion in config.h
#
######################################################################

set(CXX ${CMAKE_CXX_COMPILER})
set(CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CBT_UC}}")
set(CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CBT_UC}}")
set(HAVE_CMAKE 1 CACHE STRING "Whether built with CMake")

if(ENABLE_TESTING OR ENABLE_BENCHMARKS)
  enable_testing()
endif()

######################################################################
#
# Create the configure file
#
######################################################################

configure_file(config-cmake.h.in config.h)

configure_file(puffin/lib/common/puffProvenance.f90.in puffin/lib/common/puffProvenance.f90)
# configure_file(doc/fortran.dxg.in doc/fortran.dxg)

add_subdirectory(puffin)

if (ENABLE_BENCHMARKS)
  add_subdirectory(benchmark)
endif (ENABLE_BENCHMARKS)

if (ENABLE_TESTING)
  add_subdirectory(test)
endif (ENABLE_TESTING)

if (DOXYGEN_FOUND)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/fortran.dxg
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif()

######################################################################
#
# Static lint target (opt-in, not part of the default build)
#
######################################################################

find_program(FORTITUDE_EXECUTABLE NAMES fortitude)

if (FORTITUDE_EXECUTABLE)
    add_custom_target(lint
        ${FORTITUDE_EXECUTABLE} check ${CMAKE_SOURCE_DIR}/puffin
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running fortitude static analysis on puffin/" VERBATIM
        )
else()
    add_custom_target(lint
        ${CMAKE_COMMAND} -E echo
            "fortitude not found. Install it (e.g. 'pip install fortitude-lint' or 'cargo install fortitude-lint'), re-run cmake, then 'make lint'."
        COMMENT "fortitude not found" VERBATIM
        )
endif()

install (DIRECTORY utilities/post DESTINATION bin)
install (DIRECTORY utilities/pyPlotting DESTINATION bin)
install (DIRECTORY utilities/setup DESTINATION bin)
install (DIRECTORY utilities/visit-scripts DESTINATION bin)

install (DIRECTORY inputs/ DESTINATION examples)
# install (DIRECTORY inputs/gen/ DESTINATION python)

install (DIRECTORY scripts/ DESTINATION examples/job-submission-scripts)

######################################################################
#
# Finalize (installation, packaging)
#
######################################################################

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Puffin")

# build a CPack driven installer package

include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${Puffin_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${Puffin_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${Puffin_VERSION_PATCH}")
#include (CPack)
