# makefile overrides
# OS:       macOS
# Compiler: gfortran 13.X
# OpenMP:   enabled
# BLAS:     framework acceralate
# FFTW:     brew installed
#
# NOTE for user:
#           Check gfortran version number
#

CC=gcc-13
CXX=g++-13
FC=gfortran-13
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy

ifeq ($(PREFIX),)
    FGT_INSTALL_DIR=${HOME}/lib
endif

ifeq ($(PREFIX_FINUFFT),)
    FINUFFT_INSTALL_DIR=${HOME}/lib
endif


ifeq ($(PREFIX_FFT),)
    FFTW_ROOT = $(shell brew --prefix fftw)
    FFT_INSTALL_DIR=${FFTW_ROOT}/lib
endif


ifeq ($(PREFIX_FFT_INCLUDE),)
    FFTW_ROOT = $(shell brew --prefix fftw)
    FFT_INCLUDE_DIR=${FFTW_ROOT}/include
endif

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

LBLAS=-framework accelerate
