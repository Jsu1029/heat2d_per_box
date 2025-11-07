# makefile overrides
# OS:       macOS
# Compiler: gfortran 12.X
# OpenMP:   enabled
# BLAS:     framework acceralate
#
# NOTE for user:
#           Check gfortran version number 
#

CC=gcc-12
CXX=g++-12
FC=gfortran-12
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy 

ifeq ($(PREFIX),)
    FGT_INSTALL_DIR=/usr/local/lib
endif

ifeq ($(PREFIX_FINUFFT),)
    FINUFFT_INSTALL_DIR=/usr/local/lib
endif


ifeq ($(PREFIX_FFT),)
    FFT_INSTALL_DIR=/usr/local/lib
endif


ifeq ($(PREFIX_FFT_INCLUDE),)
    FFT_INCLUDE_DIR=/usr/local/include
endif

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

LBLAS=-framework accelerate



