# Makefile for fgt
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create
# the file make.inc, which overrides the defaults below (which are
# for ubunutu linux/gcc system).

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
CLINK = -lstdc++
FLINK = $(CLINK)

FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy -w
# -pg -no-pie is for profiling
#FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy -fcx-limited-range -pg -no-pie

# put this in your make.inc if you have FFTW>=3.3.5 and want thread-safe use...
#CXXFLAGS += -DFFTW_PLAN_SAFE
# FFTW base name, and math linking...
FFTWNAME = fftw3
# linux default is fftw3_omp, since 10% faster than fftw3_threads...
#FFTWOMPSUFFIX = omp
FFTWOMPSUFFIX = threads
LIBS := -lm $(CLINK)

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp
#OMP = OFF

LBLAS = -lblas -llapack

FGT_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	FGT_INSTALL_DIR = ${HOME}/lib
endif

FINUFFT_INSTALL_DIR=$(PREFIX_FINUFFT)/lib
FINUFFT_INCLUDE_DIR=$(PREFIX_FINUFFT)/include
ifeq ($(PREFIX_FINUFFT),)
	FINUFFT_INSTALL_DIR = ${HOME}/lib
	FINUFFT_INCLUDE_DIR = ${HOME}/include
endif

FFT_INSTALL_DIR=$(PREFIX_FFT)
ifeq ($(PREFIX_FFT),)
	FFT_INSTALL_DIR = /usr/lib
endif



FFT_INCLUDE_DIR=$(PREFIX_FFT_INCLUDE)
ifeq ($(PREFIX_FFT_INCLUDE),)
	FFT_INCLUDE_DIR = /usr/include
endif



# absolute path of this makefile, ie FGT's top-level directory...
FGT = $(dir $(realpath $(firstword $(MAKEFILE_LIST))))

# For your OS, override the above by placing make variables in make.inc
-include make.inc

INCL = -I$(FINUFFT_INCLUDE_DIR)
# here /usr/include needed for fftw3.f "fortran header"... (JiriK: no longer)
FFLAGS := $(FFLAGS) -I${FFT_INCLUDE_DIR} $(INCL)

DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

# single-thread total list of math and FFTW libs (now both precisions)...
# (Note: finufft tests use LIBSFFT; spread & util tests only need LIBS)
LIBSFFT := -l$(FFTWNAME) -l$(FFTWNAME)f -L$(FFT_INSTALL_DIR) -I$(FFT_INCLUDE_DIR) $(LIBS)

# multi-threaded libs & flags, and req'd flags (OO for new interface)...
ifneq ($(OMP),OFF)
  CXXFLAGS += $(OMPFLAGS)
  CFLAGS += $(OMPFLAGS)
  FFLAGS += $(OMPFLAGS)
  MFLAGS += $(MOMPFLAGS) -DR2008OO
  OFLAGS += $(OOMPFLAGS) -DR2008OO
  LIBS += $(OMPLIBS)
  ifneq ($(MINGW),ON)
    ifneq ($(MSYS),ON)
# omp override for total list of math and FFTW libs (now both precisions)...
      LIBSFFT := -l$(FFTWNAME) -l$(FFTWNAME)f -L$(FFT_INSTALL_DIR) -I$(FFT_INCLUDE_DIR) $(LIBS)
    endif
  endif
endif

LIBNAME=$(PREFIX_LIBNAME)
ifeq ($(LIBNAME),)
	LIBNAME=libfgt
endif
ifeq ($(MINGW),ON)
  DYNLIB = lib/$(LIBNAME).dll
else
  DYNLIB = lib/$(LIBNAME).so
endif

STATICLIB = lib-static/$(LIBNAME).a
# absolute path to the .so, useful for linking so executables portable...
ABSDYNLIB = $(FGT)$(DYNLIB)

FINUFFTLIBNAME = libfinufft
LFINUFFTLINKLIB = -lfinufft

LLINKLIB = $(subst lib, -l, $(LIBNAME))

FINUFFTSTATICLIB = $(FINUFFT_INSTALL_DIR)/$(FINUFFTLIBNAME)_static.a

#
# Note: the static library is used for DYLIBS, so that fmm3d
# does not get bundled in with the fmm3dbie dynamic library
#
LIBS := $(LIBSFFT)
DYLIBS := -L$(FINUFFT_INSTALL_DIR) $(LFINUFFTLINKLIB) $(LIBSFFT)
F2PYDYLIBS += -L$(FINUFFT_INSTALL_DIR) $(LFINUFFTLINKLIB) $(LIBSFFT)

LIBS += $(LBLAS) $(LDBLASINC)
DYLIBS += $(LBLAS) $(LDBLASINC)

#
# objects to compile
#
# Common objects
COM = common
COMOBJS = $(COM)/prini_new.o \
	$(COM)/hkrand.o \
	$(COM)/dlaran.o \
	$(COM)/cumsum.o \
	$(COM)/fmmcommon2d.o \
	$(COM)/pts_tree_nd.o \
	$(COM)/tree_routs_nd.o \
	$(COM)/besseljs3d.o \
	$(COM)/legeexps.o \
	$(COM)/chebexps.o \
	$(COM)/polytens.o \
	$(COM)/voltab2d.o \
	$(COM)/voltab3d.o \
	$(COM)/tree_data_routs_nd.o \
	$(COM)/tensor_prod_routs_nd.o \
	$(COM)/lapack_f77.o \
	$(COM)/tree_vol_coeffs_nd.o \
	$(COM)/fgtterms.o

# point Gauss transform objects
PFGT =  pfgt
PFGTOBJS = $(PFGT)/pfgt.o \
	$(PFGT)/pfgt_direct.o \
	$(PFGT)/pfgt_nufftrouts.o

# box Gauss transform objects
BFGT =  bfgt
BFGTOBJS = $(BFGT)/boxfgt.o \
	$(BFGT)/bfgt_volrouts.o \
	$(BFGT)/bfgt_pwrouts.o \
	$(BFGT)/bfgt_local.o


# Test objects
OBJS = $(COMOBJS) $(PFGTOBJS) $(BFGTOBJS)



.PHONY: usage lib install test-static test-dyn python

default: usage

usage:
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for fgt. Specify what to make:"
	@echo "  make install - compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR) - compile and install the main library at custom location given by PREFIX"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test-static - compile and run validation tests"
	@echo "  make test-dyn - test successful installation by validation tests linked to dynamic library"
	@echo "  make python - compile and test python interfaces using python"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=ON' for multi-threaded"
	@echo "-------------------------------------------------------------------------"
    

#
# implicit rules for objects (note -o ensures writes to correct dir)
#
%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

#
# build the library...
#
lib: $(STATICLIB) $(DYNLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNLIB) built, single-threaded versions"
endif

$(STATICLIB): $(OBJS)
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/
$(DYNLIB): $(OBJS)
	$(FC) -shared -fPIC $(OBJS) -o $(ABSDYNLIB) $(DYLIBS)
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/

install: $(STATICLIB) $(DYNLIB)
	echo $(FGT_INSTALL_DIR)
	mkdir -p $(FGT_INSTALL_DIR)
	cp -f $(DYNLIB) $(FGT_INSTALL_DIR)/
	cp -f $(STATICLIB) $(FGT_INSTALL_DIR)/
	@echo "Make sure to include " $(FGT_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(FGT_INSTALL_DIR)  " "$(LLINKLIB) " -L"$(FINUFFT_INSTALL_DIR)  " "$(LFINUFFTLINKLIB)



#inhomogeneous heat eqn solver, adaptive version
testheat1:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_forcing_adap1.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
testheat2:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_forcing_adap2.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
testheat3:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_forcing_adap3.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
	
# reaction diffusion tests
#test scant method
testsemi1:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_secant.f hpots/secant.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
# semilinear heat eqn solver, 2nd-4th order adams-moulton,adaptive
testsemi2:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_semilin_am24_init_adap1.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
testsemi3:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_semilin_am24_init_adap2.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
## self convergence --no exact solution semilinear heat eqn solver
testsemi4:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_semilin_am24_init_adap3.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2

testsemi5:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_semilin_am24_init_adap4.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2

#reaction-diffusion system, 2nd 3rd,4th,5th and 6th order adams-moulton,adaptive
# 2nd solver
testsys1:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am2.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
testsys2:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am2_adap1.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
testsys3:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am2_adap2.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
testsys4:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am2_adap3.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
testsys5:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am2_adap4.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
	
# am24
testsys6:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am24_adap1.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
testsys7:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am24_adap2.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
testsys8:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am24_adap3.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
#am6
testsys9:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am6_adap1.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
testsys10:  $(STATICLIB)
	$(FC) $(FFLAGS) hpots/test_readiff_am6_adap2.f hpots/heat_box_per.f hpots/hpots2dadap.f hpots/secant.f common/tree_ext.f common/ns_routs.f -o hpots/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd hpots; ./int2
	
	
# unsteady stokes equation solver,adaptive version

tests1:  $(STATICLIB)
	$(FC) $(FFLAGS)  unsteady_stokes_box2d_per/test_stokes_forcing1.f  unsteady_stokes_box2d_per/unsteady_stokes2d_box_per_order4.f  unsteady_stokes_box2d_per/fmm2dmk.f  unsteady_stokes_box2d_per/adaptfmm2d8.o  common/fmmrouts8.o  common/treerouts.o  unsteady_stokes_box2d_per/ns_routs.f  common/tree_ext.f  utils/tables_pois8.o  utils/chebex.o  utils/dfft.o -o   unsteady_stokes_box2d_per/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd  unsteady_stokes_box2d_per; ./int2
	
	
tests2:  $(STATICLIB)
	$(FC) $(FFLAGS)  unsteady_stokes_box2d_per/test_stokes_forcing2.f  unsteady_stokes_box2d_per/unsteady_stokes2d_box_per_order4.f  unsteady_stokes_box2d_per/fmm2dmk.f  unsteady_stokes_box2d_per/adaptfmm2d8.o  common/fmmrouts8.o  common/treerouts.o  unsteady_stokes_box2d_per/ns_routs.f  unsteady_stokes_box2d_per/tree_ext.f  utils/tables_pois8.o  utils/chebex.o  utils/dfft.o -o   unsteady_stokes_box2d_per/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd  unsteady_stokes_box2d_per; ./int2

tests3:  $(STATICLIB)
	$(FC) $(FFLAGS)  unsteady_stokes_box2d_per/test_stokes_forcing3.f  unsteady_stokes_box2d_per/unsteady_stokes2d_box_per_order4.f  unsteady_stokes_box2d_per/fmm2dmk.f  unsteady_stokes_box2d_per/adaptfmm2d8.o  common/fmmrouts8.o  common/treerouts.o  unsteady_stokes_box2d_per/ns_routs.f  common/tree_ext.f  utils/tables_pois8.o  utils/chebex.o  utils/dfft.o -o   unsteady_stokes_box2d_per/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd  unsteady_stokes_box2d_per; ./int2
	
tests4:  $(STATICLIB)
	$(FC) $(FFLAGS)  unsteady_stokes_box2d_per/test_stokes_forcing4.f  unsteady_stokes_box2d_per/unsteady_stokes2d_box_per_order4.f  unsteady_stokes_box2d_per/fmm2dmk.f  unsteady_stokes_box2d_per/adaptfmm2d8.o  common/fmmrouts8.o  common/treerouts.o  unsteady_stokes_box2d_per/ns_routs.f  common/tree_ext.f  utils/tables_pois8.o  utils/chebex.o  utils/dfft.o -o   unsteady_stokes_box2d_per/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd  unsteady_stokes_box2d_per; ./int2
#Navier-Stokes eqs solver,adaptive version,predict-correction
testnsadap1:  $(STATICLIB)
	$(FC) $(FFLAGS)  ns_solver/testnsadap1.f  ns_solver/ns2d_box_per_order4.f  ns_solver/fmm2dmk.f  ns_solver/adaptfmm2d8.o  common/fmmrouts8.o  common/treerouts.o  ns_solver/ns_routs.f  common/tree_ext.f  utils/tables_pois8.o  utils/chebex.o  utils/dfft.o -o   ns_solver/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd  ns_solver; ./int2

testnsadap2:  $(STATICLIB)
	$(FC) $(FFLAGS)  ns_solver/testnsadap2_video.f  ns_solver/ns2d_box_per_order4.f  ns_solver/fmm2dmk.f  ns_solver/adaptfmm2d8.o  common/fmmrouts8.o  common/treerouts.o  ns_solver/ns_routs.f  common/tree_ext.f  utils/tables_pois8.o  utils/chebex.o  utils/dfft.o -o   ns_solver/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd  ns_solver; ./int2
	
testnsadap3:  $(STATICLIB)
	$(FC) $(FFLAGS)  ns_solver/testnsadap3_video.f  ns_solver/ns2d_box_per_order4.f  ns_solver/fmm2dmk.f  ns_solver/adaptfmm2d8.o  common/fmmrouts8.o  common/treerouts.o  ns_solver/ns_routs.f  common/tree_ext.f  utils/tables_pois8.o  utils/chebex.o  utils/dfft.o -o   ns_solver/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd  ns_solver; ./int2
	
testnsadap4:  $(STATICLIB)
	$(FC) $(FFLAGS)  ns_solver/testnsadap4.f  ns_solver/ns2d_box_per_order4.f  ns_solver/fmm2dmk.f  ns_solver/adaptfmm2d8.o  common/fmmrouts8.o  common/treerouts.o  ns_solver/ns_routs.f  common/tree_ext.f  utils/tables_pois8.o  utils/chebex.o  utils/dfft.o -o   ns_solver/int2 $(STATICLIB) $(LFINUFFTLINKLIB) $(LIBS)
	cd  ns_solver; ./int2

	


	

#
# build the python bindings/interface
#
python: $(STATICLIB)
	cd python && export FGTND_LIBS='$(LIBS)' && pip install -e .

#
# housekeeping routines
#
clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f test/pfgt/int2-pfgt
	rm -f test/bfgt/int2-bfgt

objclean:
	rm -f $(OBJS) $(TOBJS)
	rm -f test/pfgt/*.o
	rm -f test/bfgt/*.o
