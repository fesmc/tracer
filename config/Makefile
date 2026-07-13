.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = libtracer/include
bindir = libtracer/bin
libdir = libs

# Command-line options at make call
debug  ?= 0
openmp ?= 0

## COMPILER CONFIGURATION ##
# Assembled by `configme`. The fragment defines the compiler (FC, FFLAGS,
# FFLAGS_OPENMP, DFLAGS_*), netCDF wiring (INC_NC, LIB_NC), and any
# host-specific overrides (e.g. LFLAGS_EXTRA).

<COMPILER_CONFIGURATION>

## SHARED DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS,
# FFLAGS_OPENMP, LIB_NC and LFLAGS_EXTRA set there. Defines the FESMUTILS
# paths, the openmp= build switch, and the final LFLAGS.

include config/common.mk

# Determine whether to use normal flags or debugging flags
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

###############################################
##
## List of tracer rules and source files
##
###############################################

## bspline-fortran ##
# Vendored: fesm-utils does not provide a b-spline interpolator.

$(objdir)/bspline_sub_module.o: $(libdir)/bspline-fortran/bspline_sub_module.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/bspline_oo_module.o: $(libdir)/bspline-fortran/bspline_oo_module.f90 $(objdir)/bspline_sub_module.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/bspline_module.o: $(libdir)/bspline-fortran/bspline_module.f90 $(objdir)/bspline_sub_module.o $(objdir)/bspline_oo_module.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

obj_bspline =	$(objdir)/bspline_sub_module.o \
				$(objdir)/bspline_oo_module.o \
				$(objdir)/bspline_module.o

## tracer ##

$(objdir)/tracer_precision.o: $(srcdir)/tracer_precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/tracer_constants.o: $(srcdir)/tracer_constants.f90 $(objdir)/tracer_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/tracer_interp.o: $(srcdir)/tracer_interp.f90 $(objdir)/tracer_precision.o $(objdir)/tracer_constants.o $(objdir)/bspline_module.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/tracer_stats.o: $(srcdir)/tracer_stats.f90 $(objdir)/tracer_precision.o $(objdir)/tracer_constants.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/tracer3D.o: $(srcdir)/tracer3D.f90 $(objdir)/tracer_precision.o $(objdir)/tracer_constants.o $(objdir)/tracer_interp.o \
					$(objdir)/tracer_stats.o $(objdir)/bspline_module.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/tracer2D.o: $(srcdir)/tracer2D.f90 $(objdir)/tracer3D.o $(objdir)/tracer_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/tracer_io.o: $(srcdir)/tracer_io.f90 $(objdir)/tracer3D.o $(objdir)/tracer_precision.o $(objdir)/tracer_constants.o $(objdir)/tracer_interp.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/tracer.o: $(srcdir)/tracer.f90 $(objdir)/tracer_precision.o $(objdir)/tracer_constants.o $(objdir)/tracer_interp.o \
					$(objdir)/tracer_stats.o $(objdir)/tracer3D.o $(objdir)/tracer2D.o $(objdir)/tracer_io.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

obj_tracer =	$(objdir)/tracer_precision.o \
				$(objdir)/tracer_constants.o \
				$(objdir)/tracer_interp.o \
				$(objdir)/tracer_stats.o \
				$(objdir)/tracer3D.o \
				$(objdir)/tracer2D.o \
				$(objdir)/tracer_io.o \
				$(objdir)/tracer.o

###############################################
##
## Compilation of complete programs
##
###############################################

# Static library compilation
tracer-static: $(obj_bspline) $(obj_tracer)
	ar rc $(objdir)/libtracer.a $(obj_bspline) $(obj_tracer)
	ranlib $(objdir)/libtracer.a
	@echo " "
	@echo "    $(objdir)/libtracer.a is ready."
	@echo " "

# Benchmark: Rybak & Huybrechts (2003) idealized ice-divide profile.
# Run with RH2003.nml (velocity-weighted deposition) or RH2003_analytic.nml
# (thickness-weighted, single tracer at the divide, compared to the analytic age).
profile: tracer-static
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_profile.x $(srcdir)/test_profile.f90 \
		$(LIB_TRACER) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_profile.x is ready."
	@echo " "

# Benchmark: 3D Greenland, forced offline from a Yelmo restart at 16 km.
greenland: tracer-static
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_greenland.x $(srcdir)/test_greenland.f90 \
		$(LIB_TRACER) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_greenland.x is ready."
	@echo " "

.PHONY: all usage clean check run-profile run-profile-analytic run-greenland
all: profile greenland

###############################################
##
## Running the benchmarks
##
## Each rule creates the output directory the program writes into; netCDF will
## not create it, so a fresh clone cannot run the binary directly.
##
###############################################

run-profile: profile
	@mkdir -p output/RH2003
	./$(bindir)/test_profile.x RH2003.nml

run-profile-analytic: profile
	@mkdir -p output/RH2003_analytic
	./$(bindir)/test_profile.x RH2003_analytic.nml

run-greenland: greenland
	@mkdir -p output/GRL-16KM
	./$(bindir)/test_greenland.x

# Run every benchmark. Any of them failing (a missing namelist parameter, an
# unreadable forcing file, zero tracers deposited) exits nonzero and fails here.
check: run-profile-analytic run-profile run-greenland
	@echo " "
	@echo "    All benchmarks passed."
	@echo " "

usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make tracer-static         : compiles the static library libtracer.a"
	@echo " make profile               : compiles test_profile.x, the RH2003 idealized profile benchmark"
	@echo " make greenland             : compiles test_greenland.x, the 3D Greenland benchmark"
	@echo " make all                   : compiles both benchmark programs"
	@echo ""
	@echo " make run-profile           : runs the RH2003 benchmark (velocity-weighted deposition)"
	@echo " make run-profile-analytic  : runs the RH2003 analytic case (one tracer at the divide)"
	@echo " make run-greenland         : runs the 3D Greenland benchmark"
	@echo " make check                 : runs all benchmarks; nonzero exit on any failure"
	@echo ""
	@echo " make clean                 : cleans object files and programs"
	@echo " Add debug=1 to any build for bounds checking and floating-point traps."
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM $(bindir)/*.x.dSYM
