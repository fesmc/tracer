# Shared build configuration for tracer (dependency wiring).
#
# Loaded after the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). References FFLAGS /
# FFLAGS_OPENMP (compiler) and LIB_NC (machine or auto-detected netCDF).

# fesm-utils provides the ncio and nml modules that tracer uses. tracer carries
# no netCDF calls of its own, so INC_NC is not needed at compile time; only
# LIB_NC at link time, via ncio inside libfesmutils.
FESMUTILSROOT = fesm-utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

# tracer needs the fesm-utils includes in its compile flags.
FFLAGS += $(INC_FESMUTILS)

ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    FFLAGS += $(FFLAGS_OPENMP)
endif

# tracer build outputs, for downstream linking by a consumer (e.g. yelmo, which
# nests this checkout at yelmo/tracer).
TRACERROOT = ${CURDIR}
INC_TRACER = -I${TRACERROOT}/libtracer/include
LIB_TRACER = -L${TRACERROOT}/libtracer/include -ltracer

# Extra link flags. -Wl,-zmuldefs works around duplicate symbols in the static
# deps (the default on Linux). A machine fragment can disable it by setting
# `LFLAGS_EXTRA =` (macOS ld rejects -zmuldefs, so the macbook fragment does).
LFLAGS_EXTRA ?= -Wl,-zmuldefs

LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LFLAGS_EXTRA)
