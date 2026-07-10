# Tracer

A module to insert Lagrangian tracer particles in a modeling domain
and follow them over the course of a simulation.

Tracers are deposited at the ice surface, advected by a 3D velocity field, and
carry the time, location and (optionally) the climate conditions of their
deposition. `tracer` builds as a static library, `libtracer.a`, and ships two
stand-alone benchmarks.

## Install

`tracer` is a [configme](https://github.com/fesmc/configme) package. It depends
on [fesm-utils](https://github.com/fesmc/fesm-utils) for its `ncio` and `nml`
modules, and vendors `bspline-fortran` for spline interpolation.

```bash
configme install tracer --only     # tracer + fesm-utils
```

It is also nested inside a `yelmo` checkout (`yelmo/tracer`), the same way
`FastHydrology` is, so `configme install yelmo` brings it along.

To configure an existing checkout for your machine:

```bash
cd tracer
configme -m macbook -c gfortran    # writes the repo-root Makefile
```

## Build

```bash
make tracer-static     # libtracer/include/libtracer.a
make all               # both benchmark programs
make usage             # all targets
```

Add `debug=1` to any build for bounds checking and floating-point traps.

## Benchmarks

All three are bit-reproducible: deposition draws on the random number generator
(for the location jitter and the `rand` weight), which is seeded from the `seed`
namelist parameter. Set `seed <= 0` for a nondeterministic run.

```bash
make check                    # runs all three cases below; nonzero exit on any failure

make run-profile              # RH2003, velocity-weighted deposition
make run-profile-analytic     # RH2003, one tracer at the divide, vs. the analytic age
make run-greenland            # 3D Greenland, forced from a 16 km Yelmo restart
```

**RH2003** is the idealized ice-divide profile of Rybak & Huybrechts (2003)
(`docs/rybak_and_huybrechts_2003.pdf`). The analytic case deposits a single
tracer at the divide, where the age profile has a closed-form solution — this is
where quantitative accuracy is checked.

To reproduce the paper's figures, run the sweep and plot it:

```bash
./run_RH2003.sh               # build, run the sweep, write analysis/figures/
./run_RH2003.sh --run-only    # runs only
./run_RH2003.sh --plot-only   # plots only, from existing output/
```

This sweeps `dt` at the divide, and grid resolution and `interp_method` on the
full profile, then writes three figures with `analysis/plot_RH2003.jl`
(Julia; `julia --project=analysis`):

| figure | compares against |
|---|---|
| `fig1_streamlines_age.png` | RH2003 Fig. 1: streamlines and the Lagrangian date field |
| `fig2_divide_age.png` | RH2003 Fig. 2 and eq. 11: age at the divide vs. the Nye–Haefeli solution |
| `fig3_age_error.png` | sensitivity of the age field to resolution and to linear vs. spline |

Precision is a compile-time choice: `prec` in `src/tracer_precision.f90` is `sp`,
and it is stamped into every output filename. For a double-precision sweep set
`prec = dp`, `make clean && make profile`, and rerun — the `dp` outputs land
beside the `sp` ones rather than over them.

**Greenland** advects tracers through a Yelmo ice-sheet state on the 16 km
grid, read from `data/initmip-grl-16km/yelmo_restart.nc`. The restart carries no
spun-up age field, so this case is a structural check — particles stay inside
the ice column, ages stay bounded by the run length — rather than a comparison
against a known answer.

The forcing is read at an explicit index along the restart's time dimension, so
the same driver serves an offline transient run forced by a sequence of
ice-sheet states. Inside a coupled model the fields are passed directly instead.

## Using the library

```fortran
use tracer

type(tracer_class) :: trc

call tracer_init(trc,"Greenland.nml",time=time,x=xc,y=yc,is_sigma=.TRUE.)

call tracer_update(trc,time=time,x=xc,y=yc,z=zeta,z_srf=z_srf,H=H_ice, &
                   ux=ux,uy=uy,uz=uz,dep_now=dep_now,stats_now=stats_now)
```

`ux`, `uy` and `uz` must be cell-centred on a single ascending sigma axis with
`sigma(nz) = 1` at the surface, and `uz` positive upward. A model that staggers
its velocities (Yelmo does) must average them onto cell centres first; see
`src/test_greenland.f90`.

The deposition tagging fields — `lon`, `lat`, `t2m_ann`, `t2m_sum`, `pr_ann`,
`pr_sum`, `d18O_ann` — are all optional. Any field omitted is recorded as
missing, so advecting particles requires no climate forcing.

## Parameters

Set in the `&tracer_par` namelist group:

| parameter | meaning |
|---|---|
| `n`, `n_max_dep` | tracer slots; max deposited per deposition step |
| `dt`, `dt_dep`, `dt_write` | timestep; deposition and write intervals [a] |
| `H_min`, `depth_max`, `U_max` | bounds beyond which a tracer is deactivated |
| `H_min_dep`, `U_max_dep` | ice must be thicker / slower than this to deposit |
| `weight` | deposition priority: `vel`, `linear`, `quadratic`, `rand` |
| `alpha` | slope of the priority function (`linear`, `quadratic` only) |
| `noise` | jitter the deposition location within its grid cell |
| `seed` | RNG seed; positive is reproducible, `<= 0` defers to the OS |
| `interp_method` | `linear` or `spline` |
| `par_trans_file` | transient parameter table, or `"None"` |

`weight` ranks cells rather than being sampled from: each deposition step takes
the `n_max_dep` highest-ranked cells outright, so `vel` deposits into the
slowest ice available and only reaches faster ice once every eligible cell holds
a tracer. (`rand` is the exception — taking the maximum of a random field is
equivalent to sampling.) Cells are eligible only where `H > H_min_dep` and the
surface speed is below `U_max_dep`.

`H_min_dep`, `dt_dep`, `dt_write` and `n_max_dep` are re-read from
`par_trans_file` on every timestep when one is given, which silently overrides
whatever the namelist set them to.

A parameter named in the namelist but not read by the code is ignored; a
parameter the code reads but the namelist omits is a hard error.

## Roadmap

Marked `== TO DO ==` in the source:

- Tracer thickness and the destruction of over-thinned tracers
  (`tracer3D.f90`), which would restore a `thk_min` parameter.
- Capping the density of tracers near the surface at deposition
  (`tracer3D.f90`), which would restore `dens_z_lim` / `dens_max`.
- Precipitation-weighted deposition temperature, `t2m_prann`.

These parameters were removed rather than left as namelist knobs that silently
do nothing.
