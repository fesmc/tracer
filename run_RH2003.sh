#!/usr/bin/env bash
#
# Run the Rybak & Huybrechts (2003) benchmark sweep and plot the results.
#
# Two cases, each swept over the parameters the comparison needs:
#
#   RH2003_analytic  one tracer at the divide (ux = 0), swept over dt.
#                    The divide age has the closed-form Nye-Haefeli solution
#                    (RH2003 eq. 11), so this is where accuracy is measured.
#                    Reproduces RH2003 Fig. 2.
#
#   RH2003           the full velocity-weighted profile, swept over grid
#                    resolution and interpolation method. The 501x201 run is
#                    the convergence reference for the coarse one; the spline
#                    run is the reference for the linear one, since RH2003
#                    Fig. 4a shows spline is orders of magnitude more accurate.
#                    Reproduces RH2003 Fig. 1 (streamlines, Lagrangian date).
#
# Spline is swept only on the full profile. The analytic case sets ux = 0, and
# RH2003 sec. 3.5 notes the interpolation error lives almost entirely in the
# horizontal velocity -- with no ux there is nothing for spline to improve.
# It also has nz = 2, below the 3-point minimum a b-spline needs.
#
# Precision is a compile-time choice: `prec` in src/tracer_precision.f90 is
# `sp`, and gen_filename stamps that into every output name. For a double
# precision sweep, set `prec = dp` there, `make clean && make profile`, and
# rerun -- the outputs land beside the sp ones rather than over them.
#
# Usage:
#   ./run_RH2003.sh              # build, run the sweep, plot
#   ./run_RH2003.sh --run-only   # build and run, no plots
#   ./run_RH2003.sh --plot-only  # plot whatever is already in output/

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

do_run=1
do_plot=1
case "${1:-}" in
    --run-only)  do_plot=0 ;;
    --plot-only) do_run=0 ;;
    "")          ;;
    *) echo "usage: $0 [--run-only|--plot-only]" >&2; exit 2 ;;
esac

log_dir=logs
mkdir -p "$log_dir"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Overwrite `key = value` assignments in a copy of a committed namelist.
# The stem of the generated file picks the output directory (test_profile.x
# writes to output/<stem>), so the basename is preserved and the copy lives in
# a temp dir. Any trailing `! comment` on the line is kept.
#
#   derive_nml <template.nml> <dest.nml> [key=value ...]
derive_nml() {
    local template=$1 out=$2; shift 2

    cp "$template" "$out"
    local kv key val
    for kv in "$@"; do
        key=${kv%%=*}
        val=${kv#*=}
        if ! grep -qE "^[[:space:]]*${key}[[:space:]]*=" "$out"; then
            echo "derive_nml: '$key' not found in $template" >&2
            return 1
        fi
        sed -i.bak -E "s|^([[:space:]]*${key}[[:space:]]*=[[:space:]]*)[^!]*|\1${val} |" "$out"
        rm -f "$out.bak"
    done
}

# Scale the transient parameter table.
#
# n_max_dep, dt_dep, dt_write and H_min_dep are re-read from this table on every
# timestep (tracer_par_update), so a namelist override of any of them is read at
# init and then silently discarded. Deposition density therefore has to be
# scaled here, not in the namelist.
#
# Columns: time  H_min_dep  dt_dep  n_max_dep  dt_write
#
#   scale_par_trans <table.txt> <out.txt> <n_max_dep_factor> <dt_write>
scale_par_trans() {
    local table=$1 out=$2 factor=$3 dt_write=$4

    awk -v f="$factor" -v dtw="$dt_write" '
        NR == 1                 { print; next }             # header
        /^[[:space:]]*[#!]/     { print; next }             # comments
        /^[[:space:]]*$/        { print; next }
        { printf "%12s %10s %8s %9d %10s\n", $1, $2, $3, $4 * f, dtw }
    ' "$table" > "$out"
}

# run_case <label> <template.nml> [key=value ...]
run_case() {
    local label=$1 template=$2; shift 2
    local nml="$tmp/$(basename "$template")"

    derive_nml "$template" "$nml" "$@"
    mkdir -p "output/$(basename "$template" .nml)"

    local log="$log_dir/RH2003_${label}.log"
    printf '  %-28s -> %s\n' "$label" "$log"
    if ! ./libtracer/bin/test_profile.x "$nml" > "$log" 2>&1; then
        echo "FAILED: $label (see $log)" >&2
        grep -v '^ ncio' "$log" | tail -5 >&2
        return 1
    fi
}

if (( do_run )); then
    echo "==> Building"
    make profile

    echo "==> Analytic case: one tracer at the divide, swept over dt"
    for dt in 1.0 5.0 10.0; do
        run_case "analytic_dt${dt}" RH2003_analytic.nml "dt=$dt"
    done

    echo "==> Full profile: swept over resolution and interpolation method"
    run_case profile_51x101_linear RH2003.nml nx=51 nz=101 interp_method='"linear"'
    run_case profile_51x101_spline RH2003.nml nx=51 nz=101 interp_method='"spline"'

    # The 501x201 reference needs deposition scaled with resolution. At nx=501
    # there are 375 columns eligible to deposit into (H > H_min_dep, u < U_max_dep)
    # against 38 at nx=51; leaving n_max_dep=50 would confine every deposition
    # event to the 50 slowest columns, i.e. the innermost 98 km, and the
    # "reference" would sample a different domain than the run it references.
    # Scale n_max_dep and the tracer budget n by the grid refinement factor, 10.
    # dt_write is relaxed to 5 ka because only the final slice is plotted and the
    # table's near-present 50 a cadence would write ~300 MB at 30000 tracers.
    scale_par_trans tracer_par_transient.txt "$tmp/par_trans_x10.txt" 10 5000.0
    run_case profile_501x201_linear RH2003.nml nx=501 nz=201 interp_method='"linear"' \
             n=30000 par_trans_file="\"$tmp/par_trans_x10.txt\""

    echo "==> Runs complete"
    ls -1sh output/RH2003/*.nc output/RH2003_analytic/*.nc
fi

if (( do_plot )); then
    echo "==> Plotting"
    julia --project=analysis analysis/plot_RH2003.jl
    echo "==> Figures in analysis/figures/"
fi
