#!/usr/bin/env julia
#
# Plot the RH2003 benchmark sweep produced by run_RH2003.sh.
#
# Three figures land in analysis/figures/:
#
#   fig1_streamlines_age.png  streamlines and Lagrangian date of deposition,
#                             against RH2003 Fig. 1 (panels A-1 and A-2)
#   fig2_divide_age.png       age at the divide vs the Nye-Haefeli analytic
#                             solution, against RH2003 Fig. 2 (Lagrangian curve)
#   fig3_age_error.png        sensitivity of the age field to grid resolution
#                             and to linear vs spline interpolation
#
# Run via ./run_RH2003.sh, or directly:  julia --project=analysis analysis/plot_RH2003.jl

using CairoMakie
using Printf
using Statistics: mean

include(joinpath(@__DIR__, "rh2003.jl"))
using .RH2003

CairoMakie.activate!(type="png", px_per_unit=2)

const OUT = joinpath(@__DIR__, "figures")
const PROF = joinpath(@__DIR__, "..", "output", "RH2003")
const ANLY = joinpath(@__DIR__, "..", "output", "RH2003_analytic")

mkpath(OUT)

# Common grid for putting scattered tracer ages onto a field. 40 x 20 cells over
# 3000 tracers averages ~5 tracers per occupied cell; finer than that and the
# maps show the scatter between individual particles rather than the age field.
const XEDGES = range(0, 1000, length=41)
const SEDGES = range(0, 1, length=21)

midpoints(e) = [(e[i] + e[i+1]) / 2 for i in 1:length(e)-1]

prof_file(nx, nz, interp, dt) =
    joinpath(PROF, @sprintf("RH2003_%d_%d_0.10_0.00_sp_%s_%.1f.nc", nx, nz, interp, dt))
anly_file(dt) =
    joinpath(ANLY, @sprintf("RH2003_51_2_0.10_0.00_sp_linear_%.1f.nc", dt))

function require(path)
    isfile(path) || error("missing output: $path\nRun ./run_RH2003.sh --run-only first.")
    return path
end

# ---------------------------------------------------------------------------
# Consistency check: the Fortran profile against the paper's equations.
# ---------------------------------------------------------------------------

let p = load_profile(require(joinpath(ANLY, "profile_RH2003.nc")))
    ref = vialov_profile(nx=51, nz=101)
    @printf("H0 (model)     = %.1f m\n", p.H0)
    @printf("H0 (eq. 7)     = %.1f m\n", ref.H0)
    @printf("H0 (RH2003 p.152) = 3598.4 m\n")
    abs(p.H0 - ref.H0) < 1.0 ||
        error("model H0 disagrees with eq. 7 by $(abs(p.H0 - ref.H0)) m")
end

# ---------------------------------------------------------------------------
# Figure 1 -- streamlines and the Lagrangian date field (cf. RH2003 Fig. 1)
# ---------------------------------------------------------------------------

function fig_streamlines_age()
    trc = load_tracer(require(prof_file(51, 101, "linear", 10.0)))
    prof = load_profile(require(joinpath(PROF, "profile_RH2003.nc")))

    fig = Figure(size=(760, 620))

    # -- streamlines: in a stationary field, particle paths are streamlines
    ax1 = Axis(fig[1, 1], ylabel="Height (km)", title="Streamlines",
               limits=(0, 1000, 0, 4))
    # The tracers deposited in the very first event, i.e. one per eligible column
    ids0 = sort(unique(filter(!isnan, trc.id[:, 1])))
    for (xs, zs) in trajectories(trc, Int.(ids0))
        lines!(ax1, xs, zs ./ 1e3, color=(:steelblue, 0.55), linewidth=0.8)
    end
    lines!(ax1, prof.xc ./ 1e3, prof.H ./ 1e3, color=:black, linewidth=2.5)
    hidexdecorations!(ax1, grid=false)

    # -- Lagrangian date of deposition, contoured in dimensionless height
    x, sigma, age = final_slice(trc)
    grid = bin_mean(x, sigma, age, XEDGES, SEDGES)
    xc, sc = midpoints(XEDGES), midpoints(SEDGES)

    ax2 = Axis(fig[2, 1], xlabel="Distance from the divide (km)",
               ylabel="Dimensionless height", title="Lagrangian date (kyr)",
               limits=(0, 1000, 0, 1))
    hm = heatmap!(ax2, xc, sc, grid, colormap=:viridis, colorrange=(0, 120))
    # The paper's contour levels (its dates are negative, i.e. before present)
    levels = [2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 100.0]
    contour!(ax2, xc, sc, grid; levels, color=:white, linewidth=1.0,
             labels=true, labelsize=9, labelcolor=:white)
    Colorbar(fig[2, 2], hm, label="Age (kyr)")

    rowsize!(fig.layout, 1, Relative(0.42))
    save(joinpath(OUT, "fig1_streamlines_age.png"), fig)
    println("wrote fig1_streamlines_age.png  ($(length(ids0)) streamlines, $(length(age)) tracers)")
end

# ---------------------------------------------------------------------------
# Figure 2 -- age at the divide vs the analytic solution (cf. RH2003 Fig. 2)
# ---------------------------------------------------------------------------

function fig_divide_age()
    prof = load_profile(require(joinpath(ANLY, "profile_RH2003.nc")))
    dts = [1.0, 5.0, 10.0]
    runs = [(dt, load_tracer(require(anly_file(dt)))) for dt in dts]

    fig = Figure(size=(820, 480))

    ax1 = Axis(fig[1, 1], xlabel="Age (ka)", ylabel="Dimensionless height at the divide",
               limits=(-165, 5, -0.01, 1.02))
    ax2 = Axis(fig[1, 2], xlabel="Relative error on Lagrangian date (%)",
               xscale=log10, limits=(1e-6, 1e1, -0.01, 1.02))
    hideydecorations!(ax2, grid=false)

    # The analytic curve, drawn on the model's own sigma axis
    sig = range(1e-4, 1.0, length=400)
    lines!(ax1, .-analytic_age(sig, prof.H0, prof.G) ./ 1e3, sig,
           color=(:grey30, 0.9), linewidth=6, label="Analytical solution (eq. 11)")

    colors = [:magenta, :darkorange, :seagreen]
    for ((dt, trc), col) in zip(runs, colors)
        # One tracer, tracked over time: each output step gives a (sigma, age) pair
        sigma = vec(trc.sigma[1, :])
        age = vec(trc.age[1, :])          # ka
        keep = @. !isnan(sigma) & !isnan(age)
        sigma, age = sigma[keep], age[keep]

        lines!(ax1, -age, sigma, color=col, linewidth=1.4,
               label=@sprintf("Lagrangian tracer, dt = %g a", dt))

        an = analytic_age(sigma, prof.H0, prof.G) ./ 1e3     # ka
        err = @. 100 * (age - an) / an
        # Near the surface the age passes through zero and the relative error
        # is meaningless; the R original masked |age| < 0.1 ka the same way.
        vis = @. !isnan(err) & (abs(age) > 0.1) & (abs(err) > 0)
        lines!(ax2, abs.(err[vis]), sigma[vis], color=col, linewidth=1.4,
               label=@sprintf("dt = %g a", dt))

        finite = err[@. !isnan(err) & (abs(age) > 0.1)]
        @printf("dt = %5.1f a:  max |err| = %.4f %%   median |err| = %.5f %%\n",
                dt, maximum(abs, finite), median_abs(finite))
    end

    # RH2003 reports the Lagrangian error staying below 0.1% along the vertical
    vlines!(ax2, [0.1], color=:red, linestyle=:dash, linewidth=1)
    text!(ax2, 0.08, 0.03, text="RH2003: Lagrangian < 0.1%", color=:red,
          fontsize=9, align=(:right, :bottom))

    axislegend(ax1, position=:lt, framevisible=false, labelsize=9)
    axislegend(ax2, position=:lt, framevisible=false, labelsize=9)

    save(joinpath(OUT, "fig2_divide_age.png"), fig)
    println("wrote fig2_divide_age.png")
end

median_abs(v) = (s = sort(abs.(v)); isempty(s) ? NaN : s[cld(length(s), 2)])

# ---------------------------------------------------------------------------
# Figure 3 -- age error from grid resolution and interpolation method
# ---------------------------------------------------------------------------

function fig_age_error()
    coarse = load_tracer(require(prof_file(51, 101, "linear", 10.0)))
    spline = load_tracer(require(prof_file(51, 101, "spline", 10.0)))
    fine = load_tracer(require(prof_file(501, 201, "linear", 10.0)))

    xc, sc = midpoints(XEDGES), midpoints(SEDGES)
    # Differencing two runs cell by cell: demand a few tracers in each, or the
    # map shows particle scatter instead of the difference between age fields.
    grid(trc; mincount) = (s = final_slice(trc);
                           bin_mean(s[1], s[2], s[3], XEDGES, SEDGES; mincount))

    g_coarse = grid(coarse, mincount=3)
    g_spline = grid(spline, mincount=3)
    g_fine = grid(fine, mincount=3)

    fig = Figure(size=(1180, 380))

    ax1 = Axis(fig[1, 1], xlabel="Distance from the divide (km)",
               ylabel="Dimensionless height", title="Age, 51x101 linear (ka)")
    hm1 = heatmap!(ax1, xc, sc, grid(coarse, mincount=1),
                   colormap=:viridis, colorrange=(0, 120))
    Colorbar(fig[1, 2], hm1)

    # vs the high-resolution reference
    d_res = g_coarse .- g_fine
    ax2 = Axis(fig[1, 3], xlabel="Distance from the divide (km)",
               title="Age error vs 501x201 (ka)")
    hm2 = heatmap!(ax2, xc, sc, d_res, colormap=:balance, colorrange=(-5, 5))
    Colorbar(fig[1, 4], hm2)
    hideydecorations!(ax2, grid=false)

    # vs the spline reference at the same resolution (cf. RH2003 Fig. 4)
    d_int = g_coarse .- g_spline
    ax3 = Axis(fig[1, 5], xlabel="Distance from the divide (km)",
               title="Linear - spline (ka)")
    hm3 = heatmap!(ax3, xc, sc, d_int, colormap=:balance, colorrange=(-5, 5))
    Colorbar(fig[1, 6], hm3)
    hideydecorations!(ax3, grid=false)

    save(joinpath(OUT, "fig3_age_error.png"), fig)

    stat(d, lab) = begin
        v = filter(!isnan, vec(d))
        @printf("%-28s  mean = %+7.3f ka  rms = %6.3f ka  max|.| = %6.3f ka  (n=%d)\n",
                lab, mean(v), sqrt(mean(v .^ 2)), maximum(abs, v), length(v))
    end
    stat(d_res, "51x101 linear - 501x201")
    stat(d_int, "51x101 linear - spline")
    println("wrote fig3_age_error.png")
end

fig_streamlines_age()
fig_divide_age()
fig_age_error()
println("\nFigures written to $(normpath(OUT))")
