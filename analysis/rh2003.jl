"""
    RH2003

Analytic solutions and output loaders for the Rybak & Huybrechts (2003)
idealized ice-divide benchmark.

Reference:
  Rybak, O. and Huybrechts, P. (2003). A comparison of Eulerian and Lagrangian
  methods for dating in numerical ice-sheet models. Ann. Glaciol. 37, 150-158.
  (docs/rybak_and_huybrechts_2003.pdf)

Equation numbers below refer to that paper.
"""
module RH2003

using NCDatasets
using Statistics: mean

export Profile, vialov_profile, analytic_age, load_profile, load_tracer,
       final_slice, trajectories, bin_mean

"The prescribed Nye-Vialov ice sheet: geometry and stationary velocity field."
struct Profile
    xc::Vector{Float64}      # distance from divide [m]
    sigma::Vector{Float64}   # dimensionless height, 0 at bed, 1 at surface
    H0::Float64              # thickness at the divide [m]
    G::Float64               # surface accumulation [m/a]
    B::Float64               # basal melt [m/a]
    L::Float64               # half-width [m]
    H::Vector{Float64}       # thickness [m]
    ux::Matrix{Float64}      # horizontal velocity [m/a], (x, sigma)
    uz::Matrix{Float64}      # vertical velocity [m/a], (x, sigma)
end

"""
    vialov_profile(; nx=51, nz=101, L=1e6, G=0.10, B=0.0, A=1e-16, ng=3,
                     rho=910.0, g=9.81, ux_fac=1.0)

The standard experiment A of RH2003: flat bed, no sliding, no basal melting.

Surface profile is the Nye-Vialov solution (eq. 6-7), the horizontal velocity
follows from it along the flowline (eq. 8), and the vertical velocity uses the
linear form of eq. 9 -- which is inconsistent with eq. 8, but is what makes the
closed-form divide age of eq. 11 exact.

This mirrors `calc_profile_RH2003` in src/test_profile.f90; agreement between
the two is checked in plot_RH2003.jl.
"""
function vialov_profile(; nx=51, nz=101, L=1e6, G=0.10, B=0.0, A=1e-16, ng=3,
                          rho=910.0, g=9.81, ux_fac=1.0)
    M = G - B

    xc = collect(range(0.0, 1000.0e3, length=nx))
    sigma = collect(range(0.0, 1.0, length=nz))
    nz > 1 && (sigma[1] = 1e-8)          # avoid the log singularity at the bed

    # eq. 7 -- with the standard parameters this gives H0 = 3598.4 m
    H0 = (20M / A)^(1 / (2 * (ng + 1))) * (1 / (rho * g))^(ng / (2 * (ng + 1))) * sqrt(L)

    # eq. 6
    H = @. H0 * (1 - (abs(xc) / L)^((ng + 1) / ng))^(ng / (2 * (ng + 1)))

    # Centred differences, one-sided ends left at zero (as in test_profile.f90)
    dHdx = zeros(nx)
    for i in 2:nx-1
        dHdx[i] = (H[i+1] - H[i-1]) / (xc[i+1] - xc[i-1])
    end

    ux = zeros(nx, nz)
    uz = zeros(nx, nz)
    for i in 1:nx, j in 1:nz
        # eq. 8
        ux[i, j] = -(2A) / (ng + 1) * (rho * g)^ng * abs(dHdx[i])^(ng - 1) * dHdx[i] *
                   (H[i]^(ng + 1) - (H[i] - sigma[j] * H[i])^(ng + 1))
        ux[i, j] *= ux_fac
        # eq. 9
        uz[i, j] = sigma[j] * (-G + B + ux[i, j] * dHdx[i]) - B
    end

    return Profile(xc, sigma, H0, G, B, L, H, ux, uz)
end

"""
    analytic_age(sigma, H0, G)

Nye-Haefeli date of deposition at the ice divide (eq. 11), in years before
present, returned as a positive age. `sigma` is dimensionless height.
Non-finite results (sigma <= 0) come back as `NaN`.
"""
function analytic_age(sigma, H0, G)
    age = @. -(H0 / G) * log(sigma)
    return @. ifelse(isfinite(age), age, NaN)
end

analytic_age(p::Profile) = analytic_age(p.sigma, p.H0, p.G)

_nan(x) = ismissing(x) ? NaN : Float64(x)
_nanmat(a) = _nan.(a)

"""
    load_profile(path)

Read the prescribed profile that `test_profile.x` wrote (`profile_RH2003.nc`).
This is the velocity field the tracers actually saw, so plots of streamlines and
of the analytic age should use it rather than recomputing from `vialov_profile`.
"""
function load_profile(path)
    NCDataset(path) do ds
        xc = _nanmat(ds["xc"][:]) .* 1e3      # file stores km
        sigma = _nanmat(ds["sigma"][:])
        Profile(xc, sigma,
                _nan(ds["H0"][1]), _nan(ds["G"][1]), _nan(ds["B"][1]), _nan(ds["L"][1]),
                _nanmat(ds["H"][:]), _nanmat(ds["ux"][:, :]), _nanmat(ds["uz"][:, :]))
    end
end

"""
    load_tracer(path)

Read a tracer output file. Fields come back as `(pt, time)` matrices with the
missing value replaced by `NaN`; `time`, `age` and `dep_time` are converted from
years to ka. `sigma` is the dimensionless height `z / H`.
"""
function load_tracer(path)
    NCDataset(path) do ds
        get2d(v) = _nanmat(ds[v][:, :])
        x = get2d("x")            # km
        z = get2d("z")            # m
        H = get2d("H")            # m
        (; path,
           time = _nanmat(ds["time"][:]) ./ 1e3,
           n_active = ds["n_active"][:],
           id = get2d("id"),
           x, z, H,
           sigma = z ./ H,
           ux = get2d("ux"),
           uz = get2d("uz"),
           age = get2d("age") ./ 1e3,
           dep_time = get2d("dep_time") ./ 1e3,
           dep_x = get2d("dep_x"))
    end
end

"""
    final_slice(trc)

The active tracers at the last output time, as `(x, sigma, age)` vectors.
A slot with no tracer in it carries `NaN` and is dropped.
"""
function final_slice(trc)
    x = trc.x[:, end]
    sigma = trc.sigma[:, end]
    age = trc.age[:, end]
    keep = @. !isnan(x) & !isnan(sigma) & !isnan(age)
    return x[keep], sigma[keep], age[keep]
end

"""
    trajectories(trc, ids)

Path of each tracer id through `(x [km], z [m])`, ordered in time. In a
stationary velocity field these coincide with streamlines (RH2003 Fig. 1a).
"""
function trajectories(trc, ids)
    want = Set(ids)
    # Time is the outer loop so each path accumulates already ordered.
    paths = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
    for t in axes(trc.id, 2), p in axes(trc.id, 1)
        i = trc.id[p, t]
        (isnan(i) || i ∉ want || isnan(trc.x[p, t])) && continue
        xs, zs = get!(paths, Int(i)) do
            (Float64[], Float64[])
        end
        push!(xs, trc.x[p, t])
        push!(zs, trc.z[p, t])
    end
    return [paths[i] for i in sort(collect(keys(paths))) if length(paths[i][1]) > 1]
end

"""
    bin_mean(x, y, z, xedges, yedges; mincount=1)

Mean of `z` over each cell of the `xedges` x `yedges` grid; cells holding fewer
than `mincount` samples are `NaN`. Used to put scattered tracer ages onto a
common grid so that two runs can be differenced.

Raise `mincount` when differencing two runs: a cell holding one tracer from each
reports the scatter between two individual particles, not the difference between
the two age fields.
"""
function bin_mean(x, y, z, xedges, yedges; mincount::Int=1)
    nx, ny = length(xedges) - 1, length(yedges) - 1
    tot = zeros(nx, ny)
    cnt = zeros(Int, nx, ny)
    for k in eachindex(x)
        (isnan(x[k]) || isnan(y[k]) || isnan(z[k])) && continue
        i = searchsortedlast(xedges, x[k])
        j = searchsortedlast(yedges, y[k])
        i == nx + 1 && (i = nx)      # include the right edge
        j == ny + 1 && (j = ny)
        (1 <= i <= nx && 1 <= j <= ny) || continue
        tot[i, j] += z[k]
        cnt[i, j] += 1
    end
    out = fill(NaN, nx, ny)
    ok = cnt .>= mincount
    @. out[ok] = tot[ok] / cnt[ok]
    return out
end

end # module
