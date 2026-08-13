### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Ray Theory and The Eikonal Equation"
#> tags = ["raytheory"]
#> layout = "layout.jlhtml"
#> description = "Paint velocity perturbations onto a background Earth model and watch how seismic rays bend, and how the first-arrival wavefront sweeps through the heterogeneity."

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 05d2c8b6-41a3-445b-8988-34b5020087d5
begin
    using LinearAlgebra
    using Interpolations
    using Eikonal
    using PlutoUI
end

# ╔═╡ ef66e643-bc46-4a44-bdb9-d641a26f4e0e
TableOfContents()

# ╔═╡ 864f274b-4f28-4cf9-83f7-69a1f03ca870
md"""
# Ray Theory and The Eikonal Equation

Seismic ray theory is the high-frequency limit of wave propagation: when a wave's wavelength is
much shorter than the scale over which the medium changes, energy travels along well-defined
paths — rays — that bend smoothly through velocity heterogeneity exactly the way light bends
through a lens. This notebook derives the Eikonal equation that governs those rays from the
elastic wave equation, then hands you a widget to see it in action: **paint velocity
perturbations directly onto a background Earth model** and watch the rays bend and the
first-arrival wavefront sweep through whatever structure you've drawn.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 28f9d159-170e-4371-85c4-98566e90d331
md"""
## Deriving the Eikonal Equation

We look for a high-frequency (WKB) solution of the scalar wave equation
`` \rho\,\partial_t^2\phi = \nabla\cdot(\mu\nabla\phi) `` of the form

```math
\phi(x,z,t) = A(x,z)\,\exp\!\big[i\omega\big(t - T(x,z)\big)\big],
```

an amplitude `` A `` that varies slowly in space compared to the rapidly-oscillating phase
`` \omega T ``. Substituting `` \phi `` into the wave equation and sorting the resulting terms by
their power of `` \omega `` gives, at leading order `` O(\omega^2) `` (the **high-frequency
approximation**: every term of lower order in `` \omega `` is discarded, exactly the terms that
matter close to the source — see "Near Field vs. Far Field" in `lame-theorem.jl` for the point-
force version of this same tradeoff), the **Eikonal equation**

```math
\big(\partial_x T\big)^2 + \big(\partial_z T\big)^2 = \frac{1}{\beta(x,z)^2}, \qquad (\text{eik})
```

and at the next order, `` O(\omega^1) ``, the **transport equation** that governs how the
amplitude `` A `` evolves:

```math
2\,\nabla A\cdot\nabla T + A\,\nabla^2 T = 0. \qquad (\text{eikA})
```

!!! tip "Dispersion relation, revisited"
	The Eikonal equation looks just like the dispersion relation from a plane-wave solution of
	the homogeneous scalar wave equation: the magnitude of the wavenumber vector `` \vec k ``
	equals `` \omega^2/\alpha^2 ``, which also constrains the magnitude of the slowness vector
	`` \vec s `` to `` 1/\alpha^2 ``.

!!! tip "A local plane wave"
	When solving the Eikonal equation, think of a *local plane wave* at `` (x,z) `` with its own
	slowness vector `` \vec s(x,z) ``. `` \partial_x T `` is exactly the x-component of that
	slowness vector.
"""

# ╔═╡ 0f652df5-37a3-40fe-a8c9-5d873d7b491a
md"""
## Wavefronts and Rays

The surfaces `` T(x,y,z) = \text{const.} `` are **wavefronts**. In a homogeneous medium they are
spheres; in a heterogeneous one they take whatever shape the velocity structure forces on them —
exactly what the widget's animated sweep shows.

**Rays** are lines perpendicular to the wavefronts, i.e. parallel to `` \nabla T(x,y,z) ``. If
`` \hat s `` is the unit vector along `` \nabla T ``, then

```math
\nabla T(x,y,z) = \vec{s}(x,z) = |\vec{s}(x,z)|\,\hat{s}(x,z), \qquad (1)
```

where `` |\vec s| = 1/\alpha `` is the local slowness.
"""

# ╔═╡ 8dd31a81-a418-4164-90aa-b65e48ed9b2c
md"""
## Tracing the Ray Path

A ray begins at the source, where we choose the direction of the outgoing wavefront
`` \hat s_0 ``. Denoting the source position `` \vec p_0 ``, we move along `` \hat s_0 `` by an
incremental length `` dl ``:

```math
\vec{p}_1 = \vec{p}_0 + \hat{s}_0 \,dl. \qquad(2)
```

To continue the ray, we need the change in `` \vec s `` between `` \vec p_0 `` and `` \vec p_1 ``
— the plane wave we're riding gets refracted as it propagates through changing slowness
`` s=|\vec s| ``. Differentiating the Eikonal equation `` |\vec s|^2=s^2 `` along the ray's own
arclength (`` d/dl = \hat s\cdot\nabla ``, since the ray direction *is* `` \hat s ``), and using
that mixed partial derivatives of the smooth scalar field `` T `` commute, gives the compact
result

```math
\frac{d\vec s}{dl} = \nabla s. \qquad (3)
```

This is the fundamental ray-tracing equation: the slowness vector bends toward regions of
*higher* slowness (lower velocity) at a rate set by the local slowness gradient — exactly how
light bends toward the normal when entering a denser medium. This notebook solves equations (2)
and (3) numerically, stepping along the ray, to trace ray paths through the 2D heterogeneous
medium you paint below.

!!! tip "What equation (3) tells us"
	`` d\vec s/dl `` determines how the horizontal and vertical components of the slowness vector
	change along the ray path. For example, if `` \partial_x s `` is zero in a layered Earth
	medium, the horizontal slowness stays exactly constant — the classic statement of Snell's
	law along a ray in a layered medium.
"""

# ╔═╡ bc87032f-8941-4380-97f8-33dbebd9afe1
md"""
## Amplitudes Along a Ray

Writing `` A=\exp(\tilde A) `` turns the transport equation (eikA) into a linear equation for the
*log*-amplitude `` \tilde A ``. Projecting it along the ray direction `` \hat s `` (the same
`` d/dl=\hat s\cdot\nabla `` used above) gives

```math
\frac{d\tilde A}{dl} = -\frac{1}{2s}\,\nabla^2 T = -\frac{1}{2s}\,\nabla\cdot\vec s.
```

The right-hand side is the divergence of the slowness field along the ray: where rays converge
(`` \nabla\cdot\vec s<0 ``), amplitude grows; where they diverge, it decays. Integrating this
along the ray gives the amplitude via geometrical spreading — exactly the `` 1/r `` decay of a
spherical wavefront in a homogeneous medium, generalized to whatever raypath geometry the local
structure produces.
"""

# ╔═╡ 75c23a27-4253-474c-97bb-925ba236f362
md"""
## Fermat's Principle

Rays are exactly the paths that extremize the total travel time between two points,

```math
\mathbb{T} = \int_A^B s(x)\,dl,
```

among *all* possible paths connecting them — not just the true ray, which is why ray tracing
(shooting a path and integrating equations (2)-(3)) and travel-time extremization are two views
of the same physics. A full variational-calculus proof is a nice exercise: perturb the path by
`` \epsilon\,\eta(x) `` with `` \eta(A)=\eta(B)=0 ``, and set `` d\mathbb{T}/d\epsilon|_{\epsilon=0}=0 ``
to recover the ray equations derived above.
"""

# ╔═╡ 56cd6049-2266-4c7a-9e5f-1749f562034a
md"""
## Appendix
"""

# ╔═╡ 322e72d8-73aa-4386-af67-547ff74eec5e
md"""
### Background Velocity Models
"""

# ╔═╡ 075b648c-e275-44f4-94ea-73a8958f15dd
md"""
### Grid & Painted Perturbation
"""

# ╔═╡ 01e1393f-5b2a-423b-8a6a-0b8163cb1a1f
begin
    const XMAX = 500.0 # km -- must match the widget's world extent below
    const ZMAX = 150.0 # km
    const DX = 5.0      # km -- grid spacing (same in both directions, required by Eikonal.jl)
    const xgrid = range(0.0, XMAX; step=DX)
    const zgrid = range(0.0, ZMAX; step=DX)
    const NX = length(xgrid)
    const NZ = length(zgrid)
    const LVZ_WIDTH = 15.0 # km -- fixed Gaussian half-width of the "waveguide" low-velocity zone
    const FAN_SPREAD_DEG = 140.0 # deg -- fixed total angular width of the ray fan around the aim direction
end

# ╔═╡ 2d7556b5-27d7-40ec-b461-0945de83e5f7
begin
    """
    	velocity_homogeneous(z, x; v0)

    Constant background velocity `v0` (km/s), independent of position.
    """
    velocity_homogeneous(z, x; v0=6.0) = v0

    """
    	velocity_gradient(z, x; v0, grad)

    Background velocity increasing linearly with depth: `v0 + grad*z` (km/s), with `grad` in
    `(km/s)/km` -- the classic first-order approximation to a real Earth velocity profile.
    """
    velocity_gradient(z, x; v0=6.0, grad=0.02) = v0 + grad * z

    """
    	velocity_layered(z, x; v0, v1, zint)

    Two-layer background: constant `v0` (km/s) above the interface at depth `zint` (km),
    constant `v1` below it. A sharp velocity discontinuity is the classic setting for
    critically-refracted head waves -- rays bend abruptly instead of curving smoothly.
    """
    velocity_layered(z, x; v0=6.0, v1=8.0, zint=60.0) = z < zint ? v0 : v1

    """
    	velocity_waveguide(z, x; v0, grad, lvzdepth, lvzamp)

    Linear-gradient background `v0 + grad*z` with a Gaussian low-velocity zone of depth
    `lvzamp` (km/s) subtracted, centered at `lvzdepth` km with fixed half-width
    [`LVZ_WIDTH`](@ref). Rays that graze the channel get bent back into it and travel as
    trapped, guided paths -- the ray-theory picture of a waveguide.
    """
    velocity_waveguide(z, x; v0=6.0, grad=0.02, lvzdepth=60.0, lvzamp=1.5) =
        v0 + grad * z - lvzamp * exp(-((z - lvzdepth) / LVZ_WIDTH)^2)
end

# ╔═╡ 4097ad6c-d573-456c-822e-3eb922ecf0b8
"""
	build_velocity_grid(bgType, v0, grad, pert_flat; v1, zint, lvzdepth, lvzamp)

Combine a smooth background velocity model (`bgType` one of `"homogeneous"`, `"gradient"`,
`"layered"`, `"waveguide"` -- see [`velocity_homogeneous`](@ref), [`velocity_gradient`](@ref),
[`velocity_layered`](@ref), [`velocity_waveguide`](@ref)) with the user-painted perturbation
`pert_flat` (a flat vector, `NZ*NX` long, in row-major `(z outer, x inner)` order matching the
widget's own JS array layout) into one velocity grid, clamped to a small positive floor so the
Eikonal solve stays well-posed.

Returns the `(NZ, NX)` velocity matrix.
"""
function build_velocity_grid(bgType, v0, grad, pert_flat; v1=8.0, zint=60.0, lvzdepth=60.0, lvzamp=1.5)
    pert_grid = permutedims(reshape(pert_flat, NX, NZ))
    bgfun = if bgType == "gradient"
        (z, x) -> velocity_gradient(z, x; v0, grad)
    elseif bgType == "layered"
        (z, x) -> velocity_layered(z, x; v0, v1, zint)
    elseif bgType == "waveguide"
        (z, x) -> velocity_waveguide(z, x; v0, grad, lvzdepth, lvzamp)
    else
        (z, x) -> velocity_homogeneous(z, x; v0)
    end
    return [max(0.3, bgfun(z, x) + pert_grid[iz, ix]) for (iz, z) in enumerate(zgrid), (ix, x) in enumerate(xgrid)]
end

# ╔═╡ 08196372-2cc4-401f-87d7-b41d32364e5d
"""
	slowness_grid_and_gradients(v_grid)

Slowness grid `1 ./ v_grid` plus its spatial derivatives via centered finite differences
(one-sided at the domain edges) -- used instead of automatic differentiation since the field now
includes a grid-interpolated painted perturbation rather than a pure closed-form function.

Returns `(s_grid, sx_grid, sz_grid)`, each `(NZ, NX)`.
"""
function slowness_grid_and_gradients(v_grid)
    s_grid = inv.(v_grid)
    sx_grid = similar(s_grid)
    sz_grid = similar(s_grid)
    for ix in 1:NX, iz in 1:NZ
        ixm, ixp = max(1, ix - 1), min(NX, ix + 1)
        izm, izp = max(1, iz - 1), min(NZ, iz + 1)
        sx_grid[iz, ix] = (s_grid[iz, ixp] - s_grid[iz, ixm]) / ((ixp - ixm) * DX)
        sz_grid[iz, ix] = (s_grid[izp, ix] - s_grid[izm, ix]) / ((izp - izm) * DX)
    end
    return s_grid, sx_grid, sz_grid
end

# ╔═╡ f862769d-7fac-4408-80fd-aa3caf90c940
md"""
### Ray Tracing
"""

# ╔═╡ e8d8f418-7d68-4004-b294-b48fe82c7da9
"""
	get_raypaths(N, ds, Xsource, initial_slowness_vectors, pa)

Trace `length(initial_slowness_vectors)` rays outward from `Xsource = [z,x]`, each `N` steps of
arclength `ds`, by numerically integrating equations (2)-(3) above: `pa` bundles the interpolated
`slowness`/`slowness_x`/`slowness_z` functions and the grid extents used to detect when a ray
has left the domain. Each initial slowness vector's *direction* is kept but its magnitude is
rescaled to the true local slowness at the source (so only takeoff angle matters, not the
caller's chosen units). Travel time is accumulated along the way (`` \\int s\\,dl ``), since it's
needed both to animate the wavefront and to cross-check against the Eikonal solution below.

Returns a vector of `(path, traveltime)` named tuples: `path` a `2×N` matrix (rows `[z,x]`),
`traveltime` an `N`-vector, both `missing`/`NaN`-padded once a ray exits the domain.
"""
function get_raypaths(N, ds, Xsource, initial_slowness_vectors, pa)
    (; slowness, slowness_x, slowness_z, xgrid, zgrid) = pa

    map(initial_slowness_vectors) do slowness_vector
        slowness_vector = (slowness_vector ./ norm(slowness_vector)) .* slowness(Xsource[1], Xsource[2])

        Xraysave = Array{Any}(missing, 2, N)
        Tsave = fill(NaN, N)
        # Xsource may arrive as Int64 (whole-number values pushed through the widget's JS bond
        # deserialize as Int, not Float64) -- force Float64 so the position update below (which
        # writes a genuinely fractional value each step) never hits an InexactError.
        Xray = Float64.(deepcopy(Xsource))
        S = deepcopy(slowness_vector)
        t = 0.0
        for i = 1:N
            Xs = view(Xraysave, :, i)
            copyto!(Xs, Xray)
            Tsave[i] = t
            s_here = slowness(Xray[1], Xray[2])
            t += s_here * ds
            # equation 2
            Xray[1] = Xray[1] + (S[1] / s_here) * ds
            Xray[2] = Xray[2] + (S[2] / s_here) * ds
            # equation 3
            S[1] = S[1] + ds * slowness_z(Xray[1], Xray[2])
            S[2] = S[2] + ds * slowness_x(Xray[1], Xray[2])

            !(Xray[2] >= xgrid[1] && Xray[2] <= xgrid[end] && Xray[1] >= zgrid[1] && Xray[1] <= zgrid[end]) && break
        end
        return (path=Xraysave, traveltime=Tsave)
    end
end

# ╔═╡ 432388ad-4852-426c-956b-cc54710cb619
md"""
!!! correct "Checking the ray tracer against a straight-line analytic solution"
	In a homogeneous medium, rays are straight lines and the accumulated travel time to any point
	is just `` \text{distance}/v ``. `get_raypaths`'s accumulated `traveltime` is checked below
	against this closed form for an arbitrary interior point along a straight ray.
"""

# ╔═╡ 439b27c2-fd03-4391-b93c-4c58f3edd3ba
let
    v = 6.0 # homogeneous velocity, km/s
    pa = (; slowness=(z, x) -> 1 / v, slowness_x=(z, x) -> 0.0, slowness_z=(z, x) -> 0.0, xgrid, zgrid)
    srcZ, srcX = 0.0, 0.0
    θ = 30.0 # takeoff angle, degrees from horizontal
    ray = only(get_raypaths(200, 1.0, [srcZ, srcX], [[sind(θ), cosd(θ)]], pa))
    i = 100 # an arbitrary interior point along the ray
    z, x = ray.path[1, i], ray.path[2, i]
    predicted_t = hypot(z - srcZ, x - srcX) / v # straight ray in a homogeneous medium ⇒ t = distance/v
    @assert isapprox(ray.traveltime[i], predicted_t; rtol=1e-6)
    "get_raypaths' accumulated travel time matches distance/velocity for a straight ray ✓"
end

# ╔═╡ 05505f75-9720-4474-9388-250bae5f2fd2
"""
	ray_surface_return(path, traveltime)

Find where a traced ray (from [`get_raypaths`](@ref)) first returns to within 1 km of the
surface (`z<1`) after leaving the source, and the accumulated travel time there -- used to
compare against the Eikonal equation's own, independently-computed first-arrival time at that
same surface point (see the widget's "Travel-time check" panel).

Returns `(x, t)`, or `nothing` if the ray never returns to the surface within its traced length
(e.g. it exits the domain sideways, or the medium doesn't bend it back up). Skips the very first
point of the path so a source placed right at the surface doesn't trivially "return" immediately.
"""
function ray_surface_return(path, traveltime)
    for i in 2:length(traveltime)
        z = path[1, i]
        (ismissing(z) || isnan(traveltime[i])) && break
        z < 1.0 && return (x=path[2, i], t=traveltime[i])
    end
    return nothing
end

# ╔═╡ a9d965b7-0e29-4552-84e9-46b3bb8ace89
md"""
### Eikonal Fast-Sweep
"""

# ╔═╡ 6260a705-c8a8-4805-bd87-382dbc8a1d22
"""
	first_arrival_traveltimes(s_grid, srcZ, srcX)

Solve the Eikonal equation for the first-arrival travel time everywhere on the grid using the
fast sweeping method (`Eikonal.jl`), given a slowness grid and a source position in km.
"""
function first_arrival_traveltimes(s_grid, srcZ, srcX)
    fastsweep = FastSweeping(s_grid)
    init!(fastsweep, (argmin(abs.(srcZ .- zgrid)), argmin(abs.(srcX .- xgrid))))
    # `sweep!`'s default epsilon is eps(Float64) (~2.2e-16) -- convergence to full machine
    # precision, which costs ~200ms/~130MB of Gauss-Seidel sweeps on this grid for a value
    # only ever used to draw a wavefront/contour plot. epsilon=1e-6 (grid units, so ~2.7e-5s
    # of real travel-time error at this grid's 5km spacing) is still ~10 orders of magnitude
    # below anything visible, and cuts the cost to ~85ms.
    sweep!(fastsweep, verbose=false, epsilon=1e-6)
    # Eikonal.jl pads fastsweep.t to (NZ+1, NX+1) -- one row/col larger than the input
    # s_grid -- so it must be trimmed back to (NZ, NX) here, or every downstream consumer
    # (the widget's `ttFlat[iz*NX+ix]` indexing in particular) silently reads off-by-one
    # rows, which shows up as a sheared, seemingly-disconnected wavefront band.
    NZ_in, NX_in = size(s_grid)
    return (fastsweep.t .* DX)[1:NZ_in, 1:NX_in]
end

# ╔═╡ 450f8df0-4270-49f0-a910-c96c6bd89882
md"""
### Field Sampling for the Widget
"""

# ╔═╡ 94ad2aee-758d-4aca-8020-006ed515ebf9
_ray_flatten_rowmajor(M) = join(vec(permutedims(M)), ",")

# ╔═╡ 94853ecf-a2f0-4b9b-83ed-faf9869c84fd
"""
	ray_fan(N_RAYS, srcZ, srcX, angleMin, angleMax, pa)

Trace a fan of `N_RAYS` rays evenly spaced in takeoff angle between `angleMin` and `angleMax`
(degrees, measured from the `+x` axis toward `+z`/down -- `0°` horizontal, `90°` straight down),
each `500` steps of `1` km arclength. This is the single source of truth the widget draws from:
[`ray_fan_paths`](@ref) and [`ray_fan_surface_returns`](@ref) both extract what they need from
the *same* traced fan rather than re-tracing it.
"""
function ray_fan(N_RAYS, srcZ, srcX, angleMin, angleMax, pa)
    angles = N_RAYS == 1 ? [angleMin] : range(angleMin, angleMax; length=N_RAYS)
    dirs = [[sind(θ), cosd(θ)] for θ in angles]
    return get_raypaths(500, 1.0, [srcZ, srcX], dirs, pa)
end

# ╔═╡ c497fae9-3cf2-4697-9eff-7e4da5024adf
begin
    """
    	ray_fan_paths(rays)

    Extract just the `(z,x)` coordinates of each ray in a fan (from [`ray_fan`](@ref)), dropping
    the `missing`-padded tail once a ray exits the domain -- ready to push to the widget for
    drawing the ray-path polylines.
    """
    function ray_fan_paths(rays)
        map(rays) do ray
            valid = filter(x -> !any(ismissing.(x)), eachslice(ray.path, dims=2))
            (z=[p[1] for p in valid], x=[p[2] for p in valid])
        end
    end

    """
    	ray_fan_surface_returns(rays)

    The [`ray_surface_return`](@ref) point of every ray in a fan (from [`ray_fan`](@ref)) that
    actually returns to the surface, dropping the ones that don't (`nothing`).
    """
    ray_fan_surface_returns(rays) = filter(!isnothing, [ray_surface_return(r.path, r.traveltime) for r in rays])
end

# ╔═╡ 35124aad-70e8-42bf-a985-ebcb23d39a67
md"""
### The Interactive Widget
"""

# ╔═╡ eeaeaec0-75f0-4c57-9ce8-394d90b0ec91
begin
    struct RayPaintInput
        bgType::String # "homogeneous" | "gradient" | "layered" | "waveguide"
        v0::Float64
        grad::Float64
        v1::Float64      # "layered": velocity below the interface
        zint::Float64    # "layered": interface depth
        lvzdepth::Float64 # "waveguide": low-velocity-zone center depth
        lvzamp::Float64   # "waveguide": low-velocity-zone amplitude
        pert::Vector{Float64} # flat NZ*NX painted perturbation grid, km/s
        srcX::Float64
        srcZ::Float64
        angleAim::Float64 # takeoff angle of the fan's center, dragged directly from the source
        nRays::Int        # ray density: how many rays populate the (fixed-width) fan
    end
    RayPaintInput(; bgType="gradient", v0=6.0, grad=0.02, v1=8.0, zint=60.0, lvzdepth=60.0, lvzamp=1.5,
        pert=zeros(NZ * NX), srcX=20.0, srcZ=5.0, angleAim=45.0, nRays=15) =
        RayPaintInput(bgType, v0, grad, v1, zint, lvzdepth, lvzamp, pert, srcX, srcZ, angleAim, nRays)

    Base.get(w::RayPaintInput) = Dict{String,Any}(
        "bgType" => w.bgType, "v0" => w.v0, "grad" => w.grad, "v1" => w.v1, "zint" => w.zint,
        "lvzdepth" => w.lvzdepth, "lvzamp" => w.lvzamp, "pert" => w.pert,
        "srcX" => w.srcX, "srcZ" => w.srcZ, "angleAim" => w.angleAim, "nRays" => w.nRays)

    function Base.show(io::IO, ::MIME"text/html", w::RayPaintInput)
        write(io, """
        <div id="raywidget">
        <style>
        pluto-cell:has(#raywidget) { width: min(80vw, 1500px) !important;
          margin-left: calc((100% - min(80vw, 1500px)) / 2) !important; }
        #raywidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #raywidget .ray-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #raywidget .ray-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #raywidget .ray-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #raywidget .ray-workspace{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start}
        #raywidget .ray-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #raywidget .ray-caption{font-size:13px;color:#9ca3af;text-align:center;margin-top:4px}
        #raywidget canvas{display:block;cursor:crosshair}
        #raywidget .ray-controls{flex:0 0 540px;width:540px;display:grid;
          grid-template-columns:repeat(2, minmax(0,1fr));gap:8px;font:14px sans-serif;align-content:start}
        #raywidget .ray-control-group{background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
        #raywidget .ray-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #raywidget .ray-control-row{display:grid;grid-template-columns:70px minmax(0,1fr) 58px;gap:6px;align-items:center;margin:5px 0}
        #raywidget .ray-control-row label{font-size:13px;color:#9ca3af}
        #raywidget .ray-control-row input[type=range]{width:100%;min-width:0}
        #raywidget .ray-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #raywidget .ray-actions{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
        #raywidget .ray-readout{font-size:13px;color:#d1d5db;line-height:1.6}
        #raywidget .ray-readout b{color:#e5e7eb}
        #raywidget select{background:#0b0b0b;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:3px;width:100%}
        #raywidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:13px;cursor:pointer}
        #raywidget button.active{background:#2563eb;border-color:#93c5fd}
        #raywidget .ray-tt-panel{width:100%;box-sizing:border-box;margin-top:14px}
        #raywidget .ray-cbar-label{font-size:11px;color:#9ca3af;margin:8px 0 2px;text-align:center}
        #raywidget .ray-cbar{cursor:default;display:block;margin:0 auto}
        </style>

        <div class="ray-title">
          <div class="ray-title-desc">Paint velocity perturbations and watch how rays and wavefronts bend through them.</div>
          <div class="ray-title-hint">drag on the map to paint &middot; drag the star to move the source &middot; drag its spoke to aim and shoot a ray &middot; press Play to sweep the wavefront</div>
        </div>

        <div class="ray-workspace">
          <div>
            <div class="ray-panel"><canvas id="ray-map"></canvas></div>
            <div class="ray-caption" id="ray-caption"></div>
            <div class="ray-cbar-label" id="ray-cbar-label">velocity perturbation from background, %</div>
            <canvas id="ray-cbar" class="ray-cbar"></canvas>
          </div>
          <div class="ray-controls">
            <div class="ray-control-group">
              <div class="ray-control-title">Map shows</div>
              <div class="ray-actions">
                <button id="ray-view-pert" type="button">Perturbation</button>
                <button id="ray-view-bg" type="button">Background</button>
              </div>
            </div>
            <div class="ray-control-group">
              <div class="ray-control-title">Rays</div>
              <div class="ray-control-row"><label>density</label><input type="range" id="ray-density" min="1" max="41" step="2" value="$(w.nRays)"><span class="ray-value" id="ray-density-v"></span></div>
            </div>
            <div class="ray-control-group">
              <div class="ray-control-title">Paint</div>
              <div class="ray-actions">
                <button id="ray-paint-slow" type="button">Slower (red)</button>
                <button id="ray-paint-fast" type="button">Faster (blue)</button>
              </div>
              <div class="ray-control-row"><label>brush</label><input type="range" id="ray-brush" min="10" max="80" step="5" value="30"><span class="ray-value" id="ray-brush-v"></span></div>
              <div class="ray-actions"><button id="ray-clear" type="button">Clear paint</button></div>
            </div>
            <div class="ray-control-group">
              <div class="ray-control-title">Background</div>
              <div class="ray-control-row"><label>model</label>
                <select id="ray-bgtype">
                  <option value="gradient">Linear gradient</option>
                  <option value="homogeneous">Homogeneous</option>
                  <option value="layered">Layered (interface)</option>
                  <option value="waveguide">Low-velocity waveguide</option>
                </select><span></span>
              </div>
              <div class="ray-control-row"><label>v&#8320; (km/s)</label><input type="range" id="ray-v0" min="2" max="8" step="0.1" value="$(w.v0)"><span class="ray-value" id="ray-v0-v"></span></div>
              <div class="ray-control-row" id="ray-row-grad"><label>gradient</label><input type="range" id="ray-grad" min="0" max="0.06" step="0.002" value="$(w.grad)"><span class="ray-value" id="ray-grad-v"></span></div>
              <div class="ray-control-row" id="ray-row-v1"><label>v&#8321; (km/s)</label><input type="range" id="ray-v1" min="2" max="10" step="0.1" value="$(w.v1)"><span class="ray-value" id="ray-v1-v"></span></div>
              <div class="ray-control-row" id="ray-row-zint"><label>interface z</label><input type="range" id="ray-zint" min="10" max="$(ZMAX-10)" step="5" value="$(w.zint)"><span class="ray-value" id="ray-zint-v"></span></div>
              <div class="ray-control-row" id="ray-row-lvzdepth"><label>LVZ depth</label><input type="range" id="ray-lvzdepth" min="10" max="$(ZMAX-10)" step="5" value="$(w.lvzdepth)"><span class="ray-value" id="ray-lvzdepth-v"></span></div>
              <div class="ray-control-row" id="ray-row-lvzamp"><label>LVZ amp</label><input type="range" id="ray-lvzamp" min="0.2" max="3" step="0.1" value="$(w.lvzamp)"><span class="ray-value" id="ray-lvzamp-v"></span></div>
            </div>
            <div class="ray-control-group">
              <div class="ray-control-title">Wavefront</div>
              <div class="ray-actions"><button id="ray-play" type="button">Play</button><button id="ray-reset" type="button">Reset</button></div>
            </div>
            <div class="ray-control-group">
              <div class="ray-control-title">Readouts</div>
              <div class="ray-readout" id="ray-readout"></div>
            </div>
          </div>
        </div>

        <div class="ray-tt-panel">
          <div class="ray-panel"><canvas id="ray-tt-compare"></canvas></div>
          <div class="ray-caption">travel time at the surface: eikonal first-arrival (line) vs. ray-traced (dots) &mdash; two independent methods, same physics</div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        const NX=$(NX), NZ=$(NZ), XMAX=$(XMAX), ZMAX=$(ZMAX), DXKM=$(DX), LVZWIDTH=$(LVZ_WIDTH), FANSPREAD=$(FAN_SPREAD_DEG);

        let bgType="$(w.bgType)", v0=$(w.v0), grad=$(w.grad), v1=$(w.v1), zint=$(w.zint), lvzdepth=$(w.lvzdepth), lvzamp=$(w.lvzamp);
        let pert = $(w.pert == zeros(NZ*NX) ? "new Array(NX*NZ).fill(0)" : "[" * join(w.pert, ",") * "]");
        let srcX=$(w.srcX), srcZ=$(w.srcZ), angleAim=$(w.angleAim), nRays=$(w.nRays);

        const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1500);
        const CONTROLS_W = 540+16; // matches .ray-controls' 2-column grid width -- keep in sync
        const SEC_W = Math.max(420, Math.min(availW - CONTROLS_W, 760));
        const SEC_H = Math.round(SEC_W * ZMAX/XMAX);
        const DPR = window.devicePixelRatio || 1;

        const cv = par.querySelector('#ray-map');
        const ctx = cv.getContext('2d');
        function hidpi(canvas, context, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          context.setTransform(DPR,0,0,DPR,0,0);
        }
        hidpi(cv, ctx, SEC_W, SEC_H);

        const ttCv = par.querySelector('#ray-tt-compare');
        const ttCtx = ttCv.getContext('2d');
        const TT_W = SEC_W + CONTROLS_W - 16, TT_H = 170;
        hidpi(ttCv, ttCtx, TT_W, TT_H);

        const cbarCv = par.querySelector('#ray-cbar');
        const cbarCtx = cbarCv.getContext('2d');
        const CBAR_W = SEC_W, CBAR_H = 34;
        hidpi(cbarCv, cbarCtx, CBAR_W, CBAR_H);

        function toScreen(x, z){ return [x/XMAX*SEC_W, z/ZMAX*SEC_H]; }
        function toWorld(px, pz){ return [px/SEC_W*XMAX, pz/SEC_H*ZMAX]; }

        function velBg(z, x){
          if(bgType==='gradient') return v0 + grad*z;
          if(bgType==='layered') return z < zint ? v0 : v1;
          if(bgType==='waveguide') return v0 + grad*z - lvzamp*Math.exp(-(((z-lvzdepth)/LVZWIDTH)**2));
          return v0;
        }
        function pertAt(iz, ix){ return pert[iz*NX + ix]; }
        function xOf(ix){ return ix/(NX-1)*XMAX; }
        function zOf(iz){ return iz/(NZ-1)*ZMAX; }

        // perturbation as a percentage of the LOCAL background velocity at that grid point --
        // works for all four background models (a flat % clamp wouldn't, since "faster" means
        // something different in km/s near a 2 km/s layer than near an 8 km/s one).
        function pctAt(iz, ix){ return 100 * pertAt(iz, ix) / velBg(zOf(iz), xOf(ix)); }
        const FIELD_PCT_MAX = 15; // fixed colorbar scale: painted perturbation saturates at +-15%
        let fieldView = 'pert'; // 'pert' (painted perturbation, %) | 'background' (absolute v, km/s)

        // diverging colormap: `mx` sets the +-range mapped to full color saturation, `mid` sets
        // what value reads as the neutral center -- blue=faster/higher, red=slower/lower (this
        // repo's usual tomography convention). Reused for both the %-perturbation view (mid=0,
        // mx=FIELD_PCT_MAX) and the absolute-background view (mid/mx = the background's own
        // domain mean/half-range), so both views and their colorbars share one implementation.
        function velColor(v, mx, mid=0){
          const t = Math.max(-1, Math.min(1, (v-mid)/mx));
          if(t >= 0) return [Math.round(255*(1-t)), Math.round(255*(1-t)), 255];
          const s = -t;
          return [255, Math.round(255*(1-s)), Math.round(255*(1-s))];
        }

        // Sampled directly from velColor -- what the bar shows always matches what the map
        // actually paints, including at the clamped extremes.
        function drawFieldColorbar(mx, mid){
          const w = CBAR_W, h = CBAR_H;
          cbarCtx.clearRect(0,0,w,h);
          const barH = 14, barY = 2;
          for(let i=0;i<w;i++){
            const t = i/(w-1);
            const v = mid + (-mx + t*2*mx);
            const [r,g,b] = velColor(v, mx, mid);
            cbarCtx.fillStyle = 'rgb('+r+','+g+','+b+')';
            cbarCtx.fillRect(i, barY, 1, barH);
          }
          cbarCtx.strokeStyle = '#4b5563'; cbarCtx.lineWidth = 1;
          cbarCtx.strokeRect(0.5, barY+0.5, w-1, barH-1);
          cbarCtx.fillStyle = '#9ca3af'; cbarCtx.font = '10px sans-serif';
          const fmt = fieldView==='background' ? (v => v.toFixed(2)+' km/s') :
            (v => (Math.abs(v)<1e-9 ? '0' : (v>=0?'+':'')+v.toFixed(1))+'%');
          cbarCtx.textAlign = 'left'; cbarCtx.fillText(fmt(mid-mx), 0, barY+barH+12);
          cbarCtx.textAlign = 'right'; cbarCtx.fillText(fmt(mid+mx), w, barY+barH+12);
          cbarCtx.textAlign = 'center'; cbarCtx.fillText(fmt(mid), w/2, barY+barH+12);
          cbarCtx.textAlign = 'left';
        }

        let waveData = null;    // {ttFlat} pushed from Julia -- full (NZ,NX) eikonal grid, flat
        let rayPaths = null;    // [{z:[...], x:[...]}] pushed from Julia -- the current fan, always live
        let surfaceReturns = null; // {x:[...], t:[...]} pushed from Julia (ray-traced check points)
        let playing = false, rafId = null, tPhase = 0, lastTs = null;
        let paintMode = 'slow'; // 'slow' (red, subtract velocity) | 'fast' (blue, add velocity)

        function drawField(){
          let mx = FIELD_PCT_MAX, mid = 0;
          if(fieldView==='background'){
            let bgMin=Infinity, bgMax=-Infinity;
            for(let iz=0; iz<NZ; iz++) for(let ix=0; ix<NX; ix++){
              const v = velBg(zOf(iz), xOf(ix));
              bgMin = Math.min(bgMin, v); bgMax = Math.max(bgMax, v);
            }
            mid = (bgMin+bgMax)/2;
            mx = Math.max(0.1, (bgMax-bgMin)/2); // floor avoids a degenerate 0-width bar for homogeneous
          }
          const img = ctx.createImageData(Math.round(SEC_W*DPR), Math.round(SEC_H*DPR));
          const wpx = img.width, hpx = img.height;
          for(let py=0; py<hpx; py++){
            const zkm = py/hpx*ZMAX;
            const iz = Math.min(NZ-1, Math.round(zkm/ZMAX*(NZ-1)));
            for(let px2=0; px2<wpx; px2++){
              const xkm = px2/wpx*XMAX;
              const ix = Math.min(NX-1, Math.round(xkm/XMAX*(NX-1)));
              const v = fieldView==='background' ? velBg(zOf(iz), xOf(ix)) : pctAt(iz, ix);
              const [r,g,b] = velColor(v, mx, mid);
              const idx = (py*wpx+px2)*4;
              img.data[idx]=r; img.data[idx+1]=g; img.data[idx+2]=b; img.data[idx+3]=255;
            }
          }
          ctx.putImageData(img, 0, 0);
          drawFieldColorbar(mx, mid);
          par.querySelector('#ray-cbar-label').textContent = fieldView==='background' ?
            'background velocity, km/s' : 'velocity perturbation from background, %';
        }

        function drawWavefrontSweep(){
          if(!waveData) return;
          ctx.save();
          for(let iz=0; iz<NZ; iz++){
            for(let ix=0; ix<NX; ix++){
              const tt = waveData.ttFlat[iz*NX+ix];
              const d = Math.abs(tt - tPhase);
              if(d < 2.2){
                const a = Math.max(0, 1 - d/2.2);
                const [x0,z0] = toScreen(xOf(ix)-DXKM/2, zOf(iz)-DXKM/2);
                const [x1,z1] = toScreen(xOf(ix)+DXKM/2, zOf(iz)+DXKM/2);
                ctx.fillStyle = 'rgba(250,204,21,' + (0.75*a) + ')';
                ctx.fillRect(x0, z0, x1-x0, z1-z0);
              }
            }
          }
          ctx.restore();
        }

        function drawRays(){
          if(!rayPaths) return;
          ctx.strokeStyle = 'rgba(17,24,39,0.9)'; ctx.lineWidth = 3.2;
          for(const ray of rayPaths){
            ctx.beginPath();
            for(let i=0;i<ray.x.length;i++){
              const [px,pz] = toScreen(ray.x[i], ray.z[i]);
              i===0 ? ctx.moveTo(px,pz) : ctx.lineTo(px,pz);
            }
            ctx.stroke();
          }
          ctx.strokeStyle = '#f5f3ef'; ctx.lineWidth = 1.4;
          for(const ray of rayPaths){
            ctx.beginPath();
            for(let i=0;i<ray.x.length;i++){
              const [px,pz] = toScreen(ray.x[i], ray.z[i]);
              i===0 ? ctx.moveTo(px,pz) : ctx.lineTo(px,pz);
            }
            ctx.stroke();
          }
        }

        // Per this repo's marker-shape convention (pluto-widget-style skill, "Marker shapes"):
        // a seismic source is drawn as a filled star, copied verbatim from
        // lame-theorem.jl's PointForceRadiationInput rather than a plain circle.
        function drawStarMarker(cx, cy, r, fill, stroke){
          const spikes = 5, rOuter = r, rInner = r * 0.45;
          ctx.beginPath();
          for(let i=0; i<spikes*2; i++){
            const rad = i % 2 === 0 ? rOuter : rInner;
            const ang = -Math.PI/2 + i*Math.PI/spikes;
            const x = cx + rad*Math.cos(ang), y = cy + rad*Math.sin(ang);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.closePath();
          ctx.fillStyle = fill; ctx.fill();
          ctx.strokeStyle = stroke; ctx.lineWidth = 1; ctx.stroke();
        }

        // Same arrowhead helper as lame-theorem.jl's PointForceRadiationInput, copied verbatim.
        function drawHandles(){
          const [spx,spz] = toScreen(srcX, srcZ);
          const SPOKE_LEN = 40;

          // faint dashed guide lines mark the fan's outer edges (fixed width, not user-adjustable)
          // whenever more than one ray is shown -- not draggable, just context for the density slider.
          if(nRays > 1){
            ctx.strokeStyle = 'rgba(250,204,21,0.35)'; ctx.lineWidth = 1.5; ctx.setLineDash([4,3]);
            for(const ang of [(angleAim-FANSPREAD/2)*Math.PI/180, (angleAim+FANSPREAD/2)*Math.PI/180]){
              const hx = spx + SPOKE_LEN*Math.cos(ang), hz = spz + SPOKE_LEN*Math.sin(ang);
              ctx.beginPath(); ctx.moveTo(spx,spz); ctx.lineTo(hx,hz); ctx.stroke();
            }
            ctx.setLineDash([]);
          }

          // the one draggable spoke: drag its tip to aim the fan in any direction.
          const aimR = angleAim*Math.PI/180;
          const hx = spx + SPOKE_LEN*Math.cos(aimR), hz = spz + SPOKE_LEN*Math.sin(aimR);
          ctx.strokeStyle = 'rgba(250,204,21,0.7)'; ctx.lineWidth = 2;
          ctx.beginPath(); ctx.moveTo(spx,spz); ctx.lineTo(hx,hz); ctx.stroke();
          ctx.beginPath(); ctx.arc(hx,hz,6,0,2*Math.PI);
          ctx.fillStyle = '#facc15'; ctx.fill(); ctx.strokeStyle='#000'; ctx.lineWidth=1; ctx.stroke();

          drawStarMarker(spx, spz, 9, '#f5f3ef', '#0a0f18');
          ctx.fillStyle = '#e5e7eb'; ctx.font='12px sans-serif'; ctx.textAlign='left';
          ctx.fillText('source', spx+12, spz-9);
        }

        function draw(){
          drawField();
          drawWavefrontSweep();
          drawRays();
          drawHandles();
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1; ctx.strokeRect(0.5,0.5,SEC_W-1,SEC_H-1);
          par.querySelector('#ray-caption').textContent = 'x: 0–' + XMAX.toFixed(0) + ' km  ·  z (depth): 0–' + ZMAX.toFixed(0) + ' km';
          updateReadout();
          drawTTCompare();
        }

        function drawTTCompare(){
          ttCtx.clearRect(0,0,TT_W,TT_H);
          ttCtx.strokeStyle = '#374151'; ttCtx.lineWidth = 1; ttCtx.strokeRect(0.5,0.5,TT_W-1,TT_H-1);
          if(!waveData){
            ttCtx.fillStyle = '#6b7280'; ttCtx.font = '12px sans-serif'; ttCtx.fillText('computing...', 10, 18);
            return;
          }
          // eikonal first-arrival time along the surface (z=0 row) is exactly the first NX
          // entries of ttFlat, since ttFlat is iz-outer/ix-inner -- a plain array slice, no
          // new physics computed here.
          const ttSurface = waveData.ttFlat.slice(0, NX);
          let tmax = 1e-6;
          for(let i=0;i<NX;i++) tmax = Math.max(tmax, ttSurface[i]);
          if(surfaceReturns) for(const t of surfaceReturns.t) tmax = Math.max(tmax, t);

          const PAD = 34;
          function xOf2(x){ return PAD + (x/XMAX)*(TT_W-PAD-10); }
          function yOf2(t){ return TT_H-20 - (t/tmax)*(TT_H-34); }

          ttCtx.strokeStyle = '#f97316'; ttCtx.lineWidth = 1.8;
          ttCtx.beginPath();
          for(let ix=0; ix<NX; ix++){
            const px = xOf2(xOf(ix)), py = yOf2(ttSurface[ix]);
            ix===0 ? ttCtx.moveTo(px,py) : ttCtx.lineTo(px,py);
          }
          ttCtx.stroke();

          if(surfaceReturns){
            ttCtx.fillStyle = '#38bdf8';
            for(let i=0;i<surfaceReturns.x.length;i++){
              const px = xOf2(surfaceReturns.x[i]), py = yOf2(surfaceReturns.t[i]);
              ttCtx.beginPath(); ttCtx.arc(px,py,3.2,0,2*Math.PI); ttCtx.fill();
            }
          }

          ttCtx.fillStyle = '#6b7280'; ttCtx.font = '11px sans-serif'; ttCtx.textAlign='left';
          ttCtx.fillText('x=0 km', PAD, TT_H-6);
          ttCtx.textAlign='right'; ttCtx.fillText('x=' + XMAX.toFixed(0) + ' km', TT_W-6, TT_H-6);
          ttCtx.textAlign='left'; ttCtx.fillStyle='#f97316'; ttCtx.fillText('eikonal', PAD, 12);
          ttCtx.fillStyle='#38bdf8'; ttCtx.fillText('ray-traced', PAD+55, 12);
        }

        function updateReadout(){
          par.querySelector('#ray-readout').innerHTML =
            'source <b>(' + srcX.toFixed(0) + ', ' + srcZ.toFixed(0) + ')</b> km<br>' +
            'aim <b>' + angleAim.toFixed(0) + '&deg;</b>' +
            (nRays > 1 ? ' &middot; <b>' + nRays + '</b> rays' : ' (1 ray)') + '<br>' +
            (waveData ? 'sweep t <b>' + tPhase.toFixed(1) + '</b> s' : 'sweep: press Play');
        }

        function syncControls(){
          par.querySelector('#ray-bgtype').value = bgType;
          par.querySelector('#ray-v0').value = v0; par.querySelector('#ray-v0-v').textContent = v0.toFixed(1);
          par.querySelector('#ray-grad').value = grad; par.querySelector('#ray-grad-v').textContent = grad.toFixed(3);
          par.querySelector('#ray-v1').value = v1; par.querySelector('#ray-v1-v').textContent = v1.toFixed(1);
          par.querySelector('#ray-zint').value = zint; par.querySelector('#ray-zint-v').textContent = zint.toFixed(0)+' km';
          par.querySelector('#ray-lvzdepth').value = lvzdepth; par.querySelector('#ray-lvzdepth-v').textContent = lvzdepth.toFixed(0)+' km';
          par.querySelector('#ray-lvzamp').value = lvzamp; par.querySelector('#ray-lvzamp-v').textContent = lvzamp.toFixed(1);
          par.querySelector('#ray-density').value = nRays;
          par.querySelector('#ray-density-v').textContent = nRays;
          par.querySelector('#ray-row-grad').style.display = bgType==='gradient' ? '' : 'none';
          par.querySelector('#ray-row-v1').style.display = bgType==='layered' ? '' : 'none';
          par.querySelector('#ray-row-zint').style.display = bgType==='layered' ? '' : 'none';
          par.querySelector('#ray-row-lvzdepth').style.display = bgType==='waveguide' ? '' : 'none';
          par.querySelector('#ray-row-lvzamp').style.display = bgType==='waveguide' ? '' : 'none';
          par.querySelector('#ray-paint-slow').classList.toggle('active', paintMode==='slow');
          par.querySelector('#ray-paint-fast').classList.toggle('active', paintMode==='fast');
          par.querySelector('#ray-view-pert').classList.toggle('active', fieldView==='pert');
          par.querySelector('#ray-view-bg').classList.toggle('active', fieldView==='background');
        }

        let commitInFlight = false;
        function commit(){
          commitInFlight = true;
          par.value = { bgType, v0, grad, v1, zint, lvzdepth, lvzamp, pert, srcX, srcZ, angleAim, nRays };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledCommit(){ if(!commitInFlight) commit(); }

        par.addEventListener('ray-update', e=>{
          rayPaths = e.detail.rayPaths;
          waveData = { ttFlat: e.detail.ttFlat };
          surfaceReturns = { x: e.detail.surfaceReturnsX, t: e.detail.surfaceReturnsT };
          commitInFlight = false;
          draw();
        });

        // ---- painting ----
        let painting = false, dragMode = null, lastX=0, lastZ=0;
        const PAINT_STEP_PCT = 0.01; // each brush touch adds/removes 1% of local velocity at the brush center
        function paintAt(xkm, zkm, brushKm){
          const sign = paintMode==='fast' ? 1 : -1;
          for(let iz=0; iz<NZ; iz++){
            const dz = zOf(iz)-zkm;
            for(let ix=0; ix<NX; ix++){
              const dx = xOf(ix)-xkm;
              const d2 = dx*dx+dz*dz;
              const r2 = brushKm*brushKm;
              if(d2 < r2*4){
                const vLocal = velBg(zOf(iz), xOf(ix));
                const falloff = Math.exp(-d2/(2*r2/4));
                const amp = falloff * PAINT_STEP_PCT * vLocal * sign;
                const cap = (FIELD_PCT_MAX/100) * vLocal; // saturate accumulated paint at +-FIELD_PCT_MAX
                const idx = iz*NX+ix;
                pert[idx] = Math.max(-cap, Math.min(cap, pert[idx] + amp));
              }
            }
          }
        }

        cv.addEventListener('mousedown', e=>{
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left, pz = e.clientY-rect.top;
          const [spx,spz] = toScreen(srcX, srcZ);
          const SPOKE_LEN = 40;
          const aimR = angleAim*Math.PI/180;
          const aimH = [spx+SPOKE_LEN*Math.cos(aimR), spz+SPOKE_LEN*Math.sin(aimR)];
          const dSrc = Math.hypot(px-spx, pz-spz);
          const dAim = Math.hypot(px-aimH[0], pz-aimH[1]);
          let best=null, bestD=14;
          if(dSrc<bestD){bestD=dSrc; best='source';}
          if(dAim<bestD){bestD=dAim; best='aim';}
          if(best){ dragMode = best; return; }
          painting = true;
          const [xkm,zkm] = toWorld(px,pz);
          const brushKm = +par.querySelector('#ray-brush').value;
          paintAt(xkm, zkm, brushKm);
          draw();
        });
        window.addEventListener('mousemove', e=>{
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left, pz = e.clientY-rect.top;
          if(dragMode){
            const [xkm,zkm] = toWorld(Math.max(0,Math.min(SEC_W,px)), Math.max(0,Math.min(SEC_H,pz)));
            if(dragMode==='source'){
              srcX = Math.max(0, Math.min(XMAX, xkm)); srcZ = Math.max(0, Math.min(ZMAX, zkm));
            } else {
              const [spx,spz] = toScreen(srcX, srcZ);
              let ang = Math.atan2(pz-spz, px-spx)*180/Math.PI;
              if(ang < -90) ang += 360;
              angleAim = Math.max(-90, Math.min(270, ang));
            }
            draw(); throttledCommit();
            return;
          }
          if(!painting) return;
          const [xkm,zkm] = toWorld(px,pz);
          const brushKm = +par.querySelector('#ray-brush').value;
          paintAt(xkm, zkm, brushKm);
          draw();
        });
        window.addEventListener('mouseup', ()=>{
          if(painting || dragMode) commit();
          painting = false; dragMode = null;
        });

        par.querySelector('#ray-paint-slow').addEventListener('click', ()=>{ paintMode='slow'; syncControls(); });
        par.querySelector('#ray-paint-fast').addEventListener('click', ()=>{ paintMode='fast'; syncControls(); });
        par.querySelector('#ray-view-pert').addEventListener('click', ()=>{ fieldView='pert'; syncControls(); draw(); });
        par.querySelector('#ray-view-bg').addEventListener('click', ()=>{ fieldView='background'; syncControls(); draw(); });
        par.querySelector('#ray-clear').addEventListener('click', ()=>{
          pert = new Array(NX*NZ).fill(0); draw(); commit();
        });

        par.addEventListener('input', e=>{
          if(e.target===par) return;
          e.stopImmediatePropagation();
          const id = e.target.id, v = e.target.value;
          if(id==='ray-v0'){ v0=+v; par.querySelector('#ray-v0-v').textContent=v0.toFixed(1); }
          else if(id==='ray-grad'){ grad=+v; par.querySelector('#ray-grad-v').textContent=grad.toFixed(3); }
          else if(id==='ray-v1'){ v1=+v; par.querySelector('#ray-v1-v').textContent=v1.toFixed(1); }
          else if(id==='ray-zint'){ zint=+v; par.querySelector('#ray-zint-v').textContent=zint.toFixed(0)+' km'; }
          else if(id==='ray-lvzdepth'){ lvzdepth=+v; par.querySelector('#ray-lvzdepth-v').textContent=lvzdepth.toFixed(0)+' km'; }
          else if(id==='ray-lvzamp'){ lvzamp=+v; par.querySelector('#ray-lvzamp-v').textContent=lvzamp.toFixed(1); }
          else if(id==='ray-density'){ nRays=+v; par.querySelector('#ray-density-v').textContent = nRays; }
          else if(id==='ray-brush'){ par.querySelector('#ray-brush-v').textContent=v+' km'; return; }
          else return;
          draw(); throttledCommit();
        }, true);

        par.addEventListener('change', e=>{
          if(e.target===par) return;
          e.stopImmediatePropagation();
          const id = e.target.id;
          if(id==='ray-bgtype'){ bgType = e.target.value; syncControls(); draw(); commit(); return; }
          if(id==='ray-v0'||id==='ray-grad'||id==='ray-v1'||id==='ray-zint'||id==='ray-lvzdepth'||id==='ray-lvzamp'||id==='ray-density'){ draw(); commit(); return; }
        }, true);

        function tLoopMax(){ return waveData ? Math.max(...waveData.ttFlat)*1.15 : 60; }
        const SIM_SPEED = 6;
        const playBtn = par.querySelector('#ray-play');
        function stepAnim(ts){
          if(lastTs===null) lastTs = ts;
          const dt = Math.min(0.1, (ts-lastTs)/1000);
          lastTs = ts;
          tPhase += dt*SIM_SPEED;
          if(tPhase > tLoopMax()) tPhase = 0;
          draw();
          rafId = requestAnimationFrame(stepAnim);
        }
        playBtn.addEventListener('click', ()=>{
          playing = !playing;
          playBtn.textContent = playing ? 'Pause' : 'Play';
          if(playing){ lastTs=null; rafId = requestAnimationFrame(stepAnim); }
          else if(rafId){ cancelAnimationFrame(rafId); rafId=null; }
        });

        par.querySelector('#ray-reset').addEventListener('click', ()=>{
          bgType="$(w.bgType)"; v0=$(w.v0); grad=$(w.grad);
          v1=$(w.v1); zint=$(w.zint); lvzdepth=$(w.lvzdepth); lvzamp=$(w.lvzamp);
          pert = new Array(NX*NZ).fill(0);
          srcX=$(w.srcX); srcZ=$(w.srcZ); angleAim=$(w.angleAim); nRays=$(w.nRays);
          tPhase=0;
          syncControls(); draw(); commit();
        });

        par.querySelector('#ray-brush-v').textContent = par.querySelector('#ray-brush').value+' km';
        syncControls(); draw();
        }
        </script>

        """)
    end

    const _ray_ready = true
end

# ╔═╡ e6cd1e51-7666-49e9-8d8c-523abf6fe6c9
begin
    _ray_ready
    WideCell(@bind ray RayPaintInput(); max_width=1500)
end

# ╔═╡ cddc22f1-85fc-489d-85d3-d406d795e5a3
# The bond starts as `nothing` until the widget's first real interaction in a live browser
# reports back -- fall back to the same defaults the widget itself opens with. A previously
# connected browser tab can also hand back a dict from *before* v1/zint/lvzdepth/lvzamp were
# added to the widget's own commit() -- read those four with `get(...)` + a default everywhere
# below, not `ray_safe["..."]`, so a stale dict from an older session can't KeyError.
ray_safe = ray isa AbstractDict ? ray : Dict{String,Any}(
    "bgType" => "gradient", "v0" => 6.0, "grad" => 0.02, "v1" => 8.0, "zint" => 60.0,
    "lvzdepth" => 60.0, "lvzamp" => 1.5, "pert" => zeros(NZ * NX),
    "srcX" => 20.0, "srcZ" => 5.0, "angleAim" => 45.0, "nRays" => 15)

# ╔═╡ 496f5399-68e2-4a83-a546-324a9de63c09
# Pluto's static analysis of an `md"""..."""` cell chokes when it contains several `$(...)`
# interpolations alongside non-ASCII characters (`·`, `°`) -- it grabs the wrong-length source
# substring and re-parses it, throwing a spurious ParseError. Sidestep the whole class of bug by
# building one plain Julia string first and handing the md block a single bare `$readout`.
let
    srcX_r = round(ray_safe["srcX"]; digits=0)
    srcZ_r = round(ray_safe["srcZ"]; digits=0)
    aim_r = round(get(ray_safe, "angleAim", 45.0); digits=0)
    nrays_r = max(1, round(Int, get(ray_safe, "nRays", 15)))
    bgType = ray_safe["bgType"]
    bgText = if bgType == "gradient"
        "**linear gradient**, **$(ray_safe["v0"])** km/s, gradient **$(ray_safe["grad"])** (km/s)/km"
    elseif bgType == "layered"
        "**layered**, **$(ray_safe["v0"])** km/s above / **$(get(ray_safe, "v1", 8.0))** km/s below " *
        "an interface at **$(get(ray_safe, "zint", 60.0))** km"
    elseif bgType == "waveguide"
        "**low-velocity waveguide**, **$(ray_safe["v0"])** km/s gradient with a **$(get(ray_safe, "lvzamp", 1.5))** km/s dip " *
        "at **$(get(ray_safe, "lvzdepth", 60.0))** km"
    else
        "**homogeneous**, **$(ray_safe["v0"])** km/s"
    end
    rayText = nrays_r > 1 ? "aim **$(aim_r)°**, **$(nrays_r)** rays" : "aim **$(aim_r)°** (1 ray)"
    readout = "Background: " * bgText * " · " *
              "source at **($(srcX_r), $(srcZ_r))** km · " * rayText
    md"$readout"
end

# ╔═╡ 1dc9add4-83d4-4516-855b-4531208d8dee
md"""
`RayPush` does no physics -- it takes the already-traced ray paths, the eikonal traveltime grid,
and the ray-traced surface-return check points, and hands them to the *already-rendered*
[`RayPaintInput`](@ref) widget by dispatching a browser `CustomEvent`, the same pattern this
repo's other widgets use (see e.g. `FieldPush` in `fault-dislocation.jl`).
"""

# ╔═╡ e7f5e18c-0434-4336-a552-264c58f92c24
begin
    struct RayPush
        rayPaths::Any
        tt_grid::Any
        surface_returns::Any
    end
    function Base.show(io::IO, ::MIME"text/html", p::RayPush)
        rayjs = join(map(p.rayPaths) do r
            "{z:[" * join(r.z, ",") * "],x:[" * join(r.x, ",") * "]}"
        end, ",")
        srx = join([r.x for r in p.surface_returns], ",")
        srt = join([r.t for r in p.surface_returns], ",")
        write(io, """
        <script>
        {
        const w = document.getElementById('raywidget');
        if(w){
          w.dispatchEvent(new CustomEvent('ray-update', { detail: {
            rayPaths: [$(rayjs)],
            ttFlat: [$(_ray_flatten_rowmajor(p.tt_grid))],
            surfaceReturnsX: [$srx],
            surfaceReturnsT: [$srt],
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 93b663f4-46b0-434f-a0de-a99cae6c14bf
_ray_pa = let
    v_grid = build_velocity_grid(ray_safe["bgType"], ray_safe["v0"], ray_safe["grad"], convert(Vector{Float64}, ray_safe["pert"]);
        v1=get(ray_safe, "v1", 8.0), zint=get(ray_safe, "zint", 60.0),
        lvzdepth=get(ray_safe, "lvzdepth", 60.0), lvzamp=get(ray_safe, "lvzamp", 1.5))
    s_grid, sx_grid, sz_grid = slowness_grid_and_gradients(v_grid)
    s_itp = extrapolate(interpolate((zgrid, xgrid), s_grid, Gridded(Interpolations.Linear())), Flat())
    sx_itp = extrapolate(interpolate((zgrid, xgrid), sx_grid, Gridded(Interpolations.Linear())), Flat())
    sz_itp = extrapolate(interpolate((zgrid, xgrid), sz_grid, Gridded(Interpolations.Linear())), Flat())
    (; slowness=(z, x) -> s_itp(z, x), slowness_x=(z, x) -> sx_itp(z, x), slowness_z=(z, x) -> sz_itp(z, x), xgrid, zgrid, s_grid)
end

# ╔═╡ 0d3ca438-5f67-4b75-a3aa-7ff8729ff1bf
_ray_rays = let
    aim = get(ray_safe, "angleAim", 45.0)
    n = max(1, round(Int, get(ray_safe, "nRays", 15)))
    # ray_fan's own N_RAYS==1 special case uses angleMin as the sole angle -- pass `aim` as
    # angleMin in that case so a 1-ray view shows exactly the dragged direction, not the
    # fan's edge.
    angleMin = n == 1 ? aim : aim - FAN_SPREAD_DEG / 2
    angleMax = n == 1 ? aim : aim + FAN_SPREAD_DEG / 2
    ray_fan(n, ray_safe["srcZ"], ray_safe["srcX"], angleMin, angleMax, _ray_pa)
end

# ╔═╡ 51a1345d-3e16-4bf4-879a-1da5e6ed5f2c
_ray_paths = ray_fan_paths(_ray_rays)

# ╔═╡ cddc57c7-2d46-4982-b3c5-22a417e2ce3b
_ray_tt = first_arrival_traveltimes(_ray_pa.s_grid, ray_safe["srcZ"], ray_safe["srcX"])

# ╔═╡ f5d81d25-b61a-4afb-b976-328e36f1aea2
_ray_surface_returns = ray_fan_surface_returns(_ray_rays)

# ╔═╡ 94947518-d508-4e31-837a-dc492753384b
RayPush(_ray_paths, _ray_tt, _ray_surface_returns)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Eikonal = "a6aab1ba-8f88-4217-b671-4d0788596809"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Eikonal = "~0.1.1"
Interpolations = "~0.16.2"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "884d130e167b68a435fbfc6edb6d38588b98bbd0"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "daa72978cd7a624246e894a4f4f067706d4e17e2"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "60f11b38ebeabd984f5535838d91e197d97202f0"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.28.1"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceFillArraysExt = "FillArrays"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChunkCodecCore]]
git-tree-sha1 = "1a3ad7e16a321667698a19e77362b35a1e94c544"
uuid = "0b6fb165-00bc-4d37-ab8b-79f91016dbe1"
version = "1.0.1"

[[deps.ChunkCodecLibZlib]]
deps = ["ChunkCodecCore", "Zlib_jll"]
git-tree-sha1 = "cee8104904c53d39eb94fd06cbe60cb5acde7177"
uuid = "4c0bbee4-addc-4d73-81a0-b6caacae83c8"
version = "1.0.0"

[[deps.ChunkCodecLibZstd]]
deps = ["ChunkCodecCore", "Zstd_jll"]
git-tree-sha1 = "34d9873079e4cb3d0c62926a225136824677073f"
uuid = "55437552-ac27-4d47-9aa3-63184e8fd398"
version = "1.0.0"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ef2022bff55342a8c9846cdf218f62e475f0444d"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.1.2"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "a692f5e257d332de1e554e4566a4e5a8a72de2b2"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.4"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.Eikonal]]
deps = ["DataStructures", "Images", "LinearAlgebra", "PrecompileTools", "Printf"]
git-tree-sha1 = "ac89a6cf8c89a741448deb8692aaacba745ecee0"
uuid = "a6aab1ba-8f88-4217-b671-4d0788596809"
version = "0.1.1"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6866aec60ef98e3164cd8d6855225684207e9dff"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.12+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "8e9c059d6857607253e837730dbf780b6b151acd"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.19.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7a98c6502f4632dbe9fb1973a4244eaa3324e84d"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.1"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.HistogramThresholding]]
deps = ["ImageBase", "LinearAlgebra", "MappedArrays"]
git-tree-sha1 = "7194dfbb2f8d945abdaf68fa9480a965d6661e69"
uuid = "2c695a8d-9458-5d45-9878-1b8a99cf7853"
version = "0.3.1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageBinarization]]
deps = ["HistogramThresholding", "ImageCore", "LinearAlgebra", "Polynomials", "Reexport", "Statistics"]
git-tree-sha1 = "33485b4e40d1df46c806498c73ea32dc17475c59"
uuid = "cbc4b850-ae4b-5111-9e64-df94c024a13d"
version = "0.3.1"

[[deps.ImageContrastAdjustment]]
deps = ["ImageBase", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "eb3d4365a10e3f3ecb3b115e9d12db131d28a386"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.12"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageCorners]]
deps = ["ImageCore", "ImageFiltering", "PrecompileTools", "StaticArrays", "StatsBase"]
git-tree-sha1 = "24c52de051293745a9bad7d73497708954562b79"
uuid = "89d5987c-236e-4e32-acd0-25bd6bd87b70"
version = "0.1.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "08b0e6354b21ef5dd5e49026028e41831401aca8"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.17"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "52116260a234af5f69969c5286e6a5f8dc3feab8"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.12"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils"]
git-tree-sha1 = "8e64ab2f0da7b928c8ae889c514a52741debc1c2"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.4.2"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Bzip2_jll", "FFTW_jll", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "Zstd_jll", "libpng_jll", "libwebp_jll", "libzip_jll"]
git-tree-sha1 = "d670e8e3adf0332f57054955422e85a4aec6d0b0"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "7.1.2005+0"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.ImageMorphology]]
deps = ["DataStructures", "ImageCore", "LinearAlgebra", "LoopVectorization", "OffsetArrays", "Requires", "TiledIteration"]
git-tree-sha1 = "cffa21df12f00ca1a365eb8ed107614b40e8c6da"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.4.6"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "PrecompileTools", "Statistics"]
git-tree-sha1 = "783b70725ed326340adf225be4889906c96b8fd1"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.7"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "7196039573b6f312864547eb7a74360d6c0ab8e6"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.9.0"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "dfde81fafbe5d6516fb864dc79362c5c6b973c82"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.10.2"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageBinarization", "ImageContrastAdjustment", "ImageCore", "ImageCorners", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "a49b96fd4a8d1a9a718dfd9cde34c154fc84fcd5"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.26.2"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "b842cbff3f44804a84fda409745cc8f04c029a20"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.6"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalSets]]
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLD2]]
deps = ["ChunkCodecLibZlib", "ChunkCodecLibZstd", "FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues"]
git-tree-sha1 = "941f87a0ae1b14d1ac2fa57245425b23a9d7a516"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.6.4"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1dae3057da6f2b9c857afef03177bbdc7c4afe92"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "38928f7999753af13d4e13966ae15958ff3a917a"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.19.1+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "bba2d9aa057d8f126415de240573e86a8f39d2a1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.1"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "a9fc7883eb9b5f04f46efb9a540833d1fad974b3"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.173"

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    ForwardDiffNNlibExt = ["ForwardDiff", "NNlib"]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.LoopVectorization.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "3a8f462a180a9d735e340f4e8d5f364d411da3a4"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.8.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ca7e18198a166a1f3eb92a3650d53d94ed8ca8a1"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.22"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "215a6666fee6d6b3a6e75f2cc22cb767e2dd393a"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.5+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "972089912ba299fba87671b025cd0da74f5f54f7"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.0"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "994cc27cdacca10e68feb291673ec3a76aa2fae9"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"
weakdeps = ["RecipesBase"]

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLPublic]]
git-tree-sha1 = "cf9aaf8b9ed5db993259ea8b24cf2b7ba9bd3b79"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.2.4"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "67a144433c4ce877ee6d1ada69a124d6b1ecf7be"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.6.2"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "3e5f165e58b18204aed03158664c4982d691f454"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.5.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "1e44e7b1dbb5249876d84c32466f8988a6b41bbb"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "e4d7a1a0edc20af42689ea6f4f3587a2175d50ee"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.12"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "98b9352a24cb6a2066f9ababcc6802de9aed8ad8"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.6"

[[deps.TiledIteration]]
deps = ["OffsetArrays", "StaticArrayInterface"]
git-tree-sha1 = "1176cc31e867217b06928e2f140c90bd1bc88283"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.5.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "3b0738bd7c5645641845da25cbd99800b8718689"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "d1d9a935a26c475ebffd54e9c7ad11627c43ea85"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.72"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.libzip_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "86addc139bca85fdf9e7741e10977c45785727b7"
uuid = "337d8026-41b4-5cde-a456-74a10e5b31d1"
version = "1.11.3+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "da8c1f6eee04831f14edcfa5dae611d309807e57"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.3.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═ef66e643-bc46-4a44-bdb9-d641a26f4e0e
# ╟─864f274b-4f28-4cf9-83f7-69a1f03ca870
# ╟─e6cd1e51-7666-49e9-8d8c-523abf6fe6c9
# ╟─cddc22f1-85fc-489d-85d3-d406d795e5a3
# ╟─496f5399-68e2-4a83-a546-324a9de63c09
# ╟─28f9d159-170e-4371-85c4-98566e90d331
# ╟─0f652df5-37a3-40fe-a8c9-5d873d7b491a
# ╟─8dd31a81-a418-4164-90aa-b65e48ed9b2c
# ╟─bc87032f-8941-4380-97f8-33dbebd9afe1
# ╟─75c23a27-4253-474c-97bb-925ba236f362
# ╟─56cd6049-2266-4c7a-9e5f-1749f562034a
# ╟─322e72d8-73aa-4386-af67-547ff74eec5e
# ╠═2d7556b5-27d7-40ec-b461-0945de83e5f7
# ╟─075b648c-e275-44f4-94ea-73a8958f15dd
# ╠═01e1393f-5b2a-423b-8a6a-0b8163cb1a1f
# ╠═4097ad6c-d573-456c-822e-3eb922ecf0b8
# ╠═08196372-2cc4-401f-87d7-b41d32364e5d
# ╟─f862769d-7fac-4408-80fd-aa3caf90c940
# ╠═e8d8f418-7d68-4004-b294-b48fe82c7da9
# ╟─432388ad-4852-426c-956b-cc54710cb619
# ╠═439b27c2-fd03-4391-b93c-4c58f3edd3ba
# ╠═05505f75-9720-4474-9388-250bae5f2fd2
# ╟─a9d965b7-0e29-4552-84e9-46b3bb8ace89
# ╠═6260a705-c8a8-4805-bd87-382dbc8a1d22
# ╟─450f8df0-4270-49f0-a910-c96c6bd89882
# ╠═94ad2aee-758d-4aca-8020-006ed515ebf9
# ╠═94853ecf-a2f0-4b9b-83ed-faf9869c84fd
# ╠═c497fae9-3cf2-4697-9eff-7e4da5024adf
# ╟─35124aad-70e8-42bf-a985-ebcb23d39a67
# ╠═eeaeaec0-75f0-4c57-9ce8-394d90b0ec91
# ╟─1dc9add4-83d4-4516-855b-4531208d8dee
# ╠═e7f5e18c-0434-4336-a552-264c58f92c24
# ╠═93b663f4-46b0-434f-a0de-a99cae6c14bf
# ╠═0d3ca438-5f67-4b75-a3aa-7ff8729ff1bf
# ╠═51a1345d-3e16-4bf4-879a-1da5e6ed5f2c
# ╠═cddc57c7-2d46-4982-b3c5-22a417e2ce3b
# ╠═f5d81d25-b61a-4afb-b976-328e36f1aea2
# ╠═94947518-d508-4e31-837a-dc492753384b
# ╠═05d2c8b6-41a3-445b-8988-34b5020087d5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
