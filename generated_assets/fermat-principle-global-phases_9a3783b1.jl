### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> tags = ["raytheory"]
#> title = "Fermat's Principle and Global Seismic Phases"
#> description = "Drag the real PREM velocity structure and watch shadow zones, triplications, and mode-converted phases like SKS appear."
#> layout = "layout.jlhtml"

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

# ╔═╡ 875d305b-6fa0-441f-b4ea-673e1053c669
begin
    using LinearAlgebra
    using Eikonal
    using PlutoUI
end

# ╔═╡ 95c66aa1-b555-4f38-a3a6-c79746262c87
TableOfContents()

# ╔═╡ e785f801-df99-4fb7-ad1e-34861320c9bf
md"""
# Fermat's Principle and Global Seismic Phases

A point source's wavefront is a clean expanding circle. A shadow behind the core, three separate P
arrivals from one earthquake (a triplication), an S wave that somehow overtakes direct S (SKS) —
none of that looks like "least time." This notebook is about the resolution: Fermat's principle
says a ray's travel time is *stationary*, not necessarily minimal, and every one of those
phenomena is a direct, visible consequence of that one idea applied to the real Earth.

**[ray-theory-eikonal-equation.jl](../Ray%20Theory/ray-theory-eikonal-equation.jl)** already shoots
rays through a paintable velocity model and derives this same principle — but flat and local
(500×150 km), so it can't show a shadow zone or a global triplication. **[global-seismic-
arrivals.jl](../Ray%20Theory/global-seismic-arrivals.jl)** already shows real phases by name on a
real Earth — but that Earth is fixed, so it can't show *why* those phases behave the way they do.
This notebook sits between the two: the real PREM Earth, curved and global, but with its velocity
structure in your hands.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 5be71c78-0c91-4d72-aa5d-da8763ea02c1
md"""
## The Exact PREM Model
"""

# ╔═╡ fa8a8574-165b-42a4-a14f-7ad7fb4b32f1
const R_EARTH = 6371.0 # km

# ╔═╡ 11042f31-3dbb-4a2a-89b4-598f6e72df47
"""
	load_prem_profile()

Parse the exact PREM model (Dziewonski & Anderson, 1981, *Phys. Earth Planet. Inter.* 25(4):297-356)
from its AXISEM card-deck file, `src/assets/data/specnm_models/prem_ani` — the same file
[`earth-internal-structure.jl`](../Earth%20Structure%20and%20Internal%20Layers/earth-internal-structure.jl)
and `specnm.jl` use. Unlike that notebook's own parser, real discontinuities (two consecutive rows
sharing one radius, each a genuine jump in Vp/Vs — including `vs=0` across the whole fluid outer
core) are kept exactly, not nudged apart: [`trace_ray`](@ref) below needs the true jump, not a
smoothed one.

Returns `(r, vp, vs)`, three vectors of equal length ordered by strictly non-increasing radius `r`
(km), exactly as the file lists them (surface to center), with a discontinuity appearing as two
consecutive equal `r` entries.
"""
function load_prem_profile()
    r_m = Float64[]
    vp = Float64[]
    vs = Float64[]
    path = joinpath(@__DIR__, "..", "assets", "data", "specnm_models", "prem_ani")
    open(path) do io
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, "#") || !occursin(r"^[0-9.]", s)) && continue
            cols = split(s)
            push!(r_m, parse(Float64, cols[1]))
            push!(vp, parse(Float64, cols[3]))
            push!(vs, parse(Float64, cols[4]))
        end
    end
    return r_m ./ 1000.0, vp ./ 1000.0, vs ./ 1000.0 # km, km/s
end

# ╔═╡ 6c0c8399-320f-4381-a007-dfccc0f7eddd
"""
	prem_discontinuities(r, vp, vs)

Every genuine jump in a radial profile `(r, vp, vs)` (as returned by [`load_prem_profile`](@ref)
or a student-edited copy of it): a place where two consecutive rows share the same radius. Returns
a vector of named tuples `(i, r, depth, vp_outer, vp_inner, vs_outer, vs_inner, label)` — `i` the
array index of the outer row (so `r[i]==r[i+1]`, letting a caller edit that jump in place), `outer`
meaning the value approached from the surface side, `inner` the value continuing toward the
center. `label` names the five best-known discontinuities (Moho, 400, 670, CMB, ICB — the same
depths [`global-seismic-arrivals.jl`](../Ray%20Theory/global-seismic-arrivals.jl)'s `DISCS` uses)
when the jump's depth is within 5 km of one of them, else a plain `"<depth> km"` label.
"""
function prem_discontinuities(r, vp, vs)
    named = [(24.4, "Moho"), (400.0, "400"), (670.0, "670"), (2891.0, "CMB"), (5149.5, "ICB")]
    discs = NamedTuple[]
    for i in 1:length(r)-1
        if r[i] == r[i+1]
            depth = R_EARTH - r[i]
            hit = findfirst(nd -> abs(nd[1] - depth) < 5.0, named)
            label = isnothing(hit) ? "$(round(depth, digits=1)) km" : named[hit][2]
            push!(discs, (i=i, r=r[i], depth=depth, vp_outer=vp[i], vp_inner=vp[i+1],
                vs_outer=vs[i], vs_inner=vs[i+1], label=label))
        end
    end
    return discs
end

# ╔═╡ 87825c53-c340-47a1-9544-7be10960fe8f
prem_r, prem_vp, prem_vs = load_prem_profile();

# ╔═╡ 738f4b56-13d1-459f-a8ac-4b22aab71c8a
md"""
## Circular Ray Shooting

A fast-sweep Eikonal solve (as used in the flat notebook) gives only the single first-arrival time
everywhere — it cannot show a genuine shadow zone (a range of angular distance nothing lands in) or
a triplication (several rays landing at the same distance). Shooting a fan of rays at many takeoff
angles is the only way to show these, and it's the direct computational embodiment of "trying every
path" that Fermat's principle is about.

`get_raypaths` in the flat notebook bends a ray using a grid-based finite-difference gradient of
slowness, which is fine for a smooth background but would produce a spurious numerical kick across
one of PREM's real jump discontinuities. Since the background here is exactly radial, `dv/dr` is
taken analytically from the profile's own piecewise-linear segments instead, and every discontinuity
crossing is handled as an explicit Snell's-law event.
"""

# ╔═╡ 6e3900a8-b1f4-458e-b7cf-8940ab465279
"""
	segment_velocity(r, v, r0)

Linearly interpolate `v(r0)` and its exact local slope `dv/dr` from a piecewise-linear profile
`(r, v)`, for `r0` strictly inside one smooth segment. [`trace_ray`](@ref) never calls this exactly
on a discontinuity radius — those are handled as explicit events instead, using the discontinuity's
own `_outer`/`_inner` values, so there is no ambiguity about which side is meant.
"""
function segment_velocity(r, v, r0)
    r0 = clamp(r0, r[end], r[1])
    n = length(r)
    for i in 1:n-1
        if r[i] >= r0 >= r[i+1] && r[i] != r[i+1]
            t = (r0 - r[i+1]) / (r[i] - r[i+1])
            return v[i+1] + t * (v[i] - v[i+1]), (v[i] - v[i+1]) / (r[i] - r[i+1])
        end
    end
    return v[end], 0.0
end

# ╔═╡ ccbed1cc-ce29-41b8-8bb1-d2f203735696
"""
	circle_crossing_fraction(x, z, dx, dz, ds, r_target)

Exact fractional step `α ∈ (0, ds]` at which the straight segment from `(x,z)` in unit direction
`(dx,dz)` first reaches radius `r_target` — the ray-circle intersection used to land a step exactly
on the Earth's surface or a discontinuity radius instead of stepping past it. Solving
`‖(x,z)+α(dx,dz)‖ = r_target` for `α` is an exact quadratic (the trial step is a straight line),
so this needs no iteration. Returns `nothing` if the segment never reaches `r_target` within
`(0, ds]`.
"""
function circle_crossing_fraction(x, z, dx, dz, ds, r_target)
    b = x * dx + z * dz
    c = x^2 + z^2 - r_target^2
    disc = b^2 - c
    disc < 0 && return nothing
    sq = sqrt(disc)
    for alpha in (-b - sq, -b + sq)
        1e-9 < alpha <= ds && return alpha
    end
    return nothing
end

# ╔═╡ 47f4eb99-0684-49cf-89eb-7b63c9d4407a
"""
	refract_or_reflect(dx, dz, x, z, v_from, v_to)

Apply Snell's law at a spherical interface crossed by a ray moving in unit direction `(dx, dz)` at
position `(x, z)` (radius `r = hypot(x,z)`), where the wave speed changes from `v_from` to `v_to`.
The ray parameter `p = r sin i / v` (`i` from the local radial direction) is conserved, so the
transmitted angle follows directly: `sin i₂ = (r sin i₁ / v_from) · v_to / r`. Returns the new unit
direction and whether the ray transmitted (`true`) or underwent total internal reflection (`false`,
direction reflected, medium/velocity unchanged).
"""
function refract_or_reflect(dx, dz, x, z, v_from, v_to)
    r = hypot(x, z)
    nx, nz = x / r, z / r
    tx, tz = -nz, nx
    d_r = dx * nx + dz * nz
    d_t = dx * tx + dz * tz
    sin_i1 = clamp(d_t, -1.0, 1.0)
    sin_i2 = (sin_i1 / v_from) * v_to
    if abs(sin_i2) <= 1.0
        cos_i2 = (d_r < 0 ? -1.0 : 1.0) * sqrt(max(0.0, 1.0 - sin_i2^2))
        return sin_i2 * tx + cos_i2 * nx, sin_i2 * tz + cos_i2 * nz, true
    else
        return d_t * tx - d_r * nx, d_t * tz - d_r * nz, false
    end
end

# ╔═╡ 5ec8e794-8b2e-4355-820a-db771c8cdece
"""
	trace_ray(profile, discs, wavetype0, takeoff_deg; ds=2.0, max_steps=8000, src_r=R_EARTH,
	          src_theta_deg=0.0, convert_on_core_exit=false, max_bounces=0)

Shoot one ray from `(src_r, src_theta_deg)` at `takeoff_deg` from the local downward radial
direction, through the radial velocity structure `profile = (r=, vp=, vs=)` (as returned by
[`load_prem_profile`](@ref) or the widget's edited copy), starting as wave type `wavetype0`
(`:P`/`:S`). `discs` is the discontinuity list from [`prem_discontinuities`](@ref), computed once
per shot and shared across a whole [`ray_fan`](@ref).

Integrates the ray equations analytically between discontinuities and treats each crossing as an
explicit Snell's-law event conserving `p = r sin i / v(r)`. Because real PREM has `vs ≈ 0` in the
fluid outer core, an `:S` ray meeting the CMB is automatically converted to `:P` to cross — this
falls straight out of the same event handling once "S cannot transmit into `vs≈0`" is checked,
it isn't a hand-coded special case. If `convert_on_core_exit`, a `:P` ray leaving a `vs≈0` layer
outward is converted back to `:S` (builds the classic SKS phase instead of PKS).

If `max_bounces > 0`, a ray that reaches the free surface `r=R_EARTH` before its budget of bounces
is used up specularly reflects (radial component flips, tangential component and wave type are
preserved — the same approximation already used for total-internal-reflection at an internal
discontinuity, just applied at the outer boundary) and keeps tracing, instead of stopping there.
This is what turns a single-leg `P`/`S` shot into a multi-leg `PP`/`PPP`-style path; only the
*final* surface arrival (once the bounce budget is exhausted) counts as `landed`.

Returns `(x, z, wavetype, t, landed, delta_deg, p0, nbounces)`: the traced path (km, km, `:P`/`:S`
per point, s); `landed` is `true` only if the ray returned to `r=R_EARTH` with no bounces left to
spend, in which case `delta_deg` is its total angular distance from the source (`nothing`
otherwise — exactly how the shadow zone shows up); `p0` is the ray's initial parameter; `nbounces`
is how many free-surface reflections it actually used.
"""
function trace_ray(profile, discs, wavetype0, takeoff_deg; ds=2.0, max_steps=8000,
    src_r=R_EARTH, src_theta_deg=0.0, convert_on_core_exit=false, max_bounces=0)
    vel_slope(wt, r) = wt == :P ? segment_velocity(profile.r, profile.vp, r) :
                        segment_velocity(profile.r, profile.vs, r)

    src_theta = deg2rad(src_theta_deg)
    x, z = src_r * sin(src_theta), src_r * cos(src_theta)
    nx0, nz0 = -x / src_r, -z / src_r # inward radial
    tx0, tz0 = -nz0, nx0
    θ = deg2rad(takeoff_deg)
    dxu, dzu = cos(θ) * nx0 + sin(θ) * tx0, cos(θ) * nz0 + sin(θ) * tz0

    wavetype = wavetype0
    v0, = vel_slope(wavetype, src_r)
    p0 = src_r * sin(θ) / v0
    S = (dxu / v0, dzu / v0)

    xs, zs, wts, ts = [x], [z], Symbol[wavetype], [0.0]
    t = 0.0
    landed = false
    delta_deg = nothing
    nbounces = 0

    for _ in 1:max_steps
        r = hypot(x, z)
        v, dvdr = vel_slope(wavetype, r)
        s = 1.0 / v
        dxu, dzu = S[1] * v, S[2] * v

        alpha_surf = circle_crossing_fraction(x, z, dxu, dzu, ds, R_EARTH)
        best_alpha = alpha_surf === nothing ? Inf : alpha_surf
        best_disc = nothing
        for d in discs
            a = circle_crossing_fraction(x, z, dxu, dzu, ds, d.r)
            if a !== nothing && a < best_alpha
                best_alpha, best_disc = a, d
            end
        end

        if best_alpha <= ds
            x, z = x + dxu * best_alpha, z + dzu * best_alpha
            t += s * best_alpha
            push!(xs, x); push!(zs, z); push!(ts, t)
            if best_disc === nothing
                if nbounces < max_bounces
                    # Specular free-surface reflection: flip the radial component of the unit
                    # direction, keep the tangential one and the wave type (same approximation as
                    # a total-internal-reflection event above, just at the outer boundary instead
                    # of an internal one).
                    r_surf = hypot(x, z)
                    nxr, nzr = x / r_surf, z / r_surf
                    d_r = dxu * nxr + dzu * nzr
                    dx2, dz2 = dxu - 2 * d_r * nxr, dzu - 2 * d_r * nzr
                    v_here, = vel_slope(wavetype, r_surf)
                    S = (dx2 / v_here, dz2 / v_here)
                    nbounces += 1
                    # Same nudge-past-the-boundary trick used for interface crossings below, so
                    # the very next step doesn't immediately re-detect this same surface point.
                    x += dx2 * 1e-3
                    z += dz2 * 1e-3
                    t += 1e-3 / v_here
                    push!(xs, x); push!(zs, z); push!(ts, t); push!(wts, wavetype)
                    continue
                end
                landed = true
                push!(wts, wavetype)
                delta_deg = mod(rad2deg(atan(x, z)) - src_theta_deg, 360.0)
                break
            end
            going_inward = (x * dxu + z * dzu) < 0
            v_from = wavetype == :P ? (going_inward ? best_disc.vp_outer : best_disc.vp_inner) :
                     (going_inward ? best_disc.vs_outer : best_disc.vs_inner)
            if wavetype == :S
                vs_far = going_inward ? best_disc.vs_inner : best_disc.vs_outer
                v_to, new_type = vs_far < 1e-6 ?
                                  ((going_inward ? best_disc.vp_inner : best_disc.vp_outer), :P) :
                                  (vs_far, :S)
            else
                vs_from = going_inward ? best_disc.vs_outer : best_disc.vs_inner
                v_to, new_type = (!going_inward && convert_on_core_exit && vs_from < 1e-6) ?
                                  (best_disc.vs_outer, :S) :
                                  ((going_inward ? best_disc.vp_inner : best_disc.vp_outer), :P)
            end
            dx2, dz2, transmitted = refract_or_reflect(dxu, dzu, x, z, v_from, v_to)
            if !transmitted
                v_to, new_type = v_from, wavetype
                dx2, dz2 = dxu, dzu # refract_or_reflect already returned the reflected direction here
            end
            wavetype = new_type
            S = (dx2 / v_to, dz2 / v_to)
            # Nudge a hair past the interface: `x,z` sit at r == the disc radius up to ~1e-7 km
            # of floating-point slop from the circle-intersection solve, and `segment_velocity`
            # can resolve that ambiguous boundary to the side just *left* -- catastrophic if that
            # side is the fluid core (v≈0), producing s=1/v=Inf and a NaN path from then on. A
            # 1 m step in the new direction (1e-3 km) is ~10,000x that slop, so it decisively
            # lands strictly inside the new medium, and is far below any visible/physical scale.
            x += dx2 * 1e-3
            z += dz2 * 1e-3
            t += 1e-3 / v_to
            push!(wts, wavetype)
        else
            dsdx = -dvdr / v^2 * (x / r)
            dsdz = -dvdr / v^2 * (z / r)
            x += dxu * ds
            z += dzu * ds
            S = (S[1] + dsdx * ds, S[2] + dsdz * ds)
            t += s * ds
            push!(xs, x); push!(zs, z); push!(wts, wavetype); push!(ts, t)
        end
    end

    return (x=xs, z=zs, wavetype=wts, t=ts, landed=landed, delta_deg=delta_deg, p0=p0, nbounces=nbounces)
end

# ╔═╡ c48a015d-66dd-407a-92c9-ddc2088612b3
"""
	ray_fan(profile, discs, N, wavetype0, angle_min, angle_max; kwargs...)

Shoot `N` rays evenly spaced in takeoff angle between `angle_min` and `angle_max` (degrees from the
local downward radial direction), forwarding `kwargs` to [`trace_ray`](@ref). This is the single
source of truth the widget draws from and the travel-time curve is built from.
"""
function ray_fan(profile, discs, N, wavetype0, angle_min, angle_max; kwargs...)
    angles = N == 1 ? [angle_min] : range(angle_min, angle_max; length=N)
    return [trace_ray(profile, discs, wavetype0, θ; kwargs...) for θ in angles]
end

# ╔═╡ 83953fdc-7750-4da6-97a8-00d08b3c7f7d
md"""
## Eikonal Cross-Check

A shot ray's own accumulated travel time can be checked against a completely independent
numerical method: `Eikonal.jl`'s fast-sweeping solver (the same package and call pattern the
flat notebook uses for its own first-arrival wavefront), run once on a Cartesian sampling of the
disk. This has nothing to do with the ray-shooting physics above — it's purely a second, unrelated
way of computing a travel time, used only to confirm the two agree.
"""

# ╔═╡ 3b85d627-6ba9-4bda-b571-489a774e1e81
"""
	eikonal_crosscheck_time(profile, wavetype, src_r, src_theta_deg, x, z; n=301)

Independently solve the Eikonal equation with `Eikonal.jl`'s fast-sweeping method (same call
pattern as `first_arrival_traveltimes` in the flat notebook) on an `n×n` Cartesian grid covering
the Earth disk, for `wavetype`'s velocity field, and return the first-arrival travel time at
`(x, z)`. Used only by the self-check below, to cross-check a [`trace_ray`](@ref) ray's own
accumulated time against a method that shares none of its code.
"""
function eikonal_crosscheck_time(profile, wavetype, src_r, src_theta_deg, x, z; n=301)
    xs = range(-R_EARTH, R_EARTH; length=n)
    zs = range(-R_EARTH, R_EARTH; length=n)
    dcell = xs[2] - xs[1]
    v_at(xx, zz) = wavetype == :P ? segment_velocity(profile.r, profile.vp, hypot(xx, zz))[1] :
                   segment_velocity(profile.r, profile.vs, hypot(xx, zz))[1]
    # Outside the disk is given a tiny velocity (a near-impassable "wall"), not a fast one --
    # a fast exterior would let the fast-sweep shortcut through it and hug the boundary,
    # systematically underestimating travel time for any point near the surface.
    vgrid = [hypot(xx, zz) <= R_EARTH ? v_at(xx, zz) : 0.01 for zz in zs, xx in xs]
    sgrid = 1.0 ./ vgrid
    src_x, src_z = src_r * sin(deg2rad(src_theta_deg)), src_r * cos(deg2rad(src_theta_deg))
    iz0, ix0 = argmin(abs.(zs .- src_z)), argmin(abs.(xs .- src_x))
    fastsweep = FastSweeping(sgrid)
    init!(fastsweep, (iz0, ix0))
    sweep!(fastsweep, verbose=false, epsilon=1e-6)
    tgrid = (fastsweep.t .* dcell)[1:n, 1:n]
    iz, ix = argmin(abs.(zs .- z)), argmin(abs.(xs .- x))
    return tgrid[iz, ix]
end

# ╔═╡ 199e63bd-a95e-4fb9-bc26-1ca0c05dd12e
md"""
### Verifying the Ray Tracer
"""

# ╔═╡ 86bf927a-b3c9-4edb-888d-8ca18afa299a
let
    r_flat = [R_EARTH, 0.0]
    v_flat = [6.0, 6.0]
    profile = (r=r_flat, vp=v_flat, vs=v_flat ./ sqrt(3))
    discs = prem_discontinuities(profile.r, profile.vp, profile.vs)
    @assert isempty(discs)
    ray = trace_ray(profile, discs, :P, 30.0)
    @assert ray.landed
    i = length(ray.x) ÷ 2
    predicted_t = hypot(ray.x[i], ray.z[i] - R_EARTH) / 6.0
    @assert isapprox(ray.t[i], predicted_t; rtol=1e-3)
    "Homogeneous-medium straight-ray check passed"
end

# ╔═╡ 35793516-2618-4552-8dee-cb23a59296de
"""
	local_ray_parameter(profile, ray, i)

Estimate the ray parameter `p = r sin i / v` at point `i` of a traced `ray` (from
[`trace_ray`](@ref)) directly from the path geometry — used only by the self-checks below to
verify `p` independently of the value `trace_ray` itself reports as `p0`.
"""
function local_ray_parameter(profile, ray, i)
    r = hypot(ray.x[i], ray.z[i])
    dx, dz = ray.x[i+1] - ray.x[i], ray.z[i+1] - ray.z[i]
    len = hypot(dx, dz)
    nx, nz = ray.x[i] / r, ray.z[i] / r
    tx, tz = -nz, nx
    sinI = (dx * tx + dz * tz) / len
    v = ray.wavetype[i] == :P ? segment_velocity(profile.r, profile.vp, r)[1] :
        segment_velocity(profile.r, profile.vs, r)[1]
    return r * abs(sinI) / v
end

# ╔═╡ 0dd771df-6970-4eec-a26a-975319f0d697
let
    profile = (r=prem_r, vp=prem_vp, vs=prem_vs)
    discs = prem_discontinuities(prem_r, prem_vp, prem_vs)
    ray = trace_ray(profile, discs, :P, 25.0)
    idxs = round.(Int, range(3, length(ray.x) - 3; length=8))
    ps = [local_ray_parameter(profile, ray, i) for i in idxs]
    @assert all(p -> isapprox(p, ray.p0; rtol=2e-2), ps) "measured p = $ps vs p0 = $(ray.p0)"
    "Ray-parameter conservation through real PREM passed (p ≈ $(round(ray.p0, digits=4)))"
end

# ╔═╡ 4d775817-c78d-4e2a-8d5f-f951daba5ce3
let
    profile = (r=prem_r, vp=prem_vp, vs=prem_vs)
    discs = prem_discontinuities(prem_r, prem_vp, prem_vs)
    rays = ray_fan(profile, discs, 60, :P, 3.0, 80.0)
    landed = filter(r -> r.landed, rays)
    order = sortperm([r.delta_deg for r in landed])
    landed = landed[order]
    deltas = [r.delta_deg for r in landed]
    times = [r.t[end] for r in landed]
    ps = [r.p0 for r in landed]
    # Δ≈50° sits on the smooth single-branch mantle-P curve, clear of both the upper-mantle
    # (400/670) triplication range (roughly 15-25° here) and the core shadow zone (>~100°) --
    # see the notebook text for how those ranges were located.
    i = argmin(abs.(deltas .- 50.0))
    i = clamp(i, 2, length(deltas) - 1)
    dTdDelta = (times[i+1] - times[i-1]) / deg2rad(deltas[i+1] - deltas[i-1])
    @assert isapprox(abs(dTdDelta), ps[i]; rtol=0.05) "dT/dΔ=$dTdDelta vs p=$(ps[i]) at Δ=$(deltas[i])°"
    "dT/dΔ ≈ p verified numerically (dT/dΔ = $(round(dTdDelta, digits=3)), p = $(round(ps[i], digits=3)) at Δ ≈ $(round(deltas[i], digits=1))°)"
end

# ╔═╡ ec425e80-1077-4f29-ba4a-befc34a94ea4
let
    profile = (r=prem_r, vp=prem_vp, vs=prem_vs)
    discs = prem_discontinuities(prem_r, prem_vp, prem_vs)
    ray = trace_ray(profile, discs, :P, 30.0)
    @assert ray.landed
    # The turning point (minimum radius) sits well below PREM's steep near-surface layer --
    # the fairest place to compare against a coarse Cartesian fast-sweep grid, whose cell size
    # (tens of km) cannot resolve that thin layer well.
    i = argmin(hypot.(ray.x, ray.z))
    t_eikonal = eikonal_crosscheck_time(profile, :P, R_EARTH, 0.0, ray.x[i], ray.z[i]; n=401)
    @assert isapprox(ray.t[i], t_eikonal; rtol=0.05) "trace_ray t=$(ray.t[i]) vs Eikonal.jl t=$t_eikonal"
    "Eikonal.jl cross-check passed (trace_ray t=$(round(ray.t[i], digits=2))s vs fast-sweep t=$(round(t_eikonal, digits=2))s)"
end

# ╔═╡ b97f0fa4-eb4d-425b-bdc4-b63a983ca8ed
let
    profile = (r=prem_r, vp=prem_vp, vs=prem_vs)
    discs = prem_discontinuities(prem_r, prem_vp, prem_vs)
    cmb = only(filter(d -> d.label == "CMB", discs))
    ray = trace_ray(profile, discs, :S, 8.0; convert_on_core_exit=true)
    switch_in = findfirst(i -> ray.wavetype[i] == :S && ray.wavetype[i+1] == :P, 1:length(ray.wavetype)-1)
    @assert switch_in !== nothing "ray never reached the fluid core -- steepen the takeoff angle"
    r_at = hypot(ray.x[switch_in], ray.z[switch_in])
    @assert isapprox(r_at, cmb.r; atol=5.0)
    p_before = local_ray_parameter(profile, ray, switch_in - 1)
    p_after = local_ray_parameter(profile, ray, switch_in + 1)
    @assert isapprox(p_before, p_after; rtol=3e-2) "p_before=$p_before vs p_after=$p_after across S→P conversion"
    "Ray parameter survives the S→P mode conversion at the CMB (SKS) ✓"
end

# ╔═╡ c2c98c49-6b83-4907-8e66-a769d35e137d
let
    profile = (r=prem_r, vp=prem_vp, vs=prem_vs)
    discs = prem_discontinuities(prem_r, prem_vp, prem_vs)
    # A full SKS round trip: S down through the mantle, forced S->P at the CMB, P->S back on
    # exit. Position sits exactly on the CMB radius (up to circle-intersection floating slop)
    # right after each conversion, and `segment_velocity` resolving that ambiguous boundary to
    # the just-exited (fluid, v≈0) side would blow the path up to NaN -- exactly what happened
    # before `trace_ray` nudged a hair past the interface after every crossing.
    ray = trace_ray(profile, discs, :S, 8.0; convert_on_core_exit=true)
    @assert !any(isnan, ray.x) && !any(isnan, ray.z) "path went NaN -- interface nudge regressed"
    @assert ray.landed "SKS ray failed to return to the surface"
    "Full SKS round trip (S→P→S) returns to the surface cleanly at Δ ≈ $(round(ray.delta_deg, digits=1))°"
end

# ╔═╡ 05a9c088-4976-4034-89df-ddbcdfe46b00
md"""
## The Interactive Widget

Every drag happens directly on the velocity-model disk itself:

- **Dashed rings** mark the five real PREM discontinuities (Moho, 400, 670, CMB, ICB). Drag one to
  change Vp or Vs right at that jump — outside the ring is the shallow (outer) side, inside is the
  deep (inner) side.
- **Dotted rings** mark three smooth, thick layers (upper mantle, lower mantle, outer core). Drag
  one to reshape the *gradient* inside that layer — the discontinuities bounding it stay exactly
  fixed, by construction (see `tent` in the Appendix), so you can steepen or flatten a stretch of
  mantle without ever touching a real jump.
- **The yellow star** is the source, fixed at the surface. Drag from it to *aim*: a faint ghost fan
  previews the spread you're about to shoot, centered on a bright arrow that tracks your drag.
  Release, and the real (curved) rays shoot along it — a bright head races along each one's actual
  computed path, leaving a trail, so what you see land is exactly what `trace_ray` found, not the
  straight ghost lines.

Whichever field the P/S toggle currently has on screen is the one a ring drag edits, so what's
plotted is what's editable. Nothing recomputes until you release. This is enough to reproduce every
phenomenon below — drag the CMB ring together to close the shadow zone, drag the 670 ring apart to
sharpen the triplication, pull the lower-mantle bump down to exaggerate the mantle's continuous
refraction, flatten everything to see straight chords.
"""

# ╔═╡ 4a03c3a7-ec91-4f9a-9b2a-d8c34c7e2af5
prem_named_discs = filter(d -> d.label in ("Moho", "400", "670", "CMB", "ICB"),
    prem_discontinuities(prem_r, prem_vp, prem_vs));

# ╔═╡ d46fabe3-57ff-4d98-8eb6-266da1a58c2f
"""
	tent(depth, lo, center, hi)

Triangular bump: `0` at and beyond `lo`/`hi` (km depth), rising linearly to `1` at `center`.
Multiplying a velocity offset by this and adding it to a smooth stretch of the profile reshapes
that stretch while leaving its value at `lo` and `hi` **exactly** unchanged — used to let the
smooth interior of a layer (not one of PREM's real discontinuities) be dragged without disturbing
the discontinuities bounding it.
"""
function tent(depth, lo, center, hi)
    if depth <= lo || depth >= hi
        0.0
    elseif depth <= center
        (depth - lo) / (center - lo)
    else
        (hi - depth) / (hi - center)
    end
end

# ╔═╡ dd31d8ea-11e5-4f84-8286-e5ef65c5ca96
"""
	bump_specs

Three smooth, thick layers picked for interior reshaping -- upper mantle (between the Moho and
400 discontinuities), lower mantle (between 670 and the CMB), and outer core (between the CMB
and ICB) -- each a `(lo, center, hi, label)` named tuple for [`tent`](@ref), with `lo`/`hi` taken
directly from the real discontinuity depths in [`prem_named_discs`](@ref) so a bump can never
leak past the boundary that must stay fixed.
"""
bump_specs = let
    byname(l) = only(filter(d -> d.label == l, prem_named_discs))
    [
        (lo=byname("Moho").depth, hi=byname("400").depth, label="Upper Mantle"),
        (lo=byname("670").depth, hi=byname("CMB").depth, label="Lower Mantle"),
        (lo=byname("CMB").depth, hi=byname("ICB").depth, label="Outer Core"),
    ] .|> spec -> (lo=spec.lo, center=(spec.lo + spec.hi) / 2, hi=spec.hi, label=spec.label)
end;

# ╔═╡ a8fa560a-01ee-4d42-89e6-5feaab505c54
begin
    """A circular Earth cross-section built from the exact PREM model, with its Vp and Vs directly
    draggable at the five real discontinuities (Moho, 400, 670, CMB, ICB) and at three smooth
    interior layers (see `bump_specs`) -- drag a ring to reshape that wave type's structure;
    everything else keeps its exact PREM value. The source sits fixed at the surface; drag from it
    to aim a fan of rays (center angle `aimAngle`, spread `fanWidth`) and shows the resulting ray
    family and travel-time-vs-distance curve."""
    struct FermatGlobalPhasesInput
        discVpOuter::Vector{Float64}
        discVpInner::Vector{Float64}
        discVsOuter::Vector{Float64}
        discVsInner::Vector{Float64}
        bumpVp::Vector{Float64}
        bumpVs::Vector{Float64}
        wavetype0::String
        convertOnExit::Bool
        nRays::Int
        aimAngle::Float64  # takeoff angle of the fan's center, dragged from the (fixed) source
        fanWidth::Float64  # total angular spread of the fan around aimAngle
        maxBounces::Int    # free-surface reflections allowed before a ray must land (0 = direct)
    end

    FermatGlobalPhasesInput(;
        discVpOuter=[d.vp_outer for d in prem_named_discs],
        discVpInner=[d.vp_inner for d in prem_named_discs],
        discVsOuter=[d.vs_outer for d in prem_named_discs],
        discVsInner=[d.vs_inner for d in prem_named_discs],
        bumpVp=zeros(length(bump_specs)), bumpVs=zeros(length(bump_specs)),
        wavetype0="P", convertOnExit=false, nRays=41, aimAngle=43.0, fanWidth=84.0, maxBounces=0) =
        FermatGlobalPhasesInput(Float64.(discVpOuter), Float64.(discVpInner), Float64.(discVsOuter),
            Float64.(discVsInner), Float64.(bumpVp), Float64.(bumpVs), wavetype0,
            convertOnExit, Int(nRays), Float64(aimAngle), Float64(fanWidth), Int(maxBounces))

    Base.get(w::FermatGlobalPhasesInput) = Dict{String,Any}(
        "discVpOuter" => w.discVpOuter, "discVpInner" => w.discVpInner,
        "discVsOuter" => w.discVsOuter, "discVsInner" => w.discVsInner,
        "bumpVp" => w.bumpVp, "bumpVs" => w.bumpVs,
        "wavetype0" => w.wavetype0, "convertOnExit" => w.convertOnExit,
        "nRays" => w.nRays, "aimAngle" => w.aimAngle, "fanWidth" => w.fanWidth,
        "maxBounces" => w.maxBounces)

    function Base.show(io::IO, ::MIME"text/html", w::FermatGlobalPhasesInput)
        disc_depths = [d.depth for d in prem_named_discs]
        disc_labels = [d.label for d in prem_named_discs]
        bump_depths = [b.center for b in bump_specs]
        bump_labels = [b.label for b in bump_specs]
        write(io, """
        <div id="fermatwidget">
        <style>
        #fermatwidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #fermatwidget .fw-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #fermatwidget .fw-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #fermatwidget .fw-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #fermatwidget .fw-workspace{display:flex;flex-wrap:nowrap;gap:16px;justify-content:center;align-items:flex-start}
        #fermatwidget .fw-col{display:flex;flex-direction:column;align-items:center}
        #fermatwidget .fw-panel-title{font-size:14px;font-weight:700;margin-bottom:4px}
        #fermatwidget .fw-panel{background:#000;border:1px solid #374151;border-radius:6px}
        #fermatwidget .fw-caption{font-size:12px;color:#9ca3af;text-align:center;margin-top:4px;max-width:340px}
        #fermatwidget canvas{display:block}
        #fermatwidget .fw-controls{margin-top:14px;display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));
          gap:10px;width:100%;box-sizing:border-box}
        #fermatwidget .fw-control-group{background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
        #fermatwidget .fw-control-title{font-size:14px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #fermatwidget .fw-control-row{display:grid;grid-template-columns:64px minmax(0,1fr) 54px;gap:6px;align-items:center;margin:5px 0}
        #fermatwidget .fw-control-row label{font-size:12px;color:#9ca3af}
        #fermatwidget .fw-control-row input[type=range]{width:100%;min-width:0}
        #fermatwidget .fw-value{font-size:12px;color:#e5e7eb;text-align:right}
        #fermatwidget select{background:#0b0b0b;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:3px;width:100%}
        #fermatwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:5px 10px;font-size:12px;cursor:pointer}
        #fermatwidget button.active{background:#2563eb;border-color:#93c5fd}
        #fermatwidget .fw-actions{display:flex;gap:8px;flex-wrap:wrap;align-items:center}
        #fermatwidget .fw-legend{font-size:12px;color:#9ca3af;margin-top:4px}
        #fermatwidget .fw-legend b{color:#e5e7eb}
        #fermatwidget .fw-view-controls{display:flex;align-items:center;gap:6px;margin-top:6px}
        #fermatwidget .fw-zoom-level{min-width:3.2rem;color:#d1d5db;font-size:12px;text-align:center}
        </style>

        <div class="fw-title">
          <div class="fw-title-desc">Drag directly on the velocity model to reshape the real PREM structure and aim the ray fan, and watch shadow zones, triplications, and mode-converted phases appear.</div>
          <div class="fw-title-hint">drag a dashed ring = edit a real discontinuity &middot; drag a dotted ring = reshape a smooth layer, boundaries stay fixed &middot; drag the star = aim the fan &middot; editing Vp or Vs follows the P/S toggle below &middot; release to shoot &middot; zoom and drag empty space to pan</div>
        </div>

        <div class="fw-workspace">
          <div class="fw-col">
            <div class="fw-panel-title">Velocity Model + Ray Fan (drag a ring to edit)</div>
            <div class="fw-panel"><canvas id="fw-model"></canvas></div>
            <div class="fw-caption" id="fw-model-caption"></div>
            <div class="fw-view-controls" aria-label="Plot zoom controls">
              <button id="fw-zoomout" type="button" aria-label="Zoom out">&minus;</button>
              <span id="fw-zoomlevel" class="fw-zoom-level" aria-live="polite">100%</span>
              <button id="fw-zoomin" type="button" aria-label="Zoom in">+</button>
              <button id="fw-zoomreset" type="button">Reset view</button>
            </div>
            <canvas id="fw-colorbar" style="display:block;margin-top:4px"></canvas>
          </div>
          <div class="fw-col">
            <div class="fw-panel-title">Travel-Time Curve T(Δ)</div>
            <div class="fw-panel"><canvas id="fw-tt"></canvas></div>
            <div class="fw-caption" id="fw-tt-caption"></div>
          </div>
        </div>

        <div class="fw-controls">
          <div class="fw-control-group">
            <div class="fw-control-title">Source</div>
            <div class="fw-control-row"><label>Wave type</label>
              <div class="fw-actions"><button id="fw-wtp" type="button">P</button><button id="fw-wts" type="button">S</button></div>
              <span></span></div>
            <div style="font-size:12px;color:#9ca3af;margin:4px 0">the source is fixed at the surface -- drag from the star to aim the fan</div>
            <label style="font-size:12px;color:#9ca3af;display:flex;gap:6px;align-items:center;margin-top:4px">
              <input id="fw-convert" type="checkbox"> convert back to S on leaving the core (SKS, not PKS)
            </label>
          </div>
          <div class="fw-control-group">
            <div class="fw-control-title">Ray Fan</div>
            <div class="fw-control-row"><label>N rays</label><input id="fw-nrays" type="range" min="5" max="120" step="1"><span class="fw-value" id="fw-nrays-v"></span></div>
            <div class="fw-control-row"><label>Fan width</label><input id="fw-fanwidth" type="range" min="2" max="90" step="1"><span class="fw-value" id="fw-fanwidth-v"></span></div>
            <div class="fw-control-row"><label>Bounces</label>
              <select id="fw-bounces">
                <option value="0">0 (direct P/S)</option>
                <option value="1">1 (PP-style)</option>
                <option value="2">2</option>
                <option value="3">3</option>
              </select>
              <span></span>
            </div>
            <div style="font-size:12px;color:#9ca3af;margin-top:4px" id="fw-aim-readout"></div>
          </div>
          <div class="fw-control-group">
            <div class="fw-control-title">Scenario Presets</div>
            <select id="fw-preset">
              <option value="prem">Exact PREM</option>
              <option value="uniform">Uniform Earth (no mirage)</option>
              <option value="nocore">No core contrast (shadow zone closes)</option>
              <option value="sharp660">Sharper 670 (triplication)</option>
              <option value="sks">SKS setup</option>
            </select>
            <div class="fw-actions" style="margin-top:8px"><button id="fw-reset" type="button">Reset to PREM</button></div>
            <div class="fw-legend"><b>Ray colors:</b> red legs = P, blue legs = S</div>
          </div>
        </div>
        </div>

        <script>
        const par = currentScript.previousElementSibling;
        // WideCell resizes `par`'s wrapper to its full (up to max_width) width via its own
        // ResizeObserver, which fires asynchronously after mount -- reading par.clientWidth
        // synchronously on script insertion races it and reliably gets the pre-resize (narrow,
        // default-cell-width) value instead, permanently undersizing both canvases. main() (the
        // sizing-sensitive setup) is deliberately started a few frames late, below, to let that
        // resize land first. The 'fermat-push' listener can't wait on that too, though -- the
        // FermatPush cell dispatches its first push as soon as the page loads, and a CustomEvent
        // fired before a listener exists is simply lost -- so it's attached here, immediately,
        // and buffers one push if main() (which replaces this with the real handler) hasn't run
        // yet.
        par.addEventListener('fermat-push', (ev)=>{
          if(par._fwApplyPush){ par._fwApplyPush(ev.detail); } else { par._fwPendingPush = ev.detail; }
        });
        function main(){
        const R_EARTH = $(R_EARTH);
        const DISC_DEPTHS = $(disc_depths);
        const DISC_LABELS = $(disc_labels);
        const BUMP_DEPTHS = $(bump_depths);
        const BUMP_LABELS = $(bump_labels);
        const PREM_DEFAULT = {
          discVpOuter: $([d.vp_outer for d in prem_named_discs]),
          discVpInner: $([d.vp_inner for d in prem_named_discs]),
          discVsOuter: $([d.vs_outer for d in prem_named_discs]),
          discVsInner: $([d.vs_inner for d in prem_named_discs]),
        };

        const SRC_DEPTH = 0.0; // fixed -- the source no longer moves; dragging from it aims instead
        let state = {
          discVpOuter: $(w.discVpOuter), discVpInner: $(w.discVpInner),
          discVsOuter: $(w.discVsOuter), discVsInner: $(w.discVsInner),
          bumpVp: $(w.bumpVp), bumpVs: $(w.bumpVs),
          wavetype0: "$(w.wavetype0)", convertOnExit: $(w.convertOnExit ? "true" : "false"),
          nRays: $(w.nRays), aimAngle: $(w.aimAngle), fanWidth: $(w.fanWidth),
          maxBounces: $(w.maxBounces),
        };

        function publish(){
          par.value = { ...state };
          par.dispatchEvent(new CustomEvent('input'));
        }

        const DPR = Math.min(window.devicePixelRatio || 1, 2);
        function hidpi(cv, cx, w, h){
          cv.width = Math.round(w*DPR); cv.height = Math.round(h*DPR);
          cv.style.width = w+'px'; cv.style.height = h+'px';
          cx.setTransform(DPR,0,0,DPR,0,0);
        }

        // Size the two square panels to actually fill the widget's width (up to WideCell's
        // max_width) instead of a fixed pixel size, the same availW pattern the flat notebook's
        // own widget uses -- computed once at render time, not with a live ResizeObserver.
        const availW = Math.min(window.innerWidth*0.88, par.clientWidth || window.innerWidth*0.88, 1500);
        // Two panels + the 16px gap + each panel's own 2px border must fit with room to spare --
        // flex-wrap is off (see .fw-workspace) so the panels never stack, only ever shrink.
        const SEC = Math.max(300, Math.floor((availW - 16 - 4) / 2));
        const mCanvas = par.querySelector('#fw-model'), mCtx = mCanvas.getContext('2d');
        const tCanvas = par.querySelector('#fw-tt'), tCtx = tCanvas.getContext('2d');
        const cbCanvas = par.querySelector('#fw-colorbar'), cbCtx = cbCanvas.getContext('2d');
        const CB_H = 34;
        hidpi(mCanvas, mCtx, SEC, SEC); hidpi(tCanvas, tCtx, SEC, SEC); hidpi(cbCanvas, cbCtx, SEC, CB_H);
        // Offscreen copy of the heatmap at the model canvas's own device-pixel size: putImageData
        // ignores the canvas transform (see computeHeatmap below), so it can't be zoomed/panned
        // directly. Rasterizing it here once and blitting it with drawImage on the visible canvas
        // instead lets the same zoom/pan transform that moves the rings and rays move the heatmap
        // too, for free.
        const heatCanvas = document.createElement('canvas');
        heatCanvas.width = mCanvas.width; heatCanvas.height = mCanvas.height;
        const heatCtx = heatCanvas.getContext('2d');

        // ---- Panels: filled in by the 'fermat-push' event from Julia ----
        const REARTH2 = R_EARTH;
        function polarXY(thetaDeg, canvasR, cx, cy){
          const th = thetaDeg*Math.PI/180;
          return [cx + canvasR*Math.sin(th), cy - canvasR*Math.cos(th)];
        }
        function toXY(x, z, canvasScale, cx, cy){
          return [cx + x*canvasScale, cy - z*canvasScale];
        }

        function drawStarMarker(ctx, cx, cy, r, fill, stroke){
          ctx.beginPath();
          for(let i=0;i<10;i++){
            const ang = -Math.PI/2 + i*Math.PI/5;
            const rr = i%2===0 ? r : r*0.42;
            const x = cx+rr*Math.cos(ang), y = cy+rr*Math.sin(ang);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.closePath();
          ctx.fillStyle = fill; ctx.fill();
          ctx.strokeStyle = stroke; ctx.lineWidth = 1.4; ctx.stroke();
        }

        const MCX=SEC/2, MCY=SEC/2, MR=SEC*0.46, MSCALE=MR/REARTH2;

        // ---- Model-panel zoom/pan (mirrors the interstation-pair notebook's ray-geometry widget) ----
        let zoom = 1, panX = 0, panY = 0;
        function constrainPan(){
          const limit = Math.max(0, (zoom-1)*SEC*0.42);
          panX = Math.max(-limit, Math.min(limit, panX));
          panY = Math.max(-limit, Math.min(limit, panY));
        }
        // Screen-space (mx,my) -> the fixed model-space coordinates everything below is drawn in,
        // so hit-testing and dragging stay accurate at any zoom/pan.
        function viewPoint(mx, my){
          return [MCX + (mx-MCX-panX)/zoom, MCY + (my-MCY-panY)/zoom];
        }
        // Clears the whole canvas in its own unscaled frame (a scaled clearRect would only clear
        // part of the canvas), then pushes the zoom/pan transform for `fn`'s draw calls, which are
        // otherwise unchanged -- they still draw in fixed model-space coordinates like MCX/MCY.
        function withZoom(fn){
          mCtx.clearRect(0,0,SEC,SEC);
          mCtx.save();
          mCtx.translate(MCX+panX, MCY+panY);
          mCtx.scale(zoom, zoom);
          mCtx.translate(-MCX, -MCY);
          fn();
          mCtx.restore();
        }
        function updateCursor(){
          if(!dragging) mCanvas.style.cursor = zoom > 1 ? 'grab' : 'default';
        }
        const zoomLevelEl = par.querySelector('#fw-zoomlevel');
        function setZoom(nextZoom){
          zoom = Math.max(0.6, Math.min(3, nextZoom));
          constrainPan();
          zoomLevelEl.textContent = Math.round(zoom*100) + '%';
          updateCursor();
          drawOverlay();
        }
        function resetView(){ panX = 0; panY = 0; setZoom(1); }

        // Velocity-model heatmap (grayscale: faster = brighter) with the ray fan overlaid directly
        // on top in red (P) / blue (S), colors that never occur in a grayscale ramp, so a ray
        // reads clearly against whatever structure it's crossing. The heatmap pixels are
        // cached (`hasHeatmap`) so dragging a ring can cheaply redraw the vector overlay (rings,
        // rays, star, drag marker) on every mousemove without re-running the per-pixel loop below,
        // which only real Julia recomputes (on release) actually need to redo.
        let hasHeatmap = false;

        // putImageData writes raw device pixels and ignores any canvas transform (unlike every
        // other draw call here, which goes through hidpi()'s setTransform plus the zoom/pan
        // transform from withZoom) -- so it's written to the offscreen `heatCanvas` at that
        // canvas's own device-pixel size, then blitted onto the visible canvas with drawImage
        // (which DOES respect the current transform) so it zooms/pans along with everything else.
        let lastVmin = 0, lastVmax = 1;
        // Gray level for a velocity `v` given the current field's [vmin,vmax] -- shared by the
        // heatmap and its colorbar so the bar is always an exact legend for what's on screen,
        // not a fixed scale that drifts out of sync with whichever field/edit is showing.
        function velToGray(v, vmin, vmax){
          const s = Math.max(0, Math.min(1, (v-vmin)/(vmax-vmin+1e-9)));
          return Math.round(25+205*s);
        }

        function computeHeatmap(lutR, lutV){
          const vmin = Math.min(...lutV.filter(v=>!isNaN(v))), vmax = Math.max(...lutV.filter(v=>!isNaN(v)));
          lastVmin = vmin; lastVmax = vmax;
          const W = mCanvas.width, H = mCanvas.height;
          const cx = W/2, cy = H/2, R = Math.min(W,H)*0.46;
          const img = heatCtx.createImageData(W,H);
          for(let py=0; py<H; py++) for(let px=0; px<W; px++){
            const x = (px-cx)/R*REARTH2, z = (cy-py)/R*REARTH2;
            const r = Math.hypot(x,z);
            const idx4 = (py*W+px)*4;
            if(r>REARTH2){ img.data[idx4+3]=0; continue; }
            const frac = r/REARTH2*(lutR.length-1);
            const i0 = Math.max(0, Math.min(lutR.length-2, Math.floor(frac)));
            const t = frac-i0;
            const v = lutV[i0]*(1-t) + lutV[i0+1]*t;
            const gray = velToGray(v, vmin, vmax); // grayscale: faster = brighter
            img.data[idx4] = gray; img.data[idx4+1] = gray; img.data[idx4+2] = gray; img.data[idx4+3] = 255;
          }
          heatCtx.putImageData(img, 0, 0);
        }

        // The source is fixed at the surface, straight up from center -- a constant screen point.
        const SRC_SX = MCX, SRC_SY = MCY - MR;

        function drawBase(){
          if(!hasHeatmap) return;
          mCtx.drawImage(heatCanvas, 0, 0, heatCanvas.width, heatCanvas.height, 0, 0, SEC, SEC);
          mCtx.strokeStyle = '#374151'; mCtx.beginPath(); mCtx.arc(MCX,MCY,MR,0,2*Math.PI); mCtx.stroke();
          for(const d of DISC_DEPTHS){
            const rr = (REARTH2-d)/REARTH2*MR;
            mCtx.strokeStyle = 'rgba(255,255,255,0.35)'; mCtx.setLineDash([2,3]);
            mCtx.beginPath(); mCtx.arc(MCX,MCY,rr,0,2*Math.PI); mCtx.stroke(); mCtx.setLineDash([]);
          }
          for(const d of BUMP_DEPTHS){
            const rr = (REARTH2-d)/REARTH2*MR;
            mCtx.strokeStyle = 'rgba(250,204,21,0.35)'; mCtx.setLineDash([1,4]);
            mCtx.beginPath(); mCtx.arc(MCX,MCY,rr,0,2*Math.PI); mCtx.stroke(); mCtx.setLineDash([]);
          }
        }

        function drawStar(){ drawStarMarker(mCtx, SRC_SX, SRC_SY, 6, '#facc15', '#4b5563'); }

        // Splits the NaN-separated flat arrays pushed from Julia back into one array per ray, so
        // the reveal animation can track each ray's own progress independently.
        function splitRays(flatX, flatZ, flatType){
          const rays = [];
          let cx_=[], cz_=[], ct_=[];
          for(let i=0;i<flatX.length;i++){
            if(isNaN(flatX[i])){ if(cx_.length>1) rays.push({x:cx_,z:cz_,t:ct_}); cx_=[]; cz_=[]; ct_=[]; continue; }
            cx_.push(flatX[i]); cz_.push(flatZ[i]); ct_.push(flatType[i]);
          }
          if(cx_.length>1) rays.push({x:cx_,z:cz_,t:ct_});
          return rays;
        }

        // Draws one ray's trail up to `frac` of its length (1 = the whole path), with a bright
        // glowing head at the current tip while frac<1 -- the "shooting" look during the reveal
        // animation. The path itself is exactly what trace_ray computed; this only decides how
        // much of it is visible yet, which is display timing, not physics.
        function drawRayPath(ray, frac){
          const n = ray.x.length;
          if(n<2) return;
          const upto = Math.max(1, Math.min(n, Math.ceil(n*frac)));
          mCtx.globalAlpha = 0.85; mCtx.lineWidth = 1.2;
          for(let i=1;i<upto;i++){
            const legColor = ray.t[i]>0.5 ? '#3b82f6' : '#ef4444'; // S = blue, P = red
            mCtx.strokeStyle = legColor;
            mCtx.beginPath();
            mCtx.moveTo(...toXY(ray.x[i-1], ray.z[i-1], MSCALE, MCX, MCY));
            mCtx.lineTo(...toXY(ray.x[i], ray.z[i], MSCALE, MCX, MCY));
            mCtx.stroke();
          }
          mCtx.globalAlpha = 1;
          if(frac<1 && upto<n){
            const [hx,hy] = toXY(ray.x[upto-1], ray.z[upto-1], MSCALE, MCX, MCY);
            const headColor = ray.t[upto-1]>0.5 ? '#3b82f6' : '#ef4444';
            mCtx.save();
            mCtx.shadowColor = headColor; mCtx.shadowBlur = 10;
            mCtx.beginPath(); mCtx.arc(hx,hy,2.6,0,2*Math.PI);
            mCtx.fillStyle = '#ffffff'; mCtx.fill();
            mCtx.restore();
          }
        }

        let lastRays = [], lastLabel = '';
        let animId = null;

        // Instant, fully-revealed redraw -- used while dragging a disc/bump ring (the ray fan
        // itself hasn't changed yet, only the drag marker/caption have) and as the animation's
        // settled end state.
        function drawOverlay(){
          withZoom(()=>{
            drawBase();
            for(const ray of lastRays) drawRayPath(ray, 1.0);
            drawStar();
            if(dragMarker){
              mCtx.beginPath(); mCtx.arc(dragMarker.x, dragMarker.y, 5, 0, 2*Math.PI);
              mCtx.fillStyle = '#ffffff'; mCtx.fill();
              mCtx.strokeStyle = '#000000'; mCtx.lineWidth = 1; mCtx.stroke();
            }
          });
          par.querySelector('#fw-model-caption').textContent = dragCaption || lastLabel;
        }

        // On every real recompute, "shoot" the new fan: each ray's head races along its own
        // already-computed path, leaving a trail, over a fixed duration -- animation timing only,
        // no physics recomputed per frame.
        function animateRayReveal(){
          if(animId) cancelAnimationFrame(animId);
          const t0 = performance.now(), DUR = 800;
          function step(now){
            const p = Math.min(1, (now-t0)/DUR);
            withZoom(()=>{
              drawBase();
              for(const ray of lastRays) drawRayPath(ray, p);
              drawStar();
            });
            par.querySelector('#fw-model-caption').textContent = lastLabel;
            if(p<1){ animId = requestAnimationFrame(step); }
            else { animId = null; drawOverlay(); }
          }
          animId = requestAnimationFrame(step);
        }

        // A legend for the grayscale heatmap: a horizontal gray ramp with numeric ticks spanning
        // the CURRENT field's own [vmin,vmax] (recomputed on every push, since editing the
        // structure changes that range) and a Vp/Vs label matching the P/S toggle.
        function drawColorbar(fieldName){
          cbCtx.clearRect(0,0,SEC,CB_H);
          const padL=46, padR=14, barY=4, barH=12;
          const barW = SEC-padL-padR;
          for(let x=0; x<barW; x++){
            const v = lastVmin + (x/barW)*(lastVmax-lastVmin);
            const gray = velToGray(v, lastVmin, lastVmax);
            cbCtx.fillStyle = 'rgb('+gray+','+gray+','+gray+')';
            cbCtx.fillRect(padL+x, barY, 1, barH);
          }
          cbCtx.strokeStyle = '#374151'; cbCtx.strokeRect(padL, barY, barW, barH);
          cbCtx.font = '10px sans-serif'; cbCtx.fillStyle = '#9ca3af';
          cbCtx.textAlign = 'left'; cbCtx.textBaseline = 'left';
          cbCtx.fillText(fieldName + ' (km/s)', 0, barY+barH-2);
          const step = niceStep(lastVmax-lastVmin, 4);
          cbCtx.textAlign = 'center'; cbCtx.textBaseline = 'top';
          const start = Math.ceil(lastVmin/step)*step;
          for(let v=start; v<=lastVmax+1e-9; v+=step){
            const x = padL + (v-lastVmin)/(lastVmax-lastVmin)*barW;
            cbCtx.beginPath(); cbCtx.moveTo(x, barY+barH); cbCtx.lineTo(x, barY+barH+3); cbCtx.stroke();
            cbCtx.fillText(v.toFixed(1), x, barY+barH+4);
          }
        }

        function drawModel(lutR, lutV, label, flatX, flatZ, flatType, fieldName){
          computeHeatmap(lutR, lutV);
          hasHeatmap = true;
          lastRays = splitRays(flatX, flatZ, flatType);
          lastLabel = label;
          drawColorbar(fieldName);
          animateRayReveal();
        }

        // The angles this fan will actually shoot at, evenly spaced -- used for both the ghost
        // preview (below) and could be cross-checked against ray_fan's own N (kept independent on
        // purpose: this is pure preview geometry, not a physics computation).
        function ghostFanAngles(){
          const n = Math.max(2, Math.min(state.nRays, 21));
          const lo = state.aimAngle - state.fanWidth/2, hi = state.aimAngle + state.fanWidth/2;
          const arr = [];
          for(let i=0;i<n;i++) arr.push(lo + (hi-lo)*i/(n-1));
          return arr;
        }

        // Faint straight-line fan around the current aim direction, drawn only while dragging the
        // source star -- an aiming guide, not a physics preview (real rays curve; these don't).
        function drawGhostFan(){
          const len = MR*1.05;
          mCtx.strokeStyle = 'rgba(255,255,255,0.22)'; mCtx.lineWidth = 1;
          for(const a of ghostFanAngles()){
            const rad = a*Math.PI/180;
            mCtx.beginPath(); mCtx.moveTo(SRC_SX,SRC_SY);
            mCtx.lineTo(SRC_SX+Math.sin(rad)*len, SRC_SY+Math.cos(rad)*len);
            mCtx.stroke();
          }
          const rad0 = state.aimAngle*Math.PI/180;
          const ax = SRC_SX+Math.sin(rad0)*len*0.72, ay = SRC_SY+Math.cos(rad0)*len*0.72;
          mCtx.strokeStyle = '#facc15'; mCtx.lineWidth = 2;
          mCtx.beginPath(); mCtx.moveTo(SRC_SX,SRC_SY); mCtx.lineTo(ax,ay); mCtx.stroke();
          const backAng = Math.atan2(ay-SRC_SY, ax-SRC_SX), headLen=8, headAng=0.4;
          mCtx.beginPath();
          mCtx.moveTo(ax,ay); mCtx.lineTo(ax-headLen*Math.cos(backAng-headAng), ay-headLen*Math.sin(backAng-headAng));
          mCtx.moveTo(ax,ay); mCtx.lineTo(ax-headLen*Math.cos(backAng+headAng), ay-headLen*Math.sin(backAng+headAng));
          mCtx.stroke();
        }

        function drawAimPreview(){
          if(animId){ cancelAnimationFrame(animId); animId=null; }
          withZoom(()=>{
            drawBase();
            drawGhostFan();
            drawStar();
          });
          par.querySelector('#fw-model-caption').textContent =
            'aim ' + state.aimAngle.toFixed(1) + '° · width ' + state.fanWidth.toFixed(1) + '° (release to shoot)';
        }

        // ---- Direct-manipulation editing, all on the model panel itself ----
        // Four kinds of drag target, checked in this priority order:
        //  1. the source star -- fixed in place; dragging from it aims the fan (shows the ghost
        //     preview above) instead of moving anything;
        //  2. a dashed discontinuity ring -- outside it (larger canvas radius, shallower depth)
        //     is that jump's "outer" side, inside is "inner", vertical drag distance is a
        //     velocity delta;
        //  3. a dotted bump ring (a smooth interior layer, see bump_specs in the Appendix) --
        //     vertical drag distance is an offset added on top of the exact PREM value, vanishing
        //     at that layer's two real discontinuities by construction (tent), so those stay put.
        // Which field -- Vp or Vs -- a ring drag edits is exactly whichever the P/S toggle has on
        // screen, so what's plotted is what's draggable. Nothing round-trips through Julia until
        // release, matching this repo's commit-on-release widget convention.
        let dragging = null, dragMarker = null, dragCaption = null, panStart = null;
        const DRAG_SENSITIVITY = 0.04; // km/s per pixel of vertical drag

        mCanvas.addEventListener('mousedown', (e)=>{
          const rect = mCanvas.getBoundingClientRect();
          const rawMx = e.clientX-rect.left, rawMy = e.clientY-rect.top;
          const [mx, my] = viewPoint(rawMx, rawMy);

          if(Math.hypot(mx-SRC_SX, my-SRC_SY) < 14){
            dragging = {type:'aim'};
            drawAimPreview();
            return;
          }

          const rho = Math.hypot(mx-MCX, my-MCY);
          if(rho <= MR+20){
            let bestK=-1, bestD=18;
            for(let k=0;k<DISC_DEPTHS.length;k++){
              const d = Math.abs(rho-(REARTH2-DISC_DEPTHS[k])/REARTH2*MR);
              if(d<bestD){ bestD=d; bestK=k; }
            }
            if(bestK>=0){
              const rr = (REARTH2-DISC_DEPTHS[bestK])/REARTH2*MR;
              const side = rho > rr ? 'Outer' : 'Inner';
              const field = state.wavetype0==='S' ? 'discVs' : 'discVp';
              const key = field+side;
              dragging = {type:'disc', k:bestK, key, startY:my, startVal: state[key][bestK]};
              dragMarker = {x:mx, y:my};
              drawOverlay();
              return;
            }

            let bestJ=-1, bestD2=18;
            for(let j=0;j<BUMP_DEPTHS.length;j++){
              const d = Math.abs(rho-(REARTH2-BUMP_DEPTHS[j])/REARTH2*MR);
              if(d<bestD2){ bestD2=d; bestJ=j; }
            }
            if(bestJ>=0){
              const field = state.wavetype0==='S' ? 'bumpVs' : 'bumpVp';
              dragging = {type:'bump', j:bestJ, field, startY:my, startVal: state[field][bestJ]};
              dragMarker = {x:mx, y:my};
              drawOverlay();
              return;
            }
          }

          // Nothing to edit under the cursor -- if zoomed in, empty space pans the view instead.
          // Pan tracks raw (unconverted) screen coordinates, unlike every drag above, since it's
          // the zoom/pan transform itself that's being adjusted, not something drawn inside it.
          if(zoom > 1){
            dragging = {type:'pan'};
            panStart = {x: rawMx, y: rawMy, panX, panY};
            mCanvas.style.cursor = 'grabbing';
          }
        });
        window.addEventListener('mousemove', (e)=>{
          if(!dragging) return;
          const rect = mCanvas.getBoundingClientRect();
          const rawMx = e.clientX-rect.left, rawMy = e.clientY-rect.top;

          if(dragging.type==='pan'){
            panX = panStart.panX + rawMx - panStart.x;
            panY = panStart.panY + rawMy - panStart.y;
            constrainPan();
            drawOverlay();
            return;
          }

          const [mx, my] = viewPoint(rawMx, rawMy);

          if(dragging.type==='aim'){
            const dx = mx-SRC_SX, dy = my-SRC_SY;
            if(Math.hypot(dx,dy) < 4) return; // avoid atan2 jitter right at the source
            const angle = Math.atan2(dx,dy)*180/Math.PI;
            state.aimAngle = Math.max(0.5, Math.min(89.5, Math.abs(angle)));
            drawAimPreview();
            return;
          }

          const dv = (dragging.startY-my)*DRAG_SENSITIVITY;
          if(dragging.type==='disc'){
            const v = Math.max(0.2, Math.min(14.5, dragging.startVal+dv));
            state[dragging.key][dragging.k] = v;
            dragCaption = DISC_LABELS[dragging.k] + ' ' + dragging.key + ': ' + v.toFixed(2) + ' km/s (release to recompute)';
          } else if(dragging.type==='bump'){
            const v = Math.max(-3, Math.min(3, dragging.startVal+dv));
            state[dragging.field][dragging.j] = v;
            dragCaption = BUMP_LABELS[dragging.j] + ' Δ' + (dragging.field==='bumpVs'?'Vs':'Vp') + ': ' +
              (v>=0?'+':'') + v.toFixed(2) + ' km/s (release to recompute)';
          }
          dragMarker = {x:dragMarker.x, y:my};
          drawOverlay();
        });
        window.addEventListener('mouseup', ()=>{
          if(dragging){
            const wasAim = dragging.type==='aim', wasPan = dragging.type==='pan';
            dragging=null; dragMarker=null; dragCaption=null; panStart=null;
            updateCursor();
            if(wasAim) syncControls();
            if(!wasPan) publish(); // panning is a view-only change, nothing for Julia to recompute
          }
        });

        // "Nice" tick step (1/2/5 x a power of ten) for a given axis range and target tick count --
        // the standard trick behind almost every plotting library's default axis ticks.
        function niceStep(range, targetCount){
          if(!(range>0)) return 1;
          const rough = range/targetCount;
          const mag = Math.pow(10, Math.floor(Math.log10(rough)));
          const norm = rough/mag;
          const step = norm<1.5 ? 1 : norm<3 ? 2 : norm<7 ? 5 : 10;
          return step*mag;
        }

        function drawTT(deltas, times){
          tCtx.clearRect(0,0,SEC,SEC);
          const padL=46, padR=14, padT=10, padB=32;
          const plotW = SEC-padL-padR, plotH = SEC-padT-padB;
          const valid = [];
          for(let i=0;i<deltas.length;i++) if(!isNaN(deltas[i]) && !isNaN(times[i])) valid.push([deltas[i], times[i]]);
          tCtx.strokeStyle = '#374151'; tCtx.strokeRect(padL, padT, plotW, plotH);
          if(valid.length<2){ par.querySelector('#fw-tt-caption').textContent = 'no rays returned to the surface (shadow zone)'; return; }
          const dmax = Math.max(...valid.map(p=>p[0])), tmax = Math.max(...valid.map(p=>p[1]));
          function toPx(d,t){ return [padL + d/dmax*plotW, padT+plotH - t/tmax*plotH]; }

          tCtx.font = '10px sans-serif'; tCtx.fillStyle = '#9ca3af'; tCtx.strokeStyle = '#374151';
          const dStep = niceStep(dmax, 6);
          tCtx.textAlign = 'center'; tCtx.textBaseline = 'top';
          for(let d=0; d<=dmax+1e-9; d+=dStep){
            const [x] = toPx(d,0);
            tCtx.beginPath(); tCtx.moveTo(x, padT+plotH); tCtx.lineTo(x, padT+plotH+4); tCtx.stroke();
            tCtx.fillText(d.toFixed(0), x, padT+plotH+6);
          }
          const tStep = niceStep(tmax, 6);
          tCtx.textAlign = 'right'; tCtx.textBaseline = 'middle';
          for(let t=0; t<=tmax+1e-9; t+=tStep){
            const [, y] = toPx(0,t);
            tCtx.beginPath(); tCtx.moveTo(padL-4, y); tCtx.lineTo(padL, y); tCtx.stroke();
            tCtx.fillText(t.toFixed(0), padL-6, y);
          }

          // Connect consecutive LANDED rays in takeoff-angle order (the order Julia shot them
          // in, not sorted by Δ) -- a triplication fold or a shadow-zone gap is a specific,
          // load-bearing shape in *that* order, and is easy to miss in a bare scatter of dots.
          // Break the line at an unlanded ray or an abnormally large Δ jump between neighbors
          // (a real branch discontinuity, e.g. crossing the shadow zone) rather than drawing a
          // misleading diagonal across it.
          const JUMP_BREAK_DEG = 15;
          tCtx.strokeStyle = '#38bdf8'; tCtx.lineWidth = 1; tCtx.globalAlpha = 0.6;
          let prev = null;
          for(let i=0;i<deltas.length;i++){
            const d = deltas[i], t = times[i];
            if(isNaN(d) || isNaN(t)){ prev = null; continue; }
            if(prev && Math.abs(d-prev[0]) <= JUMP_BREAK_DEG){
              const [x0,y0] = toPx(prev[0],prev[1]), [x1,y1] = toPx(d,t);
              tCtx.beginPath(); tCtx.moveTo(x0,y0); tCtx.lineTo(x1,y1); tCtx.stroke();
            }
            prev = [d,t];
          }
          tCtx.globalAlpha = 1;

          tCtx.fillStyle = '#38bdf8';
          for(const [d,t] of valid){
            const [x,y] = toPx(d,t);
            tCtx.beginPath(); tCtx.arc(x,y,2,0,2*Math.PI); tCtx.fill();
          }

          tCtx.fillStyle = '#9ca3af'; tCtx.textAlign = 'center'; tCtx.textBaseline = 'alphabetic';
          tCtx.fillText('epicentral distance Δ (deg)', padL+plotW/2, SEC-4);
          tCtx.save(); tCtx.translate(12, padT+plotH/2); tCtx.rotate(-Math.PI/2);
          tCtx.fillText('travel time T (s)', 0, 0); tCtx.restore();
          tCtx.textAlign = 'left';
          par.querySelector('#fw-tt-caption').textContent = valid.length + ' of ' + deltas.length + ' rays reached the surface';
        }

        // Replaces the buffering stub registered immediately at script start (see the top of
        // this script for why): from here on, a push is applied directly, and any push that
        // arrived before main() got this far is replayed once, right below.
        par._fwApplyPush = (d) => {
          drawModel(d.lutR, d.wavetype0==='S' ? d.lutVs : d.lutVp, 'field: ' + d.wavetype0,
            d.flatX, d.flatZ, d.flatType, d.wavetype0==='S' ? 'Vs' : 'Vp');
          drawTT(d.deltas, d.times);
        };
        if(par._fwPendingPush){ par._fwApplyPush(par._fwPendingPush); par._fwPendingPush = null; }

        // ---- Controls wiring ----
        function syncControls(){
          par.querySelector('#fw-nrays').value = state.nRays;
          par.querySelector('#fw-nrays-v').textContent = state.nRays;
          par.querySelector('#fw-fanwidth').value = state.fanWidth;
          par.querySelector('#fw-fanwidth-v').textContent = state.fanWidth.toFixed(0)+'°';
          par.querySelector('#fw-aim-readout').textContent =
            'aim angle: ' + state.aimAngle.toFixed(1) + '° from straight down';
          par.querySelector('#fw-convert').checked = state.convertOnExit;
          par.querySelector('#fw-wtp').classList.toggle('active', state.wavetype0==='P');
          par.querySelector('#fw-wts').classList.toggle('active', state.wavetype0==='S');
          par.querySelector('#fw-bounces').value = String(state.maxBounces);
        }
        syncControls();

        par.querySelector('#fw-nrays').addEventListener('change', (e)=>{ state.nRays=parseInt(e.target.value); syncControls(); publish(); });
        par.querySelector('#fw-fanwidth').addEventListener('change', (e)=>{ state.fanWidth=parseFloat(e.target.value); syncControls(); publish(); });
        par.querySelector('#fw-bounces').addEventListener('change', (e)=>{ state.maxBounces=parseInt(e.target.value); publish(); });
        par.querySelector('#fw-convert').addEventListener('change', (e)=>{ state.convertOnExit=e.target.checked; publish(); });
        par.querySelector('#fw-wtp').addEventListener('click', ()=>{ state.wavetype0='P'; syncControls(); publish(); });
        par.querySelector('#fw-wts').addEventListener('click', ()=>{ state.wavetype0='S'; syncControls(); publish(); });

        par.querySelector('#fw-zoomin').addEventListener('click', ()=> setZoom(zoom*1.25));
        par.querySelector('#fw-zoomout').addEventListener('click', ()=> setZoom(zoom/1.25));
        par.querySelector('#fw-zoomreset').addEventListener('click', resetView);

        // Each preset builds a COMPLETE, explicit state from scratch -- wave type, conversion
        // flag, ray-fan density/aim/width, everything -- not a partial merge on top of whatever
        // was left over from the previously selected preset. Cycling through presets in any order
        // always lands in exactly the state described here, never a hybrid of two. aimAngle/
        // fanWidth values are the midpoint/span of the angleMin/angleMax ranges validated earlier
        // (e.g. sks: 2-12° -> aim 7°, width 10°) -- same shot, reparameterized around the source.
        function freshDefaults(){
          return { discVpOuter: [...PREM_DEFAULT.discVpOuter], discVpInner: [...PREM_DEFAULT.discVpInner],
                   discVsOuter: [...PREM_DEFAULT.discVsOuter], discVsInner: [...PREM_DEFAULT.discVsInner],
                   bumpVp: BUMP_DEPTHS.map(()=>0.0), bumpVs: BUMP_DEPTHS.map(()=>0.0),
                   wavetype0: 'P', convertOnExit: false, nRays: 41, aimAngle: 43.0, fanWidth: 84.0,
                   maxBounces: 0 };
        }
        function buildPresetState(name){
          const s = freshDefaults();
          if(name==='uniform'){
            s.discVpOuter = DISC_DEPTHS.map(()=>8.0); s.discVpInner = DISC_DEPTHS.map(()=>8.0);
            s.discVsOuter = DISC_DEPTHS.map(()=>4.6); s.discVsInner = DISC_DEPTHS.map(()=>4.6);
          } else if(name==='nocore'){
            // The shadow zone turns out to be controlled by the *bulk* slowness of the outer
            // core, not just the jump at the CMB itself -- flattening only the boundary value
            // barely moved it (a ray that dips even slightly past the very edge of the core
            // still crosses the real, unedited, genuinely slow bulk of the layer). Raising the
            // whole layer via its bump -- which tapers to 0 exactly at the CMB/ICB, so those
            // stay real PREM -- is what actually closes it.
            const k = DISC_LABELS.indexOf('CMB');
            const mid = (s.discVpOuter[k]+s.discVpInner[k])/2;
            s.discVpOuter[k] = mid; s.discVpInner[k] = mid;
            const j = BUMP_LABELS.indexOf('Outer Core');
            s.bumpVp[j] = 4.5;
            s.nRays = 60; s.aimAngle = 45.0; s.fanWidth = 88.0;
          } else if(name==='sharp660'){
            const k = DISC_LABELS.indexOf('670');
            s.discVpOuter[k] -= 0.4; s.discVpInner[k] += 0.4;
            s.nRays = 100; s.aimAngle = 21.5; s.fanWidth = 37.0;
          } else if(name==='sks'){
            s.wavetype0 = 'S'; s.convertOnExit = true;
            s.nRays = 25; s.aimAngle = 7.0; s.fanWidth = 10.0;
          }
          return s;
        }
        function applyPreset(name){
          state = buildPresetState(name);
          syncControls(); publish();
        }
        par.querySelector('#fw-preset').addEventListener('change', (e)=> applyPreset(e.target.value));
        par.querySelector('#fw-reset').addEventListener('click', ()=>{ par.querySelector('#fw-preset').value='prem'; applyPreset('prem'); });
        } // end main()
        let _fwBootFrames = 0;
        (function fwBoot(){
          _fwBootFrames++;
          if(_fwBootFrames < 5){ requestAnimationFrame(fwBoot); return; }
          main();
        })();
        </script>
        """)
    end

    const _fermat_ready = true
end

# ╔═╡ 6b807f1d-2b0c-423e-9074-3df06dfe5ae4
begin
    _fermat_ready
    WideCell(@bind fermat_bond FermatGlobalPhasesInput(); max_width=1500)
end

# ╔═╡ 31be046e-f288-4376-932b-50ce5817204c
begin
    _fb = fermat_bond isa AbstractDict ? fermat_bond : Base.get(FermatGlobalPhasesInput())
    _fb_discVpOuter = Float64.(_fb["discVpOuter"])
    _fb_discVpInner = Float64.(_fb["discVpInner"])
    _fb_discVsOuter = Float64.(_fb["discVsOuter"])
    _fb_discVsInner = Float64.(_fb["discVsInner"])
    _fb_bumpVp = Float64.(_fb["bumpVp"])
    _fb_bumpVs = Float64.(_fb["bumpVs"])
    _fb_wavetype0 = _fb["wavetype0"] == "S" ? :S : :P
    _fb_convertOnExit = _fb["convertOnExit"] isa Bool ? _fb["convertOnExit"] : _fb["convertOnExit"] == true
    _fb_nRays = clamp(round(Int, _fb["nRays"]), 3, 200)
    _fb_aimAngle = Float64(_fb["aimAngle"])
    _fb_fanWidth = Float64(_fb["fanWidth"])
    _fb_angleMin = clamp(_fb_aimAngle - _fb_fanWidth / 2, 0.5, 89.5)
    _fb_angleMax = clamp(_fb_aimAngle + _fb_fanWidth / 2, 0.5, 89.5)
    _fb_maxBounces = clamp(round(Int, get(_fb, "maxBounces", 0)), 0, 3)

    _fb_r = copy(prem_r)
    _fb_vp = copy(prem_vp)
    _fb_vs = copy(prem_vs)
    for (k, d) in enumerate(prem_named_discs)
        _fb_vp[d.i] = _fb_discVpOuter[k]
        _fb_vp[d.i+1] = _fb_discVpInner[k]
        _fb_vs[d.i] = _fb_discVsOuter[k]
        _fb_vs[d.i+1] = _fb_discVsInner[k]
        # Fade each edit into the real profile over a short depth range so a shrunk (or widened)
        # jump doesn't leave a hidden near-cliff of un-edited original data sitting right beside
        # it -- see blend_edited_boundary!'s docstring for why that alone would silently undo
        # the point of the edit.
        blend_edited_boundary!(_fb_vp, prem_vp, d.i, -1)
        blend_edited_boundary!(_fb_vp, prem_vp, d.i + 1, 1)
        blend_edited_boundary!(_fb_vs, prem_vs, d.i, -1)
        blend_edited_boundary!(_fb_vs, prem_vs, d.i + 1, 1)
    end
    # Smooth-interior bumps: added to every sample point in a layer via `tent`, which is exactly
    # 0 at that layer's two bounding discontinuities -- so this can never disturb the boundary
    # values just set above, only reshape the gradient strictly between them.
    for (j, spec) in enumerate(bump_specs)
        for i in eachindex(_fb_r)
            w = tent(R_EARTH - _fb_r[i], spec.lo, spec.center, spec.hi)
            w == 0.0 && continue
            _fb_vp[i] += _fb_bumpVp[j] * w
            _fb_vs[i] += _fb_bumpVs[j] * w
        end
    end
    _fb_profile = (r=_fb_r, vp=_fb_vp, vs=_fb_vs)
    _fb_discs = prem_discontinuities(_fb_r, _fb_vp, _fb_vs)
    _fb_src_r = R_EARTH  # source is fixed at the surface
    # Each bounce roughly re-traverses the mantle again, so the step budget needs to grow with
    # the number of bounces allowed or a multi-bounce ray can silently hit max_steps mid-flight.
    _fb_rays = ray_fan(_fb_profile, _fb_discs, _fb_nRays, _fb_wavetype0, _fb_angleMin, _fb_angleMax;
        src_r=_fb_src_r, convert_on_core_exit=_fb_convertOnExit, max_bounces=_fb_maxBounces,
        max_steps=8000 * (_fb_maxBounces + 1))
    "ray fan recomputed"
end

# ╔═╡ 15c5334a-f850-43b9-a25e-69b2eb6ff2b0
begin
    """
    	FermatPush(rays, profile, wavetype0)

    A push-only widget: its `Base.show` emits nothing but a `<script>` that finds the already-
    mounted `FermatGlobalPhasesInput` widget by id and dispatches a `fermat-push` `CustomEvent`
    carrying the freshly computed ray fan and velocity-vs-radius lookup tables. Mirrors the
    `RayPaintInput`/`CwPush` pattern used elsewhere in this repo -- the widget is never
    re-rendered, only its already-drawn panels are updated in place. The source is fixed at the
    surface, so unlike earlier versions of this push there is no source depth/radius to carry --
    the JS side already knows `R_EARTH`.
    """
    struct FermatPush
        rays::Vector{Any}
        lutR::Vector{Float64}
        lutVp::Vector{Float64}
        lutVs::Vector{Float64}
        wavetype0::String
    end

    function Base.show(io::IO, ::MIME"text/html", p::FermatPush)
        flat_x = Float64[]
        flat_z = Float64[]
        flat_t = Float64[]
        for ray in p.rays
            append!(flat_x, ray.x)
            push!(flat_x, NaN)
            append!(flat_z, ray.z)
            push!(flat_z, NaN)
            append!(flat_t, [wt == :P ? 0.0 : 1.0 for wt in ray.wavetype])
            push!(flat_t, NaN)
        end
        deltas = [r.delta_deg === nothing ? NaN : r.delta_deg for r in p.rays]
        times = [r.landed ? r.t[end] : NaN for r in p.rays]
        write(io, """
        <script>
        {
        const w = document.getElementById('fermatwidget');
        if(w){
          w.dispatchEvent(new CustomEvent('fermat-push', { detail: {
            flatX: $(flat_x), flatZ: $(flat_z), flatType: $(flat_t),
            deltas: $(deltas), times: $(times),
            lutR: $(p.lutR), lutVp: $(p.lutVp), lutVs: $(p.lutVs),
            wavetype0: "$(p.wavetype0)",
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ cfeabbec-d52e-44ef-8cd8-ca9429ab3a1c
let
    lut_r = collect(range(0.0, R_EARTH; length=200))
    lut_vp = [segment_velocity(_fb_r, _fb_vp, r)[1] for r in lut_r]
    lut_vs = [segment_velocity(_fb_r, _fb_vs, r)[1] for r in lut_r]
    FermatPush(_fb_rays, lut_r, lut_vp, lut_vs, string(_fb_wavetype0))
end

# ╔═╡ 2b54656d-3150-4b15-8f37-d6cfbd13b810
"""
	blend_edited_boundary!(v, orig_v, idx, direction; nblend=8)

After overwriting `v[idx]` with a new discontinuity value, the very next raw PREM sample points
in `direction` (`+1` deeper, `-1` shallower) still hold their *original* values -- for a jump this
notebook has deliberately shrunk (e.g. "no core contrast"), that leaves an unrealistic near-cliff
sitting immediately beside the edited boundary, which bends rays almost as much as the original
jump did and defeats the whole point of editing it. This linearly tapers the next `nblend` samples
from the new value back to their own original PREM value, so the edit fades in over a short but
genuine depth range instead of leaving that hidden cliff.
"""
function blend_edited_boundary!(v, orig_v, idx, direction; nblend=8)
    for step in 1:nblend
        j = idx + direction * step
        (j < 1 || j > length(v)) && break
        t = step / (nblend + 1)
        v[j] = (1 - t) * v[idx] + t * orig_v[j]
    end
end

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

julia_version = "1.13.1"
manifest_format = "2.1"
project_hash = "884d130e167b68a435fbfc6edb6d38588b98bbd0"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "daa72978cd7a624246e894a4f4f067706d4e17e2"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
registries = "General"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "60f11b38ebeabd984f5535838d91e197d97202f0"
registries = "General"
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
registries = "General"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
registries = "General"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
registries = "General"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
registries = "General"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
registries = "General"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
registries = "General"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
registries = "General"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
registries = "General"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChunkCodecCore]]
git-tree-sha1 = "1a3ad7e16a321667698a19e77362b35a1e94c544"
registries = "General"
uuid = "0b6fb165-00bc-4d37-ab8b-79f91016dbe1"
version = "1.0.1"

[[deps.ChunkCodecLibZlib]]
deps = ["ChunkCodecCore", "Zlib_jll"]
git-tree-sha1 = "cee8104904c53d39eb94fd06cbe60cb5acde7177"
registries = "General"
uuid = "4c0bbee4-addc-4d73-81a0-b6caacae83c8"
version = "1.0.0"

[[deps.ChunkCodecLibZstd]]
deps = ["ChunkCodecCore", "Zstd_jll"]
git-tree-sha1 = "34d9873079e4cb3d0c62926a225136824677073f"
registries = "General"
uuid = "55437552-ac27-4d47-9aa3-63184e8fd398"
version = "1.0.0"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
registries = "General"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
registries = "General"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
registries = "General"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
registries = "General"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
registries = "General"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ef2022bff55342a8c9846cdf218f62e475f0444d"
registries = "General"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.1.2"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
registries = "General"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
registries = "General"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
registries = "General"
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
registries = "General"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.4"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
registries = "General"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
registries = "General"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
registries = "General"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
registries = "General"
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
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.Eikonal]]
deps = ["DataStructures", "Images", "LinearAlgebra", "PrecompileTools", "Printf"]
git-tree-sha1 = "ac89a6cf8c89a741448deb8692aaacba745ecee0"
registries = "General"
uuid = "a6aab1ba-8f88-4217-b671-4d0788596809"
version = "0.1.1"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
registries = "General"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
registries = "General"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6866aec60ef98e3164cd8d6855225684207e9dff"
registries = "General"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.12+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "8e9c059d6857607253e837730dbf780b6b151acd"
registries = "General"
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
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
registries = "General"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
registries = "General"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
registries = "General"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7a98c6502f4632dbe9fb1973a4244eaa3324e84d"
registries = "General"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.1"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
registries = "General"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.HistogramThresholding]]
deps = ["ImageBase", "LinearAlgebra", "MappedArrays"]
git-tree-sha1 = "7194dfbb2f8d945abdaf68fa9480a965d6661e69"
registries = "General"
uuid = "2c695a8d-9458-5d45-9878-1b8a99cf7853"
version = "0.3.1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
registries = "General"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
registries = "General"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
registries = "General"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
registries = "General"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageBinarization]]
deps = ["HistogramThresholding", "ImageCore", "LinearAlgebra", "Polynomials", "Reexport", "Statistics"]
git-tree-sha1 = "33485b4e40d1df46c806498c73ea32dc17475c59"
registries = "General"
uuid = "cbc4b850-ae4b-5111-9e64-df94c024a13d"
version = "0.3.1"

[[deps.ImageContrastAdjustment]]
deps = ["ImageBase", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "eb3d4365a10e3f3ecb3b115e9d12db131d28a386"
registries = "General"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.12"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
registries = "General"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageCorners]]
deps = ["ImageCore", "ImageFiltering", "PrecompileTools", "StaticArrays", "StatsBase"]
git-tree-sha1 = "24c52de051293745a9bad7d73497708954562b79"
registries = "General"
uuid = "89d5987c-236e-4e32-acd0-25bd6bd87b70"
version = "0.1.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "08b0e6354b21ef5dd5e49026028e41831401aca8"
registries = "General"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.17"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "52116260a234af5f69969c5286e6a5f8dc3feab8"
registries = "General"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.12"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
registries = "General"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils"]
git-tree-sha1 = "8e64ab2f0da7b928c8ae889c514a52741debc1c2"
registries = "General"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.4.2"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Bzip2_jll", "FFTW_jll", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "Zstd_jll", "libpng_jll", "libwebp_jll", "libzip_jll"]
git-tree-sha1 = "d670e8e3adf0332f57054955422e85a4aec6d0b0"
registries = "General"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "7.1.2005+0"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
registries = "General"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.ImageMorphology]]
deps = ["DataStructures", "ImageCore", "LinearAlgebra", "LoopVectorization", "OffsetArrays", "Requires", "TiledIteration"]
git-tree-sha1 = "cffa21df12f00ca1a365eb8ed107614b40e8c6da"
registries = "General"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.4.6"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "PrecompileTools", "Statistics"]
git-tree-sha1 = "783b70725ed326340adf225be4889906c96b8fd1"
registries = "General"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.7"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "7196039573b6f312864547eb7a74360d6c0ab8e6"
registries = "General"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.9.0"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
registries = "General"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "dfde81fafbe5d6516fb864dc79362c5c6b973c82"
registries = "General"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.10.2"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageBinarization", "ImageContrastAdjustment", "ImageCore", "ImageCorners", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "a49b96fd4a8d1a9a718dfd9cde34c154fc84fcd5"
registries = "General"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.26.2"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
registries = "General"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
registries = "General"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
registries = "General"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "b842cbff3f44804a84fda409745cc8f04c029a20"
registries = "General"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.6"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
registries = "General"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
registries = "General"
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
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
registries = "General"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLD2]]
deps = ["ChunkCodecLibZlib", "ChunkCodecLibZstd", "FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues"]
git-tree-sha1 = "941f87a0ae1b14d1ac2fa57245425b23a9d7a516"
registries = "General"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.6.4"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
registries = "General"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1dae3057da6f2b9c857afef03177bbdc7c4afe92"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
registries = "General"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
registries = "General"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
registries = "General"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.18.0+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.1+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
registries = "General"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
registries = "General"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "38928f7999753af13d4e13966ae15958ff3a917a"
registries = "General"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.19.1+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "bba2d9aa057d8f126415de240573e86a8f39d2a1"
registries = "General"
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
registries = "General"
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
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
registries = "General"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
registries = "General"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
registries = "General"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "3a8f462a180a9d735e340f4e8d5f364d411da3a4"
registries = "General"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.8.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
registries = "General"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ca7e18198a166a1f3eb92a3650d53d94ed8ca8a1"
registries = "General"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.22"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
registries = "General"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
registries = "General"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
registries = "General"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
registries = "General"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "215a6666fee6d6b3a6e75f2cc22cb767e2dd393a"
registries = "General"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.5+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
registries = "General"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
registries = "General"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
registries = "General"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
registries = "General"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
registries = "General"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "972089912ba299fba87671b025cd0da74f5f54f7"
registries = "General"
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
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
registries = "General"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
registries = "General"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "994cc27cdacca10e68feb291673ec3a76aa2fae9"
registries = "General"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.6"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
registries = "General"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
registries = "General"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
registries = "General"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
registries = "General"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
registries = "General"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
registries = "General"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"
weakdeps = ["RecipesBase"]

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
registries = "General"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
registries = "General"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
registries = "General"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLPublic]]
git-tree-sha1 = "cf9aaf8b9ed5db993259ea8b24cf2b7ba9bd3b79"
registries = "General"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.2.4"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "67a144433c4ce877ee6d1ada69a124d6b1ecf7be"
registries = "General"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.6.2"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
registries = "General"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
registries = "General"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "3e5f165e58b18204aed03158664c4982d691f454"
registries = "General"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.5.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
registries = "General"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
registries = "General"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
registries = "General"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "1e44e7b1dbb5249876d84c32466f8988a6b41bbb"
registries = "General"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
registries = "General"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
registries = "General"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
registries = "General"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "e4d7a1a0edc20af42689ea6f4f3587a2175d50ee"
registries = "General"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.12"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

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
registries = "General"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
registries = "General"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "98b9352a24cb6a2066f9ababcc6802de9aed8ad8"
registries = "General"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.6"

[[deps.TiledIteration]]
deps = ["OffsetArrays", "StaticArrayInterface"]
git-tree-sha1 = "1176cc31e867217b06928e2f140c90bd1bc88283"
registries = "General"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.5.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "3b0738bd7c5645641845da25cbd99800b8718689"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
registries = "General"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "d1d9a935a26c475ebffd54e9c7ad11627c43ea85"
registries = "General"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.72"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
registries = "General"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
registries = "General"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
registries = "General"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
registries = "General"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
registries = "General"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
registries = "General"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
registries = "General"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
registries = "General"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
registries = "General"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
registries = "General"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
registries = "General"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
registries = "General"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.libzip_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "86addc139bca85fdf9e7741e10977c45785727b7"
registries = "General"
uuid = "337d8026-41b4-5cde-a456-74a10e5b31d1"
version = "1.11.3+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "da8c1f6eee04831f14edcfa5dae611d309807e57"
registries = "General"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.3.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.2+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╠═875d305b-6fa0-441f-b4ea-673e1053c669
# ╠═95c66aa1-b555-4f38-a3a6-c79746262c87
# ╟─e785f801-df99-4fb7-ad1e-34861320c9bf
# ╟─6b807f1d-2b0c-423e-9074-3df06dfe5ae4
# ╟─5be71c78-0c91-4d72-aa5d-da8763ea02c1
# ╠═fa8a8574-165b-42a4-a14f-7ad7fb4b32f1
# ╠═11042f31-3dbb-4a2a-89b4-598f6e72df47
# ╠═6c0c8399-320f-4381-a007-dfccc0f7eddd
# ╠═87825c53-c340-47a1-9544-7be10960fe8f
# ╟─738f4b56-13d1-459f-a8ac-4b22aab71c8a
# ╠═6e3900a8-b1f4-458e-b7cf-8940ab465279
# ╠═ccbed1cc-ce29-41b8-8bb1-d2f203735696
# ╠═47f4eb99-0684-49cf-89eb-7b63c9d4407a
# ╠═5ec8e794-8b2e-4355-820a-db771c8cdece
# ╠═c48a015d-66dd-407a-92c9-ddc2088612b3
# ╟─83953fdc-7750-4da6-97a8-00d08b3c7f7d
# ╠═3b85d627-6ba9-4bda-b571-489a774e1e81
# ╟─199e63bd-a95e-4fb9-bc26-1ca0c05dd12e
# ╠═86bf927a-b3c9-4edb-888d-8ca18afa299a
# ╠═35793516-2618-4552-8dee-cb23a59296de
# ╠═0dd771df-6970-4eec-a26a-975319f0d697
# ╠═4d775817-c78d-4e2a-8d5f-f951daba5ce3
# ╠═ec425e80-1077-4f29-ba4a-befc34a94ea4
# ╠═b97f0fa4-eb4d-425b-bdc4-b63a983ca8ed
# ╠═c2c98c49-6b83-4907-8e66-a769d35e137d
# ╟─05a9c088-4976-4034-89df-ddbcdfe46b00
# ╠═4a03c3a7-ec91-4f9a-9b2a-d8c34c7e2af5
# ╠═d46fabe3-57ff-4d98-8eb6-266da1a58c2f
# ╠═dd31d8ea-11e5-4f84-8286-e5ef65c5ca96
# ╠═a8fa560a-01ee-4d42-89e6-5feaab505c54
# ╠═31be046e-f288-4376-932b-50ce5817204c
# ╠═15c5334a-f850-43b9-a25e-69b2eb6ff2b0
# ╠═cfeabbec-d52e-44ef-8cd8-ca9429ab3a1c
# ╠═2b54656d-3150-4b15-8f37-d6cfbd13b810
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
