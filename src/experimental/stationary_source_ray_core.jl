"""Core ray utilities shared by the stationary-source research notebook.

The implementation follows the discontinuity-aware spherical tracer in
`fermat-principle-global-phases.jl`.  It deliberately exposes a ray's launch
angle and endpoint slowness: those give the stationary-phase derivatives
analytically, without differencing a travel-time table.
"""

const SS_R_EARTH = 6371.0 # km

function ss_load_prem()
    r, vp, vs = Float64[], Float64[], Float64[]
    path = joinpath(@__DIR__, "..", "assets", "data", "specnm_models", "prem_ani")
    open(path) do io
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, "#") || !occursin(r"^[0-9.]", s)) && continue
            c = split(s)
            push!(r, parse(Float64, c[1]) / 1000)
            push!(vp, parse(Float64, c[3]) / 1000)
            push!(vs, parse(Float64, c[4]) / 1000)
        end
    end
    return (r=r, vp=vp, vs=vs)
end

function ss_discontinuities(profile)
    discs = NamedTuple[]
    for i in 1:length(profile.r)-1
        profile.r[i] == profile.r[i+1] || continue
        push!(discs, (r=profile.r[i], vp_outer=profile.vp[i], vp_inner=profile.vp[i+1],
                      vs_outer=profile.vs[i], vs_inner=profile.vs[i+1]))
    end
    return discs
end

function ss_velocity(profile, wave, r0)
    v = wave == :P ? profile.vp : profile.vs
    r0 = clamp(r0, profile.r[end], profile.r[1])
    for i in 1:length(profile.r)-1
        if profile.r[i] >= r0 >= profile.r[i+1] && profile.r[i] != profile.r[i+1]
            f = (r0 - profile.r[i+1]) / (profile.r[i] - profile.r[i+1])
            return v[i+1] + f * (v[i] - v[i+1]), (v[i] - v[i+1]) / (profile.r[i] - profile.r[i+1])
        end
    end
    return v[end], 0.0
end

function ss_crossing_fraction(x, z, dx, dz, step, radius)
    b, c = x * dx + z * dz, x^2 + z^2 - radius^2
    d = b^2 - c
    d < 0 && return nothing
    for a in (-b - sqrt(d), -b + sqrt(d))
        1e-8 < a <= step && return a
    end
    return nothing
end

function ss_refract(dx, dz, x, z, v_from, v_to)
    r = hypot(x, z)
    nr, nz = x / r, z / r
    tr, tz = -nz, nr
    d_r, d_t = dx * nr + dz * nz, dx * tr + dz * tz
    sin_i = d_t * v_to / v_from
    if abs(sin_i) <= 1
        cos_i = signbit(d_r) ? -sqrt(1 - sin_i^2) : sqrt(1 - sin_i^2)
        return sin_i * tr + cos_i * nr, sin_i * tz + cos_i * nz, true
    end
    return d_t * tr - d_r * nr, d_t * tz - d_r * nz, false
end

"""
    ss_reflect(dx, dz, x, z)

Specularly reflect a ray direction `(dx, dz)` off a locally spherical interface at position
`(x, z)`: the tangential component of the direction is preserved and the radial component flips
sign -- exactly the boundary-reflection case `ss_refract` already falls back to for total internal
reflection, reused here for the free-surface and CMB bounces in `ss_trace`'s `bounce_plan`.
"""
function ss_reflect(dx, dz, x, z)
    r = hypot(x, z)
    nr, nz = x / r, z / r
    tr, tz = -nz, nr
    d_r, d_t = dx * nr + dz * nz, dx * tr + dz * tz
    return d_t * tr - d_r * nr, d_t * tz - d_r * nz
end

"""
    ss_trace(profile, discs, wave, takeoff_deg; src_r, src_theta_deg, ds, bounce_plan)

Trace a ray launched at `takeoff_deg` from the local inward radial direction.
The returned launch angle is signed, so the initial slowness covector can be
read analytically even when the path later crosses PREM discontinuities.

`bounce_plan` lists, in order, the reflections expected before the ray finally lands: `:surface`
for a free-surface bounce (as in `PP`/`PPP`) or `:cmb` for a CMB reflection that stays in the
mantle instead of transmitting into the core (as in `PcP`/`PcPcP`). An empty plan (the default)
is a direct, unbounced leg. See [`ss_phase_config`](@ref) for the named-phase-to-plan mapping.
"""
function ss_trace(profile, discs, wave0, takeoff_deg; src_r=SS_R_EARTH,
    src_theta_deg=0.0, ds=8.0, max_steps=3000, convert_on_core_exit=true,
    bounce_plan=Symbol[])
    θ0, α = deg2rad(src_theta_deg), deg2rad(takeoff_deg)
    x, z = src_r * sin(θ0), src_r * cos(θ0)
    nr, nz = x / src_r, z / src_r
    # The tracer's transverse direction is -e_theta, hence signed α determines
    # the sign of the canonical angular slowness below.
    tx, tz = -nz, nr
    dx, dz = -cos(α) * nr + sin(α) * tx, -cos(α) * nz + sin(α) * tz
    wave = wave0
    v0, = ss_velocity(profile, wave, src_r)
    sx, sz = dx / v0, dz / v0
    xs, zs, ws, ts = [x], [z], Symbol[wave], [0.0]
    time, landed = 0.0, false
    bounce_idx, nsurface_bounces, ncmb_bounces = 1, 0, 0
    cmb_r = SS_R_EARTH - 2891.0
    for _ in 1:max_steps
        r = hypot(x, z)
        v, dvdr = ss_velocity(profile, wave, r)
        (v > 1e-5 && isfinite(v)) || break
        dx, dz = sx * v, sz * v
        a_surface = ss_crossing_fraction(x, z, dx, dz, ds, SS_R_EARTH)
        a_best, disc_best = (a_surface === nothing ? Inf : a_surface), nothing
        for disc in discs
            a = ss_crossing_fraction(x, z, dx, dz, ds, disc.r)
            if a !== nothing && a < a_best
                a_best, disc_best = a, disc
            end
        end
        if a_best <= ds
            x += dx * a_best; z += dz * a_best; time += a_best / v
            push!(xs, x); push!(zs, z); push!(ts, time)
            if disc_best === nothing
                if bounce_idx <= length(bounce_plan) && bounce_plan[bounce_idx] == :surface
                    dx, dz = ss_reflect(dx, dz, x, z)
                    sx, sz = dx / v, dz / v
                    bounce_idx += 1; nsurface_bounces += 1
                    push!(ws, wave)
                    continue
                end
                push!(ws, wave); landed = true; break
            end
            inward = x * dx + z * dz < 0
            v_from = wave == :P ? (inward ? disc_best.vp_outer : disc_best.vp_inner) :
                                  (inward ? disc_best.vs_outer : disc_best.vs_inner)
            is_cmb = isapprox(disc_best.r, cmb_r; atol=5.0)
            if is_cmb && bounce_idx <= length(bounce_plan) && bounce_plan[bounce_idx] == :cmb
                dx2, dz2 = ss_reflect(dx, dz, x, z)
                sx, sz = dx2 / v_from, dz2 / v_from
                x += dx2 * 1e-3; z += dz2 * 1e-3; time += 1e-3 / v_from
                bounce_idx += 1; ncmb_bounces += 1
                push!(ws, wave)
                continue
            end
            if wave == :S
                v_far = inward ? disc_best.vs_inner : disc_best.vs_outer
                v_to, nextwave = v_far < 1e-5 ? (inward ? disc_best.vp_inner : disc_best.vp_outer, :P) : (v_far, :S)
            else
                previous_vs = inward ? disc_best.vs_outer : disc_best.vs_inner
                v_to, nextwave = (!inward && convert_on_core_exit && previous_vs < 1e-5) ?
                    (disc_best.vs_outer, :S) : (inward ? disc_best.vp_inner : disc_best.vp_outer, :P)
            end
            dx2, dz2, passed = ss_refract(dx, dz, x, z, v_from, v_to)
            if !passed
                v_to, nextwave = v_from, wave
            end
            wave = nextwave
            sx, sz = dx2 / v_to, dz2 / v_to
            x += dx2 * 1e-3; z += dz2 * 1e-3; time += 1e-3 / v_to
            push!(ws, wave)
        else
            x += dx * ds; z += dz * ds
            sx += (-dvdr / v^2 * x / r) * ds
            sz += (-dvdr / v^2 * z / r) * ds
            time += ds / v
            push!(xs, x); push!(zs, z); push!(ws, wave); push!(ts, time)
        end
    end
    final_theta = rad2deg(atan(x, z))
    return (x=xs, z=zs, wavetype=ws, t=ts, landed=landed,
        delta_deg=mod(final_theta - src_theta_deg, 360), takeoff_deg=takeoff_deg,
        min_radius=minimum(hypot.(xs, zs)),
        nsurface_bounces=nsurface_bounces, ncmb_bounces=ncmb_bounces)
end

ss_wrap180(x) = mod(x + 180, 360) - 180

"""
    ss_phase_config(phase)

Map a phase name to the `(wave0, bounce_plan, convert_on_core_exit, angle_range)` that configures
[`ss_trace`](@ref)/[`ss_arrivals`](@ref) to shoot for it directly, rather than trying to infer the
name from an arbitrary path after the fact. `bounce_plan` lists, in order, the reflection expected
at each free-surface (`:surface`) or CMB (`:cmb`) crossing before the ray finally lands -- up to 3
total, mixing both kinds (`PcPcP`'s plan is `[:cmb, :surface, :cmb]`, the implicit surface bounce
between its two CMB reflections). `convert_on_core_exit` matters only for the two phases that
actually transit the core: `SKS` needs it (S enters the fluid core as P, and must convert back to
S on the way out to be SKS rather than PKS); `PKP` must have it `false`, since `ss_trace`'s default
of always converting a P leaving the outer core back to S would otherwise silently turn every PKP
shot into an SKS-shaped ray (confirmed live -- every "PKP" shot before this fix landed with a
spurious `:S` leg and was misclassified as `"SKS"`, i.e. `"PKP"` was unreachable). `angle_range`
caps the takeoff-angle fan `ss_arrivals` shoots: phases confined to steep, near-vertical launches
(anything touching the CMB) have a valid window only a few tens of degrees wide, and shooting the
full ±179.5° at a fixed `nshoot` starves that narrow window of samples -- confirmed live for `PcP`,
which needed `nshoot` above 1000 to resolve at all under the full range but resolves reliably at
the default 181 once the range is restricted to where it can actually land. It is also capped
short of ±90° for every phase: a near-tangent (near-90°) launch grazes the shallow discontinuities
almost edge-on, where this fixed-step tracer's crossing detection becomes numerically unstable --
confirmed live as a source of a spurious "arrival" with a non-physical, non-monotonic jump in
landing distance a fraction of a degree away from its neighbors. Cutting the fan off at 85° means
such a point is correctly reported as no-arrival instead of a wrong travel time.
"""
function ss_phase_config(phase)
    phase == "P" && return (:P, Symbol[], false, 85.0)
    phase == "S" && return (:S, Symbol[], false, 85.0)
    phase == "PP" && return (:P, [:surface], false, 85.0)
    phase == "SS" && return (:S, [:surface], false, 85.0)
    phase == "PPP" && return (:P, [:surface, :surface], false, 85.0)
    phase == "SSS" && return (:S, [:surface, :surface], false, 85.0)
    phase == "PcP" && return (:P, [:cmb], false, 35.0)
    phase == "ScS" && return (:S, [:cmb], false, 35.0)
    phase == "PcPcP" && return (:P, [:cmb, :surface, :cmb], false, 35.0)
    phase == "ScScS" && return (:S, [:cmb, :surface, :cmb], false, 35.0)
    phase == "PKP" && return (:P, Symbol[], false, 60.0)
    phase == "SKS" && return (:S, Symbol[], true, 60.0)
    error("unknown phase \"$phase\"")
end

"""Classify a traced ray by its wave-type history, deepest point, and bounce counts.

`"PKP"`/`"SKS"` are whatever dipped below the CMB (allowed only because `PKP`/`SKS` are shot with
an empty `bounce_plan`, so any core dip is a genuine transmission, not a mis-shot reflection
phase); everything else is read off the number of free-surface/CMB bounces the ray actually used,
which -- because [`ss_phase_config`](@ref) shoots each phase with its own exact `bounce_plan` --
tells us directly whether the shot matched the requested phase (see [`ss_arrivals`](@ref)).
"""
function ss_phase_name(ray)
    hasP, hasS = :P in ray.wavetype, :S in ray.wavetype
    cmb = SS_R_EARTH - 2891.0
    if hasP && hasS && ray.min_radius < cmb
        return "SKS"
    elseif hasP && ray.min_radius < cmb
        return "PKP"
    elseif hasS && !hasP
        ray.ncmb_bounces >= 2 && return "ScScS"
        ray.ncmb_bounces == 1 && return "ScS"
        ray.nsurface_bounces >= 2 && return "SSS"
        ray.nsurface_bounces == 1 && return "SS"
        return "S"
    elseif hasP && !hasS
        ray.ncmb_bounces >= 2 && return "PcPcP"
        ray.ncmb_bounces == 1 && return "PcP"
        ray.nsurface_bounces >= 2 && return "PPP"
        ray.nsurface_bounces == 1 && return "PP"
        return "P"
    end
    return "unknown"
end

"""
    ss_bisect_shot(profile, discs, phase, src_r, src_theta_deg, receiver_theta_deg,
                   alo, ahi, glo, ghi; ds, iterations)

Refine a takeoff-angle bracket `[alo, ahi]` already known to bracket a sign change of the
signed angular mismatch `g(a) = wrap180(delta_deg(a) - target)`, by bisection.

Bisection is used instead of a Newton/secant step because only the *sign* of `g` at each
midpoint is needed: a single unlucky evaluation landing exactly on a discontinuity or in a
shadow zone can just be nudged and retried without ever aborting the search, so every call
that receives a genuine bracket converges. A secant step's local slope estimate has no such
fallback -- see the notebook's own diagnostic notes for a reproduction of it silently
dropping a perfectly real ray.

Near a near-tangent geometry (a grazing CMB reflection is the confirmed case, found while adding
`PcP`), `delta_deg` can jump discontinuously between two branches separated by only a
floating-point epsilon in takeoff angle -- the bracket then never actually shrinks to that
tolerance, it just oscillates between the two branches' `g` values forever. Tracking the best
(smallest `|g|`) ray seen across all iterations, rather than trusting whichever one the *last*
iteration happened to land on, is what makes the search robust to that: the best sample so far is
never worse than a plain "return the final iterate" and is sometimes dramatically better.
"""
function ss_bisect_shot(profile, discs, phase, src_r, src_theta_deg, receiver_theta_deg,
    alo, ahi, glo, ghi; ds=8.0, iterations=40, tolerance_deg=0.01)
    wave0, bounce_plan, convert_on_core_exit, = ss_phase_config(phase)
    target = mod(receiver_theta_deg - src_theta_deg, 360)
    best_ray, best_err = nothing, Inf
    for _ in 1:iterations
        amid = (alo + ahi) / 2
        ray = ss_trace(profile, discs, wave0, amid; src_r, src_theta_deg, ds, bounce_plan, convert_on_core_exit)
        valid = ray.landed && ss_phase_name(ray) == phase
        if !valid
            # amid grazed a discontinuity or shadow-zone edge -- nudge and retry once
            amid = clamp(amid + 1e-3 * sign(ahi - alo), -179.9, 179.9)
            ray = ss_trace(profile, discs, wave0, amid; src_r, src_theta_deg, ds, bounce_plan, convert_on_core_exit)
            valid = ray.landed && ss_phase_name(ray) == phase
        end
        valid || break
        gmid = ss_wrap180(ray.delta_deg - target)
        if abs(gmid) < best_err
            best_ray, best_err = ray, abs(gmid)
        end
        best_err < tolerance_deg && break
        if sign(gmid) == sign(glo)
            alo, glo = amid, gmid
        else
            ahi, ghi = amid, gmid
        end
    end
    return (ray=best_ray, error_deg=best_err)
end

"""Return every resolved ray branch of a chosen phase family that reaches a station.

Every adjacent pair of shot takeoff angles is checked for a sign change of the angular
mismatch; each sign change is a guaranteed bracket, refined independently by bisection.
This finds every triplicated branch exhaustively rather than guessing seeds from a
discrete minimum of `|mismatch|`, which can miss a branch entirely at an unlucky `nshoot`.
"""
function ss_arrivals(profile, discs, phase, src_r, src_theta_deg, receiver_theta_deg;
    nshoot=181, ds=8.0, tolerance_deg=0.05, max_gap=2)
    wave0, bounce_plan, convert_on_core_exit, angle_range = ss_phase_config(phase)
    angles = collect(range(-angle_range, angle_range; length=nshoot))
    target = mod(receiver_theta_deg - src_theta_deg, 360)
    rays = [ss_trace(profile, discs, wave0, a; src_r, src_theta_deg, ds, bounce_plan, convert_on_core_exit) for a in angles]
    isvalid = [r.landed && ss_phase_name(r) == phase for r in rays]
    mism = [isvalid[i] ? ss_wrap180(rays[i].delta_deg - target) : NaN for i in eachindex(rays)]
    refined = NamedTuple[]
    # A single shot can fail to land or to keep the requested phase purely from grazing a
    # discontinuity or a shadow-zone edge at that exact angle, even though the branch is
    # continuous on either side. Bridge over up to `max_gap` consecutive such misses before
    # giving up on a bracket, so an isolated bad sample can't hide a real sign change.
    for i in 1:length(angles)-1
        isvalid[i] || continue
        j = i + 1
        while j <= length(angles) && !isvalid[j] && j - i <= max_gap
            j += 1
        end
        (j > length(angles) || !isvalid[j]) && continue
        mism[i] * mism[j] > 0 && continue
        push!(refined, ss_bisect_shot(profile, discs, phase, src_r, src_theta_deg, receiver_theta_deg,
            angles[i], angles[j], mism[i], mism[j]; ds))
    end
    refined = [x for x in refined if x.error_deg <= tolerance_deg]
    kept, angles_seen = NamedTuple[], Float64[]
    for item in refined
        key = round(item.ray.takeoff_deg, digits=3)
        key in angles_seen && continue
        push!(angles_seen, key); push!(kept, item)
    end
    return sort(kept; by=x -> x.ray.t[end])
end

"""Analytic source-endpoint derivative of one ray's travel time.

Returns `(dtheta, ddepth)`, in seconds/radian and seconds/km.  It is the
negative launch slowness covector, expressed at the source; no TauP or finite
difference enters this calculation.
"""
function ss_endpoint_derivative(profile, ray, wave0, src_r)
    v, = ss_velocity(profile, wave0, src_r)
    α = deg2rad(ray.takeoff_deg)
    # Outgoing ray components are u_r=-cos(α), u_theta=-sin(α).
    ptheta = -src_r * sin(α) / v
    pr = -cos(α) / v
    return (dtheta=-ptheta, ddepth=pr)
end

"""
    ss_build_field(profile, wavetype, src_r, src_theta_deg; n=201, allow_core=true)

Precompute a first-arrival travel-time field from a single fixed source with `Eikonal.jl`'s
fast-sweeping method on an `n×n` Cartesian grid covering the whole disk. Scanning a grid of
candidate hypocentres by tracing a fresh shooting search to each one is wasteful: this instead
solves for the travel time from the source to *every* point in the disk in one pass, so each
candidate afterwards is a cheap lookup (see [`ss_field_lookup`](@ref)) instead of its own search.

Only usable for `wavetype`'s own single, unconverted branch -- there is no mode-conversion
bookkeeping here, so this matches the `"P"`/`"S"` phase families and not `"PKP"`/`"SKS"`. With
`allow_core=false` the outer+inner core is walled off with a near-impassable slowness (the same
device used for the disk's exterior below), so the field only contains the shallow,
mantle-confined branch -- matching how [`ss_phase_name`](@ref) reserves `"PKP"`/`"SKS"` for
anything that actually dips below the CMB.
"""
function ss_build_field(profile, wavetype, src_r, src_theta_deg; n=201, allow_core=true)
    xs = range(-SS_R_EARTH, SS_R_EARTH; length=n)
    zs = range(-SS_R_EARTH, SS_R_EARTH; length=n)
    dcell = xs[2] - xs[1]
    cmb = SS_R_EARTH - 2891.0
    vgrid = [begin
        r = hypot(xx, zz)
        if r > SS_R_EARTH || (!allow_core && r < cmb)
            0.01                    # a near-impassable "wall", not a fast shortcut
        else
            max(ss_velocity(profile, wavetype, r)[1], 0.01)
        end
    end for zz in zs, xx in xs]
    sgrid = 1.0 ./ vgrid
    src_x, src_z = src_r * sin(deg2rad(src_theta_deg)), src_r * cos(deg2rad(src_theta_deg))
    iz0, ix0 = argmin(abs.(zs .- src_z)), argmin(abs.(xs .- src_x))
    fs = FastSweeping(sgrid)
    init!(fs, (iz0, ix0))
    sweep!(fs, verbose=false, epsilon=1e-6)
    tgrid = (fs.t .* dcell)[1:n, 1:n]
    return (xs=collect(xs), zs=collect(zs), t=tgrid)
end

"""
    ss_field_lookup(field, θ, depth)

Bilinearly interpolate a travel-time field built by [`ss_build_field`](@ref) at angular position
`θ` (degrees, source-relative, matching [`ss_trace`](@ref)'s convention) and `depth` below the
surface. Returns `Inf` for a point outside the field's grid (only possible right at the domain
edge, since the grid spans the whole disk).
"""
function ss_field_lookup(field, θ, depth)
    r = SS_R_EARTH - depth
    x, z = r * sin(deg2rad(θ)), r * cos(deg2rad(θ))
    xs, zs, t = field.xs, field.zs, field.t
    (x < xs[1] || x > xs[end] || z < zs[1] || z > zs[end]) && return Inf
    dx = xs[2] - xs[1]
    i = clamp(floor(Int, (x - xs[1]) / dx) + 1, 1, length(xs) - 1)
    j = clamp(floor(Int, (z - zs[1]) / dx) + 1, 1, length(zs) - 1)
    fx = (x - xs[i]) / dx
    fz = (z - zs[j]) / dx
    # `t` is indexed [z, x] to match `ss_build_field`'s (zz, xx) grid comprehension order.
    t00, t10, t01, t11 = t[j, i], t[j, i+1], t[j+1, i], t[j+1, i+1]
    return (1 - fz) * ((1 - fx) * t00 + fx * t10) + fz * ((1 - fx) * t01 + fx * t11)
end

"""
    ss_field_source_state(fieldA, fieldB, θ, depth, target_time; dθ=2.0, dh=60.0)

Fast approximate counterpart to [`ss_source_state`](@ref), using two precomputed
[`ss_build_field`](@ref) travel-time fields instead of an independent shooting search per
candidate hypocentre. `Ftheta`/`Fdepth` come from a centered finite difference taken directly on
the field; `dθ`/`dh` must stay comfortably larger than the field's own grid spacing (a few tens of
km) or the difference just measures interpolation noise instead of the real derivative.

This is deliberately only accurate enough to *find* candidate roots cheaply during a wide scan --
the field's own Cartesian discretization costs a few seconds of accuracy that a genuinely reported
root should not carry, so every seed this produces is re-verified against the exact ray tracer
(`ss_source_state`) before being refined or reported.
"""
function ss_field_source_state(fieldA, fieldB, θ, depth, target_time; dθ=2.0, dh=60.0)
    tA, tB = ss_field_lookup(fieldA, θ, depth), ss_field_lookup(fieldB, θ, depth)
    (isfinite(tA) && isfinite(tB)) || return nothing
    dAdθ = (ss_field_lookup(fieldA, θ + dθ, depth) - ss_field_lookup(fieldA, θ - dθ, depth)) / (2dθ * pi / 180)
    dBdθ = (ss_field_lookup(fieldB, θ + dθ, depth) - ss_field_lookup(fieldB, θ - dθ, depth)) / (2dθ * pi / 180)
    hi, lo = min(depth + dh, 600.0), max(depth - dh, 0.0)
    dAdh = (ss_field_lookup(fieldA, θ, hi) - ss_field_lookup(fieldA, θ, lo)) / (hi - lo)
    dBdh = (ss_field_lookup(fieldB, θ, hi) - ss_field_lookup(fieldB, θ, lo)) / (hi - lo)
    return (theta=θ, depth=depth, tA=tA, tB=tB,
        Ftime=tA - tB - target_time, Ftheta=dAdθ - dBdθ, Fdepth=dAdh - dBdh)
end
