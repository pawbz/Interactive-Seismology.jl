### A Pluto.jl notebook ###
# v0.2.6

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

# ╔═╡ 9a802663-864a-454c-9e36-e5ee2cf60e87
using PlutoUI

# ╔═╡ f32a29dc-ed47-4d1e-8ea7-efb287f82c0b
using Roots

# ╔═╡ 3b7fbabb-c4c6-4b3e-9d3a-13f0ac2a7d67
using LinearAlgebra

# ╔═╡ 6e3d880e-82d5-45b0-a823-b2df7d8be6f0
TableOfContents(include_definitions=true)

# ╔═╡ c3a8bd5e-41c3-46f4-bfb5-773eadfb67ee
md"""
# Rayleigh Wave Dispersion Curves
This notebook provides an interactive environment to explore the dispersion of Rayleigh waves, a fundamental type of surface seismic wave used in seismology and Earth structure studies.

Rayleigh waves involve both vertical and horizontal ground motion and are sensitive to the elastic properties and layering of the crust and upper mantle. Their dispersion (variation of phase velocity with period) reveals information about subsurface structure.

Here, you can:

- Build a layered Earth model interactively by adjusting thickness, P-wave velocity (Vp), S-wave velocity (Vs), and density (ρ) for each layer;
- Compute and visualize Rayleigh-wave dispersion curves, observing how changes in the model affect the shape of the curve and the inferred Earth properties.


##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ fa019722-6bca-11f1-9c8d-9ffc56627241
md"""
### Observed dispersion overlay

Phase file: $(@bind observed_phase_file FilePicker())

Group file: $(@bind observed_group_file FilePicker())
"""

# ╔═╡ b21471fd-f2f7-4715-b354-48ddd7aa0207
md"---"

# ╔═╡ 5c136747-07b4-454e-9a75-0998a24560df
struct Layer
    thickness::Float64   # km
    vp::Float64          # km/s
    vs::Float64          # km/s
    rho::Float64         # g/cm³
    Qp::Float64          # P-wave quality factor
    Qs::Float64          # S-wave quality factor
end

# ╔═╡ 48f94773-cfea-46d5-b341-ae57fbe43e55
struct EarthModel
    layers::Vector{Layer}
end

# ╔═╡ 021afa55-3f7d-4c5e-8768-c91ae697c653
"""
    varmat(p, q, ra, rb, wvno, xka, xkb, dpth)

Evaluate the trigonometric/hyperbolic building blocks (`cos(p)`, `sin(p)/ra`, ...) that
[`dnka`](@ref) combines into one layer's Dunkin compound matrix, switching automatically
between the propagating (`wvno < xka`/`xkb`) and evanescent forms and tracking the
exponent `exa` separately so large `p`/`q` don't overflow before [`normc`](@ref) rescales.
"""
function varmat(p, q, ra, rb, wvno, xka, xkb, dpth)
    # Returns a NamedTuple of all needed products for Dunkin matrix, plus exponent exa
    a0 = 1.0
    pex = 0.0
    sex = 0.0

    # P-wave
    if wvno < xka
        sinp = sin(p)
        w = sinp / ra
        x = -ra * sinp
        cosp = cos(p)
    elseif wvno == xka
        cosp = 1.0
        w = dpth
        x = 0.0
    else
        pex = p
        fac = p < 16 ? exp(-2 * p) : 0.0
        cosp = (1.0 + fac) * 0.5
        sinp = (1.0 - fac) * 0.5
        w = sinp / ra
        x = ra * sinp
    end

    # S-wave
    if wvno < xkb
        sinq = sin(q)
        y = sinq / rb
        z = -rb * sinq
        cosq = cos(q)
    elseif wvno == xkb
        cosq = 1.0
        y = dpth
        z = 0.0
    else
        sex = q
        fac = q < 16 ? exp(-2 * q) : 0.0
        cosq = (1.0 + fac) * 0.5
        sinq = (1.0 - fac) * 0.5
        y = sinq / rb
        z = rb * sinq
    end

    # Exponent for normalization
    exa = pex + sex
    a0 = exa < 60.0 ? exp(-exa) : 0.0

    # Products
    cpcq = cosp * cosq
    cpy = cosp * y
    cpz = cosp * z
    cqw = cosq * w
    cqx = cosq * x
    xy = x * y
    xz = x * z
    wy = w * y
    wz = w * z

    # S-wave normalization for evanescent case
    qmp = sex - pex
    fac = qmp > -40.0 ? exp(qmp) : 0.0
    cosq *= fac
    y *= fac
    z *= fac

    return (; a0, cpcq, cpy, cpz, cqw, cqx, xy, xz, wy, wz, exa)
end

# ╔═╡ c6e30458-2eff-4a66-88de-f75cdbfcd9a4
"""
    normc(ee)

Rescale the propagator vector `ee` by its own peak magnitude, returning the rescaled
vector and `log` of the scale factor removed. [`rayleigh_secular`](@ref) calls this after
every layer to keep the Thomson-Haskell matrix product from overflowing through very thick
or high-frequency layers — the same technique this repo's Love-wave eigenfunction shooting
uses (shooting from the far boundary) to sidestep the same instability from the other end.
"""
function normc(ee)
    t1 = maximum(abs.(ee))
    if t1 < 1e-40
        t1 = 1.0
    end
    ee ./= t1
    ex = log(t1)
    return ee, ex
end

# ╔═╡ 6930149b-0a3a-453f-bf15-8d581f62c433
"""
    dnka(wvno2, gam, gammk, rho, var)

Build one layer's 5×5 Dunkin compound (minor) matrix from the trigonometric products in
`var` (see [`varmat`](@ref)). Working with the compound matrix — minors of the raw P-SV
propagator, rather than the propagator itself — is the classic Dunkin (1965) fix for the
same Thomson-Haskell overflow problem [`normc`](@ref) also guards against.
"""
function dnka(wvno2, gam, gammk, rho, var)
    ca = zeros(5, 5)
    one = 1.0
    two = 2.0

    a0 = var[:a0]
    cpcq = var[:cpcq]
    cpy = var[:cpy]
    cpz = var[:cpz]
    cqw = var[:cqw]
    cqx = var[:cqx]
    xy = var[:xy]
    xz = var[:xz]
    wy = var[:wy]
    wz = var[:wz]

    gamm1 = gam - one
    twgm1 = gam + gamm1
    gmgmk = gam * gammk
    gmgm1 = gam * gamm1
    gm1sq = gamm1 * gamm1
    rho2 = rho * rho
    a0pq = a0 - cpcq

    ca[1, 1] = cpcq - two * gmgm1 * a0pq - gmgmk * xz - wvno2 * gm1sq * wy
    ca[1, 2] = (wvno2 * cpy - cqx) / rho
    ca[1, 3] = -(twgm1 * a0pq + gammk * xz + wvno2 * gamm1 * wy) / rho
    ca[1, 4] = (cpz - wvno2 * cqw) / rho
    ca[1, 5] = -(two * wvno2 * a0pq + xz + wvno2 * wvno2 * wy) / rho2

    ca[2, 1] = (gmgmk * cpz - gm1sq * cqw) * rho
    ca[2, 2] = cpcq
    ca[2, 3] = gammk * cpz - gamm1 * cqw
    ca[2, 4] = -wz
    ca[2, 5] = ca[1, 4]

    ca[4, 1] = (gm1sq * cpy - gmgmk * cqx) * rho
    ca[4, 2] = -xy
    ca[4, 3] = gamm1 * cpy - gammk * cqx
    ca[4, 4] = ca[2, 2]
    ca[4, 5] = ca[1, 2]

    ca[5, 1] = -(two * gmgmk * gm1sq * a0pq + gmgmk * gmgmk * xz +
                 gm1sq * gm1sq * wy) * rho2
    ca[5, 2] = ca[4, 1]
    ca[5, 3] = -(gammk * gamm1 * twgm1 * a0pq + gam * gammk * gammk * xz +
                 gamm1 * gm1sq * wy) * rho
    ca[5, 4] = ca[2, 1]
    ca[5, 5] = ca[1, 1]

    t = -two * wvno2
    ca[3, 1] = t * ca[5, 3]
    ca[3, 2] = t * ca[4, 3]
    ca[3, 3] = a0 + two * (cpcq - ca[1, 1])
    ca[3, 4] = t * ca[2, 3]
    ca[3, 5] = t * ca[1, 3]

    return ca
end

# ╔═╡ 369f3691-78d3-464c-a491-a293f58df913
"""
    rayleigh_secular(layers::Vector{Layer}, period::Float64, c::Float64)

Evaluate the Rayleigh secular function for a stack of layers at a given period and trial phase velocity c.
Returns the value whose zero crossing gives the Rayleigh phase velocity.
"""
function rayleigh_secular(layers::Vector{Layer}, ω::Float64, c::Float64)

    wvno = ω / c
    wvno2 = wvno^2

    # Prepare bottom half-space parameters
    bottom = layers[end]
    α = bottom.vp
    β = bottom.vs
    ρ = bottom.rho

    xka = ω / α
    xkb = ω / β

    # Compute vertical slownesses for P and S in half-space
    wvnop_a = wvno + xka
    wvnom_a = abs(wvno - xka)
    ra = sqrt(wvnop_a * wvnom_a)

    wvnop_b = wvno + xkb
    wvnom_b = abs(wvno - xkb)
    rb = sqrt(wvnop_b * wvnom_b)

    t = β / ω
    gammk = 2.0 * t^2
    gam = gammk * wvno2
    gamm1 = gam - 1.0

    # E vector for bottom half-space
    e = zeros(5)
    e[1] = ρ^2 * (gamm1^2 - gam * gammk * ra * rb)
    e[2] = -ρ * ra
    e[3] = ρ * (gamm1 - gammk * ra * rb)
    e[4] = ρ * rb
    e[5] = wvno2 - ra * rb
    # @show e
    # Upward matrix multiplication through layers
    for m in length(layers)-1:-1:1
        lay = layers[m]
        α = lay.vp
        β = lay.vs
        ρ = lay.rho
        dpth = lay.thickness

        xka = ω / α
        xkb = ω / β

        wvnop_a = wvno + xka
        wvnom_a = abs(wvno - xka)
        ra = sqrt(wvnop_a * wvnom_a)

        wvnop_b = wvno + xkb
        wvnom_b = abs(wvno - xkb)
        rb = sqrt(wvnop_b * wvnom_b)

        p = ra * dpth
        q = rb * dpth

        # Compute extended trigonometric/hyperbolic functions
        var = varmat(p, q, ra, rb, wvno, xka, xkb, dpth)
        t = β / ω
        gammk = 2.0 * t^2
        gam = gammk * wvno2
        ca = dnka(wvno2, gam, gammk, ρ, var)

        # Matrix multiplication: e_new = e * ca
        ee = zeros(5)
        for i in 1:5
            for j in 1:5
                ee[i] += e[j] * ca[j, i]
            end
        end
        # Normalize to avoid overflow/underflow
        e, _ = normc(ee)
    end
    # If top layer is water, apply water layer correction (not shown here, but implemented in CPS)
    # For most crustal models, just return e[1]
    return e[1]
end

# ╔═╡ 87a1ea4f-e41c-44f0-8712-6dad46b1f2a3
"""
    DispersionResult(period, phase_velocities)

All Rayleigh phase-velocity roots found at one `period` (there can be more than one —
higher modes as well as the fundamental).
"""
struct DispersionResult
    period::Float64
    phase_velocities::Vector{Float64}
end

# ╔═╡ c6187aa6-a890-4b35-8fc8-ef66373a4d2e
"""
    solve_phase_velocity_rayleigh(layers, T; cmin=nothing, cmax=nothing)

Find every Rayleigh phase-velocity root at period `T` by scanning [`rayleigh_secular`](@ref)
for sign changes between `cmin` (half the slowest Vs) and `cmax` (1.2× the fastest Vp).
"""
function solve_phase_velocity_rayleigh(layers::Vector{Layer}, T::Float64; cmin=nothing, cmax=nothing)
    ω = 2π / T
    vs_vals = [L.vs for L in layers]
    vp_vals = [L.vp for L in layers]

    if cmin === nothing
        cmin = minimum(vs_vals) * 0.5
    end
    if cmax === nothing
        cmax = maximum(vp_vals) * 1.2
    end

    f(c) = rayleigh_secular(layers, ω, c)

    roots = find_zeros(f, cmin, cmax)
    return DispersionResult(T, roots)
end

# ╔═╡ 54effe7a-7d35-40e1-a861-4fa7f6ec9c06
"""
    compute_rayleigh_dispersion(layers, periods) -> Vector{DispersionResult}

Solve [`solve_phase_velocity_rayleigh`](@ref) at every period, independently — this is the
per-period loop `res` (and the widget's live dispersion-curve panel) is built from.
"""
function compute_rayleigh_dispersion(layers::Vector{Layer}, periods::Vector{Float64})
    R = []
    for period in periods
        r = solve_phase_velocity_rayleigh(layers, period)
        push!(R, r)
    end
    return R
end

# ╔═╡ e55ee346-de7c-4b88-88f2-cd1afbbc73a2
md"""
Maximum period: $(@bind period_max Slider(20.0:5.0:200.0, default=100.0, show_value=true)) s
"""

# ╔═╡ 3e02c803-9bc3-4123-a4d8-ab07021a0062
periods = collect(range(5.0, stop = period_max, length=64))

# ╔═╡ 7a1c3608-a871-424f-9455-62350343f906
begin
	struct ObservedDispersion
	    periods::Vector{Float64}
	    phase_velocities::Vector{Float64}
	    group_velocities::Vector{Float64}
	end
	
	empty_observed_dispersion() = ObservedDispersion(Float64[], Float64[], Float64[])
	load_observed_dispersion(::Missing; kind::Symbol=:auto) = empty_observed_dispersion()
	load_observed_dispersion(::Nothing; kind::Symbol=:auto) = empty_observed_dispersion()
	
	function _parse_observed_number(token::AbstractString)
	    t = strip(token)
	    isempty(t) && return NaN
	    lowercase(t) in ("nan", "missing", "none") && return NaN
	    return tryparse(Float64, t) === nothing ? NaN : parse(Float64, t)
	end
	
	function _numeric_rows_from_string(text::AbstractString)
	    rows = Vector{Vector{Float64}}()
	    for line in split(text, '\n')
	        body = strip(first(split(first(split(line, "#")), "//")))
	        isempty(body) && continue
	        tokens = split(replace(body, ',' => ' '))
	        nums = [_parse_observed_number(tok) for tok in tokens]
	        length(nums) >= 2 && push!(rows, nums)
	    end
	    return rows
	end
	
	function _observed_rows_from_path(path::AbstractString)
		path = expanduser(strip(path))
		isempty(path) && return Vector{Vector{Float64}}()
		isfile(path) || error("Observed dispersion file not found: $path")
		return _numeric_rows_from_string(read(path, String))
	end
	
	function _observed_rows_from_filepicker(file)
		data = if haskey(file, "data")
			file["data"]
		elseif haskey(file, :data)
			file[:data]
		else
			nothing
		end
		data === nothing && return Vector{Vector{Float64}}()
		text = data isa AbstractVector{UInt8} ? String(data) : String(data)
		return _numeric_rows_from_string(text)
	end
	
	function _dispersion_from_rows(rows; kind::Symbol=:auto)
		isempty(rows) && return empty_observed_dispersion()
	    period = Float64[]
	    phase = Float64[]
	    group = Float64[]
	    for row in rows
	        T = row[1]
	        isfinite(T) && T > 0 || continue
	
	        if length(row) >= 6
	            # pDSurfTomo flat file: period lat1 lon1 lat2 lon2 velocity.
	            v = row[6]
	            if kind == :group
	                push!(period, T); push!(phase, NaN); push!(group, v)
	            else
	                push!(period, T); push!(phase, v); push!(group, NaN)
	            end
	        elseif length(row) >= 3
	            # Simple combined table: period phase_velocity group_velocity.
	            push!(period, T); push!(phase, row[2]); push!(group, row[3])
	        else
	            v = row[2]
	            if kind == :group
	                push!(period, T); push!(phase, NaN); push!(group, v)
	            else
	                push!(period, T); push!(phase, v); push!(group, NaN)
	            end
	        end
	    end
	    return ObservedDispersion(period, phase, group)
	end
	
	function load_observed_dispersion(path::AbstractString; kind::Symbol=:auto)
		return _dispersion_from_rows(_observed_rows_from_path(path); kind=kind)
	end
	
	function load_observed_dispersion(file::AbstractDict; kind::Symbol=:auto)
		return _dispersion_from_rows(_observed_rows_from_filepicker(file); kind=kind)
	end
	
	function merge_observed_dispersion(phase_obs::ObservedDispersion, group_obs::ObservedDispersion)
	    return ObservedDispersion(
	        vcat(phase_obs.periods, group_obs.periods),
	        vcat(phase_obs.phase_velocities, group_obs.phase_velocities),
	        vcat(phase_obs.group_velocities, group_obs.group_velocities),
	    )
	end
	
	"""
	    track_fundamental_mode(periods, all_roots) -> Vector{Float64}

	Follow the fundamental mode continuously across periods, instead of independently
	taking the smallest root at each period. A strong low-velocity zone (or any strong
	layer contrast) can pack many closely-spaced modes into the search interval, and
	picking "whichever root is smallest" at each period can jump between different
	physical mode branches from one period to the next — producing a jagged phase-velocity
	curve that then gets amplified into a large spurious spike by
	[`group_velocity_from_phase`](@ref)'s finite-difference.

	Starting from the longest period (modes are fewest and best separated there —
	confirmed empirically: root count grows sharply toward short period as a low-velocity
	zone traps more higher modes) and walking toward the shortest, each step picks
	whichever of the current period's roots is closest to the previous period's tracked
	value, so the returned curve follows one continuous branch instead of jumping between
	them.
	"""
	function track_fundamental_mode(periods::AbstractVector, all_roots::AbstractVector)
	    n = length(periods)
	    tracked = fill(NaN, n)
	    order = sortperm(periods; rev=true)
	    prev = NaN
	    for idx in order
	        roots = all_roots[idx]
	        isempty(roots) && continue
	        tracked[idx] = isnan(prev) ? minimum(roots) : roots[argmin(abs.(roots .- prev))]
	        prev = tracked[idx]
	    end
	    return tracked
	end

	"""
	    fundamental_phase_velocities(periods, results) -> Vector{Float64}

	Extract the fundamental-mode phase-velocity curve from `results` (one
	[`DispersionResult`](@ref) per period, each possibly holding several mode roots) via
	[`track_fundamental_mode`](@ref).
	"""
	function fundamental_phase_velocities(periods, results)
	    all_roots = [r.phase_velocities for r in results]
	    return track_fundamental_mode(periods, all_roots)
	end

	"""
	    group_velocity_from_phase(periods, phase_velocities) -> Vector{Float64}

	Group velocity `U = c / (1 + (T/c)(dc/dT))` from a two-point finite difference of an
	*already mode-tracked* phase-velocity curve (see [`track_fundamental_mode`](@ref)).

	A steep enough local slope (a strong low-velocity zone can produce one) can still push
	the denominator `1 + (T/c)(dc/dT)` through zero between two sample periods — a real
	near-singularity (physically, an Airy-phase-like extremum) that this coarse two-point
	stencil sees only as a sign flip to a huge spurious value, not the true smooth peak.
	Rather than plot that, points landing outside a generous physical band (positive, and
	no more than 1.5x the curve's own peak phase velocity) are dropped back to `NaN` — an
	honest gap in an under-resolved region instead of a wrong-looking spike.
	"""
	function group_velocity_from_phase(periods::AbstractVector, phase_velocities::AbstractVector)
	    U = fill(NaN, length(periods))
	    finite_c = filter(isfinite, phase_velocities)
	    isempty(finite_c) && return U
	    cmax = maximum(finite_c)
	    for i in eachindex(periods)
	        isfinite(phase_velocities[i]) || continue
	        left = max(1, i - 1)
	        right = min(length(periods), i + 1)
	        left == right && continue
	        isfinite(phase_velocities[left]) && isfinite(phase_velocities[right]) || continue
	        dc_dT = (phase_velocities[right] - phase_velocities[left]) / (periods[right] - periods[left])
	        denom = 1 + periods[i] * dc_dT / phase_velocities[i]
	        isfinite(denom) && denom != 0 || continue
	        u = phase_velocities[i] / denom
	        0 < u <= 1.5 * cmax || continue
	        U[i] = u
	    end
	    return U
	end
end

# ╔═╡ 6ce1ef3b-2e7f-44db-9b74-1663bb6a84af
"""
    vertical_slowness(v, p)

Vertical slowness for a wave of horizontal slowness `p` in a medium of velocity `v`:
the real, negative-root branch when propagating (`p < 1/v`), or the purely-imaginary,
positive-imaginary-part branch when evanescent (`p > 1/v`, the usual case for every layer
above a trapped Rayleigh mode's half-space). Adapted from
`src/experimental/wavenumber-integration.jl`'s function of the same name.
"""
function vertical_slowness(v, p)
    arg = (1 / abs2(v)) - abs2(p)
    return arg >= 0 ? -sqrt(arg) : 1im * sqrt(-arg)
end

# ╔═╡ e54679ad-6c92-43db-b933-8a88bbafdd89
"""
    wavefield_components(p, qP, qS, λ, μ, mode, direction)

The `[u_x, u_z, σ_xz, σ_zz]` eigenvector for one P or SV potential traveling `:up` or
`:down`, at horizontal slowness `p`. One column of [`modal_basis_psv`](@ref)'s modal
basis `V`. Adapted from `src/experimental/wavenumber-integration.jl`'s function of the
same name (dropping nothing — the formulas themselves don't involve attenuation).
"""
function wavefield_components(p, qP, qS, λ, μ, mode::Symbol, direction::Symbol)
    sign_ = direction == :down ? 1.0 : -1.0
    if mode == :P
        ux = 1im * p
        uz = 1im * sign_ * qP
        sxz = -2 * μ * p * qP * sign_
        szz = -(λ * (p^2 + qP^2) + 2 * μ * qP^2)
    else
        ux = -1im * sign_ * qS
        uz = 1im * p
        sxz = μ * (qS^2 - p^2)
        szz = -2 * μ * sign_ * p * qS
    end
    return ux, uz, sxz, szz
end

# ╔═╡ 195d2bd4-7828-4a98-980c-dd8bc3e4345d
"""
    modal_basis_psv(layer::Layer, p, ω) -> (V, qP, qS)

The 4×4 P-SV modal eigenbasis for one layer: columns are the `[u_x, u_z, σ_xz, σ_zz]`
eigenvectors of the down-going P, up-going P, down-going SV, and up-going SV plane waves
(in that order), from [`wavefield_components`](@ref). [`layer_matrix_PSV`](@ref)
diagonalizes propagation through the layer in this basis. Adapted from
`src/experimental/wavenumber-integration.jl`'s function of the same name, dropping its
`Qp`/`Qs` causal-velocity dispersion (this notebook doesn't model attenuation).
"""
function modal_basis_psv(layer::Layer, p, ω)
    vp = complex(layer.vp)
    vs = complex(layer.vs)
    rho = layer.rho
    qP = vertical_slowness(vp, p)
    qS = vertical_slowness(vs, p)
    λ, μ = rho * (vp^2 - 2vs^2), rho * vs^2
    col_Pd = collect(wavefield_components(p, qP, qS, λ, μ, :P, :down))
    col_Pu = collect(wavefield_components(p, qP, qS, λ, μ, :P, :up))
    col_Sd = collect(wavefield_components(p, qP, qS, λ, μ, :SV, :down))
    col_Su = collect(wavefield_components(p, qP, qS, λ, μ, :SV, :up))
    V = hcat(col_Pd, col_Pu, col_Sd, col_Su)
    return V, qP, qS
end

# ╔═╡ 984ff052-a12b-4718-9ede-611666ccea76
"""
    one_way_phase(q, ω, h; clip=80.0)

The one-way propagation factor `exp(iω·q·h)` for vertical slowness `q` across thickness
`h`, with the evanescent decay exponent clamped to `±clip` to avoid overflow for thick or
high-frequency layers. Adapted from `src/experimental/wavenumber-integration.jl`'s
function of the same name.
"""
function one_way_phase(q, ω::Float64, h::Float64; clip=80.0)
    if !isfinite(h)
        return one(q)
    end
    atten = clamp(-ω * imag(q) * h, -clip, clip)
    phase = exp(1im * ω * real(q) * h)
    return phase * exp(atten)
end

# ╔═╡ 838fcc3e-3c3f-4061-866b-8966e263498c
"""
    layer_matrix_PSV(layer::Layer, p, ω)

The 4×4 complex P-SV propagator `D` for one layer, mapping
`[u_x, u_z, σ_xz, σ_zz]` at the top to the same at the bottom: `D = V·diag(phase)·V⁻¹`,
built by diagonalizing propagation in [`modal_basis_psv`](@ref)'s eigenbasis. The direct
4-component analogue of [`layer_matrix_SH`](@ref) in this repo's Love-wave notebook (same
role: a per-layer transfer matrix used to shoot an eigenfunction through the stack), just
built via eigen-decomposition instead of a closed-form trig matrix, since the P-SV system
mixes two wave types where SH only has one. Passing a *negative* thickness (as
[`compute_rayleigh_eigenfunctions`](@ref) does to shoot upward) gives `D`'s own inverse
exactly, since the phase factors for `+h` and `-h` are exact reciprocals.
"""
function layer_matrix_PSV(layer::Layer, p, ω)
    h = layer.thickness
    V, qP, qS = modal_basis_psv(layer, p, ω)
    phase_Pd = one_way_phase(qP, ω, h)
    phase_Pu = one_way_phase(qP, ω, -h)
    phase_Sd = one_way_phase(qS, ω, h)
    phase_Su = one_way_phase(qS, ω, -h)
    Phase = Diagonal([phase_Pd, phase_Pu, phase_Sd, phase_Su])
    return V * Phase * inv(V)
end

# ╔═╡ ac3558dc-5882-4898-a8de-543e9483af11
"""
    compute_rayleigh_eigenfunctions(layers, T, c; samples_per_layer=25)

Compute the Rayleigh-wave radial (`ur`) and vertical (`uz`) displacement eigenfunction
depth profiles for a solved period `T` and phase velocity `c`, by shooting **upward**
from the half-space to the free surface — the P-SV generalization of this notebook's
sibling Love-wave notebook's `compute_love_eigenfunctions` (same shooting direction, same
reason: starting from the half-space's own decaying solution avoids ever admixing the
numerically-dominant growing branch that downward shooting would accumulate).

Unlike [`rayleigh_secular`](@ref)'s Dunkin compound-matrix formulation (which solves for
`c` robustly but never exposes displacement/stress directly — its `e` vector is a
5-component vector of *minors*, not the state itself), this uses the raw 4-component
complex state `[u_x, u_z, σ_xz, σ_zz]` and an explicit per-layer P-SV modal-eigenbasis
propagator (`V·diag(phase)·V⁻¹`, adapted from this repo's
`src/experimental/wavenumber-integration.jl`), which *does* expose displacement and
stress directly. This is the standard trade-off in surface-wave eigenfunction codes: the
compound-matrix method for finding `c`, then a separate direct-propagator pass for the
eigenfunctions once `c` is known.

The half-space alone has a 2-dimensional space of valid (decaying) states — any
combination of its down-going P and SV modes — so unlike Love's single-mode half-space
(pinned by one ratio), the specific combination here is found by propagating that whole
2D subspace up to the surface and solving for the linear combination whose accumulated
stress is exactly zero there (the same "find `c` such that a nontrivial solution exists"
condition [`rayleigh_secular`](@ref) already solves via its determinant/minors, arrived at
here directly instead via a 2×2 null-space solve).

A free (non-attenuating) Rayleigh mode has `ur` and `uz` in *exact quadrature* (90° phase
apart) at every depth — confirmed both by the free-surface residual and by an independent
closed-form cross-check against this notebook's sibling `Rayleigh-function.jl` for a
uniform half-space. The arbitrary overall complex phase is fixed by rotating `uz` real at
the surface (matching CPS's own `UR`/`UZ` convention); `ur` then comes out purely
imaginary, so its imaginary part — not real part, which would trivially vanish at this
reference phase — is the physically meaningful shape.

Returns a named tuple: `depths`, peak-normalized `ur`/`uz` shapes, the kinetic energy
integral `I1 = ∫ρ(ur²+uz²)dz` (needed for synthetic-seismogram amplitude scaling, the
2-component generalization of Love's `∫ρu²dz`), and the free-surface stress-residual
magnitudes `sxz_resid`/`szz_resid` (both should be tiny — a check that `c` was an accurate
root).
"""
function compute_rayleigh_eigenfunctions(layers::Vector{Layer}, T::Float64, c::Float64;
        samples_per_layer::Int=25)
    ω = 2π / T
    p = 1.0 / c
    finite_layers = layers[1:end-1]
    hs = layers[end]

    Vhs, _, _ = modal_basis_psv(hs, p, ω)
    basis = Vhs[:, [1, 3]]   # decaying (down-going) P and SV columns

    ztot = sum(l.thickness for l in finite_layers; init=0.0)

    # Per-microstep propagators, deepest layer first (nearest the half-space) --
    # built once, reused both to find the half-space's amplitude ratio and to
    # sample the eigenfunction on the way back up.
    mats = Matrix{ComplexF64}[]
    for layer in reverse(finite_layers)
        n = max(4, samples_per_layer)
        dz = layer.thickness / n
        sub_up = Layer(-dz, layer.vp, layer.vs, layer.rho, layer.Qp, layer.Qs)
        M = layer_matrix_PSV(sub_up, p, ω)
        for _ in 1:n
            push!(mats, M)
        end
    end

    Mtotal = Matrix{ComplexF64}(I, 4, 4)
    for M in mats
        Mtotal = M * Mtotal
    end

    surf_basis = Mtotal * basis
    stress_rows = surf_basis[3:4, :]   # (σ_xz, σ_zz) at the surface, as fn of (a_Pd, a_Sd)
    a1, a2 = stress_rows[1, 1], stress_rows[1, 2]
    b1, b2 = stress_rows[2, 1], stress_rows[2, 2]
    x, y = abs(a1) + abs(a2) >= abs(b1) + abs(b2) ? (a2, -a1) : (b2, -b1)
    nrm = sqrt(abs2(x) + abs2(y))
    x /= nrm
    y /= nrm

    state0 = basis * ComplexF64[x, y]

    state = copy(state0)
    depths_rev = Float64[ztot]
    ux_rev = ComplexF64[state[1]]
    uz_rev = ComplexF64[state[2]]
    sxz_rev = ComplexF64[state[3]]
    szz_rev = ComplexF64[state[4]]
    z = ztot
    idx = 1
    for layer in reverse(finite_layers)
        n = max(4, samples_per_layer)
        dz = layer.thickness / n
        for _ in 1:n
            state = mats[idx] * state
            idx += 1
            z -= dz
            push!(depths_rev, z)
            push!(ux_rev, state[1])
            push!(uz_rev, state[2])
            push!(sxz_rev, state[3])
            push!(szz_rev, state[4])
        end
    end
    depths = reverse(depths_rev)
    ux = reverse(ux_rev)
    uz = reverse(uz_rev)
    sxz = reverse(sxz_rev)
    szz = reverse(szz_rev)

    # Analytic tail into the half-space (for display only): reuses the same
    # propagator with the half-space's own properties, no special-case formula needed.
    tail_depth = max(100.0, 0.5 * ztot)
    tail_samples = 40
    tstate = copy(state0)
    tz = ztot
    for _ in 1:tail_samples
        dz = tail_depth / tail_samples
        Mdown = layer_matrix_PSV(Layer(dz, hs.vp, hs.vs, hs.rho, hs.Qp, hs.Qs), p, ω)
        tstate = Mdown * tstate
        tz += dz
        push!(depths, tz)
        push!(ux, tstate[1])
        push!(uz, tstate[2])
        push!(sxz, tstate[3])
        push!(szz, tstate[4])
    end

    # Rotate the arbitrary overall phase so uz is real at the surface (see docstring).
    phase = abs(uz[1]) > 0 ? conj(uz[1]) / abs(uz[1]) : one(ComplexF64)
    ux = ux .* phase
    uz = uz .* phase
    sxz = sxz .* phase
    szz = szz .* phase

    ur_raw = -imag.(ux)
    uz_raw = real.(uz)
    umax = max(maximum(abs.(ur_raw)), maximum(abs.(uz_raw)))
    umax = umax > 0 ? umax : 1.0
    ur_shape = ur_raw ./ umax
    uz_shape = uz_raw ./ umax

    ρ_at(zval) = begin
        zc = 0.0
        for layer in finite_layers
            zc += layer.thickness
            zval <= zc && return layer.rho
        end
        return hs.rho
    end
    I1 = 0.0
    for i in 2:length(depths)
        dz = depths[i] - depths[i-1]
        ρ_mid = ρ_at((depths[i] + depths[i-1]) / 2)
        I1 += 0.5 * ρ_mid * ((ur_raw[i-1]^2 + uz_raw[i-1]^2) + (ur_raw[i]^2 + uz_raw[i]^2)) * dz / umax^2
    end

    return (depths=depths, ur=ur_shape, uz=uz_shape, I1=I1,
        sxz_resid=abs(sxz[1]) / umax, szz_resid=abs(szz[1]) / umax)
end

# ╔═╡ 4d95603b-28d2-4321-9306-d3cebb5951bf
"""
    generate_rayleigh_synthetic_seismogram(layers, res, distance; source_depth, source_freq, duration, nt=500)

Synthetic Rayleigh-wave (vertical-component) seismogram at horizontal distance `distance`
(km) from a point source at depth `source_depth` (km), using the already-solved
fundamental-mode phase velocities in `res` and the eigenfunctions from
[`compute_rayleigh_eigenfunctions`](@ref) — the direct P-SV analogue of this notebook's
sibling Love-wave notebook's `generate_love_synthetic_seismogram`: the same far-field
normal-mode modal sum, Gaussian source-spectrum weighting, and no separate group-velocity
model layered on top (the dispersed wave packet is a direct consequence of summing many
periods' own phase velocities). Uses the vertical eigenfunction `uz` (matching what a
vertical-component seismometer records) rather than `ur`.
"""
function generate_rayleigh_synthetic_seismogram(layers::Vector{Layer}, res, distance::Float64;
        source_depth::Float64, source_freq::Float64, duration::Float64, nt::Int=500)
    valid = findall(isfinite, res.phase_velocities)
    time_vec = collect(range(0.0, duration, length=nt))
    displacement = zeros(Float64, nt)
    isempty(valid) && return (time=time_vec, displacement=displacement)

    σ_f = 0.08
    for i in valid
        T = res.periods[i]
        c = res.phase_velocities[i]
        f = 1.0 / T
        w = exp(-0.5 * ((f - source_freq) / σ_f)^2)
        w < 1e-3 && continue

        ef = compute_rayleigh_eigenfunctions(layers, T, c)
        j = argmin(abs.(ef.depths .- source_depth))
        u_source = ef.uz[j]
        u_surface = ef.uz[1]
        ω = 2π / T
        k = ω / c
        amp = w * u_source * u_surface / sqrt(ef.I1 * k * distance)
        @. displacement += amp * cos(k * distance - ω * time_vec)
    end

    dmax = maximum(abs.(displacement))
    dmax > 0 && (displacement ./= dmax)
    return (time=time_vec, displacement=displacement)
end

# ╔═╡ 8b06bf24-58b3-450e-90a0-cdcda4bf015c
observed_dispersion = merge_observed_dispersion(
	load_observed_dispersion(observed_phase_file; kind=:phase),
	load_observed_dispersion(observed_group_file; kind=:group),
)

# ╔═╡ 2c205b96-d934-4302-a107-ee6aab03c522
md"## Appendix"

# ╔═╡ 18995223-af32-4f4a-bbaf-a69c4801fb88
md"### The Layered-Medium Widget"

# ╔═╡ 7e919f31-04b3-447c-a5c0-3d03362fb11a
const GUTENBERG_MODEL = [
    Layer(19.0, 6.14, 3.55, 2.74, 1000.0, 1000.0),
    Layer(19.0, 6.58, 3.80, 3.00, 1000.0, 1000.0),
    Layer(12.0, 8.20, 4.65, 3.32, 1000.0, 1000.0),
    Layer(10.0, 8.17, 4.62, 3.34, 1000.0, 1000.0),
    Layer(10.0, 8.14, 4.57, 3.35, 1000.0, 1000.0),
    Layer(10.0, 8.10, 4.51, 3.36, 1000.0, 1000.0),
    Layer(10.0, 8.07, 4.46, 3.37, 1000.0, 1000.0),
    Layer(10.0, 8.02, 4.41, 3.38, 1000.0, 1000.0),
    Layer(25.0, 7.93, 4.37, 3.39, 1000.0, 1000.0),
    Layer(25.0, 7.85, 4.35, 3.41, 1000.0, 1000.0),
    Layer(25.0, 7.89, 4.36, 3.43, 1000.0, 1000.0),
    Layer(25.0, 7.98, 4.38, 3.46, 1000.0, 1000.0),
    Layer(25.0, 8.10, 4.42, 3.48, 1000.0, 1000.0),
    Layer(25.0, 8.21, 4.46, 3.50, 1000.0, 1000.0),
    Layer(50.0, 8.38, 4.54, 3.53, 1000.0, 1000.0),
    Layer(50.0, 8.62, 4.68, 3.58, 1000.0, 1000.0),
    Layer(50.0, 8.87, 4.85, 3.62, 1000.0, 1000.0),
    Layer(50.0, 9.15, 5.04, 3.69, 1000.0, 1000.0),
    Layer(50.0, 9.45, 5.21, 3.82, 1000.0, 1000.0),
    Layer(100.0, 9.88, 5.45, 4.01, 1000.0, 1000.0),
    Layer(100.0, 10.30, 5.76, 4.21, 1000.0, 1000.0),
    Layer(100.0, 10.71, 6.03, 4.40, 1000.0, 1000.0),
    Layer(100.0, 11.10, 6.23, 4.56, 1000.0, 1000.0),
    Layer(100.0, 11.35, 6.32, 4.63, 1000.0, 1000.0)
]

# ╔═╡ 8275ce20-7366-48be-aee1-b6a2c38f1b13
default_layers = GUTENBERG_MODEL;

# ╔═╡ c1c6a307-13b9-4e47-b391-d271820b6c02
begin
	"""
	    LayeredMediumInput(layers=<Gutenberg 5-layer default>; show_vp=true, zmax=300.0)

	A dark-canvas widget for building a layered Earth model by direct manipulation: drag a
	layer boundary up/down to resize it, drag a track's marker left/right to change that
	layer's Vp/Vs/ρ, click empty space to add a layer, or drag a boundary onto its neighbor
	to delete it. The bottom layer is always the half-space (no lower boundary).
	`show_vp=false` hides the Vp track (Love waves don't depend on it).
	"""
	struct LayeredMediumInput
	    boundaries::Vector{Float64}
	    vp::Vector{Float64}
	    vs::Vector{Float64}
	    rho::Vector{Float64}
	    show_vp::Bool
	    zmax::Float64
	end

	"""
	    layered_medium_presets(default_layers) -> Dict{String,Vector{Layer}}

	Build the three quick-load presets shown as buttons on the widget: a 5-layer and a
	2-layer thickness-weighted block average of `default_layers` (reusing the same
	averaging approach the old layer-editor table used for its Gutenberg default), and a
	uniform half-space using the deepest layer's properties.
	"""
	function layered_medium_presets(default_layers::Vector{Layer})
	    function grouped(n::Int)
	        finite_layers = default_layers[1:end-1]
	        edges = round.(Int, range(1, length(finite_layers) + 1, length=n))
	        out = Layer[]
	        for i in 1:(n - 1)
	            block = finite_layers[edges[i]:(edges[i+1]-1)]
	            h = sum(l.thickness for l in block)
	            vp = sum(l.vp * l.thickness for l in block) / h
	            vs = sum(l.vs * l.thickness for l in block) / h
	            rho = sum(l.rho * l.thickness for l in block) / h
	            push!(out, Layer(h, vp, vs, rho, 1000.0, 1000.0))
	        end
	        return vcat(out, default_layers[end])
	    end
	    return Dict(
	        "gutenberg" => grouped(5),
	        "crust2" => grouped(2),
	        "uniform" => [default_layers[end]],
	    )
	end

	function LayeredMediumInput(layers::Vector{Layer}; show_vp::Bool=true, zmax::Float64=300.0)
	    n = length(layers)
	    boundaries = Float64[]
	    d = 0.0
	    for i in 1:(n - 1)
	        d += layers[i].thickness
	        push!(boundaries, d)
	    end
	    LayeredMediumInput(boundaries, [l.vp for l in layers], [l.vs for l in layers], [l.rho for l in layers], show_vp, zmax)
	end
	LayeredMediumInput(; show_vp::Bool=true, zmax::Float64=300.0) =
	    LayeredMediumInput(layered_medium_presets(default_layers)["gutenberg"]; show_vp, zmax)

	Base.get(w::LayeredMediumInput) = Dict{String,Any}(
	    "boundaries" => w.boundaries, "vp" => w.vp, "vs" => w.vs, "rho" => w.rho,
	    "selectedPeriod" => nothing,
	    "distance" => 100.0, "sourceDepth" => 5.0, "sourceFreq" => 0.1, "duration" => 1000.0,
	    "seismogramGenerationId" => 0,
	)

	_lm_preset_js(layers::Vector{Layer}) = let
	    n = length(layers)
	    b = Float64[]; d = 0.0
	    for i in 1:(n-1)
	        d += layers[i].thickness
	        push!(b, d)
	    end
	    "{boundaries:[$(join(b, ","))],vp:[$(join([l.vp for l in layers], ","))],vs:[$(join([l.vs for l in layers], ","))],rho:[$(join([l.rho for l in layers], ","))]}"
	end

	"""
	    lm_dispersion_payload(periods, phase, group, observed::ObservedDispersion)

	Serialize the fundamental-mode dispersion curves and the (possibly empty) observed
	overlay into a JS object literal, for pushing into the `LayeredMediumInput` widget's
	results panel via a `CustomEvent`.
	"""
	function lm_dispersion_payload(periods, phase, group, observed::ObservedDispersion)
	    number(x) = isfinite(x) ? string(round(Float64(x), digits=5)) : "NaN"
	    arr(xs) = "[" * join(number.(xs), ",") * "]"
	    return string(
	        "{\"periods\":", arr(periods),
	        ",\"phase\":", arr(phase),
	        ",\"group\":", arr(group),
	        ",\"obsPeriods\":", arr(observed.periods),
	        ",\"obsPhase\":", arr(observed.phase_velocities),
	        ",\"obsGroup\":", arr(observed.group_velocities),
	        "}",
	    )
	end

	"""
	    lm_eigen_payload(ef) -> String

	Serialize a [`compute_rayleigh_eigenfunctions`](@ref) result into a JS object literal,
	for pushing into the `LayeredMediumInput` widget's eigenfunction tracks via a
	`CustomEvent`.
	"""
	function lm_eigen_payload(ef)
	    number(x) = isfinite(x) ? string(round(Float64(x), digits=5)) : "NaN"
	    arr(xs) = "[" * join(number.(xs), ",") * "]"
	    return string("{\"depths\":", arr(ef.depths), ",\"ur\":", arr(ef.ur), ",\"uz\":", arr(ef.uz), "}")
	end

	"""
	    lm_synth_payload(synth, minArrival, maxArrival) -> String

	Serialize a [`generate_rayleigh_synthetic_seismogram`](@ref) result (plus the
	fastest/slowest group-velocity arrival times, for the dashed reference lines) into a
	JS object literal, for pushing into the `LayeredMediumInput` widget's seismogram panel
	via a `CustomEvent`.
	"""
	function lm_synth_payload(synth, minArrival, maxArrival)
	    number(x) = isfinite(x) ? string(round(Float64(x), digits=5)) : "NaN"
	    arr(xs) = "[" * join(number.(xs), ",") * "]"
	    extra = isfinite(minArrival) && isfinite(maxArrival) ?
	        string(",\"minArrival\":", number(minArrival), ",\"maxArrival\":", number(maxArrival)) : ""
	    return string("{\"time\":", arr(synth.time), ",\"displacement\":", arr(synth.displacement), extra, "}")
	end

	function Base.show(io::IO, ::MIME"text/html", w::LayeredMediumInput)
	    presets = layered_medium_presets(default_layers)
	    presets_js = "{gutenberg:$(_lm_preset_js(presets["gutenberg"])),crust2:$(_lm_preset_js(presets["crust2"])),uniform:$(_lm_preset_js(presets["uniform"]))}"
	    tracks_js = w.show_vp ? "['vp','vs','rho']" : "['vs','rho']"
	    write(io, """
	    <div id="lmwidget">
	    <style>
	    #lmwidget{font-family:sans-serif;color:#d1d5db;position:relative}
	    #lmwidget .lm-titlebar{background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:12px 16px;margin-bottom:14px;text-align:center}
	    #lmwidget .lm-titlebar-headline{font-size:18px;font-weight:700;color:#f3f4f6}
	    #lmwidget .lm-titlebar-sub{font-size:13px;color:#9ca3af;margin-top:4px}
	    #lmwidget .lm-row{display:flex;gap:14px;flex-wrap:wrap;margin-bottom:14px}
	    #lmwidget .lm-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
	    #lmwidget .lm-panel-title{font-size:14px;font-weight:700;color:#f3f4f6;margin-bottom:4px}
	    #lmwidget .lm-caption{font-size:13px;color:#9ca3af;margin-top:4px}
	    #lmwidget canvas{display:block;cursor:crosshair;max-width:100%;height:auto}
	    #lmwidget .lm-mini-group{background:#0b0b0b;border:1px solid #1f2937;border-radius:6px;padding:8px 10px;margin-top:10px}
	    #lmwidget .lm-mini-title{font-size:13px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
	    #lmwidget .lm-actions{display:flex;gap:8px;flex-wrap:wrap}
	    #lmwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px;cursor:pointer}
	    #lmwidget button:hover{background:#707070}
	    #lmwidget .lm-legend{display:flex;gap:12px;flex-wrap:wrap;font-size:12px;color:#9ca3af;margin-top:6px}
	    #lmwidget .lm-swatch{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:middle}
	    #lmwidget .lm-slider-row{display:flex;align-items:center;gap:8px;margin-bottom:8px;font-size:13px}
	    #lmwidget .lm-slider-row label{flex:0 0 150px;color:#d1d5db}
	    #lmwidget .lm-slider-row input[type=range]{flex:1 1 auto}
	    #lmwidget .lm-slider-val{flex:0 0 66px;text-align:right;color:#e5e7eb;font-variant-numeric:tabular-nums}
	    #lmwidget #lm-tooltip{position:absolute;display:none;background:#0b0b0b;border:1px solid #4b5563;border-radius:4px;padding:4px 8px;font-size:12px;color:#e5e7eb;pointer-events:none;z-index:10;white-space:nowrap}
	    </style>

	    <div class="lm-titlebar">
	      <div class="lm-titlebar-headline">Build a layered Earth model by dragging directly on the depth profile.</div>
	      <div class="lm-titlebar-sub">drag a boundary line to resize a layer &middot; drag a track's marker to change Vp/Vs/&rho; &middot; click empty space to add a layer &middot; drag a boundary onto its neighbor to delete it &middot; hover a marker for its value &middot; click the phase-velocity curve to see eigenfunctions at that period</div>
	    </div>

	    <div id="lm-tooltip"></div>

	    <div class="lm-row">
	      <div class="lm-panel" style="flex:3 1 950px">
	        <div class="lm-panel-title">Layered Earth Model &amp; Eigenfunctions</div>
	        <canvas id="lm-editor" width="950" height="440"></canvas>
	        <div class="lm-caption" id="lm-editor-caption">Loading…</div>
	        <div class="lm-legend" id="lm-legend"></div>
	        <div class="lm-mini-group">
	          <div class="lm-mini-title">Presets</div>
	          <div class="lm-actions">
	            <button id="lm-preset-gutenberg">Gutenberg (5 layers)</button>
	            <button id="lm-preset-crust2">2-layer crust</button>
	            <button id="lm-preset-uniform">Uniform half-space</button>
	          </div>
	        </div>
	      </div>
	      <div class="lm-panel" style="flex:2 1 420px">
	        <div class="lm-panel-title">Dispersion Curves</div>
	        <canvas id="lm-disp" width="420" height="440"></canvas>
	        <div class="lm-caption" id="lm-disp-caption">Loading…</div>
	      </div>
	    </div>

	    <div class="lm-row">
	      <div class="lm-panel" style="flex:1 1 100%">
	        <div class="lm-panel-title">Synthetic Seismogram</div>
	        <div class="lm-slider-row">
	          <label for="lm-syn-distance">Distance (km)</label>
	          <input type="range" id="lm-syn-distance" min="10" max="500" step="10" value="100">
	          <span class="lm-slider-val" id="lm-syn-distance-val">100</span>
	        </div>
	        <div class="lm-slider-row">
	          <label for="lm-syn-depth">Source depth (km)</label>
	          <input type="range" id="lm-syn-depth" min="0" max="20" step="1" value="5">
	          <span class="lm-slider-val" id="lm-syn-depth-val">5</span>
	        </div>
	        <div class="lm-slider-row">
	          <label for="lm-syn-freq">Source frequency (Hz)</label>
	          <input type="range" id="lm-syn-freq" min="0.05" max="0.2" step="0.01" value="0.1">
	          <span class="lm-slider-val" id="lm-syn-freq-val">0.10</span>
	        </div>
	        <div class="lm-slider-row">
	          <label for="lm-syn-duration">Duration (s)</label>
	          <input type="range" id="lm-syn-duration" min="500" max="2000" step="100" value="1000">
	          <span class="lm-slider-val" id="lm-syn-duration-val">1000</span>
	        </div>
	        <div class="lm-actions" style="margin-bottom:8px">
	          <button id="lm-syn-generate">Generate Synthetic Seismogram</button>
	        </div>
	        <canvas id="lm-synth" width="1380" height="260"></canvas>
	        <div class="lm-caption" id="lm-synth-caption">Drag the sliders, then click Generate.</div>
	      </div>
	    </div>

	    <script>
	    {
	      if (window.__lmCleanup) { window.__lmCleanup() }
	      const lmInstance = {}
	      window.__lmCurrentInstance = lmInstance
	      const lmController = new AbortController()
	      window.__lmCleanup = () => lmController.abort()
	      const lmSignal = { signal: lmController.signal }

	      const root = document.getElementById('lmwidget')
	      const TRACKS = $tracks_js
	      const RANGES = {vp:[2,13], vs:[1,8], rho:[1,6]}
	      const COLORS = {vp:'#f59e0b', vs:'#3b82f6', rho:'#22c55e'}
	      const LABELS = {vp:'Vp (km/s)', vs:'Vs (km/s)', rho:'ρ (g/cm³)'}
	      const EIGEN_TRACKS = ['ur','uz']
	      const EIGEN_COLORS = {ur:'#a78bfa', uz:'#fb923c'}
	      const EIGEN_LABELS = {ur:'u_r (norm.)', uz:'u_z (norm.)'}
	      const TOTAL_COLS = TRACKS.length + EIGEN_TRACKS.length
	      const PRESETS = $presets_js
	      const MAXLAYERS = 10

	      const state = {
	        boundaries: [$(join(w.boundaries, ","))],
	        vp: [$(join(w.vp, ","))],
	        vs: [$(join(w.vs, ","))],
	        rho: [$(join(w.rho, ","))],
	        selectedPeriod: null,
	        distance: 100.0,
	        sourceDepth: 5.0,
	        sourceFreq: 0.1,
	        duration: 1000.0,
	        seismogramGenerationId: 0,
	      }
	      const zmax = $(w.zmax)

	      function nlayers(){ return state.vp.length }
	      function clamp(v,a,b){ return Math.max(a, Math.min(b, v)) }

	      function publish() {
	        root.value = {
	          boundaries: state.boundaries.slice(), vp: state.vp.slice(), vs: state.vs.slice(), rho: state.rho.slice(),
	          selectedPeriod: state.selectedPeriod,
	          distance: state.distance, sourceDepth: state.sourceDepth, sourceFreq: state.sourceFreq, duration: state.duration,
	          seismogramGenerationId: state.seismogramGenerationId,
	        }
	        root.dispatchEvent(new CustomEvent('input'))
	      }
	      root.value = {
	        boundaries: state.boundaries.slice(), vp: state.vp.slice(), vs: state.vs.slice(), rho: state.rho.slice(),
	        selectedPeriod: state.selectedPeriod,
	        distance: state.distance, sourceDepth: state.sourceDepth, sourceFreq: state.sourceFreq, duration: state.duration,
	        seismogramGenerationId: state.seismogramGenerationId,
	      }

	      // ---------- Layer editor + eigenfunction tracks ----------
	      const editorCanvas = root.querySelector('#lm-editor')
	      const ectx = editorCanvas.getContext('2d')
	      const EM = { l: 46, r: 12, t: 10, b: 26 }
	      const tooltip = root.querySelector('#lm-tooltip')

	      function trackW(){ return (editorCanvas.width - EM.l - EM.r - (TOTAL_COLS-1)*14) / TOTAL_COLS }
	      function trackX0(ti){ return EM.l + ti*(trackW()+14) }
	      function eigenTrackX0(k){ return trackX0(TRACKS.length + k) }
	      function depthMax(){
	        const b = state.boundaries
	        const dataMax = b.length ? b[b.length-1] : 0
	        return Math.max(zmax, dataMax * 1.35 + 20)
	      }
	      function yTop(){ return EM.t }
	      function plotH(){ return editorCanvas.height - EM.t - EM.b }
	      function depthToY(z){ return yTop() + plotH() * (z / depthMax()) }
	      function yToDepth(y){ return clamp((y - yTop()) / plotH(), 0, 1) * depthMax() }
	      function valToX(ti, val){
	        const [lo,hi] = RANGES[TRACKS[ti]]
	        const x0 = trackX0(ti), w = trackW()
	        return x0 + w * clamp((val-lo)/(hi-lo), 0, 1)
	      }
	      function xToVal(ti, x){
	        const [lo,hi] = RANGES[TRACKS[ti]]
	        const x0 = trackX0(ti), w = trackW()
	        return lo + (hi-lo) * clamp((x-x0)/w, 0, 1)
	      }
	      function layerTopBottom(i){
	        const top = i===0 ? 0 : state.boundaries[i-1]
	        const bot = i===nlayers()-1 ? depthMax() : state.boundaries[i]
	        return [top, bot]
	      }

	      let eigenResults = null

	      function drawEigenTracks(){
	        const dmax = depthMax()
	        EIGEN_TRACKS.forEach((tr,k)=>{
	          const x0 = eigenTrackX0(k), w = trackW()
	          ectx.strokeStyle = '#374151'; ectx.strokeRect(x0, yTop(), w, plotH())
	          ectx.fillStyle = '#9ca3af'; ectx.font='12px sans-serif'; ectx.textAlign='center'
	          ectx.fillText(EIGEN_LABELS[tr], x0+w/2, yTop()-2)
	          ectx.font='10px sans-serif'
	          ectx.textAlign='left'; ectx.fillText('-1', x0+2, yTop()+plotH()-4)
	          ectx.textAlign='right'; ectx.fillText('+1', x0+w-2, yTop()+12)
	          const xMid = x0 + w/2
	          ectx.strokeStyle = '#1f2937'; ectx.beginPath(); ectx.moveTo(xMid, yTop()); ectx.lineTo(xMid, yTop()+plotH()); ectx.stroke()

	          if (!eigenResults || !eigenResults.depths || !eigenResults.depths.length) {
	            ectx.fillStyle = '#6b7280'; ectx.font='11px sans-serif'; ectx.textAlign='center'
	            ectx.fillText('click the phase-', x0+w/2, yTop()+plotH()/2-8)
	            ectx.fillText('velocity curve', x0+w/2, yTop()+plotH()/2+4)
	            return
	          }
	          const vals = eigenResults[tr]
	          ectx.beginPath(); ectx.strokeStyle = EIGEN_COLORS[tr]; ectx.lineWidth = 2
	          let started = false
	          eigenResults.depths.forEach((z,i)=>{
	            if (z > dmax) return
	            const vx = x0 + w*0.5*(1 + clamp(vals[i], -1, 1))
	            const vy = depthToY(z)
	            if (!started) { ectx.moveTo(vx,vy); started = true } else ectx.lineTo(vx,vy)
	          })
	          ectx.stroke()
	        })
	      }

	      function drawEditor(){
	        const ctx = ectx, W = editorCanvas.width, H = editorCanvas.height
	        ctx.fillStyle = '#000'; ctx.fillRect(0,0,W,H)
	        const n = nlayers()
	        const dmax = depthMax()
	        ctx.font = '11px sans-serif'
	        const nticks = 6
	        for(let k=0;k<=nticks;k++){
	          const z = dmax*k/nticks
	          const y = depthToY(z)
	          ctx.beginPath(); ctx.moveTo(EM.l-4,y); ctx.lineTo(W-EM.r,y); ctx.strokeStyle='#1f2937'; ctx.lineWidth=1; ctx.stroke()
	          ctx.fillStyle = '#9ca3af'; ctx.textAlign='right'; ctx.fillText(z.toFixed(0), EM.l-6, y+3)
	        }
	        ctx.save(); ctx.translate(12, yTop()+plotH()/2); ctx.rotate(-Math.PI/2)
	        ctx.textAlign='center'; ctx.fillStyle='#9ca3af'; ctx.fillText('Depth (km)', 0, 0); ctx.restore()

	        TRACKS.forEach((tr,ti)=>{
	          const x0 = trackX0(ti), w = trackW()
	          ctx.strokeStyle = '#374151'; ctx.strokeRect(x0, yTop(), w, plotH())
	          ctx.fillStyle = '#9ca3af'; ctx.font='12px sans-serif'; ctx.textAlign='center'
	          ctx.fillText(LABELS[tr], x0+w/2, yTop()-2)
	          const [lo,hi] = RANGES[tr]
	          ctx.font='10px sans-serif'
	          ctx.textAlign='left'; ctx.fillText(lo.toFixed(1), x0+2, yTop()+plotH()-4)
	          ctx.textAlign='right'; ctx.fillText(hi.toFixed(1), x0+w-2, yTop()+12)

	          ctx.beginPath(); ctx.strokeStyle = COLORS[tr]; ctx.lineWidth = 2.5
	          for(let i=0;i<n;i++){
	            const [top,bot] = layerTopBottom(i)
	            const x = valToX(ti, state[tr][i])
	            const yA = depthToY(top), yB = depthToY(i===n-1 ? Math.min(bot, dmax) : bot)
	            if(i===0) ctx.moveTo(x, yA)
	            else ctx.lineTo(x, yA)
	            ctx.lineTo(x, yB)
	            if(i<n-1){
	              const xNext = valToX(ti, state[tr][i+1])
	              ctx.lineTo(xNext, yB)
	            }
	          }
	          ctx.stroke()

	          for(let i=0;i<n;i++){
	            const [top,bot] = layerTopBottom(i)
	            const midY = depthToY((top + Math.min(bot,dmax))/2)
	            const x = valToX(ti, state[tr][i])
	            ctx.beginPath(); ctx.arc(x, midY, 4.5, 0, 2*Math.PI)
	            ctx.fillStyle = COLORS[tr]; ctx.fill(); ctx.strokeStyle='#000'; ctx.lineWidth=1; ctx.stroke()
	          }

	          const [hsTop] = layerTopBottom(n-1)
	          const hsY = depthToY(hsTop)
	          ctx.strokeStyle = COLORS[tr]; ctx.setLineDash([3,3]); ctx.lineWidth=1
	          ctx.beginPath(); ctx.moveTo(x0, hsY); ctx.lineTo(x0+w, hsY); ctx.stroke()
	          ctx.setLineDash([])
	          ctx.fillStyle = '#9ca3af'; ctx.font='11px sans-serif'; ctx.textAlign='center'
	          ctx.fillText('half-space (∞)', x0+w/2, yTop()+plotH()-4)
	        })

	        drawEigenTracks()

	        ctx.strokeStyle = '#f3f4f6'; ctx.lineWidth = 1.5
	        state.boundaries.forEach((b)=>{
	          const y = depthToY(b)
	          ctx.beginPath(); ctx.moveTo(EM.l, y); ctx.lineTo(W-EM.r, y); ctx.stroke()
	        })

	        document.getElementById('lm-editor-caption').textContent =
	          n + ' layer' + (n===1?'':'s') + ' (incl. half-space) · click empty space to add a layer · drag a boundary onto its neighbor to delete'

	        const legend = document.getElementById('lm-legend')
	        legend.innerHTML = TRACKS.map(tr => '<span><span class="lm-swatch" style="background:'+COLORS[tr]+'"></span>'+LABELS[tr]+'</span>').join('')
	          + EIGEN_TRACKS.map(tr => '<span><span class="lm-swatch" style="background:'+EIGEN_COLORS[tr]+'"></span>'+EIGEN_LABELS[tr]+'</span>').join('')
	      }

	      // ---------- hit testing & interaction ----------
	      let drag = null
	      const HIT = 7

	      function canvasXY(ev){
	        const rect = editorCanvas.getBoundingClientRect()
	        const scaleX = editorCanvas.width / rect.width, scaleY = editorCanvas.height / rect.height
	        return [(ev.clientX-rect.left)*scaleX, (ev.clientY-rect.top)*scaleY]
	      }

	      function findBoundaryNear(y){
	        let best=-1, bd=1e9
	        state.boundaries.forEach((b,i)=>{
	          const d = Math.abs(depthToY(b)-y)
	          if(d<HIT && d<bd){ bd=d; best=i }
	        })
	        return best
	      }
	      function findValueNear(x,y){
	        for(let ti=0; ti<TRACKS.length; ti++){
	          const tr = TRACKS[ti]
	          const x0=trackX0(ti), w=trackW()
	          if(x < x0-2 || x > x0+w+2) continue
	          for(let i=0;i<nlayers();i++){
	            const [top,bot]=layerTopBottom(i)
	            const midY = depthToY((top+Math.min(bot,depthMax()))/2)
	            const vx = valToX(ti, state[tr][i])
	            if(Math.abs(x-vx)<HIT && Math.abs(y-midY)<12) return {track:tr, layer:i}
	          }
	        }
	        return null
	      }
	      function trackAt(x){
	        for(let ti=0; ti<TRACKS.length; ti++){
	          const x0=trackX0(ti), w=trackW()
	          if(x>=x0 && x<=x0+w) return ti
	        }
	        return -1
	      }

	      function insertBoundaryAt(depth){
	        if (nlayers() >= MAXLAYERS) return
	        let idx = state.boundaries.findIndex(b => b > depth)
	        if (idx === -1) idx = state.boundaries.length
	        const layerIdx = idx
	        state.boundaries.splice(idx, 0, depth)
	        ;['vp','vs','rho'].forEach(tr => state[tr].splice(layerIdx, 0, state[tr][layerIdx]))
	      }
	      function deleteBoundary(i){
	        if (state.boundaries.length <= 1) return
	        state.boundaries.splice(i,1)
	        ;['vp','vs','rho'].forEach(tr => state[tr].splice(i+1,1))
	      }

	      editorCanvas.addEventListener('mousedown', function(ev){
	        const [x,y] = canvasXY(ev)
	        const bi = findBoundaryNear(y)
	        if (bi >= 0) { drag = {type:'boundary', index:bi}; return }
	        const vh = findValueNear(x,y)
	        if (vh) { drag = {type:'value', track:vh.track, layer:vh.layer}; return }
	        const ti = trackAt(x)
	        if (ti >= 0) {
	          const depth = yToDepth(y)
	          if (depth > 2 && depth < depthMax()-2) {
	            insertBoundaryAt(depth)
	            drawEditor(); publish()
	          }
	        }
	      }, lmSignal)

	      editorCanvas.addEventListener('mouseleave', function(){
	        tooltip.style.display = 'none'
	      }, lmSignal)

	      window.addEventListener('mousemove', function(ev){
	        const [x,y] = canvasXY(ev)
	        if (!drag) {
	          const vh = findValueNear(x,y)
	          if (vh) {
	            const [top,bot] = layerTopBottom(vh.layer)
	            const depthLabel = vh.layer===nlayers()-1 ? top.toFixed(0)+'–∞ km' : top.toFixed(0)+'–'+bot.toFixed(0)+' km'
	            tooltip.textContent = LABELS[vh.track].split(' ')[0] + ': ' + state[vh.track][vh.layer].toFixed(2) + ' · ' + depthLabel
	            const rect = editorCanvas.getBoundingClientRect()
	            const scaleX = rect.width / editorCanvas.width, scaleY = rect.height / editorCanvas.height
	            tooltip.style.left = (rect.left - root.getBoundingClientRect().left + (x+10)*scaleX) + 'px'
	            tooltip.style.top = (rect.top - root.getBoundingClientRect().top + (y-24)*scaleY) + 'px'
	            tooltip.style.display = 'block'
	          } else {
	            tooltip.style.display = 'none'
	          }
	          return
	        }
	        tooltip.style.display = 'none'
	        if (drag.type === 'boundary') {
	          const i = drag.index
	          const lo = i===0 ? 2 : state.boundaries[i-1]+2
	          const hi = i===state.boundaries.length-1 ? depthMax()-2 : state.boundaries[i+1]-2
	          state.boundaries[i] = clamp(yToDepth(y), lo, hi)
	        } else if (drag.type === 'value') {
	          const ti = TRACKS.indexOf(drag.track)
	          let v = xToVal(ti, x)
	          if (drag.track === 'vp') v = Math.max(v, state.vs[drag.layer]*1.3)
	          if (drag.track === 'vs') v = Math.min(v, state.vp[drag.layer]/1.3)
	          state[drag.track][drag.layer] = v
	        }
	        drawEditor()
	      }, lmSignal)

	      window.addEventListener('mouseup', function(){
	        if (!drag) return
	        if (drag.type === 'boundary') {
	          const i = drag.index
	          const nb = state.boundaries.length
	          const neighborAbove = i>0 ? state.boundaries[i-1] : 0
	          const neighborBelow = i<nb-1 ? state.boundaries[i+1] : depthMax()
	          if (state.boundaries[i] - neighborAbove < 4 || neighborBelow - state.boundaries[i] < 4) {
	            deleteBoundary(i)
	          }
	        }
	        drag = null
	        drawEditor(); publish()
	      }, lmSignal)

	      // ---------- Presets ----------
	      function applyPreset(p){
	        state.boundaries = p.boundaries.slice()
	        state.vp = p.vp.slice(); state.vs = p.vs.slice(); state.rho = p.rho.slice()
	        state.selectedPeriod = null
	        eigenResults = null
	        drawEditor(); publish()
	      }
	      root.querySelector('#lm-preset-gutenberg').addEventListener('click', ()=>applyPreset(PRESETS.gutenberg), lmSignal)
	      root.querySelector('#lm-preset-crust2').addEventListener('click', ()=>applyPreset(PRESETS.crust2), lmSignal)
	      root.querySelector('#lm-preset-uniform').addEventListener('click', ()=>applyPreset(PRESETS.uniform), lmSignal)

	      // ---------- Dispersion results panel ----------
	      const dispCanvas = root.querySelector('#lm-disp')
	      const dctx = dispCanvas.getContext('2d')
	      let results = null
	      const DM = { l: 48, r: 12, t: 14, b: 42 }

	      function drawDisp(){
	        const ctx = dctx, W = dispCanvas.width, H = dispCanvas.height
	        ctx.fillStyle = '#000'; ctx.fillRect(0,0,W,H)
	        if (!results || !results.periods || !results.periods.length) {
	          ctx.fillStyle = '#9ca3af'; ctx.font = '13px sans-serif'; ctx.fillText('Computing…', 16, 30)
	          return
	        }
	        const panelH = (H - DM.t - DM.b - 24) / 2
	        const plotW = W - DM.l - DM.r
	        const pmin = Math.min(...results.periods), pmax = Math.max(...results.periods)

	        function panel(y0, ys, obsPeriods, obsYs, color, label, unit){
	          const finite = ys.filter(v=>isFinite(v))
	          let vmin = finite.length? Math.min(...finite): 2, vmax = finite.length? Math.max(...finite): 6
	          const obsFinite = (obsYs||[]).filter(v=>isFinite(v))
	          if (obsFinite.length){ vmin = Math.min(vmin, ...obsFinite); vmax = Math.max(vmax, ...obsFinite) }
	          if (vmax-vmin < 0.5) { vmin -= 0.5; vmax += 0.5 }
	          const pad = (vmax-vmin)*0.1
	          vmin -= pad; vmax += pad
	          function X(p){ return DM.l + plotW * (p-pmin)/((pmax-pmin) || 1) }
	          function Y(v){ return y0 + panelH * (1 - (v-vmin)/((vmax-vmin) || 1)) }
	          ctx.strokeStyle = '#374151'; ctx.strokeRect(DM.l, y0, plotW, panelH)
	          ctx.fillStyle = '#9ca3af'; ctx.font='11px sans-serif'
	          ctx.textAlign='right'; ctx.fillText(vmax.toFixed(1), DM.l-4, y0+10)
	          ctx.fillText(vmin.toFixed(1), DM.l-4, y0+panelH)
	          ctx.save(); ctx.translate(14, y0+panelH/2); ctx.rotate(-Math.PI/2)
	          ctx.textAlign='center'; ctx.fillText(label+' ('+unit+')', 0, 0); ctx.restore()

	          ctx.beginPath(); ctx.strokeStyle = color; ctx.lineWidth = 2.5
	          let started=false
	          results.periods.forEach((p,i)=>{
	            if(!isFinite(ys[i])) { started=false; return }
	            const x=X(p), y=Y(ys[i])
	            if(!started){ ctx.moveTo(x,y); started=true } else ctx.lineTo(x,y)
	          })
	          ctx.stroke()

	          if (obsPeriods && obsPeriods.length){
	            for(let i=0;i<obsPeriods.length;i++){
	              if(!isFinite(obsYs[i])) continue
	              const x=X(obsPeriods[i]), y=Y(obsYs[i])
	              ctx.beginPath(); ctx.moveTo(x-4,y-4); ctx.lineTo(x+4,y+4); ctx.moveTo(x-4,y+4); ctx.lineTo(x+4,y-4)
	              ctx.strokeStyle='#f3f4f6'; ctx.lineWidth=1.5; ctx.stroke()
	            }
	          }

	          if (state.selectedPeriod !== null && state.selectedPeriod >= pmin && state.selectedPeriod <= pmax) {
	            const x = X(state.selectedPeriod)
	            ctx.strokeStyle = '#f3f4f6'; ctx.lineWidth = 1; ctx.setLineDash([2,2])
	            ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y0+panelH); ctx.stroke()
	            ctx.setLineDash([])
	          }
	        }

	        const bottomY0 = DM.t+panelH+24
	        panel(DM.t, results.phase, results.obsPeriods, results.obsPhase, '#3b82f6', 'Phase vel.', 'km/s')
	        panel(bottomY0, results.group, results.obsPeriods, results.obsGroup, '#ef4444', 'Group vel.', 'km/s')

	        const nptick = 5
	        ctx.font='10px sans-serif'; ctx.textAlign='center'
	        for(let k=0;k<=nptick;k++){
	          const p = pmin + (pmax-pmin)*k/nptick
	          const x = DM.l + plotW*(p-pmin)/((pmax-pmin)||1)
	          ctx.beginPath(); ctx.moveTo(x, bottomY0+panelH); ctx.lineTo(x, bottomY0+panelH+4)
	          ctx.strokeStyle='#9ca3af'; ctx.lineWidth=1; ctx.stroke()
	          ctx.fillStyle='#9ca3af'; ctx.fillText(p.toFixed(0), x, bottomY0+panelH+15)
	        }

	        ctx.fillStyle='#9ca3af'; ctx.font='11px sans-serif'; ctx.textAlign='center'
	        ctx.fillText('Period (s)', DM.l+plotW/2, H-6)

	        document.getElementById('lm-disp-caption').textContent =
	          results.periods.length + ' periods · blue = phase velocity · red = group velocity · × = observed · click phase curve for eigenfunctions'
	      }

	      dispCanvas.addEventListener('mousedown', function(ev){
	        if (!results || !results.periods || !results.periods.length) return
	        const rect = dispCanvas.getBoundingClientRect()
	        const scaleX = dispCanvas.width / rect.width
	        const x = (ev.clientX - rect.left) * scaleX
	        const plotW = dispCanvas.width - DM.l - DM.r
	        if (x < DM.l || x > DM.l + plotW) return
	        const pmin = Math.min(...results.periods), pmax = Math.max(...results.periods)
	        const target = pmin + (x - DM.l) / plotW * (pmax - pmin)
	        let best = results.periods[0], bd = Infinity
	        results.periods.forEach(p => { const d = Math.abs(p-target); if (d < bd) { bd = d; best = p } })
	        state.selectedPeriod = best
	        drawDisp()
	        publish()
	      }, lmSignal)

	      window.addEventListener('lm-results', function(ev){
	        if (window.__lmCurrentInstance !== lmInstance) return
	        results = ev.detail
	        drawDisp()
	      }, lmSignal)

	      window.addEventListener('lm-eigen-results', function(ev){
	        if (window.__lmCurrentInstance !== lmInstance) return
	        eigenResults = ev.detail
	        drawEditor()
	      }, lmSignal)

	      // ---------- Synthetic seismogram panel ----------
	      const synthCanvas = root.querySelector('#lm-synth')
	      const sctx = synthCanvas.getContext('2d')
	      let synthResults = null

	      function wireSlider(id, valId, stateKey, fmt){
	        const el = root.querySelector('#'+id), valEl = root.querySelector('#'+valId)
	        el.addEventListener('input', () => {
	          state[stateKey] = parseFloat(el.value)
	          valEl.textContent = fmt(state[stateKey])
	        }, lmSignal)
	      }
	      wireSlider('lm-syn-distance', 'lm-syn-distance-val', 'distance', v => v.toFixed(0))
	      wireSlider('lm-syn-depth', 'lm-syn-depth-val', 'sourceDepth', v => v.toFixed(0))
	      wireSlider('lm-syn-freq', 'lm-syn-freq-val', 'sourceFreq', v => v.toFixed(2))
	      wireSlider('lm-syn-duration', 'lm-syn-duration-val', 'duration', v => v.toFixed(0))

	      root.querySelector('#lm-syn-generate').addEventListener('click', () => {
	        state.seismogramGenerationId += 1
	        publish()
	      }, lmSignal)

	      function drawSynth(){
	        const ctx = sctx, W = synthCanvas.width, H = synthCanvas.height
	        ctx.fillStyle = '#000'; ctx.fillRect(0,0,W,H)
	        if (!synthResults || !synthResults.time || !synthResults.time.length) {
	          ctx.fillStyle = '#9ca3af'; ctx.font = '13px sans-serif'; ctx.fillText('No seismogram generated yet.', 16, 30)
	          return
	        }
	        const M = { l: 50, r: 12, t: 14, b: 34 }
	        const plotW = W - M.l - M.r, plotH = H - M.t - M.b
	        const tmin = synthResults.time[0], tmax = synthResults.time[synthResults.time.length-1]
	        const dmax = Math.max(...synthResults.time.map((_,i)=>Math.abs(synthResults.displacement[i])), 1e-9)
	        function X(t){ return M.l + plotW*(t-tmin)/((tmax-tmin)||1) }
	        function Y(v){ return M.t + plotH*(1-(v/dmax+1)/2) }
	        ctx.strokeStyle = '#374151'; ctx.strokeRect(M.l, M.t, plotW, plotH)
	        ctx.strokeStyle = '#1f2937'; ctx.beginPath(); ctx.moveTo(M.l, Y(0)); ctx.lineTo(M.l+plotW, Y(0)); ctx.stroke()

	        if (synthResults.minArrival !== undefined) {
	          ctx.strokeStyle = '#ef4444'; ctx.setLineDash([4,3]); ctx.lineWidth=1
	          ctx.beginPath(); ctx.moveTo(X(synthResults.minArrival), M.t); ctx.lineTo(X(synthResults.minArrival), M.t+plotH); ctx.stroke()
	          ctx.strokeStyle = '#3b82f6'
	          ctx.beginPath(); ctx.moveTo(X(synthResults.maxArrival), M.t); ctx.lineTo(X(synthResults.maxArrival), M.t+plotH); ctx.stroke()
	          ctx.setLineDash([])
	        }

	        ctx.beginPath(); ctx.strokeStyle = '#e5e7eb'; ctx.lineWidth = 1.3
	        synthResults.time.forEach((t,i)=>{
	          const x=X(t), y=Y(synthResults.displacement[i])
	          if(i===0) ctx.moveTo(x,y); else ctx.lineTo(x,y)
	        })
	        ctx.stroke()

	        ctx.fillStyle='#9ca3af'; ctx.font='10px sans-serif'; ctx.textAlign='center'
	        for(let k=0;k<=5;k++){
	          const t = tmin + (tmax-tmin)*k/5
	          const x = X(t)
	          ctx.beginPath(); ctx.moveTo(x, M.t+plotH); ctx.lineTo(x, M.t+plotH+4); ctx.strokeStyle='#9ca3af'; ctx.stroke()
	          ctx.fillText(t.toFixed(0), x, M.t+plotH+15)
	        }
	        ctx.fillText('Time (s)', M.l+plotW/2, H-4)

	        document.getElementById('lm-synth-caption').textContent =
	          'distance = ' + state.distance.toFixed(0) + ' km · gray = synthetic Rayleigh wave · red/blue dashed = fastest/slowest group-velocity arrival window'
	      }

	      window.addEventListener('lm-synth-results', function(ev){
	        if (window.__lmCurrentInstance !== lmInstance) return
	        synthResults = ev.detail
	        drawSynth()
	      }, lmSignal)

	      drawEditor()
	      drawDisp()
	      drawSynth()
	    }
	    </script>
	    </div>
	    """)
	end

	const _lm_ready = true
end

# ╔═╡ 59e635d1-33aa-4917-9e03-7c579ef8518c
begin
	_lm_ready
	PlutoUI.WideCell(@bind lm LayeredMediumInput(); max_width=1400)
end

# ╔═╡ 863ca41a-959c-49b4-af0d-1876165a64a5
layers = let
	b = Float64.(lm["boundaries"])
	vp = Float64.(lm["vp"])
	vs = Float64.(lm["vs"])
	rho = Float64.(lm["rho"])
	n = length(vp)
	thickness = [i < n ? (i == 1 ? b[1] : b[i] - b[i-1]) : 0.0 for i in 1:n]
	[Layer(thickness[i], vp[i], vs[i], rho[i], 1000.0, 1000.0) for i in 1:n]
end

# ╔═╡ 994df54b-fab7-4a39-868d-aa97f87dec9b
res = compute_rayleigh_dispersion(layers, periods)

# ╔═╡ 0a2f9984-4cf5-4ebe-b12e-974b284256a4
begin
	fundamental_phase = fundamental_phase_velocities(periods, res)
	fundamental_group = group_velocity_from_phase(periods, fundamental_phase)
end

# ╔═╡ dc65b138-737c-44df-a866-74ff223c8489
begin
	_lm_msg = lm_dispersion_payload(periods, fundamental_phase, fundamental_group, observed_dispersion)
	HTML("<script>window.dispatchEvent(new CustomEvent('lm-results', {detail: $_lm_msg}))</script>")
end

# ╔═╡ 89e2a2b5-0d1d-435f-8cc0-8d490bb35e7d
begin
	_lm_eigen_msg = let
		sel = lm["selectedPeriod"]
		if sel === nothing
			"{\"depths\":[],\"ur\":[],\"uz\":[]}"
		else
			idx = argmin(abs.(periods .- sel))
			c = fundamental_phase[idx]
			isfinite(c) ? lm_eigen_payload(compute_rayleigh_eigenfunctions(layers, periods[idx], c)) :
				"{\"depths\":[],\"ur\":[],\"uz\":[]}"
		end
	end
	HTML("<script>window.dispatchEvent(new CustomEvent('lm-eigen-results', {detail: $_lm_eigen_msg}))</script>")
end

# ╔═╡ 49d3d0ac-7ef7-4888-ada1-dac21863ee3a
begin
	lm["seismogramGenerationId"]  # re-run only when Generate is clicked
	_lm_synth_msg = let
		distance = Float64(lm["distance"])
		synth = generate_rayleigh_synthetic_seismogram(layers, (periods=periods, phase_velocities=fundamental_phase), distance;
			source_depth=Float64(lm["sourceDepth"]), source_freq=Float64(lm["sourceFreq"]), duration=Float64(lm["duration"]))
		valid_g = findall(isfinite, fundamental_group)
		if isempty(valid_g)
			lm_synth_payload(synth, NaN, NaN)
		else
			arrivals = distance ./ fundamental_group[valid_g]
			lm_synth_payload(synth, minimum(arrivals), maximum(arrivals))
		end
	end
	HTML("<script>window.dispatchEvent(new CustomEvent('lm-synth-results', {detail: $_lm_synth_msg}))</script>")
end

# ╔═╡ dd68aded-1346-49c0-aa7d-2ed339e221db
md"""
### Reference
This notebook is inspired by the classical implementation of surface wave theory in  

Herrmann, R. B. (2013), *Computer Programs in Seismology*.

The Julia code here provides a **clean reimplementation** of some of the ideas behind CPS, designed for **interactive exploration**.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"

[compat]
PlutoUI = "~0.7.83"
Roots = "~3.0.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "0934c2d072afab596fa9cace800c30c33e633ea5"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "2eeb2c9bef11013efc6f8f97f32ee59b146b09fb"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.44"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CommonSolve]]
git-tree-sha1 = "dd91a10d8b8ae06e15706158eaf1a3e87e97b5f5"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.7"

    [deps.CommonSolve.extensions]
    CommonSolveEnzymeCoreExt = "EnzymeCore"

    [deps.CommonSolve.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "91cfb1cb4f6e27557cc2df798a31eff6089a41eb"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.0"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"
"""

# ╔═╡ Cell order:
# ╠═6e3d880e-82d5-45b0-a823-b2df7d8be6f0
# ╟─c3a8bd5e-41c3-46f4-bfb5-773eadfb67ee
# ╟─59e635d1-33aa-4917-9e03-7c579ef8518c
# ╟─863ca41a-959c-49b4-af0d-1876165a64a5
# ╟─fa019722-6bca-11f1-9c8d-9ffc56627241
# ╠═dc65b138-737c-44df-a866-74ff223c8489
# ╠═89e2a2b5-0d1d-435f-8cc0-8d490bb35e7d
# ╠═49d3d0ac-7ef7-4888-ada1-dac21863ee3a
# ╟─b21471fd-f2f7-4715-b354-48ddd7aa0207
# ╠═5c136747-07b4-454e-9a75-0998a24560df
# ╠═48f94773-cfea-46d5-b341-ae57fbe43e55
# ╠═021afa55-3f7d-4c5e-8768-c91ae697c653
# ╠═c6e30458-2eff-4a66-88de-f75cdbfcd9a4
# ╠═369f3691-78d3-464c-a491-a293f58df913
# ╠═6930149b-0a3a-453f-bf15-8d581f62c433
# ╠═c6187aa6-a890-4b35-8fc8-ef66373a4d2e
# ╠═54effe7a-7d35-40e1-a861-4fa7f6ec9c06
# ╠═87a1ea4f-e41c-44f0-8712-6dad46b1f2a3
# ╟─e55ee346-de7c-4b88-88f2-cd1afbbc73a2
# ╠═3e02c803-9bc3-4123-a4d8-ab07021a0062
# ╠═994df54b-fab7-4a39-868d-aa97f87dec9b
# ╠═7a1c3608-a871-424f-9455-62350343f906
# ╠═6ce1ef3b-2e7f-44db-9b74-1663bb6a84af
# ╠═e54679ad-6c92-43db-b933-8a88bbafdd89
# ╠═195d2bd4-7828-4a98-980c-dd8bc3e4345d
# ╠═984ff052-a12b-4718-9ede-611666ccea76
# ╠═838fcc3e-3c3f-4061-866b-8966e263498c
# ╠═ac3558dc-5882-4898-a8de-543e9483af11
# ╠═4d95603b-28d2-4321-9306-d3cebb5951bf
# ╠═8b06bf24-58b3-450e-90a0-cdcda4bf015c
# ╠═0a2f9984-4cf5-4ebe-b12e-974b284256a4
# ╟─2c205b96-d934-4302-a107-ee6aab03c522
# ╠═9a802663-864a-454c-9e36-e5ee2cf60e87
# ╠═f32a29dc-ed47-4d1e-8ea7-efb287f82c0b
# ╠═3b7fbabb-c4c6-4b3e-9d3a-13f0ac2a7d67
# ╟─18995223-af32-4f4a-bbaf-a69c4801fb88
# ╠═c1c6a307-13b9-4e47-b391-d271820b6c02
# ╠═7e919f31-04b3-447c-a5c0-3d03362fb11a
# ╠═8275ce20-7366-48be-aee1-b6a2c38f1b13
# ╟─dd68aded-1346-49c0-aa7d-2ed339e221db
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
