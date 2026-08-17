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

# ╔═╡ 8d8f440d-4b79-4835-841b-739b5171f979
using PlutoUI, Printf, Roots, LinearAlgebra

# ╔═╡ ee482657-fd34-4c92-95cd-0bf8e254676c
TableOfContents(include_definitions=true)

# ╔═╡ d3bc6646-b22f-4e22-87a8-da37a9aa1d1a
md"""
# Love Wave Dispersion Curves
This notebook provides an environment to explore dispersion of Love waves, a type of surface seismic wave that plays a key role in seismology and Earth structure studies.  

Love waves are horizontally polarized shear waves that are *trapped* in near-surface low-velocity layers. Their dispersion (variation of phase velocity with period) carries information about the elastic properties and layering of the crust and upper mantle.

Here,
- build a layered Earth model interactively, by adjusting thickness, shear velocity (`Vs`), and density (`ρ`) for each layer;
- compute dispersion curves for Love waves and observe how changes in structure affect the shape of the dispersion curve.


##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ f3653043-33f5-4155-aec4-f412a4f305eb
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
	            if kind == :phase
	                push!(period, T); push!(phase, v); push!(group, NaN)
	            elseif kind == :group
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
end

# ╔═╡ 227742dd-3cc0-4615-bde2-09dc7b890095
md"""
### Observed dispersion overlay

Phase file: $(@bind observed_phase_file FilePicker())

Group file: $(@bind observed_group_file FilePicker())
"""

# ╔═╡ 46e9bdbf-8618-4937-b789-c03fbe5396f7
observed_dispersion = merge_observed_dispersion(
	load_observed_dispersion(observed_phase_file; kind=:phase),
	load_observed_dispersion(observed_group_file; kind=:group),
)

# ╔═╡ d2c48242-92ce-11f0-01e6-cb1b8fe1fecb
struct Layer
    thickness::Float64   # km
    vp::Float64          # km/s
    vs::Float64          # km/s
    rho::Float64         # g/cm³
    Qp::Float64          # P-wave quality factor
    Qs::Float64          # S-wave quality factor
end

# ╔═╡ 6c766c4b-8a1c-4f8b-9467-fec5a91fcd05
struct Velocities
    periods::Vector{Float64}
    phase_velocities::Vector{Float64}
    group_velocities::Vector{Float64}
end

# ╔═╡ 7e9acd7c-9cf5-45e0-8892-77c77a1fd983
md"""
Maximum period: $(@bind period_max Slider(20.0:5.0:200.0, default=100.0, show_value=true)) s
"""

# ╔═╡ 756e06e3-4f25-4f55-b242-26a45a1b7dd3
periods = collect(range(5.0, stop=period_max, length=64))

# ╔═╡ 9c9ad5f6-debf-4376-aaf6-759e9b0c949b
"""
    layer_matrix_SH(layer::Layer, ω, c)

Return the 2×2 complex transfer matrix for an SH layer of thickness d:
    [ u_bottom ]   [ A  B ] [ u_top   ]
    [ t_bottom ] = [ C  D ] [ t_top   ]

Where u is horizontal displacement and t is shear traction.
"""
function layer_matrix_SH(layer::Layer, ω::Float64, c::Float64)
    # angular wavenumber horizontally
    k = ω / c

    # shear velocity and density
    vs = layer.vs
    ρ = layer.rho

    # shear modulus (units consistent if rho in mass/volume and vs in velocity units)
    μ = ρ * vs^2

    # vertical wavenumber q (can be complex)
    # q^2 = (ω/vs)^2 - k^2
    q2 = (ω / vs)^2 - k^2
    q = sqrt(complex(q2))          # complex sqrt handles evanescent/propagating

    d = layer.thickness

    # if thickness is effectively zero (half-space), return identity (no propagation)
    if isapprox(d, 0.0; atol=1e-14)
        return Matrix{ComplexF64}(I, 2, 2)
    end

    # arguments for trig/hyperbolic via complex arithmetic
    α = q * d

    # use cos and sin of complex argument (covers cosh/sinh automatically)
    ca = cos(α)
    sa = sin(α)

    # Layer transfer matrix for SH (standard form)
    # [ cos(qd)          sin(qd)/(μ*q)   ]
    # [ -μ*q*sin(qd)     cos(qd)         ]
    M = zeros(ComplexF64, 2, 2)
    M[1, 1] = ca
    M[1, 2] = sa / (μ * q)
    M[2, 1] = -μ * q * sa
    M[2, 2] = ca

    return M
end






# ╔═╡ 4f122bf3-0693-4f13-b716-c0dc2863af25
"""
    compute_love_eigenfunctions(layers, T, c; samples_per_layer=25)

Compute the Love-wave displacement (`u`) and stress (`t`) eigenfunction depth
profiles for a solved period `T` and phase velocity `c`, by shooting **upward**
from the half-space to the free surface.

Shooting upward (rather than down from the surface, as `dispersion_function_SH`
does for the scalar root search) is what keeps this numerically stable: the
state at the base of the model is constructed to satisfy the radiation
condition exactly, so it starts purely on the physically decaying solution
branch with no admixture of the (numerically dominant once amplified) growing
branch — propagating a mixed initial condition downward through hundreds of
kilometers of evanescent layers is the classic Thomson-Haskell instability
that blows up long before reaching the base.

`layers` is expected in the same shape the rest of this notebook already uses
(the finite layers, an artificial thick layer inserted only to give
`dispersion_function_SH` a numerical stand-in for infinity, then the true
half-space) — that artificial layer is dropped here since it isn't needed
once shooting from the true half-space's own analytic decay directly.

Returns a named tuple: `depths`, peak-normalized `u`/`t` shapes, the kinetic
energy integral `I1 = ∫ ρ u(z)² dz` (needed for synthetic-seismogram amplitude
scaling), and `t_surface_residual` (how well the free-surface condition ended
up satisfied — a small value confirms `c` was an accurate root).
"""
function compute_love_eigenfunctions(layers::Vector{Layer}, T::Float64, c::Float64;
    samples_per_layer::Int=25)
    ω = 2π / T
    finite_layers = vcat(layers[1:end-2], layers[end])[1:end-1]
    hs = layers[end]

    q2_hs = (ω / hs.vs)^2 - (ω / c)^2
    q_hs = sqrt(complex(q2_hs))
    decay_rate = imag(q_hs)  # trapped mode: q_hs is purely imaginary, imag > 0
    μ_hs = hs.rho * hs.vs^2
    R = -μ_hs * q_hs

    ztot = sum(l.thickness for l in finite_layers)
    tail_depth = max(100.0, 0.5 * ztot)
    tail_samples = 40

    state = ComplexF64[1.0, R]
    depths_rev = Float64[ztot]
    us_rev = ComplexF64[state[1]]
    ts_rev = ComplexF64[state[2]]
    z = ztot
    for layer in reverse(finite_layers)
        n = max(4, samples_per_layer)
        dz = layer.thickness / n
        sub_up = Layer(-dz, layer.vp, layer.vs, layer.rho, layer.Qp, layer.Qs)
        Mup = layer_matrix_SH(sub_up, ω, c)
        for _ in 1:n
            state = Mup * state
            z -= dz
            push!(depths_rev, z)
            push!(us_rev, state[1])
            push!(ts_rev, state[2])
        end
    end
    depths = reverse(depths_rev)
    us = reverse(us_rev)
    ts = reverse(ts_rev)

    u_bot, t_bot = us[end], ts[end]
    for k in 1:tail_samples
        dz = tail_depth * k / tail_samples
        decay = exp(-decay_rate * dz)
        push!(depths, ztot + dz)
        push!(us, u_bot * decay)
        push!(ts, t_bot * decay)
    end

    u_real = real.(us)
    t_real = real.(ts)
    umax = maximum(abs.(u_real))
    umax = umax > 0 ? umax : 1.0
    u_shape = u_real ./ umax
    t_shape = t_real ./ umax

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
        I1 += 0.5 * ρ_mid * (u_real[i-1]^2 + u_real[i]^2) * dz
    end

    return (depths=depths, u=u_shape, t=t_shape, I1=I1, t_surface_residual=t_real[1] / umax)
end





# ╔═╡ b2add770-6eef-41f6-b4f1-ff93e299d408
"""
    group_velocity_from_phase(periods, phase_velocities) -> Vector{Float64}

Group velocity `U = c / (1 + (T/c)(dc/dT))` from a two-point finite difference of an
*already mode-tracked* phase-velocity curve. Differentiating the raw per-period roots
directly (before tracking) would amplify any mode-jump into a large spurious spike; once
[`track_fundamental_mode`](@ref) has produced one continuous branch, this simple
finite-difference is usually well-behaved.

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

# ╔═╡ c7d305d1-3ca5-4506-8163-cfc8954b0168
"""
    generate_love_synthetic_seismogram(layers, res, distance; source_depth, source_freq, duration, nt=500)

Synthetic Love-wave seismogram at horizontal distance `distance` (km) from a point
source at depth `source_depth` (km), using the already-solved fundamental-mode phase
velocities in `res` and the eigenfunctions from [`compute_love_eigenfunctions`](@ref).

Follows the standard far-field normal-mode surface-wave Green's function (Aki & Richards,
in spirit — not a line-by-line port of CPS's internal FFT-based `hspec96`, which is a much
larger undertaking than this notebook needs): each period contributes
`U(z_source)·U(0)/√(I1·k·x) · cos(kx-ωt)`, weighted by a Gaussian source spectrum centered
on `source_freq`. Summing many periods this way, each with its own (period-dependent)
phase velocity, is what makes the dispersed wave packet emerge on its own — there's no
separate group-velocity model layered on top; it's a direct consequence of the
superposition (stationary-phase interference).
"""
function generate_love_synthetic_seismogram(layers::Vector{Layer}, res, distance::Float64;
    source_depth::Float64, source_freq::Float64, duration::Float64, nt::Int=500)

    valid = findall(isfinite, res.phase_velocities)
    time_vec = collect(range(0.0, duration, length=nt))
    displacement = zeros(Float64, nt)
    isempty(valid) && return (time=time_vec, displacement=displacement)

    # Deliberately wide relative to a single-mode source spectrum: the whole point of
    # this display is to see the packet disperse, which needs enough spectral width for
    # k(ω) = ω/c(ω)'s curvature (real, if modest, in this crustal model) to matter.
    σ_f = 0.08

    for i in valid
        T = res.periods[i]
        c = res.phase_velocities[i]
        f = 1.0 / T
        w = exp(-0.5 * ((f - source_freq) / σ_f)^2)
        w < 1e-3 && continue

        ef = compute_love_eigenfunctions(layers, T, c)
        j = argmin(abs.(ef.depths .- source_depth))
        u_source = ef.u[j]
        u_surface = ef.u[1]

        ω = 2π / T
        k = ω / c
        amp = w * u_source * u_surface / sqrt(ef.I1 * k * distance)

        @. displacement += amp * cos(k * distance - ω * time_vec)
    end

    dmax = maximum(abs.(displacement))
    dmax > 0 && (displacement ./= dmax)
    return (time=time_vec, displacement=displacement)
end

# ╔═╡ e05989a8-f2f1-49ad-a17f-8ddbf5776ab4
"""
    track_fundamental_mode(periods, all_roots) -> Vector{Float64}

Follow the fundamental mode continuously across periods, instead of independently taking
the smallest root at each period. A strong low-velocity zone (or any strong layer
contrast) can pack many closely-spaced modes into the search interval, and picking
"whichever root is smallest" at each period can jump between different physical mode
branches from one period to the next — producing a jagged phase-velocity curve that then
gets amplified into a large spurious spike by [`group_velocity_from_phase`](@ref)'s
finite-difference.

Starting from the longest period (modes are fewest and best separated there — confirmed
empirically: root count grows sharply toward short period as a low-velocity zone traps
more higher modes) and walking toward the shortest, each step picks whichever of the
current period's roots is closest to the previous period's tracked value, so the returned
curve follows one continuous branch instead of jumping between them.
"""
function track_fundamental_mode(periods::AbstractVector{Float64}, all_roots::Vector{Vector{Float64}})
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

# ╔═╡ ac92ec5d-de77-4235-a701-50de13c502b1
"""
    dispersion_function_SH(model, ω, c)

Return complex-valued F(c) whose zero corresponds to a Love-wave mode:
    F(c) = M21 - R * M11
where M = product of layer matrices from top -> bottom and
R = -μ_b * q_b (bottom half-space traction/displacement ratio for decaying solution).
"""
function dispersion_function_SH(layers, ω::Float64, c::Float64)
    N = length(layers)
    if N < 1
        error("Model must contain at least one layer (half-space).")
    end

    # Build total propagator M (top -> bottom)
    M = Matrix{ComplexF64}(I, 2, 2)
    for L in layers
        M = layer_matrix_SH(L, ω, c) * M   # multiply top->down (pre-multiply or post? we choose pre)
    end

    # bottom (half-space) properties (assumed last entry)
    Lb = layers[end]
    vs_b = Lb.vs
    ρ_b = Lb.rho
    μ_b = ρ_b * vs_b^2

    # compute bottom vertical wavenumber q_b (for half-space decaying solution)
    q2_b = (ω / vs_b)^2 - (ω / c)^2
    q_b = sqrt(complex(q2_b))

    # radiation/decay relation tb/ub = R  (we take R = -μ_b * q_b)
    R = -μ_b * q_b

    # M maps [u_top; t_top] -> [u_bot; t_bot]. For t_top=0 (free surface), ub = M11*u0, tb = M21*u0.
    # impose tb/ub = R  => M21/M11 = R  => M21 - R*M11 = 0
    F = M[2, 1] - R * M[1, 1]
    return F
end


# ╔═╡ 03e6b940-ba12-41df-a3a2-2b4e0c44d64b
"""
    solve_phase_velocities_love(layers, T; c_search_pad=0.10) -> Vector{Float64}

Find *every* Love (SH) phase-velocity root at period `T` (seconds), scanning
[`dispersion_function_SH`](@ref)'s real part across the whole search interval (set
relative to the model's own S velocities) rather than stopping at the first sign change.
Finding every root, not just one, is what [`track_fundamental_mode`](@ref) needs to
follow the true fundamental mode continuously across periods instead of picking whichever
root happens to come first.
"""
function solve_phase_velocities_love(layers, T; c_search_pad=0.10)
    ω = 2π / T
    vs_vals = [L.vs for L in layers if L.vs > 0.0]
    isempty(vs_vals) && throw(ArgumentError("Model has no positive shear velocities."))
    vmin = minimum(vs_vals)
    vmax = maximum(vs_vals)
    cmin = max(1e-5, (1 - c_search_pad) * vmin)
    cmax = (1 + c_search_pad) * vmax * 1.5
    g(c) = real(dispersion_function_SH(layers, ω, c))
    try
        return find_zeros(g, cmin, cmax)
    catch
        return Float64[]
    end
end

# ╔═╡ cca474e2-7152-43d5-be66-432847de2897
"""
    solve_love_dispersion(layers, periods) -> Velocities

Compute the Love-wave fundamental-mode dispersion curve: every root at every period via
[`solve_phase_velocities_love`](@ref), tracked into one continuous branch via
[`track_fundamental_mode`](@ref), then differentiated for group velocity via
[`group_velocity_from_phase`](@ref).

A fundamental-mode Love wave only exists when the surface is slower than the half-space
(`layers[1].vs < layers[end].vs`) — that shear-velocity contrast is the waveguide that
traps the SH energy near the surface; without it (e.g. a uniform half-space, or a top
layer no slower than the half-space) there is nothing to trap the energy, and evaluating
[`dispersion_function_SH`](@ref) anyway would just return numerical noise from a huge,
near-degenerate propagator matrix, not a real mode. That case is detected up front and
returns an all-`NaN` curve rather than a spurious near-flat one.
"""
function solve_love_dispersion(layers, periods)
    if layers[1].vs >= layers[end].vs - 1e-6
        nanvec = fill(NaN, length(periods))
        return Velocities(periods, nanvec, nanvec)
    end
    all_roots = [solve_phase_velocities_love(layers, T) for T in periods]
    phase_velocities = track_fundamental_mode(periods, all_roots)
    group_velocities = group_velocity_from_phase(periods, phase_velocities)
    return Velocities(periods, phase_velocities, group_velocities)
end


# ╔═╡ 80086cb5-7e39-444b-a85e-ddf46a67ec42
md"## Appendix"

# ╔═╡ 5a43781a-b3e4-4bac-86d1-bd822c169804
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

# ╔═╡ 05d910d2-527c-416f-9d63-c259fbe8a45d
default_layers = GUTENBERG_MODEL;

# ╔═╡ 99e9fc2c-fe13-4f94-b613-86cedb0f3653
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
	2-layer thickness-weighted block average of `default_layers`, and a uniform half-space
	using the deepest layer's properties.
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

	function Base.show(io::IO, ::MIME"text/html", w::LayeredMediumInput)
	    presets = layered_medium_presets(default_layers)
	    presets_js = "{gutenberg:$(_lm_preset_js(presets["gutenberg"])),crust2:$(_lm_preset_js(presets["crust2"])),uniform:$(_lm_preset_js(presets["uniform"]))}"
	    tracks_js = w.show_vp ? "['vp','vs','rho']" : "['vs','rho']"
	    write(io, """
	    <div id="lmwidget">
	    <style>
	    #lmwidget{font-family:sans-serif;color:#d1d5db}
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
	    </style>

	    <div class="lm-titlebar">
	      <div class="lm-titlebar-headline">Build a layered Earth model by dragging directly on the depth profile.</div>
	      <div class="lm-titlebar-sub">drag a boundary line to resize a layer &middot; drag a track's marker to change Vs/&rho; &middot; click empty space to add a layer &middot; drag a boundary onto its neighbor to delete it</div>
	    </div>

	    <div class="lm-row">
	      <div class="lm-panel" style="flex:3 1 620px">
	        <div class="lm-panel-title">Layered Earth Model</div>
	        <canvas id="lm-editor" width="620" height="440"></canvas>
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
	      const PRESETS = $presets_js
	      const MAXLAYERS = 10

	      const state = {
	        boundaries: [$(join(w.boundaries, ","))],
	        vp: [$(join(w.vp, ","))],
	        vs: [$(join(w.vs, ","))],
	        rho: [$(join(w.rho, ","))],
	      }
	      const zmax = $(w.zmax)

	      function nlayers(){ return state.vp.length }
	      function clamp(v,a,b){ return Math.max(a, Math.min(b, v)) }

	      function publish() {
	        root.value = { boundaries: state.boundaries.slice(), vp: state.vp.slice(), vs: state.vs.slice(), rho: state.rho.slice() }
	        root.dispatchEvent(new CustomEvent('input'))
	      }
	      root.value = { boundaries: state.boundaries.slice(), vp: state.vp.slice(), vs: state.vs.slice(), rho: state.rho.slice() }

	      // ---------- Layer editor ----------
	      const editorCanvas = root.querySelector('#lm-editor')
	      const ectx = editorCanvas.getContext('2d')
	      const EM = { l: 46, r: 12, t: 10, b: 26 }

	      function trackW(){ return (editorCanvas.width - EM.l - EM.r - (TRACKS.length-1)*14) / TRACKS.length }
	      function trackX0(ti){ return EM.l + ti*(trackW()+14) }
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

	        ctx.strokeStyle = '#f3f4f6'; ctx.lineWidth = 1.5
	        state.boundaries.forEach((b)=>{
	          const y = depthToY(b)
	          ctx.beginPath(); ctx.moveTo(EM.l, y); ctx.lineTo(W-EM.r, y); ctx.stroke()
	        })

	        document.getElementById('lm-editor-caption').textContent =
	          n + ' layer' + (n===1?'':'s') + ' (incl. half-space) · click empty space to add a layer · drag a boundary onto its neighbor to delete'

	        const legend = document.getElementById('lm-legend')
	        legend.innerHTML = TRACKS.map(tr => '<span><span class="lm-swatch" style="background:'+COLORS[tr]+'"></span>'+LABELS[tr]+'</span>').join('')
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

	      window.addEventListener('mousemove', function(ev){
	        if (!drag) return
	        const [x,y] = canvasXY(ev)
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
	        drawEditor(); publish()
	      }
	      root.querySelector('#lm-preset-gutenberg').addEventListener('click', ()=>applyPreset(PRESETS.gutenberg), lmSignal)
	      root.querySelector('#lm-preset-crust2').addEventListener('click', ()=>applyPreset(PRESETS.crust2), lmSignal)
	      root.querySelector('#lm-preset-uniform').addEventListener('click', ()=>applyPreset(PRESETS.uniform), lmSignal)

	      // ---------- Dispersion results panel ----------
	      const dispCanvas = root.querySelector('#lm-disp')
	      const dctx = dispCanvas.getContext('2d')
	      let results = null

	      function wrapText(ctx, text, x, y, maxWidth, lineHeight){
	        const words = text.split(' ')
	        let line = ''
	        for (const word of words) {
	          const test = line ? line + ' ' + word : word
	          if (ctx.measureText(test).width > maxWidth && line) {
	            ctx.fillText(line, x, y); line = word; y += lineHeight
	          } else line = test
	        }
	        if (line) ctx.fillText(line, x, y)
	        return y
	      }

	      function drawDisp(){
	        const ctx = dctx, W = dispCanvas.width, H = dispCanvas.height
	        ctx.fillStyle = '#000'; ctx.fillRect(0,0,W,H)
	        if (!results || !results.periods || !results.periods.length) {
	          ctx.fillStyle = '#9ca3af'; ctx.font = '13px sans-serif'; ctx.fillText('Computing…', 16, 30)
	          return
	        }
	        if (results.phase.every(v => !isFinite(v))) {
	          ctx.fillStyle = '#9ca3af'; ctx.font = '13px sans-serif'; ctx.textAlign = 'left'
	          wrapText(ctx,
	            'No fundamental-mode Love wave exists for this model. Love waves need a slower layer above a faster half-space (a shear-velocity waveguide) to trap the energy near the surface — this model has none. Try lowering the top layer Vs below the half-space Vs, or load a preset.',
	            16, 30, W - 32, 18)
	          document.getElementById('lm-disp-caption').textContent = 'no waveguide — no fundamental Love mode for this model'
	          return
	        }
	        const M = { l: 48, r: 12, t: 14, b: 42 }
	        const panelH = (H - M.t - M.b - 24) / 2
	        const plotW = W - M.l - M.r
	        const pmin = Math.min(...results.periods), pmax = Math.max(...results.periods)

	        function panel(y0, ys, obsPeriods, obsYs, color, label, unit){
	          const finite = ys.filter(v=>isFinite(v))
	          let vmin = finite.length? Math.min(...finite): 2, vmax = finite.length? Math.max(...finite): 6
	          const obsFinite = (obsYs||[]).filter(v=>isFinite(v))
	          if (obsFinite.length){ vmin = Math.min(vmin, ...obsFinite); vmax = Math.max(vmax, ...obsFinite) }
	          if (vmax-vmin < 0.5) { vmin -= 0.5; vmax += 0.5 }
	          const pad = (vmax-vmin)*0.1
	          vmin -= pad; vmax += pad
	          function X(p){ return M.l + plotW * (p-pmin)/((pmax-pmin) || 1) }
	          function Y(v){ return y0 + panelH * (1 - (v-vmin)/((vmax-vmin) || 1)) }
	          ctx.strokeStyle = '#374151'; ctx.strokeRect(M.l, y0, plotW, panelH)
	          ctx.fillStyle = '#9ca3af'; ctx.font='11px sans-serif'
	          ctx.textAlign='right'; ctx.fillText(vmax.toFixed(1), M.l-4, y0+10)
	          ctx.fillText(vmin.toFixed(1), M.l-4, y0+panelH)
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
	        }

	        const bottomY0 = M.t+panelH+24
	        panel(M.t, results.phase, results.obsPeriods, results.obsPhase, '#3b82f6', 'Phase vel.', 'km/s')
	        panel(bottomY0, results.group, results.obsPeriods, results.obsGroup, '#ef4444', 'Group vel.', 'km/s')

	        const nptick = 5
	        ctx.font='10px sans-serif'; ctx.textAlign='center'
	        for(let k=0;k<=nptick;k++){
	          const p = pmin + (pmax-pmin)*k/nptick
	          const x = M.l + plotW*(p-pmin)/((pmax-pmin)||1)
	          ctx.beginPath(); ctx.moveTo(x, bottomY0+panelH); ctx.lineTo(x, bottomY0+panelH+4)
	          ctx.strokeStyle='#9ca3af'; ctx.lineWidth=1; ctx.stroke()
	          ctx.fillStyle='#9ca3af'; ctx.fillText(p.toFixed(0), x, bottomY0+panelH+15)
	        }

	        ctx.fillStyle='#9ca3af'; ctx.font='11px sans-serif'; ctx.textAlign='center'
	        ctx.fillText('Period (s)', M.l+plotW/2, H-6)

	        document.getElementById('lm-disp-caption').textContent =
	          results.periods.length + ' periods · blue = phase velocity · red = group velocity · × = observed'
	      }

	      window.addEventListener('lm-results', function(ev){
	        if (window.__lmCurrentInstance !== lmInstance) return
	        results = ev.detail
	        drawDisp()
	      }, lmSignal)

	      drawEditor()
	      drawDisp()
	    }
	    </script>
	    </div>
	    """)
	end

	const _lm_ready = true
end

# ╔═╡ 1af153bc-e471-4608-984b-61f1116dfa16
begin
	_lm_ready
	PlutoUI.WideCell(@bind lm LayeredMediumInput(; show_vp=false); max_width=1400)
end

# ╔═╡ 7bce51e0-9191-42eb-8fe3-cc68f9b2681f
layers = let
	b = Float64.(lm["boundaries"])
	vs = Float64.(lm["vs"])
	rho = Float64.(lm["rho"])
	n = length(vs)
	thickness = [i < n ? (i == 1 ? b[1] : b[i] - b[i-1]) : 0.0 for i in 1:n]
	finite_layers = [Layer(thickness[i], 1.73 * vs[i], vs[i], rho[i], 1000.0, 1000.0) for i in 1:n]
	# Insert an artificial 1000 km layer before the half-space to impose the
	# radiation condition numerically (the half-space itself stays last).
	le = finite_layers[end]
	vcat(finite_layers[1:end-1], Layer(1000.0, le.vp, le.vs, le.rho, 1000.0, 1000.0), le)
end

# ╔═╡ 45579d35-df74-454f-91d4-03542f77e1ac
# ╠═╡ show_logs = false
res = solve_love_dispersion(layers, periods)

# ╔═╡ 6e696382-0a37-46f1-8209-0bf5d7dbc82f
begin
	_lm_msg = lm_dispersion_payload(res.periods, res.phase_velocities, res.group_velocities, observed_dispersion)
	HTML("<script>window.dispatchEvent(new CustomEvent('lm-results', {detail: $_lm_msg}))</script>")
end

# ╔═╡ fc1fd90c-020a-4f2d-aedd-66115ac6a287
md"""
### Reference
This notebook is inspired by the classical implementation of surface wave theory in  

Herrmann, R. B. (2013), *Computer Programs in Seismology*.

The Julia code here provides a **clean reimplementation** of some of the ideas behind CPS, designed for **interactive exploration**.
"""

# ╔═╡ a882824a-a4da-11f0-8360-85b7e5471cb3
md"""
### Eigenfunction Display

Select period for eigenfunction plot (s): $(@bind selected_period Slider(5:5:100, default=20, show_value=true))
"""

# ╔═╡ a8828434-a4da-11f0-ae4f-7707dd8da2ef
let
	period_idx = argmin(abs.(res.periods .- selected_period))
	T_sel = res.periods[period_idx]
	c_sel = res.phase_velocities[period_idx]

	if !isfinite(c_sel)
		md"_No valid mode at this period._"
	else
		ef = compute_love_eigenfunctions(layers, T_sel, c_sel)
		zmax_plot = max(300.0, ef.depths[end])
		boundaries = Float64[]
		zc = 0.0
		for l in layers[1:end-2]
			zc += l.thickness
			push!(boundaries, zc)
		end

		number(x) = isfinite(x) ? string(round(Float64(x), digits=6)) : "0"
		arr(xs) = "[" * join(number.(xs), ",") * "]"

		HTML("""
		<div id="ef-panel" style="font-family:sans-serif;color:#d1d5db">
		<style>
		#ef-panel canvas{display:block;max-width:100%;height:auto;background:#000;border:1px solid #374151;border-radius:6px}
		#ef-panel .ef-caption{font-size:13px;color:#9ca3af;margin-top:4px}
		</style>
		<canvas id="ef-canvas" width="620" height="420"></canvas>
		<div class="ef-caption" id="ef-caption"></div>
		<script>
		{
		  const root = document.getElementById('ef-panel')
		  const depths = $(arr(ef.depths))
		  const u = $(arr(ef.u))
		  const t = $(arr(ef.t))
		  const boundaries = $(arr(boundaries))
		  const zmax = $(zmax_plot)
		  const T_sel = $(T_sel)
		  const c_sel = $(c_sel)
		  const resid = $(ef.t_surface_residual)

		  const canvas = root.querySelector('#ef-canvas')
		  const ctx = canvas.getContext('2d')
		  const W = canvas.width, H = canvas.height
		  const M = { l: 50, r: 20, t: 24, b: 30 }
		  const panelW = (W - M.l - M.r - 30) / 2
		  const plotH = H - M.t - M.b

		  function depthToY(z){ return M.t + plotH * Math.min(z, zmax) / zmax }

		  function drawPanel(x0, values, color, label){
		    const vmax = Math.max(0.05, ...values.map(v=>Math.abs(v)))
		    function X(v){ return x0 + panelW/2 + (panelW/2 - 4) * (v / vmax) }
		    ctx.strokeStyle = '#374151'; ctx.strokeRect(x0, M.t, panelW, plotH)
		    ctx.strokeStyle = '#1f2937'; ctx.beginPath(); ctx.moveTo(X(0), M.t); ctx.lineTo(X(0), M.t+plotH); ctx.stroke()

		    boundaries.forEach(b => {
		      const y = depthToY(b)
		      ctx.strokeStyle = '#4b5563'; ctx.setLineDash([3,3]); ctx.lineWidth=1
		      ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x0+panelW, y); ctx.stroke()
		      ctx.setLineDash([])
		    })

		    ctx.beginPath(); ctx.strokeStyle = color; ctx.lineWidth = 2.2
		    for(let i=0;i<depths.length;i++){
		      const x = X(values[i]), y = depthToY(depths[i])
		      if(i===0) ctx.moveTo(x,y); else ctx.lineTo(x,y)
		    }
		    ctx.stroke()

		    ctx.fillStyle = '#e5e7eb'; ctx.font='13px sans-serif'; ctx.textAlign='center'
		    ctx.fillText(label, x0+panelW/2, M.t-8)
		  }

		  ctx.fillStyle = '#000'; ctx.fillRect(0,0,W,H)
		  drawPanel(M.l, u, '#3b82f6', 'Displacement U(z)')
		  drawPanel(M.l+panelW+30, t, '#ef4444', 'Stress T(z)')

		  ctx.fillStyle = '#9ca3af'; ctx.font='11px sans-serif'; ctx.textAlign='right'
		  const nticks = 6
		  for(let k=0;k<=nticks;k++){
		    const z = zmax*k/nticks
		    const y = depthToY(z)
		    ctx.fillText(z.toFixed(0), M.l-6, y+3)
		  }
		  ctx.save(); ctx.translate(14, M.t+plotH/2); ctx.rotate(-Math.PI/2)
		  ctx.textAlign='center'; ctx.fillText('Depth (km)', 0, 0); ctx.restore()

		  document.getElementById('ef-caption').textContent =
		    'T=' + T_sel.toFixed(1) + ' s, c=' + c_sel.toFixed(4) + ' km/s  ·  dashed lines = layer boundaries  ·  free-surface residual=' + resid.toExponential(1)
		}
		</script>
		</div>
		""")
	end
end

# ╔═╡ 3cce37c8-a4db-11f0-8445-2b16c3968493
md"""
## Synthetic Seismogram Generation

Generate synthetic Love wave seismograms using modal summation with the computed eigenfunctions:
"""

# ╔═╡ 3cce396c-a4db-11f0-bc5a-e107cd4639d6
md"""
**Distance (km):** $(@bind distance Slider(10.0:10.0:500.0, default=100.0, show_value=true))

**Source Depth (km):** $(@bind source_depth Slider(0.0:1.0:20.0, default=5.0, show_value=true))

**Source Frequency (Hz):** $(@bind source_freq Slider(0.05:0.01:0.2, default=0.1, show_value=true))

**Duration (s):** $(@bind duration Slider(500.0:100.0:2000.0, default=1000.0, show_value=true))
"""

# ╔═╡ 3cce3a82-a4db-11f0-9333-3da6f5a8cb24
@bind compute_synthetic Button("Generate Synthetic Seismogram")

# ╔═╡ 3cce3aca-a4db-11f0-aedc-edcc81a83eeb
begin
	compute_synthetic  # re-run when the button is pressed
	synth = generate_love_synthetic_seismogram(layers, res, distance;
		source_depth=source_depth, source_freq=source_freq, duration=duration)
	md"✓ Generated synthetic Love wave seismogram for distance = $(distance) km"
end

# ╔═╡ 3cce3bd8-a4db-11f0-baab-cf3610efa1c8
# # Plot synthetic seismogram
let
	valid_g = findall(isfinite, res.group_velocities)
	arrivals = isempty(valid_g) ? Float64[] : distance ./ res.group_velocities[valid_g]
	min_arrival = isempty(arrivals) ? NaN : minimum(arrivals)
	max_arrival = isempty(arrivals) ? NaN : maximum(arrivals)

	number(x) = isfinite(x) ? string(round(Float64(x), digits=6)) : "0"
	arr(xs) = "[" * join(number.(xs), ",") * "]"

	HTML("""
	<div id="synth-panel" style="font-family:sans-serif;color:#d1d5db">
	<style>
	#synth-panel canvas{display:block;max-width:100%;height:auto;background:#000;border:1px solid #374151;border-radius:6px}
	#synth-panel .synth-caption{font-size:13px;color:#9ca3af;margin-top:4px}
	</style>
	<canvas id="synth-canvas" width="800" height="280"></canvas>
	<div class="synth-caption">black = synthetic Love wave &middot; red/blue dashed = fastest/slowest group-velocity arrival window</div>
	<script>
	{
	  const root = document.getElementById('synth-panel')
	  const t = $(arr(synth.time))
	  const d = $(arr(synth.displacement))
	  const minArrival = $(number(min_arrival))
	  const maxArrival = $(number(max_arrival))
	  const distance = $(distance)

	  const canvas = root.querySelector('#synth-canvas')
	  const ctx = canvas.getContext('2d')
	  const W = canvas.width, H = canvas.height
	  const M = { l: 46, r: 12, t: 24, b: 30 }
	  const plotW = W - M.l - M.r, plotH = H - M.t - M.b
	  const tmax = t.length ? t[t.length-1] : 1

	  function X(tv){ return M.l + plotW * tv / tmax }
	  function Y(v){ return M.t + plotH/2 - v * (plotH/2 - 8) }

	  ctx.fillStyle = '#000'; ctx.fillRect(0,0,W,H)
	  ctx.strokeStyle = '#374151'; ctx.strokeRect(M.l, M.t, plotW, plotH)
	  ctx.strokeStyle = '#1f2937'; ctx.beginPath(); ctx.moveTo(M.l, Y(0)); ctx.lineTo(M.l+plotW, Y(0)); ctx.stroke()

	  if (isFinite(minArrival) && minArrival >= 0 && minArrival <= tmax){
	    ctx.strokeStyle = '#ef4444'; ctx.setLineDash([5,4]); ctx.lineWidth=1.5
	    ctx.beginPath(); ctx.moveTo(X(minArrival), M.t); ctx.lineTo(X(minArrival), M.t+plotH); ctx.stroke()
	    ctx.setLineDash([])
	  }
	  if (isFinite(maxArrival) && maxArrival >= 0 && maxArrival <= tmax){
	    ctx.strokeStyle = '#3b82f6'; ctx.setLineDash([5,4]); ctx.lineWidth=1.5
	    ctx.beginPath(); ctx.moveTo(X(maxArrival), M.t); ctx.lineTo(X(maxArrival), M.t+plotH); ctx.stroke()
	    ctx.setLineDash([])
	  }

	  ctx.beginPath(); ctx.strokeStyle = '#e5e7eb'; ctx.lineWidth = 1.2
	  for(let i=0;i<t.length;i++){
	    const x = X(t[i]), y = Y(d[i])
	    if(i===0) ctx.moveTo(x,y); else ctx.lineTo(x,y)
	  }
	  ctx.stroke()

	  ctx.fillStyle = '#9ca3af'; ctx.font='11px sans-serif'; ctx.textAlign='center'
	  ctx.fillText('Time (s)', M.l+plotW/2, H-6)
	  ctx.fillText('0', M.l, H-6)
	  ctx.fillText(tmax.toFixed(0), M.l+plotW, H-6)
	  ctx.save(); ctx.translate(14, M.t+plotH/2); ctx.rotate(-Math.PI/2)
	  ctx.textAlign='center'; ctx.fillText('Displacement (arbitrary units)', 0, 0); ctx.restore()
	  ctx.textAlign='left'; ctx.fillStyle='#e5e7eb'; ctx.font='13px sans-serif'
	  ctx.fillText('Synthetic Love-wave seismogram — distance = ' + distance.toFixed(0) + ' km', M.l, M.t-8)
	}
	</script>
	</div>
	""")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"

[compat]
PlutoUI = "~0.7.72"
Roots = "~2.2.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "fa2e0e97fe1a9d6484a9c31a05331e65b8bcab0a"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

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
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

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
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

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

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

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

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
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

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "f53232a27a8c1c836d3998ae1e17d898d4df2a46"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

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
git-tree-sha1 = "8a433b1ede5e9be9a7ba5b1cc6698daa8d718f1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.10"

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

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

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

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═ee482657-fd34-4c92-95cd-0bf8e254676c
# ╟─d3bc6646-b22f-4e22-87a8-da37a9aa1d1a
# ╟─1af153bc-e471-4608-984b-61f1116dfa16
# ╟─7bce51e0-9191-42eb-8fe3-cc68f9b2681f
# ╠═45579d35-df74-454f-91d4-03542f77e1ac
# ╠═f3653043-33f5-4155-aec4-f412a4f305eb
# ╟─227742dd-3cc0-4615-bde2-09dc7b890095
# ╠═46e9bdbf-8618-4937-b789-c03fbe5396f7
# ╟─6e696382-0a37-46f1-8209-0bf5d7dbc82f
# ╠═d2c48242-92ce-11f0-01e6-cb1b8fe1fecb
# ╠═6c766c4b-8a1c-4f8b-9467-fec5a91fcd05
# ╟─7e9acd7c-9cf5-45e0-8892-77c77a1fd983
# ╠═756e06e3-4f25-4f55-b242-26a45a1b7dd3
# ╠═9c9ad5f6-debf-4376-aaf6-759e9b0c949b
# ╠═4f122bf3-0693-4f13-b716-c0dc2863af25
# ╠═b2add770-6eef-41f6-b4f1-ff93e299d408
# ╠═c7d305d1-3ca5-4506-8163-cfc8954b0168
# ╠═e05989a8-f2f1-49ad-a17f-8ddbf5776ab4
# ╠═cca474e2-7152-43d5-be66-432847de2897
# ╠═ac92ec5d-de77-4235-a701-50de13c502b1
# ╠═03e6b940-ba12-41df-a3a2-2b4e0c44d64b
# ╟─80086cb5-7e39-444b-a85e-ddf46a67ec42
# ╠═8d8f440d-4b79-4835-841b-739b5171f979
# ╠═99e9fc2c-fe13-4f94-b613-86cedb0f3653
# ╠═5a43781a-b3e4-4bac-86d1-bd822c169804
# ╠═05d910d2-527c-416f-9d63-c259fbe8a45d
# ╟─fc1fd90c-020a-4f2d-aedd-66115ac6a287
# ╠═a882824a-a4da-11f0-8360-85b7e5471cb3
# ╠═a8828434-a4da-11f0-ae4f-7707dd8da2ef
# ╠═3cce37c8-a4db-11f0-8445-2b16c3968493
# ╠═3cce396c-a4db-11f0-bc5a-e107cd4639d6
# ╠═3cce3a82-a4db-11f0-9333-3da6f5a8cb24
# ╠═3cce3aca-a4db-11f0-aedc-edcc81a83eeb
# ╠═3cce3bd8-a4db-11f0-baab-cf3610efa1c8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
