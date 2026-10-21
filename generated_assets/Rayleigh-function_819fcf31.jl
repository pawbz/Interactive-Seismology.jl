### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> tags = ["surfacewaves"]
#> title = "Rayleigh Function"
#> description = "Rayleigh waves on a half-space and on a layer, the H/V ratio, and Scholte and Stoneley interface waves, all from one boundary-matrix recipe."
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

# ╔═╡ 5b807994-416d-11ed-22ff-77b37ff3cfac
begin
	using PlutoUI
	using PlutoTeachingTools
end

# ╔═╡ d921cf1c-ca8b-44d9-8d00-e1ed55267647
ChooseDisplayMode()

# ╔═╡ 915243cb-9619-4eab-8ad4-1d9b821e8d01
PlutoUI.TableOfContents()

# ╔═╡ 8fa36b02-b810-4cbd-ac0a-75e0596d73bb
# App shell: the Expt / Theory / Code / Help frame shared by every notebook, loaded
# from src/assets/app-shell/. Without that folder (a lone downloaded .jl) this cell just
# shows the help text and the notebook stays a plain notebook.
let
	help = md"""
	- **Setting** picks the boundary: a free surface on a half-space, a layer over a half-space, a fluid on a solid (Scholte) or two welded solids (Stoneley). Each setting has presets.
	- The **α, β, ρ** sliders, and **H** for the layer, change the media; solids keep α/β ≥ √2.
	- **f/f₀** (layer only), **z/λ** and **phase** set the frequency, the depth of the orbit and the instant. You can also click or drag in the dispersion and depth panels.
	- **Play** animates the wavefield; switch off the P or SV part to see each half of the wave.
	- The dashed **α +10 %** and **β +10 %** markers and the sensitivity bars show how much the speed follows each parameter.

	**Try this:** load *Soft sediment over bedrock* and drag through f₀ on the dispersion panel: H/V spikes and the orbit flips.
	"""
	dir = joinpath(@__DIR__, "..", "assets", "app-shell")
	if isfile(joinpath(dir, "app-shell.js"))
		HTML(string(
			"<style>", read(joinpath(dir, "app-shell.css"), String), "</style>",
			"<template class=\"app-help\">", repr(MIME"text/html"(), help), "</template>",
			"<script>", read(joinpath(dir, "app-shell.js"), String), "</script>"))
	else
		help
	end
end

# ╔═╡ bb111111-1111-4111-8111-111111111111
begin
	"""
		rayleigh_function(α, β, c)

	The classical Rayleigh function of a homogeneous half-space,
	``\\big(2 - c^2/β^2\\big)^2 - 4\\sqrt{1 - c^2/α^2}\\,\\sqrt{1 - c^2/β^2}``, whose non-zero
	root below `β` is the Rayleigh-wave speed. Used to cross-check the boundary-matrix
	solver in [`Verifying the Solver`](#Verifying-the-Solver).
	"""
	function rayleigh_function(α::Real, β::Real, c::Real)
		q, r = c / β, c / α
		return (2 - q^2)^2 - 4 * sqrt(max(0.0, 1 - q^2)) * sqrt(max(0.0, 1 - r^2))
	end

	"""
		rayleigh_speed(α, β)

	The non-zero root of [`rayleigh_function`](@ref) below `β`, by a bracketing scan and
	bisection (the trivial root ``c = 0`` is excluded).
	"""
	function rayleigh_speed(α::Real, β::Real)
		f(c) = rayleigh_function(α, β, c)
		left, fleft = 0.02β, f(0.02β)
		for right in range(0.03β, 0.999β; length = 400)
			fright = f(right)
			if fleft * fright <= 0
				lo, hi = left, right
				for _ in 1:64
					mid = (lo + hi) / 2
					f(lo) * f(mid) <= 0 ? (hi = mid) : (lo = mid)
				end
				return (lo + hi) / 2
			end
			left, fleft = right, fright
		end
		return NaN
	end
end

# ╔═╡ 66e5c7d4-a32a-4a5a-9053-dab1656cfd60
begin
	"""
		Medium(α, β, ρ)

	Isotropic medium with P speed `α`, S speed `β` (0 for a fluid) and density `ρ`, in km/s,
	km/s and g/cm³.
	"""
	struct Medium
		α::Float64
		β::Float64
		ρ::Float64
	end

	"""
		lame(m)

	Lamé parameters ``λ = ρ(α^2 - 2β^2)`` and ``μ = ρβ^2`` of medium `m`.
	"""
	lame(m::Medium) = (m.ρ * (m.α^2 - 2m.β^2), m.ρ * m.β^2)

	scale_medium(m::Medium; α=1.0, β=1.0, ρ=1.0) = Medium(m.α * α, m.β * β, m.ρ * ρ)
end

# ╔═╡ 5a8f7ef7-64c6-4655-b676-a8624ba28851
begin
	"""
	One surface- or interface-wave problem. The four concrete types differ only in which
	partial waves exist and which boundary conditions they must satisfy; everything else
	(roots, amplitudes, fields) is shared. The interface is at ``z = 0`` (``z = H`` for the
	layer) and ``z`` grows downward.
	"""
	abstract type WaveProblem end

	"""Free surface over the half-space `m`: Rayleigh waves."""
	struct HalfSpace <: WaveProblem
		m::Medium
	end

	"""Layer `m1` of thickness `H` (km) over the half-space `m2`: dispersive Rayleigh waves."""
	struct LayerOverHalfSpace <: WaveProblem
		m1::Medium
		m2::Medium
		H::Float64
	end

	"""Fluid half-space `m1` (β = 0) over solid half-space `m2`: Scholte waves."""
	struct Scholte <: WaveProblem
		m1::Medium
		m2::Medium
	end

	"""Solid half-space `m1` over solid half-space `m2`, welded: Stoneley waves."""
	struct Stoneley <: WaveProblem
		m1::Medium
		m2::Medium
	end

	media(P::HalfSpace) = (P.m,)
	media(P::Union{LayerOverHalfSpace,Scholte,Stoneley}) = (P.m1, P.m2)
	rebuild(P::HalfSpace, ms) = HalfSpace(ms[1])
	rebuild(P::LayerOverHalfSpace, ms) = LayerOverHalfSpace(ms[1], ms[2], P.H)
	rebuild(P::Scholte, ms) = Scholte(ms[1], ms[2])
	rebuild(P::Stoneley, ms) = Stoneley(ms[1], ms[2])
end

# ╔═╡ c167d1e3-6fb0-494f-b261-2dc99b29eb7a
begin
	"""
		p_state(m, p, ω, f, f′)

	Displacement–stress vector ``(i u_x,\\ u_z,\\ i σ_{xz},\\ σ_{zz})`` of a P potential
	``φ = f(z)\\,e^{iω(t - px)}`` in medium `m`, given the depth function's value `f` and
	slope `f′` at the point. All four entries are real.
	"""
	function p_state(m::Medium, p, ω, f, f′)
		λ, μ = lame(m)
		sP = p^2 - inv(m.α^2)
		return (ω * p * f, f′, 2μ * ω * p * f′, ω^2 * (2μ * sP - λ / m.α^2) * f)
	end

	"""
		s_state(m, p, ω, g, g′)

	Same as [`p_state`](@ref) for an SV potential ``ψ = i\\,g(z)\\,e^{iω(t - px)}``.
	"""
	function s_state(m::Medium, p, ω, g, g′)
		_, μ = lame(m)
		return (g′, ω * p * g, μ * ω^2 * (2p^2 - inv(m.β^2)) * g, 2μ * ω * p * g′)
	end

	"""
		decaying(s, ω, z, z0, dir)

	Value and slope of ``e^{-ω η |z - z_0|}`` with ``η = \\sqrt{s}``: an evanescent wave that
	decays away from the interface at `z0`, downward for `dir = 1`, upward for `dir = -1`.
	"""
	function decaying(s, ω, z, z0, dir)
		η = sqrt(max(s, 0.0))
		f = exp(-dir * ω * η * (z - z0))
		return f, -dir * ω * η * f
	end

	"""
		layer_basis(s, ω, z, H)

	Two real depth functions for one wave type inside a layer of thickness `H`, with their
	slopes, both solving ``f'' = ω^2 s f``. For a propagating wave (`s < 0`) they are
	``\\cos qz`` and ``\\sin(qz)/q`` with ``q = ω\\sqrt{-s}``. For an evanescent one (`s > 0`,
	``q = ω\\sqrt{s}``) they are ``e^{-qz}``, which decays from the surface, and
	``e^{-qH}\\sinh(qz)/q``, which grows towards the bottom. Neither grows with ``qH``, so the
	determinant keeps its precision at high frequency, where cosh and sinh would become
	indistinguishable. Both pairs become
	``1`` and ``z`` at `s = 0`, so nothing jumps where the wave switches from one kind to the
	other.
	"""
	function layer_basis(s, ω, z, H)
		if s > 0
			q = ω * sqrt(s)
			e = exp(-q * z)
			if q * H < 20
				g = exp(-q * H)
				return e, -q * e, (z == 0 ? 0.0 : g * sinh(q * z) / q), g * cosh(q * z)
			else
				a, b = exp(-q * (H - z)), exp(-q * (H + z))
				return e, -q * e, (a - b) / (2q), (a + b) / 2
			end
		elseif s < 0
			q = ω * sqrt(-s)
			sn, cn = sincos(q * z)
			return cn, -q * sn, (q * z == 0 ? z : sn / q), cn
		else
			return 1.0, 0.0, z, 1.0
		end
	end

	slow2(p, v) = p^2 - inv(v^2)
end

# ╔═╡ 7342a27d-2b70-43b9-a892-3eeae6db8275
begin
	const Z4 = (0.0, 0.0, 0.0, 0.0)

	"""
		partial_states(P, p, ω, z)

	The displacement–stress vectors, at depth `z`, of every partial wave of problem `P`
	(one per column of its boundary matrix), and a tuple marking each as P (`true`) or SV
	(`false`). Partial waves that don't live at depth `z` contribute zeros.
	"""
	function partial_states(P::HalfSpace, p, ω, z)
		m = P.m
		fP = decaying(slow2(p, m.α), ω, z, 0.0, 1)
		fS = decaying(slow2(p, m.β), ω, z, 0.0, 1)
		return (p_state(m, p, ω, fP...), s_state(m, p, ω, fS...)), (true, false)
	end

	function partial_states(P::LayerOverHalfSpace, p, ω, z; below::Bool = z > P.H)
		m1, m2 = P.m1, P.m2
		kinds = (true, true, false, false, true, false)
		if !below
			cP, cP′, sP, sP′ = layer_basis(slow2(p, m1.α), ω, z, P.H)
			cS, cS′, sS, sS′ = layer_basis(slow2(p, m1.β), ω, z, P.H)
			return (p_state(m1, p, ω, cP, cP′), p_state(m1, p, ω, sP, sP′),
				s_state(m1, p, ω, cS, cS′), s_state(m1, p, ω, sS, sS′), Z4, Z4), kinds
		else
			fP = decaying(slow2(p, m2.α), ω, z, P.H, 1)
			fS = decaying(slow2(p, m2.β), ω, z, P.H, 1)
			return (Z4, Z4, Z4, Z4, p_state(m2, p, ω, fP...), s_state(m2, p, ω, fS...)), kinds
		end
	end

	function partial_states(P::Scholte, p, ω, z; below::Bool = z >= 0)
		kinds = (true, true, false)
		if !below
			fA = decaying(slow2(p, P.m1.α), ω, z, 0.0, -1)
			return (p_state(P.m1, p, ω, fA...), Z4, Z4), kinds
		else
			fP = decaying(slow2(p, P.m2.α), ω, z, 0.0, 1)
			fS = decaying(slow2(p, P.m2.β), ω, z, 0.0, 1)
			return (Z4, p_state(P.m2, p, ω, fP...), s_state(P.m2, p, ω, fS...)), kinds
		end
	end

	function partial_states(P::Stoneley, p, ω, z; below::Bool = z >= 0)
		kinds = (true, false, true, false)
		if !below
			fP = decaying(slow2(p, P.m1.α), ω, z, 0.0, -1)
			fS = decaying(slow2(p, P.m1.β), ω, z, 0.0, -1)
			return (p_state(P.m1, p, ω, fP...), s_state(P.m1, p, ω, fS...), Z4, Z4), kinds
		else
			fP = decaying(slow2(p, P.m2.α), ω, z, 0.0, 1)
			fS = decaying(slow2(p, P.m2.β), ω, z, 0.0, 1)
			return (Z4, Z4, p_state(P.m2, p, ω, fP...), s_state(P.m2, p, ω, fS...)), kinds
		end
	end
end

# ╔═╡ 6dff8f48-78c1-400e-bc51-05c8a8bf4f89
begin
	"""
		boundary_matrix(P, c, ω)

	The real boundary-condition matrix of problem `P` at phase velocity `c`: one row per
	boundary condition, one column per partial-wave amplitude. A surface or interface wave
	exists exactly where its determinant vanishes.
	"""
	function boundary_matrix(P::HalfSpace, c, ω)
		b, _ = partial_states(P, inv(c), ω, 0.0)
		# free surface: σ_xz = σ_zz = 0
		return [b[1][3] b[2][3]; b[1][4] b[2][4]]
	end

	function boundary_matrix(P::LayerOverHalfSpace, c, ω)
		p = inv(c)
		top, _ = partial_states(P, p, ω, 0.0)
		up, _ = partial_states(P, p, ω, P.H; below = false)
		dn, _ = partial_states(P, p, ω, P.H; below = true)
		M = Matrix{Float64}(undef, 6, 6)
		for j in 1:6
			M[1, j] = top[j][3]           # σ_xz = 0 at the free surface
			M[2, j] = top[j][4]           # σ_zz = 0 at the free surface
			for r in 1:4                  # (u_x, u_z, σ_xz, σ_zz) continuous at z = H
				M[2+r, j] = up[j][r] - dn[j][r]
			end
		end
		return M
	end

	function boundary_matrix(P::Scholte, c, ω)
		p = inv(c)
		up, _ = partial_states(P, p, ω, 0.0; below = false)
		dn, _ = partial_states(P, p, ω, 0.0; below = true)
		M = Matrix{Float64}(undef, 3, 3)
		for j in 1:3
			M[1, j] = up[j][2] - dn[j][2]   # u_z continuous
			M[2, j] = up[j][4] - dn[j][4]   # σ_zz continuous
			M[3, j] = dn[j][3]              # σ_xz = 0 on the solid side; u_x may slip
		end
		return M
	end

	function boundary_matrix(P::Stoneley, c, ω)
		p = inv(c)
		up, _ = partial_states(P, p, ω, 0.0; below = false)
		dn, _ = partial_states(P, p, ω, 0.0; below = true)
		M = Matrix{Float64}(undef, 4, 4)
		for j in 1:4, r in 1:4             # welded: all four components continuous
			M[r, j] = up[j][r] - dn[j][r]
		end
		return M
	end
end

# ╔═╡ 3602a131-0ffd-4521-b977-bf6b4c1a6281
begin
	"""
		small_det!(A)

	Determinant of the small square matrix `A` by Gaussian elimination with partial
	pivoting, overwriting `A`.
	"""
	function small_det!(A::AbstractMatrix{Float64})
		n = size(A, 1)
		d = 1.0
		@inbounds for k in 1:n
			piv, ip = abs(A[k, k]), k
			for i in k+1:n
				abs(A[i, k]) > piv && ((piv, ip) = (abs(A[i, k]), i))
			end
			piv == 0 && return 0.0
			if ip != k
				for j in k:n
					A[k, j], A[ip, j] = A[ip, j], A[k, j]
				end
				d = -d
			end
			d *= A[k, k]
			for i in k+1:n
				f = A[i, k] / A[k, k]
				for j in k+1:n
					A[i, j] -= f * A[k, j]
				end
			end
		end
		return d
	end

	"""
		null_vector(M)

	Unit vector `a` with `M*a ≈ 0` for a matrix that is singular to rounding: Gaussian
	elimination with complete pivoting puts the vanishing pivot last; that unknown is set to
	1 and the others follow by back-substitution.
	"""
	function null_vector(M::AbstractMatrix)
		A = Matrix{Float64}(M)
		n = size(A, 1)
		cp = collect(1:n)
		for k in 1:n-1
			_, idx = findmax(abs, view(A, k:n, k:n))
			i, j = idx[1] + k - 1, idx[2] + k - 1
			A[[k, i], :] = A[[i, k], :]
			A[:, [k, j]] = A[:, [j, k]]
			cp[[k, j]] = cp[[j, k]]
			A[k, k] == 0 && continue
			for r in k+1:n
				f = A[r, k] / A[k, k]
				@views A[r, k:n] .-= f .* A[k, k:n]
			end
		end
		x = zeros(n)
		x[n] = 1.0
		for i in n-1:-1:1
			x[i] = -sum(A[i, j] * x[j] for j in i+1:n) / A[i, i]
		end
		a = zeros(n)
		a[cp] = x
		return a ./ sqrt(sum(abs2, a))
	end
end

# ╔═╡ 36599b9f-505b-4a9c-b6e2-464997aad540
begin
	"""
		secular(P, c, ω)

	The determinant of [`boundary_matrix`](@ref): the secular function whose zeros in `c`
	are the surface- or interface-wave speeds.
	"""
	secular(P::WaveProblem, c, ω) = small_det!(boundary_matrix(P, c, ω))

	"""
		normalized_secular(P, c, ω)

	[`secular`](@ref) divided by the product of the matrix's column norms. By Hadamard's
	inequality it lies in [-1, 1], so very different media plot on one axis; dividing by
	positive numbers leaves the roots and the sign pattern unchanged.
	"""
	function normalized_secular(P::WaveProblem, c, ω)
		M = boundary_matrix(P, c, ω)
		nrm = prod(sqrt(sum(abs2, view(M, :, j))) for j in 1:size(M, 2))
		return nrm > 0 ? small_det!(M) / nrm : 0.0
	end

	"""
		speed_window(P)

	The phase-velocity interval searched for roots, ending at the slowest body wave that
	must stay evanescent: the half-space shear speed for Rayleigh waves, the slower of the
	fluid sound speed and the solid shear speed for Scholte waves, the slower shear speed
	for Stoneley waves.
	"""
	speed_window(P::HalfSpace) = (0.5P.m.β, P.m.β)
	speed_window(P::LayerOverHalfSpace) = (0.5min(P.m1.β, P.m2.β), P.m2.β)
	speed_window(P::Scholte) = (0.2min(P.m1.α, P.m2.β), min(P.m1.α, P.m2.β))
	speed_window(P::Stoneley) = (0.2min(P.m1.β, P.m2.β), min(P.m1.β, P.m2.β))

	"""
		find_speeds(P, ω; n=200)

	Every phase velocity in [`speed_window`](@ref) where the secular function changes sign:
	a scan of `n` samples, each bracket refined by bisection. Sorted ascending, so for the
	layer the first entry is the fundamental mode.
	"""
	function find_speeds(P::WaveProblem, ω; n::Int = 200)
		lo, hi = speed_window(P)
		hi *= 1 - 1e-9
		roots = Float64[]
		cprev = lo
		fprev = secular(P, lo, ω)
		for i in 2:n
			ci = lo + (hi - lo) * (i - 1) / (n - 1)
			fi = secular(P, ci, ω)
			if fprev == 0
				push!(roots, cprev)
			elseif fi != 0 && signbit(fi) != signbit(fprev)
				a, b, fa = cprev, ci, fprev
				for _ in 1:55
					m = (a + b) / 2
					fm = secular(P, m, ω)
					signbit(fm) == signbit(fa) ? (a = m; fa = fm) : (b = m)
				end
				push!(roots, (a + b) / 2)
			end
			cprev, fprev = ci, fi
		end
		return roots
	end
end

# ╔═╡ 23cbc525-246b-4808-af37-2d3364857ccc
begin
	"""
		mode_amplitudes(P, c, ω)

	Partial-wave amplitudes of the wave at root `c`: the null vector of the boundary
	matrix.
	"""
	mode_amplitudes(P::WaveProblem, c, ω) = null_vector(boundary_matrix(P, c, ω))

	"""
		field_at(P, a, c, ω, z)

	The wave with amplitudes `a` at depth `z`, as the named tuple
	`(U, W, UP, WP, US, WS, T, N)`: the displacement is ``u_x = U\\sin θ``,
	``u_z = W\\cos θ`` with ``θ = ω(t - x/c)``; `UP`, `WP` and `US`, `WS` are its P and SV
	parts, and `T`, `N` the amplitudes of ``σ_{xz}`` and ``σ_{zz}``.
	"""
	function field_at(P::WaveProblem, a, c, ω, z)
		b, isP = partial_states(P, inv(c), ω, z)
		U = W = UP = WP = US = WS = T = N = 0.0
		for j in eachindex(b)
			bj = b[j]
			U += a[j] * bj[1]; W += a[j] * bj[2]; T += a[j] * bj[3]; N += a[j] * bj[4]
			if isP[j]
				UP += a[j] * bj[1]; WP += a[j] * bj[2]
			else
				US += a[j] * bj[1]; WS += a[j] * bj[2]
			end
		end
		return (; U, W, UP, WP, US, WS, T, N)
	end

	"""
		surface_hv(P, c, ω)

	Horizontal-to-vertical amplitude ratio ``|u_x/u_z|`` at ``z = 0`` (the free surface, or
	the interface for Scholte and Stoneley waves): the ellipticity of the orbit there.
	"""
	function surface_hv(P::WaveProblem, c, ω)
		f = field_at(P, mode_amplitudes(P, c, ω), c, ω, 0.0)
		return abs(f.U / f.W)
	end
end

# ╔═╡ 937eb11f-d10e-4ddb-9511-32e3d7bb69b6
begin
	f0_layer(P::LayerOverHalfSpace) = P.m1.β / (4P.H)

	"""
		layer_frequencies(P; nf=100, ratio=10.0)

	`nf` frequencies spaced logarithmically from ``f_0/\\mathrm{ratio}`` to
	``\\mathrm{ratio}\\,f_0``, with ``f_0 = β_1/4H`` the quarter-wavelength resonance of
	the layer.
	"""
	layer_frequencies(P::LayerOverHalfSpace; nf::Int = 100, ratio = 10.0) =
		f0_layer(P) .* exp10.(range(-log10(ratio), log10(ratio); length = nf))

	"""
		dispersion_curves(P, fs; nmodes=5, hv=true)

	Phase velocities of the first `nmodes` modes of the layer problem `P` at the
	frequencies `fs`, and the fundamental mode's surface H/V when `hv` is set. A mode
	below its cut-off frequency is `NaN`. Returns `(C, HV)` with `C` of size
	`nmodes × length(fs)`.
	"""
	function dispersion_curves(P::LayerOverHalfSpace, fs; nmodes::Int = 5, hv::Bool = true)
		C = fill(NaN, nmodes, length(fs))
		HV = fill(NaN, length(fs))
		for (i, f) in enumerate(fs)
			ω = 2π * f
			r = find_speeds(P, ω)
			for m in 1:min(nmodes, length(r))
				C[m, i] = r[m]
			end
			hv && !isempty(r) && (HV[i] = surface_hv(P, r[1], ω))
		end
		return C, HV
	end

	"""
		group_velocity(fs, c)

	Group velocity ``U = dω/dk`` along one phase-velocity curve `c(f)`, by centred
	differences (one-sided at the ends and beside missing samples).
	"""
	function group_velocity(fs, c)
		ω = 2π .* fs
		k = ω ./ c
		U = fill(NaN, length(fs))
		for i in eachindex(fs)
			isfinite(k[i]) || continue
			lo = i > 1 && isfinite(k[i-1]) ? i - 1 : i
			hi = i < length(fs) && isfinite(k[i+1]) ? i + 1 : i
			hi > lo && (U[i] = (ω[hi] - ω[lo]) / (k[hi] - k[lo]))
		end
		return U
	end
end

# ╔═╡ 45910096-f491-4381-aa54-4e85c1f1345b
begin
	"""
		perturbed(P, i; α=1, β=1, ρ=1)

	Problem `P` with medium `i` (or every medium, for `i = 0`) scaled by the given factors.
	"""
	function perturbed(P::WaveProblem, i::Int; α = 1.0, β = 1.0, ρ = 1.0)
		ms = media(P)
		return rebuild(P, ntuple(k -> (i == 0 || i == k) ? scale_medium(ms[k]; α, β, ρ) : ms[k], length(ms)))
	end

	"""
		nearest_speed(P, ω, c)

	The root of problem `P` closest to `c`, or `NaN` when `P` has none.
	"""
	nearest_speed(P::WaveProblem, ω, c) = (r = find_speeds(P, ω); isempty(r) ? NaN : r[argmin(abs.(r .- c))])

	"""
		log_sensitivities(P, c, ω; h=1e-3)

	``∂\\ln c/∂\\ln m`` of the wave at `c` for each `m ∈ (α, β, ρ)` of each medium, by
	centred differences of the nearest perturbed root. Returns one `(α, β, ρ)` tuple per
	medium; the fluid's β is `NaN`. The velocity entries sum to ``c/U`` (1 for a
	non-dispersive wave) and the density entries to 0.
	"""
	function log_sensitivities(P::WaveProblem, c, ω; h = 1e-3)
		ms = media(P)
		d(i, key) = begin
			key == :β && ms[i].β == 0 && return NaN
			up = perturbed(P, i; (key => 1 + h,)...)
			dn = perturbed(P, i; (key => 1 - h,)...)
			(log(nearest_speed(up, ω, c)) - log(nearest_speed(dn, ω, c))) / (log1p(h) - log1p(-h))
		end
		return ntuple(i -> (d(i, :α), d(i, :β), d(i, :ρ)), length(ms))
	end
end

# ╔═╡ 4a75f3f9-e7c6-4692-97b4-cf09e2d519fa
begin
	"""
		mode_profile(P, c, ω, zs)

	Depth profiles `U(z)`, `W(z)` (``u_x = U\\sin θ``, ``u_z = W\\cos θ``) of the wave at `c`,
	scaled so the larger of the two peaks is 1.
	"""
	function mode_profile(P::WaveProblem, c, ω, zs)
		a = mode_amplitudes(P, c, ω)
		U = Float64[]; W = Float64[]
		for z in zs
			f = field_at(P, a, c, ω, z)
			push!(U, f.U); push!(W, f.W)
		end
		s = max(maximum(abs, U), maximum(abs, W), eps())
		return U ./ s, W ./ s
	end

	"""
		ux_node(zs, U)

	The first depth at or below 0 where the horizontal profile `U` changes sign (linear
	interpolation), or `NaN`: below it the orbit turns the other way.
	"""
	function ux_node(zs, U)
		for i in 2:length(zs)
			zs[i-1] < 0 && continue
			if U[i-1] != 0 && signbit(U[i]) != signbit(U[i-1])
				return zs[i-1] + (zs[i] - zs[i-1]) * U[i-1] / (U[i-1] - U[i])
			end
		end
		return NaN
	end

	"""
		wavefield_frames(P, c, ω, xs, zs; nframes=16)

	Displacement of the wave at `c` on the grid `xs × zs` at `nframes` instants over one
	period, split into its P and SV parts, as `(frame, ix, iz)` arrays sharing one
	normalization (the largest total displacement) so either part can be switched off.
	"""
	function wavefield_frames(P::WaveProblem, c, ω, xs, zs; nframes::Int = 16)
		a = mode_amplitudes(P, c, ω)
		prof = [field_at(P, a, c, ω, z) for z in zs]
		k = ω / c
		sz = (nframes, length(xs), length(zs))
		uxP, uzP, uxS, uzS = zeros(sz), zeros(sz), zeros(sz), zeros(sz)
		for (ix, x) in enumerate(xs), it in 1:nframes
			sθ, cθ = sincos(2π * (it - 1) / nframes - k * x)
			for (iz, f) in enumerate(prof)
				uxP[it, ix, iz] = f.UP * sθ; uzP[it, ix, iz] = f.WP * cθ
				uxS[it, ix, iz] = f.US * sθ; uzS[it, ix, iz] = f.WS * cθ
			end
		end
		nrm = max(maximum(hypot.(uxP .+ uxS, uzP .+ uzS)), eps())
		return (uxP = uxP ./ nrm, uzP = uzP ./ nrm, uxS = uxS ./ nrm, uzS = uzS ./ nrm)
	end
end

# ╔═╡ 826c0d0b-2957-481c-baf4-aec6719e75a6
begin
	"""
		LastCall()

	A one-entry memo: `memo!(f, c, key)` returns `f()`, calling `f` again only when `key`
	differs from the key of the previous call.

	Every interaction with the widget re-runs every cell downstream of the bond, even when
	only the orbit phase moved. The layer's dispersion curves and the sensitivities depend
	on far fewer inputs, so they sit behind one of these each.
	"""
	mutable struct LastCall
		key::Any
		val::Any
	end
	LastCall() = LastCall(nothing, nothing)

	function memo!(f, c::LastCall, key)
		if c.key === nothing || !isequal(c.key, key)
			c.val = f()
			c.key = key
		end
		return c.val
	end

	_rf_layer_memo = LastCall()
	_rf_sens_memo = LastCall()
end

# ╔═╡ 1aac4c99-af17-4b88-97c8-d78fcb1832d5
begin
	"""
		to_json(x)

	A minimal JSON writer for the widget push: numbers (5 significant digits; `NaN` and
	`Inf` become `null`), strings, booleans, vectors, tuples, matrices (as arrays of rows)
	and named tuples.
	"""
	to_json(x::Bool) = x ? "true" : "false"
	to_json(x::Integer) = string(x)
	to_json(x::Real) = isfinite(x) ? string(round(Float64(x); sigdigits = 5)) : "null"
	to_json(::Nothing) = "null"
	to_json(s::AbstractString) = "\"" * replace(s, "\\" => "\\\\", "\"" => "\\\"", "\n" => " ") * "\""
	to_json(s::Symbol) = to_json(string(s))
	to_json(v::Union{AbstractVector,Tuple}) = "[" * join((to_json(x) for x in v), ",") * "]"
	to_json(M::AbstractMatrix) = to_json([M[i, :] for i in 1:size(M, 1)])
	to_json(t::NamedTuple) = "{" * join(("\"$(k)\":" * to_json(v) for (k, v) in pairs(t)), ",") * "}"
end

# ╔═╡ 7beb097d-9d0f-44fa-9161-5c65ad3f676c
begin
	"""
		build_problem(kind, m1, m2, H)

	The [`WaveProblem`](@ref) the widget asks for: `"half"`, `"layer"`, `"scholte"` or
	`"stoneley"`. For Scholte waves medium 1 is the fluid and its shear speed is set to 0.
	"""
	function build_problem(kind::AbstractString, m1::Medium, m2::Medium, H::Float64)
		kind == "layer" && return LayerOverHalfSpace(m1, m2, H)
		kind == "scholte" && return Scholte(Medium(m1.α, 0.0, m1.ρ), m2)
		kind == "stoneley" && return Stoneley(m1, m2)
		return HalfSpace(m2)
	end

	"""
		reference_speed(P)

	The body-wave speed the wave's speed is best compared with, and its label: the
	half-space β for Rayleigh waves, the slower bounding speed for the interface waves.
	"""
	reference_speed(P::HalfSpace) = (P.m.β, "β")
	reference_speed(P::LayerOverHalfSpace) = (P.m2.β, "β₂")
	reference_speed(P::Scholte) = P.m1.α < P.m2.β ? (P.m1.α, "α_f") : (P.m2.β, "β_s")
	reference_speed(P::Stoneley) = P.m1.β < P.m2.β ? (P.m1.β, "β₁") : (P.m2.β, "β₂")

	two_sided(P::WaveProblem) = P isa Union{Scholte,Stoneley}

	"""
		depth_window(P, λ)

	Depth range (km) shown in the eigenfunction panel, and the shallower one used for the
	wavefield: a wavelength or so below the surface, or on both sides of an interface. For
	the layer, [`rayleigh_payload`](@ref) narrows it to where the motion has faded to 3 %.
	"""
	depth_window(P::HalfSpace, λ) = ((0.0, 1.2λ), (0.0, λ))
	depth_window(P::LayerOverHalfSpace, λ) = ((0.0, max(1.2λ, 2P.H)), (0.0, max(λ, 1.5P.H)))
	depth_window(P::Union{Scholte,Stoneley}, λ) = ((-1.2λ, 1.2λ), (-λ, λ))

	"""
		scene(P, L)

	Regions and interface lines the viewer shades and draws, in units of `L` km. Pure
	display description: names, depth bounds and colours.
	"""
	function scene(P::HalfSpace, L)
		(regions = [(z0 = 0.0, z1 = 1e6, name = "half-space", color = "#0d1117")],
		 lines = [(z = 0.0, label = "free surface", color = "#ef4444", dash = false)])
	end
	function scene(P::LayerOverHalfSpace, L)
		(regions = [(z0 = 0.0, z1 = P.H / L, name = "layer", color = "#1d1606"),
		            (z0 = P.H / L, z1 = 1e6, name = "half-space", color = "#0d1117")],
		 lines = [(z = 0.0, label = "free surface", color = "#ef4444", dash = false),
		          (z = P.H / L, label = "z = H", color = "#9ca3af", dash = true)])
	end
	function scene(P::Scholte, L)
		(regions = [(z0 = -1e6, z1 = 0.0, name = "fluid", color = "#0a2036"),
		            (z0 = 0.0, z1 = 1e6, name = "solid", color = "#0d1117")],
		 lines = [(z = 0.0, label = "fluid–solid interface", color = "#94a3b8", dash = false)])
	end
	function scene(P::Stoneley, L)
		(regions = [(z0 = -1e6, z1 = 0.0, name = "upper solid", color = "#151127"),
		            (z0 = 0.0, z1 = 1e6, name = "lower solid", color = "#0d1117")],
		 lines = [(z = 0.0, label = "welded interface", color = "#94a3b8", dash = false)])
	end

	"""
		sensitivity_labels(P)

	Names of the `(α, β, ρ)` entries of [`log_sensitivities`](@ref) for each medium.
	"""
	sensitivity_labels(P::HalfSpace) = (("α", "β", "ρ"),)
	sensitivity_labels(P::LayerOverHalfSpace) = (("α₁", "β₁", "ρ₁"), ("α₂", "β₂", "ρ₂"))
	sensitivity_labels(P::Scholte) = (("α_f", "β_f", "ρ_f"), ("α_s", "β_s", "ρ_s"))
	sensitivity_labels(P::Stoneley) = (("α₁", "β₁", "ρ₁"), ("α₂", "β₂", "ρ₂"))

	no_wave_message(P::HalfSpace) = "No Rayleigh root found."
	no_wave_message(P::LayerOverHalfSpace) = "No trapped mode at this frequency: Rayleigh modes need c < β₂, so the half-space must be faster than the layer."
	no_wave_message(P::Scholte) = "No Scholte root found."
	no_wave_message(P::Stoneley) = "No Stoneley wave for these media: it exists only when the two shear speeds are close (and the densities differ enough). Try the first preset."
end

# ╔═╡ bc2008fa-ad4b-40ff-8bb3-9963d31246ed
"""
	layer_bundle(P)

Everything the dispersion panel shows that depends on the media alone: the frequency
grid, the modes' phase velocities, the fundamental mode's H/V, the same phase velocities
with every α or every β raised by 10 % (the ghost curves), and each mode's group velocity.
"""
function layer_bundle(P::LayerOverHalfSpace)
	fs = layer_frequencies(P)
	C, HV = dispersion_curves(P, fs)
	CA, _ = dispersion_curves(perturbed(P, 0; α = 1.1), fs; hv = false)
	CB, _ = dispersion_curves(perturbed(P, 0; β = 1.1), fs; hv = false)
	Ug = reduce(vcat, (group_velocity(fs, C[m, :])' for m in 1:size(C, 1)))
	return (; fs, C, HV, CA, CB, Ug)
end

# ╔═╡ 00b0672a-ff69-46aa-aa07-2a69052a14c6
begin
	frequency_of(P::LayerOverHalfSpace, logf) = f0_layer(P) * exp10(logf)
	frequency_of(P::WaveProblem, logf) = 1.0

	"""
		rayleigh_payload(P, logf, depth, phase, mode)

	Everything the widget draws for problem `P`, as one named tuple: the secular function
	(or, for the layer, the dispersion and H/V curves at ``f = f_0 10^{\\mathrm{logf}}``),
	the depth profiles, the particle orbit at depth `depth` wavelengths, the wavefield
	frames, the sensitivities and the ±10 % ghost speeds. Lengths are in wavelengths for the
	non-dispersive cases and in km for the layer. `mode` (1 = fundamental) picks the
	layer's mode.
	"""
	function rayleigh_payload(P::WaveProblem, logf, depth, phase, mode)
		isL = P isa LayerOverHalfSpace
		f = frequency_of(P, logf)
		ω = 2π * f
		roots = find_speeds(P, ω)
		cref, crefname = reference_speed(P)
		lo, hi = speed_window(P)

		curves = if isL
			B = memo!(() -> layer_bundle(P), _rf_layer_memo, P)
			(disp = (f = B.fs, C = B.C, CA = B.CA, CB = B.CB, Ug = B.Ug, HV = B.HV, f0 = f0_layer(P),
				lines = [(v = P.m1.β, label = "β₁"), (v = P.m2.β, label = "β₂")]),)
		else
			cs = collect(range(lo, hi * (1 - 1e-6); length = 240))
			lines = P isa HalfSpace ? [(v = P.m.β, label = "β")] :
				P isa Scholte ? [(v = P.m1.α, label = "α_f"), (v = P.m2.β, label = "β_s")] :
				[(v = P.m1.β, label = "β₁"), (v = P.m2.β, label = "β₂")]
			(sec = (c = cs, d = [normalized_secular(P, c, ω) for c in cs], roots = roots, lines = lines),)
		end
		base = (kind = P isa HalfSpace ? "half" : isL ? "layer" : P isa Scholte ? "scholte" : "stoneley",
			f = f, f0 = isL ? f0_layer(P) : NaN, cref = cref, crefname = crefname, nroots = length(roots))

		if isempty(roots)
			return merge(base, curves, (ok = false, msg = no_wave_message(P)))
		end

		mode = isL ? clamp(mode, 1, length(roots)) : 1
		c = roots[mode]
		λ = c / f
		L = isL ? 1.0 : λ                           # plotting length unit, km
		a = mode_amplitudes(P, c, ω)

		# depth profiles and the u_x node
		(zp0, zp1), (zw0, zw1) = depth_window(P, λ)
		if isL
			# a soft layer can hold nearly all the motion: show down to where it has faded
			zprobe = range(0, max(1.5λ, 3P.H); length = 400)
			amp = [(q = field_at(P, a, c, ω, z); max(abs(q.U), abs(q.W))) for z in zprobe]
			zlast = zprobe[findlast(>=(0.03 * maximum(amp)), amp)]
			zp1, zw1 = max(2P.H, 1.2zlast), max(1.5P.H, zlast)
		end
		zs = collect(range(zp0, zp1; length = 181))
		prof = [field_at(P, a, c, ω, z) for z in zs]
		U = [q.U for q in prof]; W = [q.W for q in prof]
		pscale = max(maximum(abs, U), maximum(abs, W), eps())
		node = ux_node(zs, U)

		# particle orbit at the chosen depth
		zsel = depth * λ
		q = field_at(P, a, c, ω, zsel)
		θs = range(0, 2π; length = 181)
		oscale = max(abs(q.U), abs(q.W), eps())
		θ = 2π * phase
		orbit = (x = [q.U * sin(t) / oscale for t in θs], z = [q.W * cos(t) / oscale for t in θs],
			sx = q.U * sin(θ) / oscale, sz = q.W * cos(θ) / oscale,
			sense = q.U * q.W > 0 ? "retrograde" : "prograde",
			amp = max(abs(q.U), abs(q.W)) / pscale)

		# wavefield frames
		xs = collect(range(0, 2λ; length = 28))
		zw = collect(range(zw0, zw1; length = two_sided(P) ? 21 : 15))
		wf = wavefield_frames(P, c, ω, xs, zw)
		# iz fastest, then ix, then frame; 3 digits is far finer than an arrow on screen
		flat(A) = round.(vec(permutedims(A, (3, 2, 1))); sigdigits = 3)

		# sensitivities, ghosts and c/U
		sens = memo!(() -> log_sensitivities(P, c, ω), _rf_sens_memo, (P, f, mode))
		labels = sensitivity_labels(P)
		srows = [(label = labels[i][j], kind = ("α", "β", "ρ")[j], value = sens[i][j])
			for i in eachindex(sens) for j in 1:3 if isfinite(sens[i][j])]
		cA = nearest_speed(perturbed(P, 0; α = 1.1), ω, c)
		cB = nearest_speed(perturbed(P, 0; β = 1.1), ω, c)
		cU = if isL
			δ = 1e-3
			cm, cp = nearest_speed(P, ω * (1 - δ), c), nearest_speed(P, ω * (1 + δ), c)
			Ugrp = 2δ * ω / (ω * (1 + δ) / cp - ω * (1 - δ) / cm)
			c / Ugrp
		else
			1.0
		end
		sc = scene(P, L)
		s0 = field_at(P, a, c, ω, 0.0)              # surface, or solid side of the interface

		return merge(base, curves, (
			ok = true, msg = "", unit = isL ? "km" : "λ", lam = λ / L, lamkm = λ, mode = mode,
			roots = roots, c = c, hv = abs(s0.U / s0.W), vh = abs(s0.W / s0.U),
			node = node / λ, depth = depth, zsel = zsel / L,
			prof = (z = zs ./ L, U = U ./ pscale, W = W ./ pscale, node = node / L),
			orbit = orbit,
			wf = (xs = xs ./ L, zs = zw ./ L, nx = length(xs), nz = length(zw), nframes = 16,
				uxP = flat(wf.uxP), uzP = flat(wf.uzP), uxS = flat(wf.uxS), uzS = flat(wf.uzS)),
			regions = sc.regions, lines = sc.lines,
			sens = srows, sumv = sum(r.value for r in srows if r.kind != "ρ"), cU = cU,
			ghost = (cA = cA, cB = cB, dA = cA / c - 1, dB = cB / c - 1),
		))
	end

end

# ╔═╡ a48f6938-2667-4401-9721-ee53a692c6ef
begin
	"""
	Starting points for each setting, `(name, a1, b1, r1, a2, b2, r2, H)`: medium 1 is the
	layer, the fluid or the upper solid, medium 2 the half-space, the solid or the lower
	solid (km/s, g/cm³, km). The first entry of each list is what selecting the setting
	loads.
	"""
	const RF_PRESETS = (
		half = [
			(name = "Poisson solid (α 6.0, β 3.46)", a1 = 6.0, b1 = 3.46, r1 = 2.7, a2 = 6.0, b2 = 3.464, r2 = 2.7, H = 35.0),   # β = α/√3
			(name = "Mantle-like (α 8.0, β 4.49)", a1 = 6.0, b1 = 3.46, r1 = 2.7, a2 = 8.0, b2 = 4.49, r2 = 3.3, H = 35.0),
			(name = "ν = 0 (α/β = √2)", a1 = 6.0, b1 = 3.46, r1 = 2.7, a2 = 6.0, b2 = 4.243, r2 = 2.7, H = 35.0),
			(name = "Nearly incompressible (ν = 0.45)", a1 = 6.0, b1 = 3.46, r1 = 2.7, a2 = 6.0, b2 = 1.809, r2 = 2.0, H = 35.0),
		],
		layer = [
			(name = "Crust over mantle (H 35 km)", a1 = 6.0, b1 = 3.46, r1 = 2.7, a2 = 8.0, b2 = 4.49, r2 = 3.3, H = 35.0),
			(name = "Soft sediment over bedrock (H 30 m)", a1 = 1.6, b1 = 0.2, r1 = 1.8, a2 = 4.0, b2 = 2.0, r2 = 2.5, H = 0.03),
			(name = "Weak contrast (H 1 km)", a1 = 5.0, b1 = 2.9, r1 = 2.6, a2 = 6.0, b2 = 3.46, r2 = 2.7, H = 1.0),
		],
		scholte = [
			(name = "Soft seafloor sediment", a1 = 1.5, b1 = 0.0, r1 = 1.0, a2 = 1.8, b2 = 0.4, r2 = 1.8, H = 35.0),
			(name = "Hard-rock seafloor", a1 = 1.5, b1 = 0.0, r1 = 1.0, a2 = 6.0, b2 = 3.46, r2 = 2.7, H = 35.0),
		],
		stoneley = [
			(name = "Equal β, ρ₂ = 4ρ₁", a1 = 6.0, b1 = 3.46, r1 = 2.0, a2 = 6.0, b2 = 3.46, r2 = 8.0, H = 35.0),
			(name = "Crust against mantle (no Stoneley wave)", a1 = 6.0, b1 = 3.46, r1 = 2.7, a2 = 8.0, b2 = 4.49, r2 = 3.3, H = 35.0),
		],
	)

	"""
		rf_preset_state(kind, i)

	The widget state (the bond's value) for preset `i` of setting `kind`, at the
	quarter-wavelength frequency, at the surface, a quarter period into the orbit, on the
	fundamental mode.
	"""
	function rf_preset_state(kind::Symbol, i::Int = 1)
		p = RF_PRESETS[kind][i]
		return Dict{String,Any}("kind" => string(kind), "a1" => p.a1, "b1" => p.b1, "r1" => p.r1,
			"a2" => p.a2, "b2" => p.b2, "r2" => p.r2, "H" => p.H,
			"logf" => 0.0, "depth" => 0.0, "phase" => 0.25, "mode" => 1)
	end
end

# ╔═╡ bb222222-2222-4222-8222-222222222222
begin
	"""
		RayleighFunctionInput(state)

	The experiment: a setting switch (half-space, layer, Scholte, Stoneley), presets, the
	two media, frequency, depth and orbit phase, plus the panels that draw what
	[`rayleigh_payload`](@ref) computes. Its bond value is the state `Dict`.
	"""
	struct RayleighFunctionInput
		state::Dict{String,Any}
	end
	RayleighFunctionInput() = RayleighFunctionInput(rf_preset_state(:half, 1))

	Base.get(w::RayleighFunctionInput) = copy(w.state)

	function Base.show(io::IO, ::MIME"text/html", w::RayleighFunctionInput)
		presets = to_json((half = RF_PRESETS.half, layer = RF_PRESETS.layer,
			scholte = RF_PRESETS.scholte, stoneley = RF_PRESETS.stoneley))
		state = "{" * join(("\"$(k)\":" * to_json(v) for (k, v) in w.state), ",") * "}"
		write(io, """
<div id="rf-widget" class="app-widget">
<style>
#rf-widget{width:100%;max-width:1400px;margin:auto;color:#e5e7eb;font:14px system-ui,sans-serif}
#rf-widget *{box-sizing:border-box}
#rf-widget .rf-title{padding:10px 14px;margin-bottom:10px;text-align:center;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px}
#rf-widget .rf-title b{font-size:17px}
#rf-widget .rf-hint{margin-top:3px;color:#9ca3af;font-size:13px}
#rf-widget .rf-workspace{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px}
#rf-widget .rf-panel,#rf-widget .rf-group{min-width:0;padding:9px;background:#050505;border:1px solid #2f3744;border-radius:6px}
#rf-widget .rf-panel{display:flex;flex-direction:column}
#rf-widget .rf-wide{grid-column:1/-1}
#rf-widget .rf-panel-title{font-size:14px;font-weight:700;color:#f3f4f6}
#rf-widget .rf-caption{margin:2px 0 6px;color:#9ca3af;font-size:13px;line-height:1.3;min-height:17px}
#rf-widget .rf-cv{position:relative;height:250px;min-height:0;overflow:hidden}
#rf-widget .rf-wide .rf-cv{height:290px}
#rf-widget canvas{position:absolute;left:0;top:0;display:block;background:#000;border-radius:4px;touch-action:none}
#rf-widget .rf-legend{display:flex;gap:12px;flex-wrap:wrap;align-items:center;margin-top:6px;font-size:13px;color:#d1d5db}
#rf-widget .rf-key{display:inline-flex;align-items:center;gap:5px;cursor:pointer;user-select:none}
#rf-widget .rf-sw{width:18px;height:3px;border-radius:2px;background:currentColor;display:inline-block}
#rf-widget .rf-dash{background:repeating-linear-gradient(90deg,currentColor 0 5px,transparent 5px 9px)}
#rf-widget .rf-controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:8px;margin-top:10px}
#rf-widget .rf-group-title{margin-bottom:6px;font-size:16px;font-weight:700;color:#f3f4f6}
#rf-widget .rf-row{display:grid;grid-template-columns:62px minmax(60px,1fr) 78px;gap:6px;align-items:center;margin:6px 0;color:#d1d5db}
#rf-widget input[type=range]{width:100%;min-width:0;accent-color:#f59e0b}
#rf-widget .rf-value{color:#fbbf24;text-align:right;font-variant-numeric:tabular-nums;font-size:13px}
#rf-widget select{width:100%;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:5px;font-size:13px}
#rf-widget .rf-seg{display:grid;grid-template-columns:1fr 1fr;gap:4px;margin-bottom:7px}
#rf-widget button.rf-btn{border:1px solid #9ca3af;border-radius:4px;padding:6px 10px;background:#606060;color:#f3f4f6;font-size:14px;cursor:pointer}
#rf-widget button.rf-btn:hover{background:#737373}
#rf-widget .rf-seg button.rf-btn{padding:6px 4px;font-size:13px}
#rf-widget .rf-seg button.active{background:#b45309;border-color:#fbbf24}
#rf-widget .rf-readout{color:#d1d5db;font-size:13px;line-height:1.55}
#rf-widget .rf-readout b{color:#fbbf24}
#rf-widget .rf-msg{color:#fca5a5}
#rf-widget .rf-bars{display:grid;grid-template-columns:34px 1fr 52px;gap:3px 6px;align-items:center;font-size:13px;margin-top:6px}
#rf-widget .rf-bar{height:10px;border-radius:2px}
#rf-widget .rf-track{position:relative;height:10px;background:#111827;border-radius:2px}
#rf-widget .rf-note{margin-top:6px;color:#9ca3af;font-size:12px;line-height:1.35}
#rf-widget .rf-hidden{display:none!important}
html[data-sa-mode="expt"] #rf-widget .rf-workspace{grid-template-rows:minmax(0,1fr) minmax(0,1.1fr);align-items:stretch!important;justify-content:stretch!important;align-content:stretch}
html[data-sa-mode="expt"] #rf-widget .rf-panel{min-height:0}
html[data-sa-mode="expt"] #rf-widget .rf-cv,html[data-sa-mode="expt"] #rf-widget .rf-wide .rf-cv{flex:1 1 auto;height:auto}
@media(max-width:820px){
  #rf-widget .rf-workspace{grid-template-columns:1fr}
  html[data-sa-mode="expt"] #rf-widget .rf-workspace{grid-template-rows:none}
  html[data-sa-mode="expt"] #rf-widget .rf-cv,html[data-sa-mode="expt"] #rf-widget .rf-wide .rf-cv{flex:none;height:250px}
}
</style>
<div class="rf-title app-title"><b>Surface and interface waves: one recipe, four settings</b><span style="display:none">. </span><div class="rf-hint">Pick a setting &middot; change the media &middot; watch the speed, the orbit and H/V respond</div></div>
<div class="rf-workspace app-stage">
  <section class="rf-panel rf-wide"><div class="rf-panel-title" id="rf-wave-title">Wavefield</div><div class="rf-caption" id="rf-wave-cap"></div><div class="rf-cv"><canvas id="rf-wave"></canvas></div>
    <div class="rf-legend"><button id="rf-play" type="button" class="rf-btn">Play</button>
      <label class="rf-key" style="color:#f59e0b"><input type="checkbox" id="rf-show-p" checked><i class="rf-sw"></i>P part</label>
      <label class="rf-key" style="color:#60a5fa"><input type="checkbox" id="rf-show-s" checked><i class="rf-sw"></i>SV part</label>
      <span style="color:#9ca3af">Switch one off: neither part alone satisfies the boundary conditions.</span></div></section>
  <section class="rf-panel"><div class="rf-panel-title" id="rf-sec-title">Secular function</div><div class="rf-caption" id="rf-sec-cap"></div><div class="rf-cv"><canvas id="rf-sec"></canvas></div><div class="rf-legend" id="rf-sec-legend"></div></section>
  <section class="rf-panel"><div class="rf-panel-title">Displacement with depth</div><div class="rf-caption" id="rf-prof-cap">Drag to choose the depth of the orbit.</div><div class="rf-cv"><canvas id="rf-prof"></canvas></div>
    <div class="rf-legend"><span class="rf-key" style="color:#e879f9"><i class="rf-sw"></i>u<sub>x</sub> amplitude</span><span class="rf-key" style="color:#22d3ee"><i class="rf-sw"></i>u<sub>z</sub> amplitude</span><span class="rf-key" style="color:#fbbf24"><i class="rf-sw rf-dash"></i>orbit depth</span></div></section>
  <section class="rf-panel"><div class="rf-panel-title">Particle orbit</div><div class="rf-caption" id="rf-orbit-cap"></div><div class="rf-cv"><canvas id="rf-orbit"></canvas></div>
    <div class="rf-legend"><span class="rf-key" style="color:#60a5fa"><i class="rf-sw"></i>path over one period</span><span class="rf-key" style="color:#fbbf24">&#9679; now</span></div></section>
</div>
<div class="rf-controls app-controls">
  <section class="rf-group"><div class="rf-group-title">Setting</div>
    <div class="rf-seg" id="rf-kind">
      <button type="button" class="rf-btn" data-kind="half">Half-space</button>
      <button type="button" class="rf-btn" data-kind="layer">Layer</button>
      <button type="button" class="rf-btn" data-kind="scholte">Scholte</button>
      <button type="button" class="rf-btn" data-kind="stoneley">Stoneley</button>
    </div>
    <select id="rf-preset"></select>
    <div class="rf-note" id="rf-kind-note"></div>
  </section>
  <section class="rf-group" id="rf-g1"><div class="rf-group-title" id="rf-g1-title">Layer</div>
    <label class="rf-row"><span>α</span><input id="rf-a1" type="range" min="0" max="1000" step="1"><span id="rf-a1-v" class="rf-value"></span></label>
    <label class="rf-row" id="rf-row-b1"><span>β</span><input id="rf-b1" type="range" min="0" max="1000" step="1"><span id="rf-b1-v" class="rf-value"></span></label>
    <label class="rf-row"><span>ρ</span><input id="rf-r1" type="range" min="0" max="1000" step="1"><span id="rf-r1-v" class="rf-value"></span></label>
    <label class="rf-row" id="rf-row-H"><span>H</span><input id="rf-H" type="range" min="0" max="1000" step="1"><span id="rf-H-v" class="rf-value"></span></label>
  </section>
  <section class="rf-group" id="rf-g2"><div class="rf-group-title" id="rf-g2-title">Half-space</div>
    <label class="rf-row"><span>α</span><input id="rf-a2" type="range" min="0" max="1000" step="1"><span id="rf-a2-v" class="rf-value"></span></label>
    <label class="rf-row"><span>β</span><input id="rf-b2" type="range" min="0" max="1000" step="1"><span id="rf-b2-v" class="rf-value"></span></label>
    <label class="rf-row"><span>ρ</span><input id="rf-r2" type="range" min="0" max="1000" step="1"><span id="rf-r2-v" class="rf-value"></span></label>
  </section>
  <section class="rf-group"><div class="rf-group-title">Observe</div>
    <label class="rf-row" id="rf-row-logf"><span>f / f₀</span><input id="rf-logf" type="range" min="0" max="1000" step="1"><span id="rf-logf-v" class="rf-value"></span></label>
    <label class="rf-row"><span>z / λ</span><input id="rf-depth" type="range" min="0" max="1000" step="1"><span id="rf-depth-v" class="rf-value"></span></label>
    <label class="rf-row"><span>phase</span><input id="rf-phase" type="range" min="0" max="1000" step="1"><span id="rf-phase-v" class="rf-value"></span></label>
  </section>
  <section class="rf-group"><div class="rf-group-title">Shear or compressional?</div>
    <label class="rf-key" style="color:#fb923c;display:flex"><input type="checkbox" id="rf-ghost-a" checked><i class="rf-sw rf-dash"></i>every α +10&nbsp;% <span id="rf-ga-v" class="rf-value" style="margin-left:auto"></span></label>
    <label class="rf-key" style="color:#a78bfa;display:flex;margin-top:4px"><input type="checkbox" id="rf-ghost-b" checked><i class="rf-sw rf-dash"></i>every β +10&nbsp;% <span id="rf-gb-v" class="rf-value" style="margin-left:auto"></span></label>
    <div class="rf-note">Sensitivity &part;ln c / &part;ln m: how many percent c changes per percent change of m.</div>
    <div class="rf-bars" id="rf-bars"></div>
    <div class="rf-note" id="rf-sum"></div>
  </section>
  <section class="rf-group"><div class="rf-group-title">Readout</div><div id="rf-readout" class="rf-readout">Computing&hellip;</div>
    <button id="rf-reset" type="button" class="rf-btn" style="margin-top:8px">Reset setting</button></section>
</div>
</div>
<script>
(function(){
if (window.__rfCleanup) window.__rfCleanup()
const ctl = new AbortController()
window.__rfCleanup = () => ctl.abort()
const SIG = {signal: ctl.signal}

const par = currentScript.previousElementSibling
const by = id => par.querySelector('#'+id)
const stage = par.querySelector('.app-stage')
const DPR = Math.min(window.devicePixelRatio || 1, 2)
const PRESETS = $(presets)
const state = $(state)
let R = window.__rfLast || null
let showP = true, showS = true, ghostA = true, ghostB = true
let playing = false, frame = 0, lastT = 0, rafId = null

// ---------- controls: every slider is 0..1000 mapped onto a linear or log range
const SPEC = {
  a1: [1, 10, true, v => v.toFixed(2) + ' km/s'], b1: [0.05, 6, true, v => v.toFixed(3) + ' km/s'],
  r1: [0.5, 10, true, v => v.toFixed(2) + ' g/cm³'], H: [0.005, 100, true, fmtLen],
  a2: [1, 10, true, v => v.toFixed(2) + ' km/s'], b2: [0.05, 6, true, v => v.toFixed(3) + ' km/s'],
  r2: [0.5, 10, true, v => v.toFixed(2) + ' g/cm³'],
  logf: [-1, 1, false, v => Math.pow(10, v).toFixed(2)],
  depth: [0, 1, false, v => v.toFixed(3) + ' λ'],
  phase: [0, 1, false, v => (360 * v).toFixed(0) + '°'],
}
function twoSided(){ return state.kind === 'scholte' || state.kind === 'stoneley' }
function range(key){ const s = SPEC[key]; return key === 'depth' && twoSided() ? [-1, 1, false] : s }
function toPos(key, v){ const [lo, hi, lg] = range(key); const t = lg ? Math.log(v / lo) / Math.log(hi / lo) : (v - lo) / (hi - lo); return Math.round(1000 * Math.min(1, Math.max(0, t))) }
function fromPos(key, p){ const [lo, hi, lg] = range(key); const t = p / 1000; return lg ? lo * Math.pow(hi / lo, t) : lo + (hi - lo) * t }
function fmtLen(v){ return v < 1 ? (1000 * v).toFixed(0) + ' m' : v.toFixed(v < 10 ? 2 : 1) + ' km' }

// same limits as the Julia boundary cell: solids keep α/β ≥ √2 (Poisson ratio ≥ 0)
function clampState(){
  for (const k of ['a1','b1','r1','a2','b2','r2','H','logf','phase']){ const [lo, hi] = SPEC[k]; if (!(k === 'b1' && state.kind === 'scholte')) state[k] = Math.min(hi, Math.max(lo, state[k])) }
  const [dlo] = range('depth'); state.depth = Math.min(1, Math.max(dlo, state.depth))
  if (state.kind !== 'scholte') state.b1 = Math.min(state.b1, state.a1 * Math.SQRT1_2)
  state.b2 = Math.min(state.b2, state.a2 * Math.SQRT1_2)
}

const TITLES = {
  half: ['', 'Half-space', 'Rayleigh waves on a free surface: non-dispersive, c = β · g(α/β).'],
  layer: ['Layer', 'Half-space', 'A layer over a half-space: dispersive Rayleigh modes, f₀ = β₁/4H.'],
  scholte: ['Fluid', 'Solid', 'Scholte waves on a fluid–solid interface (e.g. the seafloor).'],
  stoneley: ['Upper solid', 'Lower solid', 'Stoneley waves on a welded solid–solid interface: rare.'],
}
function syncControls(){
  clampState()
  for (const k of Object.keys(SPEC)){ by('rf-' + k).value = toPos(k, state[k]); by('rf-' + k + '-v').textContent = SPEC[k][3](state[k]) }
  for (const b of par.querySelectorAll('#rf-kind button')) b.classList.toggle('active', b.dataset.kind === state.kind)
  const t = TITLES[state.kind]
  by('rf-g1').classList.toggle('rf-hidden', state.kind === 'half')
  by('rf-g1-title').textContent = t[0]; by('rf-g2-title').textContent = t[1]; by('rf-kind-note').textContent = t[2]
  by('rf-row-b1').classList.toggle('rf-hidden', state.kind === 'scholte')
  by('rf-row-H').classList.toggle('rf-hidden', state.kind !== 'layer')
  by('rf-row-logf').classList.toggle('rf-hidden', state.kind !== 'layer')
  const sel = by('rf-preset')
  if (sel.dataset.kind !== state.kind){
    sel.dataset.kind = state.kind
    sel.innerHTML = PRESETS[state.kind].map((p, i) => '<option value="' + i + '">' + p.name + '</option>').join('')
  }
}
let emitTimer = null
function emit(now){
  clearTimeout(emitTimer)
  const go = () => { par.value = Object.assign({}, state); par.dispatchEvent(new CustomEvent('input')) }
  now ? go() : (emitTimer = setTimeout(go, 60))
}
function applyPreset(kind, i){
  const p = PRESETS[kind][i]
  Object.assign(state, {kind, a1: p.a1, b1: p.b1, r1: p.r1, a2: p.a2, b2: p.b2, r2: p.r2, H: p.H, mode: 1, logf: 0, depth: 0})
  syncControls(); by('rf-preset').value = String(i); emit(true)
}

// Local controls must not reach Pluto's own bond listener (capture phase, stopped here).
function onControl(e){
  if (e.target === par) return          // our own update for Pluto: let it through
  const id = e.target.id || ''
  if (id.startsWith('rf-')) e.stopImmediatePropagation()
  const key = id.slice(3)
  if (key in SPEC){
    state[key] = fromPos(key, Number(e.target.value))
    syncControls(); emit(e.type === 'change')
  } else if (id === 'rf-preset' && e.type === 'change'){
    applyPreset(state.kind, Number(e.target.value))
  } else if (id === 'rf-show-p'){ showP = e.target.checked; drawWave() }
  else if (id === 'rf-show-s'){ showS = e.target.checked; drawWave() }
  else if (id === 'rf-ghost-a'){ ghostA = e.target.checked; drawAll() }
  else if (id === 'rf-ghost-b'){ ghostB = e.target.checked; drawAll() }
}
par.addEventListener('input', onControl, true)
par.addEventListener('change', onControl, true)
for (const b of par.querySelectorAll('#rf-kind button')) b.addEventListener('click', () => applyPreset(b.dataset.kind, 0))
by('rf-reset').addEventListener('click', () => applyPreset(state.kind, Number(by('rf-preset').value) || 0))
by('rf-play').addEventListener('click', () => { playing = !playing; by('rf-play').textContent = playing ? 'Pause' : 'Play'; lastT = 0; drawOrbit() })

// ---------- canvases, refit to their boxes (the stage has a definite height in Expt)
const CV = ['rf-wave', 'rf-sec', 'rf-prof', 'rf-orbit'].map(id => ({cv: by(id), w: 0, h: 0, x: null}))
function layout(){
  let changed = false
  for (const c of CV){
    const box = c.cv.parentElement
    const w = Math.max(160, Math.floor(box.clientWidth)), h = Math.max(140, Math.floor(box.clientHeight))
    if (w !== c.w || h !== c.h){
      c.w = w; c.h = h
      c.cv.width = Math.round(w * DPR); c.cv.height = Math.round(h * DPR)
      c.cv.style.width = w + 'px'; c.cv.style.height = h + 'px'
      c.x = c.cv.getContext('2d'); c.x.setTransform(DPR, 0, 0, DPR, 0, 0)
      changed = true
    }
  }
  return changed
}
function clear(c){ c.x.fillStyle = '#000'; c.x.fillRect(0, 0, c.w, c.h); return c.x }
function note(c, text){ const x = clear(c); x.fillStyle = '#9ca3af'; x.font = '13px sans-serif'; x.textAlign = 'center'; x.fillText(text, c.w / 2, c.h / 2); x.textAlign = 'left' }
function fmt(v){ const a = Math.abs(v); return a === 0 ? '0' : a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? String(Number(v.toFixed(2))) : String(Number(v.toPrecision(2))) }
function niceStep(span, n){ const raw = span / Math.max(1, n), p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p; return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p }
function ticks(a, b, lg, n){
  const out = []
  if (lg){
    // thin the mantissas until the labels fit: 1-2-5, then 1-3, then decades only
    for (const ms of [[1, 2, 5], [1, 3], [1]]){
      out.length = 0
      for (let e = Math.floor(Math.log10(a)); e <= Math.ceil(Math.log10(b)); e++) for (const m of ms){ const v = m * Math.pow(10, e); if (v >= a * 0.999 && v <= b * 1.001) out.push(v) }
      if (out.length <= 2 * n) break
    }
    return out
  }
  const s = niceStep(b - a, n); for (let v = Math.ceil(a / s) * s; v <= b + 1e-9 * s; v += s) out.push(Math.abs(v) < 1e-9 * s ? 0 : v)
  return out
}
// axes in box B = {l, r, t, b}; returns the data-to-pixel maps
function axes(x, B, xl, xr, yl, yh, o){
  const tx = v => o.xlog ? (Math.log(v) - Math.log(xl)) / (Math.log(xr) - Math.log(xl)) : (v - xl) / (xr - xl)
  const ty = v => o.ylog ? (Math.log(v) - Math.log(yl)) / (Math.log(yh) - Math.log(yl)) : (v - yl) / (yh - yl)
  const X = v => B.l + tx(v) * (B.r - B.l)
  const Y = v => o.yflip ? B.t + ty(v) * (B.b - B.t) : B.b - ty(v) * (B.b - B.t)
  x.lineWidth = 1; x.font = '12px sans-serif'
  for (const v of ticks(xl, xr, o.xlog, o.nx || 5)){
    const p = X(v); x.strokeStyle = '#1f2937'; x.beginPath(); x.moveTo(p, B.t); x.lineTo(p, B.b); x.stroke()
    if (!o.noxlab){ const s = fmt(v), w = x.measureText(s).width; x.fillStyle = '#9ca3af'; x.fillText(s, Math.min(B.r - w, Math.max(B.l, p - w / 2)), B.b + 14) }
  }
  for (const v of ticks(yl, yh, o.ylog, o.ny || 4)){
    const p = Y(v); x.strokeStyle = '#1f2937'; x.beginPath(); x.moveTo(B.l, p); x.lineTo(B.r, p); x.stroke()
    const s = fmt(v), w = x.measureText(s).width; x.fillStyle = '#9ca3af'; x.fillText(s, Math.max(2, B.l - 5 - w), Math.min(B.b, Math.max(B.t + 9, p + 4)))
  }
  x.strokeStyle = '#374151'; x.strokeRect(B.l, B.t, B.r - B.l, B.b - B.t)
  x.fillStyle = '#9ca3af'
  if (o.xlabel){ const w = x.measureText(o.xlabel).width; x.fillText(o.xlabel, (B.l + B.r - w) / 2, B.b + 28) }
  if (o.ylabel){ x.save(); x.translate(12, (B.t + B.b) / 2); x.rotate(-Math.PI / 2); x.textAlign = 'center'; x.fillText(o.ylabel, 0, 0); x.restore() }
  return {X, Y}
}
function poly(x, xs, ys, X, Y, color, width, dash){
  x.strokeStyle = color; x.lineWidth = width; x.setLineDash(dash || []); x.beginPath(); let pen = false
  for (let i = 0; i < xs.length; i++){
    const a = xs[i], b = ys[i]
    if (a === null || b === null || !isFinite(a) || !isFinite(b)){ pen = false; continue }
    pen ? x.lineTo(X(a), Y(b)) : x.moveTo(X(a), Y(b)); pen = true
  }
  x.stroke(); x.setLineDash([])
}
function vline(x, px, B, color, dash, label){
  x.strokeStyle = color; x.lineWidth = 1.6; x.setLineDash(dash ? [5, 4] : []); x.beginPath(); x.moveTo(px, B.t); x.lineTo(px, B.b); x.stroke(); x.setLineDash([])
  if (label){ x.fillStyle = color; x.font = '12px sans-serif'; const w = x.measureText(label).width; x.fillText(label, Math.min(B.r - w - 2, px + 4), B.t + 13) }
}
function hline(x, py, B, color, dash, label){
  x.strokeStyle = color; x.lineWidth = 1.4; x.setLineDash(dash ? [5, 4] : []); x.beginPath(); x.moveTo(B.l, py); x.lineTo(B.r, py); x.stroke(); x.setLineDash([])
  if (label){ x.fillStyle = color; x.font = '12px sans-serif'; x.fillText(label, B.l + 4, py - 4) }
}
function arrow(x, sx, sy, dx, dy, color){
  const len = Math.hypot(dx, dy); if (len < 0.6) return
  const ux = dx / len, uy = dy / len, ah = Math.min(5, 0.45 * len), aw = 0.55 * ah
  x.strokeStyle = color; x.lineWidth = 1.5; x.beginPath(); x.moveTo(sx, sy); x.lineTo(sx + dx, sy + dy); x.stroke()
  const tx = sx + dx, ty = sy + dy, bx = tx - ux * ah, by2 = ty - uy * ah
  x.fillStyle = color; x.beginPath(); x.moveTo(tx, ty); x.lineTo(bx - uy * aw, by2 + ux * aw); x.lineTo(bx + uy * aw, by2 - ux * aw); x.closePath(); x.fill()
}
function shade(x, B, Y, zlo, zhi){
  for (const g of R.regions){
    const a = Math.max(g.z0, zlo), b = Math.min(g.z1, zhi); if (b <= a) continue
    const ya = Y(a), yb = Y(b); x.fillStyle = g.color; x.fillRect(B.l, Math.min(ya, yb), B.r - B.l, Math.abs(yb - ya))
    x.fillStyle = '#6b7280'; x.font = '12px sans-serif'; x.fillText(g.name, B.r - x.measureText(g.name).width - 6, Math.min(ya, yb) + 14)
  }
}
function ifaces(x, B, Y, zlo, zhi){
  for (const L of R.lines){ if (L.z < zlo || L.z > zhi) continue; x.strokeStyle = L.color; x.lineWidth = 2; x.setLineDash(L.dash ? [6, 4] : []); x.beginPath(); x.moveTo(B.l, Y(L.z)); x.lineTo(B.r, Y(L.z)); x.stroke(); x.setLineDash([]) }
}

// ---------- panels
function drawWave(){
  const c = CV[0]
  if (!R) return note(c, 'Computing…')
  if (!R.ok){ by('rf-wave-cap').textContent = ''; return note(c, R.msg) }
  const x = clear(c), wf = R.wf, B = {l: 50, r: c.w - 10, t: 8, b: c.h - 34}
  const z0 = wf.zs[0], z1 = wf.zs[wf.zs.length - 1]
  const A = axes(x, B, wf.xs[0], wf.xs[wf.xs.length - 1], z0, z1, {yflip: true, xlabel: 'distance x (' + R.unit + ')', ylabel: 'depth z (' + R.unit + ')', nx: Math.max(3, Math.floor((B.r - B.l) / 80))})
  shade(x, B, A.Y, z0, z1); ifaces(x, B, A.Y, z0, z1)
  const nx = wf.nx, nz = wf.nz, sc = 0.95 * Math.min((B.r - B.l) / (nx - 1), (B.b - B.t) / (nz - 1))
  const col = showP && showS ? '#e5e7eb' : showP ? '#f59e0b' : '#60a5fa'
  for (let ix = 0; ix < nx; ix++) for (let iz = 0; iz < nz; iz++){
    const k = (frame * nx + ix) * nz + iz
    const ux = (showP ? wf.uxP[k] : 0) + (showS ? wf.uxS[k] : 0), uz = (showP ? wf.uzP[k] : 0) + (showS ? wf.uzS[k] : 0)
    arrow(x, A.X(wf.xs[ix]), A.Y(wf.zs[iz]), ux * sc, uz * sc, col)
  }
  const what = {half: 'Rayleigh wave', layer: 'Rayleigh mode R' + (R.mode - 1), scholte: 'Scholte wave', stoneley: 'Stoneley wave'}[R.kind]
  by('rf-wave-title').textContent = 'Wavefield: ' + what + ' travelling to the right'
  by('rf-wave-cap').textContent = 'c = ' + R.c.toFixed(3) + ' km/s, wavelength λ = ' + fmtLen(R.lamkm) + (R.kind === 'layer' ? ', f = ' + fmt(R.f) + ' Hz' : '') + '. Arrows: displacement at one instant.'
}

let secGeom = null
function drawSec(){
  const c = CV[1]
  if (!R) return note(c, 'Computing…')
  const x = clear(c)
  if (R.kind !== 'layer'){
    const S = R.sec, B = {l: 46, r: c.w - 10, t: 10, b: c.h - 34}
    const lim = 1.12 * Math.max(1e-12, ...S.d.map(v => Math.abs(v || 0)))   // fit the curve, however small
    const cR = 1.15 * S.c[S.c.length - 1]   // room to the right for the β +10 % zero
    const A = axes(x, B, S.c[0], cR, -lim, lim, {xlabel: 'phase velocity c (km/s)', ylabel: 'normalised det'})
    hline(x, A.Y(0), B, '#4b5563')
    poly(x, S.c, S.d, A.X, A.Y, '#f59e0b', 2.3)
    for (const L of S.lines) if (L.v <= cR) vline(x, A.X(L.v), B, '#9ca3af', true, L.label)
    if (R.ok){
      if (ghostA && R.ghost.cA !== null) vline(x, A.X(R.ghost.cA), B, '#fb923c', true)
      if (ghostB && R.ghost.cB !== null && R.ghost.cB <= cR) vline(x, A.X(R.ghost.cB), B, '#a78bfa', true)
      for (const r of S.roots){ x.fillStyle = '#fbbf24'; x.beginPath(); x.arc(A.X(r), A.Y(0), 5, 0, 2 * Math.PI); x.fill() }
    }
    by('rf-sec-title').textContent = 'Secular function D(c)'
    by('rf-sec-cap').textContent = 'Its zero (yellow) is the wave speed; dashed: the zero with every α or every β raised by 10 %.'
    by('rf-sec-legend').innerHTML = '<span class="rf-key" style="color:#f59e0b"><i class="rf-sw"></i>D(c)</span><span class="rf-key" style="color:#fb923c"><i class="rf-sw rf-dash"></i>α +10 %</span><span class="rf-key" style="color:#a78bfa"><i class="rf-sw rf-dash"></i>β +10 %</span>'
    secGeom = null
    return
  }
  const D = R.disp, f = D.f, gap = 34
  const Ht = Math.round((c.h - 10 - 30 - gap) * 0.6)
  const B1 = {l: 46, r: c.w - 10, t: 10, b: 10 + Ht}, B2 = {l: 46, r: c.w - 10, t: 10 + Ht + gap, b: c.h - 30}
  const cs = []; for (const row of D.C.concat(D.CA, D.CB)) for (const v of row) if (v !== null) cs.push(v)
  const cmin = 0.95 * Math.min(...cs, D.lines[0].v), cmax = 1.03 * Math.max(D.lines[1].v, D.lines[0].v)
  const A = axes(x, B1, f[0], f[f.length - 1], cmin, cmax, {xlog: true, noxlab: true, ylabel: 'c (km/s)'})
  for (const L of D.lines) hline(x, A.Y(L.v), B1, '#6b7280', true, L.label)
  const m = R.ok ? R.mode - 1 : 0
  for (let k = D.C.length - 1; k >= 0; k--) poly(x, f, D.C[k], A.X, A.Y, k === m ? '#fbbf24' : '#64748b', k === m ? 2.6 : 1.4)
  if (ghostA) poly(x, f, D.CA[m], A.X, A.Y, '#fb923c', 1.6, [5, 4])
  if (ghostB) poly(x, f, D.CB[m], A.X, A.Y, '#a78bfa', 1.6, [5, 4])
  poly(x, f, D.Ug[m], A.X, A.Y, '#fde68a', 1.2, [1, 3])
  vline(x, A.X(D.f0), B1, '#34d399', true, 'f₀')
  if (R.ok) vline(x, A.X(R.f), B1, '#e5e7eb', false)
  const hv = D.HV.filter(v => v !== null && v > 0)
  const hlo = hv.length ? Math.max(0.02, Math.min(...hv) / 1.5) : 0.1, hhi = hv.length ? Math.min(200, Math.max(2, ...hv) * 1.5) : 10
  const A2 = axes(x, B2, f[0], f[f.length - 1], hlo, hhi, {xlog: true, ylog: true, xlabel: 'frequency f (Hz)', ylabel: 'H/V'})
  hline(x, A2.Y(1), B2, '#4b5563', true)
  poly(x, f, D.HV, A2.X, A2.Y, '#34d399', 2.2)
  vline(x, A2.X(D.f0), B2, '#34d399', true)
  if (R.ok) vline(x, A2.X(R.f), B2, '#e5e7eb', false)
  secGeom = {B1, B2, f, A, A2}
  by('rf-sec-title').textContent = 'Dispersion c(f) and H/V(f)'
  by('rf-sec-cap').textContent = 'Click a curve to pick a mode, drag to change the frequency. Dotted: group velocity.'
  by('rf-sec-legend').innerHTML = '<span class="rf-key" style="color:#fbbf24"><i class="rf-sw"></i>chosen mode</span><span class="rf-key" style="color:#fb923c"><i class="rf-sw rf-dash"></i>α +10 %</span><span class="rf-key" style="color:#a78bfa"><i class="rf-sw rf-dash"></i>β +10 %</span><span class="rf-key" style="color:#34d399"><i class="rf-sw"></i>H/V of R0</span>'
}
function secPointer(e, down){
  if (!secGeom || state.kind !== 'layer') return
  const r = CV[1].cv.getBoundingClientRect(), px = e.clientX - r.left, py = e.clientY - r.top
  const {B1, B2, f} = secGeom
  if (px < B1.l || px > B1.r || py < B1.t || py > B2.b) return
  const t = (px - B1.l) / (B1.r - B1.l), fsel = f[0] * Math.pow(f[f.length - 1] / f[0], t)
  state.logf = Math.log10(fsel / R.disp.f0)
  if (down && py <= B1.b){
    let best = -1, bd = 14, i = Math.round(t * (f.length - 1))
    R.disp.C.forEach((row, k) => { const v = row[i]; if (v !== null){ const d = Math.abs(secGeom.A.Y(v) - py); if (d < bd){ bd = d; best = k } } })
    if (best >= 0) state.mode = best + 1
  }
  syncControls(); emit(false)
}
let secDrag = false
CV[1].cv.addEventListener('pointerdown', e => { secDrag = true; CV[1].cv.setPointerCapture(e.pointerId); secPointer(e, true) }, SIG)
CV[1].cv.addEventListener('pointermove', e => { if (secDrag) secPointer(e, false) }, SIG)
CV[1].cv.addEventListener('pointerup', () => { secDrag = false; emit(true) }, SIG)

let profGeom = null
function drawProf(){
  const c = CV[2]
  if (!R) return note(c, 'Computing…')
  if (!R.ok) return note(c, 'No wave to show.')
  const x = clear(c), P = R.prof, B = {l: 50, r: c.w - 10, t: 8, b: c.h - 34}
  const z0 = P.z[0], z1 = P.z[P.z.length - 1]
  const A = axes(x, B, -1.1, 1.1, z0, z1, {yflip: true, xlabel: 'normalised amplitude', ylabel: 'depth z (' + R.unit + ')'})
  shade(x, B, A.Y, z0, z1); ifaces(x, B, A.Y, z0, z1)
  vline(x, A.X(0), B, '#4b5563')
  poly(x, P.U, P.z, A.X, A.Y, '#e879f9', 2.3)
  poly(x, P.W, P.z, A.X, A.Y, '#22d3ee', 2.3)
  if (P.node !== null && P.node >= z0 && P.node <= z1){
    x.fillStyle = '#e879f9'; x.beginPath(); x.arc(A.X(0), A.Y(P.node), 4.5, 0, 2 * Math.PI); x.fill()
    x.font = '12px sans-serif'; x.fillText('u_x = 0 at ' + R.node.toFixed(3) + ' λ', A.X(0) + 8, A.Y(P.node) - 6)
  }
  hline(x, A.Y(R.zsel), B, '#fbbf24', true)
  profGeom = {B, z0, z1}
}
function profPointer(e){
  if (!profGeom || !R || !R.ok) return
  const r = CV[2].cv.getBoundingClientRect(), py = e.clientY - r.top, {B, z0, z1} = profGeom
  const z = z0 + (Math.min(B.b, Math.max(B.t, py)) - B.t) / (B.b - B.t) * (z1 - z0)
  state.depth = z / R.lam
  syncControls(); emit(false)
}
let profDrag = false
CV[2].cv.addEventListener('pointerdown', e => { profDrag = true; CV[2].cv.setPointerCapture(e.pointerId); profPointer(e) }, SIG)
CV[2].cv.addEventListener('pointermove', e => { if (profDrag) profPointer(e) }, SIG)
CV[2].cv.addEventListener('pointerup', () => { profDrag = false; emit(true) }, SIG)

function drawOrbit(){
  const c = CV[3]
  if (!R) return note(c, 'Computing…')
  if (!R.ok) return note(c, 'No wave to show.')
  const x = clear(c), O = R.orbit, cx = c.w / 2, cy = c.h / 2 + 4, S = 0.38 * Math.min(c.w, c.h)
  const X = v => cx + v * S, Y = v => cy + v * S          // z grows downward on screen too
  x.strokeStyle = '#374151'; x.lineWidth = 1; x.beginPath(); x.moveTo(cx - 1.15 * S, cy); x.lineTo(cx + 1.15 * S, cy); x.moveTo(cx, cy - 1.15 * S); x.lineTo(cx, cy + 1.15 * S); x.stroke()
  poly(x, O.x, O.z, X, Y, '#60a5fa', 2.3)
  const n = O.x.length, i0 = Math.round(n / 8), i1 = i0 + 4
  arrow(x, X(O.x[i0]), Y(O.z[i0]), X(O.x[i1]) - X(O.x[i0]), Y(O.z[i1]) - Y(O.z[i0]), '#60a5fa')
  let sx = O.sx, sz = O.sz
  if (playing && R.wf){ const k = Math.round(frame / R.wf.nframes * (n - 1)); sx = O.x[k]; sz = O.z[k] }
  x.fillStyle = '#fbbf24'; x.beginPath(); x.arc(X(sx), Y(sz), 5.5, 0, 2 * Math.PI); x.fill()
  x.fillStyle = '#9ca3af'; x.font = '12px sans-serif'
  const lab = 'travel direction →'; x.fillText(lab, c.w - x.measureText(lab).width - 8, c.h - 8)
  x.fillText('↓ down', 8, c.h - 8)
  x.fillStyle = O.sense === 'retrograde' ? '#fbbf24' : '#34d399'; x.font = '600 14px sans-serif'; x.fillText(O.sense, 10, 20)
  by('rf-orbit-cap').textContent = 'At z = ' + R.depth.toFixed(3) + ' λ (' + fmtLen(Math.abs(R.depth * R.lamkm)) + (R.depth < 0 ? ' above' : '') + '), x = 0; amplitude ' + (100 * O.amp).toFixed(0) + ' % of the largest.'
}

function drawSide(){
  if (!R) return
  const ro = by('rf-readout')
  if (!R.ok){ ro.innerHTML = '<span class="rf-msg">' + R.msg + '</span>'; by('rf-bars').innerHTML = ''; by('rf-sum').textContent = ''; by('rf-ga-v').textContent = ''; by('rf-gb-v').textContent = ''; return }
  const name = {half: 'c<sub>R</sub>', layer: 'c', scholte: 'c<sub>Sch</sub>', stoneley: 'c<sub>St</sub>'}[R.kind]
  let h = name + ' = <b>' + R.c.toFixed(3) + ' km/s</b> = <b>' + (R.c / R.cref).toFixed(3) + '</b> ' + R.crefname + '<br>'
  if (R.kind === 'layer') h += 'mode <b>R' + (R.mode - 1) + '</b> of ' + R.nroots + ' at f = <b>' + fmt(R.f) + ' Hz</b> (f₀ = ' + fmt(R.f0) + ' Hz)<br>'
  const at = {half: 'at the surface', layer: 'at the surface', scholte: 'at the interface (solid side)', stoneley: 'at the interface'}[R.kind]
  h += 'V/H = <b>' + R.vh.toFixed(2) + '</b>, H/V = <b>' + R.hv.toFixed(2) + '</b> ' + at + '<br>'
  h += (R.node !== null ? 'u<sub>x</sub> changes sign at <b>' + R.node.toFixed(3) + ' λ</b>' : 'u<sub>x</sub> keeps its sign') + '; orbit at z: <b>' + R.orbit.sense + '</b>'
  ro.innerHTML = h
  const pct = v => (v >= 0 ? '+' : '') + (100 * v).toFixed(1) + ' %'
  by('rf-ga-v').textContent = R.ghost.cA !== null ? 'c ' + pct(R.ghost.dA) : 'no root'
  by('rf-gb-v').textContent = R.ghost.cB !== null ? 'c ' + pct(R.ghost.dB) : 'no root'
  const COL = {'α': '#fb923c', 'β': '#a78bfa', 'ρ': '#9ca3af'}
  const top = Math.max(1, ...R.sens.map(s => Math.abs(s.value)))
  by('rf-bars').innerHTML = R.sens.map(s => {
    const w = 50 * Math.abs(s.value) / top, left = s.value >= 0 ? 50 : 50 - w
    return '<span style="color:' + COL[s.kind] + '">' + s.label + '</span><div class="rf-track"><div class="rf-bar" style="position:absolute;left:' + left + '%;width:' + w + '%;background:' + COL[s.kind] + '"></div><div style="position:absolute;left:50%;top:-2px;bottom:-2px;width:1px;background:#4b5563"></div></div><span class="rf-value">' + s.value.toFixed(2) + '</span>'
  }).join('')
  by('rf-sum').innerHTML = 'Velocity entries sum to <b>' + R.sumv.toFixed(3) + '</b> = c/U = ' + R.cU.toFixed(3) + (R.kind === 'layer' ? ' (dispersive)' : ' (non-dispersive: exactly 1)') + '; density entries sum to 0.'
}
function drawAll(){ drawWave(); drawSec(); drawProf(); drawOrbit(); drawSide() }

window.addEventListener('rf-results', e => { R = e.detail; window.__rfLast = R; frame = Math.min(frame, R.ok ? R.wf.nframes - 1 : 0); drawAll() }, SIG)
function tick(now){
  if (!par.isConnected){ return }
  if (playing && R && R.ok && now - lastT > 80){ frame = (frame + 1) % R.wf.nframes; lastT = now; drawWave(); drawOrbit() }
  rafId = requestAnimationFrame(tick)
}
// Disconnect only once this widget has left the page: Pluto can resolve `invalidation`
// for a bond output while keeping its DOM, which would freeze the canvases' size.
const fitRO = new ResizeObserver(() => {
  if (!par.isConnected){ fitRO.disconnect(); if (rafId) cancelAnimationFrame(rafId); return }
  if (layout()) drawAll()
})
fitRO.observe(stage)
syncControls(); layout(); drawAll()
rafId = requestAnimationFrame(tick)
par.value = Object.assign({}, state)
})()
</script>
""")
	end

	const _rf_ready = true
end

# ╔═╡ 544d41b1-b6ce-4844-8269-a0251b391962
begin
	# The widget implementation lives in the Appendix. This explicit reference
	# makes a cold Pluto load wait for its definition.
	_rf_ready
	PlutoUI.WideCell(@bind _rf RayleighFunctionInput(); max_width = 1400)
end

# ╔═╡ 31a317e2-0864-4552-83ee-ba4138ab5c0b
begin
	# Convert the browser value at one boundary, with the same limits the widget applies;
	# everything downstream sees concrete numbers. The bond is `nothing` until a browser
	# has rendered the widget, so fall back to its defaults.
	rf_safe = _rf isa AbstractDict ? _rf : Base.get(RayleighFunctionInput())
	rf_in = let s = rf_safe, d = rf_preset_state(:half, 1)
		num(key, lo, hi) = clamp(Float64(get(s, key, d[key])), lo, hi)
		kind = string(get(s, "kind", "half"))
		kind in ("half", "layer", "scholte", "stoneley") || (kind = "half")
		a1, a2 = num("a1", 1.0, 10.0), num("a2", 1.0, 10.0)
		# solids keep α/β ≥ √2 (Poisson ratio ≥ 0); the Scholte fluid has no shear
		b1 = kind == "scholte" ? 0.0 : min(num("b1", 0.05, 6.0), a1 / sqrt(2))
		b2 = min(num("b2", 0.05, 6.0), a2 / sqrt(2))
		two_sided_kind = kind in ("scholte", "stoneley")
		(kind = kind,
		 m1 = Medium(a1, b1, num("r1", 0.5, 10.0)), m2 = Medium(a2, b2, num("r2", 0.5, 10.0)),
		 H = num("H", 0.005, 100.0), logf = num("logf", -1.0, 1.0),
		 depth = num("depth", two_sided_kind ? -1.0 : 0.0, 1.0), phase = num("phase", 0.0, 1.0),
		 mode = clamp(round(Int, Float64(get(s, "mode", 1))), 1, 5))
	end
	rf_problem = build_problem(rf_in.kind, rf_in.m1, rf_in.m2, rf_in.H)
end

# ╔═╡ bd7562f0-cb14-4ddb-b8c7-ffc2a37b079e
let
	payload = rayleigh_payload(rf_problem, rf_in.logf, rf_in.depth, rf_in.phase, rf_in.mode)
	# Stored on window as well, so a widget that renders after this cell still finds it.
	HTML("<script>window.__rfLast = " * to_json(payload) *
		 "; window.dispatchEvent(new CustomEvent('rf-results', {detail: window.__rfLast}))</script>")
end

# ╔═╡ c760b9b7-1896-4657-aa30-7d8ddcb8aa40
let
	ω = 2π
	rows = String[]
	r5(x) = string(round(x; sigdigits = 5))
	check(name, got, want, ok) = push!(rows, "| $name | $got | $want | $(ok ? "✓" : "✗") |")

	# 1. Poisson half-space (α = √3 β exactly) against the textbook numbers
	β = 3.46
	Hs = HalfSpace(Medium(sqrt(3) * β, β, 2.7))
	cR = find_speeds(Hs, ω)[1]
	check("Poisson ``c_R/\\beta``", r5(cR / β), "``\\sqrt{2 - 2/\\sqrt3}`` = 0.91940", abs(cR / β - sqrt(2 - 2 / sqrt(3))) < 1e-8)
	check("same root from `rayleigh_function`", r5(rayleigh_speed(sqrt(3) * β, β)), r5(cR), abs(rayleigh_speed(sqrt(3) * β, β) - cR) < 1e-8)
	a = mode_amplitudes(Hs, cR, ω)
	s0 = field_at(Hs, a, cR, ω, 0.0)
	check("Poisson V/H at the surface", r5(abs(s0.W / s0.U)), "1.47", abs(abs(s0.W / s0.U) - 1.468) < 0.005)
	zs = range(0, cR, length = 4001)
	U, _ = mode_profile(Hs, cR, ω, zs)
	check("Poisson depth of ``u_x = 0``", r5(ux_node(zs, U) / cR) * " λ", "0.19 λ", abs(ux_node(zs, U) / cR - 0.19) < 0.005)
	check("retrograde at the surface", string(s0.U * s0.W > 0), "true", s0.U * s0.W > 0)

	# 2. tractions recomputed from the displacement by finite differences
	function traction_error(P, m, c, ω, z)
		a = mode_amplitudes(P, c, ω)
		λ, μ = lame(m)
		k, h = ω / c, 1e-6 * c
		f, fp, fm = field_at(P, a, c, ω, z), field_at(P, a, c, ω, z + h), field_at(P, a, c, ω, z - h)
		U′, W′ = (fp.U - fm.U) / 2h, (fp.W - fm.W) / 2h
		T, N = μ * (U′ + k * f.W), (λ + 2μ) * W′ - λ * k * f.U
		return max(abs(T - f.T), abs(N - f.N)) / max(abs(f.T), abs(f.N), 1e-300)
	end
	check("half-space: σ from ``\\partial u`` at z = 0.1λ", r5(traction_error(Hs, Hs.m, cR, ω, 0.1cR)), "≈ 0", traction_error(Hs, Hs.m, cR, ω, 0.1cR) < 1e-6)
	Lc = LayerOverHalfSpace(Medium(6.0, 3.46, 2.7), Medium(8.0, 4.49, 3.3), 35.0)
	ωL = 2π * f0_layer(Lc)
	cL = find_speeds(Lc, ωL)[1]
	check("layer: σ from ``\\partial u`` inside the layer", r5(traction_error(Lc, Lc.m1, cL, ωL, 20.0)), "≈ 0", traction_error(Lc, Lc.m1, cL, ωL, 20.0) < 1e-5)
	sL = field_at(Lc, mode_amplitudes(Lc, cL, ωL), cL, ωL, 0.0)
	tops = max(abs(sL.T), abs(sL.N)) / max(abs(sL.U), abs(sL.W)) / (lame(Lc.m1)[2] * ωL / cL)
	check("layer: free-surface traction", r5(tops), "≈ 0", tops < 1e-8)

	# 3. limits of the layer problem
	same = LayerOverHalfSpace(Hs.m, Hs.m, 35.0)
	dev = maximum(abs(find_speeds(same, 2π * f)[1] - cR) for f in (0.005, 0.03, 0.2))
	check("layer = half-space gives ``c_R`` at every f", r5(dev), "0", dev < 1e-8)
	clo = find_speeds(Lc, 2π * f0_layer(Lc) / 100)[1]
	chi = find_speeds(Lc, 2π * f0_layer(Lc) * 100)[1]
	cR2, cR1 = rayleigh_speed(8.0, 4.49), rayleigh_speed(6.0, 3.46)
	check("low f → ``c_R`` of the half-space", r5(clo), r5(cR2), abs(clo / cR2 - 1) < 0.01)
	check("high f → ``c_R`` of the layer", r5(chi), r5(cR1), abs(chi / cR1 - 1) < 1e-4)

	# 4. interface waves
	sol = Medium(1.8, 0.4, 1.8)
	cS0 = find_speeds(Scholte(Medium(1.5, 0.0, 1e-7), sol), ω)[1]
	check("Scholte, fluid density → 0", r5(cS0), r5(rayleigh_speed(1.8, 0.4)), abs(cS0 / rayleigh_speed(1.8, 0.4) - 1) < 1e-4)
	St = Stoneley(Medium(6.0, 3.46, 2.0), Medium(6.0, 3.46, 8.0))
	cSt = find_speeds(St, ω)[1]
	aSt = mode_amplitudes(St, cSt, ω)
	u, d = field_at(St, aSt, cSt, ω, -1e-12), field_at(St, aSt, cSt, ω, 0.0)
	jump = maximum(abs.((u.U - d.U, u.W - d.W, u.T - d.T, u.N - d.N))) / max(abs(d.T), abs(d.N))
	check("Stoneley: continuity at the interface", r5(jump), "≈ 0", jump < 1e-6)

	# 5. sensitivity sum rules
	sH = log_sensitivities(Hs, cR, ω)[1]
	check("half-space: ``\\sum \\partial\\ln c/\\partial\\ln v``", r5(sH[1] + sH[2]), "1", abs(sH[1] + sH[2] - 1) < 1e-5)
	sL6 = log_sensitivities(Lc, cL, ωL)
	δ = 1e-3
	cm, cp = nearest_speed(Lc, ωL * (1 - δ), cL), nearest_speed(Lc, ωL * (1 + δ), cL)
	cU = cL / (2δ * ωL / (ωL * (1 + δ) / cp - ωL * (1 - δ) / cm))
	sumv = sL6[1][1] + sL6[1][2] + sL6[2][1] + sL6[2][2]
	check("layer: ``\\sum \\partial\\ln c/\\partial\\ln v = c/U``", r5(sumv), r5(cU), abs(sumv - cU) < 2e-3)
	check("layer: ``\\sum \\partial\\ln c/\\partial\\ln \\rho``", r5(sL6[1][3] + sL6[2][3]), "0", abs(sL6[1][3] + sL6[2][3]) < 1e-5)

	table = Markdown.parse("| check | value | expected | |\n|:--|--:|--:|:-:|\n" * join(rows, "\n"))
	# wrapped, so the app shell files this cell (code and result) under Code, not Theory
	HTML("<div class=\"rf-checks\">" * repr(MIME"text/html"(), table) * "</div>")
end

# ╔═╡ 07ee82c2-d17b-4b85-9ef6-c107b76d0cfe
md"""
# Rayleigh Function
Rayleigh waves are coupled, evanescent P and SV waves that together leave a free surface without traction. The same recipe describes every wave guided by a flat boundary, and the experiment applies it to four settings:

- **Half-space**: the classical Rayleigh wave. For a given Poisson ratio its speed, the shape of its orbit and the depth where the orbit reverses are fixed numbers.
- **Layer over a half-space**: the wave becomes dispersive, higher modes appear, and the ratio of horizontal to vertical motion (H/V) peaks near the layer's resonance frequency.
- **Scholte**: the wave on a fluid–solid interface, such as the seafloor.
- **Stoneley**: the wave on a welded interface between two solids, which exists only for a narrow range of media.

In every setting you can raise all the P speeds, or all the S speeds, by 10 % and see which of the two the wave speed actually follows.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 318ac5d5-bd8b-417b-a919-b46ee555ca1a
md"""
### What to look for

- **Half-space, Poisson solid** (``\alpha = 6.0``, ``\beta = 3.46`` km/s): ``c_R = 0.919\beta``, the vertical amplitude at the surface is 1.47 times the horizontal one, and the horizontal motion changes sign at a depth of ``0.19\lambda``, the numbers in Shearer's *Introduction to Seismology*. The orbit is retrograde above that depth and prograde below it.
- **β, not α**: raising every ``\alpha`` by 10 % changes the speed by about 1 %; raising every ``\beta`` by 10 % changes it by about 8 %. The sensitivity bars tell the same story in every setting.
- **Layer over a half-space**: at low frequency the wave samples the half-space and travels near its Rayleigh speed; at high frequency it stays in the layer and approaches the layer's Rayleigh speed. In between, near ``f_0 = \beta_1/4H``, the soft-sediment preset has a sharp H/V peak.
- **Scholte**: on soft sediment the wave travels a little below the sediment's shear speed; on hard rock it travels just below the speed of sound in water and depends almost only on the water.
- **Stoneley**: with equal shear speeds and a fourfold density contrast the wave exists, at ``0.975\beta``; for crust against mantle it does not.
"""

# ╔═╡ e89ffc5a-8252-4186-8600-dd68f91b425d
md"""
## Inhomogeneous Plane Waves
In a homogeneous medium with P speed ``\alpha``, S speed ``\beta`` and density ``\rho``, a P potential ``\phi`` and an SV potential ``\psi`` travelling along ``x`` with horizontal slowness ``p = 1/c`` have the form

```math
\phi = f(z)\,e^{i\omega(t - p x)}, \qquad \psi = i\,g(z)\,e^{i\omega(t - p x)} ,
```

and the wave equation fixes their depth dependence: ``f'' = \omega^2 s_\alpha f`` and ``g'' = \omega^2 s_\beta g``, with ``s_v = p^2 - 1/v^2``. When the wave is slower than the body wave, ``c < v``, it cannot propagate vertically: ``s_v > 0``, and the solutions are ``e^{\mp\omega\eta z}`` with ``\eta = \sqrt{s_v}``. Such an inhomogeneous wave travels along the boundary and decays away from it.

A wave guided by a boundary is built from these waves on every side that extends to infinity; `decaying` returns them. Inside a layer of finite thickness both signs are allowed, and `layer_basis` uses the pair ``C(z)``, ``S(z)`` (cosh and sinh when ``s_v > 0``, cos and sin when ``s_v < 0``), which stay finite and independent through ``s_v = 0``. The factor ``i`` on ``\psi`` is a bookkeeping choice that makes every quantity below real.
"""

# ╔═╡ 6ada314c-d48f-4d01-8a53-afc06ce649cb
md"""
## Displacement and Traction
With ``z`` pointing down, the displacement is ``u_x = \partial_x\phi - \partial_z\psi`` and ``u_z = \partial_z\phi + \partial_x\psi``, and the tractions on a horizontal plane are

```math
\sigma_{xz} = \mu\,(\partial_x u_z + \partial_z u_x), \qquad \sigma_{zz} = \lambda\,\nabla\cdot\mathbf u + 2\mu\,\partial_z u_z ,
```

with ``\mu = \rho\beta^2`` and ``\lambda = \rho(\alpha^2 - 2\beta^2)`` (`lame`). A P wave with depth function ``f`` then has the displacement–stress vector

```math
\big(i u_x,\ u_z,\ i\sigma_{xz},\ \sigma_{zz}\big) = \Big(\omega p f,\ f',\ 2\mu\omega p f',\ \omega^2\rho\,(2\beta^2p^2 - 1)\, f\Big),
```

and an SV wave with depth function ``g``

```math
\big(i u_x,\ u_z,\ i\sigma_{xz},\ \sigma_{zz}\big) = \Big(g',\ \omega p g,\ \mu\omega^2\big(2p^2 - 1/\beta^2\big)\, g,\ 2\mu\omega p g'\Big),
```

all four entries real (`p_state`, `s_state`). The physical displacement is ``u_x = U\sin\theta``, ``u_z = W\cos\theta`` with ``\theta = \omega(t - x/c)``: the two components are a quarter period apart, which is why every particle moves on an ellipse.
"""

# ╔═╡ 1b2b6012-d997-4f43-8b58-87e8bc40acf3
md"""
## One Recipe for Every Boundary
Every setting in the experiment is solved the same way.

1. List the partial waves: P and SV waves that decay away from the boundary on each side that extends to infinity, plus both depth functions of each wave type inside a layer (`partial_states`).
2. List the boundary conditions: zero traction on a free surface; on a fluid–solid interface, continuous ``u_z`` and ``\sigma_{zz}`` and ``\sigma_{xz} = 0`` (a fluid holds no shear and may slip); on a welded interface, all four components continuous.
3. Each condition is a linear equation in the partial-wave amplitudes. Together they form a square, real matrix ``\mathbf M(c, \omega)`` (`boundary_matrix`), and a non-zero wave exists only where

```math
D(c, \omega) = \det \mathbf M(c, \omega) = 0 .
```

4. `find_speeds` scans ``c`` below the slowest body wave that must stay evanescent, looking for sign changes of ``D``, and refines each one by bisection. At a root, the null vector of ``\mathbf M`` gives the amplitudes (`mode_amplitudes`), and they give the displacement anywhere (`field_at`).

| setting | partial waves | conditions | matrix |
|:--|:--|:--|:--:|
| half-space | P and SV below | ``\sigma_{xz} = \sigma_{zz} = 0`` at ``z = 0`` | 2 × 2 |
| layer over half-space | P and SV in the layer (two each), P and SV below | free surface, then all four continuous at ``z = H`` | 6 × 6 |
| Scholte | P in the fluid, P and SV in the solid | ``u_z``, ``\sigma_{zz}`` continuous, ``\sigma_{xz} = 0`` | 3 × 3 |
| Stoneley | P and SV on each side | all four continuous | 4 × 4 |

The *Secular function* panel plots ``D`` divided by the product of the lengths of the columns of ``\mathbf M``. By Hadamard's inequality that ratio lies between −1 and 1 for any media (`normalized_secular`), and dividing by positive numbers moves no root.
"""

# ╔═╡ b9ddbe66-e0cc-43ec-965f-b6ae6f3b369f
md"""
## Rayleigh Waves on a Half-Space
With the two decaying waves of the half-space, the free-surface matrix is 2 × 2 and its determinant is proportional to the classical **Rayleigh function** (`rayleigh_function`)

```math
R(c) = \Big(2 - \frac{c^2}{\beta^2}\Big)^2 - 4\sqrt{1 - \frac{c^2}{\alpha^2}}\,\sqrt{1 - \frac{c^2}{\beta^2}} .
```

Apart from the trivial ``c = 0`` it has a single root below ``\beta``: the Rayleigh speed ``c_R``. ``R`` contains no length and no frequency, so ``c_R`` is the same at every frequency. On a homogeneous half-space Rayleigh waves are **not dispersive**, and every length in the wave scales with the wavelength ``\lambda = c_R/f``; that is why the experiment measures depth in wavelengths.

For a Poisson solid (``\lambda = \mu``, so ``\alpha = \sqrt3\,\beta``) the root is ``c_R = \beta\sqrt{2 - 2/\sqrt3} = 0.9194\,\beta``. The amplitudes at that root fix the shape of the motion: at the surface the vertical amplitude is 1.47 times the horizontal one, and the horizontal motion vanishes and changes sign at a depth of ``0.19\lambda`` (`ux_node`). Above that depth the orbit is **retrograde**: at the top of the ellipse the particle moves against the direction of travel. Below it the orbit is prograde. Both amplitudes decay with depth, at the two rates ``\omega\eta_\alpha`` and ``\omega\eta_\beta``.

!!! tip "Why P and SV must be coupled"
	The P and SV waves share the horizontal slowness ``p``, and neither alone can satisfy both free-surface conditions. Switch either part off in the wavefield panel to see how different each one is from their sum.
"""

# ╔═╡ 98e4a69f-ebd5-4a79-b494-1848b456acf7
md"""
## Why Rayleigh Waves Follow β
``R(c)`` depends only on ``c/\beta`` and ``\alpha/\beta``, so

```math
c_R = \beta\, g(\alpha/\beta) .
```

The function ``g`` barely changes: from ``0.874`` at Poisson ratio ``\nu = 0`` (``\alpha/\beta = \sqrt2``) to ``0.955`` as ``\nu \to 0.5`` (incompressible). The Rayleigh speed is essentially a fixed fraction of the shear speed, and the compressional speed only nudges that fraction. In numbers (`log_sensitivities`), a 1 % change of ``\beta`` changes ``c_R`` of a Poisson solid by 0.87 %, a 1 % change of ``\alpha`` by 0.13 %, and density drops out altogether.

The sensitivities obey a sum rule. Scaling every velocity by a factor ``s`` at fixed frequency is the same as scaling ``c`` by ``s`` and the frequency by ``1/s``, which gives

```math
\sum_{\text{all } \alpha,\, \beta} \frac{\partial \ln c}{\partial \ln v} = \frac{c}{U} , \qquad \sum_{\text{all } \rho} \frac{\partial \ln c}{\partial \ln \rho} = 0 ,
```

where ``U = d\omega/dk`` is the group velocity. For a non-dispersive wave ``U = c``, and the velocity sensitivities share exactly 1 between them; the bars show how they share it. This is why surface-wave dispersion is inverted for shear velocity: it is what the data constrain.
"""

# ╔═╡ 3ff5f5f1-5198-4a1f-8788-a6007ad1f494
md"""
## A Layer Over a Half-Space
A layer of thickness ``H`` brings a length into the problem, and the wave becomes **dispersive**: its speed depends on ``\omega H/\beta_1``, and the 6 × 6 boundary matrix has roots ``c_n(\omega)``, one for each **mode** ``R_n``.

- At low frequency the wavelength is much longer than ``H``; the wave hardly notices the layer and travels near the Rayleigh speed of the half-space.
- At high frequency the wave is confined to the top of the layer and travels at the layer's Rayleigh speed.
- The higher modes ``R_1, R_2, \dots`` exist only above their **cut-off frequencies**, where they appear at ``c = \beta_2``. A guided mode needs ``c < \beta_2`` so that the half-space waves decay; a slower layer (``\beta_1 < \beta_2``) is what traps the energy.

`dispersion_curves` repeats the root search over frequency, and `group_velocity` differentiates ``\omega`` with respect to ``k = \omega/c`` along each curve (dotted in the panel). The ghost curves show the same modes with every ``\alpha`` or every ``\beta`` raised by 10 %, and the β ghost moves much further at every frequency. The sensitivity bars split the effect between the layer and the half-space: at high frequency the wave feels mostly ``\beta_1``, at low frequency mostly ``\beta_2``. That frequency-dependent split is what lets a dispersion curve be inverted for a shear-velocity profile.
"""

# ╔═╡ 80c3c6de-bb53-45c6-af71-22537145d715
md"""
## The H/V Ratio
The ellipticity of the orbit at the surface, ``H/V = |u_x/u_z|`` at ``z = 0`` (`surface_hv`), is a property of the medium. On a half-space it is a constant (``1/1.47 = 0.68`` for a Poisson solid). Over a layer it depends on frequency, and for a soft layer on a stiff half-space it has a sharp peak close to the layer's quarter-wavelength resonance

```math
f_0 = \frac{\beta_1}{4H} .
```

There the vertical motion of the fundamental mode nearly vanishes, and the orbit turns from retrograde to prograde. Try the *Soft sediment over bedrock* preset: with ``\beta_1 = 200`` m/s and ``H = 30`` m the peak sits near 1.7 Hz. Change ``H`` or ``\beta_1`` and the peak follows ``f_0``; with a weak contrast (*Weak contrast*) it flattens out.

This is the physics behind the **H/V spectral ratio method** (Nakamura, 1989). The ratio of the horizontal to the vertical spectrum of ambient noise, recorded at a single three-component station, peaks near the resonance of the sediment. That one passive measurement gives ``f_0``, and with a known ``\beta_1`` the sediment thickness ``H = \beta_1/4f_0``; it is widely used to map site amplification for seismic hazard. Real noise also contains Love waves and body waves, so a measured H/V curve only approximates the Rayleigh ellipticity drawn here.
"""

# ╔═╡ 176fcfaa-c512-4d66-b3b0-db45dc7f2922
md"""
## Interface Waves: Scholte and Stoneley
An interface between two half-spaces can also guide a wave that travels along it and decays into both sides.

**Scholte waves** live on a fluid–solid interface such as the seafloor. The fluid carries only a P wave and no shear, so its side contributes one partial wave, and there are three conditions: ``u_z`` and ``\sigma_{zz}`` continuous and ``\sigma_{xz} = 0`` in the solid, while ``u_x`` may slip. A Scholte wave always exists, slower than both the speed of sound in the fluid and the shear speed of the solid.

- On soft sediment (``\beta_s`` below the speed of sound in water) it travels a little below ``\beta_s``, is carried mainly by the solid, and its sensitivity is almost all ``\beta_s``. Ocean-bottom surveys use these waves to map shallow shear velocity.
- On hard rock it travels just below the speed of sound in water, decays slowly into the water, and its sensitivity is almost all ``\alpha_f``: it is mostly a wave in the fluid.

As the density of the fluid goes to zero, the Scholte wave becomes the Rayleigh wave of the solid.

**Stoneley waves** live on a welded interface between two solids, with two decaying waves on each side and four continuity conditions. Unlike the other cases, a Stoneley wave exists only for a narrow range of media, roughly when the two shear speeds are close, and then it travels just below the slower shear speed. With equal shear speeds and a fourfold density contrast it travels at ``0.975\beta``; for crust against mantle there is no root, and the experiment says so.

!!! note "Two different Stoneley waves"
	The "Stoneley wave" (or tube wave) of borehole logging is a guided wave in a fluid-filled hole. It is a relative of the Scholte wave, not of the solid–solid interface wave here.
"""

# ╔═╡ 379feb55-32e5-4867-b108-b10c2c0782c3
md"""
## Appendix
The solver, in the order of the recipe: the media and the four problems, the partial waves and boundary matrices, roots and null vectors, then the fields, dispersion and sensitivities the panels draw, the checks, and the widget. It uses no package beyond PlutoUI.
"""

# ╔═╡ a3db0f69-7eed-4f84-9e85-5869ba1d50f9
md"""
### Media and Problems
A `Medium` holds ``(\alpha, \beta, \rho)``; a fluid has ``\beta = 0``. Each setting is its own type, so every later function picks the right partial waves and conditions by dispatch.
"""

# ╔═╡ 159610bd-4a45-4293-96fc-230029ec55dd
md"""
### Partial Waves and Boundary Matrices
`p_state` and `s_state` are the displacement–stress vectors of *Displacement and Traction*. `partial_states` evaluates every partial wave of a problem at one depth, and `boundary_matrix` stacks the boundary conditions of *One Recipe for Every Boundary*. Inside the layer, `layer_basis` writes an evanescent wave as one function decaying from the surface and one growing towards the bottom, both bounded, rather than as cosh and sinh, which become indistinguishable at high frequency.
"""

# ╔═╡ eb87acec-3e9e-4069-9d17-464f7c4590e5
md"""
### Roots and Null Vectors
The matrices are at most 6 × 6, so plain Gaussian elimination is enough: partial pivoting for the determinant (`small_det!`), complete pivoting for the null vector (`null_vector`), whose vanishing pivot it moves to the last row.
"""

# ╔═╡ 31cbf1a6-722f-4edd-bf6c-67437800cbfc
md"""
### Fields, Dispersion and Sensitivity
`field_at` sums the partial waves into ``U``, ``W`` and the tractions at any depth. From it come the surface H/V, the depth profiles and the ``u_x`` node, and the wavefield frames. For the layer, `dispersion_curves` and `group_velocity` give the modes' speeds; `log_sensitivities` perturbs each parameter of each medium by ±0.1 % and follows the nearest root.
"""

# ╔═╡ 7eb2dd96-6853-4ff5-90da-c1d5f7c82c70
md"""
### What the Widget Is Sent
`rayleigh_payload` collects, for the current state, everything the panels draw and `to_json` writes it out for the browser. The dispersion curves of the layer and the sensitivities depend on fewer inputs than the bond carries, so `memo!` keeps the last result of each.
"""

# ╔═╡ ae264dbd-042e-4e81-9ed9-bda053bcbd50
md"""
### Verifying the Solver
The textbook numbers for a Poisson half-space, the closed-form Rayleigh function, tractions recomputed from the displacement by finite differences, the limits of the layer problem, the Scholte-to-Rayleigh limit, continuity of the Stoneley wave, and the sensitivity sum rules.
"""

# ╔═╡ aa111199-0003-4003-8003-100000000003
md"""
### The Interactive Widget
`RayleighFunctionInput` sends the setting, both media, the layer thickness, the frequency (as ``\log_{10} f/f_0``), the orbit depth (in wavelengths), the orbit phase and the chosen mode. The cell under the widget converts them to a `WaveProblem` and pushes `rayleigh_payload` to the panels. Playing the wavefield animation happens in the browser.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoTeachingTools = "~0.3.1"
PlutoUI = "~0.7.72"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.1"
manifest_format = "2.1"
project_hash = "3cea405850043ca26da5fca44cf9b16b7c2afdba"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "980f01d6d3283b3dbdfd7ed89405f96b7256ad57"
registries = "General"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "2.0.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
registries = "General"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

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
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
registries = "General"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "277779adfedf4a30d66b64edc75dc6bb6d52a16e"
registries = "General"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.6"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
registries = "General"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "e24491cb83551e44a69b9106c50666dea9d953ab"
registries = "General"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.4.4"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
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
version = "2026.8.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
registries = "General"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
registries = "General"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoLinks", "PlutoUI"]
git-tree-sha1 = "8252b5de1f81dc103eb0293523ddf917695adea1"
registries = "General"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.3.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "f53232a27a8c1c836d3998ae1e17d898d4df2a46"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "276237d98b5bfc2541d92e952bfa99d362011ac7"
registries = "General"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.11.0"

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

    [deps.Revise.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
registries = "General"
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
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
registries = "General"
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

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.2+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╠═d921cf1c-ca8b-44d9-8d00-e1ed55267647
# ╠═915243cb-9619-4eab-8ad4-1d9b821e8d01
# ╟─8fa36b02-b810-4cbd-ac0a-75e0596d73bb
# ╟─07ee82c2-d17b-4b85-9ef6-c107b76d0cfe
# ╟─544d41b1-b6ce-4844-8269-a0251b391962
# ╟─31a317e2-0864-4552-83ee-ba4138ab5c0b
# ╟─318ac5d5-bd8b-417b-a919-b46ee555ca1a
# ╟─bd7562f0-cb14-4ddb-b8c7-ffc2a37b079e
# ╟─e89ffc5a-8252-4186-8600-dd68f91b425d
# ╟─6ada314c-d48f-4d01-8a53-afc06ce649cb
# ╟─1b2b6012-d997-4f43-8b58-87e8bc40acf3
# ╟─b9ddbe66-e0cc-43ec-965f-b6ae6f3b369f
# ╟─98e4a69f-ebd5-4a79-b494-1848b456acf7
# ╟─3ff5f5f1-5198-4a1f-8788-a6007ad1f494
# ╟─80c3c6de-bb53-45c6-af71-22537145d715
# ╟─176fcfaa-c512-4d66-b3b0-db45dc7f2922
# ╟─379feb55-32e5-4867-b108-b10c2c0782c3
# ╠═5b807994-416d-11ed-22ff-77b37ff3cfac
# ╟─a3db0f69-7eed-4f84-9e85-5869ba1d50f9
# ╠═66e5c7d4-a32a-4a5a-9053-dab1656cfd60
# ╠═5a8f7ef7-64c6-4655-b676-a8624ba28851
# ╟─159610bd-4a45-4293-96fc-230029ec55dd
# ╠═c167d1e3-6fb0-494f-b261-2dc99b29eb7a
# ╠═7342a27d-2b70-43b9-a892-3eeae6db8275
# ╠═6dff8f48-78c1-400e-bc51-05c8a8bf4f89
# ╟─eb87acec-3e9e-4069-9d17-464f7c4590e5
# ╠═3602a131-0ffd-4521-b977-bf6b4c1a6281
# ╠═36599b9f-505b-4a9c-b6e2-464997aad540
# ╠═bb111111-1111-4111-8111-111111111111
# ╟─31cbf1a6-722f-4edd-bf6c-67437800cbfc
# ╠═23cbc525-246b-4808-af37-2d3364857ccc
# ╠═4a75f3f9-e7c6-4692-97b4-cf09e2d519fa
# ╠═937eb11f-d10e-4ddb-9511-32e3d7bb69b6
# ╠═45910096-f491-4381-aa54-4e85c1f1345b
# ╟─7eb2dd96-6853-4ff5-90da-c1d5f7c82c70
# ╠═826c0d0b-2957-481c-baf4-aec6719e75a6
# ╠═1aac4c99-af17-4b88-97c8-d78fcb1832d5
# ╠═7beb097d-9d0f-44fa-9161-5c65ad3f676c
# ╠═bc2008fa-ad4b-40ff-8bb3-9963d31246ed
# ╠═00b0672a-ff69-46aa-aa07-2a69052a14c6
# ╟─ae264dbd-042e-4e81-9ed9-bda053bcbd50
# ╠═c760b9b7-1896-4657-aa30-7d8ddcb8aa40
# ╟─aa111199-0003-4003-8003-100000000003
# ╠═a48f6938-2667-4401-9721-ee53a692c6ef
# ╠═bb222222-2222-4222-8222-222222222222
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
