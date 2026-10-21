### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> tags = ["surfacewaves"]
#> title = "Love Waves"
#> description = "Dive into the world of plane Love waves in a single homogeneous layer overlying a homogeneous half-space, offering a simple yet fascinating exploration of the topic."
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

# ╔═╡ 66ae626f-2543-4897-946f-8e51ecbb95e0
using PlutoUI, Roots, Symbolics, SymbolicUtils, FFTW, PlutoTeachingTools, Latexify, Peaks

# ╔═╡ 9b02e72d-be6b-485d-a186-ba2fe2dcd6fb
ChooseDisplayMode()

# ╔═╡ bdf3d004-653a-4b56-af7c-42ca1d6021e5
TableOfContents()

# ╔═╡ 58ac96df-06f8-43b1-b19c-4f16c86504e8
# App shell: the Expt / Theory / Code / Help frame shared by every notebook, loaded
# from src/assets/app-shell/. Without that folder (a lone downloaded .jl) this cell just
# shows the help text and the notebook stays a plain notebook.
let
	help = md"""
	- **Drag** the layer boundary on the wavefield to change the layer thickness; drag the frequency cursor on the dispersion panel to change the frequency.
	- The **β** and **ρ** sliders set the layer and the half-space; **click** a dispersion curve to add or remove that mode from the animation.
	- **Brush** a frequency band to synthesize the seismogram; drag its **time cursor** to see the stationary-phase frequency at that time.
	- The receiver slider moves the station along the surface.

	**Try this:** raise the frequency and watch higher modes appear one by one as they pass their cut-off frequencies.
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

# ╔═╡ 2a41b15e-1a0e-4c92-a3a5-53603faacea1
md"""
# Love Waves
This notebook studies the simplest case of plane Love waves: a single homogeneous layer over a homogeneous half-space. It is an eigenvalue problem: the free-surface and interface conditions allow only certain phase velocities and mode shapes at each frequency.

With the experiment you can

- visualize the displacement wavefield of the different modes;
- see the cut-off frequency of the ``n``th higher mode;
- follow how a dispersed wave packet builds up at a receiver.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India

"""

# ╔═╡ 4479381e-b424-4599-b692-877e4a15405e
md"""
## The Experiment
Adjust the medium (drag the layer boundary, or use the ``\beta`` and ``\rho`` sliders) and the operating frequency. The dispersion, wavefield, seismogram and stationary-phase panels react together. Click a dispersion curve to include or exclude that mode from the animated wavefield.
"""

# ╔═╡ 9f106bb0-b0c6-4c5b-bdf2-f17c15693b82
Markdown.MD(Markdown.Admonition("warning", "Observations", [md"""
* Use the frequency band to synthesize the seismogram of a wave packet with a finite band around a central frequency ``\omega_0``. Do you see an oscillation at ``\omega_0`` modulated by an envelope?
* Drag the time cursor on the seismogram: the stationary-phase panel marks the frequency whose phase is stationary at that time, which is the frequency that dominates the record there.
"""
]))

# ╔═╡ a19519d6-a96f-4a8a-8563-d5052f694aff
md"---"

# ╔═╡ 8d2fd2d9-da36-43f4-8dcb-dd931c8189f5
@syms β₁::Real β₂::Real μ₁::Real μ₂::Real ρ₁::Real ρ₂::Real

# ╔═╡ 25f32a28-4c57-4ed0-9fc4-4a40ca31e4dd
@syms H::Real

# ╔═╡ d8b52234-757b-441d-918d-585fd6593fcf
md"""
## Trial Wavefield

Love waves are SH motion: only the ``y`` component of displacement, ``u``, is nonzero. In each layer we try a pair of plane waves sharing the horizontal slowness ``p``,

```math
u = \big(A\,e^{-\omega\eta z} + B\,e^{\omega\eta z}\big)\,e^{i\omega(t - p x)},
```

implemented as `trail_soln(p, η, A, B)`. When ``p < 1/\beta`` the wave propagates vertically and ``\eta = i\sqrt{1/\beta^2 - p^2}`` is imaginary, so the two terms are down- and up-going waves. When ``p > 1/\beta``, ``\eta = \sqrt{p^2 - 1/\beta^2}`` is real and the terms grow or decay exponentially with depth: inhomogeneous waves travelling along ``x``. `getη` makes exactly this choice.
"""

# ╔═╡ b120da8e-dcd6-4280-819b-0e20b4675080
@syms p::Real η₂::Complex{Real}

# ╔═╡ c74c4a58-b47e-4509-aadd-7400c350ff6c
@syms ı::Complex{Real} # imaginary unit, going to substitute with im later

# ╔═╡ 44aa241c-bf9d-4a2a-b459-28a73a991606
@syms x::Real z::Real t::Real ω::Real

# ╔═╡ 5896fa94-ab7f-42f0-9fac-265e11faba2b
trail_soln(p, η, A, B) = (A * exp(-ω * η * z) + B * exp(ω * η * z)) * exp(ı * ω * (t - (p * x)))

# ╔═╡ 6e734951-2907-49d7-ae88-dffc463bf8a3
md"""
In the layer ``0 \le z \le H`` the field must satisfy the SH wave equation ``\partial_t^2 u = (\mu_1/\rho_1)\,\nabla^2 u``; its operator is `L1`.
"""

# ╔═╡ 2f40cb44-60af-4fbe-8f73-d8d06e99fc3b
begin
    Dz = Differential(z)
    Dx = Differential(x)
    Dt = Differential(t)
end

# ╔═╡ 4d859323-12e0-4fba-a796-48995f980f33
L1(u) = Dt(Dt(u)) - μ₁ / ρ₁ * (Dx(Dx(u)) + Dz(Dz(u)))

# ╔═╡ 79c4f3a4-65b6-4e02-9786-f168584c7781
md"""
In the layer the field `u1` has vertical slowness ``\eta_1`` and amplitudes ``A_1``, ``B_1``.
"""

# ╔═╡ dc709a17-ab03-406f-b232-d7658872a95e
@syms A₁::Real B₁::Real η₁::Complex{Real}

# ╔═╡ 83e2b2b3-f84a-4926-9332-b5570d52be8b
u1 = trail_soln(p, η₁, A₁, B₁)

# ╔═╡ d6dc2422-6be3-447e-b7c1-be655971b6f7
md"""
In the half-space ``z > H`` the field `u2` has vertical slowness ``\eta_2`` and amplitudes ``A_2``, ``B_2``.
"""

# ╔═╡ 1a13363a-9f8f-4123-9b37-c6dcad1ccdde
@syms A₂::Real B₂::Real

# ╔═╡ 2f80cfcf-9788-41a3-8797-cd1965216738
u2 = trail_soln(p, η₂, A₂, B₂)

# ╔═╡ de304d8c-253e-4cb5-8332-db79de20991c
md"""
A trapped Love wave has ``\beta_1 < c < \beta_2``, so ``p > 1/\beta_2`` and ``\eta_2`` is real. The term ``B_2 e^{\omega\eta_2 z}`` would grow without bound with depth, so we set ``B_2 = 0``: the half-space field decays away from the interface.
"""

# ╔═╡ 30219c41-a040-4056-adda-ce64ac3c6969
B₂ ~ 0

# ╔═╡ e1ab342d-84fe-4c00-9756-455436b0fad7
md"""
## Free Surface Condition
At ``z = 0`` the traction vanishes: ``\sigma_{yz} = \mu_1\,\partial u_1/\partial z = 0`` (`σyz1`). This forces ``A_1 = B_1``, so the layer field is a standing wave in depth, ``\cosh(\omega\eta_1 z)`` (a cosine when ``\eta_1`` is imaginary).
"""

# ╔═╡ ca377428-d597-47a2-b29e-f40884697bd6
σyz1 = expand_derivatives(μ₁ .* Dz(u1))

# ╔═╡ 22ef3c54-fcb8-463e-8520-4d07bbceda51
substitute(σyz1, z => 0) ~ 0

# ╔═╡ b819ce05-628a-49b6-aa57-f1ee765fb027
Symbolics.symbolic_linear_solve(substitute(σyz1, z => 0), B₁) |> simplify

# ╔═╡ 53799d32-3edc-4416-a928-868ea974c6bd
A₁ ~ B₁

# ╔═╡ f39b1e01-c16e-405e-a2fc-85408d53c762
md"""
## Dynamic Boundary Condition
Across the interface ``z = H`` the shear traction ``\sigma_{yz}`` is continuous (`σyz1H` equals `σyz2H`). With ``B_1 = A_1`` and ``B_2 = 0`` this gives one expression for ``A_2`` in terms of ``A_1`` (`A₂ex1`).
"""

# ╔═╡ b3924a35-6565-46c4-9a6e-d29b1172668e
σyz2 = expand_derivatives(μ₂ .* Dz(u2))

# ╔═╡ 30f26e47-0a37-4165-b18b-a20f1148deb1
σyz1H = substitute(σyz1, z => H)

# ╔═╡ f7b5adba-09ca-4895-8b18-5fc84b752392
σyz2H = substitute(σyz2, z => H)

# ╔═╡ 8f35d8b0-ff20-442c-b5b2-f5a2616fe345
A₂ex1 = substitute(Symbolics.symbolic_linear_solve(σyz1H ~ σyz2H, A₂), [B₁ => A₁, B₂ => 0])

# ╔═╡ dab840e2-7d7c-4f65-b4ea-289babc430c8
md"""
## Kinematic Boundary Condition
The displacement is also continuous at ``z = H`` (`u1H` equals `u2H`), which gives a second expression for ``A_2`` (`A₂ex2`).
"""

# ╔═╡ 644637c9-c32d-4335-af59-3452263a9707
u1H = substitute(u1, [z => H, B₁ => A₁])

# ╔═╡ 09a8041d-6443-4e83-8443-e0e2aa05317a
u2H = substitute(u2, [z => H, B₂ => 0])

# ╔═╡ 0fc4d1a1-25ac-49c9-b74b-a9f26928c652
A₂ex2 = Symbolics.symbolic_linear_solve(u1H ~ u2H, A₂)

# ╔═╡ 794a1f79-708d-4305-8aa0-d2673f34e31f
md"""
## Graphical Solution
A nonzero wave needs the two expressions for ``A_2`` to agree: the condition `ex1 = ex2` is the Love-wave dispersion relation. For a trapped mode (``\beta_1 < c < \beta_2``) it is the classical

```math
\tan\!\Big(\omega H\sqrt{\tfrac{1}{\beta_1^2} - p^2}\Big) = \frac{\mu_2\sqrt{p^2 - 1/\beta_2^2}}{\mu_1\sqrt{1/\beta_1^2 - p^2}} ,
```

whose roots ``p`` at each frequency are the allowed modes. Because the tangent is periodic, the number of roots grows with ``\omega H``: each new branch is a higher mode, appearing above its cut-off frequency. The cells below simplify the two expressions and set up the phase-velocity ranges for the medium chosen in the experiment.
"""

# ╔═╡ bfb54c6f-a52c-4a10-87c7-9788f2f63629
ex1 = simplify(A₂ex1 * arguments(A₂ex1)[2] / A₁)

# ╔═╡ cbbae944-5377-44a4-ae32-868a75625248
ex2 = simplify(A₂ex2 * arguments(A₂ex1)[2] / A₁)

# ╔═╡ b5c5f1d0-0075-4727-a6a9-8b20d6f65bc4
md"Note that the goal here is to find zeros of the function `F` defined below,  which outputs `ex1-ex2` for a given horizontal slowness `p`."

# ╔═╡ ad269679-e1f9-4d9d-827f-468e36f27bfc
md"""
`F`, `U1` and `U2` below are built **once**, directly from the symbolic derivation above,
with `μ₁`, `μ₂`, `H` left as *symbolic* (not substituted) parameters -- so `build_function`
compiles them a single time at notebook load, and every widget interaction below is just a
numeric call into an already-compiled function. No `substitute`/`build_function` call ever
runs again as the medium or frequency change.
"""

# ╔═╡ b2858a78-4498-4374-9ee0-c68f3dfca6e8
md"From the above plots, we can observe that there are no zeros in either (0,β₁) or (β₂, ∞). Therefore, we will now find the zeros of `F` in the interval [β₁, β₂]. The `find_zeros` function from the [`Roots.jl`](https://github.com/JuliaMath/Roots.jl) package can be used to search for all zeros in a specified interval in a derivative-free fashion. We only consider the real part of `F` as the imaginary part is zero in [β₁, β₂]. The zeros correspond to the eigenvalue phase velocities."

# ╔═╡ 8db66193-eb99-478e-bdad-44b0e4e18dad
md"""
## Dispersion Curves
`compute_dispersion_curves` repeats the root search at every frequency, giving each mode's phase velocity ``c_n(f)``. `phase_velocities` arranges them by mode (up to four), `get_k` converts them to wavenumbers ``k = \omega/c``, and `get_group_velocity` takes ``U = d\omega/dk`` by finite differences. The group velocity is what governs the arrival of energy at the receiver.
"""

# ╔═╡ 7c1e5a20-3b8d-4f61-9a0e-2d6b4c8e1f37
begin
	"""
		LastCall()

	A one-entry memo: `memo!(f, c, key)` returns `f()`, calling `f` again only when `key`
	differs from the key of the previous call.

	Every widget interaction re-runs every cell downstream of `lw`, even when only the time
	cursor moved. The dispersion curves depend on the medium alone, so the ~500 root searches
	of [`compute_dispersion_curves`](@ref) sit behind one of these.
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

	_lw_dispersion_memo = LastCall()
end


# ╔═╡ e4c1cda0-7366-4aa0-a616-47d86a1c9502
"""
	get_k(c, freqgrid)

Wavenumbers ``k = ω/c`` of one mode's phase-velocity curve `c`, sampled at `freqgrid[2:end]`.
"""
get_k(c, freqgrid) = inv.(c) .* @view(freqgrid[2:end]) .* 2π

# ╔═╡ 6ae71567-4501-4cf5-8055-47b235829f14
"""
	get_group_velocity(k, freqgrid)

Group velocity ``U = dω/dk`` of one mode by finite differences of its wavenumbers `k`
over the uniform `freqgrid`.
"""
get_group_velocity(k, freqgrid) = 2π * step(freqgrid) ./ diff(k)

# ╔═╡ df35e923-b43c-4809-a9db-32a371fc9010
md"""
## Appendix
Grids, the slowness helper `getη`, the data pipeline that feeds the dashboard, and the compiled functions `F`, `U1`, `U2` built from the symbolic derivation above.
"""

# ╔═╡ 896fb716-861a-44d1-b073-45c47773a4f8
tgrid = range(0, stop=1000, length=1000)

# ╔═╡ 4ee3be7e-7fa9-4f4a-b13a-3d20a8863f04
freqgrid = rfftfreq(length(tgrid), inv(step(tgrid)))

# ╔═╡ 840a743a-e298-4f99-ae6e-11c85b6f5bc5
# a range for phase velocity (typical shear velocity values in crust/ mantle)
cgrid = range(1, stop=7, step=0.2)

# ╔═╡ 80bb6ae1-c48b-4be1-9786-c94c6e8680a1
"""
	getη(p, β)

Vertical slowness for horizontal slowness `p` in a medium of shear velocity `β`:
``i\\sqrt{1/β^2 - p^2}`` (propagating) when ``p ≤ 1/β``, else ``\\sqrt{p^2 - 1/β^2}``
(evanescent). Both branches return a complex number, so every caller stays type-stable.
"""
function getη(p, β)
	s = inv(abs2(β)) - abs2(p)
	return s >= 0 ? complex(zero(s), sqrt(s)) : complex(sqrt(-s))
end

# ╔═╡ 039a52c9-01b8-4028-9eae-1ff5b779fd78
"""
	compute_phase_velocities(F_func, freq, medium_params)

Root-find `F_func(p, 2πfreq, η₁, η₂, μ₁, μ₂, H) = 0` over `p ∈ [1/β₁, 1/β₂]` at a
single frequency, returning the eigenvalue phase velocities `1/p` sorted ascending.
`F_func` is the generalized (medium-agnostic) `F` built once in the Appendix --
`μ₁`, `μ₂`, `H` are computed here from `medium_params` and passed in numerically,
so no symbolic substitution happens on this (reactive) path.
"""
function compute_phase_velocities(F_func, freq, medium_params)
    μ1 = medium_params.β₁^2 * medium_params.ρ₁
    μ2 = medium_params.β₂^2 * medium_params.ρ₂
    Fmedium = p -> F_func(p, 2 * pi * freq, getη(p, medium_params.β₁), getη(p, medium_params.β₂), μ1, μ2, medium_params.Hp)
    return sort(inv.(find_zeros(real ∘ Fmedium, inv(medium_params.β₁), inv(medium_params.β₂))))
end

# ╔═╡ a69a71b4-6dbc-4de7-a8d7-5cc3ad744d94
"""
	compute_dispersion_curves(F_func, freqgrid, medium_params)

[`compute_phase_velocities`](@ref) repeated over every frequency in `freqgrid[2:end]`,
giving each mode's phase-velocity dispersion curve.
"""
compute_dispersion_curves(F_func, freqgrid, medium_params) =
	[compute_phase_velocities(F_func, f, medium_params) for f in @view(freqgrid[2:end])]

# ╔═╡ 7bd19e44-b4f5-43f3-b6be-b557ca95c67f
md"""
To draw the wavefield we need the displacement in the layer and in the half-space as functions of ``x``, ``z``, ``t`` and ``p``. These are `U1` and `U2`, compiled at the end of the Appendix.
"""

# ╔═╡ 16c826bf-fa48-423c-bfec-195f7fc502f8
md"""
## TODO
- plot phase velocities as a function of period, too
"""

# ╔═╡ 28fbc1e8-2323-4d0d-8f43-e6c1041c7f8d
md"""
### Verifying the Root-Search Interval
`F` has no real zeros outside `[β₁, β₂]` -- confirming it is safe to search only that
interval for the eigenvalue phase velocities below.
"""

# ╔═╡ 361a3181-857a-4454-9fe6-64c8ae3be9cf
md"""
### Dashboard Data Pipeline
`medium` collects the layer thickness and the layer and half-space properties from the experiment, and `cn` holds the phase velocities of all modes at the operating frequency. The seismogram `record` at the receiver is a superposition over the selected frequency band of the fundamental mode, evaluated with `U1` at the surface. The stationary-phase points are the extrema over frequency of the phase ``k(\omega)x - \omega t`` at the cursor time (`phase_arr`, found with Peaks.jl). The wavefield frames sum the selected modes over one period with `U1` and `U2` (`love_wavefield_frames`); `love_seismogram` builds the record. Both evaluate the field once at ``t = 0`` and rotate it by ``e^{i\omega t}``, since time enters the trial solution only through that factor. Every interaction re-runs these cells, so the dispersion curves are memoized on the medium (`memo!`). Everything is packed into a `LovePush` for the widget.
"""

# ╔═╡ 749ed9e6-a16e-49ee-9f70-909aa0493b0d
begin
	const NX_LW = 70
	const NZ_LW = 35
	const NFRAMES_LW = 16
	const ZMAX_LW = 150.0
end

# ╔═╡ 3f9b2d71-6a4e-4c18-b5d3-8e1f0a7c2b94
"""
	love_seismogram(U1, freqs, cs, tgrid, x, β₁, fmin, fmax)

Fundamental-mode seismogram at surface distance `x`: the sum of ``\\mathrm{Re}\\,U_1(x, 0, t)``
over the frequencies of `freqs` inside `[fmin, fmax]`, each with its phase velocity from `cs`.
Time enters `U1` only through ``e^{iωt}``, so each frequency's complex amplitude is evaluated
once at ``t = 0`` and then rotated through `tgrid`.
"""
function love_seismogram(U1, freqs, cs, tgrid, x, β₁, fmin, fmax)
	rec = zeros(length(tgrid))
	for (f, c) in zip(freqs, cs)
		(fmin <= f <= fmax && !ismissing(c)) || continue
		ω, p = 2π * f, inv(c)
		a = complex(U1(x, 0.0, 0.0, p, ω, getη(p, β₁)))
		for (j, t) in enumerate(tgrid)
			rec[j] += real(a * cis(ω * t))
		end
	end
	return rec
end

# ╔═╡ a5d0c3e8-91f2-4b7a-8c6d-0e4f2b9a1d53
"""
	love_wavefield_frames(U1, U2, cs, freq, medium, xgrid, zgrid1, zgrid2, nframes)

Displacement ``\\mathrm{Re}\\,u_y(x, z, t)`` of the modes with phase velocities `cs`, all at
frequency `freq`, on the layer grid `zgrid1 × xgrid` (with `U1`) and the half-space grid
`zgrid2 × xgrid` (with `U2`), at `nframes` instants spanning one period.

Returns `(U1_frames, U2_frames, wmax)`, with `wmax` the largest absolute displacement over all
frames. Every mode shares the time factor ``e^{iωt}``, so the complex field is summed once at
``t = 0`` and each frame is its real part after a phase rotation.
"""
function love_wavefield_frames(U1, U2, cs, freq, medium, xgrid, zgrid1, zgrid2, nframes)
	ω = 2π * freq
	μ1 = medium.β₁^2 * medium.ρ₁
	μ2 = medium.β₂^2 * medium.ρ₂
	W1 = zeros(ComplexF64, length(zgrid1), length(xgrid))
	W2 = zeros(ComplexF64, length(zgrid2), length(xgrid))
	for c in cs
		p = inv(c)
		η1, η2 = getη(p, medium.β₁), getη(p, medium.β₂)
		for (ix, x) in enumerate(xgrid)
			for (iz, z) in enumerate(zgrid1)
				W1[iz, ix] += U1(x, z, 0.0, p, ω, η1)
			end
			for (iz, z) in enumerate(zgrid2)
				W2[iz, ix] += U2(x, z, 0.0, p, ω, η1, η2, μ1, μ2, medium.Hp)
			end
		end
	end
	frame_times = range(0, step=inv(freq) / nframes, length=nframes)
	U1_frames = [real.(W1 .* cis(ω * t)) for t in frame_times]
	U2_frames = [real.(W2 .* cis(ω * t)) for t in frame_times]
	wmax = max(maximum(m -> maximum(abs, m), U1_frames), maximum(m -> maximum(abs, m), U2_frames), 1.0e-9)
	return U1_frames, U2_frames, wmax
end

# ╔═╡ 648572e8-2cf4-4a1f-932d-02aa2160b6d9
begin
	jsnum(x) = ismissing(x) || x === nothing || !isfinite(x) ? "NaN" :
		x isa Integer ? string(x) : string(round(x; sigdigits=5))
	flatten_js(v) = join(jsnum.(v), ",")
	flatten_rowmajor_js(M) = join(jsnum.(vec(permutedims(M))), ",")
	flatten_frames_js(frames) = join((flatten_rowmajor_js(M) for M in frames), ",")
end

# ╔═╡ beb71023-cd94-4fd2-8108-c68ad8b58055
begin
	struct LovePush
	    freqgrid::Any
	    phase_velocities::Any
	    group_velocities::Any
	    cn::Any
	    freq::Float64
	    nx::Int
	    nz::Int
	    nframes::Int
	    Hp::Float64
	    zmax::Float64
	    xmax::Float64
	    wmax::Float64
	    U1_frames::Any
	    U2_frames::Any
	    tgrid::Any
	    record::Any
	    arrival1::Float64
	    arrival2::Float64
	    phaseArr::Any
	    minPoints::Any
	    maxPoints::Any
	end

	function Base.show(io::IO, ::MIME"text/html", p::LovePush)
	    write(io, """
	    <script>
	    {
	      const root = document.getElementById('lwwidget')
	      if (root) {
	        root.dispatchEvent(new CustomEvent('lw-data', { detail: {
	          freqgrid: [$(flatten_js(p.freqgrid))],
	          phase_velocities: [$(join(["[" * flatten_js(v) * "]" for v in p.phase_velocities], ","))],
	          group_velocities: [$(join(["[" * flatten_js(v) * "]" for v in p.group_velocities], ","))],
	          cn: [$(flatten_js(p.cn))],
	          freq: $(p.freq),
	          nx: $(p.nx), nz: $(p.nz), nframes: $(p.nframes),
	          Hp: $(p.Hp), zmax: $(p.zmax), xmax: $(p.xmax), wmax: $(p.wmax),
	          U1_frames: [$(flatten_frames_js(p.U1_frames))],
	          U2_frames: [$(flatten_frames_js(p.U2_frames))],
	          tgrid: [$(flatten_js(p.tgrid))],
	          record: [$(flatten_js(p.record))],
	          arrival1: $(p.arrival1), arrival2: $(p.arrival2),
	          phaseArr: [$(flatten_js(p.phaseArr))],
	          minPoints: [$(flatten_js(p.minPoints .- 1))],
	          maxPoints: [$(flatten_js(p.maxPoints .- 1))],
	        }}))
	      }
	    }
	    </script>
	    """)
	end
end

# ╔═╡ 42190b26-0b30-4b2d-bc1f-dbb3323d7bb4
md"""
### The Interactive Widget
`LoveWaveInput` sends the layer thickness, the four medium properties, the operating frequency, the receiver distance, the frequency band, the time cursor and the selected modes back to Julia. Each change re-evaluates the already-compiled functions on small grids; the symbolic derivation never re-runs.
"""

# ╔═╡ e80a6231-6f27-4fe4-a862-32976ba3ab92
begin
	"""
		LoveWaveInput(; Hp, β1, β2, ρ1, ρ2, freq, xrecpos, freqmin, freqmax, t, selectedModes)

	The single interactive dashboard: draggable medium boundary and frequency
	cursor, β/ρ and receiver sliders, plus the wavefield animation, dispersion
	curves with click-to-toggle mode selection, a frequency-band brush, a
	seismogram with a draggable time cursor, and a stationary-phase panel. None
	of these fields touch the symbolic pipeline -- Julia only evaluates
	already-built closed-form functions over small grids, so every interaction
	stays cheap.
	"""
	struct LoveWaveInput
	    Hp::Float64
	    β1::Float64
	    β2::Float64
	    ρ1::Float64
	    ρ2::Float64
	    freq::Float64
	    xrecpos::Float64
	    freqmin::Float64
	    freqmax::Float64
	    t::Float64
	    selectedModes::Vector{Int}
	end

	LoveWaveInput(; Hp=35.0, β1=3.5, β2=4.5, ρ1=2.6, ρ2=3.4, freq=0.08,
	    xrecpos=2000.0, freqmin=0.0, freqmax=0.4995, t=500.0, selectedModes=[1]) =
	    LoveWaveInput(Hp, β1, β2, ρ1, ρ2, freq, xrecpos, freqmin, freqmax, t, selectedModes)

	Base.get(w::LoveWaveInput) = Dict{String,Any}(
	    "Hp" => w.Hp, "β1" => w.β1, "β2" => w.β2, "ρ1" => w.ρ1, "ρ2" => w.ρ2, "freq" => w.freq,
	    "xrecpos" => w.xrecpos, "freqmin" => w.freqmin, "freqmax" => w.freqmax,
	    "t" => w.t, "selectedModes" => w.selectedModes,
	)

	function Base.show(io::IO, ::MIME"text/html", w::LoveWaveInput)
	    write(io, """
	    <div id="lwwidget" class="app-widget">
	    <style>
	    #lwwidget{font-family:sans-serif;color:#d1d5db}
	    #lwwidget .lw-titlebar{background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:12px 16px;margin-bottom:14px;text-align:center}
	    #lwwidget .lw-titlebar-headline{font-size:18px;font-weight:700;color:#f3f4f6}
	    #lwwidget .lw-titlebar-sub{font-size:13px;color:#9ca3af;margin-top:4px}
	    #lwwidget .lw-row{display:flex;gap:14px;flex-wrap:wrap;margin-bottom:14px}
	    #lwwidget .lw-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px;flex:1 1 560px}
	    #lwwidget .lw-panel-title{font-size:14px;font-weight:700;color:#f3f4f6;margin-bottom:4px}
	    #lwwidget .lw-caption{font-size:13px;color:#9ca3af;margin-top:4px}
	    #lwwidget canvas{display:block;cursor:crosshair;max-width:100%;height:auto}
	    #lwwidget .lw-mini-group{background:#0b0b0b;border:1px solid #1f2937;border-radius:6px;padding:8px 10px;margin-top:10px}
	    #lwwidget .lw-mini-title{font-size:13px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
	    #lwwidget .lw-controls{display:flex;gap:14px;flex-wrap:wrap;align-items:center}
	    #lwwidget .lw-control-item{flex:1 1 130px;min-width:110px}
	    #lwwidget .lw-label{font-size:13px;color:#9ca3af;display:flex;justify-content:space-between;margin-bottom:2px}
	    #lwwidget .lw-value{color:#f3f4f6;font-weight:600}
	    #lwwidget input[type=range]{width:100%;accent-color:#3b82f6}
	    #lwwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px;cursor:pointer}
	    #lwwidget button.active{background:#2563eb;border-color:#93c5fd}
	    #lwwidget .lw-legend{display:flex;gap:12px;flex-wrap:wrap;font-size:12px;color:#9ca3af;margin-top:4px}
	    #lwwidget .lw-swatch{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:middle}
	    </style>

	    <div class="lw-titlebar app-title">
	      <div class="lw-titlebar-headline">Adjust the medium and watch the eigenvalue phase velocities, wavefield, and seismogram respond.</div>
	      <div class="lw-titlebar-sub">drag the boundary line to set H &middot; drag &beta;/&rho; sliders &middot; drag the dispersion panel's cursor to change frequency &middot; click a curve to toggle a mode &middot; drag the band below it for the seismogram's frequency window &middot; drag the seismogram's cursor for time</div>
	    </div>

	    <div class="lw-row">
	      <div class="lw-panel" style="flex:1 1 560px">
	        <div class="lw-panel-title">Love-wave Displacement</div>
	        <canvas id="lw-wave" width="560" height="340"></canvas>
	        <div class="lw-caption" id="lw-wave-caption">Loading…</div>
	        <div style="display:flex;gap:10px;align-items:center;margin-top:6px">
	          <button id="lw-play">Play</button>
	          <div style="flex:1">
	            <div class="lw-label"><span>speed</span><span class="lw-value" id="lw-speed-v">5</span></div>
	            <input type="range" id="lw-speed" min="1" max="15" step="1" value="5">
	          </div>
	        </div>
	        <div class="lw-mini-group">
	          <div class="lw-mini-title">Medium</div>
	          <div class="lw-controls">
	            <div class="lw-control-item"><div class="lw-label"><span>β₁ (km/s)</span><span class="lw-value" id="lw-β1-v">$(w.β1)</span></div>
	              <input type="range" id="lw-β1" min="1" max="7" step="0.2" value="$(w.β1)"></div>
	            <div class="lw-control-item"><div class="lw-label"><span>ρ₁ (gm/cc)</span><span class="lw-value" id="lw-ρ1-v">$(w.ρ1)</span></div>
	              <input type="range" id="lw-ρ1" min="1" max="7" step="0.2" value="$(w.ρ1)"></div>
	            <div class="lw-control-item"><div class="lw-label"><span>β₂ (km/s)</span><span class="lw-value" id="lw-β2-v">$(w.β2)</span></div>
	              <input type="range" id="lw-β2" min="1" max="7" step="0.2" value="$(w.β2)"></div>
	            <div class="lw-control-item"><div class="lw-label"><span>ρ₂ (gm/cc)</span><span class="lw-value" id="lw-ρ2-v">$(w.ρ2)</span></div>
	              <input type="range" id="lw-ρ2" min="1" max="7" step="0.2" value="$(w.ρ2)"></div>
	          </div>
	        </div>
	      </div>
	      <div class="lw-panel" style="flex:1 1 560px">
	        <div class="lw-panel-title">Dispersion Curves &amp; Mode Selection</div>
	        <canvas id="lw-disp" width="560" height="260"></canvas>
	        <canvas id="lw-brush" width="560" height="34" style="cursor:ew-resize"></canvas>
	        <div class="lw-caption" id="lw-disp-caption">Drag the dashed cursor to change frequency. Click a curve to include/exclude that mode. Drag the shaded band below to set the seismogram's frequency window.</div>
	        <div class="lw-legend" id="lw-legend"></div>
	      </div>
	    </div>

	    <div class="lw-row">
	      <div class="lw-panel" style="flex:1 1 560px">
	        <div class="lw-panel-title">Seismogram (fundamental mode)</div>
	        <canvas id="lw-seis" width="560" height="240"></canvas>
	        <div class="lw-caption" id="lw-seis-caption">Drag to move the time cursor (blue). Red lines mark the β₁/β₂ arrival window.</div>
	        <div class="lw-mini-group">
	          <div class="lw-mini-title">Receiver</div>
	          <div class="lw-controls">
	            <div class="lw-control-item"><div class="lw-label"><span>Distance (km)</span><span class="lw-value" id="lw-xrecpos-v">$(w.xrecpos)</span></div>
	              <input type="range" id="lw-xrecpos" min="0" max="3000" step="10" value="$(w.xrecpos)"></div>
	          </div>
	        </div>
	      </div>
	      <div class="lw-panel" style="flex:1 1 560px">
	        <div class="lw-panel-title">Stationary Phase Analysis</div>
	        <canvas id="lw-phase" width="560" height="240"></canvas>
	        <div class="lw-caption" id="lw-phase-caption">Phase = k(f)·x − ωt at the current time cursor; markers show stationary points.</div>
	      </div>
	    </div>

	    <script>
	    {
	      // Pluto re-runs this whole <script> block on every re-render of this cell (e.g. a
	      // kernel restart while the tab stays connected), but window-level listeners and a
	      // running requestAnimationFrame loop survive that -- nothing else ever tears them
	      // down. Without this, each re-render stacks a fresh, duplicate set of mousemove/
	      // mouseup listeners on top of the previous ones, and an in-progress animation loop
	      // never stops, which is why the page can grind to a halt after enough re-renders.
	      if (window.__lwCleanup) { window.__lwCleanup() }
	      const lwInstance = {}
	      window.__lwCurrentInstance = lwInstance
	      const lwController = new AbortController()
	      window.__lwCleanup = () => lwController.abort()
	      const lwSignal = { signal: lwController.signal }

	      const root = document.getElementById('lwwidget')
	      const ZMAX = 150.0
	      const state = {
	        Hp: $(w.Hp), β1: $(w.β1), β2: $(w.β2), ρ1: $(w.ρ1), ρ2: $(w.ρ2), freq: $(w.freq),
	        xrecpos: $(w.xrecpos), freqmin: $(w.freqmin), freqmax: $(w.freqmax),
	        t: $(w.t), selectedModes: [$(join(w.selectedModes, ","))]
	      }
	      let userTouchedModes = false
	      let data = null

	      const MODE_COLORS = ['#3b82f6', '#ef4444', '#22c55e', '#eab308']

	      function publish() {
	        root.value = { ...state, selectedModes: state.selectedModes.slice() }
	        root.dispatchEvent(new CustomEvent('input'))
	      }
	      root.value = { ...state, selectedModes: state.selectedModes.slice() }

	      function lerp(a, b, t) { return a + (b - a) * t }
	      function clamp(v, a, b) { return Math.max(a, Math.min(b, v)) }
	      function divColorRGB(v, vmax) {
	        if (vmax <= 0) vmax = 1e-9
	        const t = clamp(v / vmax, -1, 1)
	        let r, g, b
	        if (t >= 0) { r = lerp(0, 239, t); g = lerp(0, 68, t); b = lerp(0, 68, t) }
	        else { const s = -t; r = lerp(0, 59, s); g = lerp(0, 130, s); b = lerp(0, 246, s) }
	        return [Math.round(r), Math.round(g), Math.round(b)]
	      }

	      function drawTickLabel(ctx, label, x, y, align, minX, maxX) {
	        ctx.font = '11px sans-serif'
	        ctx.fillStyle = '#9ca3af'
	        const w2 = ctx.measureText(label).width
	        let px = x
	        if (align === 'center') px = x - w2 / 2
	        else if (align === 'right') px = x - w2
	        px = clamp(px, minX, maxX - w2)
	        ctx.textAlign = 'left'
	        ctx.fillText(label, px, y)
	      }

	      // ---------- Panel A: Wavefield ----------
	      const waveCanvas = root.querySelector('#lw-wave')
	      const waveCtx = waveCanvas.getContext('2d')
	      const WAVE_M = { l: 44, r: 12, t: 12, b: 26 }
	      const offscreen = document.createElement('canvas')
	      const offCtx = offscreen.getContext('2d')

	      function waveLayout() {
	        return { plotX: WAVE_M.l, plotY: WAVE_M.t, plotW: waveCanvas.width - WAVE_M.l - WAVE_M.r, plotH: waveCanvas.height - WAVE_M.t - WAVE_M.b }
	      }
	      function waveBoundaryY() {
	        const L = waveLayout()
	        return L.plotY + L.plotH * (state.Hp / ZMAX)
	      }

	      function drawWavefield(frameIdx) {
	        const ctx = waveCtx
	        const W = waveCanvas.width, H = waveCanvas.height
	        ctx.fillStyle = '#000'; ctx.fillRect(0, 0, W, H)
	        if (!data || !data.U1_frames) {
	          ctx.fillStyle = '#9ca3af'; ctx.font = '13px sans-serif'; ctx.fillText('Loading…', 20, 40)
	          return
	        }
	        const L = waveLayout()
	        const plotX = L.plotX, plotY = L.plotY, plotW = L.plotW, plotH = L.plotH
	        const nx = data.nx, nz = data.nz
	        offscreen.width = nx; offscreen.height = nz * 2
	        const img = offCtx.createImageData(nx, nz * 2)
	        const perFrame = nz * nx
	        const base1 = frameIdx * perFrame
	        const base2 = frameIdx * perFrame
	        for (let iz = 0; iz < nz; iz++) {
	          for (let ix = 0; ix < nx; ix++) {
	            const v1 = data.U1_frames[base1 + iz * nx + ix]
	            const rgb1 = divColorRGB(v1, data.wmax)
	            let p = (iz * nx + ix) * 4
	            img.data[p] = rgb1[0]; img.data[p + 1] = rgb1[1]; img.data[p + 2] = rgb1[2]; img.data[p + 3] = 255
	          }
	        }
	        for (let iz = 0; iz < nz; iz++) {
	          for (let ix = 0; ix < nx; ix++) {
	            const v2 = data.U2_frames[base2 + iz * nx + ix]
	            const rgb2 = divColorRGB(v2, data.wmax)
	            let p = ((nz + iz) * nx + ix) * 4
	            img.data[p] = rgb2[0]; img.data[p + 1] = rgb2[1]; img.data[p + 2] = rgb2[2]; img.data[p + 3] = 255
	          }
	        }
	        offCtx.putImageData(img, 0, 0)
	        ctx.imageSmoothingEnabled = true
	        ctx.drawImage(offscreen, 0, 0, nx, nz * 2, plotX, plotY, plotW, plotH)

	        const boundaryY = waveBoundaryY()
	        ctx.strokeStyle = '#f3f4f6'; ctx.setLineDash([]); ctx.lineWidth = 3
	        ctx.beginPath(); ctx.moveTo(plotX, boundaryY); ctx.lineTo(plotX + plotW, boundaryY); ctx.stroke()
	        ctx.strokeStyle = '#9ca3af'; ctx.setLineDash([5, 4]); ctx.lineWidth = 1.5
	        ctx.beginPath(); ctx.moveTo(plotX, plotY); ctx.lineTo(plotX + plotW, plotY); ctx.stroke()
	        ctx.setLineDash([])
	        ctx.font = '12px sans-serif'; ctx.fillStyle = '#f3f4f6'
	        ctx.fillText('Free Surface', plotX + 6, plotY + 14)
	        ctx.fillText('Boundary (drag) H = ' + state.Hp.toFixed(0) + ' km', plotX + 6, boundaryY - 4)

	        ctx.fillStyle = '#9ca3af'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center'
	        drawTickLabel(ctx, '0', plotX, H - 8, 'left', 0, W)
	        drawTickLabel(ctx, data.xmax.toFixed(0), plotX + plotW, H - 8, 'right', 0, W)
	        ctx.save(); ctx.translate(12, plotY + plotH / 2); ctx.rotate(-Math.PI / 2)
	        ctx.textAlign = 'center'; ctx.fillText('Depth (km)', 0, 0); ctx.restore()
	        ctx.textAlign = 'center'; ctx.fillText('Distance (km)', plotX + plotW / 2, H - 2)
	      }

	      let waveDrag = false
	      function waveHpFromEvent(ev) {
	        const rect = waveCanvas.getBoundingClientRect()
	        const scaleY = waveCanvas.height / rect.height
	        const py = (ev.clientY - rect.top) * scaleY
	        const L = waveLayout()
	        const frac = clamp((py - L.plotY) / L.plotH, 0, 1)
	        return clamp(frac * ZMAX, 10, 100)
	      }
	      waveCanvas.addEventListener('mousedown', function (ev) {
	        const rect = waveCanvas.getBoundingClientRect()
	        const scaleY = waveCanvas.height / rect.height
	        const py = (ev.clientY - rect.top) * scaleY
	        if (Math.abs(py - waveBoundaryY()) <= 12) waveDrag = true
	      })
	      window.addEventListener('mousemove', function (ev) {
	        if (!waveDrag) return
	        state.Hp = waveHpFromEvent(ev)
	        drawWavefield(frameIdx)
	      }, lwSignal)
	      window.addEventListener('mouseup', function () {
	        if (waveDrag) { waveDrag = false; publish() }
	      }, lwSignal)

	      // ---------- Panel B: Dispersion curves ----------
	      const dispCanvas = root.querySelector('#lw-disp')
	      const dispCtx = dispCanvas.getContext('2d')
	      const DM = { l: 40, r: 12, t: 10, b: 18, split: 0.56, gap: 20 }

	      function dispLayout() {
	        const W = dispCanvas.width, H = dispCanvas.height
	        const plotW = W - DM.l - DM.r
	        const topH = (H - DM.t - DM.b - DM.gap) * DM.split
	        const botH = (H - DM.t - DM.b - DM.gap) * (1 - DM.split)
	        const topY = DM.t, botY = DM.t + topH + DM.gap
	        return { W, H, plotW, topH, botH, topY, botY, plotX: DM.l }
	      }
	      function dispFmax() {
	        return (data && data.freqgrid) ? data.freqgrid[data.freqgrid.length - 1] : 0.4995
	      }

	      function seriesRange(arrs) {
	        let mn = Infinity, mx = -Infinity
	        for (const arr of arrs) for (const v of arr) if (isFinite(v)) { if (v < mn) mn = v; if (v > mx) mx = v }
	        if (!isFinite(mn)) { mn = 0; mx = 1 }
	        const pad = (mx - mn) * 0.1 || 1
	        return [mn - pad, mx + pad]
	      }

	      function drawDispersion() {
	        const ctx = dispCtx
	        const L = dispLayout()
	        ctx.fillStyle = '#000'; ctx.fillRect(0, 0, L.W, L.H)
	        if (!data || !data.freqgrid) {
	          ctx.fillStyle = '#9ca3af'; ctx.font = '13px sans-serif'; ctx.fillText('Loading…', 20, 30)
	          return
	        }
	        const fg = data.freqgrid, fmax = fg[fg.length - 1]
	        const xOf = f => L.plotX + (f / fmax) * L.plotW
	        const [pmin, pmax] = seriesRange(data.phase_velocities)
	        const [gmin, gmax] = seriesRange(data.group_velocities)
	        const yOfTop = v => L.topY + L.topH - (v - pmin) / (pmax - pmin) * L.topH
	        const yOfBot = v => L.botY + L.botH - (v - gmin) / (gmax - gmin) * L.botH

	        ctx.strokeStyle = '#1f2937'; ctx.lineWidth = 1
	        ctx.strokeRect(L.plotX, L.topY, L.plotW, L.topH)
	        ctx.strokeRect(L.plotX, L.botY, L.plotW, L.botH)

	        ctx.font = '11px sans-serif'; ctx.fillStyle = '#9ca3af'; ctx.textAlign = 'left'
	        ctx.fillText('Phase Velocity (km/s)', L.plotX, L.topY - 2)
	        ctx.fillText('Group Velocity (km/s)', L.plotX, L.botY - 2)

	        window._lwDispCurves = []
	        for (let m = 0; m < data.phase_velocities.length; m++) {
	          const sel = state.selectedModes.indexOf(m + 1) >= 0
	          ctx.globalAlpha = sel ? 1.0 : 0.28
	          ctx.strokeStyle = MODE_COLORS[m % MODE_COLORS.length]
	          ctx.lineWidth = sel ? 2.2 : 1.4
	          const pts = []
	          ctx.beginPath(); let started = false
	          for (let i = 0; i < fg.length; i++) {
	            const v = data.phase_velocities[m][i]
	            if (!isFinite(v)) { started = false; continue }
	            const x = xOf(fg[i]), y = yOfTop(v)
	            pts.push([x, y])
	            if (!started) { ctx.moveTo(x, y); started = true } else ctx.lineTo(x, y)
	          }
	          ctx.stroke()
	          window._lwDispCurves.push({ mode: m + 1, pts })

	          const pts2 = []
	          ctx.beginPath(); started = false
	          for (let i = 0; i < fg.length; i++) {
	            const v = data.group_velocities[m][i]
	            if (!isFinite(v)) { started = false; continue }
	            const x = xOf(fg[i]), y = yOfBot(v)
	            pts2.push([x, y])
	            if (!started) { ctx.moveTo(x, y); started = true } else ctx.lineTo(x, y)
	          }
	          ctx.stroke()
	          window._lwDispCurves.push({ mode: m + 1, pts: pts2 })
	        }
	        ctx.globalAlpha = 1.0

	        const cx = xOf(state.freq)
	        ctx.strokeStyle = '#eab308'; ctx.setLineDash([4, 3]); ctx.lineWidth = 2
	        ctx.beginPath(); ctx.moveTo(cx, L.topY); ctx.lineTo(cx, L.topY + L.topH); ctx.stroke()
	        ctx.beginPath(); ctx.moveTo(cx, L.botY); ctx.lineTo(cx, L.botY + L.botH); ctx.stroke()
	        ctx.setLineDash([])
	        ctx.fillStyle = '#eab308'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center'
	        ctx.fillText(state.freq.toFixed(3) + ' Hz', cx, L.topY - 2)

	        drawTickLabel(ctx, '0', L.plotX, L.botY + L.botH + 13, 'left', 0, L.W)
	        drawTickLabel(ctx, fmax.toFixed(2) + ' Hz', L.plotX + L.plotW, L.botY + L.botH + 13, 'right', 0, L.W)

	        const legend = root.querySelector('#lw-legend')
	        let html = ''
	        for (let m = 0; m < (data.cn ? data.cn.length : 0); m++) {
	          const sel = state.selectedModes.indexOf(m + 1) >= 0
	          html += '<span data-mode="' + (m + 1) + '" style="cursor:pointer;opacity:' + (sel ? 1 : 0.45) + '">' +
	            '<span class="lw-swatch" style="background:' + MODE_COLORS[m % MODE_COLORS.length] + '"></span>' +
	            'Mode ' + (m + 1) + ': ' + data.cn[m].toFixed(2) + ' km/s</span>'
	        }
	        legend.innerHTML = html
	        legend.querySelectorAll('[data-mode]').forEach(function (el) {
	          el.addEventListener('click', function () { toggleMode(parseInt(el.getAttribute('data-mode'))) })
	        })
	      }

	      function toggleMode(m) {
	        userTouchedModes = true
	        const idx = state.selectedModes.indexOf(m)
	        if (idx >= 0) { if (state.selectedModes.length > 1) state.selectedModes.splice(idx, 1) }
	        else state.selectedModes.push(m)
	        drawDispersion()
	        publish()
	      }

	      function distToSeg(px, py, x1, y1, x2, y2) {
	        const dx = x2 - x1, dy = y2 - y1
	        const len2 = dx * dx + dy * dy
	        if (len2 === 0) return Math.hypot(px - x1, py - y1)
	        let t = ((px - x1) * dx + (py - y1) * dy) / len2
	        t = clamp(t, 0, 1)
	        return Math.hypot(px - (x1 + t * dx), py - (y1 + t * dy))
	      }

	      let freqDrag = false
	      let freqDragHappened = false
	      function dispFreqFromEvent(ev) {
	        const rect = dispCanvas.getBoundingClientRect()
	        const scaleX = dispCanvas.width / rect.width
	        const px = (ev.clientX - rect.left) * scaleX
	        const L = dispLayout()
	        const frac = clamp((px - L.plotX) / L.plotW, 0, 1)
	        return clamp(frac * dispFmax(), 0.01, 0.25)
	      }
	      dispCanvas.addEventListener('mousedown', function (ev) {
	        const rect = dispCanvas.getBoundingClientRect()
	        const scaleX = dispCanvas.width / rect.width
	        const px = (ev.clientX - rect.left) * scaleX
	        const L = dispLayout()
	        const cx = L.plotX + (state.freq / dispFmax()) * L.plotW
	        if (Math.abs(px - cx) <= 10) { freqDrag = true; freqDragHappened = false }
	      })
	      window.addEventListener('mousemove', function (ev) {
	        if (!freqDrag) return
	        freqDragHappened = true
	        state.freq = dispFreqFromEvent(ev)
	        drawDispersion()
	      }, lwSignal)
	      window.addEventListener('mouseup', function () {
	        if (freqDrag) { freqDrag = false; publish() }
	      }, lwSignal)

	      dispCanvas.addEventListener('click', function (ev) {
	        if (freqDragHappened) { freqDragHappened = false; return }
	        const rect = dispCanvas.getBoundingClientRect()
	        const scaleX = dispCanvas.width / rect.width, scaleY = dispCanvas.height / rect.height
	        const px = (ev.clientX - rect.left) * scaleX, py = (ev.clientY - rect.top) * scaleY
	        let best = null, bestD = 10
	        for (const c of (window._lwDispCurves || [])) {
	          for (let i = 0; i < c.pts.length - 1; i++) {
	            const d = distToSeg(px, py, c.pts[i][0], c.pts[i][1], c.pts[i + 1][0], c.pts[i + 1][1])
	            if (d < bestD) { bestD = d; best = c.mode }
	          }
	        }
	        if (best !== null) toggleMode(best)
	      })

	      // ---------- Frequency-band brush ----------
	      const brushCanvas = root.querySelector('#lw-brush')
	      const brushCtx = brushCanvas.getContext('2d')

	      function brushXOf(f) {
	        const fmax = dispFmax()
	        return DM.l + (f / fmax) * (brushCanvas.width - DM.l - DM.r)
	      }
	      function brushFOf(x) {
	        const fmax = dispFmax()
	        const plotW = brushCanvas.width - DM.l - DM.r
	        return clamp((x - DM.l) / plotW, 0, 1) * fmax
	      }

	      function drawBrush() {
	        const ctx = brushCtx
	        const W = brushCanvas.width, H = brushCanvas.height
	        ctx.fillStyle = '#000'; ctx.fillRect(0, 0, W, H)
	        ctx.strokeStyle = '#1f2937'; ctx.strokeRect(DM.l, 4, W - DM.l - DM.r, H - 8)
	        const x0 = brushXOf(state.freqmin), x1 = brushXOf(state.freqmax)
	        ctx.fillStyle = 'rgba(59,130,246,0.35)'
	        ctx.fillRect(x0, 4, x1 - x0, H - 8)
	        ctx.fillStyle = '#3b82f6'
	        ctx.beginPath(); ctx.arc(x0, H / 2, 6, 0, 2 * Math.PI); ctx.fill()
	        ctx.beginPath(); ctx.arc(x1, H / 2, 6, 0, 2 * Math.PI); ctx.fill()
	        ctx.fillStyle = '#e5e7eb'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center'
	        ctx.fillText(state.freqmin.toFixed(3) + '–' + state.freqmax.toFixed(3) + ' Hz', W / 2, H - 20)
	      }

	      let brushDrag = null
	      brushCanvas.addEventListener('mousedown', function (ev) {
	        const rect = brushCanvas.getBoundingClientRect()
	        const scaleX = brushCanvas.width / rect.width
	        const px = (ev.clientX - rect.left) * scaleX
	        const x0 = brushXOf(state.freqmin), x1 = brushXOf(state.freqmax)
	        brushDrag = Math.abs(px - x0) <= Math.abs(px - x1) ? 'min' : 'max'
	      })
	      window.addEventListener('mousemove', function (ev) {
	        if (!brushDrag) return
	        const rect = brushCanvas.getBoundingClientRect()
	        const scaleX = brushCanvas.width / rect.width
	        const f = brushFOf((ev.clientX - rect.left) * scaleX)
	        if (brushDrag === 'min') state.freqmin = Math.min(f, state.freqmax)
	        else state.freqmax = Math.max(f, state.freqmin)
	        drawBrush()
	      }, lwSignal)
	      window.addEventListener('mouseup', function () {
	        if (brushDrag) { brushDrag = null; publish() }
	      }, lwSignal)

	      // ---------- Panel C: Seismogram ----------
	      const seisCanvas = root.querySelector('#lw-seis')
	      const seisCtx = seisCanvas.getContext('2d')
	      const SM = { l: 48, r: 12, t: 12, b: 26 }

	      function seisLayout() {
	        const W = seisCanvas.width, H = seisCanvas.height
	        return { W, H, plotX: SM.l, plotY: SM.t, plotW: W - SM.l - SM.r, plotH: H - SM.t - SM.b }
	      }

	      function drawSeismogram() {
	        const ctx = seisCtx
	        const L = seisLayout()
	        ctx.fillStyle = '#000'; ctx.fillRect(0, 0, L.W, L.H)
	        if (!data || !data.record) {
	          ctx.fillStyle = '#9ca3af'; ctx.font = '13px sans-serif'; ctx.fillText('Loading…', 20, 30)
	          return
	        }
	        const tg = data.tgrid, tmax = tg[tg.length - 1]
	        let amax = 1e-9
	        for (const v of data.record) if (Math.abs(v) > amax) amax = Math.abs(v)
	        amax *= 1.15
	        const xOf = t => L.plotX + (t / tmax) * L.plotW
	        const yOf = v => L.plotY + L.plotH / 2 - (v / amax) * (L.plotH / 2)

	        ctx.strokeStyle = '#1f2937'; ctx.strokeRect(L.plotX, L.plotY, L.plotW, L.plotH)
	        ctx.strokeStyle = '#4b5563'; ctx.beginPath()
	        ctx.moveTo(L.plotX, yOf(0)); ctx.lineTo(L.plotX + L.plotW, yOf(0)); ctx.stroke()

	        ctx.strokeStyle = '#ef4444'; ctx.setLineDash([5, 4]); ctx.lineWidth = 1.3
	        for (const a of [data.arrival1, data.arrival2]) {
	          if (a >= 0 && a <= tmax) {
	            const x = xOf(a)
	            ctx.beginPath(); ctx.moveTo(x, L.plotY); ctx.lineTo(x, L.plotY + L.plotH); ctx.stroke()
	          }
	        }
	        ctx.setLineDash([])

	        ctx.strokeStyle = '#e5e7eb'; ctx.lineWidth = 1.4
	        ctx.beginPath()
	        for (let i = 0; i < tg.length; i++) {
	          const x = xOf(tg[i]), y = yOf(data.record[i])
	          if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y)
	        }
	        ctx.stroke()

	        const tx = xOf(state.t)
	        ctx.strokeStyle = '#3b82f6'; ctx.lineWidth = 2
	        ctx.beginPath(); ctx.moveTo(tx, L.plotY); ctx.lineTo(tx, L.plotY + L.plotH); ctx.stroke()

	        drawTickLabel(ctx, '0', L.plotX, L.H - 8, 'left', 0, L.W)
	        drawTickLabel(ctx, tmax.toFixed(0) + ' s', L.plotX + L.plotW, L.H - 8, 'right', 0, L.W)
	      }

	      let seisDrag = false
	      function seisTFromEvent(ev) {
	        const L = seisLayout()
	        const rect = seisCanvas.getBoundingClientRect()
	        const scaleX = seisCanvas.width / rect.width
	        const px = (ev.clientX - rect.left) * scaleX
	        const tg = data ? data.tgrid : null
	        const tmax = tg ? tg[tg.length - 1] : 1000
	        return clamp((px - L.plotX) / L.plotW, 0, 1) * tmax
	      }
	      seisCanvas.addEventListener('mousedown', function (ev) { seisDrag = true; state.t = seisTFromEvent(ev); drawSeismogram() })
	      window.addEventListener('mousemove', function (ev) { if (seisDrag) { state.t = seisTFromEvent(ev); drawSeismogram() } }, lwSignal)
	      window.addEventListener('mouseup', function () { if (seisDrag) { seisDrag = false; publish() } }, lwSignal)

	      // ---------- Panel D: Stationary phase ----------
	      const phaseCanvas = root.querySelector('#lw-phase')
	      const phaseCtx = phaseCanvas.getContext('2d')
	      const PM = { l: 48, r: 12, t: 12, b: 26 }

	      function drawPhase() {
	        const ctx = phaseCtx
	        const W = phaseCanvas.width, H = phaseCanvas.height
	        ctx.fillStyle = '#000'; ctx.fillRect(0, 0, W, H)
	        if (!data || !data.phaseArr) {
	          ctx.fillStyle = '#9ca3af'; ctx.font = '13px sans-serif'; ctx.fillText('Loading…', 20, 30)
	          return
	        }
	        const plotX = PM.l, plotY = PM.t, plotW = W - PM.l - PM.r, plotH = H - PM.t - PM.b
	        const fg = data.freqgrid, fmax = fg[fg.length - 1]
	        const [pmin, pmax] = seriesRange([data.phaseArr])
	        const xOf = f => plotX + (f / fmax) * plotW
	        const yOf = v => plotY + plotH - (v - pmin) / (pmax - pmin) * plotH

	        ctx.strokeStyle = '#1f2937'; ctx.strokeRect(plotX, plotY, plotW, plotH)
	        ctx.strokeStyle = '#e5e7eb'; ctx.lineWidth = 1.4
	        ctx.beginPath()
	        for (let i = 0; i < fg.length; i++) {
	          const x = xOf(fg[i]), y = yOf(data.phaseArr[i])
	          if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y)
	        }
	        ctx.stroke()

	        ctx.fillStyle = '#3b82f6'
	        for (const i of (data.minPoints || [])) {
	          const x = xOf(fg[i]), y = yOf(data.phaseArr[i])
	          ctx.beginPath(); ctx.arc(x, y, 4, 0, 2 * Math.PI); ctx.fill()
	        }
	        ctx.fillStyle = '#ef4444'
	        for (const i of (data.maxPoints || [])) {
	          const x = xOf(fg[i]), y = yOf(data.phaseArr[i])
	          ctx.beginPath(); ctx.arc(x, y, 4, 0, 2 * Math.PI); ctx.fill()
	        }

	        drawTickLabel(ctx, '0', plotX, H - 8, 'left', 0, W)
	        drawTickLabel(ctx, fmax.toFixed(2) + ' Hz', plotX + plotW, H - 8, 'right', 0, W)
	      }

	      function drawAll() {
	        const n = (data && data.nframes) ? data.nframes : 1
	        drawWavefield(Math.min(frameIdx, n - 1))
	        drawDispersion()
	        drawBrush()
	        drawSeismogram()
	        drawPhase()
	        const capt = root.querySelector('#lw-wave-caption')
	        if (data && data.cn) {
	          capt.textContent = 'Summing ' + state.selectedModes.length + ' of ' + data.cn.length + ' mode(s) at f = ' + state.freq.toFixed(3) + ' Hz'
	        }
	      }

	      root.addEventListener('lw-data', function (ev) {
	        data = ev.detail
	        const n = data.cn ? data.cn.length : 1
	        let next
	        if (!userTouchedModes) {
	          next = []
	          for (let i = 1; i <= n; i++) next.push(i)
	        } else {
	          next = state.selectedModes.filter(function (i) { return i >= 1 && i <= n })
	          if (next.length === 0) next = [1]
	        }
	        const changed = next.length !== state.selectedModes.length ||
	          next.some(function (v, i) { return v !== state.selectedModes[i] })
	        state.selectedModes = next
	        drawAll()
	        if (changed) publish()
	      })

	      // ---------- animation loop ----------
	      let playing = false, rafId = null, lastTick = 0, frameIdx = 0
	      const playBtn = root.querySelector('#lw-play')
	      playBtn.addEventListener('click', function () {
	        playing = !playing
	        playBtn.textContent = playing ? 'Pause' : 'Play'
	        playBtn.className = playing ? 'active' : ''
	        if (playing) { lastTick = 0; rafId = requestAnimationFrame(tick) }
	      })
	      const speedEl = root.querySelector('#lw-speed')
	      const speedVEl = root.querySelector('#lw-speed-v')
	      speedEl.addEventListener('input', function () { speedVEl.textContent = speedEl.value })
	      function tick(now) {
	        if (window.__lwCurrentInstance !== lwInstance) return  // a newer render replaced us; stop
	        if (!playing) return
	        const fps = parseFloat(speedEl.value)
	        if (now - lastTick > 1000 / fps) {
	          const n = (data && data.nframes) ? data.nframes : 1
	          frameIdx = (frameIdx + 1) % n
	          drawWavefield(frameIdx)
	          lastTick = now
	        }
	        rafId = requestAnimationFrame(tick)
	      }

	      // ---------- control sliders (β/ρ + receiver) ----------
	      const ids = {β1:'lw-β1', β2:'lw-β2', ρ1:'lw-ρ1', ρ2:'lw-ρ2', xrecpos:'lw-xrecpos'}
	      for (const k in ids) {
	        const el = root.querySelector('#'+ids[k])
	        const vEl = root.querySelector('#'+ids[k]+'-v')
	        el.addEventListener('input', () => { vEl.textContent = el.value })
	        el.addEventListener('change', () => {
	          state[k] = parseFloat(el.value)
	          publish()
	        })
	      }

	      drawAll()
	    }
	    </script>
	    </div>
	    """)
	end

	const _lw_ready = true
end

# ╔═╡ bec713a8-708f-42e9-a7ef-35c16f597336
begin
	_lw_ready
	WideCell(@bind lw LoveWaveInput(); max_width=1400)
end

# ╔═╡ d417ece2-0299-4a60-9299-803816aa8b1f
# The bond is `nothing` until a browser has rendered the widget: use its defaults.
lw_s = lw isa AbstractDict ? lw : Base.get(LoveWaveInput())

# ╔═╡ 60c6e69b-c00b-4ffd-87b7-237bb5d3dfcc
medium = (; Hp=Float64(lw_s["Hp"]), β₁=Float64(lw_s["β1"]), β₂=Float64(lw_s["β2"]), ρ₁=Float64(lw_s["ρ1"]), ρ₂=Float64(lw_s["ρ2"]))

# ╔═╡ b3e662bf-8d13-418b-886c-b29456d62454
crange2 = range(medium.β₁ - 0.5, stop=medium.β₁, length=100)

# ╔═╡ 3fa22358-6d61-4b9a-9aa8-39d2b1c6a917
crange3 = range(medium.β₂, stop=medium.β₂ + 0.5, length=100)

# ╔═╡ 831be90d-f2ae-4290-b195-34dbb2b2f294
begin
	freq_op = Float64(lw_s["freq"])
	xrecpos = Float64(lw_s["xrecpos"])
	freqmin = Float64(lw_s["freqmin"])
	freqmax = Float64(lw_s["freqmax"])
	tcursor = Float64(lw_s["t"])
	selected_modes = Int.(lw_s["selectedModes"])
end

# ╔═╡ 1e4a9c31-9b1a-4b64-9a49-6a7b7cb4f4e1
F = let expr = ex1 - ex2
    build_function(expr, p, ω, η₁, η₂, μ₁, μ₂, H, expression=Val{false})
end

# ╔═╡ cc173455-2ed0-42f3-a9c0-0dfdbbe982ee
cn = compute_phase_velocities(F, freq_op, medium)

# ╔═╡ 5d81104e-16a1-4c44-89a3-a492821021ed
# Every interaction re-runs this cell; the ~500 root searches only repeat when the medium changed.
cn_vec = memo!(_lw_dispersion_memo, (F, freqgrid, medium)) do
	compute_dispersion_curves(F, freqgrid, medium)
end

# ╔═╡ a37d91cc-cd80-48c0-8ecc-959df16fb925
phase_velocities = map(1:4) do i
    map(cn_vec) do c
        if (length(c) >= i)
            c[i]
        else
            missing
        end
    end
end

# ╔═╡ 6efd46bf-c421-4a29-8b8c-5a1677c16431
wavenumbers = [get_k(c, freqgrid) for c in phase_velocities];

# ╔═╡ f07922a8-e39e-42c5-876b-a291e462d583
group_velocities = [get_group_velocity(k, freqgrid) for k in wavenumbers];

# ╔═╡ 50b57477-3da3-4c67-8fc1-f19a608dbfe2
begin
	phase_arr = wavenumbers[1] .* xrecpos .- 2π .* freqgrid[2:end] .* tcursor
	min_points, _ = findminima(phase_arr)
	max_points, _ = findmaxima(phase_arr)
end

# ╔═╡ f7fcbbf9-f938-46b4-83b7-d959a0208079
let
	μ1 = medium.β₁^2 * medium.ρ₁
	μ2 = medium.β₂^2 * medium.ρ₂
	f_below = maximum(abs ∘ real, F.(inv.(crange2), 2π * freq_op, getη.(inv.(crange2), medium.β₁), getη.(inv.(crange2), medium.β₂), μ1, μ2, medium.Hp))
	f_above = maximum(abs ∘ real, F.(inv.(crange3), 2π * freq_op, getη.(inv.(crange3), medium.β₁), getη.(inv.(crange3), medium.β₂), μ1, μ2, medium.Hp))
	(min_abs_F_below_β₁ = f_below, min_abs_F_above_β₂ = f_above)
end

# ╔═╡ b0a6f5b8-9d1e-4c1a-8b8e-5e9b6e0a2c3d
U1 = let expr = substitute(u1, [ı => im, B₁ => 1, A₁ => 1])
    build_function(expr, x, z, t, p, ω, η₁, expression=Val{false})
end

# ╔═╡ 546bee6c-9faa-4f6b-a1e2-c46cae098185
begin
	record = love_seismogram(U1, @view(freqgrid[2:end]), phase_velocities[1], tgrid, xrecpos, medium.β₁, freqmin, freqmax)
	arrival1 = xrecpos / medium.β₁
	arrival2 = xrecpos / medium.β₂
end

# ╔═╡ c2f7e6c9-8a2f-4d2b-9c9f-6f8c7f1b3d4e
U2 = let
	A2p_sym = substitute(A₂ex1, [ı => im, B₁ => 1, A₁ => 1, B₂ => 0])
	expr = substitute(u2, [ı => im, B₁ => 1, A₁ => 1, B₂ => 0, A₂ => A2p_sym])
	build_function(expr, x, z, t, p, ω, η₁, η₂, μ₁, μ₂, H, expression=Val{false})
end

# ╔═╡ ac1af06b-a4cd-4112-958b-9b10e797bd5d
begin
	xgrid = range(0, stop=300, length=NX_LW)
	zgrid1 = range(0, stop=medium.Hp, length=NZ_LW)
	zgrid2 = range(medium.Hp, stop=ZMAX_LW, length=NZ_LW)

	active_cn = [cn[i] for i in selected_modes if 1 <= i <= length(cn)]
	if isempty(active_cn) && !isempty(cn)
		active_cn = [cn[1]]
	end

	U1_frames, U2_frames, wmax = love_wavefield_frames(U1, U2, active_cn, freq_op, medium, xgrid, zgrid1, zgrid2, NFRAMES_LW)
end

# ╔═╡ ea8dc149-2606-431b-85ca-a060f24cf1b3
LovePush(freqgrid[2:end], phase_velocities, group_velocities, cn, freq_op,
    NX_LW, NZ_LW, NFRAMES_LW, medium.Hp, zgrid2[end], xgrid[end], wmax,
    U1_frames, U2_frames, tgrid, record, arrival1, arrival2,
    phase_arr, min_points, max_points)

# ╔═╡ d3e8f7d0-7b3f-4e3c-8d8f-7a9d8e2c4f5a
md"""
### Verifying the Generalized Functions Against the Symbolic Derivation
`F`, `U1`, `U2` above were built with `μ₁`, `μ₂`, `H` symbolic. Here we confirm they agree,
at an arbitrary test point, with directly substituting numbers into the original symbolic
expressions (`ex1 - ex2`, `u1`, `u2`) -- i.e. that leaving those three parameters symbolic
during `build_function` didn't change the math.
"""

# ╔═╡ e4f9a8e1-6c4a-4f4d-9e9a-8b0e9f3d5a6b
let
	test_medium = (; Hp=40.0, β₁=3.2, β₂=4.8, ρ₁=2.4, ρ₂=3.6)
	μ1t = test_medium.β₁^2 * test_medium.ρ₁
	μ2t = test_medium.β₂^2 * test_medium.ρ₂
	p0, ω0 = 0.28, 1.3
	η1_0, η2_0 = getη(p0, test_medium.β₁), getη(p0, test_medium.β₂)
	x0, z0, t0 = 120.0, 15.0, 3.0

	# a fully-substituted expression (no free symbols left) doesn't auto-collapse to a
	# native number under this Symbolics version -- forcing it through build_function()
	# (a zero-arg compiled call) is what actually evaluates it, unlike Symbolics.value.
	collapse(expr) = build_function(expr, expression=Val{false})()

	F_truth = collapse(substitute(ex1 - ex2,
		[p => p0, ω => ω0, η₁ => η1_0, η₂ => η2_0, H => test_medium.Hp, μ₁ => μ1t, μ₂ => μ2t]))
	U1_truth = collapse(substitute(u1,
		[ı => im, B₁ => 1, A₁ => 1, x => x0, z => z0, t => t0, p => p0, ω => ω0, η₁ => η1_0]))
	A2p_truth = collapse(substitute(A₂ex1,
		[ı => im, B₁ => 1, A₁ => 1, B₂ => 0, H => test_medium.Hp, μ₁ => μ1t, μ₂ => μ2t, ω => ω0, η₁ => η1_0, η₂ => η2_0]))
	U2_truth = collapse(substitute(u2,
		[ı => im, B₁ => 1, A₁ => 1, B₂ => 0, A₂ => A2p_truth, x => x0, z => z0, t => t0, p => p0, ω => ω0, η₂ => η2_0]))

	(
		F_err=abs(F(p0, ω0, η1_0, η2_0, μ1t, μ2t, test_medium.Hp) - F_truth),
		U1_err=abs(U1(x0, z0, t0, p0, ω0, η1_0) - U1_truth),
		U2_err=abs(U2(x0, z0, t0, p0, ω0, η1_0, η2_0, μ1t, μ2t, test_medium.Hp) - U2_truth),
		# love_seismogram / love_wavefield_frames evaluate at t = 0 and rotate by e^{iωt}
		U1_tfactor_err=abs(U1(x0, z0, t0, p0, ω0, η1_0) - U1(x0, z0, 0.0, p0, ω0, η1_0) * cis(ω0 * t0)),
		U2_tfactor_err=abs(U2(x0, z0, t0, p0, ω0, η1_0, η2_0, μ1t, μ2t, test_medium.Hp) -
			U2(x0, z0, 0.0, p0, ω0, η1_0, η2_0, μ1t, μ2t, test_medium.Hp) * cis(ω0 * t0)),
	)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
FFTW = "~1.10.0"
Latexify = "~0.16.11"
Peaks = "~0.6.2"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.83"
Roots = "~3.0.6"
SymbolicUtils = "~4.44.0"
Symbolics = "~7.35.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.1"
manifest_format = "2.1"
project_hash = "548371942c49c878fb576946db865fcd68ee9f29"

[[deps.ADTypes]]
git-tree-sha1 = "0a81a018463de6c3f4f2c9360121c562e5add9e4"
registries = "General"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.22.3"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
registries = "General"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
registries = "General"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.45"

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

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
registries = "General"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
registries = "General"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "f54afab101687a7049833d07636418a83e9a250b"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.12"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ef2022bff55342a8c9846cdf218f62e475f0444d"
registries = "General"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.1.2"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
registries = "General"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
registries = "General"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

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

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
registries = "General"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "FunctionMaps", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c0f576ae49bd2d1bc904b9946f4783db8f0ef530"
registries = "General"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.8.1"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicPolynomials]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "StarAlgebras", "Test"]
git-tree-sha1 = "5bfabc3827dfdd164359bad6800c115a81280c00"
registries = "General"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.6"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
registries = "General"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExprTools]]
git-tree-sha1 = "d2e49e7efd29719d6f28b891b0e0e159daa9d2b4"
registries = "General"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.11"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
registries = "General"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FunctionMaps]]
deps = ["CompositeTypes", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "31bd99a57edf98990d1c21486032963955450e8d"
registries = "General"
uuid = "a85aefff-f8ca-4649-a888-c8e5398bc76c"
version = "0.1.2"

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

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7eb45fe833a5b7c51cf6d89c5a841d5967e44be3"
registries = "General"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.14.0"

    [deps.Graphs.extensions]
    GraphsSharedArraysExt = "SharedArrays"

    [deps.Graphs.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"
    SharedArrays = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
registries = "General"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "c72458f1962faeb003bf23cbdb75164fe6280906"
registries = "General"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.4"

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

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
registries = "General"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
registries = "General"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

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

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "24390f715ff0795a1c4b912d788f18c52c6abd19"
registries = "General"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.11"

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

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

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

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "60beb0717782a3bbe0f7df56decad0ef89048c23"
registries = "General"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.12"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics", "StarAlgebras"]
git-tree-sha1 = "4838893d9b035c2f6967c0d533350e1755b58a70"
registries = "General"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.19"

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

    [deps.MultivariatePolynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "dc5b2c4c111c46bc79ac4405eeb563523b39c004"
registries = "General"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.8.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
registries = "General"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.Peaks]]
deps = ["SIMD"]
git-tree-sha1 = "a9b6680fb7fb097fb6eb1210c35549218d73da84"
registries = "General"
uuid = "18e31ff7-3703-566c-8e60-38913d67486b"
version = "0.6.2"

    [deps.Peaks.extensions]
    MakieExt = "Makie"
    PlotsExt = "RecipesBase"

    [deps.Peaks.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "90b41ced6bacd8c01bd05da8aed35c5458891749"
registries = "General"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

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

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
registries = "General"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.ReadOnlyArrays]]
git-tree-sha1 = "e6f7ddf48cf141cb312b078ca21cb2d29d0dc11d"
registries = "General"
uuid = "988b38a3-91fc-5605-94a2-ee2116b3bd83"
version = "0.2.0"

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

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "7fb25a964849d90a0446366cdefca822e0e84900"
registries = "General"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.6"

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

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6368773b8b433a5c6d9f8b9427c5fbded65fc8c0"
registries = "General"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.23"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
registries = "General"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.SciMLPublic]]
git-tree-sha1 = "cf9aaf8b9ed5db993259ea8b24cf2b7ba9bd3b79"
registries = "General"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.2.4"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
registries = "General"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
registries = "General"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "97c8329c5f503d2936fb36719fe25b9f94b1ae8a"
registries = "General"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StarAlgebras]]
deps = ["LinearAlgebra", "MutableArithmetics", "SparseArrays"]
git-tree-sha1 = "235b1f9d287bbf34083b3d0829343a7942c0ad1c"
registries = "General"
uuid = "0c0c59c1-dc5f-42e9-9a8b-b5dc384a6cd1"
version = "0.3.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

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

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "73048fd086b7a169bbd7232bf60bfd43240691eb"
registries = "General"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.51"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "ab885203e8395593d65b629bd4023de089e6997b"
registries = "General"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "1.1.4"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "Graphs", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "03bbe242c7433bfca3660050d0b0cc3b4be8df71"
registries = "General"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.44.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsChainRulesCoreExt = "ChainRulesCore"
    SymbolicUtilsDistributionsExt = "Distributions"
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "AbstractPlutoDingetjes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "IntervalSets", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "Preferences", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "2a8387ae58e5c1afbd8820b5c3ddb300d15ae27f"
registries = "General"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.35.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsHypergeometricFunctionsExt = "HypergeometricFunctions"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    HypergeometricFunctions = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
registries = "General"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
registries = "General"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

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

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakCacheSets]]
git-tree-sha1 = "386050ae4353310d8ff9c228f83b1affca2f7f38"
registries = "General"
uuid = "d30d5f5c-d141-4870-aa07-aabb0f5fe7d5"
version = "0.1.0"

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
# ╠═9b02e72d-be6b-485d-a186-ba2fe2dcd6fb
# ╠═bdf3d004-653a-4b56-af7c-42ca1d6021e5
# ╟─58ac96df-06f8-43b1-b19c-4f16c86504e8
# ╟─2a41b15e-1a0e-4c92-a3a5-53603faacea1
# ╟─4479381e-b424-4599-b692-877e4a15405e
# ╟─bec713a8-708f-42e9-a7ef-35c16f597336
# ╟─d417ece2-0299-4a60-9299-803816aa8b1f
# ╟─9f106bb0-b0c6-4c5b-bdf2-f17c15693b82
# ╟─a19519d6-a96f-4a8a-8563-d5052f694aff
# ╠═8d2fd2d9-da36-43f4-8dcb-dd931c8189f5
# ╠═25f32a28-4c57-4ed0-9fc4-4a40ca31e4dd
# ╟─d8b52234-757b-441d-918d-585fd6593fcf
# ╠═b120da8e-dcd6-4280-819b-0e20b4675080
# ╠═c74c4a58-b47e-4509-aadd-7400c350ff6c
# ╠═44aa241c-bf9d-4a2a-b459-28a73a991606
# ╠═5896fa94-ab7f-42f0-9fac-265e11faba2b
# ╟─6e734951-2907-49d7-ae88-dffc463bf8a3
# ╠═2f40cb44-60af-4fbe-8f73-d8d06e99fc3b
# ╠═4d859323-12e0-4fba-a796-48995f980f33
# ╟─79c4f3a4-65b6-4e02-9786-f168584c7781
# ╠═dc709a17-ab03-406f-b232-d7658872a95e
# ╠═83e2b2b3-f84a-4926-9332-b5570d52be8b
# ╟─d6dc2422-6be3-447e-b7c1-be655971b6f7
# ╠═1a13363a-9f8f-4123-9b37-c6dcad1ccdde
# ╠═2f80cfcf-9788-41a3-8797-cd1965216738
# ╟─de304d8c-253e-4cb5-8332-db79de20991c
# ╠═30219c41-a040-4056-adda-ce64ac3c6969
# ╟─e1ab342d-84fe-4c00-9756-455436b0fad7
# ╠═ca377428-d597-47a2-b29e-f40884697bd6
# ╠═22ef3c54-fcb8-463e-8520-4d07bbceda51
# ╠═b819ce05-628a-49b6-aa57-f1ee765fb027
# ╠═53799d32-3edc-4416-a928-868ea974c6bd
# ╟─f39b1e01-c16e-405e-a2fc-85408d53c762
# ╠═b3924a35-6565-46c4-9a6e-d29b1172668e
# ╠═30f26e47-0a37-4165-b18b-a20f1148deb1
# ╠═f7b5adba-09ca-4895-8b18-5fc84b752392
# ╠═8f35d8b0-ff20-442c-b5b2-f5a2616fe345
# ╟─dab840e2-7d7c-4f65-b4ea-289babc430c8
# ╠═644637c9-c32d-4335-af59-3452263a9707
# ╠═09a8041d-6443-4e83-8443-e0e2aa05317a
# ╠═0fc4d1a1-25ac-49c9-b74b-a9f26928c652
# ╟─794a1f79-708d-4305-8aa0-d2673f34e31f
# ╠═b3e662bf-8d13-418b-886c-b29456d62454
# ╠═3fa22358-6d61-4b9a-9aa8-39d2b1c6a917
# ╠═bfb54c6f-a52c-4a10-87c7-9788f2f63629
# ╠═cbbae944-5377-44a4-ae32-868a75625248
# ╟─b5c5f1d0-0075-4727-a6a9-8b20d6f65bc4
# ╟─ad269679-e1f9-4d9d-827f-468e36f27bfc
# ╟─b2858a78-4498-4374-9ee0-c68f3dfca6e8
# ╠═039a52c9-01b8-4028-9eae-1ff5b779fd78
# ╠═cc173455-2ed0-42f3-a9c0-0dfdbbe982ee
# ╟─8db66193-eb99-478e-bdad-44b0e4e18dad
# ╠═a69a71b4-6dbc-4de7-a8d7-5cc3ad744d94
# ╠═7c1e5a20-3b8d-4f61-9a0e-2d6b4c8e1f37
# ╠═5d81104e-16a1-4c44-89a3-a492821021ed
# ╠═a37d91cc-cd80-48c0-8ecc-959df16fb925
# ╠═e4c1cda0-7366-4aa0-a616-47d86a1c9502
# ╠═6ae71567-4501-4cf5-8055-47b235829f14
# ╠═6efd46bf-c421-4a29-8b8c-5a1677c16431
# ╠═f07922a8-e39e-42c5-876b-a291e462d583
# ╟─df35e923-b43c-4809-a9db-32a371fc9010
# ╠═66ae626f-2543-4897-946f-8e51ecbb95e0
# ╠═896fb716-861a-44d1-b073-45c47773a4f8
# ╠═4ee3be7e-7fa9-4f4a-b13a-3d20a8863f04
# ╠═840a743a-e298-4f99-ae6e-11c85b6f5bc5
# ╠═80bb6ae1-c48b-4be1-9786-c94c6e8680a1
# ╟─7bd19e44-b4f5-43f3-b6be-b557ca95c67f
# ╟─16c826bf-fa48-423c-bfec-195f7fc502f8
# ╠═60c6e69b-c00b-4ffd-87b7-237bb5d3dfcc
# ╟─28fbc1e8-2323-4d0d-8f43-e6c1041c7f8d
# ╠═f7fcbbf9-f938-46b4-83b7-d959a0208079
# ╟─361a3181-857a-4454-9fe6-64c8ae3be9cf
# ╠═749ed9e6-a16e-49ee-9f70-909aa0493b0d
# ╠═3f9b2d71-6a4e-4c18-b5d3-8e1f0a7c2b94
# ╠═a5d0c3e8-91f2-4b7a-8c6d-0e4f2b9a1d53
# ╠═831be90d-f2ae-4290-b195-34dbb2b2f294
# ╠═546bee6c-9faa-4f6b-a1e2-c46cae098185
# ╠═50b57477-3da3-4c67-8fc1-f19a608dbfe2
# ╠═ac1af06b-a4cd-4112-958b-9b10e797bd5d
# ╠═648572e8-2cf4-4a1f-932d-02aa2160b6d9
# ╠═beb71023-cd94-4fd2-8108-c68ad8b58055
# ╠═ea8dc149-2606-431b-85ca-a060f24cf1b3
# ╟─42190b26-0b30-4b2d-bc1f-dbb3323d7bb4
# ╠═e80a6231-6f27-4fe4-a862-32976ba3ab92
# ╠═1e4a9c31-9b1a-4b64-9a49-6a7b7cb4f4e1
# ╠═b0a6f5b8-9d1e-4c1a-8b8e-5e9b6e0a2c3d
# ╠═c2f7e6c9-8a2f-4d2b-9c9f-6f8c7f1b3d4e
# ╟─d3e8f7d0-7b3f-4e3c-8d8f-7a9d8e2c4f5a
# ╠═e4f9a8e1-6c4a-4f4d-9e9a-8b0e9f3d5a6b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
