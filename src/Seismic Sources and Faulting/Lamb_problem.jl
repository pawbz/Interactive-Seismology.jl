### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Lamb's Problem"
#> layout = "layout.jlhtml"
#> tags = ["pointsource"]
#> description = "Superposition of planewaves"

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

# ╔═╡ 897afffa-77e8-11ef-1a54-c73d1df8f6a4
using Bessels, PlutoPlotly, FFTW, PlutoUI

# ╔═╡ aa38e7f6-b8d2-4270-9142-1aa688041eb4
using LinearAlgebra

# ╔═╡ 4621f804-6a44-4c46-af37-3d364ba74cfe
using Roots

# ╔═╡ 9b488729-aadc-4c1d-b971-6931d4bc9a08
using HypertextLiteral: @htl

# ╔═╡ 6188ff3e-cfd5-4c9c-aa39-619dc280d494
TableOfContents()

# ╔═╡ 1bbedd43-9e69-4ba0-a251-855f1f63dcf8
md"""# Lamb's Problem — Reflectivity Method

Line up receivers at increasing distance from a source and record them together — a **record
section** — and two arrivals separate out on their own: a **direct** wave whose travel time
grows linearly with offset, and a **reflection** off the interface below, curving along a
hyperbola instead. Push the receivers far enough and a *third* thing appears where none of the
above predicts one: a **head wave**, riding along the interface at the faster medium's own
velocity. This notebook is about where that third arrival actually comes from.

The tool for answering that is the same trick underlying the whole notebook: any point source's
spherical wavefield can be written as a superposition (an integral) of *plane* waves, each one a
genuine ray leaving the source at its own take-off angle, related to how far each plane wave's
horizontal slowness `p` is from the interface's own **critical slowness**. Below that critical
value every plane wave reflects and transmits like an ordinary Zoeppritz interface, contributing
to the direct and reflected arrivals. At and beyond it, transmission goes evanescent — total
internal reflection — and that whole extra range of plane waves is what builds the head wave.
The hero widget below lets you switch that range on and off directly and watch the head wave
appear and disappear.

The **Reflectivity Method** (`` R(p) `` below, via a numerically stable Kennett-style recursion,
not the classical textbook's single-interface case) computes the exact reflection response of an
arbitrarily thick stack of layers, combining every interface's Zoeppritz coefficients with the
phase delay of propagating through each layer. Combined with the point-source decomposition (the
**Sommerfeld/Weyl integral**, below), it synthesizes real time-domain seismograms: superpose
enough plane waves, each correctly reflected and phase-delayed, and inverse-Fourier-transform
the sum.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 72aca1bf-3997-4be3-b316-921b45480c1a
md"""## What You Just Saw — Reading the Widget

Three panels in the widget above tell one connected story. The **record section** is
[`get_wavefield`](@ref)'s frequency-domain output, inverse-FFT'd back to the time domain —
moveout (the diagonal slant of an arrival across offsets) is the direct visual signature of a
finite wave speed. Toggle `Direct` / `Reflected P` / `Reflected P→S` to see each arrival's own
slope, then toggle the head wave off and compare the **trace at the farthest receiver**: with
`Reflected P only` isolated, switching the head wave off drops that trace's amplitude by roughly
two orders of magnitude — most of what makes a post-critical reflection strong is exactly the
continuum of plane waves the toggle removes.

The **`` R(p) `` panel** shows why, directly, instead of only through its effect on a seismogram.
`` \mathrm{Im}(R_{pp})`` and `` \mathrm{Im}(R_{ps})``, swept across ray parameter `` p `` at a
representative frequency, sit at zero below each critical slowness (dashed lines) and kink
sharply upward right at it: below that point every plane wave reflects and transmits like an
ordinary Zoeppritz interface; at and beyond it, the transmitted wave in `layers[2]` goes
evanescent — total internal reflection — and reflectivity's whole analytic character changes.
That kink, not a single ray at one exact angle, is the head wave's mathematical signature: the
head wave is built from the *entire* post-critical continuum of plane waves, not from one plane
wave arriving at exactly the critical angle. [`get_wavefield`](@ref)'s `pmax` argument switches
it on and off by literally cutting the wavenumber integral off at this same slowness — exactly
what the widget's head-wave toggle does.

An elastic (P-SV) interface has **two** critical slownesses, `` p_c^P `` and `` p_c^S`` (an
acoustic interface only ever has one), so there are, in general, two separate head waves here,
one per converted mode, each announced by its own kink in the panel above.
"""

# ╔═╡ d1f5e223-cc49-4fe8-9b2e-c1ccecf2315a
md"""Reference frequency for the causal-Q dispersion each layer's `Qp`/`Qs` introduces
(velocities are only exactly `vp`/`vs` at this frequency; elsewhere they drift slightly,
per [`causal_velocity`](@ref)): $(@bind f_ref_hz PlutoUI.Slider(0.05:0.05:5.0, default=1.0, show_value=true)) Hz"""

# ╔═╡ 614a3832-7bc7-4bd3-a09d-50b76afdd7cf
ω_ref = 2π * f_ref_hz

# ╔═╡ 7a3bd5df-0f9f-489c-b522-4098c325c0c0
md"""## Conical Wave

A point source radiates a spherical wave, and a spherical wave has no simple plane-wave
reflection coefficient of its own. The way around this — the **Sommerfeld/Weyl integral
representation** — is to write the spherical wave as a *superposition* of plane waves instead,
one for every horizontal wavenumber `` k ``:

``\frac{e^{i\omega R/v}}{R} = i\int_0^\infty \frac{k}{k_z}\,J_0(kr)\,e^{ik_z|\delta z|}\,dk``,
with `` k_z^2 = (\omega/v)^2 - k^2 ``,

where `` R=\sqrt{r^2+\delta z^2} `` is the straight-line source-receiver distance, `` r `` the
horizontal range, `` \delta z `` the vertical offset, and `` J_0 `` the zeroth-order Bessel
function. Each term in this integral — one value of `` k `` — genuinely *is* a plane wave (a
"conical" wave, since `` J_0(kr) `` itself is a superposition of waves conical about the source
axis), so it has a perfectly ordinary reflection coefficient once it reaches a layer interface.
[`conical_wave`](@ref) below evaluates the integrand at one `` k ``; [`get_wavefield`](@ref)
sums it (numerically, as a Riemann sum — not with `QuadGK`, despite what an earlier draft of
this notebook claimed) over the whole `` k `` range to reconstruct the full field.
"""

# ╔═╡ e2fea7d5-00a3-4796-ad51-26f2dcffa55b
md"## Medium"

# ╔═╡ 5ea371a8-8035-4ba3-ab83-a9ffa2e3d504
begin
	struct Layer
	    thickness::Float64   # km
	    vp::Float64          # km/s
	    vs::Float64          # km/s
	    rho::Float64         # g/cm³
	    Qp::Float64          # P-wave quality factor
	    Qs::Float64          # S-wave quality factor
	end
	Layer(th, vp, vs, rho) = Layer(th, vp, vs, rho, 1000.0, 1000.0)
end

# ╔═╡ beef0000-0000-4000-8000-000000000000
# GUTENBERG_MODEL is a whole-Earth model -- even block-averaged down to a few layers its
# shallowest interface is tens of km down, which pushes the critical/crossover distances for
# a head wave well past any receiver array worth plotting. This notebook's default is instead
# a simple, shallow, two-layer crustal setup sized specifically so that a source at
# LP_SOURCE_DEPTH and receivers out to a realistic array width both see reflection moveout
# AND the head wave comfortably: with this velocity contrast and layer thickness, the
# reflected/head-wave crossover offset works out to ~46 km, comfortably inside the ±100 km
# receiver array below -- confirmed live, not just estimated (see the Seismograms section).
# Nothing stops a reader from adding more layers or loading the Gutenberg/uniform presets
# themselves; this is only the notebook's own starting point.
const LAMB_DEFAULT_LAYERS = [
    Layer(12.0, 4.0, 2.3, 2.4, 1000.0, 1000.0),
    Layer(Inf, 7.0, 4.0, 2.9, 1000.0, 1000.0),
]

# ╔═╡ f047c190-0bf3-4763-9f31-b4272e837dd2
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

# ╔═╡ 360de109-d7db-4402-bca8-5b39c6f17da9
md"## Derive Reflection Coefficients"

# ╔═╡ 687b7062-ee25-4498-948f-43b187e0ccfa
"""
    one_way_phase(q, ω, h; clip=80.0)

The phase (and, for an evanescent `q`, the decay) accumulated propagating one-way across a
layer of thickness `h` at vertical slowness `q` and angular frequency `ω`:
`` \\exp(i\\omega\\,\\mathrm{Re}(q)\\,h) \\cdot \\exp(-\\omega\\,\\mathrm{Im}(q)\\,h) ``.
`h < 0` (used for an *upgoing* wave crossing the same layer) flips both signs. The `clip`
guards the exponent against overflow for a very thick/high-frequency evanescent layer — see
[`compound_matrix_psv`](@ref), whose Dunkin-style propagator this is the building block of.
`h = Inf` (a half-space) returns `1` unconditionally: nothing propagates "across" an infinite
layer inside this function; the half-space's radiation condition is handled separately, inside
[`reflectivity_matrix`](@ref).
"""
function one_way_phase(q, ω::Real, h::Real; clip=80.0)
    if !isfinite(h)
        return one(q)
    end
    atten = clamp(-ω * imag(q) * h, -clip, clip)
    phase = exp(1im * ω * real(q) * h)
    return phase * exp(atten)
end

# ╔═╡ 8902f19a-6b93-458a-9198-f2322b2f5ab5
"""
    causal_velocity(v, Q, ω, ω_ref)

Frequency-dependent, causal (Kramers-Kronig-consistent) correction to a nominal velocity `v`
with quality factor `Q`, following the standard nearly-constant-Q dispersion law: velocity
increases logarithmically away from the reference frequency `ω_ref` (where it equals `v`
exactly), and a small constant imaginary part `-v/(2Q)` encodes the accompanying intrinsic
attenuation. `Q ≤ 0` or non-finite `Q` is treated as lossless/non-dispersive (`v` returned
exactly, as a `Complex` for type stability with the dispersive branch).
"""
function causal_velocity(v, Q, ω, ω_ref)
    if Q <= 0 || !isfinite(Q)
        return complex(v)
    end
    ωref = max(ω_ref, 1e-5)
    α = atan(1 / Q) / π
    mag = (abs(ω) / ωref)^α
    return v * mag * (1 - 1im / (2Q))
end

# ╔═╡ 7ca06983-11c1-4203-9bd0-0cc1a461b74f
"""
    vertical_slowness(v, p)

Vertical slowness `` q `` for a plane wave of velocity `v` and horizontal slowness (ray
parameter) `p`, from `` q^2 = 1/v^2 - p^2 ``, with the branch chosen for a physically decaying
solution: **propagating** (`` p < 1/v ``, sub-critical) returns the negative real root — the
downward-propagation sign convention every wavefield/propagator function in this notebook
shares; **evanescent** (`` p > 1/v ``, post-critical or trapped) returns the positive-imaginary
root, so `` \\exp(i\\omega q z) `` decays (not grows) with increasing depth `z`. Always returns a
`Complex` (the propagating branch has zero imaginary part) rather than switching between `Real`
and `Complex` by value — broadcasting this over a `k`/`p` grid that straddles the critical
slowness needs one concrete element type, or every caller up the chain (`conical_wave`,
`modal_basis_psv`, `solve_interface`, ...) degrades to boxed `Number`/`Any` arithmetic; see this
notebook's Performance self-check for the measured cost of getting this wrong.
"""
function vertical_slowness(v, p)
    arg = (1 / abs2(v)) - abs2(p)
    return arg >= 0 ? complex(-sqrt(arg)) : 1im * sqrt(-arg)
end

# ╔═╡ 05ec38ea-3431-490c-bc38-24f8c1b2d54f
"""
    conical_wave(k, ω, r, δz, layer::Layer, ωref)

One term (one horizontal wavenumber `k`) of the Sommerfeld/Weyl integral above: the
contribution of a single conical wave to the point-source field in `layer`, at range `r` and
vertical offset `δz` from the source. `ω`/`r`/`δz` are meant to broadcast together (typically
`ω` down rows, `r`/`δz` across columns, matching [`get_wavefield`](@ref)'s `param.Ω`/`param.R`/
`param.Z` grids), so one call fills an entire frequency-range grid for this `k`.

Uses [`causal_velocity`](@ref) and [`vertical_slowness`](@ref) exactly as the reflectivity
machinery below does, so the direct wave and every reflected/converted phase share one
consistent, causally-dispersive medium.
"""
function conical_wave(k, ω, r, δz, layer::Layer, ωref)
    c = causal_velocity.(layer.vp, layer.Qp, ω, ωref)
    q = vertical_slowness.(c, k ./ ω)
    kz = ω .* q
    return Bessels.besselj0.(k .* r) .* exp.(im .* kz .* δz) .* k ./ kz ./ im
end

# ╔═╡ c9bf8e13-eee2-45ae-9ce5-4901c344c8f6
"""
    wavefield_components(p, qP, qS, λ, μ, mode::Symbol, direction::Symbol)

Displacement-stress polarization vector `(ux, uz, sxz, szz)` for a plane P or SV wave
(`mode ∈ (:P, :SV)`) traveling `:down` or `:up` (`direction`), at horizontal slowness `p` and
vertical slownesses `qP`/`qS`, in a medium with Lamé parameters `λ`, `μ`. This is the elementary
building block every modal/interface calculation below assembles: [`modal_basis_psv`](@ref)
collects all four (P↓, P↑, SV↓, SV↑) into one eigenvector matrix, and [`solve_interface`](@ref)
matches these components across an interface to get Zoeppritz's reflection/transmission
coefficients. Does not itself include any `z`-dependent propagation phase — that's
[`one_way_phase`](@ref)'s job.
"""
function wavefield_components(p, qP, qS, λ, μ, mode::Symbol, direction::Symbol)
    sign = direction == :down ? 1.0 : -1.0
    if mode == :P
        ux = 1im * p
        uz = 1im * sign * qP
        sxz = -2 * μ * p * qP * sign
        szz = -(λ * (p^2 + qP^2) + 2 * μ * qP^2)
    else
        ux = -1im * sign * qS
        uz = 1im * p
        sxz = μ * (qS^2 - p^2)
        szz = -2 * μ * sign * p * qS
    end
    return ux, uz, sxz, szz
end

# ╔═╡ e6bbca73-e59c-4b6a-854e-438222778110
# Build modal basis matrix for P-SV (columns: P↓, P↑, SV↓, SV↑)
function modal_basis_psv(layer::Layer, p, ω, ω_ref)
    vp = causal_velocity(layer.vp, layer.Qp, ω, ω_ref)
    vs = causal_velocity(layer.vs, layer.Qs, ω, ω_ref)
    rho = layer.rho
    qP = vertical_slowness(vp, p)
    qS = vertical_slowness(vs, p)
    λ, μ = rho * (vp^2 - 2vs^2), rho * vs^2

    # Columns correspond to modal eigenvectors in displacement-stress space
    col_Pd = collect(wavefield_components(p, qP, qS, λ, μ, :P, :down))
    col_Pu = collect(wavefield_components(p, qP, qS, λ, μ, :P, :up))
    col_Sd = collect(wavefield_components(p, qP, qS, λ, μ, :SV, :down))
    col_Su = collect(wavefield_components(p, qP, qS, λ, μ, :SV, :up))

    V = hcat(col_Pd, col_Pu, col_Sd, col_Su)
    return V, qP, qS
end

# ╔═╡ d60d4d65-8896-40ee-b52f-17d66999c727
"""
    compound_matrix_psv(layer::Layer, p, ω, ω_ref; is_halfspace=false)

Exact Dunkin-style propagator for P–SV using modal basis (P↓, P↑, SV↓, SV↑).
Returns exponent `exa` (for potential scaling) and 4×4 propagator `D` such that
state_out = D * state_in. Uses one-way phases with overflow clamping.
"""
function compound_matrix_psv(layer::Layer, p, ω, ω_ref; is_halfspace=false)
    if is_halfspace || !isfinite(layer.thickness)
        return 0.0, Matrix{ComplexF64}(I, 4, 4)
    end

    h = layer.thickness
    V, qP, qS = modal_basis_psv(layer, p, ω, ω_ref)

    # One-way phases for each mode (down/up)
    # Downgoing: exp(+iqωh), Upgoing: exp(-iqωh) = conjugate for real q
    phase_Pd = one_way_phase(qP, ω, h)
    phase_Pu = one_way_phase(qP, ω, -h)  # upgoing = negative distance
    phase_Sd = one_way_phase(qS, ω, h)
    phase_Su = one_way_phase(qS, ω, -h)  # upgoing = negative distance

    Phase = Diagonal([phase_Pd, phase_Pu, phase_Sd, phase_Su])

    # Full propagator via eigenbasis similarity transform
    D = V * Phase * inv(V)

    # Exponent not pulled out explicitly (already in phases)
    return 0.0, D
end

# ╔═╡ da2a2e26-f217-41e0-8485-6517af621d2a
"""
    solve_interface(l1::Layer, l2::Layer, p, incident_mode, incident_side, ω, ω_ref)

Solve for the 4 unknown reflection/transmission amplitudes (`R = [RP, RS]`, `T = [TP, TS]`) of
a single plane wave (`incident_mode ∈ (:P, :SV)`) hitting the interface between `l1` (above) and
`l2` (below) from `incident_side ∈ (:top, :bottom)`, by matching all four
[`wavefield_components`](@ref) (continuity of `ux`, `uz`, `sxz`, `szz`) across the boundary —
the standard 4×4 linear-algebra form of the Zoeppritz equations. [`zoeppritz_interface`](@ref)
calls this 4 times (2 incident modes × 2 sides) to assemble the full 2×2 reflection/transmission
matrices for one interface.
"""
function solve_interface(l1::Layer, l2::Layer, p, incident_mode::Symbol, incident_side::Symbol, ω, ω_ref)
    α1 = causal_velocity(l1.vp, l1.Qp, ω, ω_ref)
    β1 = causal_velocity(l1.vs, l1.Qs, ω, ω_ref)
    ρ1 = l1.rho
    α2 = causal_velocity(l2.vp, l2.Qp, ω, ω_ref)
    β2 = causal_velocity(l2.vs, l2.Qs, ω, ω_ref)
    ρ2 = l2.rho
    qP1, qS1 = vertical_slowness(α1, p), vertical_slowness(β1, p)
    qP2, qS2 = vertical_slowness(α2, p), vertical_slowness(β2, p)
    λ1, μ1 = ρ1 * (α1^2 - 2β1^2), ρ1 * β1^2
    λ2, μ2 = ρ2 * (α2^2 - 2β2^2), ρ2 * β2^2

    # incident side setup
    inc_layer, ref_layer, tran_layer = incident_side == :top ? (l1, l1, l2) : (l2, l2, l1)
    qP_inc, qS_inc = incident_side == :top ? (qP1, qS1) : (qP2, qS2)
    λ_inc, μ_inc = incident_side == :top ? (λ1, μ1) : (λ2, μ2)
    qP_tr, qS_tr = incident_side == :top ? (qP2, qS2) : (qP1, qS1)
    λ_tr, μ_tr = incident_side == :top ? (λ2, μ2) : (λ1, μ1)
    dir_inc = incident_side == :top ? :down : :up
    dir_ref = incident_side == :top ? :up : :down
    dir_tr = incident_side == :top ? :down : :up

    inc = wavefield_components(p, qP_inc, qS_inc, λ_inc, μ_inc, incident_mode, dir_inc)
    refP = wavefield_components(p, qP_inc, qS_inc, λ_inc, μ_inc, :P, dir_ref)
    refS = wavefield_components(p, qP_inc, qS_inc, λ_inc, μ_inc, :SV, dir_ref)
    trP = wavefield_components(p, qP_tr, qS_tr, λ_tr, μ_tr, :P, dir_tr)
    trS = wavefield_components(p, qP_tr, qS_tr, λ_tr, μ_tr, :SV, dir_tr)

    A = zeros(ComplexF64, 4, 4)
    b = -collect(inc)
    A[:, 1] = collect(refP)
    A[:, 2] = collect(refS)
    A[:, 3] = -collect(trP)
    A[:, 4] = -collect(trS)

    x = A \ b
    R = x[1:2]
    T = x[3:4]
    return R, T
end

# ╔═╡ 43c34152-fc5f-4497-882f-7f42bd4e6b99
"""
    zoeppritz_interface(l1::Layer, l2::Layer, p, ω, ω_ref)

Full 2×2 reflection/transmission matrices for the interface between `l1` and `l2`, both
directions at once, via 4 calls to [`solve_interface`](@ref). Matrix convention throughout this
notebook: **columns index the incident mode, rows the reflected/transmitted mode** (column 1 = P
incident, column 2 = SV incident; row 1 = P, row 2 = SV) — so element `[2,1]` of `Rdown` is the
P→S (mode-converted) reflection for a downgoing incident wave, `[1,1]` is P→P, and so on.
Returns `(Rdown, Tdown, Rup, Tup)`: reflection/transmission for a wave incident from *above*
(`:top`, into `l2`) and from *below* (`:bottom`, into `l1`) respectively — exactly what
[`reflectivity_matrix`](@ref)'s layer-by-layer recursion needs at each interface.
"""
function zoeppritz_interface(l1::Layer, l2::Layer, p, ω, ω_ref)
    Rdown = zeros(ComplexF64, 2, 2)
    Tdown = zeros(ComplexF64, 2, 2)
    Rup = zeros(ComplexF64, 2, 2)
    Tup = zeros(ComplexF64, 2, 2)
    for (j, mode) in enumerate((:P, :SV))
        r1, t1 = solve_interface(l1, l2, p, mode, :top, ω, ω_ref)
        r2, t2 = solve_interface(l1, l2, p, mode, :bottom, ω, ω_ref)
        Rdown[:, j] = r1
        Tdown[:, j] = t1
        Rup[:, j] = r2
        Tup[:, j] = t2
    end
    return Rdown, Tdown, Rup, Tup
end

# ╔═╡ f25f798f-0ecc-4931-b8b5-9c958b83850f
begin
"""
    reflectivity_matrix(layers::Vector{Layer}, p, ω, ω_ref)

Full 2×2 reflectivity matrix seen at the top of layer 1 (just below the free surface — there is
no free surface in this recursion yet, see the To-Do), for a unit downgoing wave in layer 1.
Column `j`/row `i` convention matches [`zoeppritz_interface`](@ref): element `[1,1]` is P→P,
`[2,1]` is P→S, etc.

Uses Kennett's R/T recursion, working from the half-space upward: at the top of the half-space
there is (by definition) nothing left to reflect off, so `Rplus = 0` there; each step up adds
one more layer by propagating the current effective reflectivity `Rbar` through that layer's own
thickness ([`one_way_phase`](@ref)) and combining it with that layer's own interface response
([`zoeppritz_interface`](@ref)) via `Rplus = Rd + Td*Rbar*inv(I - Ru*Rbar)*Tu` — the standard
formula for the *effective* reflectivity of an interface backed by more layers, accounting for
every internal multiple exactly (not just a first-order approximation), since the recursion
itself already contains lower interfaces' full effect via `Rbar`.
"""
function reflectivity_matrix(layers::Vector{Layer}, p, ω, ω_ref)
    N = length(layers)
    if N <= 1 || iszero(ω)
        return zeros(ComplexF64, 2, 2)  # no interface -> no reflection
    end
    Rplus = zeros(ComplexF64, 2, 2)  # at top of halfspace: no reflection

    for ℓ in (N - 1):-1:1
        lay_dn = layers[ℓ + 1]
        E = Matrix{ComplexF64}(I, 2, 2)
        if isfinite(lay_dn.thickness)
            vp = causal_velocity(lay_dn.vp, lay_dn.Qp, ω, ω_ref)
            vs = causal_velocity(lay_dn.vs, lay_dn.Qs, ω, ω_ref)
            qP = vertical_slowness(vp, p)
            qS = vertical_slowness(vs, p)
            eP = one_way_phase(qP, ω, lay_dn.thickness)
            eS = one_way_phase(qS, ω, lay_dn.thickness)
            E = Diagonal([eP, eS]) |> Matrix
        end
        Rbar = E * Rplus * E

        Rd, Td, Ru, Tu = zoeppritz_interface(layers[ℓ], layers[ℓ + 1], p, ω, ω_ref)
        Den = I - Ru * Rbar
        Rplus = Rd + Td * (Rbar * (Den \ Tu))
    end

    return Rplus
end

"""
    reflectivity_pp(layers, p, ω, ω_ref)

P→P reflectivity at the top of layer 1 — `reflectivity_matrix(layers, p, ω, ω_ref)[1, 1]`.
"""
reflectivity_pp(layers::Vector{Layer}, p, ω, ω_ref) = reflectivity_matrix(layers, p, ω, ω_ref)[1, 1]

"""
    reflectivity_ps(layers, p, ω, ω_ref)

P→S (mode-converted) reflectivity at the top of layer 1, for a unit downgoing P in layer 1 —
`reflectivity_matrix(layers, p, ω, ω_ref)[2, 1]`. This is what drives the `Psreflect` phase in
[`get_phase`](@ref): the same P-source Sommerfeld integral, weighted by the P→S coefficient
instead of P→P.
"""
reflectivity_ps(layers::Vector{Layer}, p, ω, ω_ref) = reflectivity_matrix(layers, p, ω, ω_ref)[2, 1]
end

# ╔═╡ dea0645d-cb7c-4488-913b-ba225595aceb
let
    # A homogeneous half-space has no interface at all, so a downgoing P wave should see
    # exactly zero reflection back up -- regardless of ray parameter.
    test_layers_homo = [Layer(Inf, 6.0, 3.5, 2.7)]
    p_test = [0.0, 0.05, 0.1, 0.15, 0.19]
    ω_test = 2π * 1.0
    Rpp_vals = [reflectivity_pp(test_layers_homo, p, ω_test, ω_ref) for p in p_test]
    max_abs = maximum(abs.(Rpp_vals))
    @assert max_abs < 1e-10

    md"""
    !!! correct "Self-check"
        A homogeneous half-space has no interface to reflect off: `reflectivity_pp` returns
        $(round(max_abs, sigdigits=2)) (machine precision) at every ray parameter tested, for
        `` p`` = $(join(p_test, ", ")) s/km.
    """
end

# ╔═╡ 579cfd9f-e872-4ec5-b913-c90f8e183247
let
    # NOTE: this test's middle layer originally read Layer(20.0, 5.0, 5.0, 2.5) -- vp==vs,
    # which makes lambda = rho*(vp^2-2vs^2) = -62.5, a non-positive-definite (physically
    # invalid) elastic tensor. Fixed to vs=3.0, matching the vp/vs ratio used everywhere else
    # in this notebook's test layers -- with a genuinely valid medium, |Rpp| stays <= 1 for
    # every ray parameter tested below (see the next self-check), which was NOT true of the
    # original typo'd layer (max|Rpp| ~ 1.044 there, confirmed directly while fixing this).
    test_layers_2 = [
        Layer(40.0, 4.0, 3.0, 2.5),
        Layer(10.0, 5.0, 3.0, 2.5),
        Layer(Inf, 6.0, 3.5, 2.7),
    ]
    p_crit = 1.0 / 6.0  # critical P-wave ray parameter: p_crit = 1/vp of the half-space
    p_range = range(0, p_crit * 1.2, length=50)
    Rpp_phase = [angle(reflectivity_pp(test_layers_2, p, 2π * 1.0, ω_ref)) for p in p_range]
    plot(p_range, Rpp_phase, xlabel="Ray Parameter p (s/km)", ylabel="phase(Rpp) (rad)")
end

# ╔═╡ a1e5796b-6f2f-44c6-b5bf-df6c9ee20960
let
    # a genuinely valid medium (see the note on the previous self-check): amplitude
    # reflectivity should stay <= 1 below the critical angle for these realistic velocities.
    test_layers = [
        Layer(20.0, 5.0, 3.0, 2.5),
        Layer(20.0, 5.0, 3.0, 2.5),
        Layer(Inf, 7.0, 4.0, 3.0),
    ]
    p_test = range(0, 0.19, length=100)
    Rpp_test = [reflectivity_pp(test_layers, p, 2π * 1.0, ω_ref) for p in p_test]
    max_abs = maximum(abs.(Rpp_test))
    @assert max_abs <= 1.0 + 1e-8

    md"""
    !!! correct "Self-check"
        Over the full sub-critical ray-parameter range, `` \max|R_{pp}| = ``
        $(round(max_abs, digits=4)) — at or below 1, as physically expected for a real
        (non-negative-definite-violating) elastic medium.
    """
end

# ╔═╡ 55e64939-9298-4a15-a2b0-0d13c23d03dc
let
    # normal incidence (p=0): the acoustic-impedance reflection formula, an exact,
    # independent analytic cross-check unrelated to the P-SV machinery above.
    normal_incidence_pp(l1, l2) = (l2.rho * l2.vp - l1.rho * l1.vp) / (l2.rho * l2.vp + l1.rho * l1.vp)
    test_layers = [Layer(10.0, 5.0, 3.0, 2.5), Layer(Inf, 6.0, 3.5, 2.7)]
    R_analytical = normal_incidence_pp(test_layers[1], test_layers[2])
    R_numerical = reflectivity_pp(test_layers, 0.0, 2π * 1.0, ω_ref)
    err = abs(R_analytical - real(R_numerical))
    @assert isapprox(R_analytical, real(R_numerical); atol=1e-8)

    # and: no mode conversion is possible at normal incidence -- P->S reflectivity must be
    # exactly zero there, a clean, independent validation of reflectivity_ps specifically.
    Rps0 = reflectivity_ps(test_layers, 0.0, 2π * 1.0, ω_ref)
    @assert abs(Rps0) < 1e-8

    md"""
    !!! correct "Self-check"
        At normal incidence, `reflectivity_pp` matches the classical acoustic-impedance formula
        `` (Z_2-Z_1)/(Z_2+Z_1)=`` $(round(R_analytical, digits=4)) to
        $(round(err, sigdigits=2)) absolute error. And `reflectivity_ps(p=0)` =
        $(round(abs(Rps0), sigdigits=2)) — no mode conversion at normal incidence, exactly as
        physically required.
    """
end

# ╔═╡ 1d49ebac-04a4-44a5-890c-f565d34246c1
let
    # isolate the causal-Q dispersion correction specifically, not the much larger ordinary
    # frequency dependence any finite-thickness layer stack shows from pure interference
    # alone (confirmed directly: comparing Rpp at different frequencies to EACH OTHER, with
    # dispersion included throughout, swings by ~0.19 here -- dominated by interference
    # across the 10 km layer over a 50x frequency range, not by dispersion). Instead,
    # compare the SAME frequencies with a realistic Q=1000 against the same layers at
    # Q=Inf (lossless, non-dispersive) -- the difference isolates what causal_velocity's
    # own frequency dependence contributes, on its own.
    layers_Q = [Layer(10.0, 5.0, 3.0, 2.5, 1000.0, 1000.0), Layer(10.0, 4.0, 3.0, 2.5, 1000.0, 1000.0), Layer(Inf, 6.0, 3.5, 2.7, 1000.0, 1000.0)]
    layers_elastic = [Layer(10.0, 5.0, 3.0, 2.5, Inf, Inf), Layer(10.0, 4.0, 3.0, 2.5, Inf, Inf), Layer(Inf, 6.0, 3.5, 2.7, Inf, Inf)]
    p_test = 0.1
    freq_range = [0.1, 0.5, 1.0, 2.0, 5.0]
    Rpp_Q = [reflectivity_pp(layers_Q, p_test, 2π * f, ω_ref) for f in freq_range]
    Rpp_elastic = [reflectivity_pp(layers_elastic, p_test, 2π * f, ω_ref) for f in freq_range]
    @assert all(abs.(Rpp_Q) .<= 1.0 + 1e-8)
    dispersion_effect = maximum(abs.(Rpp_Q .- Rpp_elastic))
    @assert dispersion_effect < 0.08

    md"""
    !!! correct "Self-check"
        Across the $(round(freq_range[1], digits=2)) Hz to $(round(freq_range[end], digits=2)) Hz
        sweep, a realistic `` Q=1000`` moves `` |R_{pp}|`` by at most
        $(round(dispersion_effect, sigdigits=2)) relative to the same layers with infinite Q
        (lossless) — a genuinely small, isolated dispersion effect, once separated from the
        much larger ordinary frequency dependence the same layer stack shows from interference
        alone.
    """
end

# ╔═╡ f0a4c439-ee9b-4001-97db-66f6ac5afd5a
md"""## Sommerfeld Integral

Three phase types, each a dispatch on one of the singleton types below, selecting which term of
the point-source response [`get_wavefield`](@ref) sums over `` k ``: the **direct** wave
(straight from source to receiver, no reflection at all), the **reflected P** wave (weighted by
[`reflectivity_pp`](@ref), the layered stack's own P→P response), and the **reflected,
mode-converted P→S** wave (weighted by [`reflectivity_ps`](@ref) instead). All three reuse the
exact same [`conical_wave`](@ref) kernel for the source's own layer — only the reflectivity
factor differs.
"""

# ╔═╡ 0e8564a1-d300-4ca1-84f3-93ebea6ad08a
struct Direct end

# ╔═╡ e1908cc5-a4c1-4802-9f23-d9905a68cc36
struct Preflect end

# ╔═╡ 281a5cac-4e00-42d6-b331-ea5d80888f27
struct Psreflect end

# ╔═╡ 089b0602-077a-4074-a424-0a6afab8bcf4
"""
    smooth_taper(p, p_max; width=0.3)

Smooth cosine (Tukey-style) window that tapers the wavenumber integrand to zero near the upper
cutoff `p_max` instead of truncating it sharply: full amplitude for `p ≤ p_max(1-width)`, then a
raised-cosine roll-off to exactly `0` at `p_max`. Truncating a Fourier-type integral sharply
rings (Gibbs phenomenon); this taper is what keeps [`get_wavefield`](@ref)'s finite `kmax` cutoff
from showing up as spurious oscillation in the resulting seismogram.
"""
function smooth_taper(p, p_max; width=0.3)
    if p <= p_max * (1 - width)
        return 1.0
    elseif p < p_max
        t = π * (p - p_max * (1 - width)) / (width * p_max)
        return 0.5 * (1.0 + cos(t))
    else
        return 0.0
    end
end

# ╔═╡ 951a5946-ab54-4aff-be89-b8d5ea90fa1a
begin
    """
        get_phase(phase, k, kmax, param, layers::Vector{Layer}; pmax=Inf)

    One wavenumber `k`'s contribution to the point-source response, for whichever `phase` is
    selected ([`Direct`](@ref), [`Preflect`](@ref), or [`Psreflect`](@ref)): the tapered direct
    conical wave ([`conical_wave`](@ref), always evaluated in `layers[1]`, the source's own
    layer), times `1` (Direct), `reflectivity_pp` (Preflect), or `reflectivity_ps` (Psreflect).
    [`get_wavefield`](@ref) calls this once per `k` and sums the results.

    `pmax` (default `Inf`, no effect) lets a caller exclude every plane wave with horizontal
    slowness `p = k/ω` above it — reusing [`smooth_taper`](@ref) as a smooth roll-off at `pmax`
    itself, applied to `p` this time instead of `k`. Note `p` depends on *both* `k` and `ω`, so
    this is a per-`(k,ω)` mask, not a single scalar cutoff on `k` alone: the same `k` is
    post-critical at low frequency and sub-critical at high frequency. Only `Preflect`/
    `Psreflect` honor it — the direct wave has no reflection coefficient, hence no critical
    angle or head wave of its own to switch off.
    """
    function get_phase(::Direct, k, kmax, param, layers::Vector{Layer}; pmax=Inf)
        taper = smooth_taper(k, kmax; width=0.2)
        return taper .* conical_wave(k, param.Ω, param.R, param.Z, layers[1], ω_ref)
    end
    function get_phase(::Preflect, k, kmax, param, layers::Vector{Layer}; pmax=Inf)
        taper = smooth_taper(k, kmax; width=0.1)
        A = reflectivity_pp.(Ref(layers), k, param.ωgrid, Ref(ω_ref))
        if isfinite(pmax)
            A = A .* smooth_taper.(k ./ param.ωgrid, pmax; width=0.05)
        end
        C = conical_wave(k, param.Ω, param.R, param.Z, layers[1], ω_ref)
        return taper .* A .* C
    end
    function get_phase(::Psreflect, k, kmax, param, layers::Vector{Layer}; pmax=Inf)
        taper = smooth_taper(k, kmax; width=0.1)
        A = reflectivity_ps.(Ref(layers), k, param.ωgrid, Ref(ω_ref))
        if isfinite(pmax)
            A = A .* smooth_taper.(k ./ param.ωgrid, pmax; width=0.05)
        end
        C = conical_wave(k, param.Ω, param.R, param.Z, layers[1], ω_ref)
        return taper .* A .* C
    end
end

# ╔═╡ c246fdfc-321c-4856-8a4b-cf6cb4ba1594
"""
    get_wavefield(param, layers, phase; np=1024, pmax=Inf)

The Sommerfeld/Weyl integral itself, evaluated as a discrete Riemann sum over horizontal
wavenumber `k ∈ [0, kmax]` (`kmax` set generously above the fastest wave any receiver in
`param` could need, via the slowest S velocity present) — not `QuadGK` or any adaptive
quadrature, just `np` evenly-spaced samples weighted by [`get_phase`](@ref) and
[`smooth_taper`](@ref)'s roll-off. Returns the frequency-domain response on `param`'s own
`(ω, r)` grid, already windowed by the source spectrum `param.W`; [`remove_zero_frequency!`](@ref)
plus an inverse real FFT turns this into the time-domain seismogram.

`pmax` is forwarded to [`get_phase`](@ref) unchanged — see [`critical_slowness`](@ref) for the
physically meaningful values (the first interface's own critical slownesses) to pass here.
"""
function get_wavefield(param, layers::Vector{Layer}, phase; np=1024, pmax=Inf)
    vmax = maximum(l.vs for l in layers)
    kmax = maximum(param.ωgrid) / vmax * 1.2
    kgrid = range(0, kmax, length=np)
    dk = step(kgrid)

    integral_sum = zeros(ComplexF64, length(param.ωgrid), length(param.rgrid))
    for k in kgrid
        integral_sum .+= get_phase(phase, k, kmax, param, layers; pmax) .* dk
    end

    return param.W .* integral_sum
end

# ╔═╡ ef000001-0000-4000-8000-ef0000000001
"""
    get_wavefield_pp_ps(param, layers::Vector{Layer}; np=1024, pmax=Inf)

`Preflect` and `Psreflect` together, sharing each `` (k,\\omega) ``'s [`reflectivity_matrix`](@ref)
evaluation between both instead of computing it twice. `reflectivity_pp`/`reflectivity_ps` are
just that same matrix's `[1,1]`/`[2,1]` entries, so [`get_wavefield`](@ref) called once per phase
(the natural thing to do for `Direct`, which never needs a full reflectivity matrix at all) runs
the whole Kennett recursion from scratch a second time for `Psreflect` when `Preflect` already
computed it — the single largest cost in "sum" mode (the widget's own default view), confirmed
directly: `` \\mathrm{reflectivity\\_pp} `` and `` \\mathrm{reflectivity\\_ps} `` cost the same
`` \\sim 16\\,\\mu s `` per call, and there are `np \\times \\mathrm{length}(\\text{ωgrid}) ``
calls of each per `get_wavefield` call. Returns `(preflect, psreflect)`, each numerically
identical to what `get_wavefield(param, layers, Preflect(); pmax)`/`Psreflect()` would give.
"""
function get_wavefield_pp_ps(param, layers::Vector{Layer}; np=1024, pmax=Inf)
    vmax = maximum(l.vs for l in layers)
    kmax = maximum(param.ωgrid) / vmax * 1.2
    kgrid = range(0, kmax, length=np)
    dk = step(kgrid)

    sum_pp = zeros(ComplexF64, length(param.ωgrid), length(param.rgrid))
    sum_ps = zeros(ComplexF64, length(param.ωgrid), length(param.rgrid))
    for k in kgrid
        taper = smooth_taper(k, kmax; width=0.1)
        M = reflectivity_matrix.(Ref(layers), k, param.ωgrid, Ref(ω_ref))
        App = getindex.(M, 1, 1)
        Aps = getindex.(M, 2, 1)
        if isfinite(pmax)
            ptaper = smooth_taper.(k ./ param.ωgrid, pmax; width=0.05)
            App = App .* ptaper
            Aps = Aps .* ptaper
        end
        C = conical_wave(k, param.Ω, param.R, param.Z, layers[1], ω_ref)
        sum_pp .+= taper .* App .* C .* dk
        sum_ps .+= taper .* Aps .* C .* dk
    end
    return param.W .* sum_pp, param.W .* sum_ps
end

# ╔═╡ 3c1a2b4e-0001-4000-8000-100000000001
# fixed source depth (km), shared between the physics below and the widget's own drawing so
# the picture and the computation always agree
const LP_SOURCE_DEPTH = 5.0

# ╔═╡ 12f7733a-edd6-481f-8166-ef967520b35a
seismograms_param = let
    Nt = 512
    freq_snapshot = 1.0
    tgrid = range(0, 100, length=Nt)
    ωgrid = collect(rfftfreq(Nt, inv(step(tgrid)))) * 2.0 * pi
    f0 = freq_snapshot
    ω0 = 2.0 * pi * f0
    rgrid = collect(range(-100, 100, length=100))
    Ω = reshape(ωgrid, :, 1)
    R = reshape(rgrid, 1, :)
    Z = fill(LP_SOURCE_DEPTH, size(R))  # receiver at the surface; source at LP_SOURCE_DEPTH -> |δz| = LP_SOURCE_DEPTH
    W = @. exp(-0.1(Ω - ω0)^2)
    (; Nt, Ω, R, Z, W, rgrid, tgrid, ωgrid, ω0, f0)
end;

# ╔═╡ 6923e13d-a4c5-45c9-a1b2-f0104932d709
"""
    remove_zero_frequency!(C)

Zero the DC (zero-frequency) row of a frequency-domain field `C` in place. The Sommerfeld
integral's own zero-frequency response is not physically meaningful here (no static
displacement field is being modeled) and left alone it would inject a constant offset into
every trace after the inverse FFT.
"""
function remove_zero_frequency!(C)
    C[1, :] .= 0.0
    return C
end

# ╔═╡ debbbbbb-0000-4000-8000-dbdbdbdbdbdb
# Performance self-check, two independent fixes found by measuring on this live kernel (never
# guess -- see @elapsed/@allocated/@code_warntype in `git log` for the actual before/after
# session). Deliberately NOT re-measured live on every widget interaction: doing so would mean
# re-running the full expensive computation (and, for a fair before/after, the OLD slow path
# too) on every phase/head-wave toggle -- exactly the kind of overhead this fix is about
# removing, so a self-check that "verifies" the fix by re-paying its cost every click would be
# self-defeating. These numbers were measured once, on LAMB_DEFAULT_LAYERS, and are static text.
#
# 1. vertical_slowness used to return a bare Float64 on its propagating branch and a
#    ComplexF64 on its evanescent one -- broadcasting it over a k/ω grid straddling the
#    critical slowness left every downstream caller (conical_wave, wavefield_components,
#    solve_interface, compound_matrix_psv, reflectivity_matrix) doing fully-boxed Any
#    arithmetic (conical_wave's own code_warntype body was ::ANY). Fixed by wrapping the real
#    branch in complex(...) so both branches share one concrete type: get_wavefield(...,
#    Direct()) went from 6.56s/7.58GB to 1.55s/0.89GB on this model (~4x faster, ~8x less
#    allocation) from that one function alone.
# 2. In "sum" mode (the widget's own default view), Preflect and Psreflect each called
#    get_wavefield separately, and reflectivity_pp/reflectivity_ps -- despite being the exact
#    same reflectivity_matrix recursion's [1,1] and [2,1] entries -- each triggered their own
#    full recursion from scratch: the same expensive calculation, computed twice.
#    get_wavefield_pp_ps computes it once per (k,ω) and extracts both, verified bit-for-bit
#    identical to the old two-separate-calls result (`max|Δ| = 0.0`, checked during
#    development), not just "close enough".
#
# Combined, full "sum" mode on this model went from ~28.3s / ~39.6GB (neither fix) to
# 14.2s / 19.5GB (fix 1 only) to 7.8s / 9.8GB (both fixes) -- a ~3.6x speedup, ~4x less
# allocation, with the answer itself unchanged.
md"""
!!! correct "Self-check: performance"
    Full "sum" mode (`LAMB_DEFAULT_LAYERS`, measured once during development): **~28.3 s /
    ~39.6 GB → 7.8 s / 9.8 GB** (~3.6x faster, ~4x less allocation) from the two fixes
    documented in this cell's own source comment — `vertical_slowness`'s type instability, and
    `Preflect`/`Psreflect` redundantly recomputing the same `reflectivity_matrix` recursion.
    Both `get_wavefield` code paths agree bit-for-bit (`max|Δ| = 0.0`); only speed changed.
"""

# ╔═╡ c695e4d3-c49d-4587-8adc-cdd4c001325e
md"## Appendix"

# ╔═╡ fbd2b939-c816-42dd-a65e-94abce3dc91a
default_plotly_template(:plotly_dark)

# ╔═╡ 5a000001-0000-4000-8000-300000000001
md"## The Layered-Medium Widget"

# ╔═╡ 5a000002-0000-4000-8000-300000000002
default_layers = GUTENBERG_MODEL;

# ╔═╡ 5a000003-0000-4000-8000-300000000003
begin
    """
        LayeredMediumInput(layers=<Gutenberg 5-layer default>; show_vp=true, zmax=300.0)

    A dark-canvas widget for building a layered Earth model by direct manipulation: drag a
    layer boundary up/down to resize it, drag a track's marker left/right to change that
    layer's Vp/Vs/ρ, click empty space to add a layer, or drag a boundary onto its neighbor
    to delete it. The bottom layer is always the half-space (no lower boundary). Copied and
    adapted from `Love-wave-dispersion-curves.jl` (`show_vp=true` here, since Lamb's problem
    is P-SV and needs Vp, unlike Love's SH-only case) — see `pluto-widget-style` for the full
    design writeup of this pattern.
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
            for i in 1:(n-1)
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
        for i in 1:(n-1)
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
        b = Float64[]
        d = 0.0
        for i in 1:(n-1)
            d += layers[i].thickness
            push!(b, d)
        end
        "{boundaries:[$(join(b, ","))],vp:[$(join([l.vp for l in layers], ","))],vs:[$(join([l.vs for l in layers], ","))],rho:[$(join([l.rho for l in layers], ","))]}"
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
          <div class="lm-titlebar-sub">drag a boundary line to resize a layer &middot; drag a track's marker to change Vp/Vs/&rho; &middot; click empty space to add a layer &middot; drag a boundary onto its neighbor to delete it</div>
        </div>

        <div class="lm-row">
          <div class="lm-panel" style="flex:3 1 500px">
            <div class="lm-panel-title">Layered Earth Model</div>
            <canvas id="lm-editor" width="620" height="380"></canvas>
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
        </div>

        <script>
        {
          if (window.__lmCleanup) { window.__lmCleanup() }
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

          function applyPreset(p){
            state.boundaries = p.boundaries.slice()
            state.vp = p.vp.slice(); state.vs = p.vs.slice(); state.rho = p.rho.slice()
            drawEditor(); publish()
          }
          root.querySelector('#lm-preset-gutenberg').addEventListener('click', ()=>applyPreset(PRESETS.gutenberg), lmSignal)
          root.querySelector('#lm-preset-crust2').addEventListener('click', ()=>applyPreset(PRESETS.crust2), lmSignal)
          root.querySelector('#lm-preset-uniform').addEventListener('click', ()=>applyPreset(PRESETS.uniform), lmSignal)

          drawEditor()
        }
        </script>
        </div>
        """)
    end

    const _lm_ready = true
end

# ╔═╡ beef0001-0000-4000-8000-000000000001
begin
    _lm_ready
    LAMB_DEFAULT_LAYERS  # bare reference: forces Pluto to track this bind cell's dependency on
                         # LAMB_DEFAULT_LAYERS, since it's otherwise buried as a nested argument
                         # inside the @bind expression below (see pluto-widget-style's
                         # bind-cell-ordering note)
    PlutoUI.WideCell(@bind lm LayeredMediumInput(LAMB_DEFAULT_LAYERS; show_vp=true, zmax=80.0); max_width=1400)
end

# ╔═╡ beef0002-0000-4000-8000-000000000002
layers = let
    b = Float64.(lm["boundaries"])
    vp = Float64.(lm["vp"])
    vs = Float64.(lm["vs"])
    rho = Float64.(lm["rho"])
    n = length(vs)
    thickness = [i < n ? (i == 1 ? b[1] : b[i] - b[i-1]) : Inf for i in 1:n]
    [Layer(thickness[i], vp[i], vs[i], rho[i], 1000.0, 1000.0) for i in 1:n]
end

# ╔═╡ 4a000001-0000-4000-8000-400000000001
md"## The Controls Widget"

# ╔═╡ 4a000002-0000-4000-8000-400000000002
begin
    """
        LambControlsInput(layers::Vector{Layer}; phase0="sum", head_wave0=true, rmax=100.0)

    A dark-canvas widget combining the source/receiver-array picture with the two controls
    that drive everything downstream: which phase(s) the record section shows, and whether
    post-critical (head-wave-generating) plane waves are included at all. A **star** marks the
    fixed source (range 0, depth `LP_SOURCE_DEPTH`); small **downward triangles**, evenly
    spaced, mark a fixed receiver array along the free surface — nothing here is draggable, the
    array is the whole point (a record section, not one receiver at a time). Layer bands are
    shaded by Vp. The bound value is two plain values (a phase-mode string, a boolean), never
    the dispatch structs or a `layers`-derived quantity directly, so this widget's own `@bind`
    cell needs no dependency workaround (see `pluto-widget-style`'s bind-cell-ordering note).
    """
    struct LambControlsInput
        layers::Vector{Layer}
        phase0::String
        head_wave0::Bool
        rmax::Float64
    end
    LambControlsInput(layers::Vector{Layer}; phase0::String="sum", head_wave0::Bool=true, rmax::Float64=100.0) =
        LambControlsInput(layers, phase0, head_wave0, rmax)

    Base.get(w::LambControlsInput) = Dict{String,Any}("phase" => w.phase0, "head_wave" => w.head_wave0)

    _lp_layers_js(layers::Vector{Layer}) = let
        n = length(layers)
        b = Float64[]
        d = 0.0
        for i in 1:(n-1)
            d += layers[i].thickness
            push!(b, d)
        end
        "{boundaries:[$(join(b, ","))],vp:[$(join([l.vp for l in layers], ","))]}"
    end

    function Base.show(io::IO, ::MIME"text/html", w::LambControlsInput)
        layers_js = _lp_layers_js(w.layers)
        write(io, """
        <div id="lpwidget">
        <style>
        #lpwidget{font-family:sans-serif;color:#d1d5db}
        #lpwidget .lp-titlebar{background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px;margin-bottom:10px;text-align:center}
        #lpwidget .lp-titlebar-headline{font-size:16px;font-weight:700;color:#f3f4f6}
        #lpwidget .lp-row{display:flex;justify-content:center;gap:8px;flex-wrap:wrap;margin-bottom:10px}
        #lpwidget .lp-mode-btn{border-radius:4px;border:1px solid #6b7280;background:#0b0b0b;color:#e5e7eb;padding:6px 10px;font-size:13px;cursor:pointer}
        #lpwidget .lp-mode-btn:hover{background:#1f2937}
        #lpwidget .lp-mode-btn.active{border-color:#4ade80;color:#4ade80}
        #lpwidget .lp-headwave-btn{border-radius:4px;border:1px solid #6b7280;background:#0b0b0b;color:#e5e7eb;padding:6px 10px;font-size:13px;cursor:pointer}
        #lpwidget .lp-headwave-btn.active{border-color:#f59e0b;color:#f59e0b}
        #lpwidget .lp-primary{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start;margin-bottom:14px}
        #lpwidget .lp-secondary{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start}
        #lpwidget .lp-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #lpwidget .lp-panel-title{font-size:14px;font-weight:700;color:#f3f4f6;margin-bottom:4px;text-align:center}
        #lpwidget canvas{display:block;width:100%;height:auto}
        #lpwidget .lp-caption{font-size:12px;color:#9ca3af;margin-top:4px;text-align:center}
        </style>

        <div class="lp-titlebar">
          <div class="lp-titlebar-headline">Choose what the panels below show.</div>
        </div>
        <div class="lp-row">
          <button class="lp-mode-btn" data-phase="sum" type="button">Sum (realistic)</button>
          <button class="lp-mode-btn" data-phase="direct" type="button">Direct only</button>
          <button class="lp-mode-btn" data-phase="preflect" type="button">Reflected P only</button>
          <button class="lp-mode-btn" data-phase="psreflect" type="button">Reflected P→S only</button>
        </div>
        <div class="lp-row">
          <button class="lp-headwave-btn" data-hw="on" type="button">Head wave: on</button>
          <button class="lp-headwave-btn" data-hw="off" type="button">Head wave: off</button>
        </div>

        <div class="lp-primary">
          <div>
            <div class="lp-panel-title">Source and Receiver Array</div>
            <div class="lp-panel"><canvas id="lp-geom"></canvas></div>
            <div class="lp-caption" id="lp-caption">Loading…</div>
          </div>
          <div>
            <div class="lp-panel-title" id="lp-record-title">Record Section</div>
            <div class="lp-panel"><canvas id="lp-record"></canvas></div>
            <div class="lp-caption" id="lp-record-caption">computing…</div>
          </div>
        </div>

        <div class="lp-secondary">
          <div>
            <div class="lp-panel-title" id="lp-trace-title">Trace at Farthest Receiver</div>
            <div class="lp-panel"><canvas id="lp-trace"></canvas></div>
          </div>
          <div>
            <div class="lp-panel-title">Reflectivity vs Ray Parameter — the Critical-Slowness Kink</div>
            <div class="lp-panel"><canvas id="lp-kink"></canvas></div>
            <div class="lp-caption">Im(R_pp) blue, Im(R_ps) red; dashed lines mark each critical slowness</div>
          </div>
        </div>

        <script>
        {
          if (window.__lpCleanup) { window.__lpCleanup() }
          const lpController = new AbortController()
          window.__lpCleanup = () => lpController.abort()
          const lpSignal = { signal: lpController.signal }

          const root = document.getElementById('lpwidget')
          const LAYERS = $layers_js
          const SOURCE_DEPTH = $(LP_SOURCE_DEPTH)
          const RMAX = $(w.rmax)
          const NRECEIVERS = 15
          const state = { phase: "$(w.phase0)", head_wave: $(w.head_wave0) }
          const DPR = window.devicePixelRatio || 1
          let pushed = null // set by the lp-update CustomEvent, {seismo,nt,nr,umax,ttime,rmax,trace,tracemax,rfar,pgrid,imPP,imPS,kinkymax,pcrit,pcritLabels,phase,headwave} from Julia

          function publish(){
            root.value = { phase: state.phase, head_wave: state.head_wave }
            root.dispatchEvent(new CustomEvent('input'))
          }
          root.value = { phase: state.phase, head_wave: state.head_wave }

          function hidpi(canvas, context, w, h){
            canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR)
            canvas.style.width = w+'px'; canvas.style.height = h+'px'
            context.setTransform(DPR,0,0,DPR,0,0)
          }

          const availW = Math.min(window.innerWidth*0.85, root.clientWidth || window.innerWidth*0.85, 1300)
          const GEOM_W = Math.max(260, Math.floor(availW*0.32)), GEOM_H = 200
          const REC_W = Math.max(360, Math.floor(availW - GEOM_W - 40)), REC_H = 300
          const SEC_W = Math.max(280, Math.floor((availW - 40)/2)), SEC_H = 220

          const cv = root.querySelector('#lp-geom')
          const ctx = cv.getContext('2d')
          hidpi(cv, ctx, GEOM_W, GEOM_H)
          const EM = { l: 16, r: 16, t: 20, b: 16 }

          const recCv = root.querySelector('#lp-record')
          const recCtx = recCv.getContext('2d')
          hidpi(recCv, recCtx, REC_W, REC_H)

          const traceCv = root.querySelector('#lp-trace')
          const traceCtx = traceCv.getContext('2d')
          hidpi(traceCv, traceCtx, SEC_W, SEC_H)

          const kinkCv = root.querySelector('#lp-kink')
          const kinkCtx = kinkCv.getContext('2d')
          hidpi(kinkCv, kinkCtx, SEC_W, SEC_H)

          function depthMax(){
            const lastB = LAYERS.boundaries.length ? LAYERS.boundaries[LAYERS.boundaries.length - 1] : 20
            return Math.max(40, lastB * 1.6, SOURCE_DEPTH * 2)
          }
          function rToX(r){ return EM.l + (r + RMAX) / (2 * RMAX) * (GEOM_W - EM.l - EM.r) }
          function depthToY(z){ return EM.t + (z / depthMax()) * (GEOM_H - EM.t - EM.b) }
          function vpColor(vp){
            const t = Math.max(0, Math.min(1, (vp - 2) / 11))
            return 'rgb(' + Math.round(30 + 40 * t) + ',' + Math.round(60 + 80 * t) + ',' + Math.round(120 + 120 * t) + ')'
          }

          // diverging colormap for the record-section heatmap -- blue=positive, red=negative,
          // same convention as Born-approximation.jl's velColor.
          function seismicColor(v, mx){
            const t = Math.max(-1, Math.min(1, v/mx))
            if(t >= 0) return [Math.round(255*(1-t)), Math.round(255*(1-t)), 255]
            const s = -t
            return [255, Math.round(255*(1-s)), Math.round(255*(1-s))]
          }

          function drawStarMarker(cx, cy, r, fill, stroke){
            ctx.beginPath()
            const spikes = 5, rOuter = r, rInner = r * 0.45
            for(let i=0; i<spikes*2; i++){
              const rad = i % 2 === 0 ? rOuter : rInner
              const ang = -Math.PI/2 + i*Math.PI/spikes
              const x = cx + rad*Math.cos(ang), y = cy + rad*Math.sin(ang)
              i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y)
            }
            ctx.closePath()
            ctx.fillStyle = fill; ctx.fill()
            ctx.strokeStyle = stroke; ctx.lineWidth = 1; ctx.stroke()
          }
          function drawTriangleDownMarker(cx, cy, r, fill, stroke){
            ctx.beginPath()
            for(let i=0; i<3; i++){
              const ang = Math.PI/2 + i*2*Math.PI/3
              const x = cx + r*Math.cos(ang), y = cy + r*Math.sin(ang)
              i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y)
            }
            ctx.closePath()
            ctx.fillStyle = fill; ctx.fill()
            ctx.strokeStyle = stroke; ctx.lineWidth = 1.5; ctx.stroke()
          }

          // shared tick-drawing for any x/y-linear canvas panel below -- xmin/xmax (or
          // ymin/ymax) define the data range, ticks is an explicit array of data values to
          // label (avoids float step edge cases), matching the axis-band convention from
          // Born-approximation.jl's drawXAxis/drawYAxis but generalized to an arbitrary
          // (possibly negative) min instead of always starting at 0.
          function drawTicksX(c, W, H, xmin, xmax, ticks, unit){
            const bandH = 16
            c.fillStyle = 'rgba(0,0,0,0.55)'; c.fillRect(0, H-bandH, W, bandH)
            c.strokeStyle = '#4b5563'; c.beginPath(); c.moveTo(0,H-bandH+0.5); c.lineTo(W,H-bandH+0.5); c.stroke()
            c.fillStyle = '#9ca3af'; c.font = '9px sans-serif'; c.textBaseline = 'top'
            for(const x of ticks){
              const px = (x-xmin)/(xmax-xmin)*W
              c.strokeStyle = '#6b7280'; c.beginPath(); c.moveTo(px,H-bandH); c.lineTo(px,H-bandH+4); c.stroke()
              const label = Math.round(x)+(unit||'')
              const tw = c.measureText(label).width
              const tx = Math.max(2, Math.min(W - tw - 2, px - tw/2))
              c.textAlign = 'left'; c.fillText(label, tx, H-bandH+5)
            }
          }
          function drawTicksY(c, W, H, ymin, ymax, ticks, unit, digits, bandW){
            bandW = bandW || 34
            c.fillStyle = 'rgba(0,0,0,0.55)'; c.fillRect(0,0,bandW,H)
            c.strokeStyle = '#4b5563'; c.beginPath(); c.moveTo(bandW+0.5,0); c.lineTo(bandW+0.5,H); c.stroke()
            c.fillStyle = '#9ca3af'; c.font = '9px sans-serif'; c.textAlign = 'left'
            for(const y of ticks){
              const py = (y-ymin)/(ymax-ymin)*H
              c.strokeStyle = '#6b7280'; c.beginPath(); c.moveTo(bandW-4,py); c.lineTo(bandW,py); c.stroke()
              const label = y.toFixed(digits||0)+(unit||'')
              const tw = c.measureText(label).width
              const tx = Math.max(2, bandW-4-tw)
              const ty = Math.max(1, Math.min(H-1, py))
              c.textBaseline = py<=1 ? 'top' : (py>=H-1 ? 'bottom' : 'middle')
              c.fillText(label, tx, ty)
            }
          }
          function linTicks(lo, hi, n){
            const t = []
            for(let i=0;i<=n;i++) t.push(lo + (hi-lo)*i/n)
            return t
          }

          function draw(){
            ctx.clearRect(0,0,GEOM_W,GEOM_H)
            ctx.fillStyle = '#000'; ctx.fillRect(0,0,GEOM_W,GEOM_H)
            const n = LAYERS.vp.length
            for(let i=0;i<n;i++){
              const top = i===0 ? 0 : LAYERS.boundaries[i-1]
              const bot = i===n-1 ? depthMax() : LAYERS.boundaries[i]
              ctx.fillStyle = vpColor(LAYERS.vp[i])
              ctx.fillRect(EM.l, depthToY(top), GEOM_W-EM.l-EM.r, Math.max(1,depthToY(bot)-depthToY(top)))
            }
            ctx.strokeStyle = '#1f2937'; ctx.lineWidth = 1
            for(let i=0;i<n-1;i++){
              const y = depthToY(LAYERS.boundaries[i])
              ctx.beginPath(); ctx.moveTo(EM.l,y); ctx.lineTo(GEOM_W-EM.r,y); ctx.stroke()
            }
            ctx.strokeStyle = '#374151'; ctx.strokeRect(EM.l, EM.t, GEOM_W-EM.l-EM.r, GEOM_H-EM.t-EM.b)

            drawStarMarker(rToX(0), depthToY(SOURCE_DEPTH), 9, '#facc15', '#000')
            for(let i=0;i<NRECEIVERS;i++){
              const r = -RMAX + (2*RMAX)*i/(NRECEIVERS-1)
              drawTriangleDownMarker(rToX(r), depthToY(0)+2, 5, '#f5f3ef', '#0a0f18')
            }

            document.getElementById('lp-caption').textContent =
              NRECEIVERS + ' receivers, ' + (-RMAX).toFixed(0) + ' to ' + RMAX.toFixed(0) +
              ' km · source depth ' + SOURCE_DEPTH.toFixed(0) + ' km (fixed)'
          }

          function drawRecordSection(){
            recCtx.clearRect(0,0,REC_W,REC_H)
            const hwLabel = state.head_wave ? 'head wave on' : 'head wave off'
            document.getElementById('lp-record-title').textContent = 'Record Section — ' + state.phase + ' (' + hwLabel + ')'
            if(!pushed){
              recCtx.strokeStyle='#374151'; recCtx.lineWidth=1; recCtx.strokeRect(0.5,0.5,REC_W-1,REC_H-1)
              recCtx.fillStyle='#6b7280'; recCtx.font='12px sans-serif'; recCtx.fillText('computing...',10,18)
              return
            }
            const NT = pushed.nt, NR = pushed.nr, U = pushed.seismo, mx = pushed.umax
            const bandL = 34, bandB = 16
            const plotW = REC_W - bandL, plotH = REC_H - bandB
            const img = recCtx.createImageData(Math.round(plotW*DPR), Math.round(plotH*DPR))
            const wpx = img.width, hpx = img.height
            for(let py=0; py<hpx; py++){
              const it = Math.min(NT-1, Math.round(py/hpx*(NT-1)))
              for(let px2=0; px2<wpx; px2++){
                const ir = Math.min(NR-1, Math.round(px2/wpx*(NR-1)))
                const v = U[it + ir*NT]
                const [r,g,b] = seismicColor(v, mx)
                const idx = (py*wpx+px2)*4
                img.data[idx]=r; img.data[idx+1]=g; img.data[idx+2]=b; img.data[idx+3]=255
              }
            }
            recCtx.putImageData(img, Math.round(bandL*DPR), 0)
            recCtx.strokeStyle = '#374151'; recCtx.lineWidth = 1; recCtx.strokeRect(bandL+0.5, 0.5, plotW-1, plotH-1)
            drawTicksX(recCtx, REC_W, REC_H, -pushed.rmax, pushed.rmax, linTicks(-pushed.rmax, pushed.rmax, 4), '')
            drawTicksY(recCtx, REC_W, REC_H, 0, pushed.ttime[pushed.ttime.length-1], linTicks(0, pushed.ttime[pushed.ttime.length-1], 5), '', 0, bandL)
            document.getElementById('lp-record-caption').textContent =
              'distance to receiver (km) vs. time (s) · ' + NR + ' receivers'
          }

          function drawTrace(){
            traceCtx.clearRect(0,0,SEC_W,SEC_H)
            if(!pushed){
              traceCtx.strokeStyle='#374151'; traceCtx.lineWidth=1; traceCtx.strokeRect(0.5,0.5,SEC_W-1,SEC_H-1)
              traceCtx.fillStyle='#6b7280'; traceCtx.font='12px sans-serif'; traceCtx.fillText('computing...',10,18)
              return
            }
            document.getElementById('lp-trace-title').textContent =
              'Trace at Farthest Receiver (r = ' + pushed.rfar.toFixed(1) + ' km)'
            const t = pushed.ttime, y = pushed.trace, ymax = Math.max(pushed.tracemax, 1e-12)
            const bandL = 34, bandB = 16
            const plotW = SEC_W - bandL, plotH = SEC_H - bandB
            const tmax = t[t.length-1]
            function toPx(ti, yi){ return [bandL + ti/tmax*plotW, plotH/2 - (yi/ymax)*(plotH/2*0.9)] }
            traceCtx.strokeStyle = '#374151'; traceCtx.beginPath(); traceCtx.moveTo(bandL, plotH/2); traceCtx.lineTo(SEC_W, plotH/2); traceCtx.stroke()
            traceCtx.strokeStyle = '#38bdf8'; traceCtx.lineWidth = 1.2; traceCtx.beginPath()
            for(let i=0;i<t.length;i++){ const [px,py]=toPx(t[i],y[i]); i===0?traceCtx.moveTo(px,py):traceCtx.lineTo(px,py) }
            traceCtx.stroke()
            traceCtx.strokeStyle = '#374151'; traceCtx.strokeRect(bandL+0.5, 0.5, plotW-1, plotH-1)
            drawTicksX(traceCtx, SEC_W, SEC_H, 0, tmax, linTicks(0, tmax, 4), 's')
          }

          function drawKink(){
            kinkCtx.clearRect(0,0,SEC_W,SEC_H)
            if(!pushed){
              kinkCtx.strokeStyle='#374151'; kinkCtx.lineWidth=1; kinkCtx.strokeRect(0.5,0.5,SEC_W-1,SEC_H-1)
              kinkCtx.fillStyle='#6b7280'; kinkCtx.font='12px sans-serif'; kinkCtx.fillText('computing...',10,18)
              return
            }
            const p = pushed.pgrid, imPP = pushed.imPP, imPS = pushed.imPS, ymax = Math.max(pushed.kinkymax, 1e-12)
            const bandL = 34, bandB = 16
            const plotW = SEC_W - bandL, plotH = SEC_H - bandB
            const pmax = p[p.length-1]
            function toPx(pi, yi){ return [bandL + pi/pmax*plotW, plotH/2 - (yi/ymax)*(plotH/2*0.9)] }
            function drawLine(arr, color){
              kinkCtx.strokeStyle = color; kinkCtx.lineWidth = 1.3; kinkCtx.beginPath()
              for(let i=0;i<p.length;i++){ const [px,py]=toPx(p[i],arr[i]); i===0?kinkCtx.moveTo(px,py):kinkCtx.lineTo(px,py) }
              kinkCtx.stroke()
            }
            kinkCtx.strokeStyle = '#374151'; kinkCtx.beginPath(); kinkCtx.moveTo(bandL, plotH/2); kinkCtx.lineTo(SEC_W, plotH/2); kinkCtx.stroke()
            drawLine(imPP, '#3b82f6')
            drawLine(imPS, '#ef4444')
            kinkCtx.setLineDash([4,3]); kinkCtx.strokeStyle = '#9ca3af'
            for(const pc of pushed.pcrit){
              const [px] = toPx(pc, 0)
              kinkCtx.beginPath(); kinkCtx.moveTo(px, 0); kinkCtx.lineTo(px, plotH); kinkCtx.stroke()
            }
            kinkCtx.setLineDash([])
            kinkCtx.fillStyle = '#9ca3af'; kinkCtx.font = '10px sans-serif'; kinkCtx.textAlign = 'center'
            for(let i=0;i<pushed.pcrit.length;i++){
              const [px] = toPx(pushed.pcrit[i], 0)
              kinkCtx.fillText(pushed.pcritLabels[i], px, 10)
            }
            kinkCtx.strokeStyle = '#374151'; kinkCtx.strokeRect(bandL+0.5, 0.5, plotW-1, plotH-1)
            drawTicksX(kinkCtx, SEC_W, SEC_H, 0, pmax, linTicks(0, pmax, 4), '')
          }

          function setButtons(){
            root.querySelectorAll('.lp-mode-btn').forEach(b => b.classList.toggle('active', b.dataset.phase === state.phase))
            root.querySelectorAll('.lp-headwave-btn').forEach(b => b.classList.toggle('active', (b.dataset.hw==='on') === state.head_wave))
          }
          root.querySelectorAll('.lp-mode-btn').forEach(b => {
            b.addEventListener('click', () => { state.phase = b.dataset.phase; setButtons(); publish(); drawRecordSection() }, lpSignal)
          })
          root.querySelectorAll('.lp-headwave-btn').forEach(b => {
            b.addEventListener('click', () => { state.head_wave = (b.dataset.hw==='on'); setButtons(); publish(); drawRecordSection() }, lpSignal)
          })
          root.addEventListener('lp-update', e => {
            pushed = e.detail
            drawRecordSection(); drawTrace(); drawKink()
          }, lpSignal)

          setButtons()
          draw()
          drawRecordSection(); drawTrace(); drawKink()
        }
        </script>
        </div>
        """)
    end

    const _lp_ready = true
end

# ╔═╡ beef0003-0000-4000-8000-000000000003
begin
    _lp_ready
    layers  # bare reference: forces Pluto to track this bind cell's dependency on `layers`,
            # which a bind expression's own nested-argument reference alone would not (see
            # pluto-widget-style's bind-cell-ordering note)
    PlutoUI.WideCell(@bind ctrl LambControlsInput(layers); max_width=1400)
end

# ╔═╡ e7c1cbf4-7939-45a2-a754-b36ae9bdcab4
md"""### To Do
- **A free surface.** Everything here is reflection off interfaces *below* the source; there's
  no free surface above it yet, so no Rayleigh wave. Surface-wave dispersion already has a
  dedicated home in this repo — see `Rayleigh-dispersion-curves.jl` and `Rayleigh-function.jl`
  under Surface Waves — so this notebook deliberately doesn't duplicate that; a free surface
  here would mainly be about how it changes the *reflectivity* seen from below, not a second
  dispersion-curve panel.
- **Remaining performance headroom.** After fixing `vertical_slowness`'s type instability and
  the redundant `Preflect`/`Psreflect` recomputation (see the Performance self-check above,
  ~3.6x combined speedup), `reflectivity_pp`/`reflectivity_ps`'s own per-call cost (~16.7 μs,
  ~32 KB) didn't move — `solve_interface` builds several small `Matrix{ComplexF64}`s and calls
  `inv(...)`/`\\` per interface, per `(k,ω)` pair, and this is called `np × length(ωgrid)` times
  per `get_wavefield` call. Preallocating buffers across calls (in-place `A \\ b`, avoiding
  repeated `zeros(ComplexF64,4,4)`) is the next concrete place to look, not yet attempted.
- Multiply-refracted head waves from deeper interfaces (only the first interface's is shown).
- Plot the displacement field itself, not just potentials.
"""

# ╔═╡ 9d10ab94-38d5-406e-8697-7f5bc0c11df2
md"""### References
- Fuchs, K., and Gerhard Müller. "Computation of synthetic seismograms with the reflectivity method and comparison with observations." Geophysical Journal International 23.4 (1971): 417-433.
- Aki, K., and Richards, P.G. *Quantitative Seismology*, 2nd ed. University Science Books, 2002 — the reflectivity method (ch. 9) and head waves/critical refraction (ch. 6) this notebook follows.
- Lamb, H. "On the propagation of tremors over the surface of an elastic solid." *Philosophical Transactions of the Royal Society A* 203 (1904): 1-42 — the original problem this notebook is named after.
- `Rayleigh-dispersion-curves.jl`, `Rayleigh-function.jl` (this repo, Surface Waves) — the surface-wave half of the classical Lamb's-problem story, covered there rather than duplicated here.
"""

# ╔═╡ 6c000001-0000-4000-8000-600000000001
"""
    critical_slowness(layers::Vector{Layer})

The two critical horizontal slownesses at the first interface (`layers[1]` into `layers[2]`):
`p_c_P = 1/layers[2].vp` for the transmitted P wave, `p_c_S = 1/layers[2].vs` for the
transmitted S wave. Below either value the corresponding transmitted wave in `layers[2]`
still propagates and the interface behaves like an ordinary Zoeppritz boundary; at and beyond
it, [`vertical_slowness`](@ref) makes that transmission evanescent — total internal
reflection — which is exactly the range of plane waves [`get_wavefield`](@ref)'s
`pmax` argument can exclude to switch the corresponding head wave off.
"""
function critical_slowness(layers::Vector{Layer})
    length(layers) < 2 && return (p_c_P=Inf, p_c_S=Inf)
    l2 = layers[2]
    return (p_c_P=1 / l2.vp, p_c_S=1 / l2.vs)
end

# ╔═╡ 3c1a2b4e-0002-4000-8000-100000000002
# pmax=Inf includes every plane wave (head wave present); pmax=p_c_P excludes everything at or
# beyond the first interface's own P critical slowness (head wave switched off) -- see
# critical_slowness and get_phase/get_wavefield's own docstrings for why this specific value.
lamb_pmax = ctrl["head_wave"] ? Inf : critical_slowness(layers).p_c_P

# ╔═╡ 53a45db3-8571-4fc7-a855-5173030b7cb9
seismograms = let
    phase = ctrl["phase"]
    C = if phase == "direct"
        get_wavefield(seismograms_param, layers, Direct(); pmax=lamb_pmax)
    elseif phase == "preflect"
        get_wavefield(seismograms_param, layers, Preflect(); pmax=lamb_pmax)
    elseif phase == "psreflect"
        get_wavefield(seismograms_param, layers, Psreflect(); pmax=lamb_pmax)
    else # "sum" -- what a real seismogram actually records: every phase, superposed. Preflect
        # and Psreflect are computed TOGETHER here (get_wavefield_pp_ps), not as two separate
        # get_wavefield calls, since they'd otherwise redo the same expensive reflectivity_matrix
        # recursion twice for no reason -- see that function's own docstring.
        pp, ps = get_wavefield_pp_ps(seismograms_param, layers; pmax=lamb_pmax)
        get_wavefield(seismograms_param, layers, Direct(); pmax=lamb_pmax) .+ pp .+ ps
    end
    irfft(remove_zero_frequency!(C), seismograms_param.Nt, 1)
end;

# ╔═╡ 4a000003-0000-4000-8000-400000000003
begin
    _lp_flatten(v) = join(v, ",")

    """
        LambPush(...)

    `BsPush` for `LambControlsInput` (see `Born-approximation.jl`) -- does no physics itself, just
    takes the already-computed seismogram/farthest-trace/`` R(p)`` arrays and hands them to the
    widget's canvases via a `lp-update` `CustomEvent`, so the record section, trace, and kink
    panels live inside the SAME hero widget the reader already sees, instead of as separate
    downstream plot cells. All fields are pre-flattened/pre-formatted strings ready to splice
    into a JS array or object literal -- see [`_lp_flatten`](@ref).
    """
    struct LambPush
        seismo::String
        nt::Int
        nr::Int
        umax::Float64
        ttime::String
        rmax::Float64
        trace::String
        tracemax::Float64
        rfar::Float64
        pgrid::String
        imPP::String
        imPS::String
        kinkymax::Float64
        pcrit::String
        pcritLabels::String
    end

    function Base.show(io::IO, ::MIME"text/html", p::LambPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('lpwidget');
        if(w){
          w.dispatchEvent(new CustomEvent('lp-update', { detail: {
            seismo: [$(p.seismo)],
            nt: $(p.nt),
            nr: $(p.nr),
            umax: $(p.umax),
            ttime: [$(p.ttime)],
            rmax: $(p.rmax),
            trace: [$(p.trace)],
            tracemax: $(p.tracemax),
            rfar: $(p.rfar),
            pgrid: [$(p.pgrid)],
            imPP: [$(p.imPP)],
            imPS: [$(p.imPS)],
            kinkymax: $(p.kinkymax),
            pcrit: [$(p.pcrit)],
            pcritLabels: $(p.pcritLabels),
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 96556fec-4bea-4132-946a-92d17372df54
# computes the record-section/trace/R(p)-kink data and pushes it into LambControlsInput's own
# canvases (see LambPush, in the Appendix) -- this is what replaced three separate PlutoPlotly
# cells, so the reader sees everything inside the one hero widget instead of scattered below it.
let
    U = seismograms
    j = length(seismograms_param.rgrid)

    pc = critical_slowness(layers)
    crit_ps = unique(filter(isfinite, [pc.p_c_P, pc.p_c_S]))
    p_max_plot = isempty(crit_ps) ? 0.3 : 1.3 * maximum(crit_ps)
    p_range_kink = collect(range(1e-4, p_max_plot, length=400))
    ω_rp = seismograms_param.ω0
    Im_pp = [imag(reflectivity_pp(layers, p, ω_rp, ω_ref)) for p in p_range_kink]
    Im_ps = [imag(reflectivity_ps(layers, p, ω_rp, ω_ref)) for p in p_range_kink]
    pcrit_labels = "[" * join(["\"" * (pcrit == pc.p_c_P ? "p_c^P" : "p_c^S") * "\"" for pcrit in crit_ps], ",") * "]"

    LambPush(
        _lp_flatten(vec(U)), seismograms_param.Nt, length(seismograms_param.rgrid),
        maximum(abs, U) * 0.1,
        _lp_flatten(seismograms_param.tgrid), seismograms_param.rgrid[j],
        _lp_flatten(U[:, j]), maximum(abs, U[:, j]), seismograms_param.rgrid[j],
        _lp_flatten(p_range_kink), _lp_flatten(Im_pp), _lp_flatten(Im_ps),
        max(maximum(abs, Im_pp), maximum(abs, Im_ps), 1e-6),
        _lp_flatten(crit_ps), pcrit_labels,
    )
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bessels = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"

[compat]
Bessels = "~0.2.8"
FFTW = "~1.10.0"
HypertextLiteral = "~1.0.0"
PlutoPlotly = "~0.6.6"
PlutoUI = "~0.7.83"
Roots = "~3.0.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "1576736fbf9894f318d08e689061d2ed32178379"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
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
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
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

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

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

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

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

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

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

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

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

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c7345ab1a7ca4dc8a02c9f6510da0d9857bbe513"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.7.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

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

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

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

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "3de8f5e6e90ebfa8d6d1f86997d6cdcd6a912ff3"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.7"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "6256ab3ee24ef079b3afa310593817e069925eeb"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.23"

    [deps.PlotlyBase.extensions]
    DataFramesExt = "DataFrames"
    DistributionsExt = "Distributions"
    IJuliaExt = "IJulia"
    JSON3Ext = "JSON3"

    [deps.PlotlyBase.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "ColorSchemes", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "PrecompileTools", "Reexport", "ScopedValues", "Scratch", "TOML"]
git-tree-sha1 = "2b9e3d771adfe535a4fdda855f4741fdaacd3f7f"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.6"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "13a9e0164267bb9ebc55c1c5d72fceea8f09555b"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.7"

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

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "67a144433c4ce877ee6d1ada69a124d6b1ecf7be"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.6.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
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

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "2d0fc55c61321ba245c47be599570d11bac50303"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.5"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

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

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.0"

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
# ╠═897afffa-77e8-11ef-1a54-c73d1df8f6a4
# ╠═aa38e7f6-b8d2-4270-9142-1aa688041eb4
# ╠═4621f804-6a44-4c46-af37-3d364ba74cfe
# ╠═9b488729-aadc-4c1d-b971-6931d4bc9a08
# ╠═6188ff3e-cfd5-4c9c-aa39-619dc280d494
# ╟─1bbedd43-9e69-4ba0-a251-855f1f63dcf8
# ╠═beef0000-0000-4000-8000-000000000000
# ╠═beef0001-0000-4000-8000-000000000001
# ╠═beef0002-0000-4000-8000-000000000002
# ╠═beef0003-0000-4000-8000-000000000003
# ╟─72aca1bf-3997-4be3-b316-921b45480c1a
# ╟─d1f5e223-cc49-4fe8-9b2e-c1ccecf2315a
# ╠═614a3832-7bc7-4bd3-a09d-50b76afdd7cf
# ╟─7a3bd5df-0f9f-489c-b522-4098c325c0c0
# ╠═05ec38ea-3431-490c-bc38-24f8c1b2d54f
# ╟─e2fea7d5-00a3-4796-ad51-26f2dcffa55b
# ╠═5ea371a8-8035-4ba3-ab83-a9ffa2e3d504
# ╠═f047c190-0bf3-4763-9f31-b4272e837dd2
# ╟─360de109-d7db-4402-bca8-5b39c6f17da9
# ╠═f25f798f-0ecc-4931-b8b5-9c958b83850f
# ╠═687b7062-ee25-4498-948f-43b187e0ccfa
# ╠═8902f19a-6b93-458a-9198-f2322b2f5ab5
# ╠═7ca06983-11c1-4203-9bd0-0cc1a461b74f
# ╠═c9bf8e13-eee2-45ae-9ce5-4901c344c8f6
# ╠═e6bbca73-e59c-4b6a-854e-438222778110
# ╠═d60d4d65-8896-40ee-b52f-17d66999c727
# ╠═da2a2e26-f217-41e0-8485-6517af621d2a
# ╠═43c34152-fc5f-4497-882f-7f42bd4e6b99
# ╠═dea0645d-cb7c-4488-913b-ba225595aceb
# ╠═579cfd9f-e872-4ec5-b913-c90f8e183247
# ╠═a1e5796b-6f2f-44c6-b5bf-df6c9ee20960
# ╠═55e64939-9298-4a15-a2b0-0d13c23d03dc
# ╠═1d49ebac-04a4-44a5-890c-f565d34246c1
# ╟─f0a4c439-ee9b-4001-97db-66f6ac5afd5a
# ╠═0e8564a1-d300-4ca1-84f3-93ebea6ad08a
# ╠═e1908cc5-a4c1-4802-9f23-d9905a68cc36
# ╠═281a5cac-4e00-42d6-b331-ea5d80888f27
# ╠═089b0602-077a-4074-a424-0a6afab8bcf4
# ╠═951a5946-ab54-4aff-be89-b8d5ea90fa1a
# ╠═c246fdfc-321c-4856-8a4b-cf6cb4ba1594
# ╠═ef000001-0000-4000-8000-ef0000000001
# ╠═3c1a2b4e-0001-4000-8000-100000000001
# ╠═12f7733a-edd6-481f-8166-ef967520b35a
# ╠═3c1a2b4e-0002-4000-8000-100000000002
# ╠═6923e13d-a4c5-45c9-a1b2-f0104932d709
# ╠═53a45db3-8571-4fc7-a855-5173030b7cb9
# ╠═debbbbbb-0000-4000-8000-dbdbdbdbdbdb
# ╠═96556fec-4bea-4132-946a-92d17372df54
# ╟─c695e4d3-c49d-4587-8adc-cdd4c001325e
# ╠═fbd2b939-c816-42dd-a65e-94abce3dc91a
# ╟─5a000001-0000-4000-8000-300000000001
# ╠═5a000002-0000-4000-8000-300000000002
# ╠═5a000003-0000-4000-8000-300000000003
# ╟─4a000001-0000-4000-8000-400000000001
# ╠═4a000002-0000-4000-8000-400000000002
# ╟─e7c1cbf4-7939-45a2-a754-b36ae9bdcab4
# ╟─9d10ab94-38d5-406e-8697-7f5bc0c11df2
# ╠═6c000001-0000-4000-8000-600000000001
# ╠═4a000003-0000-4000-8000-400000000003
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
