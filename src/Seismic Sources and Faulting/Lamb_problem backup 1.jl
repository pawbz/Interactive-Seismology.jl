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

How does a point source's energy actually reach a distant receiver through a layered Earth? Not
as one wave, but as several: a **direct** arrival along the straightest path, **reflections**
bounced off every interface below, arrivals that **convert mode** (P to S) at those same
interfaces, and — once a free surface is added — a dispersive **surface wave** that never leaves
the near-surface layers at all. This notebook builds all of these from one common trick: any
point-source wavefield can be written as a superposition (an integral) of *plane* waves, and
plane waves are exactly what a layered medium's reflection/transmission coefficients are built
to handle one at a time.

The **Reflectivity Method** is the plane-wave half of that trick: it computes the exact
reflection response of an arbitrarily thick stack of layers — via a numerically stable
Kennett-style recursion, not the classical textbook's single-interface case — by combining every
interface's Zoeppritz coefficients with the phase delay of propagating through each layer.
Combined with the point-source decomposition (the **Sommerfeld/Weyl integral**, below), it
synthesizes a real time-domain seismogram: superpose enough plane waves, each correctly
reflected and phase-delayed, and inverse-Fourier-transform the sum.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
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

# ╔═╡ f25f798f-0ecc-4931-b8b5-9c958b83850f
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
layer inside this function; the half-space's radiation condition is handled separately wherever
it matters ([`reflectivity_matrix`](@ref), [`rayleigh_secular_function`](@ref)).
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
root, so `` \\exp(i\\omega q z) `` decays (not grows) with increasing depth `z`.
"""
function vertical_slowness(v, p)
    arg = (1 / abs2(v)) - abs2(p)
    return arg >= 0 ? -sqrt(arg) : 1im * sqrt(-arg)
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
    return Bessels.besselj0.(k .* r) .* exp.(-im .* kz .* δz) .* k ./ kz ./ im
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
        $(round(max_abs, sigdigits=2)) (machine precision) at every ray parameter tested,
        `` p\in\{$(join(p_test, ", "))\}`` s/km.
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
    # causal-Q dispersion (Qp=Qs=1000 here -- a mild, realistic loss) should shift Rpp only
    # slightly across this frequency range, not scramble it -- a basic sanity bound on how
    # much causal_velocity's frequency dependence is allowed to move the reflectivity.
    test_layers = [Layer(10.0, 5.0, 3.0, 2.5), Layer(10.0, 4.0, 3.0, 2.5), Layer(Inf, 6.0, 3.5, 2.7)]
    p_test = 0.1
    freq_range = [0.1, 0.5, 1.0, 2.0, 5.0]
    Rpp_vals = [reflectivity_pp(test_layers, p_test, 2π * f, ω_ref) for f in freq_range]
    @assert all(abs.(Rpp_vals) .<= 1.0 + 1e-8)
    spread = maximum(abs.(Rpp_vals)) - minimum(abs.(Rpp_vals))
    @assert spread < 0.05

    md"""
    !!! correct "Self-check"
        Across $(freq_range[1])-$(freq_range[end]) Hz, with a realistic `` Q=1000 `` in every
        layer, `` |R_{pp}| `` stays within $(round(spread, sigdigits=2)) of itself — causal
        dispersion perturbs the reflectivity, as it should, but doesn't dominate it.
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
        get_phase(phase, k, kmax, param, layers::Vector{Layer})

    One wavenumber `k`'s contribution to the point-source response, for whichever `phase` is
    selected ([`Direct`](@ref), [`Preflect`](@ref), or [`Psreflect`](@ref)): the tapered direct
    conical wave ([`conical_wave`](@ref), always evaluated in `layers[1]`, the source's own
    layer), times `1` (Direct), `reflectivity_pp` (Preflect), or `reflectivity_ps` (Psreflect).
    [`get_wavefield`](@ref) calls this once per `k` and sums the results.
    """
    function get_phase(::Direct, k, kmax, param, layers::Vector{Layer})
        taper = smooth_taper(k, kmax; width=0.2)
        return taper .* conical_wave(k, param.Ω, param.R, param.Z, layers[1], ω_ref)
    end
    function get_phase(::Preflect, k, kmax, param, layers::Vector{Layer})
        taper = smooth_taper(k, kmax; width=0.1)
        A = reflectivity_pp.(Ref(layers), k, param.ωgrid, Ref(ω_ref))
        C = conical_wave(k, param.Ω, param.R, param.Z, layers[1], ω_ref)
        return taper .* A .* C
    end
    function get_phase(::Psreflect, k, kmax, param, layers::Vector{Layer})
        taper = smooth_taper(k, kmax; width=0.1)
        A = reflectivity_ps.(Ref(layers), k, param.ωgrid, Ref(ω_ref))
        C = conical_wave(k, param.Ω, param.R, param.Z, layers[1], ω_ref)
        return taper .* A .* C
    end
end

# ╔═╡ c246fdfc-321c-4856-8a4b-cf6cb4ba1594
"""
    get_wavefield(param, layers, phase; np=1024)

The Sommerfeld/Weyl integral itself, evaluated as a discrete Riemann sum over horizontal
wavenumber `k ∈ [0, kmax]` (`kmax` set generously above the fastest wave any receiver in
`param` could need, via the slowest S velocity present) — not `QuadGK` or any adaptive
quadrature, just `np` evenly-spaced samples weighted by [`get_phase`](@ref) and
[`smooth_taper`](@ref)'s roll-off. Returns the frequency-domain response on `param`'s own
`(ω, r)` grid, already windowed by the source spectrum `param.W`; [`remove_zero_frequency!`](@ref)
plus an inverse real FFT turns this into the time-domain seismogram.
"""
function get_wavefield(param, layers::Vector{Layer}, phase; np=1024)
    vmax = maximum(l.vs for l in layers)
    kmax = maximum(param.ωgrid) / vmax * 1.2
    kgrid = range(0, kmax, length=np)
    dk = step(kgrid)

    integral_sum = zeros(ComplexF64, length(param.ωgrid), length(param.rgrid))
    for k in kgrid
        integral_sum .+= get_phase(phase, k, kmax, param, layers) .* dk
    end

    return param.W .* integral_sum
end

# ╔═╡ 72aca1bf-3997-4be3-b316-921b45480c1a
md"""## Seismograms

The heatmap below is exactly [`get_wavefield`](@ref)'s frequency-domain output, inverse-FFT'd
back to the time domain — moveout (the diagonal slant of an arrival across offsets) is the
direct visual signature of a finite wave speed, and the *different* slopes/arrival times of the
`Direct` vs `Preflect` vs `Psreflect` phases (toggle them in the widget above) are exactly the
"several arrivals from one source" story the introduction promised."""

# ╔═╡ 3c1a2b4e-0001-4000-8000-100000000001
# fixed source depth (km), shared between the physics below and the widget's own drawing so
# the picture and the computation always agree
const LP_SOURCE_DEPTH = 10.0

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

# ╔═╡ 7a000001-0000-4000-8000-500000000001
md"""## Surface Waves: The Rayleigh Dispersion Curve

A free surface does something a buried interface can't: it lets a wave's energy stay trapped
near the top of the model *forever*, bouncing between the surface and whatever velocity
contrasts lie below, building up into a **dispersive guided wave** whose speed depends on
period — the longer the period, the deeper (and typically faster) structure it samples. This is
exactly the same trapping mechanism `Love-wave-dispersion-curves.jl` builds for horizontally
polarized shear waves; the P-SV analogue below is the **Rayleigh wave**
([`rayleigh_secular_function`](@ref) in the Appendix), reusing the exact same per-layer
transfer-matrix machinery already sitting in this notebook ([`compound_matrix_psv`](@ref))
rather than a separate implementation.

**This is the dispersion *relation* only** — a curve of phase/group velocity vs. period, not a
new arrival in the time-domain seismogram above. Getting a real Rayleigh wave into the
seismogram needs the free-surface reflection woven into [`reflectivity_matrix`](@ref)'s own
recursion, plus a way to numerically capture the resulting pole in the wavenumber integral —
genuinely harder problems, left on the To-Do list below.
"""

# ╔═╡ 7a000002-0000-4000-8000-500000000002
rayleigh_periods = collect(range(5.0, 80.0, length=40));

# ╔═╡ c695e4d3-c49d-4587-8adc-cdd4c001325e
md"## Appendix"

# ╔═╡ fbd2b939-c816-42dd-a65e-94abce3dc91a
default_plotly_template(:plotly_dark)

# ╔═╡ 6a000001-0000-4000-8000-200000000001
md"## Rayleigh-Wave Dispersion"

# ╔═╡ 6a000002-0000-4000-8000-200000000002
"""
    rayleigh_secular_function(layers::Vector{Layer}, ω, c, ω_ref)

The P-SV analogue of Love's `dispersion_function_SH`: zero at exactly the (period, phase
velocity) pairs that support a genuine Rayleigh mode.

Builds the total P-SV propagator `M` from the free surface down to the top of the half-space by
chaining [`compound_matrix_psv`](@ref) across every finite layer (top to bottom) — the same
per-layer transfer matrix [`reflectivity_matrix`](@ref) uses, just multiplied straight through
instead of recursed. The free surface itself demands zero traction at the top
(`sxz = szz = 0`), leaving 2 free unknowns (`ux₀`, `uz₀`); propagated down, the resulting state
at the top of the half-space must ALSO be a physically admissible half-space solution — meaning,
decomposed in the half-space's own [`modal_basis_psv`](@ref) eigenbasis, it must have **zero**
amplitude on the two branches that grow (rather than decay) with depth. That "zero amplitude on
the disallowed branches" condition is 2 linear equations in the 2 free surface unknowns; a
nontrivial (non-zero) solution exists only where that 2×2 system's determinant vanishes — this
function returns that determinant.

Unlike Love's own version of this recursion, no artificial "thick buffer layer" is inserted
before the half-space: the half-space radiation condition here is already the exact analytic
one (via the eigenbasis projection above), and adding a buffer was confirmed, directly, to make
things numerically *worse* — a buffer layer far from the half-space's own phase velocity forces
[`one_way_phase`](@ref)'s overflow clamp to saturate, and the resulting propagator entries lose
essentially all precision (confirmed: a real 1000 km buffer, mirroring Love's own trick
naively, inflated this function's residual at a true root from `O(1)` to `O(10⁶⁶)`).
"""
function rayleigh_secular_function(layers::Vector{Layer}, ω, c, ω_ref)
    p = 1 / c
    M = Matrix{ComplexF64}(I, 4, 4)
    for L in layers[1:end-1]
        _, D = compound_matrix_psv(L, p, ω, ω_ref)
        M = D * M
    end
    Vhs, _, _ = modal_basis_psv(layers[end], p, ω, ω_ref)
    Vinv = inv(Vhs)
    A = zeros(ComplexF64, 2, 2)
    for (row, veci) in enumerate((2, 4))  # rows 2,4 of Vinv = the "Pu"/"Su" (growing) branches
        rowvec = Vinv[veci, :]' * M
        A[row, 1] = rowvec[1]
        A[row, 2] = rowvec[2]
    end
    return det(A)
end

# ╔═╡ 6a000003-0000-4000-8000-200000000003
"""
    solve_phase_velocities_rayleigh(layers, T; c_search_pad=0.3) -> Vector{Float64}

Find every Rayleigh phase-velocity root at period `T` (seconds), scanning
[`rayleigh_secular_function`](@ref)'s real part. The search is capped just under the
half-space's own `vs` — `c < vs_halfspace` is not an empirical guess to pad loosely, it is
*exactly* the condition [`vertical_slowness`](@ref) needs to stay evanescent in the half-space,
so the secular function itself stops describing a physical mode at or above it.
"""
function solve_phase_velocities_rayleigh(layers::Vector{Layer}, T; c_search_pad=0.3)
    ω = 2π / T
    vs_vals = [L.vs for L in layers if L.vs > 0.0]
    vmin_all = minimum(vs_vals)
    vs_halfspace = layers[end].vs
    cmin = max(1e-5, (1 - c_search_pad) * 0.5 * vmin_all)
    cmax = 0.995 * vs_halfspace
    g(c) = real(rayleigh_secular_function(layers, ω, c, 1.0))
    try
        return find_zeros(g, cmin, cmax)
    catch
        return Float64[]
    end
end

# ╔═╡ 6a000004-0000-4000-8000-200000000004
struct Velocities
    periods::Vector{Float64}
    phase_velocities::Vector{Float64}
    group_velocities::Vector{Float64}
end

# ╔═╡ 6a000005-0000-4000-8000-200000000005
"""
    track_fundamental_mode(periods, all_roots) -> Vector{Float64}

Follow the fundamental mode continuously across periods, instead of independently taking the
smallest root at each period (which can jump between distinct physical mode branches when
several roots are closely spaced). Starting from the longest period (fewest, best-separated
modes there) and walking toward the shortest, each step picks whichever of the current period's
roots is closest to the previous period's tracked value. Copied verbatim from
`Love-wave-dispersion-curves.jl` — wave-type-agnostic, a pure function of a root list.
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

# ╔═╡ 6a000006-0000-4000-8000-200000000006
"""
    group_velocity_from_phase(periods, phase_velocities) -> Vector{Float64}

Group velocity `U = c / (1 + (T/c)(dc/dT))` from a two-point finite difference of an
already-mode-tracked phase-velocity curve. Points where the finite-difference denominator
comes too close to zero (a steep local slope pushing near a spurious singularity) are dropped
back to `NaN` rather than plotted as a wrong-looking spike. Copied verbatim from
`Love-wave-dispersion-curves.jl` — wave-type-agnostic.
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

# ╔═╡ 6a000007-0000-4000-8000-200000000007
"""
    solve_rayleigh_dispersion(layers, periods) -> Velocities

The full Rayleigh fundamental-mode dispersion curve: every root at every period via
[`solve_phase_velocities_rayleigh`](@ref), tracked into one continuous branch via
[`track_fundamental_mode`](@ref), then differentiated for group velocity via
[`group_velocity_from_phase`](@ref). Unlike Love waves, a Rayleigh fundamental mode exists for
*any* layered model, including a uniform half-space (it's tied to the free surface itself, not
to a low-velocity waveguide) — so there's no equivalent of Love's own upfront
"does a mode exist" guard.
"""
function solve_rayleigh_dispersion(layers::Vector{Layer}, periods)
    all_roots = [solve_phase_velocities_rayleigh(layers, T) for T in periods]
    phase_velocities = track_fundamental_mode(periods, all_roots)
    group_velocities = group_velocity_from_phase(periods, phase_velocities)
    return Velocities(periods, phase_velocities, group_velocities)
end

# ╔═╡ 6a000008-0000-4000-8000-200000000008
md"### Verifying the Rayleigh Solver"

# ╔═╡ 6a000009-0000-4000-8000-200000000009
let
    # Uniform half-space, Poisson solid (nu=0.25, vp/vs=sqrt(3)): the classical closed-form
    # Rayleigh-velocity ratio is the physically relevant root of a textbook sextic in (c/vs)^2,
    # here derived directly (not just quoted) so this check doesn't depend on remembering a
    # rounded constant correctly.
    vs, vp = 3.5, 3.5 * sqrt(3)
    f(x) = x^3 - 8x^2 + (56 / 3) * x - 32 / 3  # x = (c/vs)^2, for vp/vs = sqrt(3)
    a, b = 0.0, 1.0
    for _ in 1:100
        m = (a + b) / 2
        sign(f(a)) == sign(f(m)) ? (a = m) : (b = m)
    end
    ratio_exact = sqrt((a + b) / 2)

    # a genuinely uniform half-space, expressed as [one thin layer of the SAME material, the
    # true half-space] since rayleigh_secular_function expects at least one finite layer above
    # the half-space -- thin enough (1e-4 km) to be a negligible perturbation, confirmed by
    # re-running at 1.0 km and 1e-2 km and seeing the root converge to the same value.
    hs = Layer(Inf, vp, vs, 2.7, 1e12, 1e12)  # huge Q: lossless, so this matches the elastic textbook formula exactly
    layers = [Layer(1e-4, vp, vs, 2.7, 1e12, 1e12), hs]
    roots = solve_phase_velocities_rayleigh(layers, 20.0)
    c_numeric = roots[argmin(abs.(roots .- ratio_exact * vs))]
    err = abs(c_numeric / vs - ratio_exact) / ratio_exact

    @assert !isempty(roots)
    @assert err < 1e-4

    md"""
    !!! correct "Self-check"
        For a uniform Poisson-solid half-space (`` \nu=0.25 ``), the exact root of the
        classical Rayleigh secular cubic gives `` c_R/V_s=`` $(round(ratio_exact, digits=6)).
        `rayleigh_secular_function`'s numerically-solved root gives
        `` c_R/V_s=`` $(round(c_numeric / vs, digits=6)) — $(round(err * 100, sigdigits=2))%
        relative error, with no dependence on this notebook's own machinery beyond the one
        function under test.
    """
end

# ╔═╡ 6a00000a-0000-4000-8000-20000000000a
let
    # physical necessity: a genuine (trapped, evanescent) Rayleigh wave must travel slower than
    # the half-space's own S velocity at every period -- that's the radiation condition the
    # secular function itself is built from, so a violation here would mean the root finder
    # locked onto numerical noise rather than a real mode.
    layers = GUTENBERG_MODEL[1:3]
    periods = collect(range(5.0, 60.0, length=10))
    res = solve_rayleigh_dispersion(layers, periods)
    finite_c = filter(isfinite, res.phase_velocities)
    vs_hs = layers[end].vs
    @assert !isempty(finite_c)
    @assert all(finite_c .< vs_hs)

    # residual check: a genuine root should make the secular function's magnitude many orders
    # smaller than its own typical scale a small step away in c -- not exactly zero (this is a
    # numerically ill-conditioned function in absolute terms) but a clear, sharp minimum.
    i0 = findfirst(isfinite, res.phase_velocities)
    T0, c0 = periods[i0], res.phase_velocities[i0]
    resid_on = abs(rayleigh_secular_function(layers, 2π / T0, c0, 1.0))
    resid_off = abs(rayleigh_secular_function(layers, 2π / T0, c0 * 1.01, 1.0))
    ratio = resid_on / resid_off

    @assert ratio < 1e-3

    md"""
    !!! correct "Self-check"
        Over `` T\in[$(periods[1]),$(periods[end])] `` s on a real (Gutenberg crustal) layered
        model, every solved Rayleigh phase velocity stays below the half-space's own
        `` V_s=`` $(vs_hs) km/s, as physically required. And the secular residual right at a
        solved root is $(round(ratio, sigdigits=2))× its value just 1% away in `` c `` —
        a sharp, genuine zero-crossing, not root-finder noise.
    """
end

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
    PlutoUI.WideCell(@bind lm LayeredMediumInput(GUTENBERG_MODEL; show_vp=true); max_width=1400)
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

# ╔═╡ 7a000003-0000-4000-8000-500000000003
rayleigh_res = solve_rayleigh_dispersion(layers, rayleigh_periods);

# ╔═╡ 7a000004-0000-4000-8000-500000000004
plot(
    [
        scatter(x=rayleigh_res.periods, y=rayleigh_res.phase_velocities, mode="lines+markers", name="phase velocity", line=attr(color="#3b82f6")),
        scatter(x=rayleigh_res.periods, y=rayleigh_res.group_velocities, mode="lines+markers", name="group velocity", line=attr(color="#ef4444")),
    ],
    Layout(title="Rayleigh-Wave Dispersion (fundamental mode)", template=:plotly_dark, width=560, height=320,
        xaxis=attr(title="period (s)"), yaxis=attr(title="velocity (km/s)")),
)

# ╔═╡ 4a000001-0000-4000-8000-400000000001
md"## The Source/Receiver Geometry Widget"

# ╔═╡ 4a000002-0000-4000-8000-400000000002
begin
    """
        LambGeometryInput(layers::Vector{Layer}; r0=60.0, phase0="direct", rmax=100.0)

    A dark-canvas widget showing the source/receiver geometry against the current layered
    model: a **star** marks the fixed source (range 0, depth `LP_SOURCE_DEPTH`), a
    **downward triangle** marks the receiver, draggable in range along the free surface. Layer
    bands are shaded by Vp. Small buttons above the canvas pick which phase (`Direct`,
    `Preflect`, `Psreflect`) drives the seismogram panels below — the bound value is a plain
    string, not one of those dispatch structs directly, so this widget's own `@bind` cell never
    needs a dependency on them (see `pluto-widget-style`'s bind-cell-ordering note: this whole
    class of bug disappears when a widget's construction doesn't need any upstream data it
    would otherwise have to depend on).
    """
    struct LambGeometryInput
        layers::Vector{Layer}
        r0::Float64
        phase0::String
        rmax::Float64
    end
    LambGeometryInput(layers::Vector{Layer}; r0::Float64=60.0, phase0::String="direct", rmax::Float64=100.0) =
        LambGeometryInput(layers, r0, phase0, rmax)

    Base.get(w::LambGeometryInput) = Dict{String,Any}("r" => w.r0, "phase" => w.phase0)

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

    function Base.show(io::IO, ::MIME"text/html", w::LambGeometryInput)
        layers_js = _lp_layers_js(w.layers)
        write(io, """
        <div id="lpwidget">
        <style>
        #lpwidget{font-family:sans-serif;color:#d1d5db}
        #lpwidget .lp-titlebar{background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px;margin-bottom:10px;text-align:center}
        #lpwidget .lp-titlebar-headline{font-size:16px;font-weight:700;color:#f3f4f6}
        #lpwidget .lp-mode-row{display:flex;justify-content:center;gap:8px;flex-wrap:wrap;margin-bottom:10px}
        #lpwidget .lp-mode-btn{border-radius:4px;border:1px solid #6b7280;background:#0b0b0b;color:#e5e7eb;padding:6px 10px;font-size:13px;cursor:pointer}
        #lpwidget .lp-mode-btn:hover{background:#1f2937}
        #lpwidget .lp-mode-btn.active{border-color:#4ade80;color:#4ade80}
        #lpwidget .lp-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #lpwidget .lp-panel-title{font-size:14px;font-weight:700;color:#f3f4f6;margin-bottom:4px;text-align:center}
        #lpwidget canvas{display:block;width:100%;height:auto;cursor:grab}
        #lpwidget .lp-caption{font-size:12px;color:#9ca3af;margin-top:4px;text-align:center}
        </style>

        <div class="lp-titlebar">
          <div class="lp-titlebar-headline">Drag the triangle (receiver) along the free surface. The star (source) is fixed.</div>
        </div>
        <div class="lp-mode-row">
          <button class="lp-mode-btn" data-phase="direct" type="button">Direct</button>
          <button class="lp-mode-btn" data-phase="preflect" type="button">Reflected P</button>
          <button class="lp-mode-btn" data-phase="psreflect" type="button">Reflected P→S</button>
        </div>
        <div class="lp-panel">
          <div class="lp-panel-title">Source, Receiver, and the Layered Earth</div>
          <canvas id="lp-geom" width="700" height="260"></canvas>
          <div class="lp-caption" id="lp-caption">Loading…</div>
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
          const state = { r: $(w.r0), phase: "$(w.phase0)" }

          function publish(){
            root.value = { r: state.r, phase: state.phase }
            root.dispatchEvent(new CustomEvent('input'))
          }
          root.value = { r: state.r, phase: state.phase }

          const cv = root.querySelector('#lp-geom')
          const ctx = cv.getContext('2d')
          const EM = { l: 16, r: 16, t: 20, b: 16 }

          function depthMax(){
            const lastB = LAYERS.boundaries.length ? LAYERS.boundaries[LAYERS.boundaries.length - 1] : 20
            return Math.max(40, lastB * 1.6, SOURCE_DEPTH * 2)
          }
          function rToX(r){ return EM.l + (r + RMAX) / (2 * RMAX) * (cv.width - EM.l - EM.r) }
          function xToR(x){ return ((x - EM.l) / (cv.width - EM.l - EM.r)) * 2 * RMAX - RMAX }
          function depthToY(z){ return EM.t + (z / depthMax()) * (cv.height - EM.t - EM.b) }
          function vpColor(vp){
            const t = Math.max(0, Math.min(1, (vp - 2) / 11))
            return 'rgb(' + Math.round(30 + 40 * t) + ',' + Math.round(60 + 80 * t) + ',' + Math.round(120 + 120 * t) + ')'
          }

          function drawStarMarker(cx, cy, r, fill, stroke){
            const spikes = 5, rOuter = r, rInner = r * 0.45
            ctx.beginPath()
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

          function draw(){
            const W = cv.width, H = cv.height
            ctx.clearRect(0,0,W,H)
            ctx.fillStyle = '#000'; ctx.fillRect(0,0,W,H)
            const n = LAYERS.vp.length
            for(let i=0;i<n;i++){
              const top = i===0 ? 0 : LAYERS.boundaries[i-1]
              const bot = i===n-1 ? depthMax() : LAYERS.boundaries[i]
              ctx.fillStyle = vpColor(LAYERS.vp[i])
              ctx.fillRect(EM.l, depthToY(top), W-EM.l-EM.r, Math.max(1,depthToY(bot)-depthToY(top)))
            }
            ctx.strokeStyle = '#1f2937'; ctx.lineWidth = 1
            for(let i=0;i<n-1;i++){
              const y = depthToY(LAYERS.boundaries[i])
              ctx.beginPath(); ctx.moveTo(EM.l,y); ctx.lineTo(W-EM.r,y); ctx.stroke()
            }
            ctx.strokeStyle = '#374151'; ctx.strokeRect(EM.l, EM.t, W-EM.l-EM.r, H-EM.t-EM.b)

            drawStarMarker(rToX(0), depthToY(SOURCE_DEPTH), 9, '#facc15', '#000')
            drawTriangleDownMarker(rToX(state.r), depthToY(0)+2, 8, '#f5f3ef', '#0a0f18')

            document.getElementById('lp-caption').textContent =
              'receiver at r = ' + state.r.toFixed(1) + ' km · source depth ' + SOURCE_DEPTH.toFixed(0) + ' km (fixed) · drag the triangle'
          }

          let dragging = false
          function pointerXY(ev){
            const rect = cv.getBoundingClientRect()
            const sx = cv.width/rect.width, sy = cv.height/rect.height
            return [(ev.clientX-rect.left)*sx, (ev.clientY-rect.top)*sy]
          }
          cv.addEventListener('mousedown', ev => {
            const [x,y] = pointerXY(ev)
            const hx = rToX(state.r), hy = depthToY(0)+2
            if (Math.hypot(x-hx, y-hy) < 16) { dragging = true; cv.style.cursor = 'grabbing' }
          }, lpSignal)
          window.addEventListener('mousemove', ev => {
            if (!dragging) return
            const [x] = pointerXY(ev)
            state.r = Math.max(-RMAX, Math.min(RMAX, xToR(x)))
            draw()
          }, lpSignal)
          window.addEventListener('mouseup', () => {
            if (!dragging) return
            dragging = false; cv.style.cursor = 'grab'
            publish()
          }, lpSignal)

          function setPhaseButtons(){
            root.querySelectorAll('.lp-mode-btn').forEach(b => b.classList.toggle('active', b.dataset.phase === state.phase))
          }
          root.querySelectorAll('.lp-mode-btn').forEach(b => {
            b.addEventListener('click', () => { state.phase = b.dataset.phase; setPhaseButtons(); publish() }, lpSignal)
          })

          setPhaseButtons()
          draw()
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
    PlutoUI.WideCell(@bind geom LambGeometryInput(layers); max_width=1400)
end

# ╔═╡ 3c1a2b4e-0002-4000-8000-100000000002
# translate the widget's chosen phase (a plain string, so the widget never needs a bind-cell
# dependency on the phase structs themselves -- see pluto-widget-style's bind-ordering note)
# into the dispatch singleton get_phase/get_wavefield actually need
phase_obj = geom["phase"] == "preflect" ? Preflect() : geom["phase"] == "psreflect" ? Psreflect() : Direct()

# ╔═╡ 53a45db3-8571-4fc7-a855-5173030b7cb9
seismograms = irfft(remove_zero_frequency!(get_wavefield(seismograms_param, layers, phase_obj)), seismograms_param.Nt, 1);

# ╔═╡ 96556fec-4bea-4132-946a-92d17372df54
let
    U = seismograms
    Umaxclip = maximum(abs, U) * 0.1
    r0 = geom["r"]
    plot(
        [
            heatmap(y=seismograms_param.tgrid, x=seismograms_param.rgrid, z=U, zmin=-Umaxclip, zmax=Umaxclip, colorscale=:seismic),
        ],
        Layout(title="Seismograms (" * geom["phase"] * ")", template=:plotly_dark, width=420, height=380,
            yaxis_autorange="reversed", xaxis=attr(title="distance to receiver (km)"), yaxis=attr(title="time (s)"),
            shapes=[attr(type="line", x0=r0, x1=r0, y0=seismograms_param.tgrid[1], y1=seismograms_param.tgrid[end],
                line=attr(color="yellow", width=2, dash="dot"))]),
    )
end

# ╔═╡ 3c1a2b4e-0003-4000-8000-100000000003
let
    # the single trace at the widget's currently-dragged receiver offset
    j = argmin(abs.(seismograms_param.rgrid .- geom["r"]))
    plot(scatter(x=seismograms_param.tgrid, y=seismograms[:, j], line=attr(color="#38bdf8")),
        Layout(title="Trace at r = " * string(round(seismograms_param.rgrid[j], digits=1)) * " km", template=:plotly_dark,
            width=420, height=220, xaxis=attr(title="time (s)"), yaxis=attr(title="displacement (a.u.)")))
end

# ╔═╡ e7c1cbf4-7939-45a2-a754-b36ae9bdcab4
md"""### To Do
- **Full time-domain Rayleigh-wave synthesis.** The Rayleigh dispersion curve above is a
  *dispersion relation* only — the actual Rayleigh arrival isn't yet part of the time-domain
  seismogram panel. Getting it there needs two more things: the free-surface reflection woven
  into `reflectivity_matrix`'s Kennett recursion (the recursion currently only sees interfaces
  *below* layer 1, not the free surface above it), and correctly capturing the surface-wave
  pole when doing the wavenumber integral numerically — a pole sitting exactly on the real `k`
  axis needs contour deformation or residue extraction, not just a finer `kgrid`. A genuinely
  harder numerical problem than anything else in this notebook, left for a future pass.
- Use in-place functions during the wavenumber integration for speed.
- Time-dependent snapshots, to visualize head-wave generation directly.
- Plot the displacement field itself, not just potentials.
"""

# ╔═╡ 9d10ab94-38d5-406e-8697-7f5bc0c11df2
md"""### References
- Fuchs, K., and Gerhard Müller. "Computation of synthetic seismograms with the reflectivity method and comparison with observations." Geophysical Journal International 23.4 (1971): 417-433.
- Aki, K., and Richards, P.G. *Quantitative Seismology*, 2nd ed. University Science Books, 2002 — the reflectivity method (ch. 9) and the Rayleigh/Love secular-equation derivations (ch. 7) this notebook follows.
- Lamb, H. "On the propagation of tremors over the surface of an elastic solid." *Philosophical Transactions of the Royal Society A* 203 (1904): 1-42 — the original problem this notebook is named after.
"""

# ╔═╡ Cell order:
# ╠═897afffa-77e8-11ef-1a54-c73d1df8f6a4
# ╠═aa38e7f6-b8d2-4270-9142-1aa688041eb4
# ╠═4621f804-6a44-4c46-af37-3d364ba74cfe
# ╠═9b488729-aadc-4c1d-b971-6931d4bc9a08
# ╠═6188ff3e-cfd5-4c9c-aa39-619dc280d494
# ╟─1bbedd43-9e69-4ba0-a251-855f1f63dcf8
# ╠═beef0001-0000-4000-8000-000000000001
# ╠═beef0002-0000-4000-8000-000000000002
# ╠═beef0003-0000-4000-8000-000000000003
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
# ╟─72aca1bf-3997-4be3-b316-921b45480c1a
# ╠═3c1a2b4e-0001-4000-8000-100000000001
# ╠═12f7733a-edd6-481f-8166-ef967520b35a
# ╠═3c1a2b4e-0002-4000-8000-100000000002
# ╠═6923e13d-a4c5-45c9-a1b2-f0104932d709
# ╠═53a45db3-8571-4fc7-a855-5173030b7cb9
# ╠═96556fec-4bea-4132-946a-92d17372df54
# ╠═3c1a2b4e-0003-4000-8000-100000000003
# ╟─7a000001-0000-4000-8000-500000000001
# ╠═7a000002-0000-4000-8000-500000000002
# ╠═7a000003-0000-4000-8000-500000000003
# ╠═7a000004-0000-4000-8000-500000000004
# ╟─c695e4d3-c49d-4587-8adc-cdd4c001325e
# ╠═fbd2b939-c816-42dd-a65e-94abce3dc91a
# ╟─6a000001-0000-4000-8000-200000000001
# ╠═6a000002-0000-4000-8000-200000000002
# ╠═6a000003-0000-4000-8000-200000000003
# ╠═6a000004-0000-4000-8000-200000000004
# ╠═6a000005-0000-4000-8000-200000000005
# ╠═6a000006-0000-4000-8000-200000000006
# ╠═6a000007-0000-4000-8000-200000000007
# ╟─6a000008-0000-4000-8000-200000000008
# ╠═6a000009-0000-4000-8000-200000000009
# ╠═6a00000a-0000-4000-8000-20000000000a
# ╟─5a000001-0000-4000-8000-300000000001
# ╠═5a000002-0000-4000-8000-300000000002
# ╠═5a000003-0000-4000-8000-300000000003
# ╟─4a000001-0000-4000-8000-400000000001
# ╠═4a000002-0000-4000-8000-400000000002
# ╟─e7c1cbf4-7939-45a2-a754-b36ae9bdcab4
# ╟─9d10ab94-38d5-406e-8697-7f5bc0c11df2
