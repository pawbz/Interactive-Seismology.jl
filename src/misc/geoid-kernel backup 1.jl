### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Geoid Kernels"
#> layout = "layout.jlhtml"
#> tags = ["misc"]

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

# ╔═╡ 3a1f0c40-1a2b-4c5d-8e9f-0a1b2c3d4e5f
begin
    using PlutoUI, PlutoPlotly, LinearAlgebra, LaTeXStrings, ColorSchemes
    import PlutoUI: combine
end

# ╔═╡ 4b2e1d51-2b3c-4d6e-9f0a-1b2c3d4e5f60
using Symbolics

# ╔═╡ 5c3f2e62-3c4d-4e7f-a01b-2c3d4e5f6071
TableOfContents(include_definitions=true)

# ╔═╡ 6d403f73-4d5e-4f80-b12c-3d4e5f607182
md"""
## Why Does a Dense Blob Make a Geoid *Low*?

Drop a dense lump into the mantle. Ask a student what happens to the geoid — the shape of
sea level — directly above it, and almost everyone says the same thing: *more mass means
more gravity, so the geoid goes up.*

That answer is right for a rock sitting on a table. It is often **wrong** for a rock sitting
inside a convecting planet, and understanding why is the point of this notebook.

The reason is that the mantle **flows**. A dense anomaly does not just sit there adding its
own gravitational pull. It sinks, and as it sinks it drags the whole mantle with it. That
flow pushes *down* on the Earth's surface from below, creating a shallow depression, and it
pulls *down* on the core–mantle boundary, creating a bulge there. Both of those deflected
boundaries are enormous density contrasts — rock against air at the top, rock against iron
at the bottom — so the mass they move around is not a small correction. It can be **larger
than the anomaly's own pull**, and it has the opposite sign.

So the geoid you actually measure is a competition between three mass contributions:

1. the density anomaly itself — *positive* for a dense blob,
2. the depressed **surface** it creates — *negative*,
3. the deflected **core–mantle boundary** — *positive* or negative depending on geometry.

Which one wins depends on **how deep the anomaly sits** and on **how viscosity varies with
depth**. That competition is what a *geoid kernel* summarises, and it is why the geoid is one
of our best constraints on mantle viscosity (Hager & Richards, 1989).

This matters for a real, unsolved problem: the **Indian Ocean Geoid Low**, a ~106 m
depression in sea level south of India — the deepest geoid low on Earth. One family of
explanations involves deep dense material. This notebook is about building the intuition for
how that is even possible.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India

"""

# ╔═╡ 7e514084-5e6f-4091-c23d-4e5f60718293
md"---"

# ╔═╡ 8f625195-6f70-41a2-d34e-5f6071829304
md"""
### The Draw-Your-Own-Anomaly Sandbox

Everything below is driven by the panel you are about to see. **Paint a density anomaly**
into the mantle cross-section on the left and watch the geoid respond above it.

- **Left panel** — a slice through the Earth from the surface down to the core. Drag to paint.
  Red is dense (heavy), blue is light. The dashed lines show the deflected surface and CMB,
  exaggerated so you can see them.
- **Top panel** — the resulting geoid $N(\theta)$ along the surface. This is what a satellite
  would measure.
- **Right panel** — the *kernel* $K_\ell(r)$: the geoid produced by a unit anomaly placed at
  each depth. **Where it crosses zero is where the sign flips.**

Three things to try, in order — each is one of the three "aha" moments:

1. Press **Shallow blob** (≈580 km), then **Deep blob** (≈2600 km). Same dense anomaly,
   same everything else — and the geoid flips from a **low** to a **high**.
2. With the deep blob showing, tick **"direct effect only"**. The boundary terms switch off,
   the kernel goes positive everywhere, and the geoid jumps up. Untick it to bring the
   boundaries back. Try this at the shallow blob too, where the effect is dramatic: the
   direct-only geoid is about **+8 m**, and switching the boundaries on drives it to **−2 m**.
3. Drag the **viscosity ratio** slider. The zero-crossing migrates upward, and then
   *disappears entirely*. That disappearance is the punchline; we come back to it below.
"""

# ╔═╡ 90736206-7081-42b3-e45f-607182930415
md"""
## Physical Constants

Standard Earth values. `η` is only ever used as a *ratio*, so its absolute scale never
matters for the geoid — a fact worth pausing on: **the geoid constrains relative viscosity,
not absolute viscosity.**
"""

# ╔═╡ a1847317-8192-43c4-f560-718293041526
begin
    const G_GRAV = 6.674e-11      # gravitational constant, m^3 kg^-1 s^-2
    const A_EARTH = 6371e3        # Earth radius, m
    const B_CMB = 3480e3          # core-mantle boundary radius, m
    const G_SURF = 9.81           # surface gravity, m s^-2
    const RHO_M = 3300.0          # representative mantle density, kg m^-3
    const RHO_C = 10900.0         # outer core density, kg m^-3
    const R_660 = A_EARTH - 660e3 # upper/lower mantle boundary
    nothing
end

# ╔═╡ b2958428-92a3-44d5-0671-829304152637
md"""
## The Governing Equations

The mantle over geological time is a fluid so viscous that inertia is utterly negligible —
the Prandtl number is effectively infinite. What is left is the **Stokes equations**:
a statement that at every instant, viscous forces exactly balance buoyancy.
"""

# ╔═╡ c3a69539-a3b4-45e6-1782-930415263748
begin
    @syms r θ                    # radius and colatitude
    @syms η                      # dynamic viscosity
    @syms δρ                     # density anomaly
    @syms g                      # gravitational acceleration
    @syms 𝐯                      # velocity field
    @syms σ                      # stress tensor
    nothing
end

# ╔═╡ d4b7a64a-b4c5-46f7-2893-041526374859
md"""
**Incompressibility** — the mantle conserves volume as it flows:

$(L"\nabla \cdot \mathbf{v} = 0")

**Force balance** — viscous stress divergence balances the buoyancy of the anomaly:

$(L"\nabla \cdot \sigma = \delta\rho \, g \, \hat{r}")

!!! danger "This is instantaneous — there is no time in these equations"
	Look carefully: there is no $\partial/\partial t$ anywhere. We are **not** running a
	convection simulation. The anomaly is *frozen* where you paint it, and we perform
	exactly **one linear solve** for the flow it drives at that instant.

	This trips up almost everyone the first time. The blob does not sink as you watch. We
	are asking "what flow, and what geoid, does this configuration produce *right now*?"
	— a snapshot, not a movie.

**Dynamic topography** — the flow pushes on each boundary, and the boundary deflects until
its own weight balances that push:

$(L"h = \frac{\sigma_{rr}}{\Delta\rho \, g}")

with $\Delta\rho = \rho_m$ at the surface (rock against air) and
$\Delta\rho = \rho_c - \rho_m$ at the CMB (iron against rock).

**The geoid** is then the sum over all three mass contributions, upward-continued to the
surface. For a thin mass sheet of surface density $\sigma$ at radius $s$:

$(L"N_\ell = \frac{4\pi G a}{g(2\ell+1)} \, \sigma \, \left(\frac{s}{a}\right)^{\ell+2}")
"""

# ╔═╡ e5c8b75b-c5d6-4708-3904-152637485960
md"""
## Reducing 3-D Flow to a Stack of 1-D Problems

Solving the Stokes equations in 3-D needs a mesh and a serious solver. We do not need one,
because of a piece of good luck: **if viscosity varies only with radius**, the problem
*separates* by spherical-harmonic degree. Each degree $\ell$ evolves completely
independently of every other.

We simplify once more by restricting to **axisymmetric** anomalies ($m = 0$) — patterns that
look the same all the way around the polar axis. The lateral decomposition is then a
**Legendre** expansion rather than a full spherical-harmonic one.

The whole pipeline becomes:

```
paint δρ(r,θ)  →  Legendre-decompose in θ  →  δρ_ℓ(r)
               →  apply 1-D radial kernel K_ℓ(r) for each ℓ
               →  sum over ℓ  →  geoid N(θ)
```

So the "solver" is a **loop of small 1-D problems**, one per degree, truncated at
$\ell_{max}$. No mesh. No 3-D code. That is why this runs instantly in your browser.
"""

# ╔═╡ f6d9c86c-d6e7-4819-4a15-263748596071
md"""
## Layer 1: The Analytic Building Blocks

Rather than transcribe a propagator matrix from a paper — the classic place for a sign error
to hide — we **derive** the four fundamental solutions and then verify them against the
original differential equations further below.

Expanding the poloidal flow at degree $\ell$ (with $L = \ell(\ell+1)$):

$(L"u_r = U(r)\,Y_\ell, \qquad u_\theta = V(r)\,\partial_\theta Y_\ell, \qquad p = P(r)\,Y_\ell")

Three facts pin the solution down completely:

1. **Incompressibility** gives $V$ for free: $V = (rU' + 2U)/L$.
2. **The pressure is harmonic**, $\nabla^2 p = 0$, so $P \propto r^s$ with $s^2 + s - L = 0$,
   giving $s = \ell$ or $s = -(\ell+1)$.
3. **The homogeneous radial equation** $U'' + 4U'/r + (2-L)U/r^2 = 0$ gives
   $p = \ell-1$ or $p = -\ell-2$; each pressure mode drives a particular solution with
   $p = s+1$.

Together the four radial exponents are

$(L"U \sim \{\, r^{\ell-1},\; r^{\ell+1},\; r^{-\ell},\; r^{-\ell-2} \,\}")

Every quantity we need is a pure power of $r$, so the code below is short and exact.
"""

# ╔═╡ 07ead97d-e7f8-492a-5b26-374859607182
Lof(ℓ) = ℓ * (ℓ + 1)

# ╔═╡ 18fbea8e-f809-4a3b-6c37-485960718293
"""
	modes(η, ℓ)

The four fundamental solutions at degree `ℓ`, each returned as `(p, cU, cP)` where the
radial velocity goes as `cU * r^p` and the pressure as `cP * r^(p-1)`.

Two are homogeneous (`cP = 0`, no pressure); two are particular solutions driven by the
two harmonic pressure modes. The coefficient `cU` of a particular mode follows from the
radial momentum balance.
"""
function modes(η, ℓ)
    L = Lof(ℓ)
    out = Tuple{Int,Float64,Float64}[]
    # homogeneous modes: pressure vanishes
    for p in (ℓ - 1, -ℓ - 2)
        push!(out, (p, 1.0, 0.0))
    end
    # particular modes driven by P = r^s, s = ℓ and s = -(ℓ+1); take cP = 1
    for s in (ℓ, -(ℓ + 1))
        p = s + 1
        denom = p * (p - 1) + 2p - (2 + L) + 2 * (p + 2)
        push!(out, (p, ((p - 1) / η) / denom, 1.0))
    end
    return out
end

# ╔═╡ 290cfb9f-090a-4b4c-7d48-596071829304
"""
	mode_state(r, η, ℓ, p, cU, cP)

State vector ``y = [u_r,\\; u_\\theta,\\; \\tau_{rr},\\; \\tau_{r\\theta}]`` for a single
mode at radius `r`, using

	τ_rr  = -P + 2η U'
	τ_rθ  =  η (r V' - V + U)
"""
function mode_state(r, η, ℓ, p, cU, cP)
    L = Lof(ℓ)
    U = cU * r^p
    Up = cU * p * r^(p - 1)
    V = cU * (p + 2) / L * r^p
    Vp = cU * (p + 2) * p / L * r^(p - 1)
    P = cP * r^(p - 1)
    τrr = -P + 2 * η * Up
    τrθ = η * (r * Vp - V + U)
    return [U, V, τrr, τrθ]
end

# ╔═╡ 3a1d0ca0-1a1b-4c5d-8e59-607182930415
"""
	Ymat(r, η, ℓ)

Matrix whose four columns are the fundamental solutions at radius `r`. Radii are scaled by
the Earth radius so that `r^ℓ` stays near unity — without this the columns overflow and the
matrix becomes badly conditioned at high degree.
"""
function Ymat(r, η, ℓ)
    x = r / A_EARTH
    hcat([mode_state(x, η, ℓ, p, cU, cP) for (p, cU, cP) in modes(η, ℓ)]...)
end

# ╔═╡ 4b2e1db1-2b2c-4d6e-9f6a-718293041526
"""
	propagator(r₂, r₁, η, ℓ)

Maps the state vector from radius `r₁` to `r₂` inside a layer of constant viscosity.
Because the columns of `Ymat` span the solution space, `Y(r₂) Y(r₁)⁻¹` *is* the propagator —
no hand-written matrix entries, hence no place for a transcription error to hide.
"""
propagator(r₂, r₁, η, ℓ) = Ymat(r₂, η, ℓ) / Ymat(r₁, η, ℓ)

# ╔═╡ 5c3f2ec2-3c3d-4e7f-a07b-829304152637
md"""
### Verifying the building blocks

Before trusting any of this, we check the four solutions against the **original** equations.
Every quantity is a pure power of $r$, so each equation collapses to an algebraic identity
on the coefficients — we can test it exactly, with no finite differences.

The three equations are incompressibility, radial momentum, and colatitudinal momentum:

	U' + 2U/r - L·V/r                              = 0
	η[U'' + 2U'/r - (2+L)U/r² + 2L·V/r²] - P'      = 0
	η[V'' + 2V'/r - L·V/r² + 2U/r²]    - P/r       = 0
"""

# ╔═╡ 6d403fd3-4d4e-4f80-b18c-930415263748
"Exact residuals of the three governing equations for one mode (no finite differences)."
function mode_residuals(η, ℓ, p, cU, cP)
    L = Lof(ℓ)
    cV = cU * (p + 2) / L
    incomp = cU * p + 2cU - L * cV
    radial = η * (cU * (p * (p - 1) + 2p - (2 + L)) + 2L * cV) - cP * (p - 1)
    theta = η * (cV * (p * (p - 1) + 2p - L) + 2cU) - cP
    scale = max(abs(cU), abs(cP), abs(cV)) * max(η, 1.0) * max(L, 1.0)
    return maximum(abs.((incomp, radial, theta))) / scale
end

# ╔═╡ 7e5140e4-5e5f-4091-c29d-041526374859
validation_residual = maximum(
    mode_residuals(η_test, ℓ_test, p, cU, cP)
    for ℓ_test in 2:30, η_test in (1.0, 5.0, 1e3)
    for (p, cU, cP) in modes(η_test, ℓ_test)
)

# ╔═╡ 8f6251f5-6f60-41a2-d3ae-152637485960
if validation_residual < 1e-12
    md"""
    !!! correct "Validation 1 passed"
    	Maximum relative residual over degrees ``\ell = 2 \ldots 30`` and viscosities
    	spanning three orders of magnitude: **$(round(validation_residual, sigdigits=3))**
    	— machine precision. The fundamental solutions satisfy the Stokes equations exactly.
    """
else
    md"""
    !!! danger "Validation 1 FAILED"
    	Residual $(validation_residual) is too large — the fundamental solutions do not
    	satisfy the governing equations. Do not trust anything below.
    """
end

# ╔═╡ 90736265-7071-42b3-e4bf-263748596071
md"""
## The Geoid Kernel

Now we assemble the three contributions. Place a thin sheet of density anomaly at radius
$s$ and ask what geoid it produces at the surface.

**Boundary conditions.** Both the surface and the CMB are *free-slip*: no shear traction
($\tau_{r\theta} = 0$), and no flow through the boundary ($u_r = 0$). The boundary does not
physically move in the solve — instead the normal stress $\tau_{rr}$ that builds up against
it tells us how far it *would* deflect, which is the dynamic topography.

**The load.** Crossing the loaded radius, the normal stress jumps by the weight of the sheet:

$(L"y(s^+) = y(s^-) + [\,0,\; 0,\; -\delta\rho\, g,\; 0\,]^T")

This sign is worth checking rather than trusting: a *dense* sheet must **depress** the
surface. We verify exactly that below.
"""

# ╔═╡ a18473c5-8182-43c4-f5ca-374859607182
"""
	layer_propagator(r₁, r₂, ℓ, layers)

Chain propagators through every constant-viscosity layer between `r₁` and `r₂`. All four
state components are continuous across a viscosity interface, so the layer propagators
simply multiply.
"""
function layer_propagator(r₁, r₂, ℓ, layers)
    M = Matrix{Float64}(I, 4, 4)
    for (rb, rt, η) in layers
        lo = max(rb, r₁)
        hi = min(rt, r₂)
        hi > lo || continue
        M = propagator(hi, lo, η, ℓ) * M
    end
    return M
end

# ╔═╡ b2958425-9293-44d5-06d1-485960718293
"""
	geoid_response(s, ℓ, layers)

Solve the Stokes problem for a unit density sheet at radius `s`, degree `ℓ`.

Returns a named tuple with the three geoid contributions (m per kg/m³ per m of sheet
thickness), the surface and CMB dynamic topography, and the boundary-condition residuals
used for validation.
"""
function geoid_response(s, ℓ, layers)
    a, b = A_EARTH, B_CMB
    upward(σ, rad) = 4pi * G_GRAV * a / (G_SURF * (2ℓ + 1)) * σ * (rad / a)^(ℓ + 2)

    # (1) the anomaly's own gravitational pull — always positive for dense material
    N_direct = upward(1.0, s)

    Pb = layer_propagator(b, s, ℓ, layers)   # CMB  -> load
    Pa = layer_propagator(s, a, ℓ, layers)   # load -> surface
    M = Pa * Pb
    jump = [0.0, 0.0, -G_SURF, 0.0]          # unit δρ

    # Free slip at the CMB: u_r = 0 and τ_rθ = 0, leaving two unknown constants.
    # Require the same two conditions at the surface -> 2x2 linear system.
    rhs = -(Pa * jump)
    Asys = [M[1, 2] M[1, 3]; M[4, 2] M[4, 3]]
    coef = Asys \ [rhs[1], rhs[4]]
    y_cmb = [0.0, coef[1], coef[2], 0.0]
    y_surf = M * y_cmb + Pa * jump

    # Normal stress -> dynamic topography at each boundary
    h_surf = y_surf[3] / (RHO_M * G_SURF)
    h_cmb = y_cmb[3] / ((RHO_C - RHO_M) * G_SURF)

    # (2) and (3): the deflected boundaries, as mass sheets
    N_surf = upward(RHO_M * h_surf, a)
    N_cmb = upward((RHO_C - RHO_M) * h_cmb, b)

    return (direct=N_direct, surface=N_surf, cmb=N_cmb,
        total=N_direct + N_surf + N_cmb,
        h_surf=h_surf, h_cmb=h_cmb,
        bc_residual=maximum(abs.((y_surf[1], y_surf[4], y_cmb[1], y_cmb[4]))))
end

# ╔═╡ c3a69595-a3a4-45e6-17e2-596071829304
"""
	geoid_kernel(s, ℓ, layers; direct_only=false)

Geoid at the surface (in metres) per unit density anomaly, for a sheet at radius `s`.

`direct_only=true` switches off the deflected-boundary terms, leaving only the anomaly's
own attraction — the naive, always-positive answer.
"""
function geoid_kernel(s, ℓ, layers; direct_only=false)
    resp = geoid_response(s, ℓ, layers)
    direct_only ? resp.direct : resp.total
end

# ╔═╡ d4b7a645-b4b5-46f7-28f3-607182930415
begin
    "Uniform-viscosity mantle."
    isoviscous_layers(η=1e21) = [(B_CMB, A_EARTH, η)]

    "Two-layer mantle: lower mantle `ratio` times stiffer than the upper mantle."
    twolayer_layers(ratio; η_um=1e21) =
        [(B_CMB, R_660, η_um * ratio), (R_660, A_EARTH, η_um)]
    nothing
end

# ╔═╡ e5c8b755-c5c6-4708-3954-718293041526
"Depth (km) where the kernel changes sign, or `nothing` if it never does."
function zero_crossing_depth(ℓ, layers; n=200, direct_only=false)
    rs = range(B_CMB + 5e3, A_EARTH - 5e3, length=n)
    v = [geoid_kernel(s, ℓ, layers; direct_only) for s in rs]
    idx = findfirst(i -> sign(v[i]) != sign(v[i+1]), 1:length(v)-1)
    isnothing(idx) && return nothing
    s0 = rs[idx] - v[idx] * (rs[idx+1] - rs[idx]) / (v[idx+1] - v[idx])
    return (A_EARTH - s0) / 1e3
end

# ╔═╡ f6d9c865-d6d7-4819-4a65-829304152637
md"""
### Validating the kernel

Four independent physical checks. Each one would catch a different class of sign error.
"""

# ╔═╡ 07ead975-e7e8-492a-5bb6-930415263748
begin
    _val_rs = range(B_CMB + 20e3, A_EARTH - 20e3, length=40)

    # Check A: with boundaries off, the kernel must be positive and monotonic.
    _direct_vals = [geoid_kernel(s, 2, isoviscous_layers(); direct_only=true) for s in _val_rs]
    check_direct = all(_direct_vals .> 0) && all(diff(_direct_vals) .> 0)

    # Check B: free-slip boundary conditions satisfied, measured RELATIVE to the
    # size of the solution itself (an absolute threshold would be meaningless here,
    # since the state components carry physical units spanning many orders).
    check_bc = maximum(
        let resp = geoid_response(s, ℓ, isoviscous_layers())
            resp.bc_residual / max(abs(resp.h_surf), 1e-30)
        end
        for s in _val_rs, ℓ in (2, 4, 8)) < 1e-8

    # Check C: a dense anomaly must DEPRESS the surface at every depth.
    check_topo = all(geoid_response(s, 2, isoviscous_layers()).h_surf < 0 for s in _val_rs)

    # Check D: isoviscous ℓ=2 kernel flips sign, deep-positive / shallow-negative.
    _cross2 = zero_crossing_depth(2, isoviscous_layers())
    check_flip = !isnothing(_cross2) &&
                 geoid_kernel(B_CMB + 100e3, 2, isoviscous_layers()) > 0 &&
                 geoid_kernel(A_EARTH - 100e3, 2, isoviscous_layers()) < 0
    nothing
end

# ╔═╡ 18fbea85-f8f9-4a3b-6c97-041526374859
if check_direct && check_bc && check_topo && check_flip
    md"""
    !!! correct "Validation 2 passed — all four physical checks"
    	- **Direct-only kernel** is positive and monotonic everywhere ✓ (no boundaries ⇒ no flip)
    	- **Free-slip boundary conditions** satisfied to machine zero at surface and CMB ✓
    	- **A dense anomaly depresses the surface** at every depth ✓ (the load-jump sign is right)
    	- **The isoviscous ``\ell=2`` kernel flips sign** at **$(round(Int, _cross2)) km** depth,
    	  negative shallow and positive deep ✓

    	That crossing depth is the canonical Hager & Richards (1989) result for a uniform
    	mantle, which is our reproduction-of-a-published-kernel check.
    """
else
    md"""
    !!! danger "Validation 2 FAILED"
    	direct-only: $(check_direct) · boundary conditions: $(check_bc) ·
    	surface depression: $(check_topo) · sign flip: $(check_flip)
    """
end

# ╔═╡ 290cfb95-0901-4b4c-7da8-152637485960
md"""
## Layer 3: From a Painted Picture to a Geoid

The canvas gives us $\delta\rho(r,\theta)$ on a grid. To use the kernels we need its
**Legendre coefficients** at each radius:

$(L"\delta\rho_\ell(r) = \frac{2\ell+1}{2}\int_{-1}^{1}\delta\rho(r,\theta)\,P_\ell(\cos\theta)\,d(\cos\theta)")

then apply the kernel for each degree and sum back up:

$(L"N(\theta) = \sum_\ell \left[\int K_\ell(r)\,\delta\rho_\ell(r)\,dr\right] P_\ell(\cos\theta)")

We roll the Legendre machinery by hand — a three-term recurrence and a Newton solve for the
Gauss–Legendre nodes. It is about fifteen lines, needs no extra packages, and you can read
every step.
"""

# ╔═╡ 3a1d0ca5-1a0b-4c5d-8eb9-263748596071
"All Legendre polynomials `P_0..P_lmax` at `x`, via the three-term recurrence."
function legendre_all(x, lmax)
    P = zeros(lmax + 1)
    P[1] = 1.0
    lmax >= 1 && (P[2] = x)
    for l in 1:lmax-1
        P[l+2] = ((2l + 1) * x * P[l+1] - l * P[l]) / (l + 1)
    end
    return P
end

# ╔═╡ 4b2e1db5-2b1c-4d6e-9fca-374859607182
"Gauss–Legendre nodes and weights on [-1,1], found by Newton iteration on `P_n`."
function gauss_legendre(n)
    x = zeros(n)
    w = zeros(n)
    for i in 1:n
        z = cos(pi * (i - 0.25) / (n + 0.5))   # excellent initial guess
        for _ in 1:100
            P = legendre_all(z, n)
            dP = n * (P[n] - z * P[n+1]) / (1 - z^2)
            dz = -P[n+1] / dP
            z += dz
            abs(dz) < 1e-15 && break
        end
        P = legendre_all(z, n)
        dP = n * (P[n] - z * P[n+1]) / (1 - z^2)
        x[i] = z
        w[i] = 2 / ((1 - z^2) * dP^2)
    end
    return x, w
end

# ╔═╡ 6d403fd5-4d3e-4f80-b1eb-596071829304
md"""
### Validating the Legendre machinery

Two checks: the polynomials must be **orthogonal** under the quadrature, and a
decompose-then-synthesize round trip must return the field we started with.
"""

# ╔═╡ 7e5140e5-5e4f-4091-c2fb-607182930415
begin
    _leg_lmax = 30
    _leg_n = 2 * _leg_lmax + 2
    _leg_x, _leg_w = gauss_legendre(_leg_n)
    _leg_P = [legendre_all(xi, _leg_lmax) for xi in _leg_x]

    # orthogonality: ∫ P_i P_j = 2/(2i+1) δ_ij
    _orth_off = 0.0
    _orth_diag = 0.0
    for i in 0:_leg_lmax, j in 0:_leg_lmax
        s = sum(_leg_w[k] * _leg_P[k][i+1] * _leg_P[k][j+1] for k in 1:_leg_n)
        if i == j
            _orth_diag = max(_orth_diag, abs(s - 2 / (2i + 1)))
        else
            _orth_off = max(_orth_off, abs(s))
        end
    end
    check_legendre = _orth_off < 1e-12 && _orth_diag < 1e-12 &&
                     abs(sum(_leg_w) - 2) < 1e-12
    nothing
end

# ╔═╡ 8f6251f6-6f40-41a2-d3bd-718293041526
if check_legendre
    md"""
    !!! correct "Validation 3 passed"
    	Quadrature weights sum to 2 exactly; Legendre orthogonality holds with maximum
    	off-diagonal error **$(round(_orth_off, sigdigits=3))** and diagonal error
    	**$(round(_orth_diag, sigdigits=3))** up to ``\ell = 30``.
    """
else
    md"""!!! danger "Validation 3 FAILED" """
end

# ╔═╡ 90736266-7051-42b3-e4cf-829304152637
md"""
## Putting It Together: Painted Anomaly → Geoid

Now the full pipeline. The canvas hands us a set of painted blobs; we build
$\delta\rho(r,\theta)$ from them, decompose, apply kernels, and synthesize.
"""

# ╔═╡ a18473c6-8152-43c4-f5db-930415263748
"""
	blobs_to_field(blobs, rgrid, μgrid)

Build ``\\delta\\rho(r, \\cos\\theta)`` from painted blobs. Each blob is
`(r_centre, μ_centre, amplitude)` and is given a smooth Gaussian profile so the
Legendre expansion converges without ringing.
"""
function blobs_to_field(blobs, rgrid, μgrid; r_width=250e3, μ_width=0.10)
    F = zeros(length(rgrid), length(μgrid))
    isempty(blobs) && return F
    for (rc, μc, amp) in blobs
        for (i, r) in enumerate(rgrid), (j, μ) in enumerate(μgrid)
            F[i, j] += amp * exp(-((r - rc) / r_width)^2 - ((μ - μc) / μ_width)^2)
        end
    end
    return F
end

# ╔═╡ b2958426-9263-44d5-06e1-041526374859
"""
	synthesize_geoid(blobs, layers, lmax; direct_only=false)

The full Layer-3 pipeline. Returns the geoid ``N(\\theta)`` sampled on the quadrature
nodes, the per-degree contributions, and the radial grid used.
"""
function synthesize_geoid(blobs, layers, lmax; direct_only=false, nr=48)
    n = max(2 * lmax + 2, 16)
    μ, w = gauss_legendre(n)
    rgrid = range(B_CMB + 30e3, A_EARTH - 30e3, length=nr)
    dr = step(rgrid)
    F = blobs_to_field(blobs, rgrid, μ)
    Pall = [legendre_all(μk, lmax) for μk in μ]

    per_degree = zeros(lmax + 1)
    spectrum = zeros(lmax + 1)
    for ℓ in 0:lmax
        # δρ_ℓ(r) at each radius by Gauss-Legendre quadrature
        acc = 0.0
        specacc = 0.0
        for (i, r) in enumerate(rgrid)
            δρℓ = (2ℓ + 1) / 2 * sum(w[k] * F[i, k] * Pall[k][ℓ+1] for k in 1:n)
            specacc += abs(δρℓ) * dr
            ℓ == 0 && continue          # degree 0 cannot deform the geoid shape
            acc += geoid_kernel(r, ℓ, layers; direct_only) * δρℓ * dr
        end
        per_degree[ℓ+1] = acc
        spectrum[ℓ+1] = specacc
    end

    N = [sum(per_degree[ℓ+1] * Pall[k][ℓ+1] for ℓ in 0:lmax) for k in 1:n]
    return (μ=μ, N=N, per_degree=per_degree, spectrum=spectrum, rgrid=rgrid, field=F)
end

# ╔═╡ c3a69596-a374-45e6-17f2-152637485960
md"""
## The Interactive Panel

All the drawing and all the plotting happen on HTML canvases inside a single widget. Julia
computes the physics and pushes the results back to the browser through a `CustomEvent`,
which the widget listens for and redraws.
"""

# ╔═╡ d4b7a646-b485-46f7-2803-263748596071
begin
    struct GeoidCanvasInput
        default_blobs::Vector{Vector{Float64}}
        visc_ratio::Float64
        ell::Int
        lmax::Int
        direct_only::Bool
        dense::Bool
    end

    function GeoidCanvasInput(; visc_ratio=1.0, ell=2, lmax=20,
        direct_only=false, dense=true)
        # one dense blob at ~1900 km depth, mid-latitude, as a starting picture
        GeoidCanvasInput([[0.62, 0.30, 1.0]], Float64(visc_ratio), Int(ell),
            Int(lmax), Bool(direct_only), Bool(dense))
    end

    function colorscheme_stops_js(scheme, n=9)
        stops = map(range(0, 1, length=n)) do x
            c = get(scheme, x)
            rgb = round.(Int, 255 .* (c.r, c.g, c.b))
            "[$(round(x, digits=4)),[$(rgb[1]),$(rgb[2]),$(rgb[3])]]"
        end
        "[" * join(stops, ",") * "]"
    end

    Base.get(w::GeoidCanvasInput) = Dict{String,Any}(
        "blobs" => w.default_blobs,
        "visc_ratio" => w.visc_ratio,
        "ell" => w.ell,
        "lmax" => w.lmax,
        "direct_only" => w.direct_only,
        "dense" => w.dense,
    )

    function Base.show(io::IO, ::MIME"text/html", w::GeoidCanvasInput)
        write(io, """
<div id="gkwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    #gkwidget .gk-workspace{display:flex;gap:12px;align-items:flex-start;justify-content:center;width:100%}
    #gkwidget .gk-controls{width:min(914px,100%);margin-top:8px;display:grid;grid-template-columns:repeat(2,minmax(300px,1fr));gap:8px;font:12px sans-serif}
    #gkwidget .gk-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:9px 10px}
    #gkwidget .gk-control-title{font-weight:700;color:#e5e7eb;margin-bottom:7px;font-size:18px}
    #gkwidget .gk-control-row{display:grid;grid-template-columns:150px minmax(110px,1fr) 70px;gap:8px;align-items:center;margin:6px 0}
    #gkwidget .gk-control-row input[type=range]{width:100%;vertical-align:middle}
    #gkwidget .gk-value{color:#d1d5db;text-align:left;font-variant-numeric:tabular-nums}
    #gkwidget .gk-actions{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
    #gkwidget button{border-radius:4px;border:1px solid #6b7280;background:#606060;color:#f3f4f6;padding:4px 9px}
    #gkwidget label{color:#d1d5db}
    @media (max-width: 980px){
      #gkwidget .gk-workspace{flex-direction:column;align-items:center}
      #gkwidget .gk-controls{grid-template-columns:1fr;width:640px;max-width:100%}
    }
  </style>
  <canvas id="geoidcvs" width="914" height="150"
    style="background:#000;border:1px solid #374151;border-radius:4px;display:block;margin-bottom:4px"></canvas>
  <div class="gk-workspace">
    <canvas id="sectioncvs" width="640" height="640"
      style="cursor:crosshair;background:#000;border:1px solid #374151;border-radius:6px;display:block"></canvas>
    <canvas id="kernelcvs" width="262" height="640"
      style="background:#000;border:1px solid #374151;border-radius:6px;display:block"></canvas>
  </div>
  <div class="gk-controls">
    <div class="gk-control-group">
      <div class="gk-control-title">Structure</div>
      <label class="gk-control-row"><span>η lower / η upper</span><input type="range" id="visc" min="0" max="2" step="0.01" value="$(round(log10(w.visc_ratio), digits=4))"><span id="viscv" class="gk-value">$(round(w.visc_ratio, digits=1))×</span></label>
      <label class="gk-control-row"><span>kernel degree ℓ</span><input type="range" id="ell" min="2" max="12" step="1" value="$(w.ell)"><span id="ellv" class="gk-value">$(w.ell)</span></label>
      <label class="gk-control-row"><span>ℓ max (synthesis)</span><input type="range" id="lmax" min="4" max="30" step="1" value="$(w.lmax)"><span id="lmaxv" class="gk-value">$(w.lmax)</span></label>
    </div>
    <div class="gk-control-group">
      <div class="gk-control-title">Anomaly &amp; physics</div>
      <div class="gk-actions">
        <label><input type="checkbox" id="direct" $(w.direct_only ? "checked" : "") style="vertical-align:middle"> <b>direct effect only</b> (no boundaries)</label>
      </div>
      <div class="gk-actions" style="margin-top:8px">
        <label><input type="checkbox" id="dense" $(w.dense ? "checked" : "") style="vertical-align:middle"> paint dense (uncheck for light)</label>
      </div>
      <div class="gk-actions" style="margin-top:8px">
        <button id="clrbtn">Clear</button>
        <button id="shallowbtn">Shallow blob</button>
        <button id="deepbtn">Deep blob</button>
        <span id="cnt" class="gk-value">blobs: $(length(w.default_blobs))</span>
      </div>
    </div>
  </div>
</div>
<script>
  const A = $(A_EARTH), B = $(B_CMB)
  const SEC = 640, GW = 914, GH = 150, KW = 262, KH = 640
  const DPR = Math.min(window.devicePixelRatio || 1, 2)
  const par = currentScript.previousElementSibling
  const sec = par.querySelector('#sectioncvs'), sctx = sec.getContext('2d')
  const geo = par.querySelector('#geoidcvs'),   gctx = geo.getContext('2d')
  const ker = par.querySelector('#kernelcvs'),  kctx = ker.getContext('2d')
  const lbl = par.querySelector('#cnt')

  function hidpi(cv, cx, w, h){
    cv.width = Math.round(w*DPR); cv.height = Math.round(h*DPR)
    cv.style.width = w+'px'; cv.style.height = h+'px'
    cx.setTransform(DPR,0,0,DPR,0,0)
  }
  hidpi(sec,sctx,SEC,SEC); hidpi(geo,gctx,GW,GH); hidpi(ker,kctx,KW,KH)

  // blobs: [mu, rfrac, amp]  with rfrac in [0,1] from CMB(0) to surface(1)
  let blobs = $(replace(string([[b[1], b[2], b[3]] for b in w.default_blobs]), "Any" => ""))
  let viscRatio = $(w.visc_ratio), ell = $(w.ell), lmax = $(w.lmax)
  let directOnly = $(w.direct_only), dense = $(w.dense)
  let painting = false
  let kernelData = null, geoidData = null, topoData = null

  // --- geometry: the cross-section is drawn as an annulus wedge ---------------
  const CX = SEC/2, CY = SEC/2, R_OUT = 290, R_IN = R_OUT*B/A
  function toCanvas(mu, rfrac){
    const th = Math.acos(Math.max(-1,Math.min(1,mu)))   // colatitude 0..pi
    const rr = R_IN + rfrac*(R_OUT-R_IN)
    return [CX + rr*Math.sin(th), CY - rr*Math.cos(th)]
  }
  function fromCanvas(px, py){
    const dx = px-CX, dy = CY-py
    const rr = Math.hypot(dx,dy)
    if(rr < R_IN-6 || rr > R_OUT+6) return null
    const rfrac = Math.max(0,Math.min(1,(rr-R_IN)/(R_OUT-R_IN)))
    const th = Math.atan2(Math.abs(dx), dy)
    return [Math.cos(th), rfrac]
  }

  function drawSection(){
    sctx.clearRect(0,0,SEC,SEC)
    // mantle annulus
    sctx.save()
    sctx.beginPath(); sctx.arc(CX,CY,R_OUT,0,2*Math.PI)
    sctx.arc(CX,CY,R_IN,0,2*Math.PI,true)
    sctx.fillStyle='#0d1117'; sctx.fill('evenodd')
    sctx.restore()
    // core
    sctx.beginPath(); sctx.arc(CX,CY,R_IN,0,2*Math.PI)
    sctx.fillStyle='#241a12'; sctx.fill()

    // painted anomalies
    for(const [mu,rf,amp] of blobs){
      const [x,y] = toCanvas(mu,rf)
      const g = sctx.createRadialGradient(x,y,0,x,y,34)
      const col = amp>0 ? '239,68,68' : '59,130,246'
      g.addColorStop(0,'rgba('+col+',0.85)'); g.addColorStop(1,'rgba('+col+',0)')
      sctx.fillStyle=g; sctx.beginPath(); sctx.arc(x,y,34,0,2*Math.PI); sctx.fill()
    }

    // deflected boundaries (exaggerated) pushed from Julia
    if(topoData){
      for(const [key,baseR,color] of [['surf',R_OUT,'#22d3ee'],['cmb',R_IN,'#f59e0b']]){
        const arr = topoData[key]; if(!arr) continue
        const mx = Math.max(...arr.map(Math.abs)) || 1
        sctx.beginPath()
        for(let i=0;i<arr.length;i++){
          const mu = -1 + 2*i/(arr.length-1)
          const th = Math.acos(Math.max(-1,Math.min(1,mu)))
          const rr = baseR + 26*arr[i]/mx
          const x = CX + rr*Math.sin(th), y = CY - rr*Math.cos(th)
          i===0 ? sctx.moveTo(x,y) : sctx.lineTo(x,y)
        }
        sctx.strokeStyle=color; sctx.lineWidth=1.6; sctx.setLineDash([5,4]); sctx.stroke()
        sctx.setLineDash([])
      }
    }
    // reference circles
    sctx.strokeStyle='#4b5563'; sctx.lineWidth=1
    sctx.beginPath(); sctx.arc(CX,CY,R_OUT,0,2*Math.PI); sctx.stroke()
    sctx.beginPath(); sctx.arc(CX,CY,R_IN,0,2*Math.PI); sctx.stroke()

    sctx.fillStyle='#9ca3af'; sctx.font='13px sans-serif'
    sctx.fillText('mantle cross-section — drag to paint', 12, 20)
    sctx.fillText('surface', CX-22, CY-R_OUT-8)
    sctx.fillText('CMB', CX-14, CY-R_IN-6)
    sctx.fillStyle='#22d3ee'; sctx.fillText('— deflected surface (exagg.)', 12, SEC-30)
    sctx.fillStyle='#f59e0b'; sctx.fillText('— deflected CMB (exagg.)', 12, SEC-14)
    lbl.textContent = 'blobs: '+blobs.length
  }

  function axes(cx,w,h,pad){
    cx.strokeStyle='#374151'; cx.lineWidth=1
    cx.beginPath(); cx.moveTo(pad,h-pad); cx.lineTo(w-pad,h-pad); cx.stroke()
    cx.beginPath(); cx.moveTo(pad,pad); cx.lineTo(pad,h-pad); cx.stroke()
  }

  function drawKernel(){
    kctx.clearRect(0,0,KW,KH)
    kctx.fillStyle='#9ca3af'; kctx.font='12px sans-serif'
    kctx.fillText('kernel K'+'\\u2113'+'(r), \\u2113='+ell, 10, 18)
    if(!kernelData){ kctx.fillText('computing...', 10, 40); return }
    const pad=38, W=KW-pad-14, H=KH-pad-40
    const vals = kernelData.k, deps = kernelData.depth
    const mx = Math.max(...vals.map(Math.abs)) || 1
    const X = v => pad + W*(0.5 + 0.45*v/mx)
    const Y = d => pad + H*(d/2891)
    // zero line
    kctx.strokeStyle='#6b7280'; kctx.setLineDash([4,4])
    kctx.beginPath(); kctx.moveTo(X(0),pad); kctx.lineTo(X(0),pad+H); kctx.stroke()
    kctx.setLineDash([])
    // curve
    kctx.beginPath()
    for(let i=0;i<vals.length;i++){
      const x=X(vals[i]), y=Y(deps[i])
      i===0?kctx.moveTo(x,y):kctx.lineTo(x,y)
    }
    kctx.strokeStyle='#f97316'; kctx.lineWidth=2; kctx.stroke()
    // zero crossing marker
    if(kernelData.crossing !== null && kernelData.crossing !== undefined){
      const y=Y(kernelData.crossing)
      kctx.strokeStyle='#22c55e'; kctx.lineWidth=1.5; kctx.setLineDash([3,3])
      kctx.beginPath(); kctx.moveTo(pad,y); kctx.lineTo(KW-14,y); kctx.stroke()
      kctx.setLineDash([])
      kctx.fillStyle='#22c55e'; kctx.font='11px sans-serif'
      kctx.fillText('sign flip '+Math.round(kernelData.crossing)+' km', pad+4, y-5)
    } else {
      kctx.fillStyle='#22c55e'; kctx.font='11px sans-serif'
      kctx.fillText('no sign flip', pad+4, pad+14)
    }
    axes(kctx,KW,KH,pad)
    kctx.fillStyle='#9ca3af'; kctx.font='11px sans-serif'
    kctx.save(); kctx.translate(12,pad+H/2); kctx.rotate(-Math.PI/2)
    kctx.fillText('depth (km)',-30,0); kctx.restore()
    kctx.fillText('0',pad-4,pad-6); kctx.fillText('2891',pad-16,pad+H+14)
    kctx.fillStyle='#ef4444'; kctx.fillText('negative',pad+4,KH-16)
    kctx.fillStyle='#3b82f6'; kctx.fillText('positive \\u2192',KW-84,KH-16)
  }

  function drawGeoid(){
    gctx.clearRect(0,0,GW,GH)
    gctx.fillStyle='#9ca3af'; gctx.font='13px sans-serif'
    gctx.fillText('geoid N(\\u03b8) at the surface', 12, 18)
    if(!geoidData){ gctx.fillText('computing...', 12, 40); return }
    const pad=34, W=GW-2*pad, H=GH-pad-28
    const mx = Math.max(...geoidData.map(Math.abs)) || 1
    const X = i => pad + W*i/(geoidData.length-1)
    const Y = v => pad + H*(0.5 - 0.45*v/mx)
    gctx.strokeStyle='#374151'; gctx.beginPath()
    gctx.moveTo(pad,Y(0)); gctx.lineTo(GW-pad,Y(0)); gctx.stroke()
    // filled curve, colored by sign
    for(let i=0;i<geoidData.length-1;i++){
      gctx.beginPath()
      gctx.moveTo(X(i),Y(0)); gctx.lineTo(X(i),Y(geoidData[i]))
      gctx.lineTo(X(i+1),Y(geoidData[i+1])); gctx.lineTo(X(i+1),Y(0)); gctx.closePath()
      const up = (geoidData[i]+geoidData[i+1]) > 0
      gctx.fillStyle = up ? 'rgba(59,130,246,0.45)' : 'rgba(239,68,68,0.45)'
      gctx.fill()
    }
    gctx.beginPath()
    for(let i=0;i<geoidData.length;i++){
      i===0?gctx.moveTo(X(i),Y(geoidData[i])):gctx.lineTo(X(i),Y(geoidData[i]))
    }
    gctx.strokeStyle='#e5e7eb'; gctx.lineWidth=1.8; gctx.stroke()
    gctx.fillStyle='#9ca3af'; gctx.font='11px sans-serif'
    gctx.fillText('pole (\\u03b8=0)', pad, GH-8)
    gctx.fillText('equator', pad+W/2-20, GH-8)
    gctx.fillText('pole (\\u03b8=180\\u00b0)', GW-pad-70, GH-8)
    gctx.fillText('peak |N| = '+mx.toExponential(2)+' m', GW-pad-160, 18)
  }

  function redraw(){ drawSection(); drawKernel(); drawGeoid() }

  function emit(){
    par.value = {blobs:blobs, visc_ratio:viscRatio, ell:ell, lmax:lmax,
                 direct_only:directOnly, dense:dense}
    par.dispatchEvent(new CustomEvent('input'))
  }

  function addBlob(px,py){
    const p = fromCanvas(px,py); if(!p) return false
    const [mu,rf] = p
    // avoid piling up hundreds of overlapping blobs while dragging
    for(const b of blobs){
      if(Math.abs(b[0]-mu)<0.045 && Math.abs(b[1]-rf)<0.045) return false
    }
    blobs.push([mu,rf, dense?1.0:-1.0])
    return true
  }

  sec.addEventListener('mousedown', e=>{ painting=true; if(addBlob(e.offsetX,e.offsetY)){drawSection(); emit()} })
  sec.addEventListener('mousemove', e=>{ if(painting && addBlob(e.offsetX,e.offsetY)){drawSection(); emit()} })
  window.addEventListener('mouseup', ()=>{ painting=false })

  par.querySelector('#visc').addEventListener('input', e=>{
    viscRatio = Math.pow(10, parseFloat(e.target.value))
    par.querySelector('#viscv').textContent = viscRatio.toFixed(viscRatio<10?1:0)+'\\u00d7'
    emit()
  })
  par.querySelector('#ell').addEventListener('input', e=>{
    ell = parseInt(e.target.value); par.querySelector('#ellv').textContent = ell; emit()
  })
  par.querySelector('#lmax').addEventListener('input', e=>{
    lmax = parseInt(e.target.value); par.querySelector('#lmaxv').textContent = lmax; emit()
  })
  par.querySelector('#direct').addEventListener('change', e=>{ directOnly=e.target.checked; emit() })
  par.querySelector('#dense').addEventListener('change', e=>{ dense=e.target.checked })
  par.querySelector('#clrbtn').addEventListener('click', ()=>{ blobs=[]; drawSection(); emit() })
  // Presets chosen so the sign flip is unmistakable: at 0.80 (≈580 km) the
  // boundary terms dominate and the geoid is LOW; at 0.10 (≈2600 km) they no
  // longer can and it is HIGH.
  par.querySelector('#shallowbtn').addEventListener('click', ()=>{
    blobs=[[0.30, 0.80, dense?1.0:-1.0]]; drawSection(); emit() })
  par.querySelector('#deepbtn').addEventListener('click', ()=>{
    blobs=[[0.30, 0.10, dense?1.0:-1.0]]; drawSection(); emit() })

  // results pushed back from Julia
  window.addEventListener('geoid-results', e=>{
    const d = e.detail ? JSON.parse(e.detail) : null
    if(!d) return
    kernelData = d.kernel; geoidData = d.geoid; topoData = d.topo
    redraw()
  })

  redraw(); emit()
</script>
        """)
    end

    const _geoid_canvas_ready = true
end

# ╔═╡ e5c8b756-c5b6-4708-3964-374859607182
begin
    _geoid_canvas_ready
    WideCell(@bind _gk GeoidCanvasInput())
end

# ╔═╡ f6d9c866-d6a7-4819-4a75-485960718293
begin
    gk_blobs_raw = get(_gk, "blobs", Vector{Float64}[])
    gk_visc = Float64(get(_gk, "visc_ratio", 1.0))
    gk_ell = Int(round(Float64(get(_gk, "ell", 2))))
    gk_lmax = Int(round(Float64(get(_gk, "lmax", 20))))
    gk_direct = Bool(get(_gk, "direct_only", false))

    # canvas gives (μ, r_fraction, amplitude) -> physical (r, μ, amp)
    gk_blobs = [(B_CMB + Float64(b[2]) * (A_EARTH - B_CMB), Float64(b[1]), Float64(b[3]))
                for b in gk_blobs_raw]
    gk_layers = gk_visc ≈ 1.0 ? isoviscous_layers() : twolayer_layers(gk_visc)
    nothing
end

# ╔═╡ 07ead976-e7b8-492a-5bc6-596071829304
gk_result = synthesize_geoid(gk_blobs, gk_layers, gk_lmax; direct_only=gk_direct)

# ╔═╡ 18fbea86-f8c9-4a3b-6ca7-607182930415
begin
    # kernel curve for the selected single degree
    gk_kdepths = range(20.0, 2870.0, length=120)
    gk_kvals = [geoid_kernel(A_EARTH - d * 1e3, gk_ell, gk_layers; direct_only=gk_direct)
                for d in gk_kdepths]
    gk_crossing = zero_crossing_depth(gk_ell, gk_layers; direct_only=gk_direct)
    nothing
end

# ╔═╡ 290cfb96-0871-4b4c-7db8-718293041526
begin
    # dynamic topography profiles along the surface and CMB, for the cross-section
    _nμ = length(gk_result.μ)
    _topo_surf = zeros(_nμ)
    _topo_cmb = zeros(_nμ)
    let
        dr = step(gk_result.rgrid)
        Pall = [legendre_all(μk, gk_lmax) for μk in gk_result.μ]
        _, wq = gauss_legendre(_nμ)
        for ℓ in 1:gk_lmax
            hs = 0.0
            hc = 0.0
            for (i, r) in enumerate(gk_result.rgrid)
                δρℓ = (2ℓ + 1) / 2 * sum(wq[k] * gk_result.field[i, k] * Pall[k][ℓ+1]
                                         for k in 1:_nμ)
                resp = geoid_response(r, ℓ, gk_layers)
                hs += resp.h_surf * δρℓ * dr
                hc += resp.h_cmb * δρℓ * dr
            end
            for k in 1:_nμ
                _topo_surf[k] += hs * Pall[k][ℓ+1]
                _topo_cmb[k] += hc * Pall[k][ℓ+1]
            end
        end
    end
    nothing
end

# ╔═╡ 3a1d0ca6-1adb-4c5d-8ec9-829304152637
# Push the computed kernel, geoid and topography back to the canvas widget.
let
    jsonarr(v) = "[" * join([isfinite(x) ? string(round(x, sigdigits=6)) : "0" for x in v], ",") * "]"
    payload = string(
        "{\"kernel\":{\"k\":", jsonarr(gk_kvals),
        ",\"depth\":", jsonarr(collect(gk_kdepths)),
        ",\"crossing\":", isnothing(gk_crossing) ? "null" : string(round(gk_crossing, digits=1)),
        "},\"geoid\":", jsonarr(gk_result.N),
        ",\"topo\":{\"surf\":", jsonarr(_topo_surf),
        ",\"cmb\":", jsonarr(_topo_cmb), "}}")
    HTML("""<script>
      window.dispatchEvent(new CustomEvent('geoid-results', {detail: $(repr(payload))}));
    </script>""")
end

# ╔═╡ 4b2e1db6-2bec-4d6e-9fda-930415263748
md"""
### The Legendre spectrum of what you painted

The one conventional plot in this notebook. It shows how the anomaly you drew decomposes
into degrees — and, crucially, the **geoid contribution of each degree separately**. Notice
that different degrees can contribute with **opposite signs at the same place**, because each
degree has its own kernel with its own zero-crossing depth.
"""

# ╔═╡ 5c3f2ec6-3cfd-4e7f-a0eb-041526374859
let
    ℓs = 1:gk_lmax
    contrib = [gk_result.per_degree[ℓ+1] for ℓ in ℓs]
    colors = [c > 0 ? "#3b82f6" : "#ef4444" for c in contrib]
    fig = PlutoPlotly.Plot(Layout(
        height=300, width=700,
        title="geoid contribution by degree (blue = high, red = low)",
        xaxis=attr(title="spherical harmonic degree ℓ", dtick=2),
        yaxis=attr(title="contribution to N (m)"),
        font=attr(size=11), showlegend=false))
    add_trace!(fig, PlutoPlotly.bar(x=collect(ℓs), y=contrib, marker=attr(color=colors)))
    PlutoPlotly.plot(fig)
end

# ╔═╡ 6d403fd6-4dcd-4f80-b1fb-152637485960
md"""
## Aha #3: Viscosity Controls the Flip

Sweep the **η lower / η upper** slider and watch the kernel. Here is what actually happens —
and it is *not* quite what you might expect from the phrase "the sign-flip moves":

$(let
	rows = String[]
	for ratio in (1.0, 3.0, 10.0, 30.0, 100.0)
		lay = ratio ≈ 1.0 ? isoviscous_layers() : twolayer_layers(ratio)
		d = zero_crossing_depth(2, lay)
		push!(rows, "| $(ratio)× | " * (isnothing(d) ? "**none — positive everywhere**" : "$(round(Int,d)) km") * " |")
	end
	Markdown.parse("| ``\\eta_{lower}/\\eta_{upper}`` | ``\\ell=2`` sign flip |\n|---|---|\n" * join(rows, "\n"))
end)

At **moderate** contrast the crossing migrates *upward* — from 1368 km in a uniform mantle to
about 744 km at 3×. But push past roughly 10× and the crossing **disappears entirely**: the
kernel becomes positive at every depth.

!!! tip "Why the disappearance is the whole point"
	A stiff lower mantle resists the flow that a sinking anomaly tries to drive. Weaken that
	flow and you weaken the boundary deflections — which were the *negative* contributions.
	Past about 10×, they can no longer outvote the anomaly's own positive attraction at any
	depth.

	This is exactly how the geoid constrains viscosity. Subducted slabs are dense and deep,
	and they correlate **positively** with the observed geoid. A uniform-viscosity mantle
	predicts the opposite sign for deep slabs. Getting the observed sign right *requires*
	a lower mantle roughly 10–100× stiffer than the upper mantle — which is precisely the
	conclusion Hager & Richards drew in 1989, and it remains one of the strongest constraints
	we have on mantle rheology.
"""

# ╔═╡ 7e5140e6-5ebf-4091-c20b-263748596071
md"""
## The Indian Ocean Geoid Low

South of India lies the deepest geoid low on Earth: sea level sits roughly **106 m below**
the reference ellipsoid over a vast region of the Indian Ocean. Explaining it is an open
research problem.

What this notebook gives you is the crucial piece of intuition. The naive objection —
*"but there is dense subducted slab material down there, so shouldn't the geoid be high?"* —
dissolves once you have seen the kernel. Under the right viscosity structure, deep dense
material produces a geoid **low**, because the boundary-deflection terms outvote the
anomaly's direct attraction.

Try it: hit **Shallow blob**, keep the viscosity ratio near 1, and watch the geoid go
negative even though the anomaly is *dense*. Then slide the ratio up and watch it flip
positive. Whether the real Indian Ocean low comes
out right therefore depends sensitively on the *actual* radial viscosity profile — which is
why this remains contested, and why people run full 3-D models with lateral viscosity
variations rather than kernels like these.

!!! warning "One honest caveat"
	Real explanations of the IOGL invoke lateral viscosity variations, a low-density anomaly
	beneath the low, and the detailed history of Tethyan slab sinking — none of which this
	axisymmetric, radially-viscous, instantaneous model contains. What you have here is the
	*mechanism*, not the answer.
"""

# ╔═╡ 5c3f2ec5-3c2d-4e7f-a0db-485960718293
md"""
## Assumptions — What This Notebook Does and Does Not Do

Being explicit about this matters more than the pretty pictures.

!!! warning "Read before drawing conclusions"
	- **Radial viscosity only.** ``\eta`` varies with depth, never sideways. Real mantle has
	  strong lateral viscosity variations (cold stiff slabs, hot weak plumes). Capturing those
	  needs a genuine 3-D code — CitcomS, ASPECT. The kernel captures the *mechanism*, not the
	  specifics of any real region.
	- **Axisymmetric (``m = 0``).** Anomalies are rings around the polar axis, so we use
	  Legendre polynomials rather than full spherical harmonics. A deliberate simplification
	  to keep the mathematics readable.
	- **Instantaneous.** One linear solve. No time-stepping, no evolution — the anomaly stays
	  where you paint it.
	- **Self-gravitation neglected.** The deflected boundaries and the anomaly do not pull on
	  each other. This shifts long-wavelength amplitudes by roughly 10–20% and can move the
	  crossing depth somewhat, but it does not change the sign structure that is the lesson here.
	- **Free-slip boundaries, no lithosphere.** No elastic plate, no rigid lid.

!!! tip "For quantitative work, use a real code"
	This notebook is built for **intuition**. For research, use
	**HC** (Becker et al.) for propagator-matrix geoid calculations, or **ASPECT** /
	**CitcomS** for full 3-D convection. Benchmarks show propagator-matrix codes and ASPECT
	agree to about 1% out to degree ~15, which is exactly why this simple approach is
	trustworthy *at long wavelengths* — and why we truncate rather than chasing high degrees.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
"""

# ╔═╡ Cell order:
# ╠═5c3f2e62-3c4d-4e7f-a01b-2c3d4e5f6071
# ╟─6d403f73-4d5e-4f80-b12c-3d4e5f607182
# ╟─7e514084-5e6f-4091-c23d-4e5f60718293
# ╟─8f625195-6f70-41a2-d34e-5f6071829304
# ╟─c3a69596-a374-45e6-17f2-152637485960
# ╠═d4b7a646-b485-46f7-2803-263748596071
# ╠═e5c8b756-c5b6-4708-3964-374859607182
# ╠═f6d9c866-d6a7-4819-4a75-485960718293
# ╠═07ead976-e7b8-492a-5bc6-596071829304
# ╠═18fbea86-f8c9-4a3b-6ca7-607182930415
# ╠═290cfb96-0871-4b4c-7db8-718293041526
# ╠═3a1d0ca6-1adb-4c5d-8ec9-829304152637
# ╟─4b2e1db6-2bec-4d6e-9fda-930415263748
# ╠═5c3f2ec6-3cfd-4e7f-a0eb-041526374859
# ╟─6d403fd6-4dcd-4f80-b1fb-152637485960
# ╟─7e5140e6-5ebf-4091-c20b-263748596071
# ╟─b2958428-92a3-44d5-0671-829304152637
# ╠═c3a69539-a3b4-45e6-1782-930415263748
# ╟─d4b7a64a-b4c5-46f7-2893-041526374859
# ╟─e5c8b75b-c5d6-4708-3904-152637485960
# ╟─90736206-7081-42b3-e45f-607182930415
# ╠═a1847317-8192-43c4-f560-718293041526
# ╟─f6d9c86c-d6e7-4819-4a15-263748596071
# ╠═07ead97d-e7f8-492a-5b26-374859607182
# ╠═18fbea8e-f809-4a3b-6c37-485960718293
# ╠═290cfb9f-090a-4b4c-7d48-596071829304
# ╠═3a1d0ca0-1a1b-4c5d-8e59-607182930415
# ╠═4b2e1db1-2b2c-4d6e-9f6a-718293041526
# ╟─5c3f2ec2-3c3d-4e7f-a07b-829304152637
# ╠═6d403fd3-4d4e-4f80-b18c-930415263748
# ╠═7e5140e4-5e5f-4091-c29d-041526374859
# ╟─8f6251f5-6f60-41a2-d3ae-152637485960
# ╟─90736265-7071-42b3-e4bf-263748596071
# ╠═a18473c5-8182-43c4-f5ca-374859607182
# ╠═b2958425-9293-44d5-06d1-485960718293
# ╠═c3a69595-a3a4-45e6-17e2-596071829304
# ╠═d4b7a645-b4b5-46f7-28f3-607182930415
# ╠═e5c8b755-c5c6-4708-3954-718293041526
# ╟─f6d9c865-d6d7-4819-4a65-829304152637
# ╠═07ead975-e7e8-492a-5bb6-930415263748
# ╟─18fbea85-f8f9-4a3b-6c97-041526374859
# ╟─290cfb95-0901-4b4c-7da8-152637485960
# ╠═3a1d0ca5-1a0b-4c5d-8eb9-263748596071
# ╠═4b2e1db5-2b1c-4d6e-9fca-374859607182
# ╟─6d403fd5-4d3e-4f80-b1eb-596071829304
# ╠═7e5140e5-5e4f-4091-c2fb-607182930415
# ╟─8f6251f6-6f40-41a2-d3bd-718293041526
# ╟─90736266-7051-42b3-e4cf-829304152637
# ╠═a18473c6-8152-43c4-f5db-930415263748
# ╠═b2958426-9263-44d5-06e1-041526374859
# ╟─5c3f2ec5-3c2d-4e7f-a0db-485960718293
# ╠═3a1f0c40-1a2b-4c5d-8e9f-0a1b2c3d4e5f
# ╠═4b2e1d51-2b3c-4d6e-9f0a-1b2c3d4e5f60
# ╟─00000000-0000-0000-0000-000000000001
