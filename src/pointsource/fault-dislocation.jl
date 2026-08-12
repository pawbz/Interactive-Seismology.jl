### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Faulting and the Surface"
#> tags = ["pointsource"]
#> layout = "layout.jlhtml"
#> description = "Drag a slipping rectangular fault in an elastic half-space and watch GPS and InSAR respond."

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

# ╔═╡ 27fe8f86-04ae-4c3f-af82-80c229db5edb
using PlutoUI

# ╔═╡ cc1d2a64-0a0e-4e46-b23e-89b4eab36db2
TableOfContents()

# ╔═╡ 3262a56e-a44b-42b5-98b3-68f191cceb39
md"""
# Faulting and the Surface

An earthquake is, mechanically, nothing more exotic than a crack that slips. Put a rectangular
patch of slip inside an elastic half-space, and everything you can measure at the surface — the
GPS offsets, the InSAR fringes — is the response of the surrounding elastic medium to that one
buried event. This notebook builds the whole class around a single concrete question, then hands
you a widget to answer it by hand.

The half-space is $z\le0$ (depth positive down). The fault living inside it is a rectangle,
parameterized by its center location, depth, size, and orientation:

```math
\boxed{x_0,\; y_0,\; d,\; L,\; W,\; \phi,\; \delta}
```

for horizontal location $(x_0,y_0)$, centroid depth $d$, length $L$, width $W$, strike $\phi$
(measured clockwise from north), and dip $\delta$ (from horizontal). Slip on the fault is
prescribed, not solved for — we decompose it into three components along the fault's own local
directions: along-strike $\hat{\mathbf s}$, along-dip $\hat{\mathbf d}$, and fault-normal
$\hat{\mathbf n}$ (opening/closing),

```math
\mathbf s = s_{\rm strike}\,\hat{\mathbf s} + s_{\rm dip}\,\hat{\mathbf d} + s_{\rm tensile}\,\hat{\mathbf n}.
```

Held uniform across the patch, for now. The question this entire notebook is built to answer:

```math
\boxed{\text{Given fault geometry + slip, what does the surface do?}}
```

The answer is a displacement field at the free surface,

```math
\mathbf u(x,y,0) = \begin{pmatrix} u_E \\ u_N \\ u_Z \end{pmatrix},
```

and that's exactly what a GPS receiver and a satellite radar measure: horizontal offset,
vertical offset, and (for InSAR) the displacement projected onto the satellite's line of sight.
Everything below is aimed at that one triple.
"""

# ╔═╡ 62575fc9-6613-421e-a3d6-9d287c913fb9
md"""
## The Governing Equation

Start from momentum conservation for the elastic medium surrounding the fault,

```math
\rho\,\ddot{\mathbf u} = \nabla\cdot\boldsymbol\sigma + \mathbf f.
```

A static (or slowly post-seismic) deformation field has no acceleration, $\ddot{\mathbf u}=0$,
and there's no body force away from the source, so this collapses to elastostatic equilibrium:

```math
\nabla\cdot\boldsymbol\sigma = 0.
```

Stress is related to strain through Hooke's law for an isotropic solid,

```math
\sigma_{ij} = \lambda\,\epsilon_{kk}\,\delta_{ij} + 2\mu\,\epsilon_{ij}, \qquad
\epsilon_{ij} = \tfrac12\left(u_{i,j} + u_{j,i}\right),
```

the same infinitesimal-strain tensor built up in the
[strain tensor notebook](../introduction/strain-tensor.html). Substituting strain into stress
into the equilibrium condition gives the **static Navier equation** — the one PDE this whole
notebook is quietly solving everywhere except right on the fault:

```math
\boxed{\mu\,\nabla^2\mathbf u + (\lambda+\mu)\,\nabla(\nabla\cdot\mathbf u) = 0}
```

By the time you reach the widget below, you've seen the full chain: displacement → strain →
stress → force balance.

### Why Ignore Inertia? The Quasi-Static Limit

Dropping `` \rho\,\ddot{\mathbf u} `` above deserves a real justification, not just the word
"static." Two length scales are in play, and everything hinges on comparing them.

The source sets a *spatial* scale — the fault length `` L `` (the same `` L `` you'll be
dragging in the widget below), which is roughly the distance over which the surface
displacement field varies. Separately, any process with characteristic timescale `` T `` has an
associated *wavelength*,

```math
\lambda = V_s\,T,
```

for shear waves traveling at speed `` V_s ``. `` L `` and `` \lambda `` are set by completely
different things — one geometric, one temporal — so nothing forces them to be comparable.

Go back to the full (not-yet-static) equation of motion for shear waves,
`` \rho\,\ddot u = \mu\,\nabla^2 u ``, and scale it: `` u\sim U ``, space `` \sim L ``,
time `` \sim T ``. Then `` \ddot u \sim U/T^2 `` while `` \nabla^2 u \sim U/L^2 ``, so

```math
\boxed{\frac{\text{inertia}}{\text{elastic}} \sim \frac{\rho\,U/T^2}{\mu\,U/L^2} = \left(\frac{L}{\lambda}\right)^2}
```

using `` V_s^2 = \mu/\rho ``. That ratio is the whole story:

- **`` L\sim\lambda ``** — the *dynamic* regime. Space and time are locked together into a
  genuinely propagating wave, inertia is not negligible, and this is ordinary
  wave-propagation seismology.
- **`` \lambda\gg L ``** — the *quasi-static* regime. Elastic waves cross the entire source
  region many times within one cycle of whatever is driving the deformation, so the medium
  has time to equilibrate everywhere before the loading changes appreciably. Inertia drops
  out and you're left with `` \nabla\cdot\boldsymbol\sigma=0 `` — exactly the equilibrium used
  above.

A concrete check: `` L=10 `` km, `` V_s=4 `` km/s. A slowly-loading process with period
`` T=1000 `` s gives `` \lambda = 4000 `` km, so `` L/\lambda\approx0.0025 `` — inertia is six
orders of magnitude smaller than the elastic term, comfortably quasi-static. Shrink `` T `` to
`` 2.5 `` s and `` \lambda=10 `` km `` \sim L ``: now you're in the wave-propagation regime, and
elastostatics no longer applies.

This is exactly the limit Okada's solution lives in — it's the `` \omega\to0 `` (equivalently
`` \lambda\gg L ``) limit of the full elastodynamic problem, appropriate for coseismic-to-postseismic
deformation measured well after the seismic waves have passed, not for the seismic waves
themselves. It's the assumption this whole notebook has been quietly making since the first
equation.
"""

# ╔═╡ c1aa8808-793f-4e87-98ab-27588b832dc4
md"""
## The Interesting Part: Boundary Conditions

The Navier equation above is satisfied everywhere in the half-space — it says nothing yet about
*why* there's any deformation at all. That comes entirely from two boundary conditions.

At the Earth's surface, traction vanishes because there's nothing above to push back:

```math
\boxed{\boldsymbol\sigma\cdot\hat{\mathbf n} = 0} \qquad (z=0).
```

Across the fault itself, the two faces are prescribed to have moved relative to each other by
exactly the slip vector, while traction stays continuous (the two faces still push on each other
with equal and opposite force — they just no longer line up):

```math
\boxed{\mathbf u^+ - \mathbf u^- = \mathbf s}.
```

This is worth sitting with: **we never apply a mysterious "earthquake force."** We prescribe a
displacement *discontinuity* across an internal surface, and let the elastic medium do the rest.

It also sharpens the connection to infinitesimal-strain theory noted above. The jump
$\mathbf s$ across the fault can be meters — not infinitesimal at all — while every point of the
surrounding medium, arbitrarily close to the fault, is still governed by the *linear*,
small-strain elasticity of the Navier equation. Large slip, small strain everywhere except
exactly on the discontinuity: that's what makes this a genuinely *elastic* dislocation problem
rather than a nonlinear fracture-mechanics one.
"""

# ╔═╡ 44f34984-05fb-447c-813c-9abb9aa47863
md"""
## Okada's Solution, as a Green's Function

Solving the Navier equation subject to those two boundary conditions, for an arbitrary
rectangular fault in a half-space, is a real derivation — and it buries the physics in algebra
that doesn't teach much on a first pass. We won't do it here. What matters is the *shape* of the
answer: because the governing equation is linear, the surface displacement is a linear function
of the slip, mediated by a **Green's function** that packages the entire boundary-value problem
for a given geometry and elastic moduli:

```math
\mathbf u(\mathbf x) = \mathbf G\bigl(\mathbf x;\, \text{fault geometry},\, \lambda,\, \mu\bigr)\,\mathbf s.
```

Okada (1985) worked out $\mathbf G$ in closed form, for exactly this problem — a uniform-slip
rectangular dislocation in an elastic half-space. It's a long but entirely explicit formula
(see the Appendix if you want to see every term), and it's the analytic engine behind essentially
every geodetic earthquake-source study since. From here, we stop deriving and start computing.
"""

# ╔═╡ 0f812b1a-45d5-481c-b4bc-eaf400f2da52
md"""
## What Depth Does to the Surface

Try the **Shallow source (d = 3 km)** and **Deep source (d = 30 km)** presets back to back, with
everything else held fixed. Switch the field checkbox between **Vertical GPS** and **InSAR** and
rotate the view to look along the surface — watch how the pattern spreads, not just its peak
amplitude.

```math
\boxed{\text{deeper source} \rightarrow \text{broader, weaker surface deformation}}
```

The peak drops because you're further from the source. The pattern *widens* because a buried
Green's function is smooth — a deeper dislocation looks, from the surface, like a low-pass
filtered version of a shallow one at the same location. This is exactly the resolution problem
that makes geodetic depth estimates hard: a deep, large-slip event and a shallow, small-slip
event can produce surface patterns that are hard to tell apart without independent depth
information (seismic arrival times, aftershock depths).
"""

# ╔═╡ f3862faf-3adb-4b2d-ab23-ecc0da655da3
md"""
## Fault Size and Seismic Moment

Now try **Large slow slip (moment)** — a fault an order of magnitude bigger in area than the
strike-slip default, at a comparable depth. The surface footprint grows to match the fault size
directly (deformation stays roughly localized over the source, unlike the depth effect above),
and the moment magnitude readout jumps with it.

Seismic moment ties the geometry you're dragging directly to earthquake size:

```math
\boxed{M_0 = \mu\,A\,D}
```

with $A=LW$ the rupture area and $D$ the slip magnitude. The **Readouts** panel above computes
this from the exact same `flt` geometry driving the surface plots — this is the same $M_0$ that
converts to the moment magnitude $M_w$ reported by every seismic catalog, computed here from pure
static geodesy rather than a seismogram. Geodesy and seismic source theory are describing the
same rupture; this is where they meet.
"""

# ╔═╡ a2f7c9e1-6b3d-4a8f-9c2e-1d5f8b0e3a7c
md"""
## Reading InSAR: Line-of-Sight Projection

Unlike a GPS receiver, a radar satellite doesn't hand you $(u_E,u_N,u_Z)$ directly. It measures a
*change in range* — the distance between satellite and ground — and that range change is only
sensitive to the part of the ground motion pointing back at the satellite, the **line of sight**
(LOS), a unit vector $\hat{\mathbf l}$:

```math
\boxed{d_{\rm LOS} = \hat{\mathbf l}\cdot\mathbf u = l_E u_E + l_N u_N + l_Z u_Z}
```

Two SAR images, before and after the earthquake, interfere to give a phase difference that wraps
every time the LOS range changes by half a radar wavelength (out-and-back doubles the one-way
path):

```math
\Delta\phi_{\rm disp} = \frac{4\pi}{\lambda_{\rm radar}}\,d_{\rm LOS}
\qquad\Longrightarrow\qquad
\text{one fringe} \;=\; \frac{\lambda_{\rm radar}}{2}.
```

That's exactly what the **InSAR fringes** checkbox draws: `uLOS` computed from the same
`(uE,uN,uZ)` Okada already gave you, wrapped into `[0,1)` cycles of $\lambda_{\rm radar}/2$. The
**radar λ** slider changes $\lambda_{\rm radar}$ directly — drag it from a short X-band wavelength
(~3 cm) toward a long L-band one (~24 cm) and watch the fringes spread apart for the *exact same*
deformation field: shorter wavelength packs more fringes into the same physical displacement, it
doesn't change $\mathbf u$ at all. The colorbar's "cm LOS" label tracks whatever the slider is
currently set to.

One real limitation the widget's own geometry makes concrete: a single interferogram gives one
number ($d_{\rm LOS}$) per pixel but three unknowns ($u_E,u_N,u_Z$). Polar-orbiting satellites
look sideways, nearly east-west, so $l_N\approx0$ in practice — this notebook's `okada_surface_field`
takes that simplification literally (`lookE`/`lookU` only, no north term at all), which is why
flipping **track** between ascending and descending flips the *sign* of the east sensitivity but
never recovers north-south motion. Reconstructing the full 3-D field for real needs multiple
viewing geometries (ascending + descending, or several tracks) combined, typically with GPS to
break the remaining ambiguity — a single fringe map alone can't do it.
"""

# ╔═╡ fcddbcd9-5018-424f-b0be-8301cbf4d3a8
md"""
## Appendix
"""

# ╔═╡ ea0ea1c3-1073-49ce-aa6a-12e3011d2538
md"""
### The Okada Green's Function

Okada's (1985) closed-form solution for the surface displacement due to a uniform-slip
rectangular dislocation in an elastic half-space. The auxiliary `_okada_*` functions below
implement Okada's own notation term for term (his Chinnery-notation four-corner evaluation, and
his $I_1,\dots,I_5$ auxiliary integrals); [`okada_surface_displacement`](@ref) is the public
entry point.
"""

# ╔═╡ cd707524-349c-412d-a7a8-e4dea3ca9b42
const OKADA_EPS = eps(1.0)

# ╔═╡ e0f54460-54ab-4c2d-a761-5e1fedad45cd
_okada_A(x, R) = (2R + x) / (R^3 * (R + x)^2)

# ╔═╡ 06919fe1-6209-4bd5-a1b1-a14d26b24e4e
begin
    function _okada_I5(xi, eta, q, dip, nu, R, dtilde)
        sd, cd = sincos(dip)
        X = sqrt(xi^2 + q^2)
        if cd > OKADA_EPS
            xi == 0 && return 0.0
            return (1 - 2nu) * 2 / cd * atan((eta * (X + q * cd) + X * (R + X) * sd) / (xi * (R + X) * cd))
        end
        return -(1 - 2nu) * xi * sd / (R + dtilde)
    end

    function _okada_I4(dtilde, eta, q, dip, nu, R)
        sd, cd = sincos(dip)
        cd > OKADA_EPS && return (1 - 2nu) / cd * (log(R + dtilde) - sd * log(R + eta))
        return -(1 - 2nu) * q / (R + dtilde)
    end

    function _okada_I3(eta, q, dip, nu, R)
        sd, cd = sincos(dip)
        ytilde = eta * cd + q * sd
        dtilde = eta * sd - q * cd
        if cd > OKADA_EPS
            return (1 - 2nu) * (ytilde / (cd * (R + dtilde)) - log(R + eta)) + sd / cd * _okada_I4(dtilde, eta, q, dip, nu, R)
        end
        return (1 - 2nu) / 2 * (eta / (R + dtilde) + ytilde * q / (R + dtilde)^2 - log(R + eta))
    end

    _okada_I2(eta, q, dip, nu, R) = (1 - 2nu) * (-log(R + eta)) - _okada_I3(eta, q, dip, nu, R)

    function _okada_I1(xi, eta, q, dip, nu, R)
        sd, cd = sincos(dip)
        dtilde = eta * sd - q * cd
        if cd > OKADA_EPS
            return (1 - 2nu) * (-xi / (cd * (R + dtilde))) - sd / cd * _okada_I5(xi, eta, q, dip, nu, R, dtilde)
        end
        return -(1 - 2nu) / 2 * xi * q / (R + dtilde)^2
    end
end

# ╔═╡ 0e36d361-eb42-4d4d-84b3-e943614b3b32
begin
    function _okada_ux_ss(xi, eta, q, dip, nu)
        sd = sin(dip)
        R = sqrt(xi^2 + eta^2 + q^2)
        u = xi * q / (R * (R + eta)) + _okada_I1(xi, eta, q, dip, nu, R) * sd
        q != 0 && (u += atan(xi * eta / (q * R)))
        return u
    end
    function _okada_uy_ss(xi, eta, q, dip, nu)
        sd, cd = sincos(dip)
        R = sqrt(xi^2 + eta^2 + q^2)
        return (eta * cd + q * sd) * q / (R * (R + eta)) + q * cd / (R + eta) + _okada_I2(eta, q, dip, nu, R) * sd
    end
    function _okada_uz_ss(xi, eta, q, dip, nu)
        sd, cd = sincos(dip)
        R = sqrt(xi^2 + eta^2 + q^2)
        dtilde = eta * sd - q * cd
        return dtilde * q / (R * (R + eta)) + q * sd / (R + eta) + _okada_I4(dtilde, eta, q, dip, nu, R) * sd
    end
    function _okada_ux_ds(xi, eta, q, dip, nu)
        sd, cd = sincos(dip)
        R = sqrt(xi^2 + eta^2 + q^2)
        return q / R - _okada_I3(eta, q, dip, nu, R) * sd * cd
    end
    function _okada_uy_ds(xi, eta, q, dip, nu)
        sd, cd = sincos(dip)
        R = sqrt(xi^2 + eta^2 + q^2)
        u = (eta * cd + q * sd) * q / (R * (R + xi)) - _okada_I1(xi, eta, q, dip, nu, R) * sd * cd
        q != 0 && (u += cd * atan(xi * eta / (q * R)))
        return u
    end
    function _okada_uz_ds(xi, eta, q, dip, nu)
        sd, cd = sincos(dip)
        R = sqrt(xi^2 + eta^2 + q^2)
        dtilde = eta * sd - q * cd
        u = dtilde * q / (R * (R + xi)) - _okada_I5(xi, eta, q, dip, nu, R, dtilde) * sd * cd
        q != 0 && (u += sd * atan(xi * eta / (q * R)))
        return u
    end
    function _okada_ux_tf(xi, eta, q, dip, nu)
        sd = sin(dip)
        R = sqrt(xi^2 + eta^2 + q^2)
        return q^2 / (R * (R + eta)) - _okada_I3(eta, q, dip, nu, R) * sd^2
    end
    function _okada_uy_tf(xi, eta, q, dip, nu)
        sd, cd = sincos(dip)
        R = sqrt(xi^2 + eta^2 + q^2)
        u = -(eta * sd - q * cd) * q / (R * (R + xi)) - sd * xi * q / (R * (R + eta)) - _okada_I1(xi, eta, q, dip, nu, R) * sd^2
        q != 0 && (u += sd * atan(xi * eta / (q * R)))
        return u
    end
    function _okada_uz_tf(xi, eta, q, dip, nu)
        sd, cd = sincos(dip)
        R = sqrt(xi^2 + eta^2 + q^2)
        dtilde = eta * sd - q * cd
        u = (eta * cd + q * sd) * q / (R * (R + xi)) + cd * xi * q / (R * (R + eta)) - _okada_I5(xi, eta, q, dip, nu, R, dtilde) * sd^2
        q != 0 && (u -= cd * atan(xi * eta / (q * R)))
        return u
    end
    _okada_chinnery(f, x, p, L, W, q, dip, nu) = f(x, p, q, dip, nu) - f(x, p - W, q, dip, nu) - f(x - L, p, q, dip, nu) + f(x - L, p - W, q, dip, nu)
end

# ╔═╡ f2f6bc7e-fc23-4e71-9c05-b7d5147ec4b2
"""
	okada_surface_displacement(E, N; depth, strike_deg, dip_deg, L, W, U1, U2, U3=0.0, nu=0.25)

Surface (`z=0`) displacement `(uE, uN, uZ)` at horizontal position `(E, N)` due to uniform slip
on a rectangular fault patch buried in a Poisson-solid elastic half-space (Okada, 1985).

`depth` is the depth (positive down) of the fault-plane centroid, in the same length unit as
`E`, `N`, `L`, `W`. `strike_deg` is measured clockwise from north, `dip_deg` from horizontal.
`U1` is strike-slip, `U2` is dip-slip, and `U3` is tensile opening, all in the same length unit
(e.g. metres) — independent of the geometric length unit used for `E`, `N`, `L`, `W`, `depth`
(e.g. kilometres).
"""
function okada_surface_displacement(E, N; depth, strike_deg, dip_deg, L, W, U1, U2, U3=0.0, nu=0.25)
    strike, dip = deg2rad(strike_deg), deg2rad(dip_deg)
    sd, cd = sincos(dip)
    d = depth + sd * W / 2
    ec = E + cos(strike) * cd * W / 2
    nc = N - sin(strike) * cd * W / 2
    x = cos(strike) * nc + sin(strike) * ec + L / 2
    y = sin(strike) * nc - cos(strike) * ec + cd * W
    p = y * cd + d * sd
    q = y * sd - d * cd

    ux = -U1 / 2π * _okada_chinnery(_okada_ux_ss, x, p, L, W, q, dip, nu) -
          U2 / 2π * _okada_chinnery(_okada_ux_ds, x, p, L, W, q, dip, nu) +
          U3 / 2π * _okada_chinnery(_okada_ux_tf, x, p, L, W, q, dip, nu)
    uy = -U1 / 2π * _okada_chinnery(_okada_uy_ss, x, p, L, W, q, dip, nu) -
          U2 / 2π * _okada_chinnery(_okada_uy_ds, x, p, L, W, q, dip, nu) +
          U3 / 2π * _okada_chinnery(_okada_uy_tf, x, p, L, W, q, dip, nu)
    uz = -U1 / 2π * _okada_chinnery(_okada_uz_ss, x, p, L, W, q, dip, nu) -
          U2 / 2π * _okada_chinnery(_okada_uz_ds, x, p, L, W, q, dip, nu) +
          U3 / 2π * _okada_chinnery(_okada_uz_tf, x, p, L, W, q, dip, nu)

    uE = sin(strike) * ux - cos(strike) * uy
    uN = cos(strike) * ux + sin(strike) * uy
    return (uE=uE, uN=uN, uZ=uz)
end

# ╔═╡ 8bdb3c71-79cf-4969-862a-dc2816406765
md"""
### Verifying the Implementation

A visible check against Okada's own published verification example (his Table, reproduced in
essentially every implementation since): a fault at $x=2,\,y=3,\,d=4$ (fault-local coordinates),
dip $70°$, $L=3$, $W=2$, unit slip in each of the three components in turn.
"""

# ╔═╡ 027c5b50-732c-4b03-8584-bb696635f596
let
    x, y, d, dip, L, W, nu = 2.0, 3.0, 4.0, deg2rad(70.0), 3.0, 2.0, 0.25
    p = y * cos(dip) + d * sin(dip)
    q = y * sin(dip) - d * cos(dip)
    ss = (_okada_chinnery(_okada_ux_ss, x, p, L, W, q, dip, nu), _okada_chinnery(_okada_uy_ss, x, p, L, W, q, dip, nu), _okada_chinnery(_okada_uz_ss, x, p, L, W, q, dip, nu))
    ds = (_okada_chinnery(_okada_ux_ds, x, p, L, W, q, dip, nu), _okada_chinnery(_okada_uy_ds, x, p, L, W, q, dip, nu), _okada_chinnery(_okada_uz_ds, x, p, L, W, q, dip, nu))
    tf = (_okada_chinnery(_okada_ux_tf, x, p, L, W, q, dip, nu), _okada_chinnery(_okada_uy_tf, x, p, L, W, q, dip, nu), _okada_chinnery(_okada_uz_tf, x, p, L, W, q, dip, nu))
    got_ss = round.((-ss[1] / 2π, -ss[2] / 2π, -ss[3] / 2π); sigdigits=4)
    got_ds = round.((-ds[1] / 2π, -ds[2] / 2π, -ds[3] / 2π); sigdigits=4)
    got_tf = round.((tf[1] / 2π, tf[2] / 2π, tf[3] / 2π); sigdigits=4)
    md"""
	| component | computed $(u_x,u_y,u_z)$ | Okada (1985) Table |
	|---|---|---|
	| strike-slip | $got_ss | (-8.689e-3, -4.298e-3, -2.747e-3) |
	| dip-slip    | $got_ds | (-4.682e-3, -3.527e-2, -3.564e-2) |
	| tensile     | $got_tf | (-2.660e-4, +1.056e-2, +3.214e-3) |
	"""
end

# ╔═╡ 504186a7-eb86-4e0c-81ac-7564aa966ba9
md"""
### Seismic Moment
"""

# ╔═╡ cdb34e89-4c64-4713-b289-7511be357d91
"""
	seismic_moment(L, W, U1, U2; mu=30e9)

Scalar seismic moment `M0 = μ A D` (N·m) and moment magnitude `Mw` (Hanks & Kanamori, 1979) for
a fault of length `L` and width `W` (km), slip components `U1` (strike-slip) and `U2` (dip-slip)
in metres, and shear modulus `mu` (Pa, default a representative crustal value).
"""
function seismic_moment(L, W, U1, U2; mu=30e9)
    Area = (L * 1000) * (W * 1000)
    D = hypot(U1, U2)
    M0 = mu * Area * D
    Mw = M0 > 0 ? (2 / 3) * (log10(M0) - 9.1) : NaN
    return (M0=M0, Mw=Mw)
end

# ╔═╡ eb3d526c-19c7-43bf-8a1c-eccd82bc15ce
md"""
### Sampling the Surface
"""

# ╔═╡ c3396aa7-69a9-4ec1-ae16-b4922ff636e7
"""
	okada_surface_field(flt; n=80, world_half=90.0, wavelength=0.056)

Evaluate [`okada_surface_displacement`](@ref) on an `n`×`n` grid spanning `±world_half`
kilometres in East and North around the origin, for the fault geometry and slip in `flt` (a
`Dict` as produced by `FaultDislocationInput`, with keys `"E0"`, `"N0"`, `"strikeDeg"`, `"L"`,
`"dipDeg"`, `"depth"`, `"W"`, `"sStrike"`, `"sDip"`, `"ascending"`, `"lambdaRadar"`).

Returns the coordinate vectors and the `uE`, `uN`, `uZ` displacement fields, the horizontal
magnitude `uHoriz = hypot(uE,uN)`, and the line-of-sight displacement `uLOS` with its wrapped
InSAR fringe `phase` (in `[0,1)`) for a simplified satellite geometry — 23° incidence from
vertical, ascending tracks looking east and descending tracks looking west — at `wavelength`
metres (default: Sentinel-1 C-band, 5.6 cm; pass `flt["lambdaRadar"]`, the widget's own radar-λ
slider, to keep the picture's fringe spacing in sync with what the colorbar reports).

`world_half` must match the fixed axis half-range the [`FaultDislocationInput`](@ref) widget
draws its map at, so the two line up spatially.
"""
function okada_surface_field(flt; n=56, world_half=90.0, wavelength=0.056)
    es = range(-world_half, world_half; length=n)
    ns = range(world_half, -world_half; length=n)
    uE = Matrix{Float64}(undef, n, n)
    uN = Matrix{Float64}(undef, n, n)
    uZ = Matrix{Float64}(undef, n, n)
    for (j, N) in enumerate(ns), (i, E) in enumerate(es)
        # okada_surface_displacement takes (E,N) relative to the fault CENTROID, not absolute
        # world coordinates -- forgetting this offset silently pins the computed field to the
        # origin regardless of where the fault (E0,N0) actually is.
        Erel = E - flt["E0"]
        Nrel = N - flt["N0"]
        Erel = (Erel == 0 && Nrel == 0) ? Erel + 1e-6 : Erel
        d = okada_surface_displacement(Erel, Nrel;
            depth=flt["depth"], strike_deg=flt["strikeDeg"], dip_deg=flt["dipDeg"],
            L=flt["L"], W=flt["W"], U1=flt["sStrike"], U2=flt["sDip"])
        uE[j, i], uN[j, i], uZ[j, i] = d.uE, d.uN, d.uZ
    end
    incidence = deg2rad(23.0)
    lookE = (flt["ascending"] ? 1 : -1) * sin(incidence)
    lookU = cos(incidence)
    uLOS = uE .* lookE .+ uZ .* lookU
    phase = mod.(uLOS ./ (wavelength / 2), 1.0)
    uHoriz = hypot.(uE, uN)
    return (es=collect(es), ns=collect(ns), uE=uE, uN=uN, uZ=uZ, uHoriz=uHoriz, uLOS=uLOS, phase=phase)
end

# ╔═╡ af3f2c0e-5bb8-490b-8332-c0ef928e9d0b
md"""
### Pushing the Field into the 3D View

`FieldPush` does no physics either — it takes the already-computed `field` (from
[`okada_surface_field`](@ref)) and hands its raw numbers to the *already-rendered*
[`FaultDislocationInput`](@ref) widget below, by dispatching a browser `CustomEvent` at the
widget's own `<div>`. The widget's own script (defined once, when it was first drawn) is
listening for that event and repaints its 3D view from the new numbers. Nothing here decides
*where* the deformation is — that arithmetic already happened in Julia; this cell only delivers
the answer to the picture.
"""

# ╔═╡ b7708494-743f-4578-b6e1-64e474f2723b
_flatten_rowmajor(M) = join(vec(permutedims(M)), ",")

# ╔═╡ cf9a4232-f9d4-4538-b243-d7c599ac896d
begin
    struct FieldPush
        field::Any
    end

    function Base.show(io::IO, ::MIME"text/html", p::FieldPush)
        f = p.field
        n = length(f.es)
        write(io, """
        <script>
        {
        const w = document.getElementById('fdwidget');
        if(w){
          w.dispatchEvent(new CustomEvent('fd-field-update', { detail: {
            n: $(n),
            es: [$(join(f.es, ","))],
            ns: [$(join(f.ns, ","))],
            uE: [$(_flatten_rowmajor(f.uE))],
            uN: [$(_flatten_rowmajor(f.uN))],
            uZ: [$(_flatten_rowmajor(f.uZ))],
            uHoriz: [$(_flatten_rowmajor(f.uHoriz))],
            phase: [$(_flatten_rowmajor(f.phase))],
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 038ee4a8-5a8c-43e9-bed4-d985366e5873
md"""
### The Interactive Widget

`FaultDislocationInput` does no physics either — it draws the 3D view (a flat, colorable ground
plane plus the tilted fault rectangle below it), handles dragging the fault's two endpoints,
dragging the on-fault arrow that sets the slip *direction* (rake), dragging to rotate the
camera, the field-select checkboxes, and the geometry/magnitude sliders, and reports the
resulting fault dictionary back to Julia. The rake and magnitude are pure display/UI state —
they're converted to the `sStrike`/`sDip` components Okada's formula actually takes
(``s_{\rm strike}=D\cos(\text{rake})``, ``s_{\rm dip}=D\sin(\text{rake})``) before being
reported, so nothing downstream changes. Everything downstream (`okada_surface_field`, `seismic_moment`,
`FieldPush`) is what turns that dictionary into the surface response and paints it back onto
this same widget.
"""

# ╔═╡ 5378e71e-753f-49fe-8679-dd57cd0cc329
begin
    struct FaultDislocationInput
        E0::Float64
        N0::Float64
        strikeDeg::Float64
        L::Float64
        dipDeg::Float64
        depth::Float64
        W::Float64
        sStrike::Float64
        sDip::Float64
        ascending::Bool
        lambdaRadar::Float64
    end
    FaultDislocationInput(; E0=0.0, N0=0.0, strikeDeg=10.0, L=40.0, dipDeg=90.0, depth=8.0, W=15.0, sStrike=3.0, sDip=0.0, ascending=true, lambdaRadar=0.056) =
        FaultDislocationInput(E0, N0, strikeDeg, L, dipDeg, depth, W, sStrike, sDip, ascending, lambdaRadar)

    Base.get(w::FaultDislocationInput) = Dict{String,Any}(
        "E0" => w.E0, "N0" => w.N0, "strikeDeg" => w.strikeDeg, "L" => w.L,
        "dipDeg" => w.dipDeg, "depth" => w.depth, "W" => w.W,
        "sStrike" => w.sStrike, "sDip" => w.sDip, "ascending" => w.ascending,
        "lambdaRadar" => w.lambdaRadar)

    function Base.show(io::IO, ::MIME"text/html", w::FaultDislocationInput)
        write(io, """
        <div id="fdwidget">
        <style>
        #fdwidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #fdwidget .fd-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #fdwidget .fd-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #fdwidget .fd-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #fdwidget .fd-workspace{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start}
        #fdwidget .fd-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #fdwidget .fd-caption{font-size:13px;color:#9ca3af;text-align:center;margin-top:4px}
        #fdwidget canvas{display:block;cursor:grab}
        #fdwidget .fd-controls{flex:0 0 260px;width:260px;display:flex;flex-direction:column;
          gap:8px;font:14px sans-serif}
        #fdwidget .fd-control-group{background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
        #fdwidget .fd-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #fdwidget .fd-control-row{display:grid;grid-template-columns:70px minmax(0,1fr) 52px;gap:6px;align-items:center;margin:5px 0}
        #fdwidget .fd-control-row label{font-size:13px;color:#9ca3af}
        #fdwidget .fd-control-row input[type=range]{width:100%;min-width:0}
        #fdwidget .fd-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #fdwidget .fd-check-row{display:flex;align-items:center;gap:6px;margin:5px 0;font-size:13px;color:#d1d5db}
        #fdwidget select{background:#0b0b0b;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:3px;width:100%}
        #fdwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:13px;cursor:pointer;margin-top:4px}
        </style>

        <div class="fd-title">
          <div class="fd-title-desc">Drag a rectangular slipping fault under an elastic half-space.</div>
          <div class="fd-title-hint">drag the yellow handles to move the fault &middot; drag the yellow arrow to set slip direction &middot; drag the orange dot to change depth &middot; drag the teal diamond to change dip/width &middot; drag empty space to rotate the view</div>
        </div>

        <div class="fd-workspace">
          <div>
            <div class="fd-panel"><canvas id="fd-3d"></canvas></div>
            <div class="fd-caption" id="fd-geom"></div>
          </div>
          <div class="fd-controls">
            <div class="fd-control-group">
              <div class="fd-control-title">Fault</div>
              <div class="fd-control-row"><label>preset</label>
                <select id="fd-preset">
                  <option value="ss">Vertical strike-slip</option>
                  <option value="thrust">Thrust earthquake</option>
                  <option value="shallow">Shallow source (d=3 km)</option>
                  <option value="deep">Deep source (d=30 km)</option>
                  <option value="big">Large slow slip (moment)</option>
                  <option value="custom">Custom</option>
                </select><span></span>
              </div>
              <div class="fd-control-row"><label>dip</label><input type="range" id="fd-dip" min="0" max="90" step="1" value="$(w.dipDeg)"><span class="fd-value" id="fd-dip-v"></span></div>
              <div class="fd-control-row"><label>depth</label><input type="range" id="fd-depth" min="1" max="30" step="0.5" value="$(w.depth)"><span class="fd-value" id="fd-depth-v"></span></div>
              <div class="fd-control-row"><label>width W</label><input type="range" id="fd-width" min="5" max="50" step="1" value="$(w.W)"><span class="fd-value" id="fd-width-v"></span></div>
            </div>
            <div class="fd-control-group">
              <div class="fd-control-title">Slip &amp; InSAR</div>
              <div class="fd-control-row"><label>magnitude</label><input type="range" id="fd-slipmag" min="0" max="10" step="0.1" value="0"><span class="fd-value" id="fd-slipmag-v"></span></div>
              <div class="fd-caption" id="fd-rake" style="text-align:left;margin:2px 0 6px">rake —</div>
              <div class="fd-control-row"><label>track</label>
                <select id="fd-track"><option value="asc">ascending</option><option value="desc">descending</option></select><span></span>
              </div>
              <div class="fd-control-row"><label>radar &lambda;</label><input type="range" id="fd-lambda" min="2" max="24" step="0.1" value="$(w.lambdaRadar*100)"><span class="fd-value" id="fd-lambda-v"></span></div>
              <button id="fd-reset">reset</button>
            </div>
            <div class="fd-control-group">
              <div class="fd-control-title">Surface heatmap</div>
              <div class="fd-check-row"><input type="checkbox" id="fd-chk-vert" checked><label for="fd-chk-vert">Vertical GPS</label></div>
              <div class="fd-check-row"><input type="checkbox" id="fd-chk-horiz"><label for="fd-chk-horiz">Horizontal GPS (magnitude)</label></div>
              <div class="fd-check-row"><input type="checkbox" id="fd-chk-insar"><label for="fd-chk-insar">InSAR fringes</label></div>
            </div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        const WORLD_HALF = 90; // km -- must match okada_surface_field's world_half in the Appendix
        const DEFAULT_N = 80;  // must match okada_surface_field's default n in the Appendix

        let E0=$(w.E0), N0=$(w.N0), strikeDeg=$(w.strikeDeg), L=$(w.L), dipDeg=$(w.dipDeg),
            depth=$(w.depth), W=$(w.W), sStrike=$(w.sStrike), sDip=$(w.sDip), ascending=$(w.ascending);
        // lambdaRadar is kept in METRES internally (matching okada_surface_field's `wavelength`
        // argument and Julia's `flt_safe["lambdaRadar"]`); only the slider's own value/label are
        // in centimetres, since that's the unit satellite radar bands are normally quoted in.
        let lambdaRadar=$(w.lambdaRadar);
        // rake/slipMag are the primary UI state for slip (set by dragging the arrow / the
        // magnitude slider); sStrike & sDip (what Okada / the Julia side actually take) are
        // kept in sync from them via syncSlipFromRakeMag().
        let rakeDeg = Math.atan2(sDip, sStrike) * 180 / Math.PI;
        let slipMag = Math.hypot(sStrike, sDip);
        function syncSlipFromRakeMag(){
          const r = rakeDeg * Math.PI / 180;
          sStrike = slipMag * Math.cos(r);
          sDip = slipMag * Math.sin(r);
        }
        function syncRakeMagFromSlip(){
          rakeDeg = Math.atan2(sDip, sStrike) * 180 / Math.PI;
          slipMag = Math.hypot(sStrike, sDip);
        }
        let selectedField = 'vert';
        let az = -35*Math.PI/180, el = 35*Math.PI/180;

        const PRESETS = {
          ss:      {strikeDeg:10, L:40,  dipDeg:90, depth:8,  W:15, sStrike:3, sDip:0},
          thrust:  {strikeDeg:0,  L:50,  dipDeg:30, depth:10, W:25, sStrike:0, sDip:4},
          shallow: {strikeDeg:10, L:40,  dipDeg:90, depth:3,  W:15, sStrike:3, sDip:0},
          deep:    {strikeDeg:10, L:40,  dipDeg:90, depth:30, W:15, sStrike:3, sDip:0},
          big:     {strikeDeg:0,  L:100, dipDeg:20, depth:14, W:50, sStrike:0, sDip:3},
        };

        function makeZeroField(n){
          const es=[], ns=[];
          for(let i=0;i<n;i++) es.push(-WORLD_HALF + i*(2*WORLD_HALF)/(n-1));
          for(let j=0;j<n;j++) ns.push(WORLD_HALF - j*(2*WORLD_HALF)/(n-1));
          const z = new Array(n*n).fill(0);
          return {n, es, ns, uE:z, uN:z, uZ:z, uHoriz:z, phase:z};
        }
        let fieldData = makeZeroField(DEFAULT_N);

        function handles(){
          const sr = strikeDeg*Math.PI/180;
          const dxA = Math.sin(sr)*L/2, dyA = Math.cos(sr)*L/2;
          return { A:{x:E0+dxA,y:N0+dyA}, B:{x:E0-dxA,y:N0-dyA} };
        }

        // fault-local basis at the centroid: sHat along strike, dHat up-dip (matches Okada's
        // own U1>0/U2>0 sign convention -- U2>0 is reverse motion, i.e. hanging wall moves
        // up-dip), both unit vectors, both lying in the fault plane.
        function faultBasis(){
          const sr = strikeDeg*Math.PI/180, dr = dipDeg*Math.PI/180;
          return {
            center: {x:E0, y:N0, z:-depth},
            sHat: {x:Math.sin(sr), y:Math.cos(sr), z:0},
            dHat: {x:Math.cos(sr)*Math.cos(dr), y:-Math.sin(sr)*Math.cos(dr), z:Math.sin(dr)},
          };
        }
        // tip position of the on-fault slip-direction arrow (rake measured from sHat toward
        // dHat, exactly as U1=slipMag*cos(rake), U2=slipMag*sin(rake) in Okada's formula) --
        // its length just scales with slipMag for visual feedback, capped relative to the
        // fault's own size so it never dwarfs or disappears against it.
        function slipArrowGeom(){
          const {center, sHat, dHat} = faultBasis();
          const rr = rakeDeg*Math.PI/180;
          const dir = {
            x: Math.cos(rr)*sHat.x + Math.sin(rr)*dHat.x,
            y: Math.cos(rr)*sHat.y + Math.sin(rr)*dHat.y,
            z: Math.cos(rr)*sHat.z + Math.sin(rr)*dHat.z,
          };
          const frac = Math.min(1, 0.4 + 0.6*(slipMag/10));
          const lenWorld = frac * Math.min(L,W) * 0.8;
          const tip = {x:center.x+dir.x*lenWorld, y:center.y+dir.y*lenWorld, z:center.z+dir.z*lenWorld};
          return {center, tip, sHat, dHat};
        }
        // surface point directly above the fault centroid, joined to the centroid by a dashed
        // plumb line -- deliberately NOT placed at the centroid itself (faultBasis().center),
        // which sits under the slip arrow's tail and can be within the hit-test radius of the
        // arrow's tip whenever rake is near 0 (the default "ss" preset draws the arrow exactly
        // along the strike line, i.e. along the A-B handles, so its tip lands close to the
        // centroid). Anchoring the depth handle at the surface (z=0) instead keeps it well
        // clear of every in-plane handle regardless of rake, strike, or fault size.
        function depthHandleGeom(){
          return { surfacePos: {x:E0, y:N0, z:0}, faultPos: {x:E0, y:N0, z:-depth} };
        }
        // midpoint of the fault's down-dip (far/deep) edge -- dragging it within the vertical
        // cross-section plane (spanned by awayDir, horizontal and perpendicular to strike, and
        // downDir, straight down) sets dip and width together: dip = atan2(down, away),
        // W = 2*hypot(away, down). At the current dip/W this point sits exactly at
        // center - dHat*(W/2), which is the same corner geometry the fault slab itself uses.
        function dipWidthHandleGeom(){
          const {center, dHat} = faultBasis();
          const half = W/2;
          const pos = {x:center.x-dHat.x*half, y:center.y-dHat.y*half, z:center.z-dHat.z*half};
          const sr = strikeDeg*Math.PI/180;
          const awayDir = {x:-Math.cos(sr), y:Math.sin(sr), z:0};
          const downDir = {x:0, y:0, z:-1};
          return {center, pos, awayDir, downDir};
        }

        const cv = par.querySelector('#fd-3d');
        const CONTROLS_W = 260 + 16; // .fd-controls width + workspace gap
        const availW = Math.min(window.innerWidth*0.8, (par.clientWidth || window.innerWidth*0.8) - CONTROLS_W);
        const SEC = Math.max(320, Math.min(availW, 760));
        // size the canvas's backing store in real device pixels (not just CSS pixels) so the
        // plot stays crisp on HiDPI/Retina screens -- all drawing code below still works in
        // SEC-sized logical units, because draw3D() resets the 2D context's transform to this
        // same DPR scale factor every time it runs, before any drawing happens.
        const DPR = window.devicePixelRatio || 1;
        cv.style.width = SEC + 'px'; cv.style.height = SEC + 'px';
        cv.width = Math.round(SEC * DPR); cv.height = Math.round(SEC * DPR);
        const scale = SEC/(2*WORLD_HALF)*0.85;
        const cx = SEC/2, cy = SEC*0.48;

        function proj(e,n,z){
          // z is elevation (positive up, so a buried fault passes z<0) -- z*ce must ADD to
          // y2 so that more-negative z (deeper) pushes the point further down/away on screen,
          // i.e. visually under the z=0 ground plane, not floating above it.
          const ca=Math.cos(az), sa=Math.sin(az);
          const x1 = e*ca - n*sa, y1 = e*sa + n*ca;
          const ce=Math.cos(el), se=Math.sin(el);
          const y2 = y1*se + z*ce;
          return [cx + x1*scale, cy - y2*scale];
        }
        function unprojGround(px,py){
          const ce=Math.cos(el), se=Math.sin(el);
          const x1 = (px-cx)/scale, y1 = (cy-py)/(se*scale);
          const ca=Math.cos(az), sa=Math.sin(az);
          return [x1*ca + y1*sa, -x1*sa + y1*ca];
        }

        function colorFor(mode, v, vmax){
          if(mode==='cyc') return 'hsl('+Math.round(v*360)+',85%,55%)';
          if(mode==='div'){
            const t=Math.max(-1,Math.min(1, v/vmax));
            if(t>=0){ const g=Math.round(255*(1-t)); return 'rgb(255,'+g+','+g+')'; }
            const g=Math.round(255*(1+t)); return 'rgb('+g+','+g+',255)';
          }
          const t=Math.max(0,Math.min(1, v/vmax));
          return 'rgb('+Math.round(255*t)+','+Math.round(160*t+30)+',30)';
        }

        function drawAxisTicks(ctx){
          const step = 30; // km
          ctx.font = '10px sans-serif'; ctx.fillStyle = '#9ca3af'; ctx.strokeStyle = '#4b5563'; ctx.lineWidth = 1;
          // East axis, along the N = -WORLD_HALF edge
          ctx.textAlign = 'center'; ctx.textBaseline = 'top';
          for (let g = -WORLD_HALF; g <= WORLD_HALF; g += step) {
            const [x0, y0] = proj(g, -WORLD_HALF, 0);
            ctx.beginPath(); ctx.arc(x0, y0, 1.5, 0, 2 * Math.PI); ctx.fill();
            ctx.fillText(g.toFixed(0), x0, y0 + 5);
          }
          const [ex, ey] = proj(0, -WORLD_HALF, 0);
          ctx.fillText('East (km)', ex, ey + 17);
          // North axis, along the E = -WORLD_HALF edge
          ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
          for (let g = -WORLD_HALF; g <= WORLD_HALF; g += step) {
            const [x0, y0] = proj(-WORLD_HALF, g, 0);
            ctx.beginPath(); ctx.arc(x0, y0, 1.5, 0, 2 * Math.PI); ctx.fill();
            ctx.fillText(g.toFixed(0), x0 - 5, y0);
          }
          ctx.save();
          const [nx, ny] = proj(-WORLD_HALF, 0, 0);
          ctx.translate(nx - 24, ny); ctx.rotate(-Math.PI / 2);
          ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
          ctx.fillText('North (km)', 0, 0);
          ctx.restore();
        }

        function drawColorbar(ctx, mode, vmax){
          // reserve enough right-margin for the widest label -- the LOS-per-fringe caption
          // scales with the radar-wavelength slider, so size for its longest possible text
          // (e.g. "(12.0 cm LOS)" at lambdaRadar=24cm) rather than the old fixed "2.8 cm" value
          const bw = 14, bh = 130, labelW = 112, bx = SEC - bw - labelW, by = 16;
          const steps = 60;
          for (let k = 0; k < steps; k++) {
            const t0 = k / steps, t1 = (k + 1) / steps;
            const v0 = mode === 'cyc' ? t0 : (mode === 'div' ? vmax * (1 - 2 * t0) : vmax * (1 - t0));
            const y0 = by + bh * t0, y1 = by + bh * t1;
            ctx.fillStyle = colorFor(mode, v0, mode === 'cyc' ? 1 : vmax);
            ctx.fillRect(bx, y0, bw, y1 - y0 + 1);
          }
          ctx.strokeStyle = '#e5e7eb'; ctx.lineWidth = 1; ctx.strokeRect(bx, by, bw, bh);
          ctx.fillStyle = '#e5e7eb'; ctx.font = '11px sans-serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
          if (mode === 'cyc') {
            // one fringe is one half-wavelength of line-of-sight range change -- this label
            // tracks the widget's own radar-wavelength slider, not a fixed Sentinel-1 value
            ctx.fillText('1 fringe', bx + bw + 4, by + 4);
            ctx.fillText('(' + (lambdaRadar*50).toFixed(1) + ' cm LOS)', bx + bw + 4, by + 18);
            ctx.fillText('0', bx + bw + 4, by + bh);
            ctx.textAlign = 'right'; ctx.fillText('InSAR phase', bx + bw, by - 8);
          } else {
            ctx.fillText('+' + vmax.toFixed(2) + ' m', bx + bw + 4, by + 4);
            ctx.fillText((mode === 'div' ? '-' : '') + (mode === 'div' ? vmax.toFixed(2) : '0') + ' m', bx + bw + 4, by + bh);
            if (mode === 'div') ctx.fillText('0 m', bx + bw + 4, by + bh / 2);
            ctx.textAlign = 'right';
            ctx.fillText((mode === 'div' ? 'vertical' : 'horizontal') + ' (m)', bx + bw, by - 8);
          }
        }

        function draw3D(){
          const ctx = cv.getContext('2d');
          // absolute (not cumulative) reset every frame -- safe to call every redraw, unlike
          // ctx.scale(), which would compound and blow up the transform over repeated calls
          ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
          ctx.clearRect(0,0,SEC,SEC); ctx.fillStyle='#000'; ctx.fillRect(0,0,SEC,SEC);
          const n = fieldData.n, es=fieldData.es, ns=fieldData.ns;
          let arr, mode;
          if(selectedField==='vert'){ arr=fieldData.uZ; mode='div'; }
          else if(selectedField==='horiz'){ arr=fieldData.uHoriz; mode='seq'; }
          else { arr=fieldData.phase; mode='cyc'; }
          let vmax=1e-9;
          if(mode!=='cyc') for(let k=0;k<arr.length;k++) vmax=Math.max(vmax, mode==='div'?Math.abs(arr[k]):arr[k]);

          // each (es[i],ns[j]) is a sample point, not a quad corner -- draw a cell CENTERED on
          // it (half a grid step in each direction), otherwise the fill for a value sampled at
          // one corner ends up rendered a half-cell away from where it was actually computed,
          // which is invisible for a smooth field but visibly shifts a sharp discontinuity
          // (like the jump straight across the fault) off of the fault line drawn below.
          const stepE = n>1 ? (es[1]-es[0]) : (2*WORLD_HALF);
          const stepN = n>1 ? (ns[0]-ns[1]) : (2*WORLD_HALF);
          for(let j=0;j<n;j++) for(let i=0;i<n;i++){
            const idx=j*n+i;
            const eLo=es[i]-stepE/2, eHi=es[i]+stepE/2;
            const nHi=ns[j]+stepN/2, nLo=ns[j]-stepN/2;
            const p1=proj(eLo,nHi,0), p2=proj(eHi,nHi,0), p3=proj(eHi,nLo,0), p4=proj(eLo,nLo,0);
            ctx.beginPath(); ctx.moveTo(p1[0],p1[1]); ctx.lineTo(p2[0],p2[1]); ctx.lineTo(p3[0],p3[1]); ctx.lineTo(p4[0],p4[1]); ctx.closePath();
            ctx.fillStyle = colorFor(mode, arr[idx], vmax); ctx.fill();
          }

          ctx.strokeStyle='#374151'; ctx.lineWidth=1;
          const rim=[[-WORLD_HALF,-WORLD_HALF],[WORLD_HALF,-WORLD_HALF],[WORLD_HALF,WORLD_HALF],[-WORLD_HALF,WORLD_HALF]];
          ctx.beginPath();
          rim.forEach((c,k)=>{ const [px,py]=proj(c[0],c[1],0); k===0?ctx.moveTo(px,py):ctx.lineTo(px,py); });
          ctx.closePath(); ctx.stroke();
          drawAxisTicks(ctx);

          const sr = strikeDeg*Math.PI/180, dr = dipDeg*Math.PI/180;
          const {A,B} = handles();
          const perp = { x: Math.cos(sr)*Math.cos(dr), y: -Math.sin(sr)*Math.cos(dr) };
          const half = W/2;
          const zNear = -depth + half*Math.sin(dr), zFar = -depth - half*Math.sin(dr);
          const c1 = proj(A.x+perp.x*half, A.y+perp.y*half, zNear);
          const c2 = proj(B.x+perp.x*half, B.y+perp.y*half, zNear);
          const c3 = proj(B.x-perp.x*half, B.y-perp.y*half, zFar);
          const c4 = proj(A.x-perp.x*half, A.y-perp.y*half, zFar);
          ctx.beginPath(); ctx.moveTo(c1[0],c1[1]); ctx.lineTo(c2[0],c2[1]); ctx.lineTo(c3[0],c3[1]); ctx.lineTo(c4[0],c4[1]); ctx.closePath();
          // mostly-opaque mid-slate so the fault plane still reads as a solid slab against any
          // heatmap color underneath (blue/red/white or the InSAR rainbow), without going as
          // dark as pure charcoal
          ctx.fillStyle='rgba(90,100,115,0.85)'; ctx.fill();
          ctx.strokeStyle='#cbd5e1'; ctx.lineWidth=1.5; ctx.stroke();

          const [ax,ay]=proj(A.x,A.y,0), [bx,by]=proj(B.x,B.y,0);
          ctx.strokeStyle='#facc15'; ctx.lineWidth=2; ctx.setLineDash([4,3]);
          ctx.beginPath(); ctx.moveTo(ax,ay); ctx.lineTo(bx,by); ctx.stroke(); ctx.setLineDash([]);
          [[ax,ay],[bx,by]].forEach(([px,py])=>{
            ctx.beginPath(); ctx.arc(px,py,7,0,2*Math.PI);
            ctx.fillStyle='#facc15'; ctx.fill(); ctx.strokeStyle='#000'; ctx.lineWidth=1.5; ctx.stroke();
          });

          // slip-direction arrow: drag its tip to set rake, magnitude comes from the slider.
          // Drawn yellow (matching the position handles) with a thin dark halo underneath so
          // it stays legible against every heatmap color and against the dark fault slab.
          {
            const {center, tip} = slipArrowGeom();
            const [cxp,cyp] = proj(center.x, center.y, center.z);
            const [txp,typ] = proj(tip.x, tip.y, tip.z);
            const ang = Math.atan2(typ-cyp, txp-cxp);
            const headLen = 12;
            const head = [
              [txp, typ],
              [txp-headLen*Math.cos(ang-0.4), typ-headLen*Math.sin(ang-0.4)],
              [txp-headLen*Math.cos(ang+0.4), typ-headLen*Math.sin(ang+0.4)],
            ];
            ctx.lineCap = 'round';
            // halo
            ctx.strokeStyle='#111827'; ctx.lineWidth=6;
            ctx.beginPath(); ctx.moveTo(cxp,cyp); ctx.lineTo(txp,typ); ctx.stroke();
            ctx.beginPath(); head.forEach((p,k)=>k===0?ctx.moveTo(p[0],p[1]):ctx.lineTo(p[0],p[1])); ctx.closePath();
            ctx.fillStyle='#111827'; ctx.fill();
            // core
            ctx.strokeStyle='#facc15'; ctx.lineWidth=3;
            ctx.beginPath(); ctx.moveTo(cxp,cyp); ctx.lineTo(txp,typ); ctx.stroke();
            ctx.beginPath(); head.forEach((p,k)=>k===0?ctx.moveTo(p[0],p[1]):ctx.lineTo(p[0],p[1])); ctx.closePath();
            ctx.fillStyle='#facc15'; ctx.fill();
            ctx.lineCap = 'butt';
          }

          // depth handle: a plumb line from the surface down to the fault centroid; drag the
          // surface dot straight up/down to move the whole fault shallower/deeper.
          {
            const {surfacePos, faultPos} = depthHandleGeom();
            const [spx,spy] = proj(surfacePos.x, surfacePos.y, surfacePos.z);
            const [fpx,fpy] = proj(faultPos.x, faultPos.y, faultPos.z);
            ctx.strokeStyle='#fb923c'; ctx.lineWidth=1.5; ctx.setLineDash([3,3]);
            ctx.beginPath(); ctx.moveTo(spx,spy); ctx.lineTo(fpx,fpy); ctx.stroke(); ctx.setLineDash([]);
            ctx.beginPath(); ctx.arc(spx,spy,7,0,2*Math.PI);
            ctx.fillStyle='#fb923c'; ctx.fill(); ctx.strokeStyle='#000'; ctx.lineWidth=1.5; ctx.stroke();
          }
          // dip/width handle: drag the fault's down-dip edge to reshape both at once
          {
            const {pos} = dipWidthHandleGeom();
            const [px,py] = proj(pos.x, pos.y, pos.z);
            ctx.beginPath();
            ctx.moveTo(px,py-8); ctx.lineTo(px+8,py); ctx.lineTo(px,py+8); ctx.lineTo(px-8,py); ctx.closePath();
            ctx.fillStyle='#2dd4bf'; ctx.fill(); ctx.strokeStyle='#000'; ctx.lineWidth=1.5; ctx.stroke();
          }

          drawColorbar(ctx, mode, vmax);

          par.querySelector('#fd-geom').textContent = 'L = '+L.toFixed(0)+' km, strike = '+strikeDeg.toFixed(0)+'°';
          par.querySelector('#fd-rake').textContent = 'rake ' + Math.round(((rakeDeg%360)+360)%360) + '°';
        }

        function syncControls(){
          par.querySelector('#fd-dip').value = dipDeg;
          par.querySelector('#fd-dip-v').textContent = dipDeg.toFixed(0)+'°';
          par.querySelector('#fd-depth').value = depth;
          par.querySelector('#fd-depth-v').textContent = depth.toFixed(1)+' km';
          par.querySelector('#fd-width').value = W;
          par.querySelector('#fd-width-v').textContent = W.toFixed(0)+' km';
          par.querySelector('#fd-slipmag').value = slipMag;
          par.querySelector('#fd-slipmag-v').textContent = slipMag.toFixed(1)+' m';
          par.querySelector('#fd-track').value = ascending ? 'asc' : 'desc';
          par.querySelector('#fd-lambda').value = lambdaRadar*100;
          par.querySelector('#fd-lambda-v').textContent = (lambdaRadar*100).toFixed(1)+' cm';
        }

        let commitInFlight = false;
        function commit(){
          commitInFlight = true;
          par.value = { E0, N0, strikeDeg, L, dipDeg, depth, W, sStrike, sDip, ascending, lambdaRadar };
          par.dispatchEvent(new CustomEvent('input'));
        }

        // live updates while dragging/sliding, not just on release: paced by round-trip
        // COMPLETION (commitInFlight), not a fixed timer -- a fixed-interval throttle turned
        // out to queue requests faster than Julia could actually drain them under sustained
        // dragging (measured a ~90ms single round trip, but a steady stream of 120ms-spaced
        // commits built an 11-SECOND backlog of stale updates trickling in after the drag
        // ended). Gating on "is a request already outstanding" instead means at most one commit
        // is ever in flight, so the picture always catches up to wherever the cursor currently
        // is, as fast as the kernel can go, with no backlog possible by construction. The final
        // value is still always sent unconditionally on release/change (see below), so this
        // gate never causes the last tick to be dropped -- only intermediate ticks are skipped.
        function throttledCommit(){
          if (!commitInFlight) commit();
        }

        function applyPreset(name){
          const p = PRESETS[name]; if(!p) return;
          E0=0; N0=0; strikeDeg=p.strikeDeg; L=p.L; dipDeg=p.dipDeg; depth=p.depth; W=p.W; sStrike=p.sStrike; sDip=p.sDip;
          syncRakeMagFromSlip();
          par.querySelector('#fd-preset').value = name;
          syncControls(); draw3D(); commit();
        }
        function setCustom(){ par.querySelector('#fd-preset').value = 'custom'; }
        function setField(id){
          selectedField = id==='fd-chk-vert' ? 'vert' : (id==='fd-chk-horiz' ? 'horiz' : 'insar');
          ['fd-chk-vert','fd-chk-horiz','fd-chk-insar'].forEach(cid=>{ par.querySelector('#'+cid).checked = (cid===id); });
          draw3D();
        }

        // capture-phase: block descendant controls' native input/change from reaching Pluto's
        // own bond listener directly (it would fire with stale par.value before our state
        // updates below), then re-dispatch our own clean event once state is current.
        par.addEventListener('input', e => {
          if(e.target === par) return;
          e.stopImmediatePropagation();
          const id = e.target.id, v = e.target.value;
          // radar wavelength is a sensor property, not a fault parameter -- change it without
          // flipping the fault-geometry preset selector to "Custom" (setCustom() below).
          if(id==='fd-lambda'){ lambdaRadar=+v/100; par.querySelector('#fd-lambda-v').textContent=(+v).toFixed(1)+' cm'; draw3D(); throttledCommit(); return; }
          if(id==='fd-dip'){ dipDeg=+v; par.querySelector('#fd-dip-v').textContent=dipDeg.toFixed(0)+'°'; }
          else if(id==='fd-depth'){ depth=+v; par.querySelector('#fd-depth-v').textContent=depth.toFixed(1)+' km'; }
          else if(id==='fd-width'){ W=+v; par.querySelector('#fd-width-v').textContent=W.toFixed(0)+' km'; }
          else if(id==='fd-slipmag'){ slipMag=+v; par.querySelector('#fd-slipmag-v').textContent=slipMag.toFixed(1)+' m'; syncSlipFromRakeMag(); }
          else return;
          setCustom(); draw3D(); throttledCommit();
        }, true);

        par.addEventListener('change', e => {
          if(e.target === par) return;
          e.stopImmediatePropagation();
          const id = e.target.id, v = e.target.value;
          if(id==='fd-preset'){ if(v!=='custom') applyPreset(v); return; }
          if(id==='fd-track'){ ascending = (v==='asc'); draw3D(); commit(); return; }
          if(id==='fd-lambda'){ draw3D(); commit(); return; }
          if(id==='fd-chk-vert' || id==='fd-chk-horiz' || id==='fd-chk-insar'){
            if(!e.target.checked){ e.target.checked = true; return; } // keep exactly one selected
            setField(id);
            return; // purely local recoloring -- no Julia round trip needed
          }
          if(['fd-dip','fd-depth','fd-width','fd-slipmag'].includes(id)){ draw3D(); commit(); return; }
        }, true);

        par.querySelector('#fd-reset').addEventListener('click', ()=> applyPreset('ss'));
        par.addEventListener('fd-field-update', e=>{ fieldData = e.detail; commitInFlight = false; draw3D(); });

        let dragMode = null, lastMouse = null;
        cv.addEventListener('mousedown', e=>{
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left, py = e.clientY-rect.top;
          const {A,B} = handles();
          const [ax,ay]=proj(A.x,A.y,0), [bx,by]=proj(B.x,B.y,0);
          const {tip} = slipArrowGeom();
          const [tpx,tpy] = proj(tip.x,tip.y,tip.z);
          const {surfacePos: depthPos} = depthHandleGeom();
          const [dpx,dpy] = proj(depthPos.x,depthPos.y,depthPos.z);
          const {pos: dwPos} = dipWidthHandleGeom();
          const [dwx,dwy] = proj(dwPos.x,dwPos.y,dwPos.z);
          // pick whichever handle is genuinely CLOSEST (not a fixed priority chain) -- with
          // some fault geometries two handles' screen positions land within the hit radius of
          // each other (e.g. the dip/width corner can sit near the B endpoint when W~L), and a
          // fixed if/else order would always resolve the tie the same wrong way.
          const candidates = [
            ['A', ax, ay], ['B', bx, by], ['slip', tpx, tpy],
            ['depth', dpx, dpy], ['dipwidth', dwx, dwy],
          ];
          let best = null, bestDist = 12;
          for(const [name, hx, hy] of candidates){
            const d = Math.hypot(px-hx, py-hy);
            if(d < bestDist){ bestDist = d; best = name; }
          }
          dragMode = best || 'rotate';
          lastMouse = {x:px,y:py};
        });
        window.addEventListener('mousemove', e=>{
          if(!dragMode) return;
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left, py = e.clientY-rect.top;
          if(dragMode==='rotate'){
            const dx = px-lastMouse.x, dy = py-lastMouse.y;
            az += dx*0.008;
            el = Math.max(10*Math.PI/180, Math.min(89*Math.PI/180, el - dy*0.008));
            lastMouse = {x:px,y:py};
            draw3D();
            return;
          }
          if(dragMode==='slip'){
            // invert the (linear, orthographic) projection restricted to the fault plane's own
            // basis: [screenSHat screenDHat] * [u,v]^T = mouse offset from the projected center
            const {center, sHat, dHat} = faultBasis();
            const [cxp,cyp] = proj(center.x,center.y,center.z);
            const [sxp,syp] = proj(center.x+sHat.x, center.y+sHat.y, center.z+sHat.z);
            const [dxp,dyp] = proj(center.x+dHat.x, center.y+dHat.y, center.z+dHat.z);
            const sx=sxp-cxp, sy=syp-cyp, dx2=dxp-cxp, dy2=dyp-cyp;
            const det = sx*dy2 - dx2*sy;
            if(Math.abs(det) > 1e-6){
              const mx = px-cxp, my = py-cyp;
              const u = (mx*dy2 - dx2*my)/det;
              const v = (sx*my - mx*sy)/det;
              rakeDeg = Math.atan2(v,u)*180/Math.PI;
              syncSlipFromRakeMag();
              setCustom(); draw3D(); throttledCommit();
            }
            return;
          }
          if(dragMode==='depth'){
            // pure vertical drag: screenY changes by cos(el)*scale per unit of depth (derived
            // from proj()'s z-term, which is a fixed screen-space slope independent of position),
            // so this is a plain 1D conversion, not a 2D basis inversion like the other handles.
            const dy = py - lastMouse.y;
            depth = Math.max(1, Math.min(30, depth + dy/(Math.cos(el)*scale)));
            lastMouse = {x:px,y:py};
            syncControls(); setCustom(); draw3D(); throttledCommit();
            return;
          }
          if(dragMode==='dipwidth'){
            // same style of basis inversion as the slip arrow, but in the vertical cross-section
            // plane (awayDir horizontal-perpendicular-to-strike, downDir straight down) instead
            // of the fault plane itself: recovers dip and W together from one drag.
            const {center, awayDir, downDir} = dipWidthHandleGeom();
            const [cxp,cyp] = proj(center.x,center.y,center.z);
            const [axp,ayp] = proj(center.x+awayDir.x, center.y+awayDir.y, center.z+awayDir.z);
            const [dxp,dyp] = proj(center.x+downDir.x, center.y+downDir.y, center.z+downDir.z);
            const ax2=axp-cxp, ay2=ayp-cyp, dx2=dxp-cxp, dy2=dyp-cyp;
            const det = ax2*dy2 - dx2*ay2;
            if(Math.abs(det) > 1e-6){
              const mx = px-cxp, my = py-cyp;
              const a = (mx*dy2 - dx2*my)/det;
              const d = (ax2*my - mx*ay2)/det;
              dipDeg = Math.max(0, Math.min(90, Math.atan2(d,a)*180/Math.PI));
              W = Math.max(5, Math.min(50, 2*Math.hypot(a,d)));
              syncControls(); setCustom(); draw3D(); throttledCommit();
            }
            return;
          }
          let [wx,wy] = unprojGround(px,py);
          wx = Math.max(-WORLD_HALF, Math.min(WORLD_HALF, wx));
          wy = Math.max(-WORLD_HALF, Math.min(WORLD_HALF, wy));
          const {A,B} = handles();
          const other = dragMode==='A' ? B : A;
          const moved = {x:wx,y:wy};
          E0 = (moved.x+other.x)/2; N0 = (moved.y+other.y)/2;
          L = Math.max(2, Math.hypot(moved.x-other.x, moved.y-other.y));
          const dx = dragMode==='A' ? moved.x-other.x : other.x-moved.x;
          const dy = dragMode==='A' ? moved.y-other.y : other.y-moved.y;
          let s = Math.atan2(dx,dy)*180/Math.PI; if(s<0) s+=360;
          strikeDeg = s;
          setCustom(); draw3D(); throttledCommit();
        });
        window.addEventListener('mouseup', ()=>{
          if(dragMode==='A'||dragMode==='B'||dragMode==='slip'||dragMode==='depth'||dragMode==='dipwidth') commit();
          dragMode=null;
        });

        syncControls(); draw3D();
        }
        </script>
        """)
    end

    const _fd_ready = true
end

# ╔═╡ 79bf0c9d-4e9e-4f28-9905-76af5f60a52c
begin
    _fd_ready
    WideCell(@bind flt FaultDislocationInput(); max_width=1500)
end

# ╔═╡ cb9f6341-becd-4ed2-bda1-dd309733762c
# The bond starts as `nothing` until the widget's first real interaction in a live browser
# reports back — fall back to the same defaults the widget itself opens with.
flt_safe = flt isa AbstractDict ? flt : Dict{String,Any}(
    "E0" => 0.0, "N0" => 0.0, "strikeDeg" => 10.0, "L" => 40.0, "dipDeg" => 90.0,
    "depth" => 8.0, "W" => 15.0, "sStrike" => 3.0, "sDip" => 0.0, "ascending" => true,
    "lambdaRadar" => 0.056)

# ╔═╡ b3975a01-14a4-4de9-8781-6fb8c6439fa1
field = okada_surface_field(flt_safe; wavelength=flt_safe["lambdaRadar"])

# ╔═╡ eff5c5f6-a648-4ff0-8027-a6d5b64c60af
FieldPush(field)

# ╔═╡ 258d0573-8a98-40f2-ba4f-46573e43a865
moment = seismic_moment(flt_safe["L"], flt_safe["W"], flt_safe["sStrike"], flt_safe["sDip"])

# ╔═╡ 469ddcaa-fcd0-47ec-87d1-3bca5d2fd764
md"""
strike **$(round(flt_safe["strikeDeg"]; digits=0))°**, dip **$(round(flt_safe["dipDeg"]; digits=0))°**,
depth **$(flt_safe["depth"])** km, $L\times W$ **$(round(flt_safe["L"]; digits=0)) × $(flt_safe["W"])** km ·
slip: strike **$(flt_safe["sStrike"])** m, dip **$(flt_safe["sDip"])** m

``M_0=\mu A D`` = **$(round(moment.M0, sigdigits=3))** N·m &nbsp;→&nbsp; ``M_w`` = **$(round(moment.Mw; digits=2))**

max uplift **$(round(maximum(field.uZ); digits=2))** m · max subsidence **$(round(minimum(field.uZ); digits=2))** m ·
max horizontal **$(round(maximum(hypot.(field.uE, field.uN)); digits=2))** m
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "40c9f1cac973d64f8ca3ef3a09f769ff947e80f3"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

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

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

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
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

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
git-tree-sha1 = "3b0738bd7c5645641845da25cbd99800b8718689"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.2"

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
# ╟─cc1d2a64-0a0e-4e46-b23e-89b4eab36db2
# ╟─3262a56e-a44b-42b5-98b3-68f191cceb39
# ╟─62575fc9-6613-421e-a3d6-9d287c913fb9
# ╟─c1aa8808-793f-4e87-98ab-27588b832dc4
# ╟─44f34984-05fb-447c-813c-9abb9aa47863
# ╟─79bf0c9d-4e9e-4f28-9905-76af5f60a52c
# ╠═cb9f6341-becd-4ed2-bda1-dd309733762c
# ╟─b3975a01-14a4-4de9-8781-6fb8c6439fa1
# ╟─258d0573-8a98-40f2-ba4f-46573e43a865
# ╟─eff5c5f6-a648-4ff0-8027-a6d5b64c60af
# ╟─469ddcaa-fcd0-47ec-87d1-3bca5d2fd764
# ╟─0f812b1a-45d5-481c-b4bc-eaf400f2da52
# ╟─f3862faf-3adb-4b2d-ab23-ecc0da655da3
# ╠═a2f7c9e1-6b3d-4a8f-9c2e-1d5f8b0e3a7c
# ╟─fcddbcd9-5018-424f-b0be-8301cbf4d3a8
# ╟─ea0ea1c3-1073-49ce-aa6a-12e3011d2538
# ╠═cd707524-349c-412d-a7a8-e4dea3ca9b42
# ╠═e0f54460-54ab-4c2d-a761-5e1fedad45cd
# ╠═06919fe1-6209-4bd5-a1b1-a14d26b24e4e
# ╠═0e36d361-eb42-4d4d-84b3-e943614b3b32
# ╠═f2f6bc7e-fc23-4e71-9c05-b7d5147ec4b2
# ╟─8bdb3c71-79cf-4969-862a-dc2816406765
# ╠═027c5b50-732c-4b03-8584-bb696635f596
# ╟─504186a7-eb86-4e0c-81ac-7564aa966ba9
# ╠═cdb34e89-4c64-4713-b289-7511be357d91
# ╟─eb3d526c-19c7-43bf-8a1c-eccd82bc15ce
# ╠═c3396aa7-69a9-4ec1-ae16-b4922ff636e7
# ╟─af3f2c0e-5bb8-490b-8332-c0ef928e9d0b
# ╠═b7708494-743f-4578-b6e1-64e474f2723b
# ╠═cf9a4232-f9d4-4538-b243-d7c599ac896d
# ╟─038ee4a8-5a8c-43e9-bed4-d985366e5873
# ╠═5378e71e-753f-49fe-8679-dd57cd0cc329
# ╟─27fe8f86-04ae-4c3f-af82-80c229db5edb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
