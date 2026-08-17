### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Waves On A String"
#> tags = ["introduction"]
#> layout = "layout.jlhtml"
#> description = "A staggered-grid finite-difference tutorial: reflection/transmission at an impedance contrast, and free-vs-fixed boundary polarity"

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

# ╔═╡ dbfb7f10-8fe2-471b-8626-7a8894626caf
using PlutoUI

# ╔═╡ 46dca275-df8b-4060-b81d-841428142791
TableOfContents()

# ╔═╡ f3c52928-1450-4596-87b7-b6665a1b8f6a
md"""
# Waves On A String

Drop a velocity bump on a 1-D elastic string and watch it split into two
pulses, run into a change in the string's density, and bounce off the ends.
Every 2-D and 3-D finite-difference seismic code you'll meet later in this
course does the same three things this notebook does, just in more
dimensions: it puts different physical quantities on *different* grid points
(a "staggered grid"), it only stays numerically stable below a specific ratio
of time step to space step (the "Courant condition"), and it has to decide
what happens at a boundary -- free, like Earth's actual surface, or
artificially rigid, like a wall.

A tempting first guess is to put velocity and stress on the *same* grid
point and update them together -- it's simpler to code. It also smears the
two fields together: a spatial derivative evaluated *at* a grid point needs
neighbors on both sides, so collocating both fields throws away exactly the
symmetry that makes a centered difference accurate. Offsetting them by half a
grid cell -- a *staggered grid* -- turns every derivative in the scheme
below into a clean, centered difference for free, which is why it tracks a
sharp pulse without extra numerical smearing.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ f69bc5a5-d354-41b6-8d12-61a04ef1c43d
md"""
## Governing Equations

The transverse velocity ``v_y`` and shear stress ``\sigma`` on the string obey
a coupled pair of first-order equations -- the same velocity-stress form used
for elastic wave propagation in 2-D and 3-D finite-difference seismic codes:

```math
\rho\,\partial_t v_y = \partial_x \sigma \qquad (1)
```
```math
\partial_t \sigma = \mu\,\partial_x v_y \qquad (2)
```

Combining them eliminates ``\sigma`` and gives the familiar scalar wave
equation ``\partial_t^2 v_y = c^2\,\partial_x^2 v_y`` with speed
``c=\sqrt{\mu/\rho}`` -- but the widget below solves the *coupled first-order*
pair directly, not the second-order equation, because that's what lets
``\rho`` and ``\mu`` vary independently in space (a real seismic medium)
without ever needing a second spatial derivative of a discontinuous
quantity.
"""

# ╔═╡ 7676054a-2e9a-43f4-ac46-af162541e269
md"""
## The Staggered Grid

The widget's top panel is a schematic of the grid the solver actually uses:
``v_y`` lives on integer grid nodes ``x_i``, and ``\sigma`` lives on the
*half*-grid nodes ``x_{i+1/2}`` sitting exactly between them. One leapfrog
time step does two things, in order:

1. **Update `` \sigma``** at every half-node from its two neighboring
   ``v_y`` values (the orange arrows in the schematic) -- a centered
   difference of ``v_y`` with no interpolation needed, since the two ``v_y``
   nodes it needs are exactly ``\pm\,dx/2`` away.
2. **Update `` v_y``** at every interior node from its two neighboring
   ``\sigma`` values (the blue arrows) -- likewise a clean centered
   difference.

This is why the scheme is called "staggered": it isn't one grid with two
fields interleaved awkwardly, it's two half-cell-shifted grids, each of which
sees the *other* field sitting perfectly centered around it.
"""

# ╔═╡ d0be7679-5df5-4ad0-93d1-8df6f94ba299
md"""
## Numerical Stability: The Courant Condition

The time step ``dt`` and grid spacing ``dx`` can't be chosen independently:
for this leapfrog scheme to stay numerically stable, information can't be
allowed to cross more than one grid cell per time step. That requirement is
the Courant number,

```math
C = \frac{v_{\max}\,dt}{dx},
```

which must stay below ``1`` (using the fastest wave speed anywhere in the
medium for ``v_{\max}``). Every explicit finite-difference wave solver you'll
meet later in this course has some version of this same limit.

!!! warning "Try it"
	Push the **Courant** slider above `1` in the widget above and watch the
	string blow up within a few steps -- that's the scheme amplifying its own
	rounding error every iteration once the stability condition is violated,
	not a bug.
"""

# ╔═╡ 6feab4ba-89ed-4be0-87e7-c1cc1d78b6c2
md"""
## Reflection and Transmission at an Interface

At the density jump the widget lets you drag, the string behaves like two
different media joined at ``x=`` **Interface position**, each with its own
mechanical impedance ``Z=\rho c`` (``c=\sqrt{\mu/\rho}``, the local wave
speed). Matching velocity and stress across the boundary gives the (velocity)
reflection and transmission coefficients

```math
R = \frac{Z_1-Z_2}{Z_1+Z_2}, \qquad T = \frac{2Z_1}{Z_1+Z_2},
```

for a pulse arriving from medium 1 (the denser side) into medium 2. Both are
printed live on the widget's main panel. Edit the density-ratio number box
floating over medium 1 toward `1` and the interface disappears (``R\to0``);
push it far from `1` and almost everything reflects. This is the same
physics -- just in one dimension, at normal incidence -- as the
oblique-incidence reflection coefficients in the P-SV and SH free-surface
notebooks later in this course.
"""

# ╔═╡ 260264ff-24bb-425f-abb8-d3556a1bb19f
md"""
## Boundary Conditions: Free Surface vs. Fixed End

The two ends of the string need their own boundary condition, and the widget
lets you toggle between the two extremes:

- **Fixed (rigid)**: the end can't move, ``v_y=0`` there always -- physically,
  what an infinitely-hard wall looks like, one limit of the reflection
  formula above (the far side's impedance going to infinity).
- **Free (`` \sigma=0``)**: no wall at all -- the string's own end carries no
  stress, since there's no material beyond it to stress against. This is the
  boundary condition that actually applies at Earth's real surface.

The two behave oppositely on reflection: a pulse reflecting off a **fixed**
end comes back **inverted**, the same story as a string tied to a wall, while
a pulse reflecting off a **free** end comes back **the same sign** it arrived
with. Toggle the boundary control in the widget and watch the sign of the
returning pulse flip -- the Appendix verifies this numerically below, and
it's the reason a seismogram's first surface-reflected arrival doesn't carry
a sign flip the way a rigid-wall reflection would.
"""

# ╔═╡ 8a11db06-56c3-42cc-b786-d22b5387bbeb
md"## Appendix"

# ╔═╡ 90e40955-3aaa-4efb-b48b-21ae94a93bd5
md"""
### Medium and Grid

The distance axis is in kilometers and time in seconds throughout, so
``\rho_0`` and ``\mu_0`` below are expressed in the matching units
(``\mathrm{kg/km^3}`` and ``\mathrm{kg\,km^{-1}\,s^{-2}}``, i.e. a "Pa" built
from km instead of m) rather than plain SI -- chosen to give a shear-wave
speed ``v_{s,0}=\sqrt{\mu_0/\rho_0}\approx5\ \mathrm{km/s}``, typical of
Earth's upper mantle.
"""

# ╔═╡ fec2f28a-ecd0-4bfc-a113-448966ca023d
ρ0 = 3.22e-3 * 1.0e15   # ≈ 3220 kg/m³ (upper-mantle olivine), expressed in kg/km³

# ╔═╡ ac4eb794-7cf7-4fc2-9a9c-adf1cc7eb9ce
μ0 = 82.0e9 * 1.0e3     # ≈ 82 GPa (upper-mantle shear modulus), expressed in kg·km⁻¹·s⁻²

# ╔═╡ 858ed5e1-5881-4bb9-927f-ebe22afc3e5c
md"### Forward Modeling"

# ╔═╡ 134bee99-e691-4e4c-beba-45d3c541bd4a
"""
	step_string!(vy, sigma, mu, rho, dx, dt, boundary)

Advance the staggered-grid velocity `vy` (at integer grid nodes) and stress
`sigma` (at the half-grid nodes between them) by one leapfrog time step,
using the interior stencils for ``\\partial_t\\sigma = \\mu\\,\\partial_x
v_y`` and ``\\rho\\,\\partial_t v_y = \\partial_x \\sigma``.

The two boundary velocity nodes are handled according to `boundary`:
`"free"` mirrors the nearest interior stress across the edge (enforcing
``\\sigma=0`` exactly at the physical boundary via a ghost point), while
`"fixed"` leaves the boundary velocity untouched, i.e. a rigid, immovable
end.
"""
function step_string!(vy, sigma, mu, rho, dx, dt, boundary)
    nx = length(vy)
    @inbounds for i in 1:nx-1
        sigma[i] += dt * mu[i] * (vy[i+1] - vy[i]) / dx
    end
    @inbounds for i in 2:nx-1
        vy[i] += dt / rho[i] * (sigma[i] - sigma[i-1]) / dx
    end
    if boundary == "free"
        vy[1] += dt / rho[1] * (2 * sigma[1]) / dx
        vy[end] += dt / rho[end] * (-2 * sigma[end]) / dx
    end
    return vy, sigma
end

# ╔═╡ 806c8e93-15fb-402a-bab4-c58cae7dd6d2
"""
	simulate_string(; courant, density_ratio, interface_x, boundary, mu0, rho0,
	                  nx=600, xmin=-100.0, xmax=100.0, n_frames=220, duration=nothing)

Run the 1-D velocity-stress staggered-grid simulation of a transverse pulse
launched at `x=-50`, propagating through a medium with a density
discontinuity of `density_ratio` at `x=interface_x` (denser for `x <
interface_x`), reflecting off both ends of the domain according to
`boundary`.

Returns a named tuple with the down-sampled time series `vy_frames`
(`n_frames × nx`), `frames_t`, the spatial grid `xgrid`, the wave speeds and
impedances `vs1`/`vs2`/`Z1`/`Z2` either side of the interface, the analytic
reflection/transmission coefficients `R`/`T`, and `stable::Bool` (`false` if
the Courant condition was violated and the run was aborted early, in which
case the remaining frames repeat the last valid one).
"""
function simulate_string(; courant, density_ratio, interface_x, boundary,
    mu0, rho0, nx=600, xmin=-100.0, xmax=100.0, n_frames=220, duration=nothing)

    xg = range(xmin, xmax, length=nx)
    xgrid = collect(xg)
    dx = step(xg)

    mu = fill(mu0, nx - 1)
    rho = fill(rho0, nx)
    ileft = count(x -> x <= interface_x, xgrid)
    rho[1:ileft] .*= density_ratio

    vs1 = sqrt(mu0 / rho[1])
    vs2 = sqrt(mu0 / rho[end])
    Z1, Z2 = rho[1] * vs1, rho[end] * vs2
    R = (Z1 - Z2) / (Z1 + Z2)
    T = 2 * Z1 / (Z1 + Z2)

    vmax = max(vs1, vs2)
    dt = courant * dx / vmax
    dur = something(duration, 1.6 * (xmax - xmin) / sqrt(mu0 / rho0))
    nt = max(2, round(Int, dur / dt))
    frame_stride = max(1, nt ÷ n_frames)
    nframes_actual = length(1:frame_stride:nt)

    vy = exp.(-0.11 .* (xgrid .+ 50.0) .^ 2)
    sigma = zeros(nx - 1)

    vy_frames = zeros(nframes_actual, nx)
    frames_t = zeros(nframes_actual)
    stable = true
    fidx = 0
    for it in 1:nt
        step_string!(vy, sigma, mu, rho, dx, dt, boundary)
        if maximum(abs, vy) > 1.0e3
            stable = false
            break
        end
        if it % frame_stride == 0
            fidx += 1
            fidx > nframes_actual && break
            vy_frames[fidx, :] .= vy
            frames_t[fidx] = it * dt
        end
    end
    if fidx < nframes_actual
        for k in fidx+1:nframes_actual
            vy_frames[k, :] .= fidx == 0 ? 0.0 : vy_frames[fidx, :]
            frames_t[k] = fidx == 0 ? 0.0 : frames_t[fidx]
        end
    end

    (; xgrid, vy_frames, frames_t, rho, mu, interface_x, vs1, vs2, Z1, Z2, R, T,
        courant, boundary, stable, xmin, xmax, nx)
end

# ╔═╡ ddeef9e3-ba88-4db9-9e53-025b0e9bde36
md"""
### Verifying Reflection and Transmission

`simulate_string`'s symmetric velocity initial condition launches two
counter-propagating pulses, each carrying half the initial peak amplitude (a
standard property of a velocity-only initial condition for this first-order
system: the two characteristics ``v\pm\sigma/Z`` each inherit half of the
initial profile). The check below places the interface close enough to the
source, and stops the run early enough, that neither the outgoing twin pulse
nor either transmitted/reflected pulse has reached a domain edge yet -- so
the only physics left in the isolated snapshot is the interface interaction
itself, and the simulated reflection/transmission amplitudes can be compared
directly against the analytic ``Z``-based formulas above.
"""

# ╔═╡ 149cbe37-8c4f-4cd2-82e7-c938d6865963
let
    ratio_check = 4.0
    iface = -20.0
    vs1c = sqrt(μ0 / (ρ0 * ratio_check))
    Z1c, Z2c = ρ0 * ratio_check * vs1c, ρ0 * sqrt(μ0 / ρ0)
    R_analytic = (Z1c - Z2c) / (Z1c + Z2c)
    T_analytic = 2 * Z1c / (Z1c + Z2c)

    dur = 42.0 / vs1c
    res = simulate_string(; courant=0.3, density_ratio=ratio_check, interface_x=iface,
        boundary="fixed", mu0=μ0, rho0=ρ0, n_frames=300, duration=dur)

    vy_last = res.vy_frames[end, :]
    xg = res.xgrid
    # the symmetric IC also launches an untouched twin pulse heading left from
    # the source, in parallel with (and never interacting with) the interface --
    # it stays at full amplitude, so it must be excluded from the reflected-pulse
    # search window or it dominates the argmax below
    twin_pos = -50.0 - vs1c * dur
    i_refl = findall(x -> (x < iface - 3) && (x > twin_pos + 8), xg)
    i_trans = findall(x -> x > iface + 3, xg)
    idxR = i_refl[argmax(abs.(vy_last[i_refl]))]
    idxT = i_trans[argmax(abs.(vy_last[i_trans]))]

    peak_incident = 0.5 # the symmetric IC splits a unit-peak Gaussian into two counter-propagating waves of half the amplitude -- that's what actually reaches the interface
    sim_R = vy_last[idxR] / peak_incident
    sim_T = vy_last[idxT] / peak_incident

    @assert isapprox(sim_R, R_analytic; atol=0.05) "reflection coefficient mismatch: simulated $(round(sim_R, digits=3)), analytic $(round(R_analytic, digits=3))"
    @assert isapprox(sim_T, T_analytic; atol=0.05) "transmission coefficient mismatch: simulated $(round(sim_T, digits=3)), analytic $(round(T_analytic, digits=3))"

    (simulated_R=round(sim_R, digits=3), analytic_R=round(R_analytic, digits=3),
        simulated_T=round(sim_T, digits=3), analytic_T=round(T_analytic, digits=3))
end

# ╔═╡ 97527c8d-70f6-4665-b0e4-0e7be62a80a4
md"""
### Verifying Boundary Polarity

Both boundaries in `simulate_string` share the same `boundary` type, so a run
long enough for the first reflection to have happened -- but short enough to
avoid a second one arriving from either wall and complicating the sign --
directly tests the free-vs-fixed polarity claim made above: whichever
boundary the dominant reflected pulse actually came from, both walls behave
identically, so the sign of the largest-amplitude sample after one reflection
is an unambiguous readout of the boundary type's effect.
"""

# ╔═╡ 9831ddbe-a0c4-459f-8217-9416b297ae04
let
    check_t = 25.0
    results = map(("fixed", "free")) do bt
        res = simulate_string(; courant=0.3, density_ratio=1.0, interface_x=0.0,
            boundary=bt, mu0=μ0, rho0=ρ0, n_frames=300)
        k = argmin(abs.(res.frames_t .- check_t))
        vyk = res.vy_frames[k, :]
        xg = res.xgrid
        # by t=check_t only the nearer (left) wall has reflected -- the twin
        # heading right is still mid-flight at full, untouched amplitude, so
        # restrict the search to the side that has actually bounced
        i_left = findall(x -> x < 0, xg)
        idx = i_left[argmax(abs.(vyk[i_left]))]
        vyk[idx]
    end
    fixed_peak, free_peak = results

    @assert fixed_peak<0 "expected an inverted (negative) reflection off a fixed end, got $(round(fixed_peak, digits=4))"
    @assert free_peak>0 "expected a same-sign (positive) reflection off a free end, got $(round(free_peak, digits=4))"

    (fixed_end_peak=round(fixed_peak, digits=4), free_end_peak=round(free_peak, digits=4))
end

# ╔═╡ e091f4b2-d23d-4ad1-8820-793368159ec6
md"### The Interactive Widget"

# ╔═╡ ff461695-5ea8-4f6d-8de9-3fa5204b5794
begin
    """
    A dark-canvas widget for the 1-D velocity-stress staggered-grid string
    solver. Courant number, density ratio, interface position, and boundary
    type are all physics inputs -- dragging any of them triggers a fresh
    `simulate_string` solve and a push of the resulting time series back into
    this widget (see the `wos-results` push, right below the widget). Once
    the series is on hand, Play/pause and the time-scrub slider are pure
    client-side state: nothing about *which frame* is currently displayed
    needs another round trip to Julia.
    """
    struct WaveOnStringInput
        courant::Float64
        density_ratio::Float64
        interface_x::Float64
        boundary::String
    end

    WaveOnStringInput(; courant=0.4, density_ratio=4.0, interface_x=-25.0, boundary="free") =
        WaveOnStringInput(Float64(courant), Float64(density_ratio), Float64(interface_x), boundary)

    Base.get(w::WaveOnStringInput) = Dict{String,Any}(
        "courant" => w.courant, "density_ratio" => w.density_ratio,
        "interface_x" => w.interface_x, "boundary" => w.boundary,
    )

    function Base.show(io::IO, ::MIME"text/html", w::WaveOnStringInput)
        write(io, """
<div id="wos-root" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    #wos-root { width: 100%; box-sizing: border-box; color: #d1d5db; font: 14px sans-serif; }
    #wos-root .wos-title { width: 100%; box-sizing: border-box; text-align: center; margin-bottom: 10px;
      background: #0a0f18; border: 1px solid #3b5c85; border-radius: 6px; padding: 10px 14px; }
    #wos-root .wos-title-desc { font-size: 17px; font-weight: 700; color: #e5e7eb; }
    #wos-root .wos-title-hint { font-size: 13px; color: #9ca3af; margin-top: 3px; }
    #wos-root .wos-workspace { display: flex; flex-direction: column; gap: 8px; align-items: center; width: 100%; }
    #wos-root canvas { background: #000; border: 1px solid #374151; border-radius: 6px; display: block; }
    #wos-root .wos-panel-title { font-size: 13px; font-weight: 700; color: #e5e7eb; align-self: flex-start; margin: 4px 0 -2px 2px; }
    #wos-root .wos-controls { width: 100%; margin-top: 10px; display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); gap: 8px; font: 14px sans-serif; }
    #wos-root .wos-control-group { box-sizing: border-box; background: #050505; border: 1px solid #2f3744;
      border-radius: 6px; padding: 10px 12px; }
    #wos-root .wos-control-title { font-weight: 700; color: #e5e7eb; margin-bottom: 8px; font-size: 18px; }
    #wos-root .wos-control-row { display: grid; grid-template-columns: minmax(90px,120px) minmax(70px,1fr) minmax(50px,70px);
      gap: 6px; align-items: center; margin: 6px 0; }
    #wos-root .wos-control-row input[type=range] { width: 100%; min-width: 0; }
    #wos-root .wos-value { color: #d1d5db; text-align: right; font-variant-numeric: tabular-nums; }
    #wos-root .wos-actions { display: flex; gap: 8px; align-items: center; flex-wrap: wrap; }
    #wos-root button { border-radius: 4px; border: 1px solid #9ca3af; background: #606060; color: #f3f4f6;
      padding: 6px 12px; font-size: 14px; cursor: pointer; }
    #wos-root button.wos-active { background: #1d4ed8; border-color: #93c5fd; }
    #wos-root .wos-hint { font-size: 12px; color: #6b7280; margin-top: 6px; }
    #wos-root .wos-panel { position: relative; }
    #wos-root .wos-ratio-input { position: absolute; width: 58px; font-size: 11px; text-align: center;
      pointer-events: auto; background: #0b0b0b; color: #e5e7eb; border: 1px solid #4ade80; border-radius: 3px;
      padding: 2px; font-variant-numeric: tabular-nums; }
  </style>
  <div class="wos-title">
    <div class="wos-title-desc">Drag the interface line and edit the density-ratio box directly on the plot; Courant and boundary type control everything else.</div>
    <div class="wos-title-hint">push Courant above 1 to see it blow up &middot; press Play to watch the pulse split, reflect, and transmit</div>
  </div>
  <div class="wos-workspace">
    <div class="wos-panel-title">Staggered grid (schematic)</div>
    <canvas id="wos-schem"></canvas>
    <div class="wos-panel-title">String velocity v_y(x) -- drag the white interface line, edit the density-ratio box</div>
    <div class="wos-panel">
      <canvas id="wos-main"></canvas>
      <div id="wos-overlay" style="position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none">
        <input type="number" id="wos-ratio-input" class="wos-ratio-input" min="0.2" max="8" step="0.1" value="$(w.density_ratio)" title="density ratio, medium 1 (left of interface) / baseline">
      </div>
    </div>
  </div>
  <div class="wos-controls">
    <div class="wos-control-group">
      <div class="wos-control-title">Numerics</div>
      <label class="wos-control-row"><span>Courant C</span><input type="range" id="wos-courant" min="0.05" max="1.3" step="0.05" value="$(w.courant)"><span id="wos-courant-v" class="wos-value">$(w.courant)</span></label>
      <div class="wos-hint">Density ratio and interface position: drag directly on the plot above.</div>
    </div>
    <div class="wos-control-group">
      <div class="wos-control-title">Boundary</div>
      <div class="wos-actions">
        <button type="button" id="wos-b-free">Free (σ=0)</button>
        <button type="button" id="wos-b-fixed">Fixed (rigid)</button>
      </div>
      <div class="wos-hint">Free = Earth's actual surface. Fixed = a rigid wall.</div>
    </div>
    <div class="wos-control-group">
      <div class="wos-control-title">Playback</div>
      <div class="wos-actions"><button id="wos-play" type="button">Play</button><button id="wos-reset" type="button">Reset defaults</button></div>
      <input type="range" id="wos-time" min="0" max="1000" step="1" value="0" style="width:100%;margin-top:8px">
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.86, par.clientWidth || window.innerWidth*0.86, 1400)
  const totalW = Math.max(700, availW)
  const DPR = Math.min(window.devicePixelRatio || 1, 2)

  let courant = $(w.courant), densityRatio = $(w.density_ratio), interfaceX = $(w.interface_x)
  let boundaryType = "$(w.boundary)"
  let data = null   // filled in by the 'wos-results' push from Julia
  let playing = false, rafId = null, tPhase = 0, lastTs = null
  let draggingIface = false

  const schemCvs = par.querySelector('#wos-schem'), sctx = schemCvs.getContext('2d')
  const mainCvs = par.querySelector('#wos-main'), mctx = mainCvs.getContext('2d')
  const ratioInput = par.querySelector('#wos-ratio-input')

  function hidpi(cv, cx, w, h){
    cv.width = Math.round(w*DPR); cv.height = Math.round(h*DPR)
    cv.style.width = w+'px'; cv.style.height = h+'px'
    cx.setTransform(DPR,0,0,DPR,0,0)
  }
  const SCH = 100, MCH = 380
  hidpi(schemCvs, sctx, totalW, SCH)
  hidpi(mainCvs, mctx, totalW, MCH)

  // Shared main-panel coordinate transform -- used by both drawMain and the
  // interface-drag mouse handlers so a click and its own drawn line always
  // agree on where "the line" actually is.
  const MPAD = {l:48, r:16, t:16, b:54}
  function domainMinMax(){ return data ? [data.xmin, data.xmax] : [-100, 100] }
  function mainX(xv, xmin, xmax){
    const x0 = MPAD.l, x1 = totalW-MPAD.r
    return x0 + ((xv-xmin)/(xmax-xmin))*(x1-x0)
  }
  function mainXInv(px, xmin, xmax){
    const x0 = MPAD.l, x1 = totalW-MPAD.r
    return xmin + ((px-x0)/(x1-x0))*(xmax-xmin)
  }
  function positionRatioInput(){
    const [xmin, xmax] = domainMinMax()
    const lineX = mainX(interfaceX, xmin, xmax)
    const cx = (MPAD.l + lineX)/2
    ratioInput.style.left = Math.round(cx-29) + 'px'
    ratioInput.style.top = (MPAD.t+6) + 'px'
  }

  function drawArrow(cx0,cy0,cx1,cy1,ccx,ccy){
    sctx.beginPath()
    sctx.moveTo(cx0,cy0)
    sctx.quadraticCurveTo(ccx,ccy,cx1,cy1)
    sctx.lineWidth = 1.4
    sctx.stroke()
    const ang = Math.atan2(cy1-ccy, cx1-ccx)
    sctx.beginPath()
    sctx.moveTo(cx1,cy1)
    sctx.lineTo(cx1-7*Math.cos(ang-0.4), cy1-7*Math.sin(ang-0.4))
    sctx.lineTo(cx1-7*Math.cos(ang+0.4), cy1-7*Math.sin(ang+0.4))
    sctx.closePath(); sctx.fill()
  }

  function drawSchematic(){
    sctx.clearRect(0,0,totalW,SCH)
    const n = 5
    const padL = 60, padR = 60
    const usableW = totalW - padL - padR
    const xAt = i => padL + (i/(n-1))*usableW
    const sxAt = i => (xAt(i)+xAt(i+1))/2
    const midY = SCH*0.5

    sctx.strokeStyle = '#374151'; sctx.beginPath(); sctx.moveTo(padL-14,midY); sctx.lineTo(totalW-padR+14,midY); sctx.stroke()

    sctx.strokeStyle = '#f97316'; sctx.fillStyle = '#f97316'
    for(let i=1;i<n-2;i++){
      drawArrow(xAt(i), midY, sxAt(i), midY-24, xAt(i)+(sxAt(i)-xAt(i))*0.5, midY-24)
      drawArrow(xAt(i+1), midY, sxAt(i), midY-24, xAt(i+1)-(xAt(i+1)-sxAt(i))*0.5, midY-24)
    }
    sctx.strokeStyle = '#38bdf8'; sctx.fillStyle = '#38bdf8'
    for(let i=2;i<n-1;i++){
      drawArrow(sxAt(i-1), midY, xAt(i), midY+24, sxAt(i-1)+(xAt(i)-sxAt(i-1))*0.5, midY+24)
      drawArrow(sxAt(i), midY, xAt(i), midY+24, sxAt(i)-(sxAt(i)-xAt(i))*0.5, midY+24)
    }

    for(let i=0;i<n;i++){
      sctx.fillStyle = '#38bdf8'
      sctx.beginPath(); sctx.arc(xAt(i), midY, 5, 0, Math.PI*2); sctx.fill()
      sctx.strokeStyle = '#0b1220'; sctx.lineWidth=1; sctx.stroke()
    }
    for(let i=0;i<n-1;i++){
      sctx.save(); sctx.translate(sxAt(i), midY); sctx.rotate(Math.PI/4)
      sctx.fillStyle = '#f97316'; sctx.fillRect(-4,-4,8,8)
      sctx.restore()
    }

    sctx.fillStyle = '#9ca3af'; sctx.font = '11px sans-serif'; sctx.textAlign = 'center'
    sctx.fillText('v_y  (velocity, node i)', xAt(2), midY+40)
    sctx.fillText('σ  (stress, half-node i+1/2)', sxAt(2), midY-34)
    sctx.textAlign = 'left'
    sctx.fillStyle = '#6b7280'; sctx.font = '10px sans-serif'
    sctx.fillText('orange: σ update reads neighboring v_y   ·   blue: v_y update reads neighboring σ', padL-14, SCH-6)
  }

  function drawBoundaryMarker(px, side){
    const y0=32, y1=MCH-52
    // labelY sits just below the Z/R/T/C readout row and well above the
    // x-axis tick labels at the bottom -- both would otherwise collide with
    // this label at the domain edge, where it's drawn directly under the tick
    const labelY = y0 + 12
    const label = boundaryType === 'fixed' ? 'fixed' : 'free'
    const color = boundaryType === 'fixed' ? '#9ca3af' : '#4ade80'
    if(boundaryType === 'fixed'){
      mctx.strokeStyle = color; mctx.lineWidth = 1.3
      for(let y=y0; y<y1; y+=9){
        mctx.beginPath()
        const dx = side==='left' ? 9 : -9
        mctx.moveTo(px, y); mctx.lineTo(px+dx, y+9)
        mctx.stroke()
      }
    } else {
      mctx.strokeStyle = color; mctx.lineWidth = 1.6
      const midy = (y0+y1)/2, dx = side==='left' ? -13 : 13
      mctx.beginPath(); mctx.moveTo(px,midy-16); mctx.lineTo(px+dx,midy-16); mctx.stroke()
      mctx.beginPath(); mctx.moveTo(px,midy+16); mctx.lineTo(px+dx,midy+16); mctx.stroke()
    }
    mctx.fillStyle = color; mctx.font = '11px sans-serif'; mctx.textAlign = 'left'
    const tw = mctx.measureText(label).width
    const rawX = side==='left' ? px+4 : px-4-tw
    const clampedX = Math.max(2, Math.min(totalW-tw-2, rawX))
    mctx.fillText(label, clampedX, labelY)
  }

  function currentFrameIdx(){
    if(!data) return 0
    const nF = data.nFrames
    const frameDT = nF>1 ? (data.framesT[1]-data.framesT[0]) : 1
    let idx = Math.round(tPhase/frameDT)
    return Math.max(0, Math.min(nF-1, idx))
  }

  function drawMain(){
    mctx.clearRect(0,0,totalW,MCH)
    const x0=MPAD.l, x1=totalW-MPAD.r, y0=MPAD.t, y1=MCH-MPAD.b
    positionRatioInput()

    if(!data){
      mctx.fillStyle = '#6b7280'; mctx.font = '13px sans-serif'; mctx.fillText('computing...', 12, 20)
      return
    }

    const {xgrid, nx, nFrames, vyFlat, framesT, vs1, vs2, Z1, Z2, R, T, stable, xmin, xmax, ampMax} = data
    const ifaceX = interfaceX

    const X = xv => mainX(xv, xmin, xmax)
    const amp = Math.max(ampMax, 1e-9)
    const Y = v => (y0+y1)/2 - (v/amp)*((y1-y0)/2)*0.92

    mctx.fillStyle = 'rgba(56,189,248,0.06)'; mctx.fillRect(x0, y0, X(ifaceX)-x0, y1-y0)
    mctx.fillStyle = 'rgba(74,222,128,0.06)'; mctx.fillRect(X(ifaceX), y0, x1-X(ifaceX), y1-y0)
    mctx.strokeStyle = '#374151'; mctx.strokeRect(x0,y0,x1-x0,y1-y0)
    mctx.strokeStyle = '#1f2937'; mctx.beginPath(); mctx.moveTo(x0,Y(0)); mctx.lineTo(x1,Y(0)); mctx.stroke()

    if(!stable){
      mctx.fillStyle = 'rgba(127,29,29,0.4)'; mctx.fillRect(x0,y0,x1-x0,y1-y0)
      mctx.fillStyle = '#fca5a5'; mctx.font = 'bold 14px sans-serif'; mctx.textAlign = 'center'
      mctx.fillText('⚠ UNSTABLE  —  C = '+courant.toFixed(2)+' > 1', (x0+x1)/2, (y0+y1)/2-8)
      mctx.font = '12px sans-serif'
      mctx.fillText('lower the Courant slider below 1', (x0+x1)/2, (y0+y1)/2+12)
      mctx.textAlign = 'left'
      return
    }

    // draggable directly (mousedown/mousemove below) -- brighter while held
    mctx.strokeStyle = draggingIface ? '#ffffff' : '#e5e7eb'; mctx.lineWidth = draggingIface ? 2.6 : 1.6
    mctx.beginPath(); mctx.moveTo(X(ifaceX),y0); mctx.lineTo(X(ifaceX),y1); mctx.stroke()

    drawBoundaryMarker(x0, 'left')
    drawBoundaryMarker(x1, 'right')

    const fIdx = currentFrameIdx()
    mctx.strokeStyle = '#facc15'; mctx.lineWidth = 2
    mctx.beginPath()
    for(let i=0;i<nx;i++){
      const px = X(xgrid[i]), py = Y(vyFlat[fIdx*nx+i])
      i===0 ? mctx.moveTo(px,py) : mctx.lineTo(px,py)
    }
    mctx.stroke()

    mctx.fillStyle = '#9ca3af'; mctx.font = '11px sans-serif'; mctx.textAlign = 'center'
    for(let k=0;k<=4;k++){
      const xv = xmin + k*(xmax-xmin)/4, px = X(xv)
      mctx.strokeStyle = '#374151'; mctx.beginPath(); mctx.moveTo(px,y1); mctx.lineTo(px,y1+4); mctx.stroke()
      mctx.fillText(Math.round(xv)+'', px, y1+16)
    }
    mctx.fillText('distance (km)', (x0+x1)/2, y1+32)

    mctx.textAlign = 'left'; mctx.fillStyle = '#d1d5db'; mctx.font = '12px sans-serif'
    mctx.fillText('Z₁='+Z1.toExponential(2)+'   Z₂='+Z2.toExponential(2)+'   R='+R.toFixed(2)+'   T='+T.toFixed(2), x0, 12)
    mctx.textAlign = 'right'
    mctx.fillStyle = courant > 1 ? '#f87171' : '#9ca3af'
    mctx.fillText('C='+courant.toFixed(2)+'   t='+framesT[fIdx].toFixed(1)+'s  ('+(fIdx+1)+'/'+nFrames+')', x1, 12)
    mctx.textAlign = 'left'
  }

  function updateTimeSlider(){
    if(!data) return
    const totalDur = data.framesT[data.nFrames-1]
    par.querySelector('#wos-time').value = totalDur>0 ? Math.round(1000*tPhase/totalDur) : 0
  }

  function emit(){
    par.value = {courant, density_ratio: densityRatio, interface_x: interfaceX, boundary: boundaryType}
    par.dispatchEvent(new CustomEvent('input'))
  }

  function syncBoundaryButtons(){
    par.querySelector('#wos-b-free').classList.toggle('wos-active', boundaryType==='free')
    par.querySelector('#wos-b-fixed').classList.toggle('wos-active', boundaryType==='fixed')
  }

  par.querySelector('#wos-courant').addEventListener('input', e=>{
    courant = parseFloat(e.target.value); par.querySelector('#wos-courant-v').textContent = courant.toFixed(2); drawMain(); emit()
  })
  ratioInput.addEventListener('input', e=>{
    const v = parseFloat(e.target.value)
    if(!isFinite(v) || v<=0) return
    densityRatio = Math.max(0.2, Math.min(8, v))
    drawMain(); emit()
  })

  // Interface line: drag it directly on the canvas rather than a slider,
  // same pattern as the layer-boundary drag in wave-mode-duality-1D.jl.
  function nearIfaceLine(px){
    const [xmin,xmax] = domainMinMax()
    return Math.abs(px - mainX(interfaceX, xmin, xmax)) < 10
  }
  mainCvs.addEventListener('mousedown', e=>{
    if(!data) return
    const rect = mainCvs.getBoundingClientRect()
    if(nearIfaceLine(e.clientX-rect.left)) draggingIface = true
  })
  mainCvs.addEventListener('mousemove', e=>{
    if(draggingIface || !data) return
    const rect = mainCvs.getBoundingClientRect()
    mainCvs.style.cursor = nearIfaceLine(e.clientX-rect.left) ? 'ew-resize' : 'default'
  })
  window.addEventListener('mousemove', e=>{
    if(!draggingIface || !data) return
    const rect = mainCvs.getBoundingClientRect()
    const px = Math.max(0, Math.min(totalW, e.clientX-rect.left))
    const [xmin,xmax] = domainMinMax()
    let xv = mainXInv(px, xmin, xmax)
    interfaceX = Math.max(xmin+20, Math.min(xmax-20, xv))
    drawMain(); emit()
  })
  window.addEventListener('mouseup', ()=>{
    if(draggingIface){ draggingIface = false; drawMain() }
  })

  par.querySelector('#wos-b-free').addEventListener('click', ()=>{ boundaryType='free'; syncBoundaryButtons(); drawMain(); emit() })
  par.querySelector('#wos-b-fixed').addEventListener('click', ()=>{ boundaryType='fixed'; syncBoundaryButtons(); drawMain(); emit() })

  const playBtn = par.querySelector('#wos-play')
  function stepAnim(ts){
    if(lastTs===null) lastTs = ts
    const ddt = Math.min(0.1, (ts-lastTs)/1000)
    lastTs = ts
    if(data){
      const totalDur = data.framesT[data.nFrames-1]
      const simSpeed = totalDur > 0 ? totalDur/10 : 1
      tPhase += ddt*simSpeed
      if(tPhase > totalDur) tPhase = 0
    }
    drawMain(); updateTimeSlider()
    rafId = requestAnimationFrame(stepAnim)
  }
  playBtn.addEventListener('click', ()=>{
    playing = !playing
    playBtn.textContent = playing ? 'Pause' : 'Play'
    if(playing){ lastTs=null; rafId = requestAnimationFrame(stepAnim) }
    else if(rafId){ cancelAnimationFrame(rafId); rafId = null }
  })

  par.querySelector('#wos-time').addEventListener('input', e=>{
    if(playing){ playing=false; playBtn.textContent='Play'; if(rafId){ cancelAnimationFrame(rafId); rafId=null } }
    if(data){
      const totalDur = data.framesT[data.nFrames-1]
      tPhase = (parseFloat(e.target.value)/1000)*totalDur
    }
    drawMain()
  })

  par.querySelector('#wos-reset').addEventListener('click', ()=>{
    courant = $(w.courant); densityRatio = $(w.density_ratio); interfaceX = $(w.interface_x); boundaryType = "$(w.boundary)"
    par.querySelector('#wos-courant').value = courant; par.querySelector('#wos-courant-v').textContent = courant.toFixed(2)
    ratioInput.value = densityRatio
    syncBoundaryButtons()
    tPhase = 0; playing = false; if(rafId){ cancelAnimationFrame(rafId); rafId=null }; playBtn.textContent = 'Play'
    drawMain(); emit()
  })

  window.addEventListener('wos-results', e=>{
    const d = e.detail ? JSON.parse(e.detail) : null
    if(!d) return
    data = d
    tPhase = 0
    drawMain(); updateTimeSlider()
  })

  syncBoundaryButtons()
  drawSchematic(); drawMain(); emit()
</script>
""")
    end

    const _wos_ready = true
end

# ╔═╡ bc9056a7-6c7f-4664-8167-d6ab9792563c
begin
    # `WaveOnStringInput` is defined in the Appendix, displayed below this cell
    # -- a bare reference forces Pluto to run that cell first on a cold restart.
    # See "the one thing that will silently break on a fresh restart" in
    # pluto-widget-SKILL.md.
    _wos_ready
    WideCell(@bind _wos WaveOnStringInput(); max_width=1400)
end

# ╔═╡ c6cd50d4-b42a-4aec-a6b1-b16d88940d9f
begin
    # Courant number, density ratio, interface position, and boundary type are
    # all physics inputs -- each one triggers the Julia recompute below.
    # Which already-computed frame is on screen (Play/pause, the scrub slider)
    # is pure client-side state inside the widget; see its docstring.
    wos_courant = Float64(_wos["courant"])
    wos_density_ratio = Float64(_wos["density_ratio"])
    wos_interface_x = Float64(_wos["interface_x"])
    wos_boundary = String(_wos["boundary"])
end

# ╔═╡ b19d0fbf-d0ae-4472-aaca-89369188d6a3
wos_result = simulate_string(; courant=wos_courant, density_ratio=wos_density_ratio,
    interface_x=wos_interface_x, boundary=wos_boundary, mu0=μ0, rho0=ρ0)

# ╔═╡ 4188a0d9-101a-474a-b9d5-69eb8808e24a
# Push the down-sampled time series to the widget above -- it stays mounted
# across reruns of this cell, same CustomEvent pattern used throughout this
# repo's other widgets. This is the only place the widget's Julia side runs
# after a slider change: playback itself never touches Julia again.
let
    res = wos_result
    nx = length(res.xgrid)
    nFrames = length(res.frames_t)
    ampMax = maximum(abs, res.vy_frames)

    num(x) = isfinite(x) ? string(round(Float64(x), digits=6)) : "0"
    jsonarr(v) = "[" * join(num.(v), ",") * "]"
    jsonflat(m) = "[" * join(num.(vec(permutedims(m))), ",") * "]"

    payload = string(
        "{\"xgrid\":", jsonarr(res.xgrid),
        ",\"nx\":", nx,
        ",\"nFrames\":", nFrames,
        ",\"vyFlat\":", jsonflat(res.vy_frames),
        ",\"framesT\":", jsonarr(res.frames_t),
        ",\"vs1\":", num(res.vs1), ",\"vs2\":", num(res.vs2),
        ",\"Z1\":", num(res.Z1), ",\"Z2\":", num(res.Z2),
        ",\"R\":", num(res.R), ",\"T\":", num(res.T),
        ",\"stable\":", res.stable ? "true" : "false",
        ",\"xmin\":", num(res.xmin), ",\"xmax\":", num(res.xmax),
        ",\"ampMax\":", num(ampMax),
        "}",
    )
    HTML("""<script>
      window.dispatchEvent(new CustomEvent('wos-results', {detail: $(repr(payload))}));
    </script>""")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "4c0a816760ff43ca9ffe29e328a423fcbd2cbca4"

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
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.0"

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
# ╠═46dca275-df8b-4060-b81d-841428142791
# ╟─f3c52928-1450-4596-87b7-b6665a1b8f6a
# ╟─bc9056a7-6c7f-4664-8167-d6ab9792563c
# ╠═c6cd50d4-b42a-4aec-a6b1-b16d88940d9f
# ╠═b19d0fbf-d0ae-4472-aaca-89369188d6a3
# ╟─4188a0d9-101a-474a-b9d5-69eb8808e24a
# ╟─f69bc5a5-d354-41b6-8d12-61a04ef1c43d
# ╟─7676054a-2e9a-43f4-ac46-af162541e269
# ╟─d0be7679-5df5-4ad0-93d1-8df6f94ba299
# ╟─6feab4ba-89ed-4be0-87e7-c1cc1d78b6c2
# ╟─260264ff-24bb-425f-abb8-d3556a1bb19f
# ╟─8a11db06-56c3-42cc-b786-d22b5387bbeb
# ╠═dbfb7f10-8fe2-471b-8626-7a8894626caf
# ╟─90e40955-3aaa-4efb-b48b-21ae94a93bd5
# ╠═fec2f28a-ecd0-4bfc-a113-448966ca023d
# ╠═ac4eb794-7cf7-4fc2-9a9c-adf1cc7eb9ce
# ╟─858ed5e1-5881-4bb9-927f-ebe22afc3e5c
# ╠═134bee99-e691-4e4c-beba-45d3c541bd4a
# ╠═806c8e93-15fb-402a-bab4-c58cae7dd6d2
# ╟─ddeef9e3-ba88-4db9-9e53-025b0e9bde36
# ╠═149cbe37-8c4f-4cd2-82e7-c938d6865963
# ╟─97527c8d-70f6-4665-b0e4-0e7be62a80a4
# ╠═9831ddbe-a0c4-459f-8217-9416b297ae04
# ╟─e091f4b2-d23d-4ad1-8820-793368159ec6
# ╠═ff461695-5ea8-4f6d-8de9-3fa5204b5794
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
