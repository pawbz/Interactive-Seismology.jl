### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> tags = ["faulting", "friction", "earthquake-cycle"]
#> title = "Stick-Slip Faulting — A Loaded Spring on a Frictional Surface"
#> description = "A rate-and-state spring-slider that shows why some faults creep and others lock, load, and rupture."
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

# ╔═╡ 32526135-0ec9-4865-bc05-fa81450ff3a3
begin
    using PlutoUI
    using Printf
end

# ╔═╡ 7f48d81c-5664-459a-bf10-c265462d57a0
PlutoUI.TableOfContents(include_definitions=true)

# ╔═╡ 1103da07-e3f7-408b-937e-5ee0a1563ac7
md"""
# Stick-Slip: Why a Locked Fault Doesn't Creep Smoothly

A tectonic fault is loaded by steady plate motion, millimeters a year, forever. The naive
expectation is that it should simply creep at that same steady rate — a smooth, continuous
slip that never lets extra stress build up. Real faults mostly don't do this. Stress
accumulates for decades to centuries while the fault sits locked, then releases in seconds:
an earthquake.

Whether a patch of fault creeps or locks-and-slips turns out not to depend on how hard it is
loaded, but on a property of the friction itself, combined with how stiff the surrounding
rock is. This notebook builds the simplest system that shows the effect: one block, one
spring, one frictional surface.
"""

# ╔═╡ 8417cf97-4555-44c6-9f44-a6a32799bde8
md"""
## Drag the sliders, or try a preset

The right-hand panel shows the slip rate ``v`` (log scale — it swings over orders of
magnitude between the locked and slipping phases) and the friction coefficient ``\mu``, both
against time. Press **Play** to watch loading turn into a slip event.

!!! tip "A useful first experiment"
	Start with the **Locked fault** preset and press Play — watch the long quiet buildup,
	then the sudden spike. Then switch to **Creeping fault**: same friction law, same
	loading, but the sign of ``b-a`` alone changes the outcome completely.
"""

# ╔═╡ d723bd18-4389-445e-b747-fab5684f85a7
md"""
## From a Force Balance to Two Coupled Rates

The spring is stretched by the plate moving at rate ``V_{pl}`` and shortened by the block's
own slip ``\delta``; a small radiation-damping term ``\eta v`` stands in for the seismic
energy a real fault would radiate away during rapid slip (without it, this quasi-static
approximation would let the slip rate run away to infinity instead of a large-but-finite
peak). The force balance is

``k(V_{pl}t - \delta) - \eta v = \sigma\mu(v,\theta)``, with ``v = \dot\delta``.

Friction depends on the current slip rate *and* on a state variable ``\theta`` that evolves
more slowly — physically, something like the average age of the microscopic contacts making
up the sliding surface (see [`rate_state_friction`](@ref) in the Appendix):

``\mu(v,\theta) = \mu_* + a\ln(v/v_*) + b\ln(v_*\theta/D_c)``, evolving by the aging law
``\dot\theta = 1 - v\theta/D_c``.

Differentiating the force balance in time and solving for ``\dot v`` gives a two-variable ODE
for ``(v,\theta)`` — worked out in full in [`spring_slider_rhs`](@ref), using ``\psi=\ln v``
so the numerical solution can never accidentally go negative. Everything the widget shows
comes from integrating that ODE for the parameters you set.
"""

# ╔═╡ d54456fb-3b04-4a9a-ba51-3df430f0837f
md"""
## Two "Aha"s: A Sign, and a Ratio

**Aha #1 — the sign of ``b-a`` decides whether instability is even possible.** If ``b\le a``
(*velocity-strengthening*), speeding up strengthens the fault immediately, which is
self-correcting: a patch of fault with this property can only creep steadily, no matter how
soft the surrounding rock is. Self-check 1 in the Appendix confirms this numerically. Only
``b>a`` (*velocity-weakening*) opens the door to instability.

**Aha #2 — when ``b>a``, stiffness decides the rest.** Linearizing the ODE around steady
sliding ``(v,\theta)=(V_{pl}, D_c/V_{pl})`` gives a critical stiffness (derived in
[`critical_stiffness`](@ref)): ``k_c = (\sigma(b-a) - \eta V_{pl})/D_c``.

A stiff spring (``k>k_c``) damps perturbations back to steady sliding — the fault creeps. A
soft spring (``k<k_c``) cannot; perturbations grow into a self-sustaining stick-slip limit
cycle — the fault locks, loads, and ruptures, repeatedly, forever, under the *same* friction
law. Self-check 2 confirms both outcomes from the same friction parameters, changed only by a
factor of four in ``k``. The widget's **k / k_c** readout and regime badge track exactly this
ratio live as you drag.
"""

# ╔═╡ 20022d79-0eb9-4e08-9a07-bfde9a9a4db5
md"""
## The Same Idea, at Fault Scale

Real faults are not one block on one spring, but the logic transfers directly: ``k`` stands
in for the elastic stiffness of the rock surrounding a fault patch (roughly the shear modulus
divided by the patch size — a *larger* patch behaves like a *softer* spring), and ``a-b`` is
measured in the lab from how a rock's frictional strength responds to a step change in
sliding speed. Sections of a fault where ``a-b>0`` — often the deeper, hotter parts, or
serpentinite-rich creeping sections like part of the San Andreas Fault — creep steadily and
rarely host large earthquakes. Sections where ``a-b<0`` and the surrounding rock is soft
enough (``k<k_c``) are exactly the locked, seismogenic patches that store centuries of strain
and release it as earthquakes.
"""

# ╔═╡ 874a667c-c6cd-4797-818e-01b8dfc57e1b
md"""
## Appendix

### Stick-Slip: A Loaded Spring on a Frictional Surface

Every function below is plain Julia, documented and — where the result is a genuine physical
claim rather than bookkeeping — checked numerically against a known limit. The widget itself
(struct, `Base.get`, `Base.show`) is defined last, since everything above it is one of its
dependencies.
"""

# ╔═╡ 69a1dd35-07d4-4d7a-8fcf-6dfc00d2e092
"""
    rate_state_friction(v, theta; a, b, mustar, vstar, Dc)

Dieterich–Ruina rate-and-state friction coefficient: a logarithmic direct effect in slip rate
`v` (coefficient `a`) plus a logarithmic evolution effect in the state variable `theta`
(coefficient `b`), relative to reference values `vstar`/`mustar`. `theta` has units of time;
the aging law in [`spring_slider_rhs`](@ref) is what evolves it.
"""
rate_state_friction(v, theta; a, b, mustar, vstar, Dc) =
    mustar + a * log(v / vstar) + b * log(vstar * theta / Dc)

# ╔═╡ 4133723b-9d54-403e-84cb-92f268523f99
"""
    spring_slider_rhs(psi, theta; k, Vpl, sigma, eta, a, b, Dc)

Right-hand side of the quasi-dynamic spring-slider ODE in `(psi, theta)` coordinates, where
`psi = log(v)` (so `v = exp(psi)` stays positive no matter how many decades it swings during
a slip event). Derived by differentiating the force balance
`k*(Vpl*t - delta) - eta*v = sigma*mu(v,theta)` in `t` — using
`d(mu)/dt = a*psidot + b*thetadot/theta` — and solving for `psidot`:
``\\dot\\theta = 1 - v\\theta/D_c``, ``\\dot\\psi = (k(V_{pl}-v) - \\sigma b\\,\\dot\\theta/\\theta)/(\\sigma a + \\eta v)``.

Returns `(psidot, thetadot)`. `eta` is the radiation-damping coefficient that regularizes the
quasi-static approximation (see the "Governing Equations" section above).
"""
function spring_slider_rhs(psi, theta; k, Vpl, sigma, eta, a, b, Dc)
    v = exp(psi)
    thetadot = 1 - v * theta / Dc
    psidot = (k * (Vpl - v) - sigma * b * thetadot / theta) / (sigma * a + eta * v)
    return psidot, thetadot
end

# ╔═╡ 9b71d1c6-70b5-40fa-9e43-c3c02b4082e0
"""
    critical_stiffness(; sigma, a, b, Dc, eta, Vpl)

Linear-stability threshold for the spring-slider's steady-sliding fixed point
`(v, theta) = (Vpl, Dc/Vpl)`. Linearizing [`spring_slider_rhs`](@ref) around that point gives
a 2×2 Jacobian with `det(J) = k*Vpl^2 / (Dc*(sigma*a + eta*Vpl)) > 0` always (`k, Dc, Vpl > 0`
and `sigma*a + eta*Vpl > 0`), so stability is decided entirely by the trace:
`trace(J) < 0 ⟺ k > (sigma*(b-a) - eta*Vpl)/Dc`. That threshold is `k_c` here. When `b ≤ a`
(velocity-strengthening), `k_c ≤ 0` and the fixed point is stable for *any* positive
stiffness — there is no stick-slip instability to speak of. Radiation damping `eta` lowers
`k_c` slightly relative to the textbook `sigma*(b-a)/Dc` result (recovered when `eta=0`),
since it adds dissipation that helps stabilize the system.
"""
critical_stiffness(; sigma, a, b, Dc, eta, Vpl) = (sigma * (b - a) - eta * Vpl) / Dc

# ╔═╡ 48043f9d-d937-443d-9132-8c8c921a5d86
"""
    rk4_step(f, u, h)

One classical 4th-order Runge–Kutta step for `du/dt = f(u)` (`f` already captures the fixed
parameters via a closure). Deliberately the textbook 4-stage tableau, not a higher-order
method, so its coefficients need no citation and cannot be misremembered.
"""
function rk4_step(f, u, h)
    k1 = f(u)
    k2 = f(u .+ (h / 2) .* k1)
    k3 = f(u .+ (h / 2) .* k2)
    k4 = f(u .+ h .* k3)
    return u .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

# ╔═╡ cabca0a8-ffb6-40aa-ab08-d507b0e0e063
"""
    integrate_adaptive(f, u0, tspan; h0, tol, hmin, hmax)

Adaptive-step integration of `du/dt = f(u)` by step-doubling Richardson extrapolation: each
candidate step is taken once at size `h` and again as two steps of size `h/2`; since
[`rk4_step`](@ref) is 4th order, the leading error term scales as `h^5` and
`(u_half - u_full)/15` estimates it. A step is accepted (keeping the more accurate `u_half`)
when that estimate is within `tol` of the state's own scale, and the next step size grows or
shrinks accordingly. This adaptivity matters here because a stick-slip cycle spends almost
all of its time in slow interseismic loading and a tiny fraction in a fast slip event — a
fixed step size would either waste effort during the quiet phase or miss the event entirely.

Returns `(ts, us)`: a `Vector{Float64}` of times and a `Vector{Vector{Float64}}` of states.
"""
function integrate_adaptive(f, u0, tspan; h0=1.0e-3, tol=1.0e-7, hmin=1.0e-10, hmax=0.5)
    t = tspan[1]
    u = collect(float.(u0))
    ts = [t]
    us = [copy(u)]
    h = h0
    while t < tspan[2]
        h = min(h, tspan[2] - t)
        u_full = rk4_step(f, u, h)
        u_half = rk4_step(f, u, h / 2)
        u_half = rk4_step(f, u_half, h / 2)
        err = maximum(abs.(u_half .- u_full)) / 15
        scale = max(maximum(abs.(u_half)), 1.0e-6)
        if err <= tol * scale || h <= hmin
            t += h
            u = u_half
            push!(ts, t)
            push!(us, copy(u))
            growth = err > 0 ? 0.9 * (tol * scale / err)^0.2 : 2.0
            h = min(hmax, h * clamp(growth, 0.2, 2.0))
        else
            h = max(hmin, h / 2)
        end
    end
    return ts, us
end

# ╔═╡ 0e435cbe-662d-4b75-9530-0d0fa898399a
"""
    integrate_spring_slider(; a, b, Dc, Vpl, sigma, eta, k, mustar=0.6, vstar=Vpl, npoints=2000)

Integrate the rate-and-state spring-slider from a small perturbation off steady sliding for
`250*Dc/Vpl` time units (250 "loading times" — enough to show several stick-slip cycles when
the parameters are unstable, or clean convergence to steady creep when they are stable), then
resample onto a uniform time grid of `npoints` samples. [`integrate_adaptive`](@ref)'s own
time steps are wildly uneven (long strides through the stick phase, tiny ones through a slip
event); a uniform grid is what the widget's JS side plays back at a constant frame rate.

Returns a named tuple `(t, v, theta, mu, stretch, kc)` where `stretch = Vpl.*t .- slip` is the
spring's elastic stretch — the physically bounded, oscillating quantity the widget animates.
Slip itself grows without bound, but stretch is proportional to the spring force (hence to
stress) and saws between the stick and slip phases.
"""
function integrate_spring_slider(; a, b, Dc, Vpl, sigma, eta, k, mustar=0.6, vstar=Vpl, npoints=2000)
    theta_ss = Dc / Vpl
    u0 = (log(Vpl) + 0.05, theta_ss * 0.9)
    f(u) = collect(spring_slider_rhs(u[1], u[2]; k, Vpl, sigma, eta, a, b, Dc))
    tend = 250 * Dc / Vpl
    ts, us = integrate_adaptive(f, u0, (0.0, tend))

    tgrid = collect(range(0.0, tend; length=npoints))
    v = similar(tgrid)
    theta = similar(tgrid)
    j = 1
    for (i, tq) in enumerate(tgrid)
        while j < length(ts) && ts[j+1] <= tq
            j += 1
        end
        j2 = min(j + 1, length(ts))
        frac = ts[j2] > ts[j] ? (tq - ts[j]) / (ts[j2] - ts[j]) : 0.0
        v[i] = exp(us[j][1] + frac * (us[j2][1] - us[j][1]))
        theta[i] = us[j][2] + frac * (us[j2][2] - us[j][2])
    end

    dt = tgrid[2] - tgrid[1]
    slip = similar(tgrid)
    slip[1] = 0.0
    for i in 2:npoints
        slip[i] = slip[i-1] + 0.5 * (v[i-1] + v[i]) * dt
    end

    mu = rate_state_friction.(v, theta; a, b, mustar, vstar, Dc)
    stretch = Vpl .* tgrid .- slip
    kc = critical_stiffness(; sigma, a, b, Dc, eta, Vpl)
    return (t=tgrid, v=v, theta=theta, mu=mu, stretch=stretch, kc=kc)
end

# ╔═╡ 6df71f7e-1f23-4edb-8878-e07aba98a32d
_ss_flatten(x) = join(x, ",")

# ╔═╡ a3fd5f1a-90a3-44cd-b1d4-cdbd9b22b6c5
begin
    struct StickSlipPush
        t::Any
        v::Any
        theta::Any
        mu::Any
        stretch::Any
        kc::Any
    end
    function Base.show(io::IO, ::MIME"text/html", p::StickSlipPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('ss-widget');
        if(w){
          w.dispatchEvent(new CustomEvent('ss-push', { detail: {
            t: [$(_ss_flatten(p.t))],
            v: [$(_ss_flatten(p.v))],
            theta: [$(_ss_flatten(p.theta))],
            mu: [$(_ss_flatten(p.mu))],
            stretch: [$(_ss_flatten(p.stretch))],
            kc: $(p.kc),
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 331f1292-1115-444f-931f-ac46bf88943f
md"""
### Verifying the Spring-Slider Model
"""

# ╔═╡ e63b931a-94d9-42f1-a6c6-e8143faa0911
begin
    _check1_params = (a=0.015, b=0.010, Dc=1.0, Vpl=1.0, sigma=1.0, eta=0.002, k=0.01)
    _check1_f(u) = collect(spring_slider_rhs(u[1], u[2]; _check1_params...))
    _check1_theta_ss = _check1_params.Dc / _check1_params.Vpl
    _check1_ts, _check1_us = integrate_adaptive(_check1_f, (log(0.3), 2 * _check1_theta_ss), (0.0, 40.0))
    _check1_v_end = exp(_check1_us[end][1])
    _check1_theta_end = _check1_us[end][2]
    _check1_v_err = abs(_check1_v_end - _check1_params.Vpl)
    _check1_theta_err = abs(_check1_theta_end - _check1_theta_ss)
    _check1_pass = _check1_v_err < 1.0e-4 && _check1_theta_err < 1.0e-4
    md"""
    **Steady creep** (velocity-strengthening, ``a`` = $(_check1_params.a) > ``b`` = $(_check1_params.b)):
    starting well off equilibrium, the slip rate settles to ``v`` = $(@sprintf("%.6f", _check1_v_end))
    (target ``V_{pl}`` = $(_check1_params.Vpl), error $(@sprintf("%.1e", _check1_v_err))) and the
    state variable to ``\theta`` = $(@sprintf("%.6f", _check1_theta_end)) (target
    ``D_c/V_{pl}`` = $(@sprintf("%.6f", _check1_theta_ss)), error $(@sprintf("%.1e", _check1_theta_err)))
    — the analytically known steady-sliding fixed point.
    $(_check1_pass ? "✅ PASS" : "❌ FAIL")
    """
end

# ╔═╡ b27449d6-07fd-411b-af09-88f0524939e0
begin
    _check2_base = (a=0.010, b=0.016, Dc=1.0, Vpl=1.0, sigma=1.0, eta=0.002)
    _check2_kc = critical_stiffness(; _check2_base...)
    function _check2_late_amplitude(kfac)
        p = merge(_check2_base, (k=kfac * _check2_kc,))
        f(u) = collect(spring_slider_rhs(u[1], u[2]; p...))
        theta_ss = p.Dc / p.Vpl
        ts, us = integrate_adaptive(f, (log(p.Vpl) + 0.05, theta_ss * 0.9), (0.0, 400.0))
        vs = [exp(u[1]) for u in us]
        tail = vs[max(1, length(vs) - 200):end]
        return maximum(tail) - minimum(tail)
    end
    _check2_amp_stable = _check2_late_amplitude(2.0)
    _check2_amp_unstable = _check2_late_amplitude(0.5)
    _check2_pass = _check2_amp_stable < 1.0e-3 && _check2_amp_unstable > 1.0
    md"""
    **Stability bifurcation** (velocity-weakening, ``a`` = $(_check2_base.a) < ``b`` = $(_check2_base.b),
    ``k_c`` = $(@sprintf("%.4f", _check2_kc))): at ``k = 2k_c`` the late-time slip-rate oscillation
    amplitude is $(@sprintf("%.2e", _check2_amp_stable)) (decays to steady sliding), while at
    ``k = 0.5k_c`` it is $(@sprintf("%.2e", _check2_amp_unstable)) (a sustained limit cycle) — the
    same stiffness change flips the fault between creeping and stick-slip.
    $(_check2_pass ? "✅ PASS" : "❌ FAIL")
    """
end

# ╔═╡ 5897a5e9-3807-4fa9-bd9b-7778d5ea4c10
begin
    struct StickSlipInput
        a::Float64
        b_minus_a::Float64
        Dc::Float64
        Vpl::Float64
        sigma::Float64
        eta::Float64
        k::Float64
    end

    StickSlipInput() = StickSlipInput(0.010, 0.006, 1.0, 1.0, 1.0, 0.002, 0.0004)

    Base.get(w::StickSlipInput) = Dict{String,Any}(
        "a" => w.a,
        "b_minus_a" => w.b_minus_a,
        "Dc" => w.Dc,
        "Vpl" => w.Vpl,
        "sigma" => w.sigma,
        "eta" => w.eta,
        "k" => w.k,
    )

    """
        Base.show(io, ::MIME"text/html", w::StickSlipInput)

    Render the spring-slider cartoon, the slip-rate/friction strip chart with a phase-portrait
    inset, and the friction/stiffness control panel. The widget publishes only physical
    parameters back to Pluto; the ODE solution itself is always computed in Julia and pushed
    back in via the `ss-push` CustomEvent — JS only plays back and draws what Julia integrated.
    """
    function Base.show(io::IO, ::MIME"text/html", w::StickSlipInput)
        write(io, """
<div id="ss-widget">
  <style>
    pluto-cell:has(#ss-widget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #ss-widget{width:100%;box-sizing:border-box;font:14px sans-serif;color:#d1d5db;background:#000;padding:12px;border-radius:7px}
    #ss-widget .ss-title{box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #ss-widget .ss-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #ss-widget .ss-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #ss-widget .ss-workspace{display:grid;grid-template-columns:minmax(300px,1fr) minmax(300px,1fr);gap:9px}
    #ss-widget .ss-panel{min-width:0;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:8px}
    #ss-widget .ss-panel-title{font-size:16px;font-weight:700;color:#e5e7eb;margin:0 0 6px}
    #ss-widget .ss-caption-note{color:#ef4444;font-weight:400;font-size:13px}
    #ss-widget .ss-canvas-wrap{position:relative}
    #ss-widget canvas{display:block;width:100%;height:300px;background:#000;border:1px solid #374151;border-radius:5px;box-sizing:border-box;touch-action:none}
    #ss-widget .ss-canvas-overlay{position:absolute;bottom:10px;left:50%;transform:translateX(-50%);display:flex;gap:10px;align-items:center;background:rgba(0,0,0,.55);padding:5px 10px;border-radius:6px}
    #ss-widget .ss-canvas-overlay button{border-radius:4px;border:1px solid #9ca3af;background:#075985;color:#f3f4f6;padding:5px 10px;font-size:13px;cursor:pointer}
    #ss-widget .ss-canvas-overlay button:hover{background:#0c4a6e}
    #ss-widget #ss-phase-inset{position:absolute;top:8px;right:8px;width:132px;background:rgba(0,0,0,.75);border:1px solid #374151;border-radius:6px}
    #ss-widget #ss-phase-inset canvas{height:112px;width:132px;border:none;border-radius:0;background:transparent}
    #ss-widget .ss-phase-label{text-align:center;font-size:10px;color:#9ca3af;padding:0 3px 3px}
    #ss-widget .ss-controls{width:min(var(--ss-totalw,1100px),100%);display:grid;grid-template-columns:repeat(auto-fit,minmax(245px,1fr));gap:8px;margin-top:9px}
    #ss-widget .ss-group{min-width:0;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px}
    #ss-widget .ss-control-title{font-size:20px;font-weight:700;color:#e5e7eb;margin-bottom:8px}
    #ss-widget .ss-row{display:grid;grid-template-columns:minmax(90px,140px) minmax(70px,1fr) minmax(52px,74px);align-items:center;gap:7px;color:#9ca3af;margin:7px 0}
    #ss-widget input[type=range]{width:100%;min-width:0;accent-color:#38bdf8}
    #ss-widget .ss-value{text-align:right;color:#f3f4f6;font-variant-numeric:tabular-nums}
    #ss-widget .ss-hint{color:#9ca3af;font-size:13px;line-height:1.35;margin-top:7px}
    #ss-widget .ss-metric{display:flex;justify-content:space-between;gap:12px;margin:7px 0;color:#9ca3af}
    #ss-widget .ss-metric b{color:#f3f4f6;font-variant-numeric:tabular-nums}
    #ss-widget .ss-metric b.ss-badge-stable{color:#22c55e}
    #ss-widget .ss-metric b.ss-badge-unstable{color:#ef4444}
    #ss-widget .ss-metric b.ss-badge-near{color:#f59e0b}
    #ss-widget select{width:100%;background:#111827;color:#f3f4f6;border:1px solid #4b5563;border-radius:4px;padding:6px;font-size:14px}
    @media(max-width:700px){#ss-widget .ss-workspace{grid-template-columns:1fr}#ss-widget canvas{height:260px}}
  </style>
  <div class="ss-title">
    <div class="ss-title-desc">A spring-loaded block on a frictional surface: will it creep, or lock and slip?</div>
    <div class="ss-title-hint">Drag the friction/stiffness sliders or pick a preset &middot; press Play to watch loading turn into a slip event</div>
  </div>
  <div class="ss-workspace">
    <section class="ss-panel">
      <div class="ss-panel-title">Spring-slider model</div>
      <div class="ss-canvas-wrap">
        <canvas id="ss-scene" aria-label="Spring-slider cartoon"></canvas>
        <div class="ss-canvas-overlay">
          <button id="ss-play" type="button">Play</button>
          <button id="ss-reset" type="button">Reset</button>
        </div>
      </div>
    </section>
    <section class="ss-panel">
      <div class="ss-panel-title">Slip rate &amp; friction vs. time <span class="ss-caption-note">&mdash; red line = now</span></div>
      <div class="ss-canvas-wrap">
        <canvas id="ss-strip" aria-label="Slip rate and friction time series"></canvas>
        <div id="ss-phase-inset"><canvas id="ss-phase" aria-label="Phase portrait"></canvas><div class="ss-phase-label">phase portrait: v vs. &theta; (log&ndash;log)</div></div>
      </div>
    </section>
  </div>
  <div class="ss-controls">
    <section class="ss-group">
      <div class="ss-control-title">Rate-and-state friction</div>
      <div class="ss-row"><label for="ss-a">direct effect a</label><input id="ss-a" type="range" min="0.002" max="0.03" step="0.0005" value="$(w.a)"><span id="ss-av" class="ss-value"></span></div>
      <div class="ss-row"><label for="ss-b_minus_a">b &minus; a</label><input id="ss-b_minus_a" type="range" min="-0.02" max="0.02" step="0.0005" value="$(w.b_minus_a)"><span id="ss-b_minus_av" class="ss-value"></span></div>
      <div class="ss-row"><label for="ss-Dc">D_c</label><input id="ss-Dc" type="range" min="0.3" max="3" step="0.05" value="$(w.Dc)"><span id="ss-Dcv" class="ss-value"></span></div>
      <div class="ss-hint">b &minus; a &gt; 0 (velocity-weakening) is necessary for stick-slip at all; b &minus; a &le; 0 always creeps steadily, no matter the stiffness.</div>
    </section>
    <section class="ss-group">
      <div class="ss-control-title">Loading &amp; medium</div>
      <div class="ss-row"><label for="ss-Vpl">plate rate V_pl</label><input id="ss-Vpl" type="range" min="0.3" max="3" step="0.05" value="$(w.Vpl)"><span id="ss-Vplv" class="ss-value"></span></div>
      <div class="ss-row"><label for="ss-sigma">normal stress &sigma;</label><input id="ss-sigma" type="range" min="0.3" max="3" step="0.05" value="$(w.sigma)"><span id="ss-sigmav" class="ss-value"></span></div>
      <div class="ss-row"><label for="ss-eta">radiation damping &eta;</label><input id="ss-eta" type="range" min="0" max="0.02" step="0.0002" value="$(w.eta)"><span id="ss-etav" class="ss-value"></span></div>
    </section>
    <section class="ss-group">
      <div class="ss-control-title">Spring stiffness</div>
      <div class="ss-row"><label for="ss-k">stiffness k</label><input id="ss-k" type="range" min="0.0001" max="0.02" step="0.0001" value="$(w.k)"><span id="ss-kv" class="ss-value"></span></div>
      <div class="ss-metric"><span>k_c (critical)</span><b id="ss-kcval">&ndash;</b></div>
      <div class="ss-metric"><span>k / k_c</span><b id="ss-kratio">&ndash;</b></div>
      <div class="ss-metric"><span>regime</span><b id="ss-regime">&ndash;</b></div>
    </section>
    <section class="ss-group">
      <div class="ss-control-title">Presets</div>
      <select id="ss-preset">
        <option value="locked" selected>Locked fault (stick-slip)</option>
        <option value="creeping">Creeping fault (stable)</option>
        <option value="critical">Near-critical stiffness</option>
      </select>
      <div class="ss-hint">Same friction law, three outcomes: whether the fault ruptures or creeps depends only on stiffness relative to k_c.</div>
    </section>
  </div>
  <script>
  const script=document.currentScript, par=script.parentElement
  const scene=par.querySelector('#ss-scene'), strip=par.querySelector('#ss-strip'), phase=par.querySelector('#ss-phase')
  const sctx=scene.getContext('2d'), stx=strip.getContext('2d'), phctx=phase.getContext('2d')
  const DPR=Math.min(window.devicePixelRatio||1,2)
  const W=900,H=300,PW=132,PH=112
  const PLAYBACK_SECONDS=12
  let state={a:$(w.a),b_minus_a:$(w.b_minus_a),Dc:$(w.Dc),Vpl:$(w.Vpl),sigma:$(w.sigma),eta:$(w.eta),k:$(w.k)}
  let data=null, playIdx=0, playing=false, playStart=0, commitInFlight=false, pendingCommit=false
  const sliderIds=['a','b_minus_a','Dc','Vpl','sigma','eta','k']
  const inputs={}, values={}
  sliderIds.forEach(key=>{inputs[key]=par.querySelector('#ss-'+key);values[key]=par.querySelector('#ss-'+key+'v')})
  const playBtn=par.querySelector('#ss-play'), resetBtn=par.querySelector('#ss-reset')
  const kcValEl=par.querySelector('#ss-kcval'), kratioEl=par.querySelector('#ss-kratio'), regimeEl=par.querySelector('#ss-regime')
  const PRESETS={
    locked:{a:0.010,b_minus_a:0.006,Dc:1.0,Vpl:1.0,sigma:1.0,eta:0.002,k:0.0004},
    creeping:{a:0.016,b_minus_a:-0.006,Dc:1.0,Vpl:1.0,sigma:1.0,eta:0.002,k:0.0004},
    critical:{a:0.010,b_minus_a:0.006,Dc:1.0,Vpl:1.0,sigma:1.0,eta:0.002,k:0.0032},
  }
  function kcJS(s){return (s.sigma*s.b_minus_a - s.eta*s.Vpl)/s.Dc}
  function regimeOf(s){
    if(s.b_minus_a<=0) return 'stable (velocity-strengthening)'
    const kc=kcJS(s)
    if(kc<=0) return 'stable (velocity-strengthening)'
    const ratio=s.k/kc
    if(ratio>1.3) return 'stable (creeping)'
    if(ratio<0.7) return 'unstable (stick-slip)'
    return 'near-critical'
  }
  function regimeClass(r){return r.startsWith('stable')?'ss-badge-stable':r.startsWith('unstable')?'ss-badge-unstable':'ss-badge-near'}
  function fmt(key,val){
    if(key==='a'||key==='Dc') return val.toFixed(3)
    if(key==='b_minus_a') return (val>=0?'+':'')+val.toFixed(3)
    if(key==='eta'||key==='k') return val.toFixed(4)
    return val.toFixed(2)
  }
  function updateReadouts(){
    const kc=kcJS(state)
    kcValEl.textContent=kc.toFixed(5)
    kratioEl.textContent=(state.b_minus_a<=0||kc<=0)?'n/a':(state.k/kc).toFixed(2)
    const r=regimeOf(state)
    regimeEl.textContent=r
    regimeEl.className=regimeClass(r)
  }
  function updateControls(){sliderIds.forEach(key=>{inputs[key].value=state[key];values[key].textContent=fmt(key,state[key])});updateReadouts()}
  function emit(){commitInFlight=true;pendingCommit=false;par.value={a:state.a,b_minus_a:state.b_minus_a,Dc:state.Dc,Vpl:state.Vpl,sigma:state.sigma,eta:state.eta,k:state.k};par.dispatchEvent(new CustomEvent('input'))}
  function throttledCommit(){if(commitInFlight){pendingCommit=true;return}emit()}
  function stopPlayback(){playing=false;playBtn.textContent='Play';playIdx=0}
  par.addEventListener('input',e=>{
    const id=e.target.id
    if(id && id.indexOf('ss-')===0){
      const key=id.slice(3)
      if(sliderIds.includes(key)){
        state[key]=Number(e.target.value)
        updateControls()
        stopPlayback()
        draw()
        throttledCommit()
      }
    }
  })
  par.addEventListener('change',e=>{
    if(e.target.id==='ss-preset'){
      const p=PRESETS[e.target.value]
      if(p){Object.assign(state,p);updateControls();stopPlayback();draw();throttledCommit()}
    }
  })
  par.addEventListener('ss-push',e=>{
    const d=e.detail
    data={t:d.t,v:d.v,theta:d.theta,mu:d.mu,stretch:d.stretch,kc:d.kc}
    data.vmax=Math.max(...data.v)
    data.smin=Math.min(...data.stretch);data.smax=Math.max(...data.stretch)
    playIdx=0
    commitInFlight=false
    if(pendingCommit){pendingCommit=false;emit()}
    draw()
  })
  playBtn.addEventListener('click',()=>{
    if(!data) return
    playing=!playing
    if(playing){
      if(playIdx>=data.t.length-1) playIdx=0
      playStart=performance.now()-(playIdx/(data.t.length-1))*PLAYBACK_SECONDS*1000
      playBtn.textContent='Pause'
      requestAnimationFrame(tick)
    } else {
      playBtn.textContent='Play'
    }
  })
  resetBtn.addEventListener('click',()=>{stopPlayback();draw()})
  function tick(now){
    if(!playing) return
    const elapsed=(now-playStart)/1000
    let phase=elapsed/PLAYBACK_SECONDS
    if(phase>=1){phase=1;playing=false;playBtn.textContent='Play'}
    playIdx=Math.min(data.t.length-1,Math.floor(phase*(data.t.length-1)))
    draw()
    if(playing) requestAnimationFrame(tick)
  }
  function hidpi(canvas,ctx,w,h){const rect=canvas.getBoundingClientRect();canvas.width=Math.max(1,Math.round(rect.width*DPR));canvas.height=Math.max(1,Math.round(rect.height*DPR));const scale=Math.min(canvas.width/w,canvas.height/h),dx=(canvas.width-w*scale)/2,dy=(canvas.height-h*scale)/2;ctx.setTransform(scale,0,0,scale,dx,dy)}
  function line(ctx,x1,y1,x2,y2,color,width){ctx.save();ctx.strokeStyle=color;ctx.lineWidth=width;ctx.beginPath();ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);ctx.stroke();ctx.restore()}
  function spring(ctx,x1,x2,y){ctx.beginPath();const n=60;for(let i=0;i<=n;i++){const f=i/n,x=x1+(x2-x1)*f,yy=y+(i>0&&i<n?Math.sin(f*Math.PI*10)*9:0);if(i===0)ctx.moveTo(x,yy);else ctx.lineTo(x,yy)}ctx.stroke()}
  function drawEmpty(ctx,w,h){ctx.clearRect(0,0,w,h);ctx.fillStyle='#000';ctx.fillRect(0,0,w,h);ctx.fillStyle='#9ca3af';ctx.font='14px sans-serif';ctx.textAlign='center';ctx.fillText('computing…',w/2,h/2)}
  function drawScene(){
    if(!data) return drawEmpty(sctx,W,H)
    sctx.clearRect(0,0,W,H);sctx.fillStyle='#000';sctx.fillRect(0,0,W,H)
    const i=playIdx, stretch=data.stretch[i], v=data.v[i], mu=data.mu[i]
    const wallX=90, groundY=210, baseX=180, spanPx=560
    const frac=(stretch-data.smin)/Math.max(1e-9,(data.smax-data.smin))
    const blockX=baseX+frac*spanPx
    sctx.fillStyle='#4b5563';sctx.fillRect(wallX-14,60,14,170)
    for(let yy=64;yy<220;yy+=14){line(sctx,wallX-14,yy,wallX,yy+10,'#374151',2)}
    line(sctx,60,groundY,860,groundY,'#6b7280',2)
    for(let x=64;x<860;x+=16){line(sctx,x,groundY,x-8,groundY+10,'#374151',1.5)}
    sctx.fillStyle='#9ca3af';sctx.font='13px sans-serif';sctx.textAlign='left'
    sctx.fillText('plate loading at V_pl →',560,52)
    sctx.strokeStyle='#e5e7eb';sctx.lineWidth=2.4
    spring(sctx,wallX,blockX-34,135)
    const glow=Math.min(1,v/(data.vmax+1e-9))
    sctx.fillStyle='rgba(56,189,248,'+(0.5+0.5*glow)+')'
    sctx.fillRect(blockX-34,135-24,68,48)
    sctx.strokeStyle='#38bdf8';sctx.lineWidth=2;sctx.strokeRect(blockX-34,135-24,68,48)
    sctx.fillStyle='#0b0b0b';sctx.font='700 15px sans-serif';sctx.textAlign='center'
    sctx.fillText('block',blockX,135+5)
    sctx.fillStyle='#f3f4f6';sctx.font='14px sans-serif';sctx.textAlign='left'
    sctx.fillText('t = '+data.t[i].toFixed(2),40,34)
    sctx.fillText('v = '+v.toExponential(2),40,54)
    sctx.fillText('μ = '+mu.toFixed(4),40,74)
  }
  function drawStrip(){
    if(!data) return drawEmpty(stx,W,H)
    stx.clearRect(0,0,W,H);stx.fillStyle='#000';stx.fillRect(0,0,W,H)
    const left=66,right=860,topV=24,botV=170,topMu=192,botMu=280
    const tmax=data.t[data.t.length-1]
    const logv=data.v.map(x=>Math.log10(Math.max(x,1e-12)))
    const lvmin=Math.min(...logv),lvmax=Math.max(...logv)
    function xOf(t){return left+(right-left)*t/tmax}
    function yOfV(lv){return botV-(lv-lvmin)/Math.max(1e-9,(lvmax-lvmin))*(botV-topV)}
    const mumin=Math.min(...data.mu),mumax=Math.max(...data.mu)
    function yOfMu(m){return botMu-(m-mumin)/Math.max(1e-9,(mumax-mumin))*(botMu-topMu)}
    stx.fillStyle='#9ca3af';stx.font='12px sans-serif';stx.textAlign='right'
    for(let p=Math.ceil(lvmin);p<=Math.floor(lvmax);p++){const y=yOfV(p);line(stx,left,y,right,y,'#1f2937',1);stx.fillText('1e'+p,left-8,y+4)}
    stx.textAlign='left';stx.fillStyle='#e5e7eb';stx.font='13px sans-serif';stx.fillText('slip rate v (log scale)',left,14)
    stx.strokeStyle='#38bdf8';stx.lineWidth=2;stx.beginPath()
    for(let i=0;i<data.t.length;i++){const x=xOf(data.t[i]),y=yOfV(logv[i]);if(i===0)stx.moveTo(x,y);else stx.lineTo(x,y)}
    stx.stroke()
    stx.fillStyle='#e5e7eb';stx.font='13px sans-serif';stx.textAlign='left';stx.fillText('friction coefficient μ',left,topMu-6)
    line(stx,left,botMu,right,botMu,'#374151',1)
    stx.strokeStyle='#a78bfa';stx.lineWidth=2;stx.beginPath()
    for(let i=0;i<data.t.length;i++){const x=xOf(data.t[i]),y=yOfMu(data.mu[i]);if(i===0)stx.moveTo(x,y);else stx.lineTo(x,y)}
    stx.stroke()
    const nowX=xOf(data.t[playIdx])
    line(stx,nowX,topV,nowX,botMu,'#ef4444',2)
    stx.fillStyle='#9ca3af';stx.textAlign='center';stx.fillText('time',(left+right)/2,H-6)
  }
  function drawPhase(){
    if(!data){phctx.clearRect(0,0,PW,PH);return}
    phctx.clearRect(0,0,PW,PH)
    const pad=14,w=PW-2*pad,h=PH-2*pad
    const logv=data.v.map(x=>Math.log10(Math.max(x,1e-12)))
    const logth=data.theta.map(x=>Math.log10(Math.max(x,1e-12)))
    const lvmin=Math.min(...logv),lvmax=Math.max(...logv),ltmin=Math.min(...logth),ltmax=Math.max(...logth)
    function xOf(lt){return pad+(lt-ltmin)/Math.max(1e-9,(ltmax-ltmin))*w}
    function yOf(lv){return (PH-pad)-(lv-lvmin)/Math.max(1e-9,(lvmax-lvmin))*h}
    phctx.strokeStyle='rgba(56,189,248,0.55)';phctx.lineWidth=1.2;phctx.beginPath()
    for(let i=0;i<data.t.length;i++){const x=xOf(logth[i]),y=yOf(logv[i]);if(i===0)phctx.moveTo(x,y);else phctx.lineTo(x,y)}
    phctx.stroke()
    const x0=xOf(logth[playIdx]),y0=yOf(logv[playIdx])
    phctx.fillStyle='#ef4444';phctx.beginPath();phctx.arc(x0,y0,3,0,2*Math.PI);phctx.fill()
  }
  function draw(){drawScene();drawStrip();drawPhase()}
  const observer=new ResizeObserver(()=>{hidpi(scene,sctx,W,H);hidpi(strip,stx,W,H);hidpi(phase,phctx,PW,PH);draw()})
  observer.observe(scene);observer.observe(strip)
  par.style.setProperty('--ss-totalw',Math.round(Math.min(window.innerWidth*.8,par.clientWidth||1100))+'px')
  updateControls()
  hidpi(scene,sctx,W,H);hidpi(strip,stx,W,H);hidpi(phase,phctx,PW,PH)
  draw()
  emit()
  </script>
</div>
        """)
    end

    const _ss_ready = true
end

# ╔═╡ 45405706-2ca9-466b-b67a-8b8f387c3423
begin
    _ss_ready
    WideCell(@bind ss_state StickSlipInput(); max_width=1500)
end

# ╔═╡ e22c0ca4-43c7-4c3b-a12f-27794ac48616
"""
    bounded(raw, key, default, lo, hi)

Read one finite control value from the widget's bound state and clip it to the physical
range used by the demonstration, so a malformed or missing browser value cannot create a
non-finite or sign-flipped parameter downstream.
"""
function bounded(raw, key, default, lo, hi)
    candidate = try
        raw isa AbstractDict ? Float64(get(raw, key, default)) : default
    catch
        default
    end
    isfinite(candidate) ? clamp(candidate, lo, hi) : default
end

# ╔═╡ e446deb7-e077-435c-a178-b5a0c8b35190
begin
    _ss_a = bounded(ss_state, "a", 0.010, 0.001, 0.05)
    _ss_bma = bounded(ss_state, "b_minus_a", 0.006, -0.03, 0.03)
    _ss_Dc = bounded(ss_state, "Dc", 1.0, 0.1, 5.0)
    _ss_Vpl = bounded(ss_state, "Vpl", 1.0, 0.1, 5.0)
    _ss_sigma = bounded(ss_state, "sigma", 1.0, 0.1, 5.0)
    _ss_eta = bounded(ss_state, "eta", 0.002, 0.0, 0.05)
    _ss_k = bounded(ss_state, "k", 0.0004, 1.0e-5, 0.05)
    _ss_b = _ss_a + _ss_bma
end

# ╔═╡ c856f825-e50a-4f7c-8b2f-d99f2a5ab8b3
_ss_solution = integrate_spring_slider(a=_ss_a, b=_ss_b, Dc=_ss_Dc, Vpl=_ss_Vpl, sigma=_ss_sigma, eta=_ss_eta, k=_ss_k)

# ╔═╡ 88175629-ae25-4826-94c5-ed588807ddb6
StickSlipPush(_ss_solution.t, _ss_solution.v, _ss_solution.theta, _ss_solution.mu, _ss_solution.stretch, _ss_solution.kc)

# ╔═╡ d76775c3-8e3b-443a-ac15-1a1cd8a9173b
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
PlutoUI = "~0.7.83"
"""

# ╔═╡ 7babf4d8-c4df-446e-af75-1f819db1b454
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "e4dcbcd30a3d5ffa1505607c8b44298e5599088d"

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.0"
manifest_format = "2.1"
project_hash = "e4dcbcd30a3d5ffa1505607c8b44298e5599088d"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

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
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

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
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

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
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.18.0+1"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.103+0"

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

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
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

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
registries = "General"
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
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"

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
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
registries = "General"
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

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╠═32526135-0ec9-4865-bc05-fa81450ff3a3
# ╠═7f48d81c-5664-459a-bf10-c265462d57a0
# ╟─1103da07-e3f7-408b-937e-5ee0a1563ac7
# ╟─8417cf97-4555-44c6-9f44-a6a32799bde8
# ╟─45405706-2ca9-466b-b67a-8b8f387c3423
# ╠═e22c0ca4-43c7-4c3b-a12f-27794ac48616
# ╠═e446deb7-e077-435c-a178-b5a0c8b35190
# ╟─d723bd18-4389-445e-b747-fab5684f85a7
# ╟─d54456fb-3b04-4a9a-ba51-3df430f0837f
# ╟─20022d79-0eb9-4e08-9a07-bfde9a9a4db5
# ╟─874a667c-c6cd-4797-818e-01b8dfc57e1b
# ╠═69a1dd35-07d4-4d7a-8fcf-6dfc00d2e092
# ╠═4133723b-9d54-403e-84cb-92f268523f99
# ╠═9b71d1c6-70b5-40fa-9e43-c3c02b4082e0
# ╠═48043f9d-d937-443d-9132-8c8c921a5d86
# ╠═cabca0a8-ffb6-40aa-ab08-d507b0e0e063
# ╠═0e435cbe-662d-4b75-9530-0d0fa898399a
# ╠═c856f825-e50a-4f7c-8b2f-d99f2a5ab8b3
# ╠═6df71f7e-1f23-4edb-8878-e07aba98a32d
# ╠═a3fd5f1a-90a3-44cd-b1d4-cdbd9b22b6c5
# ╠═88175629-ae25-4826-94c5-ed588807ddb6
# ╟─331f1292-1115-444f-931f-ac46bf88943f
# ╠═e63b931a-94d9-42f1-a6c6-e8143faa0911
# ╠═b27449d6-07fd-411b-af09-88f0524939e0
# ╠═5897a5e9-3807-4fa9-bd9b-7778d5ea4c10
# ╟─d76775c3-8e3b-443a-ac15-1a1cd8a9173b
# ╟─7babf4d8-c4df-446e-af75-1f819db1b454
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
