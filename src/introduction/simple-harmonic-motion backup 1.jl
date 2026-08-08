### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Simple Harmonic Motion — Drag, Release, Explain"
#> tags = ["introduction", "oscillations", "waves"]
#> layout = "layout.jlhtml"
#> description = "A draggable mass–spring–damper that connects initial conditions to exact oscillator solutions."

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

# ╔═╡ 27ab0195-38a6-4235-a782-8a5b8a7bac01
begin
    using PlutoUI
    using Printf
end

# ╔═╡ 548b0e6a-7284-49e1-88f7-e5ef9f8a2de5
PlutoUI.TableOfContents(include_definitions=true)

# ╔═╡ 9a3af4df-76fe-4381-87b7-56a3b43925cc
md"""
# Simple harmonic motion: make the initial conditions visible

An oscillator is the smallest possible model of a vibrating Earth. Pull a mass away from
equilibrium, give it an initial velocity, and its later motion is no longer a guess: the
equation of motion determines every point on the time history.

This notebook keeps that connection in view. Your hand sets the initial condition; the
animation, graph, and exact solution all update together.
"""

# ╔═╡ c1d09d0d-2c25-4db5-9298-e4a6e9fa550a
md"""
## Pull, release, and read the answer

Drag the blue block horizontally to set the initial displacement ``x_0``. Use the initial
velocity control when you want to launch it rather than release it from rest. The red line
on the time-history plot is the instant currently shown in the mechanical model.

!!! tip "A useful first experiment"
	Start with **undamped**, drag the block, and press Play. Then change only the mass or
	stiffness. The amplitude stays fixed, but the clock runs at a different rate.
"""

# ╔═╡ 6e3ba5f5-70ec-45bf-9a56-e0a5990d453d
begin
    _shm_ready
    WideCell(@bind shm_state SHMInput(); max_width=1500)
end

# ╔═╡ 6fba3aa2-6192-4f45-8931-2e715e213786
begin
    m = bounded_state_value(shm_state, "m", 1.0, 0.25, 4.0)
    k = bounded_state_value(shm_state, "k", 9.0, 1.0, 36.0)
    ζ = bounded_state_value(shm_state, "zeta", 0.15, 0.0, 1.6)
    x₀ = bounded_state_value(shm_state, "x0", 0.70, -1.2, 1.2)
    v₀ = bounded_state_value(shm_state, "v0", 0.0, -4.0, 4.0)
    playback_speed = bounded_state_value(shm_state, "speed", 0.6, 0.2, 2.0)

    ωₙ = sqrt(k / m)
    c = 2ζ * sqrt(k * m)
    ωd = ζ < 1 ? ωₙ * sqrt(1 - ζ^2) : 0.0
    damping_regime = ζ < 1e-6 ? "undamped" : ζ < 1 - 1e-6 ? "underdamped" : abs(ζ - 1) ≤ 1e-6 ? "critically damped" : "overdamped"
    Tₙ = 2π / ωₙ
    t_end = ζ < 1e-6 ? 5Tₙ : min(20.0, max(2Tₙ, 5 / (max(ζ, 0.08) * ωₙ)))
    t_history = collect(range(0.0, t_end; length=1001))
    x_history = exact_displacement.(t_history, Ref(m), Ref(k), Ref(ζ), Ref(x₀), Ref(v₀))
end

# ╔═╡ 94ab3b93-4831-48d3-8ae8-b3f82e8738a2
begin
    solution_latex, solution_note = current_solution_latex(m, k, ζ, x₀, v₀)
    HTML("""
    <section style="margin:18px 0;padding:16px 18px;background:#050505;border:1px solid #2f3744;border-radius:7px;color:#d1d5db">
      <div style="font:700 20px sans-serif;color:#f3f4f6">Your exact $(damping_regime) solution</div>
      <div style="margin:6px 0 10px;color:#9ca3af;font:14px sans-serif">$(solution_note)</div>
      <div class="tex" style="overflow-x:auto;font-size:1.08rem">\\[ $(solution_latex) \\]</div>
      <div style="display:flex;flex-wrap:wrap;gap:8px 22px;font:14px sans-serif">
        <span><b style="color:#38bdf8">ωₙ</b> = $(@sprintf("%.3f", ωₙ)) rad s⁻¹</span>
        <span><b style="color:#a78bfa">ζ</b> = $(@sprintf("%.3f", ζ))</span>
        <span><b style="color:#f3f4f6">c</b> = $(@sprintf("%.3f", c)) N s m⁻¹</span>
        <span><b style="color:#f3f4f6">Tₙ</b> = $(@sprintf("%.3f", Tₙ)) s</span>
      </div>
    </section>
    """)
end

# ╔═╡ b4f616c5-925c-40ea-9a74-1b1f68848aed
md"""
## From a restoring force to an equation

The spring pulls toward equilibrium, so its force is ``F_s=-kx``. The dashpot opposes
velocity, ``F_d=-c\dot{x}``. Newton's law gives

```math
m\ddot{x} + c\dot{x} + kx = 0.
```

The two initial conditions are indispensable: ``x(0)=x_0`` says where the block starts;
``\dot{x}(0)=v_0`` says how it is launched. The widget's blue curve is the solution to
this equation with exactly those two numbers.
"""

# ╔═╡ 8ed80c0f-8177-4229-ad4e-a1a9f075718f
md"""
## Damping changes the *kind* of solution

The damping ratio ``\zeta=c/(2\sqrt{km})`` compares the actual damping to the critical
value. It is more informative than ``c`` alone because it stays meaningful when mass or
stiffness changes.

| Setting | What you see | Mathematical signature |
|:--|:--|:--|
| ``\zeta=0`` | perpetual oscillation | sine and cosine |
| ``0<\zeta<1`` | oscillation inside a shrinking envelope | ``e^{-\zeta\omega_n t}`` times a sinusoid |
| ``\zeta=1`` | fastest return without overshoot | exponential times a line |
| ``\zeta>1`` | slow, non-oscillatory return | sum of two exponentials |

!!! warning "Critical damping is not the most damping"
	It is the boundary between oscillatory and non-oscillatory motion. Adding more damping
	does not make the block settle faster; it makes the long, slow exponential tail stronger.
"""

# ╔═╡ 174420d3-9c28-4d3f-93b2-6b9a267a222a
md"""
## Why seismologists care

A real Earth has enormously many coupled degrees of freedom, but each normal mode has the
same temporal form as this one oscillator. An earthquake supplies initial displacement and
velocity; Earth structure sets the mode frequency; attenuation supplies weak damping.

For light damping, the seismic quality factor is approximately ``Q \approx 1/(2\zeta)``.
Large ``Q`` means many visible cycles before the envelope fades. The coupled-oscillations
notebook is the natural next step: replace this one mass with many masses and the single
frequency becomes a family of normal modes.
"""

# ╔═╡ 405af09e-aa7b-46f5-9877-e9dfc4e1427c
md"""
## Appendix

### Exact free-vibration solutions

The widget evaluates the same closed-form expressions shown above in the browser so the
block and playhead remain smooth. Julia evaluates them independently here for the bound
state, which is what makes the displayed numerical solution and the verification below
auditable.
"""

# ╔═╡ 73c5bd9b-9a8f-41a9-b2a4-d076a15bbbca
"""
    bounded_state_value(raw, key, default, lo, hi)

Read one finite control value from the widget state and clip it to the physical range used
by the demonstration. A malformed browser value therefore cannot create a negative mass,
negative stiffness, or a singular formula.
"""
function bounded_state_value(raw, key, default, lo, hi)
    candidate = try
        raw isa AbstractDict ? Float64(get(raw, key, default)) : default
    catch
        default
    end
    isfinite(candidate) ? clamp(candidate, lo, hi) : default
end

# ╔═╡ bc14b27d-af7e-499f-b169-cb4814810a5d
"""
    exact_displacement(t, m, k, ζ, x₀, v₀)

Return the exact displacement of a free mass–spring–damper system at time `t`. The branch
is chosen by the damping ratio, while `x₀` and `v₀` enforce the two initial conditions.
"""
function exact_displacement(t, m, k, ζ, x₀, v₀)
    ωₙ = sqrt(k / m)
    if ζ < 1e-6
        return x₀ * cos(ωₙ * t) + (v₀ / ωₙ) * sin(ωₙ * t)
    elseif ζ < 1 - 1e-6
        ωd = ωₙ * sqrt(1 - ζ^2)
        return exp(-ζ * ωₙ * t) * (x₀ * cos(ωd * t) + (v₀ + ζ * ωₙ * x₀) / ωd * sin(ωd * t))
    elseif ζ ≤ 1 + 1e-6
        return exp(-ωₙ * t) * (x₀ + (v₀ + ωₙ * x₀) * t)
    else
        r₁ = -ωₙ * (ζ - sqrt(ζ^2 - 1))
        r₂ = -ωₙ * (ζ + sqrt(ζ^2 - 1))
        A = (v₀ - r₂ * x₀) / (r₁ - r₂)
        B = x₀ - A
        return A * exp(r₁ * t) + B * exp(r₂ * t)
    end
end

# ╔═╡ 90b3d802-0993-4df2-8b48-3c3c462d5c8f
"""
    current_solution_latex(m, k, ζ, x₀, v₀)

Build the numerical closed-form solution displayed above the derivation. Coefficients are
rounded only for presentation; [`exact_displacement`](@ref) retains full precision.
"""
function current_solution_latex(m, k, ζ, x₀, v₀)
    ωₙ = sqrt(k / m)
    if ζ < 1e-6
        return @sprintf("x(t) = %.3f\\cos(%.3f\\,t) + %.3f\\sin(%.3f\\,t)", x₀, ωₙ, v₀ / ωₙ, ωₙ), "The amplitude does not decay."
    elseif ζ < 1 - 1e-6
        ωd = ωₙ * sqrt(1 - ζ^2)
        B = (v₀ + ζ * ωₙ * x₀) / ωd
        return @sprintf("x(t) = e^{-%.3f t}\\left[%.3f\\cos(%.3f\\,t) + %.3f\\sin(%.3f\\,t)\\right]", ζ * ωₙ, x₀, ωd, B, ωd), "The exponential envelope and oscillation have different clocks."
    elseif ζ ≤ 1 + 1e-6
        B = v₀ + ωₙ * x₀
        return @sprintf("x(t) = e^{-%.3f t}\\left[%.3f + %.3f\\,t\\right]", ωₙ, x₀, B), "There is no sinusoid, and no overshoot."
    else
        r₁ = -ωₙ * (ζ - sqrt(ζ^2 - 1))
        r₂ = -ωₙ * (ζ + sqrt(ζ^2 - 1))
        A = (v₀ - r₂ * x₀) / (r₁ - r₂)
        B = x₀ - A
        return @sprintf("x(t) = %.3f\\,e^{%.3f t} + %.3f\\,e^{%.3f t}", A, r₁, B, r₂), "Two decaying exponentials return the mass without a cycle."
    end
end

# ╔═╡ f6d03a7b-3c89-4799-bcd3-8c799711ea1e
begin
    h = min(1e-4 / ωₙ, t_end / 5000)
    check_times = range(3h, t_end - 3h; length=61)
    residuals = [
        m * (exact_displacement(t + h, m, k, ζ, x₀, v₀) - 2 * exact_displacement(t, m, k, ζ, x₀, v₀) + exact_displacement(t - h, m, k, ζ, x₀, v₀)) / h^2 +
        c * (exact_displacement(t + h, m, k, ζ, x₀, v₀) - exact_displacement(t - h, m, k, ζ, x₀, v₀)) / (2h) +
        k * exact_displacement(t, m, k, ζ, x₀, v₀)
        for t in check_times
    ]
    md"""
    ### Verifying the equation

    A centered finite-difference check of ``m\ddot{x}+c\dot{x}+kx`` over the displayed
    interval has maximum residual **$(@sprintf("%.2e", maximum(abs.(residuals))))** N.
    The small nonzero number is numerical differentiation error, not an additional force.
    """
end

# ╔═╡ 6694ed08-e429-41b8-b8ac-a2f9b417954b
md"""
### The interactive widget

The canvas keeps a fixed logical coordinate system and uses a device-pixel-ratio transform
for a sharp display. Dragging therefore remains physically calibrated while the graphic
stays clear on high-resolution laptop and projector screens.
"""

# ╔═╡ de2db126-21d0-4f72-9d26-79ba19d1b65b
begin
    struct SHMInput
        m::Float64
        k::Float64
        zeta::Float64
        x0::Float64
        v0::Float64
        speed::Float64
    end

    SHMInput() = SHMInput(1.0, 9.0, 0.15, 0.70, 0.0, 0.60)

    Base.get(w::SHMInput) = Dict{String,Any}(
        "m" => w.m,
        "k" => w.k,
        "zeta" => w.zeta,
        "x0" => w.x0,
        "v0" => w.v0,
        "speed" => w.speed,
    )

    """
        Base.show(io, ::MIME"text/html", w::SHMInput)

    Render the draggable mass–spring–damper and its browser-smooth exact time history.
    The widget publishes only physical state, never pixel coordinates, back to Pluto.
    """
    function Base.show(io::IO, ::MIME"text/html", w::SHMInput)
        write(io, """
<div id="shm-widget">
  <style>
    pluto-cell:has(#shm-widget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #shm-widget{width:100%;box-sizing:border-box;font:14px sans-serif;color:#d1d5db;background:#000;padding:12px;border-radius:7px}
    #shm-widget .shm-title{box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #shm-widget .shm-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #shm-widget .shm-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #shm-widget .shm-workspace{display:grid;grid-template-columns:minmax(300px,1fr) minmax(300px,1fr);gap:9px}
    #shm-widget .shm-panel{min-width:0;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:8px}
    #shm-widget .shm-panel-title{font-size:16px;font-weight:700;color:#e5e7eb;margin:0 0 6px}
    #shm-widget canvas{display:block;width:100%;height:270px;background:#000;border:1px solid #374151;border-radius:5px;box-sizing:border-box;touch-action:none}
    #shm-widget .shm-controls{width:min(var(--shm-totalw,960px),100%);display:grid;grid-template-columns:repeat(auto-fit,minmax(245px,1fr));gap:8px;margin-top:9px}
    #shm-widget .shm-group{min-width:0;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px}
    #shm-widget .shm-control-title{font-size:20px;font-weight:700;color:#e5e7eb;margin-bottom:8px}
    #shm-widget .shm-row{display:grid;grid-template-columns:minmax(73px,128px) minmax(70px,1fr) minmax(44px,68px);align-items:center;gap:7px;color:#9ca3af;margin:7px 0}
    #shm-widget input[type=range]{width:100%;min-width:0;accent-color:#38bdf8}
    #shm-widget .shm-value{text-align:right;color:#f3f4f6;font-variant-numeric:tabular-nums}
    #shm-widget .shm-buttons{display:flex;flex-wrap:wrap;gap:6px;margin-top:8px}
    #shm-widget button{border-radius:4px;border:1px solid #9ca3af;background:#075985;color:#f3f4f6;padding:6px 10px;font-size:14px;cursor:pointer}
    #shm-widget button:hover{background:#0c4a6e}
    #shm-widget button.shm-preset{background:#171717;border-color:#4b5563;color:#d1d5db}
    #shm-widget button.shm-preset:hover{background:#27272a}
    #shm-widget .shm-hint{color:#9ca3af;font-size:13px;line-height:1.35;margin-top:7px}
    #shm-widget .shm-metric{display:flex;justify-content:space-between;gap:12px;margin:7px 0;color:#9ca3af}
    #shm-widget .shm-metric b{color:#f3f4f6;font-variant-numeric:tabular-nums}
    @media(max-width:700px){#shm-widget .shm-workspace{grid-template-columns:1fr}#shm-widget canvas{height:245px}}
    @media(max-width:430px){#shm-widget .shm-row{grid-template-columns:78px minmax(60px,1fr) 50px}}
  </style>
  <div class="shm-title">
    <div class="shm-title-desc">One initial push determines an entire exact time history.</div>
    <div class="shm-title-hint">Drag the block to set displacement &middot; set launch velocity &middot; compare damping regimes</div>
  </div>
  <div class="shm-workspace">
    <section class="shm-panel"><div class="shm-panel-title">Mechanical model</div><canvas id="shm-scene" aria-label="Draggable mass spring damper animation"></canvas></section>
    <section class="shm-panel"><div class="shm-panel-title">Displacement history <span style="color:#ef4444;font-weight:400">— red line = now</span></div><canvas id="shm-plot" aria-label="Displacement versus time plot"></canvas></section>
  </div>
  <div class="shm-controls">
    <section class="shm-group">
      <div class="shm-control-title">Physical parameters</div>
      <div class="shm-row"><label for="shm-m">mass m</label><input id="shm-m" type="range" min="0.25" max="4" step="0.05" value="$(w.m)"><span id="shm-mv" class="shm-value"></span></div>
      <div class="shm-row"><label for="shm-k">stiffness k</label><input id="shm-k" type="range" min="1" max="36" step="0.5" value="$(w.k)"><span id="shm-kv" class="shm-value"></span></div>
      <div class="shm-row"><label for="shm-z">damping ζ</label><input id="shm-z" type="range" min="0" max="1.6" step="0.01" value="$(w.zeta)"><span id="shm-zv" class="shm-value"></span></div>
      <div class="shm-hint">ζ is the damping ratio. The actual dashpot coefficient is derived below.</div>
    </section>
    <section class="shm-group">
      <div class="shm-control-title">Initial condition &amp; playback</div>
      <div class="shm-row"><label for="shm-v">launch v₀</label><input id="shm-v" type="range" min="-4" max="4" step="0.05" value="$(w.v0)"><span id="shm-vv" class="shm-value"></span></div>
      <div class="shm-row"><label for="shm-speed">play speed</label><input id="shm-speed" type="range" min="0.2" max="2" step="0.1" value="$(w.speed)"><span id="shm-speedv" class="shm-value"></span></div>
      <div class="shm-buttons"><button id="shm-play" type="button">Play</button><button id="shm-restart" type="button">Restart</button><button id="shm-home" type="button">Reset IC</button></div>
      <div class="shm-metric"><span>current time</span><b id="shm-time">0.00 s</b></div>
      <div class="shm-metric"><span>regime</span><b id="shm-regime">underdamped</b></div>
    </section>
    <section class="shm-group">
      <div class="shm-control-title">Damping presets</div>
      <div class="shm-buttons"><button class="shm-preset" data-z="0">Undamped</button><button class="shm-preset" data-z="0.15">Underdamped</button><button class="shm-preset" data-z="1">Critical</button><button class="shm-preset" data-z="1.4">Overdamped</button></div>
      <div class="shm-hint">The graph always shows the full response. Critical damping returns fastest without crossing equilibrium.</div>
    </section>
  </div>
  <script>
  const script=document.currentScript, par=script.parentElement
  const scene=par.querySelector('#shm-scene'), plot=par.querySelector('#shm-plot')
  const sctx=scene.getContext('2d'), pctx=plot.getContext('2d')
  const inputs={m:par.querySelector('#shm-m'),k:par.querySelector('#shm-k'),zeta:par.querySelector('#shm-z'),v0:par.querySelector('#shm-v'),speed:par.querySelector('#shm-speed')}
  const values={m:par.querySelector('#shm-mv'),k:par.querySelector('#shm-kv'),zeta:par.querySelector('#shm-zv'),v0:par.querySelector('#shm-vv'),speed:par.querySelector('#shm-speedv')}
  const playButton=par.querySelector('#shm-play'), timeReadout=par.querySelector('#shm-time'), regimeReadout=par.querySelector('#shm-regime')
  const W=900,H=270,DPR=Math.min(window.devicePixelRatio||1,2)
  let state={m:$(w.m),k:$(w.k),zeta:$(w.zeta),x0:$(w.x0),v0:$(w.v0),speed:$(w.speed),t:0,playing:false,last:0,dragging:false}
  function clamp(x,a,b){return Math.max(a,Math.min(b,x))}
  function naturalFrequency(){return Math.sqrt(state.k/state.m)}
  function regime(){return state.zeta<1e-6?'undamped':state.zeta<1-1e-6?'underdamped':Math.abs(state.zeta-1)<=1e-6?'critical damping':'overdamped'}
  function displacement(t){const wn=naturalFrequency(),z=state.zeta,x0=state.x0,v0=state.v0;if(z<1e-6)return x0*Math.cos(wn*t)+(v0/wn)*Math.sin(wn*t);if(z<1-1e-6){const wd=wn*Math.sqrt(1-z*z);return Math.exp(-z*wn*t)*(x0*Math.cos(wd*t)+(v0+z*wn*x0)/wd*Math.sin(wd*t))}if(z<=1+1e-6)return Math.exp(-wn*t)*(x0+(v0+wn*x0)*t);const r1=-wn*(z-Math.sqrt(z*z-1)),r2=-wn*(z+Math.sqrt(z*z-1)),A=(v0-r2*x0)/(r1-r2),B=x0-A;return A*Math.exp(r1*t)+B*Math.exp(r2*t)}
  function timeEnd(){const wn=naturalFrequency(),T=2*Math.PI/wn;return state.zeta<1e-6?5*T:Math.min(20,Math.max(2*T,5/(Math.max(state.zeta,.08)*wn)))}
  function emit(){par.value={m:state.m,k:state.k,zeta:state.zeta,x0:state.x0,v0:state.v0,speed:state.speed};par.dispatchEvent(new CustomEvent('input'))}
  function updateControls(){Object.keys(inputs).forEach(key=>{inputs[key].value=state[key];values[key].textContent=key==='speed'?state[key].toFixed(1)+'×':key==='zeta'?state[key].toFixed(2):state[key].toFixed(2)});regimeReadout.textContent=regime()}
  function hidpi(canvas,ctx){const rect=canvas.getBoundingClientRect();canvas.width=Math.max(1,Math.round(rect.width*DPR));canvas.height=Math.max(1,Math.round(rect.height*DPR));ctx.setTransform(canvas.width/W,0,0,canvas.height/H,0,0)}
  function line(ctx,x1,y1,x2,y2,color,width,dash){ctx.save();ctx.strokeStyle=color;ctx.lineWidth=width;if(dash)ctx.setLineDash(dash);ctx.beginPath();ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);ctx.stroke();ctx.restore()}
  function spring(ctx,x1,x2,y){ctx.beginPath();const n=84;for(let i=0;i<=n;i++){const f=i/n,x=x1+(x2-x1)*f,yy=y+Math.sin(f*Math.PI*12)*10;if(i===0)ctx.moveTo(x,yy);else ctx.lineTo(x,yy)}ctx.stroke()}
  function drawScene(){sctx.clearRect(0,0,W,H);sctx.fillStyle='#000';sctx.fillRect(0,0,W,H);const eq=465,scale=205,x=displacement(state.t),bx=eq+scale*x,y=132;sctx.fillStyle='#9ca3af';sctx.font='14px sans-serif';sctx.textAlign='left';sctx.fillText('drag block → x₀',36,31);line(sctx,85,203,820,203,'#6b7280',2);line(sctx,eq,60,eq,221,'#4b5563',1,[5,5]);sctx.fillStyle='#9ca3af';sctx.textAlign='center';sctx.fillText('equilibrium',eq,244);line(sctx,76,70,76,194,'#e5e7eb',5);sctx.strokeStyle='#e5e7eb';sctx.lineWidth=2;spring(sctx,84,bx-31,y);sctx.strokeStyle='#a78bfa';sctx.lineWidth=3;line(sctx,88,172,bx-31,172,'#a78bfa',3);sctx.fillStyle='#a78bfa';sctx.fillRect(bx-49,164,18,17);sctx.fillStyle='#38bdf8';sctx.fillRect(bx-31,y-29,62,58);sctx.fillStyle='#f3f4f6';sctx.font='700 18px sans-serif';sctx.textAlign='center';sctx.fillText('m',bx,y+6);line(sctx,eq,87,bx,87,x>=0?'#38bdf8':'#f59e0b',2);sctx.fillStyle=x>=0?'#38bdf8':'#f59e0b';sctx.font='14px sans-serif';sctx.fillText('x = '+x.toFixed(3)+' m',Math.max(120,Math.min(780,(eq+bx)/2)),76);sctx.fillStyle='#9ca3af';sctx.font='13px sans-serif';sctx.fillText('spring k',220,113);sctx.fillText('dashpot c',235,192)}
  function drawPlot(){pctx.clearRect(0,0,W,H);pctx.fillStyle='#000';pctx.fillRect(0,0,W,H);const left=70,right=860,top=30,bottom=222,tw=right-left,th=bottom-top,tmax=timeEnd();let amp=.9;for(let i=0;i<=480;i++)amp=Math.max(amp,Math.abs(displacement(tmax*i/480)));amp*=1.16;for(let i=0;i<=4;i++){const y=top+i*th/4;line(pctx,left,y,right,y,'#1f2937',1);const val=amp*(1-2*i/4);pctx.fillStyle='#9ca3af';pctx.font='12px sans-serif';pctx.textAlign='right';pctx.fillText(val.toFixed(1),left-8,y+4)}for(let i=0;i<=5;i++){const x=left+i*tw/5;line(pctx,x,top,x,bottom,'#1f2937',1);pctx.fillStyle='#9ca3af';pctx.textAlign='center';pctx.fillText((tmax*i/5).toFixed(1),x,bottom+19)}const yzero=top+th/2;line(pctx,left,yzero,right,yzero,'#6b7280',1.5);pctx.strokeStyle='#38bdf8';pctx.lineWidth=3;pctx.beginPath();for(let i=0;i<=480;i++){const t=tmax*i/480,x=left+tw*i/480,y=yzero-displacement(t)/amp*th/2;if(i===0)pctx.moveTo(x,y);else pctx.lineTo(x,y)}pctx.stroke();const now=clamp(state.t,0,tmax),px=left+tw*now/tmax;line(pctx,px,top,px,bottom,'#ef4444',2);pctx.fillStyle='#ef4444';pctx.font='13px sans-serif';pctx.textAlign='left';pctx.fillText('now',Math.min(px+5,right-28),top+14);pctx.fillStyle='#9ca3af';pctx.textAlign='center';pctx.fillText('time (s)',(left+right)/2,H-12);pctx.save();pctx.translate(16,(top+bottom)/2);pctx.rotate(-Math.PI/2);pctx.fillText('displacement x (m)',0,0);pctx.restore()}
  function draw(){drawScene();drawPlot();timeReadout.textContent=state.t.toFixed(2)+' s'}
  function tick(now){if(state.playing){if(state.last)state.t+=((now-state.last)/1000)*state.speed;if(state.t>=timeEnd()){state.t=timeEnd();state.playing=false;playButton.textContent='Play'}state.last=now;draw();if(state.playing)requestAnimationFrame(tick)}}
  function resetClock(){state.t=0;state.last=0;draw()}
  Object.keys(inputs).forEach(key=>inputs[key].addEventListener('input',()=>{state[key]=Number(inputs[key].value);resetClock();updateControls();emit()}))
  playButton.addEventListener('click',()=>{state.playing=!state.playing;state.last=0;playButton.textContent=state.playing?'Pause':'Play';if(state.playing)requestAnimationFrame(tick);draw()})
  par.querySelector('#shm-restart').addEventListener('click',()=>{state.playing=false;playButton.textContent='Play';resetClock()})
  par.querySelector('#shm-home').addEventListener('click',()=>{state.playing=false;playButton.textContent='Play';state.x0=.70;state.v0=0;resetClock();updateControls();emit()})
  par.querySelectorAll('.shm-preset').forEach(button=>button.addEventListener('click',()=>{state.zeta=Number(button.dataset.z);state.playing=false;playButton.textContent='Play';resetClock();updateControls();emit()}))
  function pointerX(event){const rect=scene.getBoundingClientRect();return (event.clientX-rect.left)*W/rect.width}
  scene.addEventListener('pointerdown',event=>{const px=pointerX(event),bx=465+205*displacement(state.t);if(Math.abs(px-bx)<62){state.dragging=true;state.playing=false;playButton.textContent='Play';scene.setPointerCapture(event.pointerId)}})
  scene.addEventListener('pointermove',event=>{if(!state.dragging)return;state.t=0;state.x0=clamp((pointerX(event)-465)/205,-1.2,1.2);updateControls();draw()})
  function endDrag(event){if(state.dragging){state.dragging=false;try{scene.releasePointerCapture(event.pointerId)}catch(e){}emit();draw()}}
  scene.addEventListener('pointerup',endDrag);scene.addEventListener('pointercancel',endDrag)
  const observer=new ResizeObserver(()=>{hidpi(scene,sctx);hidpi(plot,pctx);draw()});observer.observe(scene);observer.observe(plot);par.style.setProperty('--shm-totalw',Math.round(Math.min(window.innerWidth*.8,par.clientWidth||960))+'px');updateControls();hidpi(scene,sctx);hidpi(plot,pctx);draw();emit()
  </script>
</div>
        """)
    end

    const _shm_ready = true
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7"
"""

# ╔═╡ Cell order:
# ╠═27ab0195-38a6-4235-a782-8a5b8a7bac01
# ╠═548b0e6a-7284-49e1-88f7-e5ef9f8a2de5
# ╟─9a3af4df-76fe-4381-87b7-56a3b43925cc
# ╟─c1d09d0d-2c25-4db5-9298-e4a6e9fa550a
# ╠═6e3ba5f5-70ec-45bf-9a56-e0a5990d453d
# ╠═6fba3aa2-6192-4f45-8931-2e715e213786
# ╠═94ab3b93-4831-48d3-8ae8-b3f82e8738a2
# ╟─b4f616c5-925c-40ea-9a74-1b1f68848aed
# ╟─8ed80c0f-8177-4229-ad4e-a1a9f075718f
# ╟─174420d3-9c28-4d3f-93b2-6b9a267a222a
# ╟─405af09e-aa7b-46f5-9877-e9dfc4e1427c
# ╠═73c5bd9b-9a8f-41a9-b2a4-d076a15bbbca
# ╠═bc14b27d-af7e-499f-b169-cb4814810a5d
# ╠═90b3d802-0993-4df2-8b48-3c3c462d5c8f
# ╠═f6d03a7b-3c89-4799-bcd3-8c799711ea1e
# ╟─6694ed08-e429-41b8-b8ac-a2f9b417954b
# ╠═de2db126-21d0-4f72-9d26-79ba19d1b65b
# ╟─00000000-0000-0000-0000-000000000001
