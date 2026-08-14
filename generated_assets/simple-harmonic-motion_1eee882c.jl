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
  function hidpi(canvas,ctx){const rect=canvas.getBoundingClientRect();canvas.width=Math.max(1,Math.round(rect.width*DPR));canvas.height=Math.max(1,Math.round(rect.height*DPR));const scale=Math.min(canvas.width/W,canvas.height/H),dx=(canvas.width-W*scale)/2,dy=(canvas.height-H*scale)/2;ctx.setTransform(scale,0,0,scale,dx,dy)}
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

# ╔═╡ 6e3ba5f5-70ec-45bf-9a56-e0a5990d453d
begin
    _shm_ready
    WideCell(@bind shm_state SHMInput(); max_width=1500)
end

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

# ╔═╡ Cell order:
# ╠═27ab0195-38a6-4235-a782-8a5b8a7bac01
# ╠═548b0e6a-7284-49e1-88f7-e5ef9f8a2de5
# ╟─9a3af4df-76fe-4381-87b7-56a3b43925cc
# ╟─c1d09d0d-2c25-4db5-9298-e4a6e9fa550a
# ╟─6e3ba5f5-70ec-45bf-9a56-e0a5990d453d
# ╠═6fba3aa2-6192-4f45-8931-2e715e213786
# ╟─94ab3b93-4831-48d3-8ae8-b3f82e8738a2
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
# ╟─00000000-0000-0000-0000-000000000002
