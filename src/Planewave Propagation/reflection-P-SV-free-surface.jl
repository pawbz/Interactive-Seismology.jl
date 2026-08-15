### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "P-SV Free-surface Reflection"
#> layout = "layout.jlhtml"
#> tags = ["planewaves"]
#> description = "In this notebook, we shall investigate the interaction of plane P waves with Earth's free surface. "

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

# ╔═╡ ab40f79c-3d8a-11ed-0697-a7b794dbba99
begin
    using Symbolics
    using PlutoUI
end

# ╔═╡ 08429397-3964-4600-bc14-c45d22c915ec
TableOfContents()

# ╔═╡ 32a757ff-aa7a-41d5-b8b3-6ac9e0125875
md"""
# P-SV Free-surface Reflection
In this notebook, we will delve into the analysis of the reflection coefficients of P and S waves at a free surface, which is a solid-vacuum boundary. By assuming a free surface at `z=0`, we will apply stress-free boundary conditions to derive the expressions for the reflection coefficients of both P and S waves. This analysis is crucial for understanding the behavior of seismic waves as they interact with the Earth's surface.


##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 4476bf78-39e3-4674-a152-db19fe80929a
md"Spatial coordinates, time, and angular frequency."

# ╔═╡ 927dcc43-202c-4ad2-a76c-837d41f1ed6c
@syms x::Real z::Real ω::Real t::Real

# ╔═╡ c2eacb91-38e8-428e-9247-691950668bc3
@syms ı::Complex{Real} # imaginary unit, going to substitute with im later

# ╔═╡ 0bea69f9-3ffa-4695-a5eb-6e962d2a81ce
md"Differential operators."

# ╔═╡ e8729ceb-e85f-4c21-89d2-f15ec69a840f
begin
    Dx = Differential(x)
    Dz = Differential(z)
end

# ╔═╡ 670d5198-f810-472e-8db8-b13385ca294a
md"In 2D, a harmonic plane wave with an frequency `ω`, amplitude `A`, horizontal slowness `p` and vertical slowness `η` is defined below."

# ╔═╡ 6043b904-6991-4776-bc1b-13eda3fcc936
"""
	plane(p, η, A)

A harmonic plane-wave potential with horizontal slowness `p`, vertical slowness `η`, and
amplitude `A`, sharing the notebook-wide symbolic angular frequency `ω` and coordinates
`x`, `z`, `t`.
"""
plane(p, η, A) = A * exp(ı * ω * (t - (p * x + η * z)))

# ╔═╡ c7c52926-6e2e-4c4c-a01f-09f57c2ececd
md"The horizontal slowness i.e., the ray parameter is denoted using `p`."

# ╔═╡ 33c6aa70-87cb-464a-88f7-b82d17476a2f
md"The vertical component of the slowness vector in the first and second layer are denoted using `η` and `η_t`. These plane waves satisfy the scalar wave equation only if the dispersion relation
```math
p^2 + η^2 = \frac{1}{β^2}
```
is satisfied."

# ╔═╡ fbc6cd6c-0b3b-44e9-b2a6-5edb0eeced1b
@syms A A₁ A₂

# ╔═╡ 8dad08c2-b9a8-40ec-a8ef-c9383e5ec7b1
@syms p::Real

# ╔═╡ 13b5abc9-10be-4ef1-817d-bcee1271c89e
@syms η₁ η₂

# ╔═╡ aae7b893-4b33-44ae-b01c-6345a1180d0d
md"""
## Potentials : Incident P-wave, Reflected P-wave and Reflected SV-wave 
For the P-SV system on a free surface, all the waves are traveling in the same medium. The velocity of the P-wave is denoted by `α` and the velocity of the S-wave is denoted by `β`. The plane waves representing potentials share the same horizontal slowness `p`.

- Now, we create an incident wave using the `plane` function with parameters `p`, `η₁`, and `A`, representing the horizontal slowness, verticle slowness, and amplitude, respectively. 
- Then, two reflected waves, `reflect1` and `reflect2`, are also defined using the `plane` function. The first reflected wave, `reflect1`, uses the parameters `p`, `-η₁`, and `A₁`, indicating a reflection with opposite vertical slowness and a different amplitude. 
- Similarly, the second reflected wave, `reflect2`, is defined with `p`, `-η₂`, and `A₂`, representing another reflection with a different vertical slowness and amplitude. This setup models the behavior of P and S incident and reflected waves at a free surface.
"""

# ╔═╡ f488f9f3-e73b-42ae-b3c2-2262661fd839
incident = plane(p, η₁, A)

# ╔═╡ 4dc428b0-f141-4e94-9edb-e847780333aa
reflect1 = plane(p, -η₁, A₁)

# ╔═╡ 8940ee90-2dd1-448c-92cd-09fd1ebbaea7
reflect2 = plane(p, -η₂, A₂)

# ╔═╡ a1213dc4-3ea9-4b26-aa82-a15cc3e72401
md"""
## Particle Displacements
We now write expressions for particle displacements given plane wave potentials.
Depending on whether the incident wave is P or S, these expressions change.
"""

# ╔═╡ 52d2ba9b-826a-4805-9583-f1a80b10a926
md"""
## Free-surface Condition
At the free-surface i.e. z=0 , kinematic boundary condition does not exist.
We have two dynamic boundary conditions `σ_{zz}=0` and `σ_{xz}=0` and two unknowns (`A_{1}`,`A_{2}`) to estimate (when assuming `A=1`).
"""

# ╔═╡ 7cf515f8-2496-4bde-b3b3-9d39d971764a
@syms λ::Real μ::Real

# ╔═╡ 5729b459-b283-41ce-95be-c4d33a7c28c0
md"## Example"

# ╔═╡ 281cb873-760c-411f-98b5-1c64d218e7e9
md"""
We shall begin defining a function that computes the ray parameter, given the
angle of incidence `θ`.
"""

# ╔═╡ 3d50985b-b79a-45ab-ba1f-5287936c56d9
@syms θ::Real

# ╔═╡ a7cb9cd5-7013-461b-960d-5da73db62aea
md"We now compute the vertical component of the slowness vector for both the reflected waves using the 
dispersion relation."

# ╔═╡ 8c81ddb5-bf4d-4610-bfea-3d1a27ffd61f
md"""
We can finally update the expressions of reflection coefficients using input parameters and plot them.
"""

# ╔═╡ eeee5555-5555-5555-5555-555555555555
md"""
### Verifying the coefficients
At normal incidence (`θ=0`), there is no mode conversion: a P wave striking the surface
head-on reflects only as P (and likewise for S), so the mode-converted amplitude should
vanish, and since a free surface reflects all incident energy, the direct reflection
coefficient's magnitude should be exactly 1.
"""

# ╔═╡ dba7ea14-e0dd-4dc4-ad9c-4627fd16cc62
md"## Appendix"

# ╔═╡ ea3b7089-bda8-4694-8042-98534b1739bd
Markdown.MD(Markdown.Admonition("warning", "Sign convention",
    [md"""
The sign of the vertical slowness in the transmitted field above is chosen to prevent the exponential growth of the wavefield away from the boundary.
"""]
))

# ╔═╡ 602f13f9-6d14-41fd-9183-b8255d64399b
Markdown.MD(Markdown.Admonition("note", "Observations",
    [md"""
- The angle of reflection of the S-wave is always less than that of the P-reflected wave. This is because the vertical slowness of the S-reflected wave (`η_β`) is always greater than that of the P-reflected wave.
- The wavenumber of the S-reflected wave is higher than that of the P-wave. The wavenumber of the S-reflected wave increases as `β` becomes smaller compared to `α`.
- When `β` is very small, the amplitude of the S-reflected wave is negligible.
- A negative reflection coefficient indicates a phase change of the reflected wave by π.
- The P-reflection coefficient and S-reflection coefficient follow the conservation of energy. The P-reflection coefficient is maximum when the S-reflection coefficient is minimum.
"""]
))

# ╔═╡ aaaa7777-7777-7777-7777-777777777777
md"## The Interactive Widget"

# ╔═╡ dddd4444-4444-4444-4444-444444444444
begin
    """A draggable coefficient-vs-angle canvas: the incidence angle is set by
    dragging directly on the reflectivity curves, not by a slider. Only the
    medium properties (`alpha`, `beta`, `rho`) and the incident-wave type feed
    back to Julia -- angle, frequency, and which waves are drawn are pure
    client-side state, since a full angle sweep of the coefficients is pushed to
    the widget once (see the cell above `## Appendix`) and everything else is
    cheap trigonometry over that sweep."""
    struct PSVWidgetInput
        alpha::Float64
        beta::Float64
        rho::Float64
        incident_type::String
    end

    PSVWidgetInput(; alpha=3.9, beta=2.0, rho=2.9, incident_type="P") =
        PSVWidgetInput(Float64(alpha), Float64(beta), Float64(rho), incident_type)

    Base.get(w::PSVWidgetInput) = Dict{String,Any}(
        "alpha" => w.alpha, "beta" => w.beta, "rho" => w.rho, "incident_type" => w.incident_type,
    )

    function Base.show(io::IO, ::MIME"text/html", w::PSVWidgetInput)
        write(io, """
<div id="psvwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#psvwidget) {
      width: min(80vw, 1100px) !important;
      margin-left: calc((100% - min(80vw, 1100px)) / 2) !important;
    }
    #psvwidget { width: 100%; box-sizing: border-box; color: #d1d5db; font: 14px sans-serif; }
    #psvwidget .psv-title { width: 100%; box-sizing: border-box; text-align: center; margin-bottom: 10px;
      background: #0a0f18; border: 1px solid #3b5c85; border-radius: 6px; padding: 10px 14px; }
    #psvwidget .psv-title-desc { font-size: 17px; font-weight: 700; color: #e5e7eb; }
    #psvwidget .psv-title-hint { font-size: 13px; color: #9ca3af; margin-top: 3px; }
    #psvwidget .psv-workspace { display: flex; gap: 12px; align-items: flex-start; justify-content: center;
      width: 100%; flex-wrap: wrap; }
    #psvwidget .psv-controls { width: min(var(--totalw,900px),100%); margin-top: 10px; display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 8px; font: 14px sans-serif; }
    #psvwidget .psv-control-group { box-sizing: border-box; background: #050505; border: 1px solid #2f3744;
      border-radius: 6px; padding: 10px 12px; }
    #psvwidget .psv-control-title { font-weight: 700; color: #e5e7eb; margin-bottom: 8px; font-size: 18px; }
    #psvwidget .psv-control-row { display: grid; grid-template-columns: minmax(60px,90px) minmax(70px,1fr) minmax(50px,70px);
      gap: 6px; align-items: center; margin: 6px 0; }
    #psvwidget .psv-control-row input[type=range] { width: 100%; min-width: 0; }
    #psvwidget .psv-value { color: #d1d5db; text-align: right; font-variant-numeric: tabular-nums; }
    #psvwidget .psv-actions { display: flex; gap: 8px; align-items: center; flex-wrap: wrap; }
    #psvwidget button { border-radius: 4px; border: 1px solid #9ca3af; background: #606060; color: #f3f4f6;
      padding: 6px 12px; font-size: 14px; cursor: pointer; }
    #psvwidget label { color: #d1d5db; }
  </style>
  <div class="psv-title">
    <div class="psv-title-desc">Drag the angle axis on the coefficient plot to see how the reflected P and S waves change.</div>
    <div class="psv-title-hint">drag left/right on the coefficient plot &middot; hover a curve to identify it &middot; press play to animate</div>
  </div>
  <div class="psv-workspace">
    <canvas id="psvrefl" style="background:#000;border:1px solid #374151;border-radius:6px;display:block"></canvas>
    <canvas id="psvwave" style="background:#000;border:1px solid #374151;border-radius:6px;display:block"></canvas>
  </div>
  <div class="psv-controls">
    <div class="psv-control-group">
      <div class="psv-control-title">Medium</div>
      <label class="psv-control-row"><span>α (km/s)</span><input type="range" id="psv-alpha" min="2" max="6" step="0.1" value="$(w.alpha)"><span id="psv-alpha-v" class="psv-value">$(w.alpha)</span></label>
      <label class="psv-control-row"><span>β (km/s)</span><input type="range" id="psv-beta" min="2" max="6" step="0.1" value="$(w.beta)"><span id="psv-beta-v" class="psv-value">$(w.beta)</span></label>
      <label class="psv-control-row"><span>ρ (g/cm³)</span><input type="range" id="psv-rho" min="2" max="6" step="0.1" value="$(w.rho)"><span id="psv-rho-v" class="psv-value">$(w.rho)</span></label>
    </div>
    <div class="psv-control-group">
      <div class="psv-control-title">Incidence</div>
      <div class="psv-actions">
        <label><input type="radio" name="psv-itype" id="psv-itype-p" $(w.incident_type == "P" ? "checked" : "")> incident <b style="color:#f97316">P</b></label>
        <label><input type="radio" name="psv-itype" id="psv-itype-s" $(w.incident_type == "S" ? "checked" : "")> incident <b style="color:#38bdf8">S</b></label>
      </div>
      <label class="psv-control-row"><span>ω</span><input type="range" id="psv-omega" min="0.1" max="2" step="0.1" value="0.3"><span id="psv-omega-v" class="psv-value">0.3</span></label>
      <div class="psv-actions"><button id="psv-play" type="button">Play</button></div>
    </div>
    <div class="psv-control-group">
      <div class="psv-control-title">Display</div>
      <div class="psv-actions">
        <label><input type="checkbox" id="psv-show-i" checked> Incident</label>
        <label><input type="checkbox" id="psv-show-r1" checked> Reflected 1</label>
        <label><input type="checkbox" id="psv-show-r2" checked> Reflected 2</label>
      </div>
      <div class="psv-actions" style="margin-top:8px"><button id="psv-reset" type="button">Reset defaults</button></div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1100)
  const totalW = Math.max(700, availW)
  par.style.setProperty('--totalw', Math.round(totalW)+'px')
  const SEC = Math.round(Math.min(totalW/2 - 10, 380))
  const DPR = Math.min(window.devicePixelRatio || 1, 2)

  let alpha = $(w.alpha), beta = $(w.beta), rho = $(w.rho), incidentType = "$(w.incident_type)"
  let thetaDeg = 40, omega = 0.3
  let showI = true, showR1 = true, showR2 = true
  let hoverCurve = -1, playing = false, tPhase = 0, rafId = null
  let sweep = null   // filled in by the 'psv-results' push from Julia, below

  const reflCvs = par.querySelector('#psvrefl'), rctx = reflCvs.getContext('2d')
  const waveCvs = par.querySelector('#psvwave'), wctx = waveCvs.getContext('2d')
  function hidpi(cv, cx, w, h){
    cv.width = Math.round(w*DPR); cv.height = Math.round(h*DPR)
    cv.style.width = w+'px'; cv.style.height = h+'px'
    cx.setTransform(DPR,0,0,DPR,0,0)
  }
  hidpi(reflCvs, rctx, SEC, SEC)
  hidpi(waveCvs, wctx, SEC, SEC)

  // ---- reflectivity canvas: quarter-circle polar plot (angle 0-90deg sweeping
  // clockwise from "up"), matching the notebook's own Reference Plots section --
  // the radial axis is the diverging coefficient value itself, so the pole
  // (r=RVMIN) sits at the plot's inner corner and the outer arc is r=RVMAX.
  const RPOLE = {x: 46, y: SEC-30}
  const RVMIN = -2, RVMAX = 2
  const RMAX = Math.max(60, Math.min(SEC-RPOLE.x-34, RPOLE.y-32))
  function angDir(thetaDeg){
    const rad = thetaDeg*Math.PI/180
    return [Math.sin(rad), -Math.cos(rad)]   // 0deg -> straight up, 90deg -> straight right
  }
  function polarXY(thetaDeg, v){
    const [dx,dy] = angDir(thetaDeg)
    const r = RMAX * (v-RVMIN)/(RVMAX-RVMIN)
    return [RPOLE.x + dx*r, RPOLE.y + dy*r]
  }
  function polarThetaOf(px, py){
    const dx = px-RPOLE.x, dy = py-RPOLE.y
    return Math.max(0, Math.min(90, Math.atan2(dx, -dy) * 180/Math.PI))
  }

  const CURVES = [
    {key:'A1', part:'abs', label:'|P-reflection|', color:'#f97316', dash:[]},
    {key:'A1', part:'re',  label:'Re(P-reflection)', color:'#f97316', dash:[6,3]},
    {key:'A1', part:'im',  label:'Im(P-reflection)', color:'#f97316', dash:[1,3]},
    {key:'A2', part:'abs', label:'|S-reflection|', color:'#38bdf8', dash:[]},
    {key:'A2', part:'re',  label:'Re(S-reflection)', color:'#38bdf8', dash:[6,3]},
    {key:'A2', part:'im',  label:'Im(S-reflection)', color:'#38bdf8', dash:[1,3]},
  ]

  function curveValue(c, i){
    const re = sweep[c.key+'_re'][i], im = sweep[c.key+'_im'][i]
    if(c.part === 'abs') return Math.hypot(re, im)
    if(c.part === 're') return re
    return im
  }

  function sweepIndex(theta){ return Math.max(0, Math.min(sweep.theta_deg.length-1, Math.round(theta*2))) }

  function drawReflectivity(){
    rctx.clearRect(0,0,SEC,SEC)

    // radial gridlines: concentric quarter-circle arcs at nice diverging values
    const rticks = [-2,-1,0,1,2]
    rctx.font = '10px sans-serif'; rctx.textAlign = 'center'
    rticks.forEach(v => {
      const r = RMAX*(v-RVMIN)/(RVMAX-RVMIN)
      rctx.strokeStyle = v===0 ? '#4b5563' : '#242b38'
      rctx.lineWidth = v===0 ? 1.3 : 1
      rctx.beginPath(); rctx.arc(RPOLE.x, RPOLE.y, r, -Math.PI/2, 0); rctx.stroke()
      rctx.fillStyle = '#6b7280'
      rctx.fillText(String(v), RPOLE.x+r, RPOLE.y+14)
    })

    // angular gridlines (spokes) every 15deg, with degree labels -- clamped to
    // stay inside the canvas (see pluto-widget-style's tick-label-clipping note)
    rctx.strokeStyle = '#242b38'; rctx.lineWidth = 1
    for(let d=0; d<=90; d+=15){
      const [x,y] = polarXY(d, RVMAX)
      rctx.beginPath(); rctx.moveTo(RPOLE.x,RPOLE.y); rctx.lineTo(x,y); rctx.stroke()
      const [lx,ly] = polarXY(d, RVMAX*1.1)
      const label = d+'°'
      rctx.fillStyle = '#6b7280'; rctx.font = '10px sans-serif'; rctx.textAlign = 'center'
      const tw = rctx.measureText(label).width
      const clx = Math.max(tw/2+2, Math.min(SEC-tw/2-2, lx))
      const cly = Math.max(10, Math.min(SEC-2, ly))
      rctx.fillText(label, clx, cly)
    }
    rctx.strokeStyle = '#374151'; rctx.lineWidth = 1.2
    rctx.beginPath(); rctx.arc(RPOLE.x, RPOLE.y, RMAX, -Math.PI/2, 0); rctx.stroke()
    rctx.textAlign = 'left'

    if(!sweep){
      rctx.fillStyle = '#6b7280'; rctx.font = '12px sans-serif'
      rctx.fillText('computing...', 12, 18)
      return
    }

    CURVES.forEach((c, idx) => {
      const isHover = idx === hoverCurve
      rctx.globalAlpha = (hoverCurve === -1 || isHover) ? 1 : 0.2
      rctx.strokeStyle = c.color; rctx.lineWidth = isHover ? 2.6 : 1.6
      rctx.setLineDash(c.dash)
      rctx.beginPath()
      sweep.theta_deg.forEach((th, i) => {
        const [x,y] = polarXY(th, curveValue(c, i))
        i === 0 ? rctx.moveTo(x,y) : rctx.lineTo(x,y)
      })
      rctx.stroke()
      rctx.setLineDash([])
      rctx.globalAlpha = 1
    })

    // draggable angle cursor -- a filled wedge (matching the reference plots'
    // own "Animation" marker) plus the spoke line on top
    const thetaRad = thetaDeg*Math.PI/180
    const [cx,cy] = polarXY(thetaDeg, RVMAX)
    rctx.fillStyle = 'rgba(228,255,135,0.22)'
    rctx.beginPath()
    rctx.moveTo(RPOLE.x, RPOLE.y)
    rctx.arc(RPOLE.x, RPOLE.y, RMAX, -Math.PI/2+thetaRad-0.025, -Math.PI/2+thetaRad+0.025)
    rctx.closePath(); rctx.fill()
    rctx.strokeStyle = '#e5e7eb'; rctx.lineWidth = 1.4
    rctx.beginPath(); rctx.moveTo(RPOLE.x,RPOLE.y); rctx.lineTo(cx,cy); rctx.stroke()
    rctx.fillStyle = '#e5e7eb'; rctx.font = '12px sans-serif'
    rctx.fillText(Math.round(thetaDeg)+'°', Math.min(cx+4, SEC-30), Math.max(cy-4, 12))

    if(hoverCurve >= 0){
      const c = CURVES[hoverCurve]
      const v = curveValue(c, sweepIndex(thetaDeg))
      const label = c.label + '  ' + v.toFixed(3)
      rctx.font = '12px sans-serif'
      const tw = rctx.measureText(label).width
      const tx = Math.min(cx+8, SEC-tw-14), ty = 20
      rctx.fillStyle = 'rgba(11,18,32,0.9)'; rctx.fillRect(tx-5, ty-13, tw+10, 18)
      rctx.strokeStyle = '#374151'; rctx.lineWidth = 1; rctx.strokeRect(tx-5, ty-13, tw+10, 18)
      rctx.fillStyle = '#e5e7eb'; rctx.fillText(label, tx, ty)
    }
  }

  function distToPolyline(mx, my, c){
    let best = 1e9
    sweep.theta_deg.forEach((th, i) => {
      const [x,y] = polarXY(th, curveValue(c, i))
      best = Math.min(best, Math.hypot(mx-x, my-y))
    })
    return best
  }
  function nearestCurve(mx, my){
    if(!sweep) return -1
    let best = -1, bestD = 8
    CURVES.forEach((c, idx) => {
      const d = distToPolyline(mx, my, c)
      if(d < bestD){ bestD = d; best = idx }
    })
    return best
  }

  // ---- wavefield canvas: x in [-200,200], z in [-200,0] (free surface at top) ----
  const REARTH_X = 200, REARTH_Z = 200
  const off = document.createElement('canvas')
  const OFFW = 110, OFFH = 90
  off.width = OFFW; off.height = OFFH
  const offCtx = off.getContext('2d')

  // Diverging "seismic" colormap: blue (negative) - white (zero) - red (positive).
  function seismicColor(v, vmax){
    const t = Math.max(-1, Math.min(1, vmax > 0 ? v/vmax : 0))
    if(t >= 0){
      return [255, Math.round(255*(1-t)), Math.round(255*(1-t))]
    } else {
      const s = -t
      return [Math.round(255*(1-s)), Math.round(255*(1-s)), 255]
    }
  }

  function drawWavefield(){
    wctx.clearRect(0,0,SEC,SEC)
    if(!sweep){
      wctx.fillStyle = '#6b7280'; wctx.font = '12px sans-serif'
      wctx.fillText('computing...', 12, 18)
      return
    }
    const i0 = sweepIndex(thetaDeg)
    const p = sweep.p[i0]
    const eta1 = sweep.eta1_re[i0]   // incident/reflected-1 leg is always real for a valid 0-90deg angle
    const eta2_re = sweep.eta2_re[i0], eta2_im = sweep.eta2_im[i0]
    const isReal2 = Math.abs(eta2_im) < 1e-6
    const eU2re = isReal2 ? -eta2_re : eta2_re
    const eU2im = isReal2 ? -eta2_im : eta2_im
    const a1 = sweep.A1_re[i0], b1 = sweep.A1_im[i0]
    const a2 = sweep.A2_re[i0], b2 = sweep.A2_im[i0]
    const t = tPhase

    const img = offCtx.createImageData(OFFW, OFFH)
    // first pass: compute raw values to find a robust color scale
    const vals = new Float64Array(OFFW*OFFH)
    let vmax = 1e-6
    for(let j=0;j<OFFH;j++){
      const z = -(j/(OFFH-1))*REARTH_Z
      for(let i=0;i<OFFW;i++){
        const x = -REARTH_X + (i/(OFFW-1))*2*REARTH_X
        let U = 0
        if(showI){
          const phase = omega*(t - p*x - eta1*z)
          U += Math.cos(phase)
        }
        if(showR1){
          const phase = omega*(t - p*x + eta1*z)
          U += a1*Math.cos(phase) - b1*Math.sin(phase)
        }
        if(showR2){
          const mag = Math.exp(omega*eU2im*z)
          const phase = omega*(t - p*x) - omega*eU2re*z
          U += mag*(a2*Math.cos(phase) - b2*Math.sin(phase))
        }
        vals[j*OFFW+i] = U
        vmax = Math.max(vmax, Math.abs(U))
      }
    }
    for(let k=0;k<vals.length;k++){
      const [r,g,b] = seismicColor(vals[k], vmax)
      img.data[k*4] = r; img.data[k*4+1] = g; img.data[k*4+2] = b; img.data[k*4+3] = 255
    }
    offCtx.putImageData(img, 0, 0)
    wctx.imageSmoothingEnabled = true
    wctx.drawImage(off, 0, 0, OFFW, OFFH, 0, 0, SEC, SEC)

    // free-surface line + labels
    wctx.strokeStyle = '#facc15'; wctx.lineWidth = 2
    wctx.beginPath(); wctx.moveTo(0,1); wctx.lineTo(SEC,1); wctx.stroke()
    wctx.fillStyle = '#facc15'; wctx.font = '12px sans-serif'
    wctx.fillText('Free surface', SEC*0.55, 16)
    wctx.fillStyle = '#e5e7eb'
    wctx.fillText('(α,β,ρ)', 10, SEC-10)
  }

  function redraw(){ drawReflectivity(); drawWavefield() }

  function emit(){
    par.value = {alpha, beta, rho, incident_type: incidentType}
    par.dispatchEvent(new CustomEvent('input'))
  }

  // ---- interaction: dragging the angle cursor on the reflectivity canvas ----
  let draggingTheta = false
  reflCvs.addEventListener('mousedown', e=>{ draggingTheta = true; thetaDeg = polarThetaOf(e.offsetX, e.offsetY); redraw() })
  reflCvs.addEventListener('mousemove', e=>{
    if(draggingTheta){
      thetaDeg = polarThetaOf(e.offsetX, e.offsetY)
      redraw()
    } else {
      hoverCurve = nearestCurve(e.offsetX, e.offsetY)
      redraw()
    }
  })
  window.addEventListener('mouseup', ()=>{ draggingTheta = false })
  reflCvs.addEventListener('mouseleave', ()=>{ if(hoverCurve !== -1){ hoverCurve = -1; redraw() } })

  // ---- controls ----
  par.querySelector('#psv-alpha').addEventListener('input', e=>{
    alpha = parseFloat(e.target.value); par.querySelector('#psv-alpha-v').textContent = alpha.toFixed(1); emit()
  })
  par.querySelector('#psv-beta').addEventListener('input', e=>{
    beta = parseFloat(e.target.value); par.querySelector('#psv-beta-v').textContent = beta.toFixed(1); emit()
  })
  par.querySelector('#psv-rho').addEventListener('input', e=>{
    rho = parseFloat(e.target.value); par.querySelector('#psv-rho-v').textContent = rho.toFixed(1); emit()
  })
  par.querySelector('#psv-itype-p').addEventListener('change', ()=>{ incidentType = 'P'; emit() })
  par.querySelector('#psv-itype-s').addEventListener('change', ()=>{ incidentType = 'S'; emit() })
  par.querySelector('#psv-omega').addEventListener('input', e=>{
    omega = parseFloat(e.target.value); par.querySelector('#psv-omega-v').textContent = omega.toFixed(1); drawWavefield()
  })
  par.querySelector('#psv-show-i').addEventListener('change', e=>{ showI = e.target.checked; drawWavefield() })
  par.querySelector('#psv-show-r1').addEventListener('change', e=>{ showR1 = e.target.checked; drawWavefield() })
  par.querySelector('#psv-show-r2').addEventListener('change', e=>{ showR2 = e.target.checked; drawWavefield() })

  const playBtn = par.querySelector('#psv-play')
  function stepAnim(){
    tPhase += 0.12
    drawWavefield()
    rafId = requestAnimationFrame(stepAnim)
  }
  playBtn.addEventListener('click', ()=>{
    playing = !playing
    playBtn.textContent = playing ? 'Pause' : 'Play'
    if(playing){ rafId = requestAnimationFrame(stepAnim) } else if(rafId){ cancelAnimationFrame(rafId); rafId = null }
  })

  par.querySelector('#psv-reset').addEventListener('click', ()=>{
    alpha = $(w.alpha); beta = $(w.beta); rho = $(w.rho); incidentType = "$(w.incident_type)"
    thetaDeg = 40; omega = 0.3; showI = true; showR1 = true; showR2 = true
    par.querySelector('#psv-alpha').value = alpha; par.querySelector('#psv-alpha-v').textContent = alpha.toFixed(1)
    par.querySelector('#psv-beta').value = beta; par.querySelector('#psv-beta-v').textContent = beta.toFixed(1)
    par.querySelector('#psv-rho').value = rho; par.querySelector('#psv-rho-v').textContent = rho.toFixed(1)
    par.querySelector('#psv-omega').value = omega; par.querySelector('#psv-omega-v').textContent = omega.toFixed(1)
    par.querySelector('#psv-itype-p').checked = (incidentType === 'P')
    par.querySelector('#psv-itype-s').checked = (incidentType === 'S')
    par.querySelector('#psv-show-i').checked = true
    par.querySelector('#psv-show-r1').checked = true
    par.querySelector('#psv-show-r2').checked = true
    redraw(); emit()
  })

  window.addEventListener('psv-results', e=>{
    const d = e.detail ? JSON.parse(e.detail) : null
    if(!d) return
    sweep = d
    redraw()
  })

  redraw(); emit()
</script>
""")
    end

    const _psv_ready = true
end

# ╔═╡ bbbb2222-2222-2222-2222-222222222222
begin
    # `PSVWidgetInput` is defined in the Appendix, displayed below this cell -- a
    # bare reference forces Pluto to run that cell first on a cold restart. See
    # "the one thing that will silently break on a fresh restart" in
    # pluto-widget-SKILL.md.
    _psv_ready
    @bind _psv PSVWidgetInput()
end

# ╔═╡ afdb5b7d-d670-4a98-a91d-3ff638fb0294
begin
    # Only the medium parameters and the incident-wave type need a Julia recompute
    # (they feed the direct numerical coefficients below) -- angle, frequency, and
    # display toggles live entirely client-side inside the widget (see its
    # docstring for why: everything downstream of theta is pre-swept once and
    # pushed to the widget, so dragging theta never touches Julia at all).
    αin = Float64(_psv["alpha"])
    βin = Float64(_psv["beta"])
    ρin = Float64(_psv["rho"])
    incident_waves = String(_psv["incident_type"])
end


# ╔═╡ b4a603fb-5ee7-4c97-8e1d-fd0a89692503
begin
    if (incident_waves == "P")
        UxInc = expand_derivatives.(Dx.(incident))
        UzInc = expand_derivatives.(Dz.(incident))
    else
        UxInc = expand_derivatives.(-Dz.(incident))
        UzInc = expand_derivatives.(Dz.(incident))
    end
end

# ╔═╡ 040f1440-14c4-491b-b792-570aa3c2080b
σzz_incident = expand_derivatives.(λ .* Dx.(UxInc) .+ λ .* Dz.(UzInc)) + expand_derivatives.(2 .* μ .* Dz.(UzInc))

# ╔═╡ f51dfb73-b5ef-4de8-b4c3-fc8814a31485
σxz_incident = expand_derivatives.(μ .* Dx.(UzInc) .+ μ .* Dz.(UxInc))

# ╔═╡ 74484bb8-a9dc-4b54-894e-3bf86c18278e
begin
    if (incident_waves == "P")
        UxRef1 = expand_derivatives.(Dx.(reflect1))
        UzRef1 = expand_derivatives.(Dz.(reflect1))
    else
        UxRef1 = expand_derivatives.(-Dz.(reflect1))
        UzRef1 = expand_derivatives.(Dx.(reflect1))
    end
end

# ╔═╡ 3525bf50-b0e5-4c3b-ba75-77689a5799f4
σzz_reflect1 = expand_derivatives.(λ .* Dx.(UxRef1) .+ λ .* Dz.(UzRef1)) + expand_derivatives.(2 .* μ .* Dz.(UzRef1))

# ╔═╡ 860ac07f-391b-4d62-bf9d-a32abea7cd60
σxz_reflect1 = expand_derivatives.(μ .* Dx.(UzRef1) .+ μ .* Dz.(UxRef1))

# ╔═╡ fdba06fe-e596-4122-b2d5-81b67dd78a97
begin
    if (incident_waves == "P")
        UxRef2 = expand_derivatives.(-Dz.(reflect2))
        UzRef2 = expand_derivatives.(Dx.(reflect2))
    else
        UxRef2 = expand_derivatives.(Dx.(reflect2))
        UzRef2 = expand_derivatives.(Dz.(reflect2))
    end

end

# ╔═╡ 0dfd91aa-f9e1-4d93-a471-0bf99e476af5
σzz_reflect2 = expand_derivatives.(λ .* Dx.(UxRef2) .+ λ .* Dz.(UzRef2)) + expand_derivatives.(2 .* μ .* Dz.(UzRef2))

# ╔═╡ 06a568fa-ba94-4c5f-815d-f30c6c8a4260
σzz_z0 = substitute.(σzz_incident .+ σzz_reflect1 .+ σzz_reflect2, z => 0) .~ 0

# ╔═╡ 64164ec7-a3cf-41d0-bcef-e55475eabeea
σzz = simplify(substitute(σzz_z0, Dict([A => 1, t => 0, x => 0])))

# ╔═╡ 47f2439f-28a1-4ad9-8d8d-67a4a30af5ea
σxz_reflect2 = expand_derivatives.(μ .* Dx.(UzRef2) .+ μ .* Dz.(UxRef2))

# ╔═╡ b4030cef-75ef-40d1-b75d-0e5f233a3023
σxz_z0 = substitute.(σxz_incident .+ σxz_reflect1 .+ σxz_reflect2, z => 0) .~ 0

# ╔═╡ 605195de-f4a1-4331-9b82-78e4a061f0a4
σxz = simplify(substitute(σxz_z0, Dict([A => 1, t => 0, x => 0])))

# ╔═╡ a95e361b-d195-4ec0-b2fd-cd9adb79efc5
sol = simplify.(Symbolics.symbolic_linear_solve([σzz, σxz], [A₁, A₂]))

# ╔═╡ 284e4a79-6cfe-4c4a-939b-55fc69611ecb
"""
	rp(θ)

The ray parameter (horizontal slowness) of a wave incident at angle `θ` (radians),
using whichever of `αin`/`βin` matches the current `incident_waves` type.
"""
rp(θ) = (incident_waves == "P") ? sin(θ) / αin : sin(θ) / βin

# ╔═╡ a06affae-47c3-4dfa-a997-ee75b35ab122
"""
	ηz1(θ)

Vertical slowness of the leg sharing the incident wave's type (the direct
reflection: P→P or S→S). Always real for `θ ∈ [0, π/2]`, since a wave cannot be
incident past its own critical angle.
"""
ηz1(θ) = (incident_waves == "P") ? sqrt((inv(αin)^2 - (rp(θ))^2) + 0im) : sqrt((inv(βin)^2 - (rp(θ))^2) + 0im)

# ╔═╡ edb72f1f-6c99-48fb-b32d-fc4fe7d3328e
ηz1(45)

# ╔═╡ 23d62611-aeb5-4ed3-897d-822c0d17781c
"""
	ηz2(θ)

Vertical slowness of the mode-converted leg (P→S or S→P). Complex once `θ`
exceeds this leg's critical angle -- see the sign convention noted where it's
used to build the wavefield, below.
"""
ηz2(θ) = (incident_waves == "S") ? sqrt((inv(αin)^2 - (rp(θ))^2) + 0im) : sqrt((inv(βin)^2 - (rp(θ))^2) + 0im)

# ╔═╡ 4c871d5c-1b62-4b1c-9c1e-37cae19cbe0e
ηz2(45)

# ╔═╡ a089ab5b-4703-4d4d-a7ab-11197b4b907c
begin
    """
        psv_traction(p, η, λ, μ, mode; incident_s=false)

    Return the `(σzz, σxz)` traction coefficients of one P or SV potential.
    It is the direct numerical form of the symbolic displacement/stress
    construction above. `incident_s=true` retains this notebook's established
    incident-S potential convention exactly, so this performance refactor does
    not alter the displayed physics.
    """
    function psv_traction(
        p::Float64,
        η::ComplexF64,
        λ::Float64,
        μ::Float64,
        mode::Symbol;
        incident_s::Bool = false,
    )
        ∂x = -im * p
        ∂z = -im * η

        if mode === :P
            ux, uz = ∂x, ∂z
        elseif incident_s
            ux, uz = -∂z, ∂z
        else
            ux, uz = -∂z, ∂x
        end

        σzz = λ * (∂x * ux + ∂z * uz) + 2μ * ∂z * uz
        σxz = μ * (∂x * uz + ∂z * ux)
        return ComplexF64(σzz), ComplexF64(σxz)
    end

    """
        psv_free_surface_coefficients(α, β, ρ, incident_type, θ)

    Compute the two free-surface reflection amplitudes directly from the
    traction-free 2×2 system. This preserves the coefficient order consumed by
    the widget while avoiding per-angle symbolic substitution and simplification.
    """
    function psv_free_surface_coefficients(
        α::Float64,
        β::Float64,
        ρ::Float64,
        incident_type::String,
        θ::Float64,
    )
        λ = ρ * (α^2 - 2β^2)
        μ = ρ * β^2
        p = sin(θ) / (incident_type == "P" ? α : β)
        ηP = ComplexF64(sqrt((inv(α)^2 - p^2) + 0im))
        ηS = ComplexF64(sqrt((inv(β)^2 - p^2) + 0im))

        if incident_type == "P"
            σzzᵢ, σxzᵢ = psv_traction(p, ηP, λ, μ, :P)
            σzz₁, σxz₁ = psv_traction(p, -ηP, λ, μ, :P)
            σzz₂, σxz₂ = psv_traction(p, -ηS, λ, μ, :S)
        else
            σzzᵢ, σxzᵢ = psv_traction(p, ηS, λ, μ, :S; incident_s = true)
            σzz₁, σxz₁ = psv_traction(p, -ηS, λ, μ, :S)
            σzz₂, σxz₂ = psv_traction(p, -ηP, λ, μ, :P)
        end

        determinant = σzz₁ * σxz₂ - σzz₂ * σxz₁
        A₁ = (-σzzᵢ * σxz₂ + σzz₂ * σxzᵢ) / determinant
        A₂ = (-σzz₁ * σxzᵢ + σzzᵢ * σxz₁) / determinant

        # The existing symbolic implementation reports P then S for S incidence.
        return incident_type == "P" ?
               (A1 = ComplexF64(A₁), A2 = ComplexF64(A₂)) :
               (A1 = ComplexF64(A₂), A2 = ComplexF64(A₁))
    end

    Avec = (
        θ -> psv_free_surface_coefficients(αin, βin, ρin, incident_waves, Float64(θ)).A1,
        θ -> psv_free_surface_coefficients(αin, βin, ρin, incident_waves, Float64(θ)).A2,
    )
end

# ╔═╡ ffff6666-6666-6666-6666-666666666666
let
    A1_normal = abs(Avec[1](0.0))
    A2_normal = abs(Avec[2](0.0))
    @assert A2_normal<1e-8 "expected no mode conversion at normal incidence, got |A2|=$(A2_normal)"
    @assert abs(A1_normal-1)<1e-8 "expected total reflection at normal incidence, got |A1|=$(A1_normal)"
    (direct_reflection_magnitude=A1_normal, mode_converted_magnitude=A2_normal)
end

# ╔═╡ cccc3333-3333-3333-3333-333333333333
# Push angle-swept slowness/coefficient curves to the widget above -- it stays
# mounted across reruns of this cell (its `@bind` cell doesn't depend on medium
# parameters), same CustomEvent pattern geoid-kernel uses to push computed maps
# back to an already-rendered widget. This is the ONLY place the widget's Julia
# side runs after a medium-parameter change: everything angle/frequency/display
# related happens client-side from these arrays.
let
    thetas_deg = 0:0.5:90
    thetas_rad = deg2rad.(thetas_deg)
    p_sweep = rp.(thetas_rad)
    eta1_sweep = ηz1.(thetas_rad)
    eta2_sweep = ηz2.(thetas_rad)
    A1_sweep = Vector{ComplexF64}(undef, length(thetas_rad))
    A2_sweep = Vector{ComplexF64}(undef, length(thetas_rad))
    for i in eachindex(thetas_rad)
        coefficients = psv_free_surface_coefficients(αin, βin, ρin, incident_waves, thetas_rad[i])
        A1_sweep[i] = coefficients.A1
        A2_sweep[i] = coefficients.A2
    end

    num(x) = isfinite(x) ? string(round(Float64(x), digits=6)) : "0"
    jsonarr(v) = "[" * join(num.(v), ",") * "]"

    payload = string(
        "{\"theta_deg\":", jsonarr(thetas_deg),
        ",\"p\":", jsonarr(real.(p_sweep)),
        ",\"eta1_re\":", jsonarr(real.(eta1_sweep)), ",\"eta1_im\":", jsonarr(imag.(eta1_sweep)),
        ",\"eta2_re\":", jsonarr(real.(eta2_sweep)), ",\"eta2_im\":", jsonarr(imag.(eta2_sweep)),
        ",\"A1_re\":", jsonarr(real.(A1_sweep)), ",\"A1_im\":", jsonarr(imag.(A1_sweep)),
        ",\"A2_re\":", jsonarr(real.(A2_sweep)), ",\"A2_im\":", jsonarr(imag.(A2_sweep)),
        ",\"incident_type\":\"", incident_waves, "\"",
        "}",
    )
    HTML("""<script>
      window.dispatchEvent(new CustomEvent('psv-results', {detail: $(repr(payload))}));
    </script>""")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
PlutoUI = "~0.7.83"
Symbolics = "~7.36.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "b306ae83a4f51b18c94e39299c3fce183e167072"

[[deps.ADTypes]]
git-tree-sha1 = "5970c86505ae9c07bf5bc521ef2bbbb3849e8b7b"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.23.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

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

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "daa72978cd7a624246e894a4f4f067706d4e17e2"
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
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "60f11b38ebeabd984f5535838d91e197d97202f0"
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
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ef2022bff55342a8c9846cdf218f62e475f0444d"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.1.2"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

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
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "FunctionMaps", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c0f576ae49bd2d1bc904b9946f4783db8f0ef530"
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
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.6"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExprTools]]
git-tree-sha1 = "d2e49e7efd29719d6f28b891b0e0e159daa9d2b4"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.11"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.FunctionMaps]]
deps = ["CompositeTypes", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "31bd99a57edf98990d1c21486032963955450e8d"
uuid = "a85aefff-f8ca-4649-a888-c8e5398bc76c"
version = "0.1.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7eb45fe833a5b7c51cf6d89c5a841d5967e44be3"
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

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "c72458f1962faeb003bf23cbdb75164fe6280906"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

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

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "bba2d9aa057d8f126415de240573e86a8f39d2a1"
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
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "60beb0717782a3bbe0f7df56decad0ef89048c23"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.12"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics", "StarAlgebras"]
git-tree-sha1 = "4838893d9b035c2f6967c0d533350e1755b58a70"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.19"

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

    [deps.MultivariatePolynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "dc5b2c4c111c46bc79ac4405eeb563523b39c004"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.8.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

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

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
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
uuid = "988b38a3-91fc-5605-94a2-ee2116b3bd83"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "65c9e1142f0372bfc16ba14b9edd57737fe0039f"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.24"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLPublic]]
git-tree-sha1 = "cf9aaf8b9ed5db993259ea8b24cf2b7ba9bd3b79"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.2.4"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "c3ac026e735264e9bdc6a9bcbd1b1e781b36e3bc"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.3"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StarAlgebras]]
deps = ["LinearAlgebra", "MutableArithmetics", "SparseArrays"]
git-tree-sha1 = "235b1f9d287bbf34083b3d0829343a7942c0ad1c"
uuid = "0c0c59c1-dc5f-42e9-9a8b-b5dc384a6cd1"
version = "0.3.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
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
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "2167b9913f3013a1485bdc9bb249123eb8b53cb0"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.54"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "90c1f0a9d1c65e462bdb7180c3ea21e0139e9748"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "1.1.5"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "Graphs", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "3bfccb39a7de6ebb9c834cd46909f6a7a009a967"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.45.0"

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
git-tree-sha1 = "83ae8cd3decb77ab9bcde0ed19f211a879335a3f"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.36.0"

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

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

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

[[deps.WeakCacheSets]]
git-tree-sha1 = "386050ae4353310d8ff9c228f83b1affca2f7f38"
uuid = "d30d5f5c-d141-4870-aa07-aabb0f5fe7d5"
version = "0.1.0"

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
# ╠═08429397-3964-4600-bc14-c45d22c915ec
# ╟─32a757ff-aa7a-41d5-b8b3-6ac9e0125875
# ╟─bbbb2222-2222-2222-2222-222222222222
# ╟─afdb5b7d-d670-4a98-a91d-3ff638fb0294
# ╟─4476bf78-39e3-4674-a152-db19fe80929a
# ╠═927dcc43-202c-4ad2-a76c-837d41f1ed6c
# ╠═c2eacb91-38e8-428e-9247-691950668bc3
# ╟─0bea69f9-3ffa-4695-a5eb-6e962d2a81ce
# ╠═e8729ceb-e85f-4c21-89d2-f15ec69a840f
# ╟─670d5198-f810-472e-8db8-b13385ca294a
# ╠═6043b904-6991-4776-bc1b-13eda3fcc936
# ╟─c7c52926-6e2e-4c4c-a01f-09f57c2ececd
# ╟─33c6aa70-87cb-464a-88f7-b82d17476a2f
# ╠═fbc6cd6c-0b3b-44e9-b2a6-5edb0eeced1b
# ╠═8dad08c2-b9a8-40ec-a8ef-c9383e5ec7b1
# ╠═13b5abc9-10be-4ef1-817d-bcee1271c89e
# ╟─aae7b893-4b33-44ae-b01c-6345a1180d0d
# ╠═f488f9f3-e73b-42ae-b3c2-2262661fd839
# ╠═4dc428b0-f141-4e94-9edb-e847780333aa
# ╠═8940ee90-2dd1-448c-92cd-09fd1ebbaea7
# ╟─a1213dc4-3ea9-4b26-aa82-a15cc3e72401
# ╠═b4a603fb-5ee7-4c97-8e1d-fd0a89692503
# ╠═74484bb8-a9dc-4b54-894e-3bf86c18278e
# ╠═fdba06fe-e596-4122-b2d5-81b67dd78a97
# ╟─52d2ba9b-826a-4805-9583-f1a80b10a926
# ╠═7cf515f8-2496-4bde-b3b3-9d39d971764a
# ╠═040f1440-14c4-491b-b792-570aa3c2080b
# ╠═3525bf50-b0e5-4c3b-ba75-77689a5799f4
# ╠═0dfd91aa-f9e1-4d93-a471-0bf99e476af5
# ╠═f51dfb73-b5ef-4de8-b4c3-fc8814a31485
# ╠═860ac07f-391b-4d62-bf9d-a32abea7cd60
# ╠═47f2439f-28a1-4ad9-8d8d-67a4a30af5ea
# ╠═06a568fa-ba94-4c5f-815d-f30c6c8a4260
# ╠═b4030cef-75ef-40d1-b75d-0e5f233a3023
# ╠═64164ec7-a3cf-41d0-bcef-e55475eabeea
# ╠═605195de-f4a1-4331-9b82-78e4a061f0a4
# ╠═a95e361b-d195-4ec0-b2fd-cd9adb79efc5
# ╟─5729b459-b283-41ce-95be-c4d33a7c28c0
# ╟─281cb873-760c-411f-98b5-1c64d218e7e9
# ╠═3d50985b-b79a-45ab-ba1f-5287936c56d9
# ╠═284e4a79-6cfe-4c4a-939b-55fc69611ecb
# ╟─a7cb9cd5-7013-461b-960d-5da73db62aea
# ╠═a06affae-47c3-4dfa-a997-ee75b35ab122
# ╠═23d62611-aeb5-4ed3-897d-822c0d17781c
# ╠═edb72f1f-6c99-48fb-b32d-fc4fe7d3328e
# ╠═4c871d5c-1b62-4b1c-9c1e-37cae19cbe0e
# ╟─8c81ddb5-bf4d-4610-bfea-3d1a27ffd61f
# ╠═a089ab5b-4703-4d4d-a7ab-11197b4b907c
# ╠═eeee5555-5555-5555-5555-555555555555
# ╠═ffff6666-6666-6666-6666-666666666666
# ╠═cccc3333-3333-3333-3333-333333333333
# ╟─dba7ea14-e0dd-4dc4-ad9c-4627fd16cc62
# ╟─ea3b7089-bda8-4694-8042-98534b1739bd
# ╠═ab40f79c-3d8a-11ed-0697-a7b794dbba99
# ╠═602f13f9-6d14-41fd-9183-b8255d64399b
# ╟─aaaa7777-7777-7777-7777-777777777777
# ╠═dddd4444-4444-4444-4444-444444444444
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
