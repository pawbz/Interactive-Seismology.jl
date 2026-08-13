### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Spherical Harmonics"
#> tags = ["normalmodes"]
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

# ╔═╡ 348240ee-f35f-11ef-01ab-d11599b51ca1
begin
    using PlutoUI, PlutoPlotly, AssociatedLegendrePolynomials, Bessels



	
end

# ╔═╡ 88195cf5-894b-4676-99bc-a559e0d5ebd9
using ForwardDiff

# ╔═╡ ca388aaf-f515-4c8e-8b39-a7173641dca0
TableOfContents()

# ╔═╡ 2a7c4b45-0153-4095-a67d-b48f8f343d6a
md"""
# Visualization of Spherical Harmonics

Spherical harmonics are fundamental mathematical functions used to describe **wave phenomena** on a sphere. They appear in seismology while modeling Earth's free oscillations. 

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
    """

# ╔═╡ 866c00fd-bf69-460b-90b2-79ecf1f72d66
md"""
##### What are `l`, `m`, and the R/S/T motion field?

The **degree `l`** counts how many nodal lines (zero-crossings) the pattern has on the
sphere &mdash; higher `l` means a more finely-lobed pattern. The **order `m`** (with
`-l ≤ m ≤ l`) splits those lobes between latitude and longitude: `m = 0` gives purely
latitudinal bands, `|m| = l` gives purely longitudinal (orange-slice) lobes.

For Earth's free oscillations, a mode's spatial pattern on the sphere is exactly a
spherical harmonic `Y_l^m`, and its particle motion decomposes into three vector fields:
- **Radial (`R`)**: purely outward/inward motion, like a breathing sphere &mdash;
  associated with the "spheroidal" family (also involves horizontal motion in general;
  here isolated to show the radial part alone).
- **Poloidal (`S`)**: horizontal motion aligned with the surface gradient of `Y_l^m`
  &mdash; the tangential part of spheroidal modes.
- **Toroidal (`T`)**: horizontal motion perpendicular to the gradient of `Y_l^m`
  &mdash; a purely rotational shearing motion with no radial component at all, the
  signature of toroidal modes.

The radial panel shows a spherical Bessel function `j_l(ωr/c)` for a **uniform,
non-dispersive sphere** (`c = 5 km/s` throughout, not a real depth-dependent Earth
model) &mdash; a simplification that keeps the radial shape simple enough to read at a
glance, in exchange for not being quantitatively Earth-like.
"""

# ╔═╡ 4d6aed73-9b2f-47c6-a161-38e2633cb8db
md"---"

# ╔═╡ 6750df2d-459b-4467-b405-1bd8177bee42
md"## Appendix"

# ╔═╡ fff4b588-59c3-4cd5-b333-222026b8d666
c = 5.0 # in km/s

# ╔═╡ 4c892f32-ffa3-4b45-a840-6f1a165293b3
md"### Surface Spherical Harmonics "

# ╔═╡ 186c6279-4573-4fda-8589-5fabcbc1cc0e
"""
    Y_lm(l, m, θ, φ)

The complex surface spherical harmonic `Y_l^m(θ,φ)`, fully normalized so that
`∫|Y_l^m|² dΩ = 1` over the unit sphere. `θ` is colatitude (0 at the north pole),
`φ` is longitude. Built from `AssociatedLegendrePolynomials.λlm` (which already
carries the normalization and Condon-Shortley phase) times the azimuthal factor
`exp(im*m*φ)`.
"""
function Y_lm(l, m, θ, φ)
	P_lm = λlm(l, abs(m), cos(θ))
    return P_lm * exp(im * m * φ)
end

# ╔═╡ b5be41c0-ff78-4f5d-a795-4c5a79135fb6
"""
    Y_lm(l, m, thetaphi)

Vector-argument form of [`Y_lm`](@ref) -- `thetaphi = [θ, φ]`. Exists so
`ForwardDiff.gradient` can differentiate `Y_lm` with respect to both angles at once
(see [`spherical_harmonic_gradient`](@ref)).
"""
function Y_lm(l, m, thetaphi)
	return Y_lm(l, m, thetaphi[1], thetaphi[2])
end

# ╔═╡ 7eb05e0a-12e1-4033-b081-1fc79b38218a
"""
    spherical_harmonics(l, m, θ, φ)

`Y_l^m(θ,φ)`, or `zero(θ)` when `|m| > l` (an invalid but commonly-swept-over
combination when `m` is driven by a slider that hasn't yet clamped to the new `l`).
Returns the full complex value -- callers take `real(...)` where a scalar is needed.
"""
function spherical_harmonics(l, m, θ, φ)
        if abs(m) > l
            return zero(θ)
        end
        Y = Y_lm(l, m, θ, φ)
        return Y  # complex value; callers take real(...) where a scalar is needed
end

# ╔═╡ 79fa4c7b-5f93-44f9-bf33-311f5efecf32
md"### Vector Spherical Harmonics"

# ╔═╡ ea0d2403-03ea-4c42-a452-59f830c22a16
"""
    spherical_harmonic_gradient(l::Int, m::Int, θ, φ)

The surface gradient `∇Y_l^m = (∂Y/∂θ, ∂Y/∂φ)`, computed with `ForwardDiff` by
differentiating the real and imaginary parts of [`Y_lm`](@ref) separately and
recombining as a complex vector. Feeds the poloidal (`S`) and toroidal (`T`)
components in [`vector_spherical_harmonics`](@ref).
"""
function spherical_harmonic_gradient(l::Int, m::Int, θ, φ)
    if abs(m) > l
        return (zero(θ), zero(θ))  # Edge case where |m| > l
    end
	Y = Y_lm(l, m, θ, φ)

	dY = ForwardDiff.gradient(x->real(Y_lm(l, m, x)), [θ,φ]) + im * ForwardDiff.gradient(x->imag(Y_lm(l, m, x)), [θ,φ])

    return dY
end

# ╔═╡ 38e90f9a-5cba-451f-bae3-341f03988e3c
"""
    vector_spherical_harmonics(l, m, θ, φ, mode_type)

One of the three vector spherical harmonics at `(θ,φ)`, returned as spherical
components `[r, θ, φ]`:
- `"R"` (radial): purely outward/inward motion, magnitude `Y_l^m` itself.
- `"S"` (poloidal): tangential motion along `∇Y_l^m`, normalized by `1/√(l(l+1))`.
- `"T"` (toroidal): tangential motion perpendicular to `∇Y_l^m` (a 90° rotation of
  the poloidal field), same normalization.
"""
function vector_spherical_harmonics(l, m, θ, φ, mode_type)
        if abs(m) > l
            return zero(θ)
        end
		Y = Y_lm(l, m, θ, φ)
		dY = spherical_harmonic_gradient(l, m, θ, φ)
		s = inv(sqrt(l*(l+1)))

 	if mode_type == "R"
        return [Y, 0, 0]  # Purely outward motion (gravity-like)
    elseif mode_type == "S"
        return s.*[0, dY[1], dY[2] / sin(θ)]  # Tangential movement
    elseif mode_type == "T"
        return s.*[0, dY[2] / sin(θ), - dY[1]] # Rotational motion
    else
        return [0, 0, 0]
    end
end

# ╔═╡ 433d799f-9243-4140-ab1f-006feeaf8816
md"### Grids"

# ╔═╡ 7100ad22-5762-465c-a6a4-bbb06e11ac2f
begin
	    θ_ε = 1e-3  # keep grid off the poles: dY[2]/sin(θ) blows up exactly at θ=0,π
	    θ_range = range(θ_ε, π - θ_ε, length=50)  # Latitude
	    φ_range = range(0, 2π, length=100)  # Longitude
	    θ_grid = first.(Iterators.product(θ_range, φ_range))
		φ_grid = last.(Iterators.product(θ_range, φ_range))
end;

# ╔═╡ 71c529bf-77d8-481c-94d3-86a80dfcd641
radius_max = 6000 # in km

# ╔═╡ 8ed29e6e-583c-4fb3-86d3-f34b8d025c0f
radius_grid = range(0., radius_max, length=1000);

# ╔═╡ ee3ca7fd-327b-408f-918f-55a387e9cb7c
md"### The Interactive Widget"

# ╔═╡ e2bb7ecd-bcad-4b61-9eae-020c2a507435
begin
    """
        SphericalHarmonicInput(; l=3, m=0, omega=0.05, mode="S")

    Bound widget: degree `l` (0-20), order `m` (auto-clamped client-side to `[-l,l]`),
    radial frequency `omega`, and R/S/T motion-field `mode`. Renders a drag-to-rotate 3D
    sphere (colored and radius-lobed by `real(Y_l^m)`, optionally overlaid with motion
    vectors) plus a radial Bessel-function profile. Bound value is a
    `Dict{String,Any}` with keys `"l"`, `"m"`, `"omega"`, `"mode"` (the same shape
    `Base.get` returns below, matching Pluto's default JS-object bond transport).
    """
    struct SphericalHarmonicInput
        l::Int
        m::Int
        omega::Float64
        mode::String
    end

    SphericalHarmonicInput(; l=3, m=0, omega=0.05, mode="S") =
        SphericalHarmonicInput(Int(l), Int(m), Float64(omega), mode)

    Base.get(w::SphericalHarmonicInput) = Dict{String,Any}(
        "l" => w.l, "m" => w.m, "omega" => w.omega, "mode" => w.mode,
    )

    function Base.show(io::IO, ::MIME"text/html", w::SphericalHarmonicInput)
        write(io, """
<div id="shwwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#shwwidget){width:min(80vw,1400px)!important;margin-left:calc((100% - min(80vw,1400px))/2)!important}
    #shwwidget .shw-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #shwwidget .shw-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #shwwidget .shw-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #shwwidget .shw-workspace{display:flex;gap:12px;align-items:flex-start;justify-content:center;width:100%;flex-wrap:wrap}
    #shwwidget canvas{background:#000;border:1px solid #374151;border-radius:6px;display:block}
    #shwwidget #shwsphere{cursor:grab}
    #shwwidget #shwsphere.dragging{cursor:grabbing}
    #shwwidget .shw-controls{width:min(var(--totalw,960px),100%);margin-top:8px;display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:8px;font:14px sans-serif}
    #shwwidget .shw-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
    #shwwidget .shw-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #shwwidget .shw-control-row{display:grid;grid-template-columns:minmax(60px,90px) minmax(70px,1fr) minmax(40px,64px);gap:6px;align-items:center;margin:7px 0}
    #shwwidget .shw-control-row input[type=range]{width:100%;min-width:0;vertical-align:middle}
    #shwwidget .shw-value{color:#d1d5db;text-align:left;font-variant-numeric:tabular-nums}
    #shwwidget .shw-actions{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
    #shwwidget label{color:#d1d5db}
    #shwwidget select{background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:3px 6px}
    @media (max-width: 900px){
      #shwwidget .shw-workspace{flex-direction:column;align-items:center}
      #shwwidget .shw-controls{grid-template-columns:1fr;width:660px;max-width:100%}
    }
  </style>
  <div class="shw-title">
    <div class="shw-title-desc">Degree `l`, order `m` shape the pattern; R/S/T picks which motion field is drawn on top of it.</div>
    <div class="shw-title-hint">drag the sphere to rotate &middot; slide l/m/&omega; to change the mode &middot; toggle motion vectors on/off</div>
  </div>
  <div class="shw-workspace">
    <canvas id="shwsphere" width="560" height="560"></canvas>
    <canvas id="shwradial" width="260" height="560"></canvas>
  </div>
  <div class="shw-controls">
    <div class="shw-control-group">
      <div class="shw-control-title">Spherical Harmonic</div>
      <label class="shw-control-row"><span>degree l</span><input type="range" id="shw-l" min="0" max="20" step="1" value="$(w.l)"><span id="shw-l-v" class="shw-value">$(w.l)</span></label>
      <label class="shw-control-row"><span>order m</span><input type="range" id="shw-m" min="$(-w.l)" max="$(w.l)" step="1" value="$(w.m)"><span id="shw-m-v" class="shw-value">$(w.m)</span></label>
    </div>
    <div class="shw-control-group">
      <div class="shw-control-title">Radial Frequency</div>
      <label class="shw-control-row"><span>&omega;</span><input type="range" id="shw-omega" min="0" max="0.1" step="0.0001" value="$(w.omega)"><span id="shw-omega-v" class="shw-value">$(w.omega)</span></label>
      <div style="font-size:12px;color:#6b7280;margin-top:4px">uniform sphere, c = 5 km/s</div>
    </div>
    <div class="shw-control-group">
      <div class="shw-control-title">Motion Field</div>
      <div class="shw-actions">
        <select id="shw-mode">
          <option value="R" $(w.mode == "R" ? "selected" : "")>Radial</option>
          <option value="S" $(w.mode == "S" ? "selected" : "")>Poloidal</option>
          <option value="T" $(w.mode == "T" ? "selected" : "")>Toroidal</option>
        </select>
        <label><input type="checkbox" id="shw-vec" checked> show vectors</label>
      </div>
    </div>
    <div class="shw-control-group">
      <div class="shw-control-title">Legend</div>
      <div style="font-size:13px;line-height:1.6">
        <span style="color:#ef4444">&#9632;</span> Y &gt; 0&emsp;<span style="color:#3b82f6">&#9632;</span> Y &lt; 0<br>
        <span style="color:#f5f3ef">&#8594;</span> motion vector (current R/S/T field)
      </div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1400)
  const totalW = Math.max(700, availW)
  par.style.setProperty('--totalw', Math.round(totalW)+'px')
  const SEC = Math.round(Math.min(totalW*0.62, 560))
  const RADW = Math.round(Math.min(totalW*0.30, 260))
  const DPR = Math.min(window.devicePixelRatio || 1, 2)

  let l = $(w.l), m = $(w.m), omega = $(w.omega), mode = "$(w.mode)"
  let showVectors = true
  let yaw = 0.6, pitch = 0.3, dragging = false, lastX = 0, lastY = 0
  let data = null   // filled in by the 'shw-results' push from Julia, below

  const sphCvs = par.querySelector('#shwsphere'), sctx = sphCvs.getContext('2d')
  const radCvs = par.querySelector('#shwradial'), rctx = radCvs.getContext('2d')
  function hidpi(cv, cx, w, h){ cv.width=Math.round(w*DPR); cv.height=Math.round(h*DPR); cv.style.width=w+'px'; cv.style.height=h+'px'; cx.setTransform(DPR,0,0,DPR,0,0) }
  hidpi(sphCvs, sctx, SEC, SEC)
  hidpi(radCvs, rctx, RADW, SEC)

  const CX = SEC/2, CY = SEC/2, RPIX = SEC*0.36

  function rot(p){
    const x = p[0]*Math.cos(yaw) - p[1]*Math.sin(yaw)
    const y = p[0]*Math.sin(yaw) + p[1]*Math.cos(yaw)
    const z = p[2]
    return [x, y*Math.cos(pitch) - z*Math.sin(pitch), y*Math.sin(pitch) + z*Math.cos(pitch)]
  }
  function proj(p){ return [CX + p[0]*RPIX, CY - p[2]*RPIX] }

  // Diverging ramp on a dark background: positive Y -> red, negative Y -> blue
  // (same construction as tomoColor in earth-internal-structure.jl's globe, opposite
  // sign convention since here it's a mathematical field, not a physical anomaly).
  function ylmColor(v, mx){
    const t = Math.max(-1, Math.min(1, mx > 0 ? v/mx : 0))
    const bg = [18, 22, 29]
    const hi = t >= 0 ? [239, 68, 68] : [59, 130, 246]
    const a = Math.min(1, Math.abs(t) * 1.6)
    const c = [0,1,2].map(i => Math.round(bg[i] + (hi[i]-bg[i])*a))
    return 'rgb('+c[0]+','+c[1]+','+c[2]+')'
  }

  function idx(i, j, nt){ return j*nt + i }   // 0-based, matches Julia's column-major (theta-fastest) flattening

  function drawSphere(){
    sctx.clearRect(0,0,SEC,SEC)
    if(!data){
      sctx.fillStyle = '#6b7280'; sctx.font = '13px sans-serif'
      sctx.fillText('computing...', 12, 20)
      return
    }
    const {nt, np, x, y, z, ylm} = data
    let ymax = 1e-9
    for(const v of ylm) ymax = Math.max(ymax, Math.abs(v))

    for(let i=0;i<nt-1;i++){
      for(let j=0;j<np;j++){
        const j1 = (j+1) % np
        const c00 = idx(i,j,nt), c01 = idx(i,j1,nt), c11 = idx(i+1,j1,nt), c10 = idx(i+1,j,nt)
        const pts = [c00,c01,c11,c10].map(k => rot([x[k],y[k],z[k]]))
        const avgDepth = (pts[0][1]+pts[1][1]+pts[2][1]+pts[3][1])/4
        if(avgDepth < 0) continue
        const vAvg = (ylm[c00]+ylm[c01]+ylm[c11]+ylm[c10])/4
        sctx.beginPath()
        const s0 = proj(pts[0]); sctx.moveTo(s0[0], s0[1])
        for(let k=1;k<4;k++){ const s = proj(pts[k]); sctx.lineTo(s[0], s[1]) }
        sctx.closePath()
        sctx.fillStyle = ylmColor(vAvg, ymax)
        sctx.fill()
      }
    }

    if(showVectors && data.arrow_idx && data.arrow_idx.length){
      let umax = 1e-9
      for(let k=0;k<data.arrow_ux.length;k++) umax = Math.max(umax, Math.hypot(data.arrow_ux[k], data.arrow_uy[k], data.arrow_uz[k]))
      const maxPix = 26
      sctx.strokeStyle = '#f5f3ef'; sctx.fillStyle = '#f5f3ef'; sctx.lineWidth = 1.6
      for(let a=0;a<data.arrow_idx.length;a++){
        const k = data.arrow_idx[a]
        const base = rot([x[k], y[k], z[k]])
        if(base[1] < 0) continue
        const mag = Math.hypot(data.arrow_ux[a], data.arrow_uy[a], data.arrow_uz[a])
        if(mag < 1e-9) continue
        const scale = (mag/umax) * maxPix / RPIX
        const dir = rot([data.arrow_ux[a], data.arrow_uy[a], data.arrow_uz[a]])
        const tip = [base[0]+dir[0]*scale, base[1]+dir[1]*scale, base[2]+dir[2]*scale]
        const s0 = proj(base), s1 = proj(tip)
        sctx.beginPath(); sctx.moveTo(s0[0], s0[1]); sctx.lineTo(s1[0], s1[1]); sctx.stroke()
        const ang = Math.atan2(s1[1]-s0[1], s1[0]-s0[0])
        sctx.beginPath()
        sctx.moveTo(s1[0], s1[1])
        sctx.lineTo(s1[0]-4*Math.cos(ang-0.4), s1[1]-4*Math.sin(ang-0.4))
        sctx.lineTo(s1[0]-4*Math.cos(ang+0.4), s1[1]-4*Math.sin(ang+0.4))
        sctx.closePath(); sctx.fill()
      }
    }

    sctx.beginPath(); sctx.arc(CX,CY,RPIX,0,2*Math.PI)
    sctx.strokeStyle='#4b5563'; sctx.lineWidth=1; sctx.stroke()
    sctx.fillStyle='#9ca3af'; sctx.font='14px sans-serif'
    sctx.fillText('Y_'+l+'^'+m+' colored on the sphere, radius lobed by |Y|', 12, 20)
    sctx.font='13px sans-serif'
    sctx.fillText('drag to rotate', 12, 38)
  }

  function drawRadial(){
    rctx.clearRect(0,0,RADW,SEC)
    if(!data){
      rctx.fillStyle = '#6b7280'; rctx.font = '12px sans-serif'
      rctx.fillText('computing...', 10, 18)
      return
    }
    const PAD_L = 42, PAD_R = 10, PAD_T = 30, PAD_B = 24
    const x0 = PAD_L, x1 = RADW-PAD_R, y0 = PAD_T, y1 = SEC-PAD_B
    const rMax = Math.max(...data.radial.map(Math.abs)) * 1.1 || 1
    const depthMax = data.radius_grid[data.radius_grid.length-1]
    const X = v => x0 + (x1-x0)*0.5 + (v/rMax)*(x1-x0)*0.5
    const Y = d => y0 + (d/depthMax)*(y1-y0)

    rctx.strokeStyle = '#1f2937'; rctx.fillStyle = '#9ca3af'; rctx.font='10px sans-serif'; rctx.textAlign='right'
    for(let i=0;i<=5;i++){
      const d = depthMax*i/5, y = Y(d)
      rctx.beginPath(); rctx.moveTo(x0,y); rctx.lineTo(x1,y); rctx.stroke()
      rctx.fillText(Math.round(d)+'', x0-6, y+3)
    }
    rctx.strokeStyle = '#374151'
    rctx.beginPath(); rctx.moveTo(X(0),y0); rctx.lineTo(X(0),y1); rctx.stroke()

    rctx.strokeStyle = '#38bdf8'; rctx.lineWidth = 2
    rctx.beginPath()
    data.radius_grid.forEach((r,i) => { const px=X(data.radial[i]), py=Y(r); i===0?rctx.moveTo(px,py):rctx.lineTo(px,py) })
    rctx.stroke()

    rctx.textAlign='center'; rctx.fillStyle='#e5e7eb'; rctx.font='12px sans-serif'
    rctx.fillText('Radial function', (x0+x1)/2, 16)
    rctx.save(); rctx.translate(12, (y0+y1)/2); rctx.rotate(-Math.PI/2)
    rctx.fillStyle='#9ca3af'; rctx.font='11px sans-serif'; rctx.textAlign='center'
    rctx.fillText('radius (km)', 0, 0)
    rctx.restore()
    rctx.fillStyle='#9ca3af'; rctx.font='11px sans-serif'
    rctx.fillText('period '+data.period_min+' min', 8, SEC-6)
  }

  function redraw(){ drawSphere(); drawRadial() }
  redraw()

  function emit(){
    par.value = {l, m, omega, mode}
    par.dispatchEvent(new CustomEvent('input'))
  }

  const lSlider = par.querySelector('#shw-l'), mSlider = par.querySelector('#shw-m')
  lSlider.addEventListener('input', e=>{
    l = parseInt(e.target.value)
    par.querySelector('#shw-l-v').textContent = l
    mSlider.min = -l; mSlider.max = l
    if(Math.abs(m) > l){ m = 0; mSlider.value = 0; par.querySelector('#shw-m-v').textContent = 0 }
    emit()
  })
  mSlider.addEventListener('input', e=>{
    m = parseInt(e.target.value)
    par.querySelector('#shw-m-v').textContent = m
    emit()
  })
  par.querySelector('#shw-omega').addEventListener('input', e=>{
    omega = parseFloat(e.target.value)
    par.querySelector('#shw-omega-v').textContent = omega.toFixed(4)
    emit()
  })
  par.querySelector('#shw-mode').addEventListener('change', e=>{
    mode = e.target.value
    emit()
  })
  par.querySelector('#shw-vec').addEventListener('change', e=>{
    showVectors = e.target.checked
    redraw()
  })

  sphCvs.addEventListener('mousedown', e=>{ dragging=true; lastX=e.offsetX; lastY=e.offsetY; sphCvs.classList.add('dragging') })
  sphCvs.addEventListener('mousemove', e=>{
    if(dragging){
      const dx=e.offsetX-lastX, dy=e.offsetY-lastY
      lastX=e.offsetX; lastY=e.offsetY
      yaw += dx*0.01
      pitch = Math.max(-1.3, Math.min(1.3, pitch + dy*0.01))
      drawSphere()
    }
  })
  window.addEventListener('mouseup', ()=>{ dragging=false; sphCvs.classList.remove('dragging') })

  window.addEventListener('shw-results', ev => {
    data = JSON.parse(ev.detail)
    redraw()
  })
</script>
""")
    end

    # Forces the correct execution order on a fresh restart (Pluto's static dependency
    # analysis doesn't reliably detect that the bind cell below depends on this cell
    # defining SphericalHarmonicInput) -- same pattern as `_epi_ready`/`_tgi_ready` in
    # earth-internal-structure.jl.
    const _shwi_ready = true
end

# ╔═╡ 5b7a1c40-8e91-4a2a-9a2c-2f1c0f6a9a01
begin
    _shwi_ready
    @bind shw SphericalHarmonicInput()
end

# ╔═╡ 2c9e8a11-4f2b-4a7e-9a1a-5b6c7d8e9f01
begin
    # `shw` is `nothing` until the widget's JS has fired its first `emit()` (Pluto has no
    # JS-side initial-value hook wired up for this bond) -- fall back to the same
    # defaults as `SphericalHarmonicInput()`'s keyword constructor until then.
    l_value = shw isa AbstractDict ? round(Int, shw["l"]) : 3
    m_value = shw isa AbstractDict ? round(Int, shw["m"]) : 0
    ω_value = shw isa AbstractDict ? shw["omega"] : 0.05
    mode = shw isa AbstractDict ? shw["mode"] : "S"
end;

# ╔═╡ ca84a03e-053e-47c3-aefc-e81202f0f9f3
spherical_harmonic_gradient(l_value, m_value, 1.2, -0.1)

# ╔═╡ 665f39ad-0764-4166-8872-8179eba23374
# Compute spherical harmonics
Ylm_grid = [real(spherical_harmonics(l_value, m_value, θ, φ)) for (θ, φ) in zip(θ_grid, φ_grid)];

# ╔═╡ 32fa5f32-72b4-4ce2-8d6b-0ab12acff038
begin
	 # Convert spherical to Cartesian coordinates, lobing the radius by Ylm_grid so
	 # the classic "bulging sphere" shape of a spherical harmonic is visible, not just its color
	lobe_amplitude = 0.35
	lobe_radius = 1 .+ lobe_amplitude .* Ylm_grid ./ maximum(abs, Ylm_grid)
	X = lobe_radius .* sin.(θ_grid) .* cos.(φ_grid)
	Y = lobe_radius .* sin.(θ_grid) .* sin.(φ_grid)
	Z = lobe_radius .* cos.(θ_grid)
end;

# ╔═╡ eaa2c2dd-f109-4a2b-b8b8-4f15d5576ca8
begin
	    # Compute selected VSH function in spherical coordinates
    U_spherical = [vector_spherical_harmonics(l_value, m_value, θ, φ, mode) for (θ, φ) in zip(θ_grid, φ_grid)]
    U_r = [u[1] for u in U_spherical]  # Radial component
    U_θ = [u[2] for u in U_spherical]  # Theta component
    U_φ = [u[3] for u in U_spherical]  # Phi component


    Ux = U_r .* sin.(θ_grid) .* cos.(φ_grid) .+ U_θ .* cos.(θ_grid) .* cos.(φ_grid) .+ U_φ .* -sin.(φ_grid)
    Uy = U_r .* sin.(θ_grid) .* sin.(φ_grid) + U_θ .* cos.(θ_grid) .* sin.(φ_grid) + U_φ .* cos.(φ_grid)
    Uz = U_r .* cos.(θ_grid) + U_θ .* (-sin.(θ_grid)) + U_φ .* 0

end;

# ╔═╡ 4d8bfdc4-175d-4f7a-ad65-338726ad343e
R = Bessels.sphericalbesselj.(l_value, ω_value .* radius_grid ./c);

# ╔═╡ 785b9164-9a3b-4abe-90ed-1cdf7fae7252
# Push the currently-computed sphere mesh, subsampled motion-vector field, and radial
# profile to the widget above -- same CustomEvent pattern used throughout this repo
# (e.g. geoid-kernel.jl, global-seismic-arrivals.jl) to feed an already-rendered widget
# from a downstream, non-@bind cell.
let
    num(x) = isfinite(x) ? string(round(Float64(x), digits=5)) : "0"
    jsonarr(v) = "[" * join(num.(v), ",") * "]"

    nt, np = length(θ_range), length(φ_range)
    # Sparse subgrid for motion-vector glyphs -- the full nt*np mesh is far too dense
    # to draw as individual arrows without turning into visual noise.
    arrow_i = 1:4:nt
    arrow_j = 1:6:np
    arrow_idx = vec([(j - 1) * nt + i for i in arrow_i, j in arrow_j])

    payload = string(
        "{\"nt\":", nt, ",\"np\":", np,
        ",\"x\":", jsonarr(X), ",\"y\":", jsonarr(Y), ",\"z\":", jsonarr(Z),
        ",\"ylm\":", jsonarr(Ylm_grid),
        ",\"arrow_idx\":", jsonarr(arrow_idx .- 1),
        ",\"arrow_ux\":", jsonarr(real.(Ux[arrow_idx])),
        ",\"arrow_uy\":", jsonarr(real.(Uy[arrow_idx])),
        ",\"arrow_uz\":", jsonarr(real.(Uz[arrow_idx])),
        ",\"radius_grid\":", jsonarr(radius_grid),
        ",\"radial\":", jsonarr(R),
        ",\"period_min\":", num(2 * pi / ω_value / 60),
        ",\"l\":", l_value, ",\"m\":", m_value, ",\"mode\":\"", mode, "\"",
        "}",
    )
    HTML("""<script>
      window.dispatchEvent(new CustomEvent('shw-results', {detail: $(repr(payload))}));
    </script>""")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AssociatedLegendrePolynomials = "2119f1ac-fb78-50f5-8cc0-dda848ebdb19"
Bessels = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
AssociatedLegendrePolynomials = "~1.0.2"
Bessels = "~0.2.8"
ForwardDiff = "~1.4.3"
PlutoPlotly = "~0.6.6"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "eddc843fa371db1223676d21a5c2050a3affbb73"

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

[[deps.AssociatedLegendrePolynomials]]
git-tree-sha1 = "c5b6a5ac656586d038dd04441b6e165a21c80f09"
uuid = "2119f1ac-fb78-50f5-8cc0-dda848ebdb19"
version = "1.0.2"

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
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "73d5084cae45f9d0857776ad78cf303fec09eb02"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.4.3"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c89d196f5ffb64bfbf80985b699ea913b0d2c211"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.6.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

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

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

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

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "32a4e09c5f29402573d673901778a0e03b0807b9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.6"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "97c8329c5f503d2936fb36719fe25b9f94b1ae8a"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

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
git-tree-sha1 = "82bee338d650aa515f31866c460cb7e3bcef90b8"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.2"

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
git-tree-sha1 = "3b0738bd7c5645641845da25cbd99800b8718689"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.2"

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

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═ca388aaf-f515-4c8e-8b39-a7173641dca0
# ╟─2a7c4b45-0153-4095-a67d-b48f8f343d6a
# ╟─5b7a1c40-8e91-4a2a-9a2c-2f1c0f6a9a01
# ╟─2c9e8a11-4f2b-4a7e-9a1a-5b6c7d8e9f01
# ╟─866c00fd-bf69-460b-90b2-79ecf1f72d66
# ╟─4d6aed73-9b2f-47c6-a161-38e2633cb8db
# ╟─6750df2d-459b-4467-b405-1bd8177bee42
# ╠═348240ee-f35f-11ef-01ab-d11599b51ca1
# ╠═88195cf5-894b-4676-99bc-a559e0d5ebd9
# ╠═fff4b588-59c3-4cd5-b333-222026b8d666
# ╟─4c892f32-ffa3-4b45-a840-6f1a165293b3
# ╠═186c6279-4573-4fda-8589-5fabcbc1cc0e
# ╠═b5be41c0-ff78-4f5d-a795-4c5a79135fb6
# ╠═7eb05e0a-12e1-4033-b081-1fc79b38218a
# ╟─79fa4c7b-5f93-44f9-bf33-311f5efecf32
# ╠═ea0d2403-03ea-4c42-a452-59f830c22a16
# ╠═ca84a03e-053e-47c3-aefc-e81202f0f9f3
# ╠═38e90f9a-5cba-451f-bae3-341f03988e3c
# ╟─433d799f-9243-4140-ab1f-006feeaf8816
# ╠═7100ad22-5762-465c-a6a4-bbb06e11ac2f
# ╠═665f39ad-0764-4166-8872-8179eba23374
# ╠═32fa5f32-72b4-4ce2-8d6b-0ab12acff038
# ╠═eaa2c2dd-f109-4a2b-b8b8-4f15d5576ca8
# ╠═71c529bf-77d8-481c-94d3-86a80dfcd641
# ╠═8ed29e6e-583c-4fb3-86d3-f34b8d025c0f
# ╠═4d8bfdc4-175d-4f7a-ad65-338726ad343e
# ╟─785b9164-9a3b-4abe-90ed-1cdf7fae7252
# ╟─ee3ca7fd-327b-408f-918f-55a387e9cb7c
# ╟─e2bb7ecd-bcad-4b61-9eae-020c2a507435
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
