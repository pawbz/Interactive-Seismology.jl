### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> order = "4"
#> title = "Stress Tensor"
#> tags = ["introduction"]
#> layout = "layout.jlhtml"
#> description = "Traction, principal stresses, and the Mohr-Coulomb failure criterion"

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

# ╔═╡ 13a3429e-12f6-11ed-326f-c154f5debceb
begin
	using LinearAlgebra
	using PlutoUI
	using Symbolics
	using SymbolicUtils
	using Einsum

	# Forces TableOfContents() below to depend on this cell -- a bare call to a
	# PlutoUI-exported function carries no argument for Pluto's static analysis to detect,
	# so on a fresh restart it can run before `using PlutoUI` has (same class of ordering
	# bug as the @bind-constructor case documented for the widget below).
	const _pkgs_ready = true
end

# ╔═╡ 3f1bc2d1-2327-48cc-984a-df09c936da87
begin
	_pkgs_ready
	TableOfContents()
end

# ╔═╡ fd14256d-8c69-4d8a-91b5-924a32479866
md"""
# Stress Tensor

A fault is just a plane cutting through rock. Whether it slips depends entirely on the
**stress acting on that one plane** — not the full 3×3 stress tensor ``\sigma`` by itself,
but what you get when you project ``\sigma`` onto a particular orientation: a **normal**
stress that clamps the fault shut (or pries it open) and a **shear** stress that drives it to
slip. This notebook is built around that one question — *for a given* ``\sigma``*, what is
the traction on a plane with normal* ``\mathbf{n}``*, and when is it enough to cause failure?*

!!! note "Sign convention"
	Stress here is **tension-positive** (the elasticity/seismology convention, matching the
	``s_{ij}`` sliders in the widget below): a positive normal stress pulls the plane apart, a
	negative one compresses it. The **Mohr diagram** panel of the widget flips this and plots
	**compression as positive** instead — the standard convention in rock and fault
	mechanics, where the interesting physics (friction, failure) only happens under
	compression. Watch for the flipped sign once here and the rest follows.

Drag the white handle on the sphere below to choose a plane; the widget decomposes the
traction on that plane live, plots it on a classical **Mohr circle for stress**, and overlays
a **Mohr-Coulomb failure criterion** with an adjustable **pore pressure** — the mechanism
behind injection-induced seismicity, where pumping fluid into the ground can trigger an
earthquake on a fault that was already close to failure, with no change in the tectonic
stress at all.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ b1a10001-0000-4000-8000-000000000001
md"""
## Cauchy's Traction Theorem

For any plane through a point, with unit normal ``\mathbf{n}``, the **traction** — force per
unit area exerted across that plane — is exactly the stress tensor contracted with the
normal:

```math
\mathbf{t} = \sigma \mathbf{n}, \qquad t_i = \sigma_{ij} n_j.
```

This is Cauchy's stress theorem: it says the traction on *any* orientation is completely
determined by the same six independent components of ``\sigma`` — you never need more than
the tensor itself, no matter which plane you ask about. Splitting ``\mathbf{t}`` into the
part along ``\mathbf{n}`` and the part perpendicular to it gives exactly the two quantities
that matter for fault mechanics:

```math
\sigma_n = \mathbf{n}\cdot\mathbf{t} = n_i\sigma_{ij}n_j, \qquad
\boldsymbol\tau = \mathbf{t} - \sigma_n\mathbf{n}.
```

``\sigma_n`` (a scalar) is the **normal stress** — positive means the plane is being pulled
apart, negative means squeezed shut. ``\boldsymbol\tau`` (a vector, tangent to the plane) is
the **shear traction** — the part that actually drives slip. The widget above computes
exactly these two quantities live, for whichever plane you're currently dragging: the red/blue
arrow is ``\sigma_n\mathbf{n}``, the orange arrow is ``\boldsymbol\tau``.
"""

# ╔═╡ 50eb4916-b3d7-4f2a-90c6-8661cbbd8e7a
@variables σ[1:3, 1:3]

# ╔═╡ 88b33835-d993-4396-8605-bb3456200eb1
@variables n[1:3]

# ╔═╡ b1a10001-0000-4000-8000-000000000002
@einsum t[i] := σ[i, j] * n[j]

# ╔═╡ b1a10001-0000-4000-8000-000000000003
σₙ = sum(n .* t)

# ╔═╡ b1a10001-0000-4000-8000-000000000004
τvec = Symbolics.scalarize(t .- σₙ .* n)

# ╔═╡ b1a10001-0000-4000-8000-000000000005
md"""
## Principal Stresses

Because ``\sigma`` is symmetric, it has three real eigenvalues
``\sigma_1\ge\sigma_2\ge\sigma_3`` — the **principal stresses** — with an orthonormal set of
eigenvectors, the **principal axes**. On a plane oriented exactly along a principal axis, the
shear traction vanishes identically: the traction there is *purely* normal, equal to the
corresponding eigenvalue. Every other orientation mixes some shear in. This is exactly what
the red/blue triad drawn on the sphere in the widget shows (red = tensile, blue = compressive
principal stress) — the three directions where dragging the handle would make the orange
shear arrow disappear entirely.

As a worked numeric example (independent of whatever stress state you've dragged the widget
to — this uses the notebook's original default ``\sigma``):
"""

# ╔═╡ ca99052e-ff24-4e8a-9d06-74bc24974873
σ_example = [1 0.5 0; 0.5 0 0; 0 0 0.5]

# ╔═╡ 987d9253-2f7d-4ca8-8eef-aae8e8d53bc6
eig_example = eigen(σ_example)

# ╔═╡ b1a10001-0000-4000-8000-000000000006
md"""
For this example, the principal stresses come out to ``\sigma_1,\sigma_2,\sigma_3 =
$(join(round.(reverse(eig_example.values), digits=3), ", "))`` — two positive (tensile) and
one negative (compressive) — with the corresponding eigenvectors (columns of
`eig_example.vectors`) giving the three principal directions.
"""

# ╔═╡ b1a10001-0000-4000-8000-000000000007
md"""
## Mean and Deviatoric Stress

Just like strain, the stress tensor splits into an **isotropic** part and a **deviatoric**
remainder. The mean (hydrostatic) stress

```math
p = \tfrac13\operatorname{tr}(\sigma) = \tfrac13(\sigma_{11}+\sigma_{22}+\sigma_{33})
```

is the same in every direction — pure confining pressure (negative ``p``) or pure tension
(positive ``p``) — with **zero shear on any plane at all** (try the "Hydrostatic" preset in
the widget above and watch the Mohr circles collapse to a single point on the axis).
Everything else, ``\sigma' = \sigma - p\,\mathbb{I}``, is the **deviatoric stress**: the part
that actually has a shape/shear character and is capable of driving failure. This is the same
split that separates P-wave-like volume change from S-wave-like shape change in the strain
notebook — here it separates the part of ``\sigma`` that just squeezes a rock (irrelevant to
whether it fails in shear) from the part that can actually break it.
"""

# ╔═╡ b1a10001-0000-4000-8000-000000000008
p_example = tr(σ_example) / 3

# ╔═╡ b1a10001-0000-4000-8000-000000000009
dev_example = σ_example - p_example * I

# ╔═╡ b1a10001-0000-4000-8000-00000000000a
md"""
## Mohr Circle for Stress

Plotting ``(\sigma_n,\tau)`` for *every possible* plane orientation traces out a region
bounded by three circles built from the principal stresses — the classical **Mohr diagram**
(right-hand panel of the widget, plotted with **compression positive**, the rock-mechanics
convention). The outer circle spans the two extreme principal stresses,
``C=\tfrac12(\sigma_1+\sigma_3)``, ``R=\tfrac12(\sigma_1-\sigma_3)``; two smaller circles
between each adjacent pair of principal stresses sit inside it. Every physically realizable
plane's ``(\sigma_n,\tau)`` point lies inside the big circle and outside the two small ones —
drag the handle on the sphere and watch the marker sweep across exactly that shaded region,
never outside it.

The point for the *current* dragged plane is computed directly (not just placed
approximately in the shaded region): writing the normal in the principal-axis frame as
``(n_1,n_2,n_3)``,

```math
\sigma_n = \sum_i \sigma_i n_i^2, \qquad
\tau = \sqrt{\sum_i (\sigma_i n_i)^2 - \sigma_n^2}.
```
"""

# ╔═╡ b1a10001-0000-4000-8000-00000000000b
md"""
## Mohr-Coulomb Failure & Pore Pressure

Rock (and fault gouge) resists shear up to a limit that grows with how hard the fault is
being clamped shut — the **Mohr-Coulomb failure criterion**,

```math
\tau = \mu\,\sigma_{n,\mathrm{eff}} + c_0,
```

with friction coefficient ``\mu`` (≈0.6–0.85 for most rock — Byerlee's law) and cohesion
``c_0``. Drawn on the Mohr diagram, it's the red line: any plane whose point sits *above* the
line is unstable — the widget's stability verdict flags exactly this. Because the line only
depends on the *orientation-independent* ``\mu,c_0`` while the circle depends on the full
stress state, whether *any* plane fails at all is really a question about the **big circle**:
if it's small and low enough to clear the line everywhere, the rock is stable on every
orientation; if it grows enough to touch the line, some critically-oriented plane sits right
at failure. The "Critically stressed" preset is built to sit almost exactly on that boundary
— drag the plane normal around and find the orientation where the marker meets the line.

!!! danger "Induced seismicity"
	``\sigma_{n,\mathrm{eff}} = \sigma_n - P_p`` is the **effective stress principle**: pore
	fluid pressure ``P_p`` pushes back against the confining stress, reducing how hard a fault
	is actually clamped shut *without changing the tectonic stress* ``\sigma`` *at all*.
	Pumping fluid underground (wastewater injection, geothermal stimulation, hydraulic
	fracturing) raises ``P_p`` directly. On the widget's Mohr diagram this shows up as a
	**leftward shift** of the point at fixed ``\tau`` — drag the ``P_p`` slider on an
	already-close-to-failure state (try "Critically stressed") and watch the point cross the
	line with *zero* change to the stress tensor itself. This is the accepted mechanism behind
	most documented cases of injection-induced earthquakes.
"""

# ╔═╡ 6f33c7ed-e8a8-440d-a771-84789ee1f397
md"""
#### Appendix
"""

# ╔═╡ b1a10001-0000-4000-8000-00000000000c
md"### The Interactive Widget"

# ╔═╡ b1a10001-0000-4000-8000-00000000000d
begin
	"""
	    StressTensorInput(; s11=1.0, s22=0.0, s33=0.5, s12=0.5, s23=0.0, s13=0.0,
	                        nx=0.0, ny=0.0, nz=1.0, mu=0.6, c0=0.1, pp=0.0)

	A fully self-contained stress-tensor explorer. Drag the six `sᵢⱼ` sliders to set the
	stress tensor, then drag the white handle directly on the sphere to pick which plane
	(unit normal `n`) you're asking about. The sphere shows the traction `t = σ·n`
	decomposed into its normal component (along `n`, red for tension/blue for compression)
	and shear component (tangential, orange — the part that drives a fault to slip), plus the
	three principal-stress axes as a fixed background triad. The right-hand panel is the
	classic 3-circle Mohr diagram for stress (compression plotted positive), with a live
	marker at the current plane's `(σₙ, τ)`, a draggable-by-slider Mohr-Coulomb failure line
	`τ = μσₙ + c₀`, and a pore-pressure slider `Pₚ` that shifts the marker via the
	effective-stress principle `σₙ,eff = σₙ - Pₚ`. Presets jump to canonical stress states.
	"""
	struct StressTensorInput
		s11::Float64
		s22::Float64
		s33::Float64
		s12::Float64
		s23::Float64
		s13::Float64
		nx::Float64
		ny::Float64
		nz::Float64
		mu::Float64
		c0::Float64
		pp::Float64
	end

	StressTensorInput(; s11=1.0, s22=0.0, s33=0.5, s12=0.5, s23=0.0, s13=0.0,
		nx=0.0, ny=0.0, nz=1.0, mu=0.6, c0=0.1, pp=0.0) =
		StressTensorInput(Float64(s11), Float64(s22), Float64(s33), Float64(s12), Float64(s23), Float64(s13),
			Float64(nx), Float64(ny), Float64(nz), Float64(mu), Float64(c0), Float64(pp))

	Base.get(w::StressTensorInput) = Dict{String,Any}(
		"s11" => w.s11, "s22" => w.s22, "s33" => w.s33,
		"s12" => w.s12, "s23" => w.s23, "s13" => w.s13,
		"nx" => w.nx, "ny" => w.ny, "nz" => w.nz,
		"mu" => w.mu, "c0" => w.c0, "pp" => w.pp,
	)

	function Base.show(io::IO, ::MIME"text/html", w::StressTensorInput)
		write(io, """
<div id="stswidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#stswidget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #stswidget .sts-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #stswidget .sts-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #stswidget .sts-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #stswidget .sts-panels{display:flex;gap:14px;align-items:flex-start;justify-content:center;width:100%}
    #stswidget .sts-col{flex:1 1 0;min-width:0;display:flex;flex-direction:column;align-items:center}
    #stswidget .sts-panel-label{font-size:12px;color:#6b7280;margin-bottom:4px;text-transform:uppercase;letter-spacing:.04em}
    #stswidget canvas{cursor:grab;background:#000;border:1px solid #374151;border-radius:6px;display:block;max-width:100%}
    #stswidget canvas.dragging{cursor:grabbing}
    #stswidget .sts-controls{display:flex;gap:10px;flex-wrap:wrap;width:100%;margin-top:12px}
    #stswidget .sts-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px;flex:1 1 240px;min-width:240px}
    #stswidget .sts-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #stswidget select{width:100%;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:5px 6px;font-size:14px;margin-bottom:8px}
    #stswidget .sts-control-row{display:grid;grid-template-columns:minmax(34px,44px) minmax(70px,1fr) minmax(44px,56px);gap:6px;align-items:center;margin:5px 0}
    #stswidget .sts-control-row input[type=range]{width:100%;min-width:0;vertical-align:middle}
    #stswidget .sts-value{color:#d1d5db;text-align:left;font-variant-numeric:tabular-nums;font-size:13px}
    #stswidget .sts-readout{font-size:13px;line-height:1.7}
    #stswidget .sts-mat{font-variant-numeric:tabular-nums;white-space:pre}
    @media (max-width: 900px){
      #stswidget .sts-panels{flex-direction:column;align-items:center}
      #stswidget .sts-col{width:100%;max-width:520px}
    }
  </style>
  <div class="sts-title">
    <div class="sts-title-desc">Traction on a plane splits exactly into the normal stress clamping a fault shut and the shear stress driving it to slip.</div>
    <div class="sts-title-hint">drag the white handle to set the plane normal n &middot; drag empty space to rotate &middot; tune &mu;, c&#8320;, P&#8346; below</div>
  </div>
  <div class="sts-panels">
    <div class="sts-col">
      <div class="sts-panel-label">Traction on the plane (drag n)</div>
      <canvas id="stsSphere" width="480" height="480"></canvas>
    </div>
    <div class="sts-col">
      <div class="sts-panel-label">Mohr circle for stress + Mohr-Coulomb failure</div>
      <canvas id="stsMohr" width="480" height="480"></canvas>
    </div>
  </div>
  <div class="sts-controls">
    <div class="sts-control-group">
      <div class="sts-control-title">Stress Tensor</div>
      <select id="stsPreset">
        <option value="custom">Custom</option>
        <option value="uniC">Uniaxial compression</option>
        <option value="uniT">Uniaxial tension</option>
        <option value="shear">Pure shear</option>
        <option value="hydro">Hydrostatic (confining)</option>
        <option value="crit">Critically stressed</option>
      </select>
      <label class="sts-control-row"><span>s&#8321;&#8321;</span><input type="range" id="sts-s11" min="-1" max="1" step="0.02" value="$(w.s11)"><span id="sts-s11-v" class="sts-value">$(w.s11)</span></label>
      <label class="sts-control-row"><span>s&#8322;&#8322;</span><input type="range" id="sts-s22" min="-1" max="1" step="0.02" value="$(w.s22)"><span id="sts-s22-v" class="sts-value">$(w.s22)</span></label>
      <label class="sts-control-row"><span>s&#8323;&#8323;</span><input type="range" id="sts-s33" min="-1" max="1" step="0.02" value="$(w.s33)"><span id="sts-s33-v" class="sts-value">$(w.s33)</span></label>
      <label class="sts-control-row"><span>s&#8321;&#8322;</span><input type="range" id="sts-s12" min="-1" max="1" step="0.02" value="$(w.s12)"><span id="sts-s12-v" class="sts-value">$(w.s12)</span></label>
      <label class="sts-control-row"><span>s&#8322;&#8323;</span><input type="range" id="sts-s23" min="-1" max="1" step="0.02" value="$(w.s23)"><span id="sts-s23-v" class="sts-value">$(w.s23)</span></label>
      <label class="sts-control-row"><span>s&#8321;&#8323;</span><input type="range" id="sts-s13" min="-1" max="1" step="0.02" value="$(w.s13)"><span id="sts-s13-v" class="sts-value">$(w.s13)</span></label>
    </div>
    <div class="sts-control-group">
      <div class="sts-control-title">Failure Criterion</div>
      <label class="sts-control-row"><span>&mu;</span><input type="range" id="sts-mu" min="0" max="1.2" step="0.01" value="$(w.mu)"><span id="sts-mu-v" class="sts-value">$(w.mu)</span></label>
      <label class="sts-control-row"><span>c&#8320;</span><input type="range" id="sts-c0" min="0" max="0.5" step="0.01" value="$(w.c0)"><span id="sts-c0-v" class="sts-value">$(w.c0)</span></label>
      <label class="sts-control-row"><span>P&#8346;</span><input type="range" id="sts-pp" min="0" max="1.5" step="0.01" value="$(w.pp)"><span id="sts-pp-v" class="sts-value">$(w.pp)</span></label>
      <div style="font-size:12px;color:#6b7280;margin-top:4px">&tau; = &mu;&sigma;<sub>n,eff</sub> + c&#8320;, &sigma;<sub>n,eff</sub> = &sigma;<sub>n</sub> &minus; P&#8346;</div>
    </div>
    <div class="sts-control-group" style="flex:1 1 280px">
      <div class="sts-control-title">Readouts</div>
      <div class="sts-readout" id="stsReadout"></div>
    </div>
    <div class="sts-control-group" style="flex:1 1 280px">
      <div class="sts-control-title">Legend</div>
      <div class="sts-readout" id="stsLegend"></div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8)
  const heightBudget = Math.max(320, window.innerHeight - 560)
  const SEC = Math.round(Math.min((availW-14)/2, heightBudget, 480))
  const DPR = Math.min(window.devicePixelRatio || 1, 2)

  function hidpi(cv, cx, w, h){ cv.width=Math.round(w*DPR); cv.height=Math.round(h*DPR); cv.style.width=w+'px'; cv.style.height=h+'px'; cx.setTransform(DPR,0,0,DPR,0,0) }

  const sphCvs = par.querySelector('#stsSphere'), sctx = sphCvs.getContext('2d')
  const mohrCvs = par.querySelector('#stsMohr'), mctx = mohrCvs.getContext('2d')
  hidpi(sphCvs, sctx, SEC, SEC)
  hidpi(mohrCvs, mctx, SEC, SEC)

  let s11=$(w.s11), s22=$(w.s22), s33=$(w.s33), s12=$(w.s12), s23=$(w.s23), s13=$(w.s13)
  let nx=$(w.nx), ny=$(w.ny), nz=$(w.nz)
  let mu=$(w.mu), c0=$(w.c0), pp=$(w.pp)
  let yaw = 0.7, pitch = 0.35, dragMode = null, lastX = 0, lastY = 0
  let currentPreset = 'custom'

  const PRESETS = {
    uniC:  {s:[-1,    0,     0,    0,0,0]},
    uniT:  {s:[ 1,    0,     0,    0,0,0]},
    shear: {s:[ 0,    0,     0,    1,0,0]},
    hydro: {s:[-0.6, -0.6,  -0.6,  0,0,0]},
    crit:  {s:[-0.86,-0.16, -0.16, 0,0,0]},
  }

  const CX = SEC/2, CY = SEC/2, RPIX = SEC*0.34

  function rot(p){
    const x = p[0]*Math.cos(yaw) - p[1]*Math.sin(yaw)
    const y = p[0]*Math.sin(yaw) + p[1]*Math.cos(yaw)
    const z = p[2]
    return [x, y*Math.cos(pitch) - z*Math.sin(pitch), y*Math.sin(pitch) + z*Math.cos(pitch)]
  }
  function unrot(P){
    const y1 = P[1]*Math.cos(pitch) + P[2]*Math.sin(pitch)
    const z1 = -P[1]*Math.sin(pitch) + P[2]*Math.cos(pitch)
    const px = P[0]*Math.cos(yaw) + y1*Math.sin(yaw)
    const py = -P[0]*Math.sin(yaw) + y1*Math.cos(yaw)
    return [px, py, z1]
  }
  function toScreen(p){ const r = rot(p); return [CX + r[0]*RPIX, CY - r[2]*RPIX] }
  // Inverse of toScreen(), restricted to the front hemisphere -- given a screen click, solve
  // for the point on the unit sphere facing the camera, then unrotate it back to world (model)
  // coordinates. Clamped to the silhouette instead of returning null so dragging near the
  // sphere's edge still tracks smoothly.
  function screenToXYZ(mx, my){
    let px = (mx-CX)/RPIX, pz = (CY-my)/RPIX
    let r2 = px*px + pz*pz
    if(r2 > 1){ const s = 1/Math.sqrt(r2); px *= s; pz *= s; r2 = 1 }
    const py = Math.sqrt(Math.max(0, 1-r2))
    const v = unrot([px, py, pz])
    const n = Math.hypot(v[0], v[1], v[2]) || 1
    return [v[0]/n, v[1]/n, v[2]/n]
  }

  // Jacobi eigenvalue algorithm for a symmetric 3x3 -- classic cyclic-pivot-free ("largest
  // off-diagonal") variant, a handful of iterations converge to machine precision for 3x3.
  function jacobiEigen(Ain){
    const a = Ain.map(r => r.slice())
    const v = [[1,0,0],[0,1,0],[0,0,1]]
    for(let iter=0; iter<60; iter++){
      let p=0, q=1, mx=Math.abs(a[0][1])
      if(Math.abs(a[0][2])>mx){ mx=Math.abs(a[0][2]); p=0; q=2 }
      if(Math.abs(a[1][2])>mx){ mx=Math.abs(a[1][2]); p=1; q=2 }
      if(mx < 1e-13) break
      const apq=a[p][q], app=a[p][p], aqq=a[q][q]
      const theta = (aqq-app)/(2*apq)
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta*theta+1))
      const c = 1/Math.sqrt(t*t+1), s = t*c
      const newapp = app - t*apq, newaqq = aqq + t*apq
      a[p][p]=newapp; a[q][q]=newaqq; a[p][q]=0; a[q][p]=0
      for(let k=0;k<3;k++){
        if(k!==p && k!==q){
          const akp=a[k][p], akq=a[k][q]
          a[k][p]=a[p][k]=c*akp - s*akq
          a[k][q]=a[q][k]=s*akp + c*akq
        }
      }
      for(let k=0;k<3;k++){
        const vkp=v[k][p], vkq=v[k][q]
        v[k][p]=c*vkp - s*vkq
        v[k][q]=s*vkp + c*vkq
      }
    }
    const vals=[a[0][0],a[1][1],a[2][2]]
    const order=[0,1,2].sort((i,j)=>vals[j]-vals[i])
    return { values: order.map(i=>vals[i]), vectors: order.map(i=>[v[0][i],v[1][i],v[2][i]]) }
  }

  function sigmaMat(){ return [[s11,s12,s13],[s12,s22,s23],[s13,s23,s33]] }
  function dot3(a,b){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
  function sub3(a,b){ return [a[0]-b[0],a[1]-b[1],a[2]-b[2]] }
  function scale3(a,k){ return [a[0]*k,a[1]*k,a[2]*k] }
  function norm3(a){ return Math.hypot(a[0],a[1],a[2]) }
  function traction(nvec){
    const S = sigmaMat()
    return [S[0][0]*nvec[0]+S[0][1]*nvec[1]+S[0][2]*nvec[2],
            S[1][0]*nvec[0]+S[1][1]*nvec[1]+S[1][2]*nvec[2],
            S[2][0]*nvec[0]+S[2][1]*nvec[1]+S[2][2]*nvec[2]]
  }
  function tractionDecomp(){
    const nv = [nx,ny,nz]
    const t = traction(nv)
    const sn = dot3(nv, t)
    const nvecPart = scale3(nv, sn)
    const shear = sub3(t, nvecPart)
    return { n: nv, t, sn, nvecPart, shear, tau: norm3(shear) }
  }
  function principal(){ return jacobiEigen(sigmaMat()) }

  function mohrCircles(pr){
    // pr.values: [σ1,σ2,σ3] descending, tension-positive -- flip to compression-positive
    // for the plotted diagram (rock-mechanics convention).
    const Scomp = [-pr.values[2], -pr.values[1], -pr.values[0]]
    const circ = (a,b) => ({ c: (a+b)/2, r: Math.abs(a-b)/2 })
    return { big: circ(Scomp[0], Scomp[2]), c12: circ(Scomp[0], Scomp[1]), c23: circ(Scomp[1], Scomp[2]) }
  }
  function mohrPointFor(decomp){
    const xPlot = -decomp.sn
    return { x: xPlot, xEff: xPlot - pp, y: decomp.tau }
  }
  function verdict(xEff, y){
    const lineY = mu*xEff + c0
    const margin = lineY - y
    if(margin < -1e-6) return { text: 'UNSTABLE — this plane would slip', color: '#ef4444' }
    if(margin < 0.02) return { text: 'AT FAILURE — right on the edge', color: '#facc15' }
    return { text: 'stable', color: '#4ade80' }
  }

  function drawLine3(ctx, p0, p1, color, width){
    const s0=toScreen(p0), s1=toScreen(p1)
    ctx.strokeStyle=color; ctx.lineWidth=width||1.3
    ctx.beginPath(); ctx.moveTo(s0[0],s0[1]); ctx.lineTo(s1[0],s1[1]); ctx.stroke()
  }
  function drawArrow3(ctx, p0, p1, color, width){
    const s0=toScreen(p0), s1=toScreen(p1)
    ctx.strokeStyle=color; ctx.lineWidth=width||2.2
    ctx.beginPath(); ctx.moveTo(s0[0],s0[1]); ctx.lineTo(s1[0],s1[1]); ctx.stroke()
    const ang = Math.atan2(s1[1]-s0[1], s1[0]-s0[0])
    ctx.beginPath(); ctx.fillStyle=color
    ctx.moveTo(s1[0],s1[1])
    ctx.lineTo(s1[0]-7*Math.cos(ang-0.4), s1[1]-7*Math.sin(ang-0.4))
    ctx.lineTo(s1[0]-7*Math.cos(ang+0.4), s1[1]-7*Math.sin(ang+0.4))
    ctx.closePath(); ctx.fill()
  }

  function drawSphereWire(){
    sctx.clearRect(0,0,SEC,SEC)
    sctx.strokeStyle = 'rgba(107,114,128,0.35)'; sctx.lineWidth = 1
    const NLAT=6, NLON=10, NSEG=48
    for(let i=1;i<NLAT;i++){
      const theta = Math.PI*i/NLAT
      sctx.beginPath(); let started=false
      for(let j=0;j<=NSEG;j++){
        const phi = 2*Math.PI*j/NSEG
        const p = [Math.sin(theta)*Math.cos(phi), Math.sin(theta)*Math.sin(phi), Math.cos(theta)]
        if(rot(p)[1] < 0){ started=false; continue }
        const s = toScreen(p)
        if(!started){ sctx.moveTo(s[0],s[1]); started=true } else sctx.lineTo(s[0],s[1])
      }
      sctx.stroke()
    }
    for(let j=0;j<NLON;j++){
      const phi = Math.PI*j/NLON
      sctx.beginPath(); let started=false
      for(let i=0;i<=NSEG;i++){
        const theta = Math.PI*i/NSEG
        const p = [Math.sin(theta)*Math.cos(phi), Math.sin(theta)*Math.sin(phi), Math.cos(theta)]
        if(rot(p)[1] < 0){ started=false; continue }
        const s = toScreen(p)
        if(!started){ sctx.moveTo(s[0],s[1]); started=true } else sctx.lineTo(s[0],s[1])
      }
      sctx.stroke()
    }
    sctx.beginPath(); sctx.arc(CX,CY,RPIX,0,2*Math.PI); sctx.strokeStyle='#4b5563'; sctx.lineWidth=1.2; sctx.stroke()
  }

  function drawSphere(){
    drawSphereWire()
    const pr = principal()
    for(let k=0;k<3;k++){
      const v = pr.vectors[k], val = pr.values[k]
      const color = val >= 0 ? 'rgba(239,68,68,0.5)' : 'rgba(59,130,246,0.5)'
      drawLine3(sctx, [-v[0]*1.25,-v[1]*1.25,-v[2]*1.25], [v[0]*1.25,v[1]*1.25,v[2]*1.25], color, 1.3)
    }
    const d = tractionDecomp()
    const P = d.n
    const SCALE = 0.7
    drawArrow3(sctx, P, [P[0]+d.t[0]*SCALE, P[1]+d.t[1]*SCALE, P[2]+d.t[2]*SCALE], 'rgba(229,231,235,0.55)', 1.6)
    const ncolor = d.sn >= 0 ? '#ef4444' : '#3b82f6'
    drawArrow3(sctx, P, [P[0]+d.nvecPart[0]*SCALE, P[1]+d.nvecPart[1]*SCALE, P[2]+d.nvecPart[2]*SCALE], ncolor, 2.6)
    drawArrow3(sctx, P, [P[0]+d.shear[0]*SCALE, P[1]+d.shear[1]*SCALE, P[2]+d.shear[2]*SCALE], '#fb923c', 2.6)
    drawLine3(sctx, [0,0,0], P, 'rgba(245,243,239,0.6)', 1.3)
    const hs = toScreen(P)
    sctx.beginPath(); sctx.arc(hs[0],hs[1],7,0,2*Math.PI); sctx.fillStyle='#f5f3ef'; sctx.fill()
    sctx.strokeStyle='#0a0f18'; sctx.lineWidth=1.5; sctx.stroke()
    sctx.fillStyle='#9ca3af'; sctx.font='13px sans-serif'
    sctx.fillText('drag the white handle to set the plane normal n', 10, 20)
  }

  function drawMohr(){
    const pr = principal()
    const mc = mohrCircles(pr)
    const d = tractionDecomp()
    const pt = mohrPointFor(d)

    const maxR = Math.max(mc.big.r, 0.15)
    let lo = Math.min(mc.big.c - mc.big.r, pt.xEff, 0) - maxR*0.3
    let hi = Math.max(mc.big.c + mc.big.r, pt.x, 0) + maxR*0.3
    if(hi - lo < 0.6){ const mid=(lo+hi)/2; lo=mid-0.3; hi=mid+0.3 }
    const yhi = Math.max(maxR*1.3, pt.y*1.3, 0.3)

    // A single isotropic scale for both axes -- using independent x/y scales would draw
    // the circles as ellipses, which would misrepresent tangency to the failure line (the
    // whole point of the diagram) even though the underlying (σₙ,τ) data is correct.
    const MARGIN = SEC*0.08
    const xCenter = (lo+hi)/2
    const scale = Math.min((SEC-2*MARGIN)/(hi-lo), (SEC-2*MARGIN)/yhi)
    const to = (x,y) => [ SEC/2 + (x-xCenter)*scale, (SEC-MARGIN) - y*scale ]

    mctx.clearRect(0,0,SEC,SEC)

    // Mohr-Coulomb failure line + unstable shading (above the line)
    const yL0=mu*lo+c0, yL1=mu*hi+c0
    const pL0=to(lo,yL0), pL1=to(hi,yL1)
    const top0=to(hi,yhi), top1=to(lo,yhi)
    mctx.fillStyle='rgba(239,68,68,0.10)'
    mctx.beginPath(); mctx.moveTo(pL0[0],pL0[1]); mctx.lineTo(pL1[0],pL1[1])
    mctx.lineTo(top0[0],top0[1]); mctx.lineTo(top1[0],top1[1]); mctx.closePath(); mctx.fill()
    mctx.strokeStyle='#ef4444'; mctx.lineWidth=2
    mctx.beginPath(); mctx.moveTo(pL0[0],pL0[1]); mctx.lineTo(pL1[0],pL1[1]); mctx.stroke()

    // axes
    mctx.strokeStyle='rgba(156,163,175,0.4)'; mctx.lineWidth=1
    let a0=to(lo,0), a1=to(hi,0); mctx.beginPath(); mctx.moveTo(a0[0],a0[1]); mctx.lineTo(a1[0],a1[1]); mctx.stroke()
    a0=to(0,0); a1=to(0,yhi); mctx.beginPath(); mctx.moveTo(a0[0],a0[1]); mctx.lineTo(a1[0],a1[1]); mctx.stroke()

    function drawCircle(circ, color, lw){
      mctx.beginPath()
      for(let i=0;i<=100;i++){
        const th = Math.PI*i/100
        const p = to(circ.c+circ.r*Math.cos(th), circ.r*Math.sin(th))
        i===0 ? mctx.moveTo(p[0],p[1]) : mctx.lineTo(p[0],p[1])
      }
      mctx.strokeStyle=color; mctx.lineWidth=lw; mctx.stroke()
    }
    drawCircle(mc.c12, 'rgba(56,189,248,0.55)', 1.4)
    drawCircle(mc.c23, 'rgba(56,189,248,0.55)', 1.4)
    drawCircle(mc.big, '#38bdf8', 2)

    if(Math.abs(pp) > 1e-6){
      const p0=to(pt.x, pt.y), p1=to(pt.xEff, pt.y)
      mctx.strokeStyle='rgba(196,132,252,0.85)'; mctx.lineWidth=1.6
      mctx.beginPath(); mctx.moveTo(p0[0],p0[1]); mctx.lineTo(p1[0],p1[1]); mctx.stroke()
      mctx.beginPath(); mctx.arc(p0[0],p0[1],3,0,2*Math.PI); mctx.fillStyle='rgba(196,132,252,0.85)'; mctx.fill()
    }

    const vd = verdict(pt.xEff, pt.y)
    const pm = to(pt.xEff, pt.y)
    mctx.beginPath(); mctx.arc(pm[0],pm[1],7,0,2*Math.PI)
    mctx.fillStyle=vd.color; mctx.fill()
    mctx.strokeStyle='#0a0f18'; mctx.lineWidth=1.5; mctx.stroke()

    mctx.fillStyle='#9ca3af'; mctx.font='12px sans-serif'
    mctx.fillText('compressive normal stress σₙ →', 10, SEC-10)
    mctx.save(); mctx.translate(14,20); mctx.fillText('↑ shear stress τ', 0, 0); mctx.restore()
  }

  function fmt(v){ return (v>=0?'+':'') + v.toFixed(2) }
  function matStr3(S){
    return '['+fmt(S[0][0])+'  '+fmt(S[0][1])+'  '+fmt(S[0][2])+']\\n'+
           '['+fmt(S[1][0])+'  '+fmt(S[1][1])+'  '+fmt(S[1][2])+']\\n'+
           '['+fmt(S[2][0])+'  '+fmt(S[2][1])+'  '+fmt(S[2][2])+']'
  }

  function updateReadouts(){
    const pr = principal()
    const d = tractionDecomp()
    const pt = mohrPointFor(d)
    const S = sigmaMat()
    const p = (s11+s22+s33)/3
    let devNorm = 0
    for(let i=0;i<3;i++) for(let j=0;j<3;j++){ const dv = S[i][j] - (i===j? p:0); devNorm += dv*dv }
    devNorm = Math.sqrt(devNorm)
    const vd = verdict(pt.xEff, pt.y)

    par.querySelector('#stsReadout').innerHTML =
      '<b>&sigma;</b><br><span class="sts-mat">'+matStr3(S)+'</span><br>'+
      'principal: <b>'+fmt(pr.values[0])+'</b>, <b>'+fmt(pr.values[1])+'</b>, <b>'+fmt(pr.values[2])+'</b><br>'+
      'mean stress p = <b>'+fmt(p)+'</b> &middot; &Vert;dev&Vert; = <b>'+devNorm.toFixed(2)+'</b><br>'+
      'plane: n = ('+d.n.map(v=>v.toFixed(2)).join(', ')+')<br>'+
      '&sigma;<sub>n</sub> = <b>'+fmt(d.sn)+'</b> &middot; &tau; = <b>'+d.tau.toFixed(2)+'</b><br>'+
      '<b style="color:'+vd.color+'">'+vd.text+'</b>'

    par.querySelector('#stsLegend').innerHTML =
      '<span style="color:#ef4444">red</span>/<span style="color:#3b82f6">blue</span> axes = tensile/compressive principal stress<br>'+
      '<span style="color:#ef4444">normal</span> traction (red=tension, blue=compression) &middot; '+
      '<span style="color:#fb923c">shear</span> traction &middot; gray = full traction t=&sigma;&middot;n<br>'+
      'Mohr diagram: compression plotted positive &middot; <span style="color:#ef4444">red line</span> = Mohr-Coulomb failure &middot; '+
      '<span style="color:#c084fc">purple</span> shift = pore pressure'
  }

  function drawAll(){ drawSphere(); drawMohr(); updateReadouts() }
  function setCustom(){ if(currentPreset !== 'custom'){ currentPreset='custom'; par.querySelector('#stsPreset').value='custom' } }

  drawAll()

  function bindSlider(id, setter, decimals){
    const el = par.querySelector('#'+id)
    el.addEventListener('input', e=>{
      const v = parseFloat(e.target.value)
      setter(v)
      par.querySelector('#'+id+'-v').textContent = v.toFixed(decimals)
      setCustom()
      drawAll()
    })
  }
  bindSlider('sts-s11', v=>s11=v, 2)
  bindSlider('sts-s22', v=>s22=v, 2)
  bindSlider('sts-s33', v=>s33=v, 2)
  bindSlider('sts-s12', v=>s12=v, 2)
  bindSlider('sts-s23', v=>s23=v, 2)
  bindSlider('sts-s13', v=>s13=v, 2)

  function bindPlainSlider(id, setter, decimals){
    const el = par.querySelector('#'+id)
    el.addEventListener('input', e=>{
      const v = parseFloat(e.target.value)
      setter(v)
      par.querySelector('#'+id+'-v').textContent = v.toFixed(decimals)
      drawAll()
    })
  }
  bindPlainSlider('sts-mu', v=>mu=v, 2)
  bindPlainSlider('sts-c0', v=>c0=v, 2)
  bindPlainSlider('sts-pp', v=>pp=v, 2)

  par.querySelector('#stsPreset').addEventListener('change', e=>{
    const key = e.target.value
    currentPreset = key
    if(key !== 'custom'){
      const [a,b,c,d2,f,g] = PRESETS[key].s
      s11=a; s22=b; s33=c; s12=d2; s23=f; s13=g
      const vals = {s11:a, s22:b, s33:c, s12:d2, s23:f, s13:g}
      for(const k in vals){
        par.querySelector('#sts-'+k).value = vals[k]
        par.querySelector('#sts-'+k+'-v').textContent = vals[k].toFixed(2)
      }
    }
    drawAll()
  })

  sphCvs.addEventListener('mousedown', e=>{
    const mx=e.offsetX, my=e.offsetY
    const hs = toScreen([nx,ny,nz])
    if(Math.hypot(hs[0]-mx, hs[1]-my) < 14){ dragMode='n' }
    else { dragMode='rotate'; lastX=mx; lastY=my }
  })
  sphCvs.addEventListener('mousemove', e=>{
    if(!dragMode) return
    if(dragMode === 'n'){
      const v = screenToXYZ(e.offsetX, e.offsetY)
      nx=v[0]; ny=v[1]; nz=v[2]
      drawAll()
    } else if(dragMode === 'rotate'){
      const dx=e.offsetX-lastX, dy=e.offsetY-lastY
      lastX=e.offsetX; lastY=e.offsetY
      yaw += dx*0.01
      pitch = Math.max(-1.3, Math.min(1.3, pitch+dy*0.01))
      drawAll()
    }
  })
  window.addEventListener('mouseup', ()=>{ dragMode=null })
</script>
""")
	end

	const _sts_ready = true
end

# ╔═╡ b1a10001-0000-4000-8000-00000000000e
begin
	_sts_ready
	@bind sts StressTensorInput()
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Einsum = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
Einsum = "~0.4.1"
PlutoUI = "~0.7.83"
SymbolicUtils = "~4.44.0"
Symbolics = "~7.35.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "3956ebaae45e166ac285d2101d18e357833d57bf"

[[deps.ADTypes]]
git-tree-sha1 = "9b38b82a9fe131f3d331a53b7203d9d1a2a4602c"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.22.4"

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

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

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

[[deps.Einsum]]
deps = ["Compat"]
git-tree-sha1 = "4a6b3eee0161c89700b6c1949feae8b851da5494"
uuid = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
version = "0.4.1"

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
git-tree-sha1 = "05f45c2e0de6259db764adbfd2f1dc6d3f8de13c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.1"

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
git-tree-sha1 = "cd4b115137894ced9830a92bcdb95a6bd8f38880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.2"

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
git-tree-sha1 = "6c7e35e1984e98f6e314dcd4c92cd68842af5456"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.52"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "ab885203e8395593d65b629bd4023de089e6997b"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "1.1.4"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "Graphs", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "03bbe242c7433bfca3660050d0b0cc3b4be8df71"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.44.0"

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
git-tree-sha1 = "2a8387ae58e5c1afbd8820b5c3ddb300d15ae27f"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.35.0"

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
# ╠═3f1bc2d1-2327-48cc-984a-df09c936da87
# ╟─fd14256d-8c69-4d8a-91b5-924a32479866
# ╟─b1a10001-0000-4000-8000-00000000000e
# ╟─b1a10001-0000-4000-8000-000000000001
# ╠═50eb4916-b3d7-4f2a-90c6-8661cbbd8e7a
# ╠═88b33835-d993-4396-8605-bb3456200eb1
# ╠═b1a10001-0000-4000-8000-000000000002
# ╠═b1a10001-0000-4000-8000-000000000003
# ╠═b1a10001-0000-4000-8000-000000000004
# ╟─b1a10001-0000-4000-8000-000000000005
# ╠═ca99052e-ff24-4e8a-9d06-74bc24974873
# ╠═987d9253-2f7d-4ca8-8eef-aae8e8d53bc6
# ╟─b1a10001-0000-4000-8000-000000000006
# ╟─b1a10001-0000-4000-8000-000000000007
# ╠═b1a10001-0000-4000-8000-000000000008
# ╠═b1a10001-0000-4000-8000-000000000009
# ╟─b1a10001-0000-4000-8000-00000000000a
# ╟─b1a10001-0000-4000-8000-00000000000b
# ╟─6f33c7ed-e8a8-440d-a771-84789ee1f397
# ╠═13a3429e-12f6-11ed-326f-c154f5debceb
# ╟─b1a10001-0000-4000-8000-00000000000c
# ╠═b1a10001-0000-4000-8000-00000000000d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
