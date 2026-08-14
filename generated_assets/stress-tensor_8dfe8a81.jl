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
    #stswidget .sts-mat3{display:inline-grid;grid-template-columns:repeat(3,auto);gap:5px 6px;padding:6px 12px;position:relative;margin-bottom:8px}
    #stswidget .sts-mat3::before, #stswidget .sts-mat3::after{content:'';position:absolute;top:3px;bottom:3px;width:7px;border:2px solid #9ca3af}
    #stswidget .sts-mat3::before{left:0;border-right:none}
    #stswidget .sts-mat3::after{right:0;border-left:none}
    #stswidget .sts-mat3 input{width:46px;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:4px 3px;font-size:13px;text-align:center;font-variant-numeric:tabular-nums}
    #stswidget .sts-mat3 input:focus{outline:2px solid #38bdf8;border-color:#38bdf8}
    #stswidget .sts-matro{display:inline-grid;grid-template-columns:repeat(3,auto);gap:2px 6px;padding:3px 9px;position:relative;vertical-align:middle;margin:2px 0}
    #stswidget .sts-matro::before, #stswidget .sts-matro::after{content:'';position:absolute;top:1px;bottom:1px;width:5px;border:1.5px solid currentColor}
    #stswidget .sts-matro::before{left:0;border-right:none}
    #stswidget .sts-matro::after{right:0;border-left:none}
    #stswidget .sts-matro span{min-width:34px;text-align:center;font-size:12px;font-variant-numeric:tabular-nums}
    @media (max-width: 900px){
      #stswidget .sts-panels{flex-direction:column;align-items:center}
      #stswidget .sts-col{width:100%;max-width:520px}
    }
  </style>
  <div class="sts-title">
    <div class="sts-title-desc">Traction on a plane splits exactly into the normal stress clamping a fault shut and the shear stress driving it to slip.</div>
    <div class="sts-title-hint">drag the white handle to set the plane normal n &middot; edit &sigma; directly or drag its sliders &middot; drag empty space to rotate &middot; tune &mu;, c&#8320;, P&#8346; below</div>
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
      <div class="sts-mat3">
        <input type="number" step="0.02" id="sts-mat-11" value="$(w.s11)">
        <input type="number" step="0.02" id="sts-mat-12" value="$(w.s12)">
        <input type="number" step="0.02" id="sts-mat-13" value="$(w.s13)">
        <input type="number" step="0.02" id="sts-mat-21" value="$(w.s12)">
        <input type="number" step="0.02" id="sts-mat-22" value="$(w.s22)">
        <input type="number" step="0.02" id="sts-mat-23" value="$(w.s23)">
        <input type="number" step="0.02" id="sts-mat-31" value="$(w.s13)">
        <input type="number" step="0.02" id="sts-mat-32" value="$(w.s23)">
        <input type="number" step="0.02" id="sts-mat-33" value="$(w.s33)">
      </div>
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
  function traction(nvec, S=sigmaMat()){
    return [S[0][0]*nvec[0]+S[0][1]*nvec[1]+S[0][2]*nvec[2],
            S[1][0]*nvec[0]+S[1][1]*nvec[1]+S[1][2]*nvec[2],
            S[2][0]*nvec[0]+S[2][1]*nvec[1]+S[2][2]*nvec[2]]
  }
  function tractionDecomp(S=sigmaMat()){
    const nv = [nx,ny,nz]
    const t = traction(nv, S)
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

  // fixed x1/x2/x3 reference frame the stress components (s11, s12, ...) are actually defined
  // against -- without this the wireframe sphere has no visible orientation cue at all, so
  // there's no way to tell which direction on screen is "1" vs "2" vs "3".
  function drawAxes(){
    const AXLEN = 1.18
    const axes = [ [[1,0,0],'x₁'], [[0,1,0],'x₂'], [[0,0,1],'x₃'] ]
    sctx.lineWidth = 1
    sctx.font = '13px sans-serif'; sctx.fillStyle = '#9ca3af'
    sctx.textAlign = 'center'; sctx.textBaseline = 'middle'
    for(const [v,label] of axes){
      drawLine3(sctx, [-v[0]*AXLEN,-v[1]*AXLEN,-v[2]*AXLEN], [v[0]*AXLEN,v[1]*AXLEN,v[2]*AXLEN], 'rgba(156,163,175,0.55)', 1)
      const s1 = toScreen([v[0]*AXLEN, v[1]*AXLEN, v[2]*AXLEN])
      sctx.fillText(label, s1[0], s1[1])
    }
    sctx.textAlign = 'left'; sctx.textBaseline = 'alphabetic'
  }

  function drawSphere(state){
    drawSphereWire()
    drawAxes()
    const {pr, d} = state
    for(let k=0;k<3;k++){
      const v = pr.vectors[k], val = pr.values[k]
      const color = val >= 0 ? 'rgba(239,68,68,0.5)' : 'rgba(59,130,246,0.5)'
      drawLine3(sctx, [-v[0]*1.25,-v[1]*1.25,-v[2]*1.25], [v[0]*1.25,v[1]*1.25,v[2]*1.25], color, 1.3)
    }
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

  function drawMohr(state){
    const {pr, mc, d, pt, vd} = state

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

    const pm = to(pt.xEff, pt.y)
    mctx.beginPath(); mctx.arc(pm[0],pm[1],7,0,2*Math.PI)
    mctx.fillStyle=vd.color; mctx.fill()
    mctx.strokeStyle='#0a0f18'; mctx.lineWidth=1.5; mctx.stroke()

    mctx.fillStyle='#9ca3af'; mctx.font='12px sans-serif'
    mctx.fillText('compressive normal stress σₙ →', 10, SEC-10)
    mctx.save(); mctx.translate(14,20); mctx.fillText('↑ shear stress τ', 0, 0); mctx.restore()
  }

  function fmt(v){ return (v>=0?'+':'') + v.toFixed(2) }
  // Real HTML bracket-grid matrix (same visual trick as the editable .sts-mat3 panel above,
  // read-only here) instead of a monospace ASCII-bracket string.
  function matHTML3(S, color){
    let html = '<span class="sts-matro" style="color:'+color+'">'
    for(let i=0;i<3;i++) for(let j=0;j<3;j++) html += '<span style="color:#e5e7eb">'+fmt(S[i][j])+'</span>'
    return html + '</span>'
  }

  function updateReadouts(state){
    const {pr, d, S, vd} = state
    const p = (s11+s22+s33)/3
    let devNorm = 0
    for(let i=0;i<3;i++) for(let j=0;j<3;j++){ const dv = S[i][j] - (i===j? p:0); devNorm += dv*dv }
    devNorm = Math.sqrt(devNorm)
    par.querySelector('#stsReadout').innerHTML =
      '<b style="color:#38bdf8">&sigma;</b> '+matHTML3(S, '#38bdf8')+'<br>'+
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

  // The nine σ matrix boxes mirror s11..s13 both ways: dragging a slider (or picking a preset)
  // updates the boxes here, and typing in a box updates the underlying variable + its slider
  // (registered below, next to bindSlider) -- symmetric off-diagonal pairs (12/21, 13/31,
  // 23/32) share one variable, so editing either box updates both. Skip overwriting whichever
  // box currently has focus so a live edit isn't clobbered mid-keystroke.
  const matBoxIds = {s11:['sts-mat-11'], s22:['sts-mat-22'], s33:['sts-mat-33'],
                      s12:['sts-mat-12','sts-mat-21'], s13:['sts-mat-13','sts-mat-31'], s23:['sts-mat-23','sts-mat-32']}
  function syncMatBoxes(){
    const vals = {s11,s22,s33,s12,s13,s23}
    for(const k in matBoxIds) for(const id of matBoxIds[k]){
      const el = par.querySelector('#'+id)
      if(document.activeElement !== el) el.value = Math.round(vals[k]*1000)/1000
    }
  }

  function drawAll(){
    const S = sigmaMat()
    const pr = jacobiEigen(S)
    const d = tractionDecomp(S)
    const pt = mohrPointFor(d)
    const state = {S, pr, d, pt, mc:mohrCircles(pr), vd:verdict(pt.xEff, pt.y)}
    drawSphere(state); drawMohr(state); updateReadouts(state); syncMatBoxes()
  }
  function setCustom(){ if(currentPreset !== 'custom'){ currentPreset='custom'; par.querySelector('#stsPreset').value='custom' } }

  function wireMatBoxes(ids, sliderId, setter){
    const els = ids.map(id => par.querySelector('#'+id))
    function handle(v){
      setter(v)
      const sl = par.querySelector('#'+sliderId)
      sl.value = v
      par.querySelector('#'+sliderId+'-v').textContent = v.toFixed(2)
      setCustom()
      drawAll()
    }
    els.forEach(el => el.addEventListener('input', e=>{
      const v = parseFloat(e.target.value)
      if(Number.isFinite(v)) handle(v)
    }))
  }
  wireMatBoxes(['sts-mat-11'], 'sts-s11', v=>s11=v)
  wireMatBoxes(['sts-mat-22'], 'sts-s22', v=>s22=v)
  wireMatBoxes(['sts-mat-33'], 'sts-s33', v=>s33=v)
  wireMatBoxes(['sts-mat-12','sts-mat-21'], 'sts-s12', v=>s12=v)
  wireMatBoxes(['sts-mat-13','sts-mat-31'], 'sts-s13', v=>s13=v)
  wireMatBoxes(['sts-mat-23','sts-mat-32'], 'sts-s23', v=>s23=v)

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
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "9bbb622e4cd9995606b539e0b3f2495d359cd8e4"

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
# ╠═3f1bc2d1-2327-48cc-984a-df09c936da87
# ╟─fd14256d-8c69-4d8a-91b5-924a32479866
# ╟─b1a10001-0000-4000-8000-00000000000e
# ╟─b1a10001-0000-4000-8000-000000000001
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
