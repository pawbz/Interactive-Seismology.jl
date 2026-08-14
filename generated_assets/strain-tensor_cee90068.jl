### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> order = "3"
#> title = "Strain Tensor"
#> tags = ["introduction"]
#> layout = "layout.jlhtml"
#> description = "Deformation of a 2-D square element"

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
using PlutoUI

# ╔═╡ 3f1bc2d1-2327-48cc-984a-df09c936da87
TableOfContents()

# ╔═╡ fd14256d-8c69-4d8a-91b5-924a32479866
md"""
# Strain Tensor

If you split any small piece of a deforming body into a translation, a rotation, and a
*change of shape*, that last piece — how much a small square actually stretches, squashes,
or shears, independent of where it moves to or which way it turns — is the **strain**. It is
the part of a deformation an elastic medium resists (and therefore radiates seismic energy
to relax), and the part that matters for whether rock breaks. Real strain rates are tiny.

!!! note "How tiny"
	Geodetic (GPS) measurements across the Southern California plate boundary show relative
	plate motion of about 2-7 cm/year, a strain *rate* of roughly ``3\times10^{-7}`` per
	year — about ten times faster than the strain accumulating deep in the plate interior
	(``<3\times10^{-8}`` per year). A large earthquake releases a *coseismic* shear strain
	drop of order ``3\times10^{-4}`` in one rupture — five orders of magnitude more than a
	year's worth of loading. The **Earthquake Recurrence** section near the end of this
	notebook turns those two numbers into a live estimate of how often a fault like that
	should rupture.

This notebook builds the strain tensor up from a 2-D Jacobian you can drag by hand: how a
deformation decomposes into rotation + strain (Toeplitz decomposition), how strain itself
splits further into dilatation (uniform area/volume change) and deviatoric strain (pure
shape change at constant area/volume), what principal strains and the strain ellipse are,
and the Mohr circle construction for reading off the strain on *any* cut through the
material.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 4d860dc6-4944-44da-81f2-14b4aa018694
md"""
## Configurations, Bodies, and the Lagrangian Description

Everything in this notebook describes a **body** — a piece of material made of some fixed set of
particles — as it moves and deforms. We single out one instant as the **reference
configuration**, where each particle sits at its reference position ``\mathbf{x}`` (the same
``\mathbf{x}`` used throughout this notebook); at any later time ``t`` the body occupies its
**deformed** (or *current*) **configuration**, where that same particle has moved to
``\mathbf{x}+\mathbf{u}(\mathbf{x},t)``.

Tracking a property by *which particle* carries it — rather than by *which spatial point*
currently happens to hold it — is the **Lagrangian description**, and it's what every
``\mathbf{u}(\mathbf{x},t)`` in this notebook means: the displacement of the *particle labeled*
``\mathbf{x}``, not of "whatever is currently at some fixed spatial location." (The alternative —
tracking properties at fixed spatial points as different particles pass through them — is the
*Eulerian* description: more natural for fluids in motion than for a solid Earth we mostly care
about at rest between earthquakes.)

The widget below drags **a part of the body** — a small region of four labeled particles — from
its reference configuration to any deformed configuration you like, by hand. Two things to watch
as you drag:

- The numeric labels stay attached to their own particle — dragging relabels nothing, exactly
  the Lagrangian point above.
- **Mass** is an **extensive** property: however you deform the region, the *same* particles are
  still in it, so its mass ``M=\rho_0 A_0`` never changes. **Density** ``\rho=M/A`` is not — it's
  mass per unit *current* area, so it changes the instant the region's area does, in exactly the
  way that keeps ``\rho A=\rho_0 A_0`` (the region's mass) constant. This is the 2-D statement of
  mass conservation, ``\rho_0\,dA_0=\rho\,dA``, behind the continuity equation used throughout
  continuum seismology.
"""

# ╔═╡ 174e1f52-bf7a-4201-b579-c784115d15f1
md"""
## Displacement
Particle displacement is a function of space and time
```math
𝐮 = 𝐮(𝐱, t).
```
Here, 𝐱 denotes the position the particle occupies at a reference time t=0, which means, a particle at position 𝐱 will be moved to position 𝐱+𝐮(𝐱, t) at time t — Lagrangian, since it's the *same particle* 𝐱 being tracked. The particle velocity and particle acceleration are given by
∂ₜ𝐮 and ∂ₜ²𝐮, respectively.
"""

# ╔═╡ b2cf2af1-564e-4a3f-a846-a50042883a7a
md"""
## Distortion
Distortion of the medium changes the relative position of nearby particles — it's therefore
related to the *gradients* of ``\mathbf{u}``, not ``\mathbf{u}`` itself (a spatially-uniform
``\mathbf{u}`` is just a translation, no distortion at all). For two nearby particles separated
by a small vector ``\delta\mathbf{x}`` in the reference configuration, their relative
displacement to first order is
```math
\delta u_i = \frac{\partial u_i}{\partial x_j}\,\delta x_j = J_{ij}\,\delta x_j,
\qquad J_{ij} = \frac{\partial u_i}{\partial x_j},
```
so the two particles end up separated by ``\delta\mathbf{x}+\delta\mathbf{u} =
(I+J)\,\delta\mathbf{x}`` — the deformation gradient ``I+J`` acting on their original
separation. Everything below is about what this matrix ``J`` (the Jacobian of ``\mathbf{u}``)
does to that separation.
"""

# ╔═╡ bbb43caf-32a4-4f70-a263-246953076d84
md"""
## Finite vs. Infinitesimal Strain

The **deformation gradient** ``\mathbf{F}=\mathbf{I}+\nabla_0\mathbf{u}`` (the gradient taken
with respect to the *reference* position) maps a reference line element exactly onto its
deformed image, ``d\mathbf{x}=\mathbf{F}\,d\mathbf{x}_0``, with no restriction on how large the
deformation is. Splitting it into "rotation + stretch" *exactly*, for any size deformation, needs
a **polar decomposition** ``\mathbf{F}=\mathbf{R}\mathbf{U}`` — ``\mathbf{R}`` a genuine (finite)
rotation matrix, ``\mathbf{U}`` a genuine (finite) stretch — and ``\mathbf{R}``, ``\mathbf{U}``
do *not* simply add or subtract from ``\mathbf{F}``.

Seismology almost never needs that. Real elastic strains in the Earth are tiny (the
introduction's "How tiny" note above: ``10^{-8}`` to ``10^{-4}``), so
``\|\nabla\mathbf{u}\|\ll1``, and to first order in that small number:

- The reference and current configurations are close enough that it no longer matters which one
  you differentiate with respect to — ``\nabla_0\mathbf{u}\approx\nabla\mathbf{u}`` — so this
  notebook never has to distinguish them, and just writes ``J=\nabla\mathbf{u}``.
- The exact multiplicative decomposition collapses to the simple **additive** one used
  throughout this notebook, ``J=e+\Omega``, with ``e=\tfrac12(J+J^T)`` symmetric and
  ``\Omega=\tfrac12(J-J^T)`` antisymmetric — the **infinitesimal strain tensor** and the
  **infinitesimal rotation tensor**.
- ``\mathbf{I}+\Omega`` is only the *linearized* rotation matrix (correct to first order in
  ``\Omega``), not a true finite rotation — already visible below in **Dilatation vs.
  Distortion**, where a pure rotation preserves area only to first order, not exactly.

Every result in this notebook — the strain ellipse, the Mohr circle, the dilatation/deviatoric
split — is this infinitesimal theory, valid because real seismic strains stay in that tiny
regime; it is *not* the general finite-strain theory needed for, say, mantle convection over
geologic time.
"""

# ╔═╡ a660a926-0835-4ed2-a2a6-615e3dba3088
md"""
## Strain Tensor
We use a strain tensor to analyze the distortion of the medium, whether it is solid or fluid,
elastic or inelastic. The Jacobian ``J=\nabla\mathbf{u}`` above is completely general — it holds
for any displacement field. The interactive widget below explores one important special case: a
**spatially-uniform** ``J`` (every particle in the square is stretched/rotated by the *same*
constant matrix), which is exactly what a small enough neighborhood of *any* smooth displacement
field looks like.

Any matrix splits uniquely into a symmetric and an antisymmetric part (the **Toeplitz
decomposition**) — the symmetric part is the **strain tensor**, the antisymmetric part the
**(infinitesimal) rotation tensor**:
```math
J = e + \Omega,\qquad
e = \tfrac12\!\left(J+J^T\right)\ \text{(symmetric)},\qquad
\Omega = \tfrac12\!\left(J-J^T\right)\ \text{(antisymmetric)}.
```
Acting separately on the same small separation ``\delta\mathbf{x}``, the two parts do genuinely
different things: ``e\,\delta\mathbf{x}`` stretches/shears ``\delta\mathbf{x}`` with *no*
rotation (it's the part responsible for elastic restoring forces), while
``\Omega\,\delta\mathbf{x}`` is a pure infinitesimal *rotation* of ``\delta\mathbf{x}`` with *no*
change in length at all — a rigid spin that costs a deforming solid nothing. Splitting ``J`` this
way is exactly what lets the widget below decompose a general deformation into rotation,
dilatation, and deviatoric strain.
"""

# ╔═╡ d8f1ceea-1968-4e59-93ae-a21f035b2a09
md"""
## Principal Strains

Just like any symmetric matrix, the strain tensor ``e`` has an orthonormal set of
eigenvectors — the **principal strain axes** — along which the deformation is a pure
stretch/squash with *no* shear at all. In 2-D, the two eigenvalues (principal strains) are

```math
\lambda_{1,2} = \frac{e_{11}+e_{22}}{2} \pm \sqrt{\left(\frac{e_{11}-e_{22}}{2}\right)^2 + e_{12}^2},
```

at an angle ``\theta_p = \tfrac12\operatorname{atan2}(2e_{12},\, e_{11}-e_{22})`` from the
``x``-axis. Geometrically, a unit circle of material points deforms into an **ellipse**
under ``e`` alone, and the principal axes are exactly that ellipse's own major/minor axes —
this is the classic structural-geology **strain ellipse**, drawn on the "Strain only" panel
in the widget above, with the two principal strain values reported in the Readouts panel.
Drag towards **pure rotation** and watch both principal strains go to zero (a rotation has
no eigenvector with real, non-trivial stretch); drag towards **isotropic dilation** and
watch the ellipse become a circle again, just a bigger or smaller one — the ellipse only
grows "lopsided" when the two principal strains differ.
"""

# ╔═╡ 6400c97b-9f07-4f52-9002-e8b5f11aa182
md"""
## Dilatation vs. Distortion

The strain tensor ``e`` itself splits one level further: its **trace**
``\Delta = \operatorname{tr}(e) = e_{11}+e_{22}(+e_{33}) = \nabla\cdot𝐮`` is the
**dilatation** — the fractional change in area (2-D) or volume (3-D), the *same* in every
direction, with zero shape change. Subtracting that isotropic part back out of ``e`` leaves
the **deviatoric strain**: shape change at *exactly* constant area/volume, zero dilatation
by construction. This is the split behind P waves (which involve dilatation — a volume
oscillation) versus S waves (purely deviatoric — a shape oscillation with no volume change
at all): the "Dilatation only" and "Deviatoric only" panels in the widget above show these
in isolation, and the **Δ** readout is exactly ``\operatorname{tr}(e)``.

We can check this claim numerically: deform a unit square by a purely isotropic matrix
``J_\Delta = \begin{bmatrix}j_\Delta&0\\0&j_\Delta\end{bmatrix}``, by a general symmetric
(strain) matrix, and by a pure rotation, and compare their areas.
"""

# ╔═╡ a759772e-a5b4-4f35-9150-c17627b58c4a
corners = [[0, 0], [1, 0], [1, 1], [0, 1]]

# ╔═╡ f4b45b1b-b7f8-4f2b-a580-98def9cc8fbc
md"""
## Mohr Circle for Strain

The Mohr circle is a graphical trick for reading off the normal and shear strain on a cut
plane at *any* orientation, without recomputing the tensor-rotation formula by hand each
time. Rotating the cut angle by ``\theta`` moves you *twice* as far around the circle
(``2\theta``) — a consequence of strain being a rank-2 tensor, not a vector. The circle is
centered at ``(\Delta/2,\, 0)`` with radius
``R=\sqrt{\left(\tfrac{e_{11}-e_{22}}{2}\right)^2+e_{12}^2}``; the two points where it
crosses the (zero-shear) horizontal axis are exactly the principal strains from the section
above.

Drag the marker on the Mohr-circle panel in the widget: the matching cut-plane line on the
main square updates with it — a dashed line at the *reference* (undeformed) orientation, and
a solid cyan line showing how that same material fiber looks after the full deformation
(rotation included). The Mohr circle itself deliberately plots only the **strain** part
``e``; the extra rotation visible on the main square's cyan line is *not* in the circle's
numbers — a fiber that looks like it swung around a lot on the square can still show almost
no shear strain on the circle, if most of what you're seeing is rigid rotation rather than
distortion.
"""

# ╔═╡ 31eb56b0-25dd-4335-bc17-f89431d79b22
md"""
## Earthquake Recurrence

Dividing the coseismic strain drop of a major earthquake by the tectonic strain *rate*
loading the fault back up gives a rough estimate of the **recurrence interval** — how long
it takes to build up enough strain for another rupture of similar size. Drag the sliders
below (log-scaled, since both quantities span orders of magnitude) around the real Southern
California numbers quoted in the introduction.
"""

# ╔═╡ a4ec2280-8bac-4e65-a849-909441cee0b6
md"""
Strain rate ``\dot\epsilon`` (/year): $(@bind log_rate Slider(-8:0.05:-6, default=-6.5, show_value=false))

Coseismic strain drop ``\Delta\epsilon`` (per event): $(@bind log_drop Slider(-5:0.05:-3, default=-3.5, show_value=false))
"""

# ╔═╡ 6f33c7ed-e8a8-440d-a771-84789ee1f397
md"""
#### Appendix
"""

# ╔═╡ ca079784-f5ac-4e9b-9d10-8501564c47f0
md"### Layer 1: Area-Change Helper"

# ╔═╡ 962355bc-2e2e-4bba-b293-aefdca8f3627
"""
    calculate_area(corner_vectors)

Shoelace-formula area of the quadrilateral with the 4 given corner points (in order around
the perimeter) — used above to check, numerically, exactly which part of a deformation
(dilatation, deviatoric strain, or rotation) actually changes a material element's area.
Splits the quadrilateral into the two triangles `ABC` and `ACD` sharing diagonal `AC` and
sums their (signed) areas via the cross product.
"""
function calculate_area(corner_vectors)
    A = corner_vectors[1]
    B = corner_vectors[2]
    C = corner_vectors[3]
    D = corner_vectors[4]

    area = 0.5 * abs((B[1] - A[1]) * (C[2] - A[2]) - (C[1] - A[1]) * (B[2] - A[2]) +
                     (C[1] - A[1]) * (D[2] - A[2]) - (D[1] - A[1]) * (C[2] - A[2]))
end

# ╔═╡ 9eeeac07-c2f6-4759-a328-4dee9d576723
let
    je1, je2, je, jd, jom = 0.35, -0.15, 0.22, 0.18, 0.12
    Je = [je1 je; je je2]      # a general symmetric (strain) matrix
    Jdil = [jd 0; 0 jd]        # a purely isotropic (dilatation-only) matrix
    Jom = [0 jom; -jom 0]      # a general antisymmetric (rotation) matrix
    A_strain = calculate_area(map(c -> Je * c + c, corners))
    A_dilatation = calculate_area(map(c -> Jdil * c + c, corners))
    A_rotation = calculate_area(map(c -> Jom * c + c, corners))
    md"""
    Plugging in one representative numeric case (``j_{e1}=$(je1)``, ``j_{e2}=$(je2)``,
    ``j_e=$(je)``, ``j_\Delta=$(jd)``, ``j_\omega=$(jom)``):

    - **general symmetric (strain)** ``J_e``: computed area = **$(round(A_strain, digits=5))**,
      *not* simply ``1+j_{e1}+j_{e2}=$(round(1+je1+je2, digits=5))`` — it depends on both the
      dilatation *and* the deviatoric shear ``j_e``.
    - **purely isotropic (dilatation-only)** ``J_\Delta``: computed area =
      **$(round(A_dilatation, digits=5))**, matching
      ``(1+j_\Delta)^2=$(round((1+jd)^2, digits=5))`` exactly — confirming dilatation alone
      controls the area change.
    - **pure (infinitesimal) rotation** ``J_\omega``: computed area =
      **$(round(A_rotation, digits=5))**, matching ``1+j_\omega^2=$(round(1+jom^2, digits=5))``
      — equal to `1` only to *first order* in ``j_\omega``. That leftover ``j_\omega^2`` is
      exactly the reminder from **Finite vs. Infinitesimal Strain** above that ``I+\Omega`` is
      only the *linearized* rotation matrix, not a finite one — an exact finite rotation
      ``\begin{bmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{bmatrix}`` would give
      area `1` identically. Either way, rotation contributes nothing to first order, which is
      the only order that matters for infinitesimal strain.
    """
end

# ╔═╡ 3f45f0c5-e672-4780-a2c6-006e49fecfb8
md"### Layer 2: Earthquake Recurrence"

# ╔═╡ 3f7a0102-597e-49cc-978b-39ad7fcb3015
recurrence_years(coseismic_drop, strain_rate) = coseismic_drop / strain_rate

# ╔═╡ c098aa8d-6106-4bf4-b6c4-b2c75e9225a2
let
    rate = 10.0^log_rate
    drop = 10.0^log_drop
    T = recurrence_years(drop, rate)
    md"""
    - strain rate ``\dot\epsilon \approx`` **$(round(rate, sigdigits=2))** per year
    - coseismic drop ``\Delta\epsilon \approx`` **$(round(drop, sigdigits=2))**
    - estimated recurrence interval ``T = \Delta\epsilon / \dot\epsilon \approx``
      **$(round(T, sigdigits=3)) years**
    """
end

# ╔═╡ 6b829564-d656-4c81-885f-212d1b8d6d91
md"### The Configuration Widget"

# ╔═╡ 0168e21c-7af1-491d-be1b-7b3f458f9ad3
begin
    """
        ConfigurationInput(; rho0=1.0)

    A body's **reference configuration** ℛ₀ (dashed, fixed) and its **deformed configuration**
    ℛ (solid, draggable) shown at once. The four corner handles carry a persistent numeric
    label attached to *that* material point — dragging a handle moves the deformed position
    `x` of the *same* labeled particle, never relabels which particle is which. That labeling
    is the Lagrangian description: properties are tracked by which material point they belong
    to, not by which fixed spatial location currently holds them.

    Mass is an **extensive** property: `M = ρ₀·A₀` is fixed by the reference configuration and
    the `ρ₀` slider alone, unaffected by how the region is subsequently dragged/deformed
    (that's what "extensive" means — conservation of mass under the mapping ℛ₀→ℛ). Density
    `ρ = M/A` is **intensive**: it changes as the region's *current* area `A` changes, exactly
    compensating so `ρ·A = ρ₀·A₀` always holds.
    """
    struct ConfigurationInput
        rho0::Float64
    end

    ConfigurationInput(; rho0=1.0) = ConfigurationInput(Float64(rho0))

    Base.get(w::ConfigurationInput) = Dict{String,Any}("rho0" => w.rho0)

    function Base.show(io::IO, ::MIME"text/html", w::ConfigurationInput)
        write(io, """
<div id="cfgwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#cfgwidget){width:min(70vw,1100px)!important;margin-left:calc((100% - min(70vw,1100px))/2)!important}
    #cfgwidget .cfg-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #cfgwidget .cfg-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #cfgwidget .cfg-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #cfgwidget canvas{cursor:crosshair;background:#000;border:1px solid #374151;border-radius:6px;display:block;max-width:100%}
    #cfgwidget .cfg-controls{display:flex;gap:10px;flex-wrap:wrap;width:100%;margin-top:12px}
    #cfgwidget .cfg-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px;flex:1 1 220px;min-width:220px}
    #cfgwidget .cfg-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #cfgwidget .cfg-readout{font-size:13px;line-height:1.8}
    #cfgwidget .cfg-readout b{font-variant-numeric:tabular-nums}
    #cfgwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px;cursor:pointer}
    #cfgwidget input[type=range]{width:100%}
  </style>
  <div class="cfg-title">
    <div class="cfg-title-desc">A body's part (region) has a reference configuration and a deformed one &mdash; drag the labeled particles to deform it.</div>
    <div class="cfg-title-hint">drag a numbered handle &middot; labels are Lagrangian tags, fixed to the particle &middot; ρ&#8320; slider sets the reference density</div>
  </div>
  <canvas id="cfgMain" width="420" height="420"></canvas>
  <div class="cfg-controls">
    <div class="cfg-control-group">
      <div class="cfg-control-title">Reference Density</div>
      <label>ρ&#8320; <input type="range" id="cfg-rho0" min="0.2" max="5" step="0.1" value="$(w.rho0)"><span id="cfg-rho0-v">$(w.rho0)</span></label>
      <div style="margin-top:8px"><button id="cfg-reset" type="button">Reset to reference</button></div>
    </div>
    <div class="cfg-control-group" style="flex:1 1 260px">
      <div class="cfg-control-title">Readouts</div>
      <div class="cfg-readout" id="cfg-readout"></div>
    </div>
    <div class="cfg-control-group" style="flex:1 1 260px">
      <div class="cfg-control-title">Legend</div>
      <div class="cfg-readout">
        <span style="color:#6b7280">dashed</span> = reference configuration ℛ&#8320; (fixed)<br>
        <span style="color:#38bdf8">solid</span> = deformed configuration ℛ (drag it)<br>
        <span style="color:#facc15">1..4</span> = Lagrangian labels, attached to their own particle<br>
        <b>M</b> = ρ&#8320;·A&#8320; is extensive (fixed) &middot; <b>ρ</b> = M/A is intensive (changes)
      </div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.6, par.clientWidth || window.innerWidth*0.6)
  const SEC = Math.round(Math.min(availW, 460, Math.max(280, window.innerHeight-560)))
  const DPR = Math.min(window.devicePixelRatio || 1, 2)
  const cvs = par.querySelector('#cfgMain'), ctx = cvs.getContext('2d')
  cvs.width = Math.round(SEC*DPR); cvs.height = Math.round(SEC*DPR)
  cvs.style.width = SEC+'px'; cvs.style.height = SEC+'px'
  ctx.setTransform(DPR,0,0,DPR,0,0)

  // Fixed logical data range -- the canvas shows a window onto a much larger continuum, of
  // which the labeled quadrilateral below is only "a part of the body".
  const VMIN = 0, VMAX = 4
  const to = (x,y) => [(x-VMIN)/(VMAX-VMIN)*SEC, SEC-(y-VMIN)/(VMAX-VMIN)*SEC]
  const invTo = (px,py) => [px/SEC*(VMAX-VMIN)+VMIN, (SEC-py)/SEC*(VMAX-VMIN)+VMIN]

  // Reference configuration ℛ₀ -- fixed for the lifetime of the widget. The deformed
  // configuration ℛ starts out coincident with it (undeformed) and is what dragging changes.
  const refPts = [[1,1],[3,1],[3,3],[1,3]]
  let curPts = refPts.map(p=>p.slice())
  let rho0 = $(w.rho0)
  let dragIdx = -1

  function shoelaceArea(pts){
    let a = 0
    for(let i=0;i<pts.length;i++){
      const [x1,y1] = pts[i], [x2,y2] = pts[(i+1)%pts.length]
      a += x1*y2 - x2*y1
    }
    return Math.abs(a)/2
  }
  const A0 = shoelaceArea(refPts)

  function drawLattice(){
    ctx.strokeStyle = 'rgba(156,163,175,0.18)'; ctx.lineWidth = 1
    for(let i=0;i<=8;i++){
      const t = VMIN + i*(VMAX-VMIN)/8
      let a=to(t,VMIN), b=to(t,VMAX); ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke()
      a=to(VMIN,t); b=to(VMAX,t); ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke()
    }
  }
  function drawPoly(pts, color, dashed){
    ctx.setLineDash(dashed ? [5,4] : [])
    ctx.strokeStyle = color; ctx.lineWidth = dashed ? 1.5 : 2.5
    ctx.beginPath()
    pts.forEach((p,i)=>{ const q=to(p[0],p[1]); i===0?ctx.moveTo(q[0],q[1]):ctx.lineTo(q[0],q[1]) })
    ctx.closePath(); ctx.stroke()
    ctx.setLineDash([])
  }

  function drawAll(){
    ctx.clearRect(0,0,SEC,SEC); ctx.fillStyle='#000'; ctx.fillRect(0,0,SEC,SEC)
    drawLattice()
    drawPoly(refPts, '#6b7280', true)
    drawPoly(curPts, '#38bdf8', false)
    curPts.forEach((p,i)=>{
      const q = to(p[0],p[1])
      ctx.beginPath(); ctx.arc(q[0],q[1],8,0,2*Math.PI)
      ctx.fillStyle = '#facc15'; ctx.fill()
      ctx.strokeStyle = '#0a0f18'; ctx.lineWidth = 1.5; ctx.stroke()
      ctx.fillStyle = '#0a0f18'; ctx.font = 'bold 12px sans-serif'; ctx.textAlign='center'; ctx.textBaseline='middle'
      ctx.fillText(String(i+1), q[0], q[1]+1)
      ctx.textAlign='left'; ctx.textBaseline='alphabetic'
    })
    const A = shoelaceArea(curPts)
    const M = rho0*A0
    const rho = M/A
    par.querySelector('#cfg-readout').innerHTML =
      'reference area A&#8320; = <b>'+A0.toFixed(2)+'</b><br>'+
      'current area A = <b>'+A.toFixed(2)+'</b><br>'+
      'mass M = ρ&#8320;·A&#8320; = <b>'+M.toFixed(2)+'</b> <span style="color:#6b7280">(fixed)</span><br>'+
      'density ρ = M/A = <b>'+rho.toFixed(2)+'</b> <span style="color:#38bdf8">(live)</span>'
  }
  drawAll()

  cvs.addEventListener('mousedown', e=>{
    const mx=e.offsetX, my=e.offsetY
    dragIdx = -1
    for(let i=0;i<curPts.length;i++){
      const q = to(curPts[i][0], curPts[i][1])
      if(Math.hypot(q[0]-mx,q[1]-my) < 14){ dragIdx = i; break }
    }
  })
  cvs.addEventListener('mousemove', e=>{
    if(dragIdx<0) return
    const [x,y] = invTo(e.offsetX, e.offsetY)
    curPts[dragIdx] = [x,y]
    drawAll()
  })
  window.addEventListener('mouseup', ()=>{ dragIdx = -1 })

  par.querySelector('#cfg-rho0').addEventListener('input', e=>{
    rho0 = parseFloat(e.target.value)
    par.querySelector('#cfg-rho0-v').textContent = rho0.toFixed(1)
    drawAll()
  })
  par.querySelector('#cfg-reset').addEventListener('click', ()=>{
    curPts = refPts.map(p=>p.slice())
    drawAll()
  })
</script>
""")
    end

    const _cfg_ready = true
end

# ╔═╡ 11966694-3009-4016-be5b-2a4eb1be6be2
begin
    _cfg_ready
    @bind cfg ConfigurationInput()
end

# ╔═╡ 906d8e1c-67c5-46ca-ad8d-0e087561690c
md"### The Interactive Widget"

# ╔═╡ e2c19bc6-e3df-4c93-aded-b9fda77247ec
begin
    """
        StrainTensorInput(; j11=0.9, j12=0.1, j21=0.4, j22=0.1, theta=0.0)

    A draggable 2-D Jacobian: drag either handle on the main square to set the deformation
    directly (handle A sets column 1 of `J` — where the `x`-basis vector goes — handle B sets
    column 2), and watch it decompose live into rotation, dilatation (uniform area change,
    `Δ = tr(e)`), and deviatoric strain (pure shape change at constant area) in the comparison
    strip below. The right-hand panel is the Mohr circle for the current strain `e`: drag its
    angle handle to pick a cut-plane orientation `θ` and read off the normal/shear strain on
    that plane, with the matching cut-plane line drawn on the main square. Presets jump to
    canonical cases (pure shear, simple shear, pure rotation, isotropic dilation).
    """
    struct StrainTensorInput
        j11::Float64
        j12::Float64
        j21::Float64
        j22::Float64
        theta::Float64
    end

    StrainTensorInput(; j11=0.9, j12=0.1, j21=0.4, j22=0.1, theta=0.0) =
        StrainTensorInput(Float64(j11), Float64(j12), Float64(j21), Float64(j22), Float64(theta))

    Base.get(w::StrainTensorInput) = Dict{String,Any}(
        "j11" => w.j11, "j12" => w.j12, "j21" => w.j21, "j22" => w.j22, "theta" => w.theta,
    )

    function Base.show(io::IO, ::MIME"text/html", w::StrainTensorInput)
        write(io, """
<div id="stwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#stwidget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #stwidget .st-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #stwidget .st-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #stwidget .st-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #stwidget .st-panels{display:flex;gap:14px;align-items:flex-start;justify-content:center;width:100%}
    #stwidget .st-col{flex:1 1 0;min-width:0;display:flex;flex-direction:column;align-items:center}
    #stwidget .st-panel-label{font-size:12px;color:#6b7280;margin-bottom:4px;text-transform:uppercase;letter-spacing:.04em}
    #stwidget canvas{cursor:crosshair;background:#000;border:1px solid #374151;border-radius:6px;display:block;max-width:100%}
    #stwidget .st-strip{display:flex;gap:10px;align-items:flex-start;justify-content:center;width:100%;margin-top:14px;flex-wrap:wrap}
    #stwidget .st-strip-item{display:flex;flex-direction:column;align-items:center}
    #stwidget .st-strip-label{font-size:11px;color:#9ca3af;margin-bottom:3px;text-align:center}
    #stwidget .st-controls{display:flex;gap:10px;flex-wrap:wrap;width:100%;margin-top:12px}
    #stwidget .st-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px;flex:1 1 220px;min-width:220px}
    #stwidget .st-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #stwidget select{width:100%;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:5px 6px;font-size:14px}
    #stwidget .st-readout{font-size:13px;line-height:1.7}
    #stwidget .st-mat{font-variant-numeric:tabular-nums;white-space:pre}
    #stwidget .st-jcol{flex:0 0 auto;padding-top:22px}
    #stwidget .st-jmatrix{display:inline-grid;grid-template-columns:auto auto;gap:6px 8px;padding:6px 12px;position:relative}
    #stwidget .st-jmatrix::before, #stwidget .st-jmatrix::after{content:'';position:absolute;top:3px;bottom:3px;width:7px;border:2px solid #9ca3af}
    #stwidget .st-jmatrix::before{left:0;border-right:none}
    #stwidget .st-jmatrix::after{right:0;border-left:none}
    #stwidget .st-jmatrix input{width:58px;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:5px 4px;font-size:14px;text-align:center;font-variant-numeric:tabular-nums}
    #stwidget .st-jmatrix input:focus{outline:2px solid #38bdf8;border-color:#38bdf8}
    #stwidget .st-jhint{font-size:11px;color:#6b7280;margin-top:10px;text-align:center;max-width:150px}
    #stwidget .st-mat2{display:inline-grid;grid-template-columns:auto auto;gap:2px 6px;padding:3px 9px;position:relative;vertical-align:middle;margin:2px 0}
    #stwidget .st-mat2::before, #stwidget .st-mat2::after{content:'';position:absolute;top:1px;bottom:1px;width:5px;border:1.5px solid currentColor}
    #stwidget .st-mat2::before{left:0;border-right:none}
    #stwidget .st-mat2::after{right:0;border-left:none}
    #stwidget .st-mat2 span{min-width:38px;text-align:center;font-size:13px;font-variant-numeric:tabular-nums}
    @media (max-width: 900px){
      #stwidget .st-panels{flex-direction:column;align-items:center}
      #stwidget .st-col{width:100%;max-width:520px}
      #stwidget .st-jcol{padding-top:0}
    }
  </style>
  <div class="st-title">
    <div class="st-title-desc">A uniform deformation splits exactly into rotation, dilatation, and shape-changing shear &mdash; drag to see how much of each.</div>
    <div class="st-title-hint">drag a corner handle, or edit J directly, to deform &middot; drag the Mohr circle's marker to pick a cut plane &middot; pick a preset below</div>
  </div>
  <div class="st-panels">
    <div class="st-col">
      <div class="st-panel-label">Deformation (drag A / B)</div>
      <canvas id="stMain" width="480" height="480"></canvas>
    </div>
    <div class="st-col st-jcol">
      <div class="st-panel-label">Displacement Gradient J</div>
      <div class="st-jmatrix">
        <input type="number" step="0.05" id="st-j11" value="$(w.j11)">
        <input type="number" step="0.05" id="st-j12" value="$(w.j12)">
        <input type="number" step="0.05" id="st-j21" value="$(w.j21)">
        <input type="number" step="0.05" id="st-j22" value="$(w.j22)">
      </div>
      <div class="st-jhint">edit directly, or drag A/B on the square &mdash; both stay in sync</div>
    </div>
    <div class="st-col">
      <div class="st-panel-label">Mohr circle for strain (drag the marker)</div>
      <canvas id="stMohr" width="480" height="480"></canvas>
    </div>
  </div>
  <div class="st-strip">
    <div class="st-strip-item"><div class="st-strip-label" style="color:#38bdf8">Strain only (e)</div><canvas id="stStrain" width="190" height="190"></canvas></div>
    <div class="st-strip-item"><div class="st-strip-label" style="color:#c084fc">Rotation only (&Omega;)</div><canvas id="stRot" width="190" height="190"></canvas></div>
    <div class="st-strip-item"><div class="st-strip-label" style="color:#facc15">Dilatation only</div><canvas id="stDil" width="190" height="190"></canvas></div>
    <div class="st-strip-item"><div class="st-strip-label" style="color:#4ade80">Deviatoric only</div><canvas id="stDev" width="190" height="190"></canvas></div>
  </div>
  <div class="st-controls">
    <div class="st-control-group">
      <div class="st-control-title">Preset</div>
      <select id="stPreset">
        <option value="custom">Custom</option>
        <option value="shear">Pure shear</option>
        <option value="simple">Simple shear</option>
        <option value="rotation">Pure rotation</option>
        <option value="dilation">Isotropic dilation</option>
      </select>
    </div>
    <div class="st-control-group" style="flex:1 1 260px">
      <div class="st-control-title">Readouts</div>
      <div class="st-readout" id="stReadout"></div>
    </div>
    <div class="st-control-group" style="flex:1 1 260px">
      <div class="st-control-title">Legend</div>
      <div class="st-readout" id="stLegend"></div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8)
  const heightBudget = Math.max(320, window.innerHeight - 520)
  // A third, narrow column (the editable J matrix) now sits between the two square canvases --
  // two gaps instead of one, and its own fixed width, come out of the same available space.
  const JCOLW = 170
  const SEC = Math.round(Math.min((availW-28-JCOLW)/2, heightBudget, 480))
  const STRIP = Math.round(Math.min((availW-30)/4, 190))
  const DPR = Math.min(window.devicePixelRatio || 1, 2)

  function hidpi(cv, cx, w, h){ cv.width=Math.round(w*DPR); cv.height=Math.round(h*DPR); cv.style.width=w+'px'; cv.style.height=h+'px'; cx.setTransform(DPR,0,0,DPR,0,0) }

  const mainCvs = par.querySelector('#stMain'), mainCtx = mainCvs.getContext('2d')
  const mohrCvs = par.querySelector('#stMohr'), mohrCtx = mohrCvs.getContext('2d')
  const strainCvs = par.querySelector('#stStrain'), strainCtx = strainCvs.getContext('2d')
  const rotCvs = par.querySelector('#stRot'), rotCtx = rotCvs.getContext('2d')
  const dilCvs = par.querySelector('#stDil'), dilCtx = dilCvs.getContext('2d')
  const devCvs = par.querySelector('#stDev'), devCtx = devCvs.getContext('2d')
  hidpi(mainCvs, mainCtx, SEC, SEC)
  hidpi(mohrCvs, mohrCtx, SEC, SEC)
  hidpi(strainCvs, strainCtx, STRIP, STRIP)
  hidpi(rotCvs, rotCtx, STRIP, STRIP)
  hidpi(dilCvs, dilCtx, STRIP, STRIP)
  hidpi(devCvs, devCtx, STRIP, STRIP)

  // Fixed logical data range shared by the main canvas and the four strip canvases -- a
  // generous window around the unit square so handles can be dragged well past it.
  const VMIN = -2.2, VMAX = 3.2
  function xform(sizePx){
    const scale = sizePx/(VMAX-VMIN)
    return (x,y) => [(x-VMIN)*scale, sizePx-(y-VMIN)*scale]
  }
  function invxform(sizePx){
    const scale = sizePx/(VMAX-VMIN)
    return (px,py) => [px/scale+VMIN, (sizePx-py)/scale+VMIN]
  }

  // Fixed axis range for the Mohr circle too -- the circle itself moves/grows/shrinks as J
  // changes, rather than the axes continuously rescaling to hug it. A rescaling axis makes it
  // hard to tell "did the circle actually change, or did the frame just rezoom" at a glance;
  // a fixed frame makes size/position changes directly readable. ±2.6 comfortably covers the
  // circle's reachable center/radius across the full A/B drag range (VMIN..VMAX above), with
  // headroom to spare for the four small presets; a circle from a hand-typed extreme J value
  // can legitimately run off the fixed frame, same as dragging A/B off the main canvas.
  const MOHR_HALF = 2.6
  function mohrXform(x,y){
    return [(x+MOHR_HALF)/(2*MOHR_HALF)*SEC, SEC-(y+MOHR_HALF)/(2*MOHR_HALF)*SEC]
  }
  function mohrInvXform(px,py){
    return [(px/SEC)*(2*MOHR_HALF)-MOHR_HALF, ((SEC-py)/SEC)*(2*MOHR_HALF)-MOHR_HALF]
  }

  let j11 = $(w.j11), j12 = $(w.j12), j21 = $(w.j21), j22 = $(w.j22)
  let theta = $(w.theta)
  let currentPreset = 'custom'
  let dragMode = null

  const PRESETS = {
    shear:    {label:'Pure shear',          j11:0.5,  j12:0.0,  j21:0.0,  j22:-0.5},
    simple:   {label:'Simple shear',        j11:0.0,  j12:0.6,  j21:0.0,  j22:0.0},
    rotation: {label:'Pure rotation',       j11:0.0,  j12:-0.4, j21:0.4,  j22:0.0},
    dilation: {label:'Isotropic dilation',  j11:0.35, j12:0.0,  j21:0.0,  j22:0.35},
  }

  // J = e + Ω (Toeplitz decomposition), then e = dilatation (Δ/2)·I + deviatoric remainder.
  function decompose(){
    const e11=j11, e22=j22, e12=(j12+j21)/2
    const delta = e11+e22
    return {
      e:   {m11:e11, m12:e12, m21:e12, m22:e22},
      Om:  {m11:0, m12:(j12-j21)/2, m21:(j21-j12)/2, m22:0},
      delta,
      dil: {m11:delta/2, m12:0, m21:0, m22:delta/2},
      dev: {m11:e11-delta/2, m12:e12, m21:e12, m22:e22-delta/2},
    }
  }
  function deform(M, x, y){ return [x + M.m11*x + M.m12*y, y + M.m21*x + M.m22*y] }
  function principal(e){
    const c=(e.m11+e.m22)/2, A=(e.m11-e.m22)/2, R=Math.hypot(A, e.m12)
    const angle = 0.5*Math.atan2(2*e.m12, e.m11-e.m22)
    return {l1:c+R, l2:c-R, angle, c, R}
  }

  function drawGridDeform(ctx, sizePx, M, opts){
    const to = xform(sizePx)
    const N = opts.n || 6
    ctx.clearRect(0,0,sizePx,sizePx)
    // faint reference (undeformed) unit square grid
    ctx.strokeStyle = 'rgba(156,163,175,0.35)'; ctx.lineWidth = 1
    for(let i=0;i<=N;i++){
      const t=i/N
      let a=to(t,0), b=to(t,1); ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke()
      a=to(0,t); b=to(1,t); ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke()
    }
    // deformed grid
    ctx.strokeStyle = opts.color || '#e5e7eb'; ctx.lineWidth = opts.lineWidth || 2
    for(let i=0;i<=N;i++){
      const t=i/N
      ctx.beginPath()
      for(let k=0;k<=N;k++){ const s=k/N; const [x,y]=deform(M,t,s); const p=to(x,y); k===0?ctx.moveTo(p[0],p[1]):ctx.lineTo(p[0],p[1])}
      ctx.stroke()
      ctx.beginPath()
      for(let k=0;k<=N;k++){ const s=k/N; const [x,y]=deform(M,s,t); const p=to(x,y); k===0?ctx.moveTo(p[0],p[1]):ctx.lineTo(p[0],p[1])}
      ctx.stroke()
    }
  }

  function drawOrigin(ctx, sizePx){
    const to = xform(sizePx)
    const o = to(0,0)
    ctx.beginPath(); ctx.arc(o[0],o[1],3,0,2*Math.PI); ctx.fillStyle='#9ca3af'; ctx.fill()
  }

  function handlePos(sizePx){
    const to = xform(sizePx)
    const A = deform({m11:j11,m12:j12,m21:j21,m22:j22}, 1, 0)
    const B = deform({m11:j11,m12:j12,m21:j21,m22:j22}, 0, 1)
    return { A: to(A[0],A[1]), B: to(B[0],B[1]) }
  }

  function drawHandle(ctx, pos, label, color){
    ctx.beginPath(); ctx.arc(pos[0], pos[1], 7, 0, 2*Math.PI)
    ctx.fillStyle = color; ctx.fill()
    ctx.strokeStyle = '#0a0f18'; ctx.lineWidth = 1.5; ctx.stroke()
    ctx.fillStyle = '#e5e7eb'; ctx.font = 'bold 13px sans-serif'
    ctx.fillText(label, pos[0]+10, pos[1]-8)
  }

  function drawCutLine(ctx, sizePx){
    const to = xform(sizePx)
    const cx=0.5, cy=0.5, L=0.9
    const dx=Math.cos(theta)*L, dy=Math.sin(theta)*L
    // reference (undeformed) orientation, dashed
    let p0=to(cx-dx,cy-dy), p1=to(cx+dx,cy+dy)
    ctx.setLineDash([4,3]); ctx.strokeStyle='rgba(245,243,239,0.55)'; ctx.lineWidth=1.5
    ctx.beginPath(); ctx.moveTo(p0[0],p0[1]); ctx.lineTo(p1[0],p1[1]); ctx.stroke()
    ctx.setLineDash([])
    // how that same material line looks after the full deformation J
    const M = {m11:j11,m12:j12,m21:j21,m22:j22}
    const q0c = deform(M, cx-dx, cy-dy), q1c = deform(M, cx+dx, cy+dy)
    const q0 = to(q0c[0], q0c[1]), q1 = to(q1c[0], q1c[1])
    ctx.strokeStyle = '#22d3ee'; ctx.lineWidth = 2.5
    ctx.beginPath(); ctx.moveTo(q0[0],q0[1]); ctx.lineTo(q1[0],q1[1]); ctx.stroke()
  }

  function drawEllipseAndAxes(ctx, sizePx, e, pr=principal(e)){
    const to = xform(sizePx)
    ctx.beginPath()
    for(let i=0;i<=64;i++){
      const a = 2*Math.PI*i/64
      const [x,y] = deform(e, Math.cos(a)*0.5+0.5, Math.sin(a)*0.5+0.5)
      const p = to(x,y)
      i===0 ? ctx.moveTo(p[0],p[1]) : ctx.lineTo(p[0],p[1])
    }
    ctx.closePath()
    ctx.strokeStyle = 'rgba(56,189,248,0.9)'; ctx.lineWidth = 2; ctx.stroke()
    // principal axes through the ellipse's own center (0.5,0.5)
    const cx=0.5, cy=0.5
    for(const [lam, sign] of [[pr.l1,1],[pr.l2,-1]]){
      const ang = pr.angle + (sign<0 ? Math.PI/2 : 0)
      const r = 0.5*(1+lam)
      const p0 = to(cx-Math.cos(ang)*r, cy-Math.sin(ang)*r)
      const p1 = to(cx+Math.cos(ang)*r, cy+Math.sin(ang)*r)
      ctx.beginPath(); ctx.moveTo(p0[0],p0[1]); ctx.lineTo(p1[0],p1[1])
      ctx.strokeStyle = 'rgba(245,243,239,0.7)'; ctx.lineWidth = 1.5; ctx.stroke()
    }
  }

  function drawMain(){
    const M = {m11:j11,m12:j12,m21:j21,m22:j22}
    drawGridDeform(mainCtx, SEC, M, {color:'#e5e7eb', lineWidth:2.2})
    drawCutLine(mainCtx, SEC)
    drawOrigin(mainCtx, SEC)
    const h = handlePos(SEC)
    drawHandle(mainCtx, h.A, 'A', '#facc15')
    drawHandle(mainCtx, h.B, 'B', '#f472b6')
  }

  function drawStrip(state){
    const {d, pr} = state
    drawGridDeform(strainCtx, STRIP, d.e,  {color:'#38bdf8', lineWidth:1.6, n:5})
    drawEllipseAndAxes(strainCtx, STRIP, d.e, pr)
    drawGridDeform(rotCtx,    STRIP, d.Om, {color:'#c084fc', lineWidth:1.6, n:5})
    drawGridDeform(dilCtx,    STRIP, d.dil,{color:'#facc15', lineWidth:1.6, n:5})
    drawGridDeform(devCtx,    STRIP, d.dev,{color:'#4ade80', lineWidth:1.6, n:5})
  }

  function mohrGeom(e){
    const A=(e.m11-e.m22)/2, B=e.m12, c=(e.m11+e.m22)/2, R=Math.hypot(A,B)
    return {A,B,c,R}
  }
  function mohrPoint(e, th){
    const {A,B,c} = mohrGeom(e)
    const x = c + A*Math.cos(2*th) + B*Math.sin(2*th)
    const y = -A*Math.sin(2*th) + B*Math.cos(2*th)
    return [x,y]
  }

  function drawMohr(state){
    const {d} = state
    const {A,B,c,R} = mohrGeom(d.e)
    mohrCtx.clearRect(0,0,SEC,SEC)
    const to = (x,y) => mohrXform(x,y)

    // axes -- fixed at the frame's own origin, independent of the circle's own center/radius
    mohrCtx.strokeStyle = 'rgba(156,163,175,0.4)'; mohrCtx.lineWidth = 1
    let a0=to(-MOHR_HALF,0), a1=to(MOHR_HALF,0); mohrCtx.beginPath(); mohrCtx.moveTo(a0[0],a0[1]); mohrCtx.lineTo(a1[0],a1[1]); mohrCtx.stroke()
    a0=to(0,-MOHR_HALF); a1=to(0,MOHR_HALF); mohrCtx.beginPath(); mohrCtx.moveTo(a0[0],a0[1]); mohrCtx.lineTo(a1[0],a1[1]); mohrCtx.stroke()

    // the circle itself
    mohrCtx.beginPath()
    for(let i=0;i<=100;i++){
      const th = 2*Math.PI*i/100
      const p = to(c+R*Math.cos(th), R*Math.sin(th))
      i===0 ? mohrCtx.moveTo(p[0],p[1]) : mohrCtx.lineTo(p[0],p[1])
    }
    mohrCtx.strokeStyle = '#38bdf8'; mohrCtx.lineWidth = 2; mohrCtx.stroke()

    // center + principal-strain intercepts
    const ctr = to(c,0)
    mohrCtx.beginPath(); mohrCtx.arc(ctr[0],ctr[1],3,0,2*Math.PI); mohrCtx.fillStyle='#9ca3af'; mohrCtx.fill()
    for(const lam of [c+R, c-R]){
      const p = to(lam,0)
      mohrCtx.beginPath(); mohrCtx.arc(p[0],p[1],4,0,2*Math.PI); mohrCtx.fillStyle='#f5f3ef'; mohrCtx.fill()
    }

    // current marker at angle theta, radius line from center
    const mp = mohrPoint(d.e, theta)
    const mpx = to(mp[0], mp[1])
    mohrCtx.beginPath(); mohrCtx.moveTo(ctr[0],ctr[1]); mohrCtx.lineTo(mpx[0],mpx[1])
    mohrCtx.strokeStyle = 'rgba(34,211,238,0.7)'; mohrCtx.lineWidth = 1.5; mohrCtx.stroke()
    mohrCtx.beginPath(); mohrCtx.arc(mpx[0],mpx[1],7,0,2*Math.PI)
    mohrCtx.fillStyle = '#22d3ee'; mohrCtx.fill()
    mohrCtx.strokeStyle = '#0a0f18'; mohrCtx.lineWidth = 1.5; mohrCtx.stroke()

    mohrCtx.fillStyle = '#9ca3af'; mohrCtx.font = '12px sans-serif'
    mohrCtx.fillText('normal strain →', 10, SEC-10)
    mohrCtx.save(); mohrCtx.translate(14,20); mohrCtx.fillText('↑ shear strain', 0, 0); mohrCtx.restore()
  }

  function fmt(v){ return (v>=0?'+':'') + v.toFixed(2) }
  // Real HTML bracket-grid matrix (reuses the same visual trick as the editable .st-jmatrix
  // panel above, read-only here) instead of a monospace ASCII-bracket string -- `color` tints
  // the brackets to match the matrix's own strip-panel accent, so e/Ω stay visually associated
  // with the "Strain only"/"Rotation only" panels at a glance.
  function matHTML(M, color){
    return '<span class="st-mat2" style="color:'+color+'">'+
      '<span style="color:#e5e7eb">'+fmt(M.m11)+'</span><span style="color:#e5e7eb">'+fmt(M.m12)+'</span>'+
      '<span style="color:#e5e7eb">'+fmt(M.m21)+'</span><span style="color:#e5e7eb">'+fmt(M.m22)+'</span>'+
      '</span>'
  }

  function updateReadouts(state){
    const {d, pr} = state
    par.querySelector('#stReadout').innerHTML =
      '<b style="color:#38bdf8">e</b> (strain) '+matHTML(d.e, '#38bdf8')+'<br>'+
      '<b style="color:#c084fc">Ω</b> (rotation) '+matHTML(d.Om, '#c084fc')+'<br>'+
      'dilatation &Delta; = tr(e) = <b>'+fmt(d.delta)+'</b><br>'+
      'principal strains: <b>'+fmt(pr.l1)+'</b>, <b>'+fmt(pr.l2)+'</b> at '+Math.round(pr.angle*180/Math.PI)+'&deg;<br>'+
      'cut-plane &theta; = '+Math.round(theta*180/Math.PI)+'&deg;'
    par.querySelector('#stLegend').innerHTML =
      '<span style="color:#facc15">A</span>/<span style="color:#f472b6">B</span> handles, or the <b>J</b> boxes above, set the same matrix &mdash; both stay in sync<br>'+
      '<span style="color:#38bdf8">strain</span> &middot; <span style="color:#c084fc">rotation</span> &middot; '+
      '<span style="color:#facc15">dilatation</span> &middot; <span style="color:#4ade80">deviatoric</span><br>'+
      'dashed = reference cut &middot; <span style="color:#22d3ee">cyan</span> = same fiber after deformation'
  }

  // The four J entry boxes mirror j11..j22 both ways: dragging A/B (or picking a preset) updates
  // the boxes here, and typing in a box updates j11..j22 (registered below, next to the drag
  // handlers) -- skip overwriting whichever box currently has focus so a live edit's cursor/
  // selection isn't clobbered out from under the person typing.
  const jBoxes = {j11:par.querySelector('#st-j11'), j12:par.querySelector('#st-j12'),
                  j21:par.querySelector('#st-j21'), j22:par.querySelector('#st-j22')}
  function syncJBoxes(){
    const vals = {j11,j12,j21,j22}
    for(const k in jBoxes){
      if(document.activeElement !== jBoxes[k]) jBoxes[k].value = Math.round(vals[k]*1000)/1000
    }
  }

  function drawAll(){
    const d = decompose()
    const state = {d, pr:principal(d.e)}
    drawMain(); drawStrip(state); drawMohr(state); updateReadouts(state); syncJBoxes()
  }

  function setCustom(){
    if(currentPreset !== 'custom'){ currentPreset = 'custom'; par.querySelector('#stPreset').value = 'custom' }
  }

  drawAll()

  mainCvs.addEventListener('mousedown', e=>{
    const mx=e.offsetX, my=e.offsetY
    const h = handlePos(SEC)
    const near = p => Math.hypot(p[0]-mx, p[1]-my) < 14
    if(near(h.A)) dragMode = 'A'
    else if(near(h.B)) dragMode = 'B'
  })
  mainCvs.addEventListener('mousemove', e=>{
    if(!dragMode) return
    const inv = invxform(SEC)
    const [x,y] = inv(e.offsetX, e.offsetY)
    if(dragMode === 'A'){ j11 = x-1; j21 = y } else { j12 = x; j22 = y-1 }
    setCustom()
    drawAll()
  })
  window.addEventListener('mouseup', ()=>{ dragMode = null; mohrDrag = false })

  jBoxes.j11.addEventListener('input', e=>{ const v=parseFloat(e.target.value); if(Number.isFinite(v)){ j11=v; setCustom(); drawAll() } })
  jBoxes.j12.addEventListener('input', e=>{ const v=parseFloat(e.target.value); if(Number.isFinite(v)){ j12=v; setCustom(); drawAll() } })
  jBoxes.j21.addEventListener('input', e=>{ const v=parseFloat(e.target.value); if(Number.isFinite(v)){ j21=v; setCustom(); drawAll() } })
  jBoxes.j22.addEventListener('input', e=>{ const v=parseFloat(e.target.value); if(Number.isFinite(v)){ j22=v; setCustom(); drawAll() } })

  let mohrDrag = false
  mohrCvs.addEventListener('mousedown', ()=>{ mohrDrag = true })
  mohrCvs.addEventListener('mousemove', e=>{
    if(!mohrDrag) return
    const d = decompose()
    const {A,B,c} = mohrGeom(d.e)
    const [dataX,dataY] = mohrInvXform(e.offsetX, e.offsetY)
    const arg0 = Math.atan2(B, A)
    const phi = Math.atan2(dataY, dataX - c)
    theta = (arg0 - phi)/2
    drawAll()
  })

  par.querySelector('#stPreset').addEventListener('change', e=>{
    const key = e.target.value
    currentPreset = key
    if(key !== 'custom'){
      const pr = PRESETS[key]
      j11=pr.j11; j12=pr.j12; j21=pr.j21; j22=pr.j22
    }
    drawAll()
  })
</script>
""")
    end

    const _st_ready = true
end

# ╔═╡ 8b707bb1-f496-4d94-b81d-8a5fbb69a7f7
begin
    _st_ready
    @bind stw StrainTensorInput()
end

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
# ╠═3f1bc2d1-2327-48cc-984a-df09c936da87
# ╟─fd14256d-8c69-4d8a-91b5-924a32479866
# ╟─4d860dc6-4944-44da-81f2-14b4aa018694
# ╟─11966694-3009-4016-be5b-2a4eb1be6be2
# ╟─174e1f52-bf7a-4201-b579-c784115d15f1
# ╟─b2cf2af1-564e-4a3f-a846-a50042883a7a
# ╟─bbb43caf-32a4-4f70-a263-246953076d84
# ╟─a660a926-0835-4ed2-a2a6-615e3dba3088
# ╟─8b707bb1-f496-4d94-b81d-8a5fbb69a7f7
# ╟─d8f1ceea-1968-4e59-93ae-a21f035b2a09
# ╟─6400c97b-9f07-4f52-9002-e8b5f11aa182
# ╠═a759772e-a5b4-4f35-9150-c17627b58c4a
# ╠═9eeeac07-c2f6-4759-a328-4dee9d576723
# ╟─f4b45b1b-b7f8-4f2b-a580-98def9cc8fbc
# ╟─31eb56b0-25dd-4335-bc17-f89431d79b22
# ╟─a4ec2280-8bac-4e65-a849-909441cee0b6
# ╠═c098aa8d-6106-4bf4-b6c4-b2c75e9225a2
# ╟─6f33c7ed-e8a8-440d-a771-84789ee1f397
# ╠═13a3429e-12f6-11ed-326f-c154f5debceb
# ╟─ca079784-f5ac-4e9b-9d10-8501564c47f0
# ╠═962355bc-2e2e-4bba-b293-aefdca8f3627
# ╟─3f45f0c5-e672-4780-a2c6-006e49fecfb8
# ╠═3f7a0102-597e-49cc-978b-39ad7fcb3015
# ╟─6b829564-d656-4c81-885f-212d1b8d6d91
# ╠═0168e21c-7af1-491d-be1b-7b3f458f9ad3
# ╟─906d8e1c-67c5-46ca-ad8d-0e087561690c
# ╠═e2c19bc6-e3df-4c93-aded-b9fda77247ec
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
