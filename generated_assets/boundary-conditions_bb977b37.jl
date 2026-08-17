### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> layout = "layout.jlhtml"
#> title = "Boundary Conditions"
#> tags = ["planewaves"]
#> description = "Boundary conditions at an interface between two different media"

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

# ╔═╡ a08acdbe-9e7f-4404-8a29-aab612984839
using PlutoUI

# ╔═╡ 5a5bbfc8-6ff7-11ef-39cb-b12385eb48ed
TableOfContents()

# ╔═╡ 486cb239-e68d-491c-a158-5320344caf25
md"""
# Boundary Conditions

A seismic wave doesn't just keep propagating obliviously when it hits an interface between
two different rock types, or the ocean floor, or the free surface of the Earth itself —
part of it reflects, part transmits, and the amplitudes of each are entirely dictated by
what has to stay continuous (and what's allowed to jump) right at that interface. Those
rules are the **boundary conditions**, and they're not an afterthought bolted onto the wave
equation: they're the reason reflection/transmission coefficients, head waves, and even
surface waves exist at all. Get them wrong and every downstream calculation — an AVO
analysis, a receiver function, a dispersion curve — is built on the wrong physics.

Every interface in this notebook is described by two things at a single point on it: its
unit normal vector `` \mathbf{n} ``, and the stress tensors `` \sigma^{(1)} `` (just above
the interface) and `` \sigma^{(2)} `` (just below), producing displacement fields
`` \mathbf{u}^{(1)} `` and `` \mathbf{u}^{(2)} ``. Two families of condition apply at that
point:

- **Kinematic conditions** — what the *displacement* field is allowed to do across the
  interface (stay glued together, or slip).
- **Dynamic conditions** — what the *traction* `` \sigma\mathbf{n} `` (the force per unit
  area transmitted across the interface) is allowed to do — this is just Newton's third
  law applied to an infinitesimal slab straddling the boundary.

Three physically distinct interfaces recur throughout seismology, and each imposes a
different pair of conditions:
- **Solid–Solid** (welded) — e.g. two rock layers bonded together.
- **Solid–Fluid** — e.g. the seafloor, or a magma chamber boundary.
- **Solid–Vacuum** (free surface) — e.g. the Earth's actual surface.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 0b2dc3ed-0491-485d-9571-29d5759f81d7
md"""
## Explore It

The widget below puts a single point on an interface under your control. Drag the two
dots to set the displacement `` \mathbf{u}^{(1)}, \mathbf{u}^{(2)} `` on each side; scrub
any stress-tensor entry left/right to change it; switch the interface type above the
canvas and watch which quantities are even allowed to differ. The two colored arrows are
the traction vectors `` \mathbf{t}^{(1)}=\sigma^{(1)}\mathbf{n} `` (cyan, medium 1) and
`` \mathbf{t}^{(2)}=\sigma^{(2)}\mathbf{n} `` (amber, medium 2) — when they land on top of
each other, traction is continuous. The zig-zag "spring" between the two displacement dots
turns red the instant the current interface type's kinematic condition is violated, and
green when it's satisfied — a direct, literal picture of the springs-not-broken language
used to describe a welded interface below.

The interface's own local frame is fixed for clarity: `` \mathbf{n} `` points straight up
out of medium 2 into medium 1, and the tangential direction is the one other direction in
this 2D cross-section. The **Detailed Look** section below restates every condition the
widget enforces in the fully general 3D form, for an arbitrary `` \mathbf{n} ``.
"""

# ╔═╡ adbec4f1-f939-40da-a0c1-054ada342c21
md"""
## Detailed Look at Each Interface Type

The widget above works in a fixed local (tangential, normal) frame. The equations below
restate the same three conditions in the fully general 3D form, for an arbitrary unit
normal `` \mathbf{n} `` — exactly what a real, tilted interface needs.
"""

# ╔═╡ 9cbee290-2093-4741-a217-fa715d57ea7b
md"""
### Solid–Solid (Welded Interface)

Consider a welded interface between two solids, where none of the springs connecting
particles of one medium to the other are "broken" — nothing can separate, open a gap, or
slide.

**Kinematic condition.** A welded interface means *all* components of the displacement
field are continuous across the boundary:

```math
\mathbf{u}^{(1)} = \mathbf{u}^{(2)}
```

**Dynamic condition.** This follows from considering a Gaussian pillbox of negligible
volume straddling the interface: it can have zero net force (and hence zero acceleration)
only if the traction on the interface is continuous across it:

```math
\sigma^{(1)}\mathbf{n} = \sigma^{(2)}\mathbf{n}
```

In the widget, switch to **Welded** and drag either dot away from the other — the spring
immediately goes red, and the readout panel flags whichever component (normal or
tangential) you broke. Click **Snap medium 2 to satisfy** to see both conditions restored
at once.
"""

# ╔═╡ cac42b11-d044-4e9e-be26-806654c1172c
md"""
### Solid–Fluid

Unlike the welded interface, some of the springs are broken: particles can slide
tangentially to the interface, because a fluid cannot apply a shear force to the solid
particles next to it. Take medium 1 to be the fluid.

A fluid's stress state is necessarily **isotropic**, `` \sigma^{(1)} = -p\,I `` for some
pressure `` p `` — there is no preferred direction, which is exactly *why* it cannot
support shear. (This is the same rule the widget enforces: switching to **Solid–Fluid**
locks `` \sigma^{(1)} `` to a single pressure slider instead of a free tensor.)

**Kinematic condition.** The normal component of displacement is continuous — the fluid
and solid can't separate or interpenetrate — but the tangential components are free to
differ (sliding):

```math
\mathbf{u}^{(1)}\cdot\mathbf{n} = \mathbf{u}^{(2)}\cdot\mathbf{n}
```

**Dynamic condition.** The normal component of traction (the pressure balance) is
continuous, and the tangential (shear) component of traction vanishes on *both* sides —
automatically on the fluid side because `` \sigma^{(1)} `` has no shear at all, and as a
consequence on the solid side too, since traction must match:

```math
\mathbf{n}\cdot\bigl(\sigma^{(1)}\mathbf{n}\bigr) = \mathbf{n}\cdot\bigl(\sigma^{(2)}\mathbf{n}\bigr)
```
```math
\sigma^{(2)}\mathbf{n} - \Bigl[\mathbf{n}\cdot\bigl(\sigma^{(2)}\mathbf{n}\bigr)\Bigr]\mathbf{n} = \mathbf{0}
```

In the widget, drag the two dots sideways past each other in **Solid–Fluid** mode — the
spring stays green and is labeled "slip free," since that tangential offset is exactly
what this condition permits. Now try to give medium 2 (the solid) a shear traction by
dragging its own stress entries — the tangential-traction readout turns red, even though
nothing about the displacement changed.
"""

# ╔═╡ 9b776562-f534-4a27-8433-cbeb8f7e4133
md"""
### Solid–Vacuum (Free Surface)

A vacuum has no particles at all to exert a force on the solid next to it, so the solid's
displacement is completely unconstrained — there is no kinematic condition to write down.
Taking medium 2 to be the solid, the only requirement is that the *entire* traction vector
on the interface vanishes:

```math
\sigma^{(2)}\mathbf{n} = \mathbf{0}
```

This single vector equation is why a free surface reflects seismic energy so efficiently
(nothing is there to carry any of it onward) and why surface waves — which need exactly
this zero-traction condition at `` z=0 `` — exist at all.

In the widget, switch to **Free surface**: medium 1 disappears entirely (there's nothing
there), and the only thing being checked is whether medium 2's own traction arrow shrinks
to a point at the interface.
"""

# ╔═╡ bb6dddb4-19b8-420e-805c-4a143fa8ebed
md"""
## Appendix
"""

# ╔═╡ c514fac3-7874-4a6c-8848-56c66ca80260
md"""
### Traction and Continuity Physics

Everything the widget draws or flags red/green is computed here, in Julia, from the raw
values the widget reports — the widget itself only handles dragging, scrubbing, and
drawing (see the `pluto-widget-style` skill's "physics lives in Julia" rule).
"""

# ╔═╡ 491b8789-e1b3-441c-8cac-f9ec24d2e653
"""
	traction_2d(σtt, σtn, σnn, n̂)

Cauchy traction ``t = \\sigma \\hat n`` for a 2D symmetric stress state given in the local
(tangential, normal) basis, evaluated on a surface with unit normal `n̂` (also expressed in
that basis). Returns `(t_t, t_n)`, the traction resolved back onto the same (tangential,
normal) axes.

!!! note "Why this stays general"
	The widget always calls this with `n̂ = (0,1)` (its fixed local frame), which makes
	`t_t = σtn` and `t_n = σnn` — but the function itself performs the real dot product for
	an arbitrary `n̂`, so it's the same formula the Detailed Look section's 3D
	`` \\sigma\\mathbf{n} `` equations describe, not a shortcut specific to one frame.
"""
function traction_2d(σtt::Float64, σtn::Float64, σnn::Float64, n̂::NTuple{2,Float64})
	nt, nn = n̂
	tt = σtt * nt + σtn * nn
	tn = σtn * nt + σnn * nn
	return (tt, tn)
end

# ╔═╡ c20ee476-b26c-470f-819d-6d413bf8ab0c
"""
	fluid_stress(p)

Isotropic stress state ``\\sigma=-p\\,I`` for a fluid at pressure `p` (tension-positive
convention throughout this notebook, so a positive `p` gives a compressive normal stress).
Returns `(σtt, σtn, σnn)` in the local (tangential, normal) basis — the same in *any*
basis, since an isotropic tensor has no preferred direction, which is exactly why a fluid
cannot support shear traction.
"""
fluid_stress(p::Float64) = (-p, 0.0, -p)

# ╔═╡ ce088892-7d2f-482e-88db-578548954f1d
"""
	interface_check(itype, t1t, t1n, t2t, t2n, u1t, u1n, u2t, u2n; tol=1e-3)

Evaluate the kinematic (displacement) and dynamic (traction) boundary conditions at a
welded (`:welded`), solid–fluid (`:fluid`), or free-surface (`:vacuum`) interface, given
the already-resolved tractions `` \\mathbf{t}^{(1)}=(t1t,t1n) ``,
`` \\mathbf{t}^{(2)}=(t2t,t2n) `` and displacements `` \\mathbf{u}^{(1)}=(u1t,u1n) ``,
`` \\mathbf{u}^{(2)}=(u2t,u2n) ``.

Returns a named tuple of `Bool`s: `kin_normal_ok`, `kin_tangential_ok`,
`kin_tangential_free` (tangential slip is unconstrained by this interface type, e.g. a
fluid or free surface), `dyn_normal_ok`, `dyn_tangential_ok`, plus the raw residuals
`dun`, `dut`, `dtn`, `dtt` those checks are built from.

!!! note "Solid–fluid asymmetry"
	For `:fluid`, the tangential *dynamic* check tests `t2t` (the solid side) against
	zero directly, rather than comparing `t1t` to `t2t` — a fluid's own shear traction is
	already zero by construction (see [`fluid_stress`](@ref)), so the real physical
	content of "shear traction vanishes at a fluid interface" is a constraint on the
	*solid* side alone.
"""
function interface_check(itype::Symbol, t1t::Float64, t1n::Float64, t2t::Float64, t2n::Float64,
	u1t::Float64, u1n::Float64, u2t::Float64, u2n::Float64; tol::Float64=1e-3)
	dun = u1n - u2n
	dut = u1t - u2t
	dtn = t1n - t2n
	dtt = t1t - t2t
	if itype == :welded
		return (kin_normal_ok=abs(dun) < tol, kin_tangential_ok=abs(dut) < tol, kin_tangential_free=false,
			dyn_normal_ok=abs(dtn) < tol, dyn_tangential_ok=abs(dtt) < tol, dun=dun, dut=dut, dtn=dtn, dtt=dtt)
	elseif itype == :fluid
		return (kin_normal_ok=abs(dun) < tol, kin_tangential_ok=true, kin_tangential_free=true,
			dyn_normal_ok=abs(dtn) < tol, dyn_tangential_ok=abs(t2t) < tol, dun=dun, dut=dut, dtn=dtn, dtt=dtt)
	else # :vacuum -- medium 1 is vacuum, only medium 2's own traction is constrained
		return (kin_normal_ok=true, kin_tangential_ok=true, kin_tangential_free=true,
			dyn_normal_ok=abs(t2n) < tol, dyn_tangential_ok=abs(t2t) < tol, dun=dun, dut=dut, dtn=dtn, dtt=dtt)
	end
end

# ╔═╡ b2c3cc61-69d7-4f3d-a394-4010437960b8
md"""
#### Verifying the Continuity Checks
"""

# ╔═╡ dd0e1670-4704-41ac-854f-5e00479aa09d
let
	n̂ = (0.0, 1.0)

	# traction_2d reduces to (σtn, σnn) for the widget's own n̂=(0,1)...
	@assert traction_2d(3.0, 1.5, -2.0, n̂) == (1.5, -2.0)
	# ...but is a genuine dot product for a general n̂, not a hardcoded component pick
	n45 = (sqrt(2) / 2, sqrt(2) / 2)
	tt45, tn45 = traction_2d(1.0, 0.0, -1.0, n45)
	@assert isapprox(tt45, sqrt(2) / 2; atol=1e-12) && isapprox(tn45, -sqrt(2) / 2; atol=1e-12)

	@assert fluid_stress(4.0) == (-4.0, 0.0, -4.0)

	# welded: matching state passes both conditions; a tangential-slip violation is caught
	t1 = traction_2d(2.0, 0.5, -1.0, n̂)
	t2 = traction_2d(99.0, 0.5, -1.0, n̂)  # σtt irrelevant to traction, deliberately different
	r_ok = interface_check(:welded, t1..., t2..., 0.3, -0.1, 0.3, -0.1)
	r_bad = interface_check(:welded, t1..., t2..., 0.3, -0.1, 0.9, -0.1)
	@assert r_ok.kin_normal_ok && r_ok.kin_tangential_ok && r_ok.dyn_normal_ok && r_ok.dyn_tangential_ok
	@assert r_bad.dyn_normal_ok && r_bad.dyn_tangential_ok && !r_bad.kin_tangential_ok

	# solid-fluid: tangential slip never flags; nonzero solid shear traction does
	p = 4.0
	t1f = traction_2d(fluid_stress(p)..., n̂)
	t2f = traction_2d(1.0, 0.0, -p, n̂)
	t2f_bad = traction_2d(1.0, 0.6, -p, n̂)
	r_fluid = interface_check(:fluid, t1f..., t2f..., 0.0, 0.2, 0.7, 0.2)
	r_fluid_bad = interface_check(:fluid, t1f..., t2f_bad..., 0.0, 0.2, 0.2, 0.2)
	@assert r_fluid.kin_normal_ok && r_fluid.kin_tangential_free && r_fluid.dyn_normal_ok && r_fluid.dyn_tangential_ok
	@assert !r_fluid_bad.dyn_tangential_ok && r_fluid_bad.dyn_normal_ok

	# free surface: only the solid side's own traction is checked, against zero
	t2v_ok = traction_2d(1.0, 0.0, 0.0, n̂)
	t2v_bad = traction_2d(1.0, 0.4, -2.0, n̂)
	r_vac_ok = interface_check(:vacuum, 0.0, 0.0, t2v_ok..., 0.0, 0.0, 5.0, -3.0)
	r_vac_bad = interface_check(:vacuum, 0.0, 0.0, t2v_bad..., 0.0, 0.0, 5.0, -3.0)
	@assert r_vac_ok.dyn_normal_ok && r_vac_ok.dyn_tangential_ok
	@assert !r_vac_bad.dyn_normal_ok && !r_vac_bad.dyn_tangential_ok

	md"All nine continuity-check assertions above (welded / solid–fluid / free-surface, satisfied and violated) pass. ✅"
end

# ╔═╡ 18152216-ae92-41f5-b5eb-37c3421b7840
md"""
### From the Widget's Raw State to the Physics

`bc` is the widget's raw bound dict (interface type, stress-tensor entries, displacement
vectors). This section turns that into the traction vectors and pass/fail verdicts pushed
back to the widget as a `CustomEvent`.
"""

# ╔═╡ b06f333a-ac4b-4bc9-8deb-271578cd4d6d
# the widget's own fixed local frame -- n points out of medium 2, into medium 1
_bc_n̂ = (0.0, 1.0)

# ╔═╡ 605c8c08-ce05-4c5c-adbc-9fb19398026e
md"""
`BcPush` carries no physics of its own — just the already-computed tractions and pass/fail
verdicts, handed to the *already-rendered* `BoundaryConditionInput` widget via a
`CustomEvent`, exactly like `MtPush`/`PfrPush` elsewhere in this repo.
"""

# ╔═╡ c56c134d-8a38-41e2-bace-f6016f3d6038
begin
	struct BcPush
		t1t::Float64
		t1n::Float64
		t2t::Float64
		t2n::Float64
		u1t::Float64
		u1n::Float64
		u2t::Float64
		u2n::Float64
		kinNormalOk::Bool
		kinTangentialOk::Bool
		kinTangentialFree::Bool
		dynNormalOk::Bool
		dynTangentialOk::Bool
	end
	function Base.show(io::IO, ::MIME"text/html", p::BcPush)
		write(io, """
		<script>
		{
		const w = document.getElementById('bcwidget');
		if(w){
		  w.dispatchEvent(new CustomEvent('bc-update', { detail: {
		    t1t: $(p.t1t), t1n: $(p.t1n), t2t: $(p.t2t), t2n: $(p.t2n),
		    u1t: $(p.u1t), u1n: $(p.u1n), u2t: $(p.u2t), u2n: $(p.u2n),
		    kinNormalOk: $(p.kinNormalOk), kinTangentialOk: $(p.kinTangentialOk), kinTangentialFree: $(p.kinTangentialFree),
		    dynNormalOk: $(p.dynNormalOk), dynTangentialOk: $(p.dynTangentialOk)
		  }}));
		}
		}
		</script>
		""")
	end
end

# ╔═╡ 3e1bb6d0-9c2b-4af2-9e18-c187ae808ccf
md"""
### The Interactive Widget

`BoundaryConditionInput` only ever handles UI/display: dragging the two displacement dots,
scrubbing stress-tensor entries, toggling the interface type, and drawing whatever the
physics cells above last pushed to it via the `bc-update` event. It never recomputes a
traction or a pass/fail verdict itself.
"""

# ╔═╡ 236b5432-5064-4207-a68f-e2c651fee95c
begin
	struct BoundaryConditionInput
		itype::String
		s1tt::Float64
		s1tn::Float64
		s1nn::Float64
		p1::Float64
		s2tt::Float64
		s2tn::Float64
		s2nn::Float64
		u1t::Float64
		u1n::Float64
		u2t::Float64
		u2n::Float64
	end

	BoundaryConditionInput(; itype="welded", s1tt=0.3, s1tn=0.4, s1nn=1.0, p1=1.0,
		s2tt=-0.3, s2tn=0.4, s2nn=1.0, u1t=0.3, u1n=-0.2, u2t=0.3, u2n=-0.2) =
		BoundaryConditionInput(itype, Float64(s1tt), Float64(s1tn), Float64(s1nn), Float64(p1),
			Float64(s2tt), Float64(s2tn), Float64(s2nn),
			Float64(u1t), Float64(u1n), Float64(u2t), Float64(u2n))

	Base.get(w::BoundaryConditionInput) = Dict{String,Any}(
		"itype" => w.itype, "s1tt" => w.s1tt, "s1tn" => w.s1tn, "s1nn" => w.s1nn, "p1" => w.p1,
		"s2tt" => w.s2tt, "s2tn" => w.s2tn, "s2nn" => w.s2nn,
		"u1t" => w.u1t, "u1n" => w.u1n, "u2t" => w.u2t, "u2n" => w.u2n,
	)

	function Base.show(io::IO, ::MIME"text/html", w::BoundaryConditionInput)
		write(io, """
<div id="bcwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#bcwidget){width:min(78vw,1250px)!important;margin-left:calc((100% - min(78vw,1250px))/2)!important}
    #bcwidget .bc-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #bcwidget .bc-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #bcwidget .bc-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #bcwidget .bc-itype{display:flex;gap:8px;margin-bottom:12px}
    #bcwidget .bc-itype button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:7px 14px;font-size:14px;cursor:pointer}
    #bcwidget .bc-itype button.active{background:#3b5c85;border-color:#38bdf8}
    #bcwidget .bc-workspace{display:flex;gap:16px;align-items:flex-start;width:100%;justify-content:center}
    #bcwidget .bc-canvas-col{display:flex;flex-direction:column;align-items:center;flex:0 0 auto}
    #bcwidget canvas{background:#000;border:1px solid #374151;border-radius:6px;display:block}
    #bcwidget .bc-readout{width:100%;box-sizing:border-box;margin-top:10px;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:8px 12px;font-size:13px}
    #bcwidget .bc-readout-row{display:grid;grid-template-columns:1fr auto;gap:8px;align-items:center;padding:4px 0;border-bottom:1px solid #1f2937}
    #bcwidget .bc-readout-row:last-child{border-bottom:none}
    #bcwidget .bc-readout-label{color:#d1d5db}
    #bcwidget .bc-readout-val{font-variant-numeric:tabular-nums;color:#9ca3af;font-size:12px}
    #bcwidget .bc-badge{border-radius:10px;padding:2px 9px;font-size:12px;font-weight:700;white-space:nowrap}
    #bcwidget .bc-badge-ok{background:#0f2e1a;color:#4ade80;border:1px solid #16a34a}
    #bcwidget .bc-badge-bad{background:#3a0f12;color:#f87171;border:1px solid #dc2626}
    #bcwidget .bc-badge-free{background:#0f1e2e;color:#7dd3fc;border:1px solid #0369a1}
    #bcwidget .bc-controls{display:flex;flex-direction:column;gap:10px;flex:0 0 300px;width:300px}
    #bcwidget .bc-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
    #bcwidget .bc-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:15px}
    #bcwidget .bc-control-hint{font-size:12px;color:#9ca3af;margin-top:6px;line-height:1.5}
    #bcwidget .bc-mat2{display:inline-grid;grid-template-columns:repeat(2,auto);gap:5px 6px;padding:6px 10px;position:relative}
    #bcwidget .bc-mat2::before, #bcwidget .bc-mat2::after{content:'';position:absolute;top:3px;bottom:3px;width:6px;border:2px solid #9ca3af}
    #bcwidget .bc-mat2::before{left:0;border-right:none}
    #bcwidget .bc-mat2::after{right:0;border-left:none}
    #bcwidget .bc-mat2 input{width:52px;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:4px 3px;font-size:13px;text-align:center;font-variant-numeric:tabular-nums;cursor:ew-resize;-moz-appearance:textfield}
    #bcwidget .bc-mat2 input::-webkit-inner-spin-button, #bcwidget .bc-mat2 input::-webkit-outer-spin-button{-webkit-appearance:none;margin:0}
    #bcwidget .bc-mat2 input:focus{outline:2px solid #38bdf8;border-color:#38bdf8}
    #bcwidget .bc-mat2 input:disabled{opacity:0.55;cursor:default}
    #bcwidget .bc-pgroup{display:flex;align-items:center;gap:8px;margin-bottom:6px}
    #bcwidget .bc-pgroup input[type=range]{flex:1 1 auto;min-width:0}
    #bcwidget .bc-value{color:#d1d5db;font-variant-numeric:tabular-nums;font-size:13px;min-width:40px;text-align:right}
    #bcwidget .bc-snap{width:100%;border-radius:4px;border:1px solid #9ca3af;background:#3b5c85;color:#f3f4f6;padding:8px 12px;font-size:14px;cursor:pointer;font-weight:700}
    #bcwidget .bc-snap:hover{background:#4b7099}
    #bcwidget .bc-dim{opacity:0.45}
  </style>

  <div class="bc-title">
    <div class="bc-title-desc">Kinematic and dynamic conditions at an interface between two media</div>
    <div class="bc-title-hint">drag the two dots to set u&sup1;, u&sup2; &middot; drag a stress value to scrub it &middot; switch interface type above</div>
  </div>

  <div class="bc-itype">
    <button type="button" id="bc-it-welded" class="active">Welded (solid&ndash;solid)</button>
    <button type="button" id="bc-it-fluid">Solid&ndash;fluid</button>
    <button type="button" id="bc-it-vacuum">Free surface (solid&ndash;vacuum)</button>
  </div>

  <div class="bc-workspace">
    <div class="bc-canvas-col">
      <canvas id="bc-canvas" width="520" height="380"></canvas>
      <div class="bc-readout" id="bc-readout"></div>
    </div>
    <div class="bc-controls">
      <div class="bc-control-group" id="bc-group-1">
        <div class="bc-control-title" id="bc-g1-title">Medium 1</div>
        <div class="bc-mat2" id="bc-s1-mat">
          <input type="text" id="bc-s1-tt" value="0.30">
          <input type="text" id="bc-s1-tn-a" value="0.40">
          <input type="text" id="bc-s1-tn-b" value="0.40">
          <input type="text" id="bc-s1-nn" value="1.00">
        </div>
        <div class="bc-pgroup" id="bc-p1-row" style="display:none">
          <span>p</span><input type="range" id="bc-p1" min="0" max="3" step="0.05" value="1.0"><span class="bc-value" id="bc-p1-v">1.00</span>
        </div>
        <div class="bc-control-hint" id="bc-g1-hint">σ¹ is a free symmetric tensor in this welded case.</div>
      </div>
      <div class="bc-control-group">
        <div class="bc-control-title">Medium 2 (solid)</div>
        <div class="bc-mat2" id="bc-s2-mat">
          <input type="text" id="bc-s2-tt" value="-0.30">
          <input type="text" id="bc-s2-tn-a" value="0.40">
          <input type="text" id="bc-s2-tn-b" value="0.40">
          <input type="text" id="bc-s2-nn" value="1.00">
        </div>
      </div>
      <div class="bc-control-group">
        <div class="bc-control-title">Actions</div>
        <button type="button" class="bc-snap" id="bc-snap">Snap medium 2 to satisfy</button>
        <div class="bc-control-hint">Copies whatever medium&nbsp;2 quantities the current interface type requires to match medium&nbsp;1 (or zero, for a free surface).</div>
      </div>
    </div>
  </div>
</div>

<script>
  {
    const par = currentScript.previousElementSibling;

    let itype = '$(w.itype)';
    let s1tt = $(w.s1tt), s1tn = $(w.s1tn), s1nn = $(w.s1nn), p1 = $(w.p1);
    let s2tt = $(w.s2tt), s2tn = $(w.s2tn), s2nn = $(w.s2nn);
    let u1t = $(w.u1t), u1n = $(w.u1n), u2t = $(w.u2t), u2n = $(w.u2n);

    // physics results, filled in once the first 'bc-update' arrives from Julia
    let phys = null;

    const cv = par.querySelector('#bc-canvas');
    const ctx = cv.getContext('2d');
    const SEC_W = 520, SEC_H = 380;
    const PX = 260, PY = 190;           // the interface point, in canvas pixels
    const IF_X0 = 30, IF_X1 = 490;      // interface line extent
    const STRESS_SCALE = 46, DISP_SCALE = 90;
    const MAXLEN = 150;
    const COL1 = '#38bdf8', COL2 = '#f59e0b';

    function clampLen(dx, dy){
      const len = Math.hypot(dx, dy);
      if(len <= MAXLEN || len === 0) return [dx, dy];
      const s = MAXLEN / len;
      return [dx*s, dy*s];
    }

    function drawArrow(x0, y0, dxRaw, dyRaw, color, label){
      const [dx, dy] = clampLen(dxRaw, dyRaw);
      const x1 = x0 + dx, y1 = y0 - dy;   // +n is up on screen
      ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = 2.4;
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
      const ang = Math.atan2(y1-y0, x1-x0), hs = 8;
      ctx.beginPath();
      ctx.moveTo(x1, y1);
      ctx.lineTo(x1 - hs*Math.cos(ang - 0.4), y1 - hs*Math.sin(ang - 0.4));
      ctx.lineTo(x1 - hs*Math.cos(ang + 0.4), y1 - hs*Math.sin(ang + 0.4));
      ctx.closePath(); ctx.fill();
      if(label){
        ctx.font = '13px sans-serif'; ctx.textAlign = 'left';
        ctx.fillText(label, x1 + 6, y1);
      }
      return [x1, y1];
    }

    function particlePos(ut, un){
      const [dx, dy] = clampLen(ut*DISP_SCALE, un*DISP_SCALE);
      return [PX + dx, PY - dy];
    }

    function drawSpring(x0, y0, x1, y1, color, dashed){
      const n = 7, dx = (x1-x0)/n, dy = (y1-y0)/n;
      const nx = -(y1-y0), ny = (x1-x0);
      const nlen = Math.hypot(nx, ny) || 1;
      const amp = dashed ? 0 : 4.5;
      ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.setLineDash(dashed ? [4,3] : []);
      ctx.beginPath(); ctx.moveTo(x0, y0);
      for(let i=1; i<n; i++){
        const px = x0 + dx*i, py = y0 + dy*i;
        const side = (i % 2 === 0) ? 1 : -1;
        ctx.lineTo(px + (nx/nlen)*amp*side, py + (ny/nlen)*amp*side);
      }
      ctx.lineTo(x1, y1); ctx.stroke();
      ctx.setLineDash([]);
    }

    function draw(){
      ctx.clearRect(0, 0, SEC_W, SEC_H);

      // media regions
      ctx.fillStyle = 'rgba(56,189,248,0.06)'; ctx.fillRect(0, 0, SEC_W, PY);
      ctx.fillStyle = 'rgba(245,158,11,0.06)'; ctx.fillRect(0, PY, SEC_W, SEC_H-PY);
      ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(IF_X0, PY); ctx.lineTo(IF_X1, PY); ctx.stroke();

      ctx.font = '13px sans-serif'; ctx.fillStyle = '#7dd3fc'; ctx.textAlign = 'left';
      const label1 = itype === 'fluid' ? 'Medium 1 (fluid)' : (itype === 'vacuum' ? 'Medium 1 (vacuum)' : 'Medium 1');
      ctx.fillText(label1, 12, 20);
      ctx.fillStyle = '#fbbf24';
      ctx.fillText('Medium 2 (solid)', 12, SEC_H - 12);

      // axis compass
      const ax = 50, ay = 60;
      ctx.strokeStyle = '#6b7280'; ctx.fillStyle = '#6b7280'; ctx.lineWidth = 1.4;
      drawArrow(ax, ay, 0, 26, '#9ca3af', 'n');
      drawArrow(ax, ay, 26, 0, '#9ca3af', 't');

      const vacuum = (itype === 'vacuum');

      // particles at rest + guide lines
      const [p1x, p1y] = particlePos(u1t, u1n);
      const [p2x, p2y] = particlePos(u2t, u2n);
      ctx.strokeStyle = 'rgba(156,163,175,0.35)'; ctx.setLineDash([3,3]); ctx.lineWidth = 1;
      if(!vacuum){ ctx.beginPath(); ctx.moveTo(PX, PY); ctx.lineTo(p1x, p1y); ctx.stroke(); }
      ctx.beginPath(); ctx.moveTo(PX, PY); ctx.lineTo(p2x, p2y); ctx.stroke();
      ctx.setLineDash([]);

      // spring / roller between the two particles
      if(!vacuum && phys){
        const kinOK = itype === 'welded' ? (phys.kinNormalOk && phys.kinTangentialOk) : phys.kinNormalOk;
        drawSpring(p1x, p1y, p2x, p2y, kinOK ? '#4ade80' : '#f87171', false);
        if(phys.kinTangentialFree){
          const mx = (p1x+p2x)/2, my = (p1y+p2y)/2;
          ctx.fillStyle = '#7dd3fc'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
          ctx.fillText('slip free', mx, my - 8);
        }
      }

      // particle markers
      if(!vacuum){
        ctx.beginPath(); ctx.arc(p1x, p1y, 6, 0, 2*Math.PI); ctx.fillStyle = COL1; ctx.fill();
        ctx.strokeStyle = '#0a0f18'; ctx.lineWidth = 1.4; ctx.stroke();
      }
      ctx.beginPath(); ctx.arc(p2x, p2y, 6, 0, 2*Math.PI); ctx.fillStyle = COL2; ctx.fill();
      ctx.strokeStyle = '#0a0f18'; ctx.lineWidth = 1.4; ctx.stroke();

      // traction arrows, anchored at the interface point P -- t=σn is physics, so these are only
      // ever drawn from Julia-pushed values (phys), never recomputed from the raw σ boxes here
      if(phys){
        if(!vacuum) drawArrow(PX, PY, phys.t1t*STRESS_SCALE, phys.t1n*STRESS_SCALE, COL1, 't¹');
        drawArrow(PX, PY, phys.t2t*STRESS_SCALE, phys.t2n*STRESS_SCALE, COL2, 't²');
      }

      ctx.beginPath(); ctx.arc(PX, PY, 3, 0, 2*Math.PI); ctx.fillStyle = '#e5e7eb'; ctx.fill();

      syncMatBoxes();
      updateReadout();
    }

    function badge(ok, free){
      if(free) return '<span class="bc-badge bc-badge-free">free</span>';
      return ok ? '<span class="bc-badge bc-badge-ok">match</span>' : '<span class="bc-badge bc-badge-bad">mismatch</span>';
    }

    function updateReadout(){
      const rd = par.querySelector('#bc-readout');
      const vacuum = (itype === 'vacuum');
      if(!phys){ rd.innerHTML = '<div class="bc-readout-row"><span class="bc-readout-label">computing…</span></div>'; return; }
      let rows = '';
      if(vacuum){
        rows += '<div class="bc-readout-row"><span class="bc-readout-label">t&sup2;<sub>n</sub> = ' + phys.t2n.toFixed(2) + ' (must be 0)</span>' + badge(phys.dynNormalOk) + '</div>';
        rows += '<div class="bc-readout-row"><span class="bc-readout-label">t&sup2;<sub>t</sub> = ' + phys.t2t.toFixed(2) + ' (must be 0)</span>' + badge(phys.dynTangentialOk) + '</div>';
        rows += '<div class="bc-readout-row"><span class="bc-readout-label">kinematic: u&sup2; unconstrained</span><span class="bc-readout-val">&mdash;</span></div>';
      } else {
        rows += '<div class="bc-readout-row"><span class="bc-readout-label">u&sup1;&middot;n = ' + phys.u1n.toFixed(2) + ',  u&sup2;&middot;n = ' + phys.u2n.toFixed(2) + '</span>' + badge(phys.kinNormalOk) + '</div>';
        rows += '<div class="bc-readout-row"><span class="bc-readout-label">u&sup1;&middot;t = ' + phys.u1t.toFixed(2) + ',  u&sup2;&middot;t = ' + phys.u2t.toFixed(2) + '</span>' + badge(phys.kinTangentialOk, phys.kinTangentialFree) + '</div>';
        rows += '<div class="bc-readout-row"><span class="bc-readout-label">t&sup1;<sub>n</sub> = ' + phys.t1n.toFixed(2) + ',  t&sup2;<sub>n</sub> = ' + phys.t2n.toFixed(2) + '</span>' + badge(phys.dynNormalOk) + '</div>';
        const tanLabel = itype === 'fluid' ? ('t&sup2;<sub>t</sub> = ' + phys.t2t.toFixed(2) + ' (must be 0)') : ('t&sup1;<sub>t</sub> = ' + phys.t1t.toFixed(2) + ',  t&sup2;<sub>t</sub> = ' + phys.t2t.toFixed(2));
        rows += '<div class="bc-readout-row"><span class="bc-readout-label">' + tanLabel + '</span>' + badge(phys.dynTangentialOk) + '</div>';
      }
      rd.innerHTML = rows;
    }

    function syncMatBoxes(){
      par.querySelector('#bc-s1-tt').value = s1tt.toFixed(2);
      par.querySelector('#bc-s1-tn-a').value = s1tn.toFixed(2);
      par.querySelector('#bc-s1-tn-b').value = s1tn.toFixed(2);
      par.querySelector('#bc-s1-nn').value = s1nn.toFixed(2);
      par.querySelector('#bc-s2-tt').value = s2tt.toFixed(2);
      par.querySelector('#bc-s2-tn-a').value = s2tn.toFixed(2);
      par.querySelector('#bc-s2-tn-b').value = s2tn.toFixed(2);
      par.querySelector('#bc-s2-nn').value = s2nn.toFixed(2);
      par.querySelector('#bc-p1').value = p1;
      par.querySelector('#bc-p1-v').textContent = p1.toFixed(2);
    }

    function syncModeUI(){
      par.querySelector('#bc-it-welded').classList.toggle('active', itype === 'welded');
      par.querySelector('#bc-it-fluid').classList.toggle('active', itype === 'fluid');
      par.querySelector('#bc-it-vacuum').classList.toggle('active', itype === 'vacuum');
      const g1 = par.querySelector('#bc-group-1');
      const s1mat = par.querySelector('#bc-s1-mat');
      const p1row = par.querySelector('#bc-p1-row');
      const g1title = par.querySelector('#bc-g1-title');
      const g1hint = par.querySelector('#bc-g1-hint');
      g1.classList.toggle('bc-dim', itype === 'vacuum');
      [...s1mat.querySelectorAll('input')].forEach(el => el.disabled = (itype !== 'welded'));
      if(itype === 'welded'){
        g1title.textContent = 'Medium 1';
        s1mat.style.display = ''; p1row.style.display = 'none';
        g1hint.textContent = 'σ¹ is a free symmetric tensor in this welded case.';
      } else if(itype === 'fluid'){
        g1title.textContent = 'Medium 1 (fluid)';
        s1mat.style.display = ''; p1row.style.display = '';
        g1hint.textContent = 'A fluid is isotropic: σ¹ = −pI, no shear. Drag the pressure slider instead of the tensor.';
      } else {
        g1title.textContent = 'Medium 1 (vacuum)';
        s1mat.style.display = 'none'; p1row.style.display = 'none';
        g1hint.textContent = 'Vacuum has no stress state at all.';
      }
    }

    let commitInFlight = false;
    function commit(){
      commitInFlight = true;
      par.value = { itype, s1tt, s1tn, s1nn, p1, s2tt, s2tn, s2nn, u1t, u1n, u2t, u2n };
      par.dispatchEvent(new CustomEvent('input'));
    }
    function throttledCommit(){ if(!commitInFlight) commit(); }

    par.addEventListener('bc-update', e => {
      phys = e.detail;
      commitInFlight = false;
      draw();
    });

    // ---- interface-type toggle ----
    function setItype(v){
      if(itype === v) return;
      itype = v; syncModeUI(); draw(); commit();
    }
    par.querySelector('#bc-it-welded').addEventListener('click', () => setItype('welded'));
    par.querySelector('#bc-it-fluid').addEventListener('click', () => setItype('fluid'));
    par.querySelector('#bc-it-vacuum').addEventListener('click', () => setItype('vacuum'));

    // ---- pressure slider ----
    par.querySelector('#bc-p1').addEventListener('input', e => {
      p1 = parseFloat(e.target.value);
      par.querySelector('#bc-p1-v').textContent = p1.toFixed(2);
      draw(); throttledCommit();
    });

    // ---- scrub stress-tensor entries (drag left/right; click to type) ----
    function wireMatBoxes(ids, setter){
      const els = ids.map(id => par.querySelector('#'+id));
      function handle(v){ setter(v); draw(); throttledCommit(); }
      els.forEach(el => {
        el.addEventListener('input', e => {
          const v = parseFloat(e.target.value);
          if(Number.isFinite(v)) handle(v);
        });
        let dragging = false, moved = false, startX = 0, startVal = 0;
        el.addEventListener('mousedown', e => {
          if(el.disabled) return;
          startX = e.clientX; startVal = parseFloat(el.value) || 0;
          dragging = true; moved = false;
          e.preventDefault();
        });
        window.addEventListener('mousemove', e => {
          if(!dragging) return;
          const dx = e.clientX - startX;
          if(!moved && Math.abs(dx) > 3){ moved = true; el.blur(); }
          if(!moved) return;
          let v = startVal + dx*0.02;
          v = Math.max(-3, Math.min(3, Math.round(v/0.05)*0.05));
          v = Math.round(v*1000)/1000;
          el.value = v.toFixed(2);
          handle(v);
        });
        window.addEventListener('mouseup', () => {
          if(dragging && !moved){ el.focus(); el.select(); }
          dragging = false;
        });
      });
    }
    wireMatBoxes(['bc-s1-tt'], v => s1tt = v);
    wireMatBoxes(['bc-s1-tn-a', 'bc-s1-tn-b'], v => s1tn = v);
    wireMatBoxes(['bc-s1-nn'], v => s1nn = v);
    wireMatBoxes(['bc-s2-tt'], v => s2tt = v);
    wireMatBoxes(['bc-s2-tn-a', 'bc-s2-tn-b'], v => s2tn = v);
    wireMatBoxes(['bc-s2-nn'], v => s2nn = v);

    // ---- drag the particle markers to set u1 / u2 ----
    let dragU = null;
    cv.addEventListener('mousedown', e => {
      const rect = cv.getBoundingClientRect();
      const mx = e.clientX - rect.left, my = e.clientY - rect.top;
      const vacuum = (itype === 'vacuum');
      const [p1x, p1y] = particlePos(u1t, u1n);
      const [p2x, p2y] = particlePos(u2t, u2n);
      const d1 = vacuum ? Infinity : Math.hypot(mx-p1x, my-p1y);
      const d2 = Math.hypot(mx-p2x, my-p2y);
      const best = Math.min(d1, d2);
      if(best > 16) return;
      // medium 2's dot is drawn last (on top), so an exact tie should grab it, not u1
      dragU = (d1 < d2) ? 'u1' : 'u2';
    });
    window.addEventListener('mousemove', e => {
      if(!dragU) return;
      const rect = cv.getBoundingClientRect();
      const mx = e.clientX - rect.left, my = e.clientY - rect.top;
      const ut = (mx - PX) / DISP_SCALE, un = (PY - my) / DISP_SCALE;
      if(dragU === 'u1'){ u1t = ut; u1n = un; } else { u2t = ut; u2n = un; }
      draw(); throttledCommit();
    });
    window.addEventListener('mouseup', () => {
      if(dragU) commit();
      dragU = null;
    });

    // ---- snap medium 2 to satisfy the current interface type ----
    par.querySelector('#bc-snap').addEventListener('click', () => {
      if(!phys) return;
      if(itype === 'vacuum'){
        s2nn = 0; s2tn = 0;
      } else if(itype === 'fluid'){
        s2nn = phys.t1n; s2tn = 0; u2n = u1n;
      } else {
        s2nn = phys.t1n; s2tn = phys.t1t; u2n = u1n; u2t = u1t;
      }
      draw(); commit();
    });

    syncModeUI(); draw();
  }
  </script>
""")
	end

	const _bc_ready = true
end

# ╔═╡ 159c6430-0478-4b9d-83a4-6ee8ba3edd11
begin
	_bc_ready
	WideCell(@bind bc BoundaryConditionInput(); max_width=1250)
end

# ╔═╡ 7eea3fc3-dbc8-42c2-a850-f3d392844306
bc_safe = bc isa AbstractDict ? bc : Dict{String,Any}(
	"itype" => "welded", "s1tt" => 0.3, "s1tn" => 0.4, "s1nn" => 1.0, "p1" => 1.0,
	"s2tt" => -0.3, "s2tn" => 0.4, "s2nn" => 1.0,
	"u1t" => 0.3, "u1n" => -0.2, "u2t" => 0.3, "u2n" => -0.2)

# ╔═╡ 0dbd33d9-886a-4e77-8b28-01bd10ee3efe
_bc_itype = Symbol(bc_safe["itype"])

# ╔═╡ ba569904-7fdb-4e92-860c-f66d78e963a4
# medium 1 is forced isotropic (a real fluid) whenever the interface type says it's a fluid
_bc_σ1 = _bc_itype == :fluid ? fluid_stress(Float64(bc_safe["p1"])) :
	(Float64(bc_safe["s1tt"]), Float64(bc_safe["s1tn"]), Float64(bc_safe["s1nn"]))

# ╔═╡ 42d92b46-016a-4e3c-823e-4e1639d82131
_bc_t1 = traction_2d(_bc_σ1..., _bc_n̂)

# ╔═╡ f084c2ce-3aa2-49f4-9fe0-806f2e9f6909
_bc_σ2 = (Float64(bc_safe["s2tt"]), Float64(bc_safe["s2tn"]), Float64(bc_safe["s2nn"]))

# ╔═╡ 176d3c95-15ee-40f8-8c36-f096df519daa
_bc_t2 = traction_2d(_bc_σ2..., _bc_n̂)

# ╔═╡ 6f22682b-8967-4c35-968f-9522ea24c6a3
_bc_u1t, _bc_u1n, _bc_u2t, _bc_u2n = Float64(bc_safe["u1t"]), Float64(bc_safe["u1n"]),
Float64(bc_safe["u2t"]), Float64(bc_safe["u2n"])

# ╔═╡ 5c8d43be-e06d-4523-9055-5b1cae7d8c10
_bc_result = interface_check(_bc_itype, _bc_t1..., _bc_t2..., _bc_u1t, _bc_u1n, _bc_u2t, _bc_u2n)

# ╔═╡ d4792547-8d8b-4268-a6ab-90a041a5a032
BcPush(_bc_t1[1], _bc_t1[2], _bc_t2[1], _bc_t2[2], _bc_u1t, _bc_u1n, _bc_u2t, _bc_u2n,
	_bc_result.kin_normal_ok, _bc_result.kin_tangential_ok, _bc_result.kin_tangential_free,
	_bc_result.dyn_normal_ok, _bc_result.dyn_tangential_ok)

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
# ╠═a08acdbe-9e7f-4404-8a29-aab612984839
# ╟─5a5bbfc8-6ff7-11ef-39cb-b12385eb48ed
# ╟─486cb239-e68d-491c-a158-5320344caf25
# ╟─0b2dc3ed-0491-485d-9571-29d5759f81d7
# ╟─159c6430-0478-4b9d-83a4-6ee8ba3edd11
# ╟─adbec4f1-f939-40da-a0c1-054ada342c21
# ╟─9cbee290-2093-4741-a217-fa715d57ea7b
# ╟─cac42b11-d044-4e9e-be26-806654c1172c
# ╟─9b776562-f534-4a27-8433-cbeb8f7e4133
# ╟─bb6dddb4-19b8-420e-805c-4a143fa8ebed
# ╟─c514fac3-7874-4a6c-8848-56c66ca80260
# ╠═491b8789-e1b3-441c-8cac-f9ec24d2e653
# ╠═c20ee476-b26c-470f-819d-6d413bf8ab0c
# ╠═ce088892-7d2f-482e-88db-578548954f1d
# ╟─b2c3cc61-69d7-4f3d-a394-4010437960b8
# ╠═dd0e1670-4704-41ac-854f-5e00479aa09d
# ╟─18152216-ae92-41f5-b5eb-37c3421b7840
# ╠═7eea3fc3-dbc8-42c2-a850-f3d392844306
# ╠═0dbd33d9-886a-4e77-8b28-01bd10ee3efe
# ╠═b06f333a-ac4b-4bc9-8deb-271578cd4d6d
# ╠═ba569904-7fdb-4e92-860c-f66d78e963a4
# ╠═f084c2ce-3aa2-49f4-9fe0-806f2e9f6909
# ╠═42d92b46-016a-4e3c-823e-4e1639d82131
# ╠═176d3c95-15ee-40f8-8c36-f096df519daa
# ╠═6f22682b-8967-4c35-968f-9522ea24c6a3
# ╠═5c8d43be-e06d-4523-9055-5b1cae7d8c10
# ╟─605c8c08-ce05-4c5c-adbc-9fb19398026e
# ╠═c56c134d-8a38-41e2-bace-f6016f3d6038
# ╠═d4792547-8d8b-4268-a6ab-90a041a5a032
# ╟─3e1bb6d0-9c2b-4af2-9e18-c187ae808ccf
# ╠═236b5432-5064-4207-a68f-e2c651fee95c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
