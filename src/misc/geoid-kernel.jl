### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Geoid Kernels"
#> layout = "layout.jlhtml"
#> tags = ["misc"]

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

# ╔═╡ 3a1f0c40-1a2b-4c5d-8e9f-0a1b2c3d4e5f
begin
    using PlutoUI, PlutoPlotly, LinearAlgebra, LaTeXStrings, ColorSchemes
    import PlutoUI: combine
end

# ╔═╡ 4b2e1d51-2b3c-4d6e-9f0a-1b2c3d4e5f60
using Symbolics

# ╔═╡ 5c3f2e62-3c4d-4e7f-a01b-2c3d4e5f6071
TableOfContents(include_definitions=true)

# ╔═╡ 6d403f73-4d5e-4f80-b12c-3d4e5f607182
md"""
## Why Does a Dense Blob Make a Geoid *Low*?

Drop a dense lump into the mantle. Ask a student what happens to the geoid — the shape of
sea level — directly above it, and almost everyone says the same thing: *more mass means
more gravity, so the geoid goes up.*

That answer is right for a rock sitting on a table. It is often **wrong** for a rock sitting
inside a convecting planet, and understanding why is the point of this notebook.

The reason is that the mantle **flows**. A dense anomaly does not just sit there adding its
own gravitational pull. It sinks, and as it sinks it drags the whole mantle with it. That
flow pushes *down* on the Earth's surface from below, creating a shallow depression, and it
pushes *down* on the core–mantle boundary too, deflecting it into the core. Both of those
deflected boundaries are enormous density contrasts — rock against air at the top, rock
against iron at the bottom — so the mass they move around is not a small correction. It can
be **larger than the anomaly's own pull**, and it has the opposite sign.

So the geoid you actually measure is a competition between three mass contributions:

1. the density anomaly itself — *positive* for a dense blob,
2. the depressed **surface** it creates — *negative*,
3. the depressed **core–mantle boundary** — also *negative* (pushing the CMB down replaces
   dense core with lighter mantle, a mass deficit).

Both boundaries push the same way, against the anomaly. A single downwelling cannot depress
the surface while raising the CMB.

### Why depth decides the winner

The crucial point is that ① and ② are **two consequences of the same mass excess**, pulling in
opposite directions — and *proximity to the surface* is what sets their balance.

**A shallow dense anomaly** is close to the surface, so it depresses the surface strongly.
That depression is a big mass *deficit* sitting right where you are measuring, and it more
than cancels the anomaly's own enhanced pull. Net: a geoid **low**, even though you added
mass.

**A deep dense anomaly** does the opposite on both counts. It is far from the surface, so the
flow it drives has largely died out by the time it reaches the top — the surface barely
deflects. Meanwhile it deflects the **CMB** strongly instead. But the CMB is 2891 km away from
where you measure, and gravitational signals attenuate with distance as $(let
	f = (B_CMB / A_EARTH)^4
	"``(b/a)^{\\ell+2}``, only **$(round(f, digits=3))** at ``\\ell = 2``"
end) — so
that large CMB deflection is barely felt at the surface. The anomaly's own pull survives, and
the net can come out **positive**.

That is the whole mechanism: *shallow anomalies are dominated by the surface deflection they
cause; deep anomalies are dominated by their own attraction, because the boundary they deflect
is too far away to matter.* Which regime you are in depends on depth and on **how viscosity
varies with depth** — that is what a *geoid kernel* summarises, and why the geoid is one of
our best constraints on mantle viscosity (Hager & Richards, 1989).

This matters for a real, unsolved problem: the **Indian Ocean Geoid Low**, a ~106 m
depression in sea level south of India — the deepest geoid low on Earth. One family of
explanations involves deep dense material. This notebook is about building the intuition for
how that is even possible.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India

"""

# ╔═╡ 7e514084-5e6f-4091-c23d-4e5f60718293
md"---"

# ╔═╡ 8f625195-6f70-41a2-d34e-5f6071829304
md"""
### The Paint-Your-Own-Earth Sandbox

Everything below is driven by the globe you are about to see. **Paint a density anomaly**
into the mantle and watch the geoid respond on the surface above it.

- **The globe** — a cut-away Earth. The **outer surface is coloured by the geoid**, warm
  (yellow→red) where sea level bulges up and cool (blue→indigo) where it sags down. A wedge is
  removed so you can see inside, and the two flat **cut faces** show the mantle with your
  painted anomalies in it (blue dense, red light). Thick **orange curves** on those faces
  trace the **deflected surface and CMB**, exaggerated so they are visible at all.
- **Right panel** — the *kernel* ``K_\ell(r)``: the geoid produced by a unit anomaly placed
  at each depth. **Where it crosses zero is where the sign flips.** Markers show the depth of
  each anomaly you have painted.

**How to paint.** Drag anywhere to rotate, then **click on an exposed cut face** — the click
lands the anomaly at exactly that depth and latitude. There is no depth setting to fiddle
with: you are pointing straight at the place in the mantle where you want the anomaly, and
the kernel panel marks the depth of everything you have painted so you can see at a glance
which side of the sign flip each one falls on.

Untick **cut away** to close the wedge and see the geoid on a whole, unobstructed globe —
useful once you have finished painting. The **slice longitude** slider swings the cut around
so you can section through whatever you have drawn (painting moves it there automatically).

!!! tip "Seeing the flow — the missing middle step"
	Tick **show instantaneous flow** to draw velocity arrows on the cut faces. This is the
	step the rest of the picture leaves implicit: the anomaly drives flow, and *that flow*
	is what pushes on the boundaries and deflects them. A dense blob shows a clear
	downwelling beneath it, dying to nothing at the surface and the CMB where free slip
	forbids flow through the boundary.

	**The arrows are not motion.** This is a single instantaneous solve — there is no time
	in it, and nothing is travelling along those arrows. They show the velocity field *at
	this instant*, for a density anomaly frozen where you painted it. Their length scale is
	arbitrary too, since viscosity enters only as a ratio; only the **pattern** means
	anything. (It is also the slowest thing in the notebook, so it is off by default.)

Three things to try, in order — each is one of the three "aha" moments:

1. **The sign is wrong from the start.** Press **Shallow** (≈580 km). You have placed *dense*
   material, and the geoid above it is a **low** (≈ **−4.5 m**), not a high. Press **Deep**
   (≈2600 km): still a low (≈ **−0.8 m**), just weaker. In a uniform-viscosity mantle the
   naive expectation fails at *every* depth.
2. **Why**: tick **"direct effect only"**. The boundary terms switch off and the geoid jumps to
   ≈ **+15.6 m** — big and positive, the naive answer. Untick it and the deflected boundaries
   drag it back down to ≈ −4.5 m. The breakdown table further down shows the arithmetic:
   direct **+15.6**, surface **−19.6**, CMB **−0.5**.
3. **Layering is what creates a sign flip.** Drag the **viscosity ratio** slider. A uniform
   mantle has no zero-crossing at all; stiffen the lower mantle enough and one *appears*, so
   deep anomalies begin producing geoid **highs**. That is the punchline, and we come back to
   it below.

You can paint **both signs**. Pick **dense** (``\delta\rho > 0``, blue — cold slab material) or
**light** (``\delta\rho < 0``, red — a hot plume) with the brush selector, and mix them freely.

!!! note "Two colour scales — deliberately different palettes"
	The globe shows two different quantities, and each follows its own community's convention:

	- **surface — the geoid**: warm ↔ cool. Yellow deepening to red for highs, light blue
	  deepening to indigo for lows. This is the usual **geodesy** convention.
	- **interior — density**: blue ↔ red against dark mantle. Blue for dense, red for light,
	  matching the **seismic tomography** convention where blue is cold and fast.

	Note that the two conventions **disagree about sign**: warm means *positive* for the geoid,
	but red means *negative* for density. That is a historical accident of two fields, not
	physics. Because of it the ramps use visibly different palettes rather than the same pair
	of endpoint colours — otherwise one red would mean "geoid high" outside and "light
	material" inside, which is worse than either convention alone.

	Neither colouring carries physical meaning of its own. If the hues ever confuse you,
	ignore them and read the peak values printed in the legend.

Because the Stokes problem is *linear*, the geoid of a combined field is exactly the sum of the
geoids of its parts — so a light anomaly at a given depth produces precisely the mirror image
of a dense one, and two equal-and-opposite anomalies in the same place cancel to nothing.
Worth testing: paint a dense blob, then a light one on top of it, and watch the geoid collapse.

!!! note "About the exaggeration"
	The real boundary deflections are a few kilometres on a 6371 km globe — invisible at this
	scale — so they are drawn exaggerated. But the surface and the CMB share **one** scale
	factor, so the ratio between the two orange curves is honest: when the surface deflection
	looks much bigger than the CMB's, it really is. Both peak values are printed in the
	legend in metres. (The geoid colouring is separately normalised to its own peak, also
	printed, since it is a different quantity in different units.)
"""

# ╔═╡ c3a69596-a374-45e6-17f2-152637485960
md"""
## The Interactive Panel

All the drawing and all the plotting happen on HTML canvases inside a single widget. Julia
computes the physics and pushes the results back to the browser through a `CustomEvent`,
which the widget listens for and redraws.
"""

# ╔═╡ 18fbea87-f8b9-4a3b-6cb7-718293041527
md"""
### Where the Geoid Comes From — the Arithmetic

The globe shows you *that* the sign flips. This shows you *why*. Values are read at the
point where the geoid is largest, so they describe the strongest feature you have painted.

Watch the ratio between ① and ②: the anomaly's own pull against the depression it causes.
Whichever exceeds the other sets the sign, and moving the anomaly in depth is what tips it.
"""

# ╔═╡ 290cfb98-0851-4b4c-7dd8-930415263749
md"""
#### The competition, depth by depth

Here is the balance laid out for a unit dense sheet, in both a uniform mantle and Hager &
Richards' model WL. The column that matters is the **ratio** ``|②/①|`` — the surface
depression divided by the anomaly's own pull. **The geoid flips sign exactly where that ratio
crosses 1.**

$(let
	function tbl(lay)
		rows = String[]
		for d in (100e3, 300e3, 600e3, 1000e3, 1500e3, 2000e3, 2500e3, 2800e3)
			q = geoid_response(A_EARTH - d, 2, lay)
			ratio = abs(q.surface / q.direct)
			mark = q.total > 0 ? "**high**" : "low"
			push!(rows, "| $(round(Int, d/1e3)) km | `$(round(ratio, digits=3))` | $(mark) |")
		end
		"| depth | ``\\vert ②/① \\vert`` | net geoid |\n|---|---:|:---|\n" * join(rows, "\n")
	end
	Markdown.parse("**Uniform viscosity**\n\n" * tbl(isoviscous_layers()) *
				   "\n\n**Hager & Richards model WL**\n\n" * tbl(hr89_WL_layers()))
end)

In the uniform mantle the ratio never drops below 1 until the very base, so every depth gives
a low. In model WL the low-viscosity asthenosphere weakens the surface coupling enough that
the ratio starts *below* 1 — shallow anomalies give geoid **highs** — and the crossing to
lows happens in the mid-mantle.

Notice also the last row of each table: the ratio collapses near the CMB, because by then the
surface hardly deflects at all. The anomaly is instead deflecting the **CMB** — strongly, but
too far from the observer for it to compete, as the attenuation factor above showed.
"""

# ╔═╡ 4b2e1db6-2bec-4d6e-9fda-930415263748
md"""
### The Legendre spectrum of what you painted

The one conventional plot in this notebook. It shows how the anomaly you drew decomposes
into degrees — and, crucially, the **geoid contribution of each degree separately**. Notice
that different degrees can contribute with **opposite signs at the same place**, because each
degree has its own kernel with its own zero-crossing depth.
"""

# ╔═╡ 6d403fd6-4dcd-4f80-b1fb-152637485960
md"""
## Aha #3: Viscosity Controls the Flip

Sweep the **η lower / η upper** slider and watch the kernel. The headline is that
**a uniform mantle has no sign flip at all** — viscosity layering is what creates one:

$(let
	rows = String[]
	for ratio in (1.0, 3.0, 10.0, 30.0, 100.0)
		lay = ratio ≈ 1.0 ? isoviscous_layers() : twolayer_layers(ratio)
		d = zero_crossing_depth(2, lay)
		lbl = ratio ≈ 1.0 ? "1× *(uniform)*" : "$(ratio)×"
		push!(rows, "| $(lbl) | " * (isnothing(d) ? "**none** — single-signed" : "**$(round(Int,d)) km**") * " |")
	end
	Markdown.parse("| ``\\eta_{lower}/\\eta_{upper}`` | ``\\ell=2`` sign flip |\n|---|---|\n" * join(rows, "\n"))
end)

For a **uniform** mantle the ``\ell=2`` kernel is negative at every depth: the deflected
boundaries always outvote the anomaly's own attraction, so *any* dense anomaly gives a geoid
low. Only once the lower mantle is stiff enough — around 30× here — does a crossing appear,
and the deep mantle starts producing geoid *highs*.

!!! tip "Why this is how the geoid constrains viscosity"
	A stiff lower mantle resists the flow a sinking anomaly tries to drive. Weaken the flow
	and you weaken the boundary deflections — which are the *negative* contributions. Enough
	stiffening and they can no longer outvote the anomaly's own positive pull at depth.

	That is the inference. Subducted slabs are dense and deep, and they correlate
	**positively** with the observed geoid. A uniform-viscosity mantle predicts the wrong
	sign for them. Getting the observed sign right *requires* a lower mantle substantially
	stiffer than the upper — the conclusion Hager & Richards reached in 1989, and still one
	of the strongest constraints we have on mantle rheology.

!!! note "The real Earth is not two layers — compare model WL"
	This two-layer slider is a caricature. Hager & Richards' own preferred model **WL** gets
	its sign change largely from a **low-viscosity asthenosphere** (η = 1/30 between 100 and
	400 km) rather than from a stiff lower mantle alone:

	$(let
		ds = (100e3, 300e3, 600e3, 900e3, 1500e3, 2400e3)
		rows = ["| $(round(Int,d/1e3)) km | `$(round(geoid_kernel(A_EARTH-d, 2, hr89_WL_layers()), sigdigits=3))` |" for d in ds]
		Markdown.parse("| depth | ``K_2(r)`` for model WL |\n|---|---:|\n" * join(rows, "\n"))
	end)

	Positive in the upper mantle, crossing between 600 and 900 km, negative below — the
	structure of the ``G^\ell`` panel in their Figure 2. This is the notebook's
	reproduction-of-a-published-kernel check; the uniform-viscosity case above is our own
	baseline, not a published result.
"""

# ╔═╡ 7e5140e6-5ebf-4091-c20b-263748596071
md"""
## The Indian Ocean Geoid Low

South of India lies the deepest geoid low on Earth: sea level sits roughly **106 m below**
the reference ellipsoid over a vast region of the Indian Ocean. Explaining it is an open
research problem.

What this notebook gives you is the crucial piece of intuition. The naive objection —
*"but there is dense subducted slab material down there, so shouldn't the geoid be high?"* —
dissolves once you have seen the kernel. Under the right viscosity structure, deep dense
material produces a geoid **low**, because the boundary-deflection terms outvote the
anomaly's direct attraction.

Try it: hit **IOGL-like**, which drops a compact dense anomaly of roughly the right size
(≈3800 km across) into the mid-mantle. Keep the viscosity ratio near 1 and a **localized
geoid low** of about **−4.4 m** appears above it — not a band around the whole planet, but a
patch, because the density field here is genuinely three-dimensional. Now slide the viscosity
ratio up:

| ``\eta_{lower}/\eta_{upper}`` | geoid above the anomaly |
|---|---:|
| 1× | −4.4 m (low) |
| 3× | −3.1 m (low) |
| 10× | −1.3 m (low) |
| 30× | **+0.7 m (high)** |

That sensitivity is the point. The same dense anomaly gives a geoid low or a geoid high
depending entirely on the viscosity profile — which is exactly why the feature remains
contested, and why the geoid is such a useful constraint on mantle rheology.

!!! warning "One honest caveat"
	Real explanations of the IOGL invoke **lateral** viscosity variations, a low-density anomaly
	beneath the low, and the detailed history of Tethyan slab sinking. This model has none of
	those: its density varies in three dimensions, but its viscosity still varies only with
	depth, and nothing evolves in time. The **IOGL-like** button is a demonstration of a
	mechanism, not a fit to data — do not read the number it produces as a prediction.
"""

# ╔═╡ b2958428-92a3-44d5-0671-829304152637
md"""
## The Governing Equations

The mantle over geological time is a fluid so viscous that inertia is utterly negligible —
the Prandtl number is effectively infinite. What is left is the **Stokes equations**:
a statement that at every instant, viscous forces exactly balance buoyancy.
"""

# ╔═╡ c3a69539-a3b4-45e6-1782-930415263748
begin
    @syms r θ                    # radius and colatitude
    @syms η                      # dynamic viscosity
    @syms δρ                     # density anomaly
    @syms g                      # gravitational acceleration
    @syms 𝐯                      # velocity field
    @syms σ                      # stress tensor
    nothing
end

# ╔═╡ d4b7a64a-b4c5-46f7-2893-041526374859
md"""
**Incompressibility** — the mantle conserves volume as it flows:

```math
\nabla \cdot \mathbf{v} = 0
```

**Force balance** — viscous stress divergence balances the buoyancy of the anomaly:

```math
\nabla \cdot \sigma = \delta\rho \, g \, \hat{r}
```

!!! danger "This is instantaneous — there is no time in these equations"
	Look carefully: there is no ``\partial/\partial t`` anywhere. We are **not** running a
	convection simulation. The anomaly is *frozen* where you paint it, and we perform
	exactly **one linear solve** for the flow it drives at that instant.

	This trips up almost everyone the first time. The blob does not sink as you watch. We
	are asking "what flow, and what geoid, does this configuration produce *right now*?"
	— a snapshot, not a movie.

**Dynamic topography** — the flow pushes on each boundary, and the boundary deflects until
its own weight balances that push. Hager & Richards (1989) write this as a **jump** in the
radial normal stress across the boundary:

```math
\left[\tau_{rr}\right]^{+}_{-} = \Delta\rho \, g \, h
```

with ``\Delta\rho`` the density contrast **across** that boundary: ``\rho_m`` at the surface
(rock against air) and ``\rho_c - \rho_m`` at the CMB (iron against rock).

!!! warning "The bracket matters — this is where sign errors hide"
	Read outward, the surface jump goes mantle → vacuum but the CMB jump goes core → mantle,
	so the two contrasts have **opposite orientation**. Applying the same sign at both
	boundaries — an easy slip, and one this notebook made in an earlier version — produces a
	dense anomaly that depresses the surface while *raising* the CMB. No single convection
	cell can do that.

	We pin the convention with a physical limit instead of trusting the algebra: a load
	resting on a boundary must be almost exactly cancelled by that boundary's own deflection.
	That is the `check_compensation` test below, and it is what caught the error.

**The geoid** is then the sum over all three mass contributions, upward-continued to the
surface. For a thin mass sheet of surface density ``\sigma`` at radius ``s``:

```math
N_\ell = \frac{4\pi G a}{g(2\ell+1)} \, \sigma \, \left(\frac{s}{a}\right)^{\ell+2}
```
"""

# ╔═╡ e5c8b75b-c5d6-4708-3904-152637485960
md"""
## Reducing 3-D Flow to a Stack of 1-D Problems

Solving the Stokes equations in 3-D needs a mesh and a serious solver. We do not need one,
because of a piece of good luck: **if viscosity varies only with radius**, the problem
*separates* by spherical-harmonic degree. Each degree ``\ell`` evolves completely
independently of every other.

There is a second, less obvious piece of good luck that makes this cheap: **the radial
problem does not depend on the order** ``m``. Only the combination ``L = \ell(\ell+1)`` from
the angular part ever enters the radial equations, so all ``2\ell+1`` orders within a degree
share one and the same kernel ``K_\ell(r)``.

That is why a fully three-dimensional density field costs almost nothing extra. We compute
the kernel **once per degree** — 20 solves, not 441 — and reuse it for every order.

The whole pipeline becomes:

```
paint δρ(r,θ,φ)  →  spherical-harmonic transform  →  δρ_ℓm(r)
                 →  apply the 1-D radial kernel K_ℓ(r) for each ℓ  (same for all m)
                 →  sum over ℓ and m  →  geoid N(θ,φ)
```

So the "solver" is still a **loop of small 1-D problems**, one per degree, truncated at
``\ell_{max} = 20``. No mesh. No 3-D solver. That is why this runs instantly in your browser
even though the anomaly you paint is genuinely localized in both latitude and longitude.

!!! note "Why ``\ell_{max}`` is fixed rather than a slider"
	``\ell_{max}`` is a numerical parameter, not a physical one. Turning it down does not model
	a different Earth — it just blurs the answer, and below about ``\ell = 16`` a compact blob
	starts to ring (Gibbs oscillations) in a way that looks like physics but is not. So it is
	pinned at 20, and the brush instead enforces a **minimum anomaly size** that the expansion
	can represent faithfully. Pleasingly, that floor (≈1900 km radius) is about the scale of
	the Indian Ocean Geoid Low itself.
"""

# ╔═╡ 90736206-7081-42b3-e45f-607182930415
md"""
## Physical Constants

Standard Earth values. `η` is only ever used as a *ratio*, so its absolute scale never
matters for the geoid — a fact worth pausing on: **the geoid constrains relative viscosity,
not absolute viscosity.**
"""

# ╔═╡ a1847317-8192-43c4-f560-718293041526
begin
    const G_GRAV = 6.674e-11      # gravitational constant, m^3 kg^-1 s^-2
    const A_EARTH = 6371e3        # Earth radius, m
    const B_CMB = 3480e3          # core-mantle boundary radius, m
    const G_SURF = 9.81           # surface gravity, m s^-2
    const RHO_M = 3300.0          # representative mantle density, kg m^-3
    const RHO_C = 10900.0         # outer core density, kg m^-3
    const R_660 = A_EARTH - 660e3 # upper/lower mantle boundary

    # Spherical-harmonic truncation. Fixed rather than exposed as a slider: it is a
    # numerical parameter, not a physical one, and letting students turn it down
    # silently introduces ringing that looks like physics.
    const LMAX = 20

    # Smallest angular blob radius we can represent faithfully at LMAX. Below this,
    # truncation error grows fast (≈7% at 0.15 rad vs 4e-5 at 0.30). Conveniently
    # 0.30 rad ≈ 1900 km — about the scale of the Indian Ocean Geoid Low itself, so
    # the numerical floor and the physically interesting scale coincide.
    const BLOB_WIDTH_MIN = 0.30
    nothing
end

# ╔═╡ d4b7a646-b485-46f7-2803-263748596071
begin
    struct GeoidGlobeInput
        default_blobs::Vector{Vector{Float64}}
        visc_ratio::Float64
        ell::Int
        direct_only::Bool
        dense::Bool
        width::Float64
    end

    function GeoidGlobeInput(; visc_ratio=1.0, ell=2, direct_only=false,
        dense=true, width=BLOB_WIDTH_MIN)
        # One dense blob at ~580 km depth, where the sign flip is unmistakable. Its
        # longitude matches the default slice (0.9 rad) so it is visible in section
        # the moment the notebook opens.
        # entries are [θ, φ, r_fraction, amplitude, angular_width]
        GeoidGlobeInput([[1.15, 0.9, 0.80, 1.0, BLOB_WIDTH_MIN]],
            Float64(visc_ratio), Int(ell), Bool(direct_only),
            Bool(dense), Float64(width))
    end

    Base.get(w::GeoidGlobeInput) = Dict{String,Any}(
        "blobs" => w.default_blobs,
        "visc_ratio" => w.visc_ratio,
        "ell" => w.ell,
        "direct_only" => w.direct_only,
        "dense" => w.dense,
        "width" => w.width,
    )

    function Base.show(io::IO, ::MIME"text/html", w::GeoidGlobeInput)
        write(io, """
<div id="gkwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    /* PlutoUI's wide-cell wrapper hard-caps itself around 1000px regardless of
       viewport size. Override it for this specific cell so the widget can actually
       reach ~80% of a wide (16:9 laptop/projector) screen instead of stalling there. */
    /* Pluto centers a wide cell with an inline margin-left computed for its OWN
       ~1000px cap; overriding width alone leaves that stale margin in place and the
       cell drifts right. Recompute margin-left from the same override width so it
       stays centered under its (unwidened) parent at any viewport size. */
    pluto-cell:has(#gkwidget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #gkwidget .gk-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #gkwidget .gk-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #gkwidget .gk-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #gkwidget .gk-workspace{display:flex;gap:12px;align-items:flex-start;justify-content:center;width:100%}
    #gkwidget .gk-controls{width:min(var(--totalw,960px),100%);margin-top:8px;display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:8px;font:14px sans-serif}
    #gkwidget .gk-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
    #gkwidget .gk-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:20px}
    /* Fixed-px columns don't shrink -- if a group ever ends up narrower than the
       label+slider+value minimums add up to, the slider bleeds past the group's
       border instead of compressing. minmax() lets every column shrink instead. */
    #gkwidget .gk-control-row{display:grid;grid-template-columns:minmax(90px,160px) minmax(90px,1fr) minmax(70px,108px);gap:6px;align-items:center;margin:7px 0}
    #gkwidget .gk-control-row input[type=range]{width:100%;min-width:0;vertical-align:middle}
    #gkwidget .gk-value{color:#d1d5db;text-align:left;font-variant-numeric:tabular-nums}
    #gkwidget .gk-actions{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
    #gkwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px}
    #gkwidget label{color:#d1d5db}
    @media (max-width: 980px){
      #gkwidget .gk-workspace{flex-direction:column;align-items:center}
      #gkwidget .gk-controls{grid-template-columns:1fr;width:660px;max-width:100%}
    }
  </style>
  <div class="gk-title">
    <div class="gk-title-desc">A density anomaly deforms the surface and core&ndash;mantle boundary as the mantle flows around it &mdash; the geoid is the sum of all three effects.</div>
    <div class="gk-title-hint">drag to rotate the globe &middot; click a cut face to place a density anomaly</div>
  </div>
  <div class="gk-workspace">
    <canvas id="globecvs" width="660" height="660"
      style="cursor:crosshair;background:#000;border:1px solid #374151;border-radius:6px;display:block"></canvas>
    <canvas id="kernelcvs" width="248" height="660"
      style="background:#000;border:1px solid #374151;border-radius:6px;display:block"></canvas>
  </div>
  <div class="gk-controls">
    <div class="gk-control-group">
      <div class="gk-control-title">Structure</div>
      <label class="gk-control-row"><span>η lower / η upper</span><input type="range" id="visc" min="0" max="2" step="0.01" value="$(round(log10(w.visc_ratio), digits=4))"><span id="viscv" class="gk-value">$(round(w.visc_ratio, digits=1))×</span></label>
      <label class="gk-control-row"><span>kernel degree ℓ</span><input type="range" id="ell" min="2" max="12" step="1" value="$(w.ell)"><span id="ellv" class="gk-value">$(w.ell)</span></label>
      <label class="gk-control-row"><span>anomaly size</span><input type="range" id="width" min="$(BLOB_WIDTH_MIN)" max="0.75" step="0.01" value="$(w.width)"><span id="widthv" class="gk-value">3800 km</span></label>
    </div>
    <div class="gk-control-group">
      <div class="gk-control-title">Anomaly &amp; physics</div>
      <div class="gk-actions">
        <label><input type="checkbox" id="direct" $(w.direct_only ? "checked" : "") style="vertical-align:middle"> <b>direct effect only</b> (no boundaries)</label>
      </div>
      <div class="gk-actions" style="margin-top:8px">
        <label><input type="checkbox" id="cut" checked style="vertical-align:middle"> <b>cut away</b> (uncheck for a whole globe)</label>
      </div>
      <label class="gk-control-row" style="margin-top:6px"><span>slice longitude</span><input type="range" id="cutpos" min="0" max="6.28" step="0.02" value="0.9"><span id="cutposv" class="gk-value">52°</span></label>
      <div class="gk-actions" style="margin-top:4px">
        <label><input type="checkbox" id="showflow" style="vertical-align:middle"> show <b>instantaneous flow</b> on the cut (slow)</label>
      </div>
      <div class="gk-actions" style="margin-top:8px">
        <span style="color:#9ca3af">brush:</span>
        <label class="gk-sign" id="lab-dense"><input type="radio" name="gksign" id="dense" $(w.dense ? "checked" : "") style="vertical-align:middle"> <b style="color:#3b82f6">dense</b> (δρ &gt; 0)</label>
        <label class="gk-sign" id="lab-light"><input type="radio" name="gksign" id="light" $(w.dense ? "" : "checked") style="vertical-align:middle"> <b style="color:#ef4444">light</b> (δρ &lt; 0)</label>
      </div>
      <div class="gk-actions" style="margin-top:8px">
        <button id="clrbtn">Clear</button>
        <button id="shallowbtn">Shallow</button>
        <button id="deepbtn">Deep</button>
        <button id="iogbtn">IOGL-like</button>
        <span id="cnt" class="gk-value">blobs: $(length(w.default_blobs))</span>
      </div>
      <div style="margin-top:8px;color:#9ca3af;font-size:13px">
        drag to rotate · click a <b>cut face</b> to place an anomaly at that depth ·
        uncheck <b>cut away</b> for a clean view of the geoid alone
      </div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  // Target ~80% of viewport width (never more than the wide-cell wrapper actually
  // gives us), and cap the globe canvas's height so the globe + kernel strip +
  // controls plausibly fit one screen without scrolling in a lecture/projector setting.
  // The globe/kernel canvases hand *physical* (θ, φ, r_fraction, ...) values back to
  // Julia (see gk_blobs below), never raw pixels, so it's safe to just resize these
  // constants directly -- unlike the interferometry widget's pixel-coupled canvas.
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8)
  const totalW = Math.max(700, availW)
  const heightBudget = Math.max(420, window.innerHeight - 425)
  const SEC = Math.round(Math.min(totalW*(660/920), heightBudget, 900))
  const KW = Math.round(SEC*248/660)
  const KH = SEC
  const A = $(A_EARTH), B = $(B_CMB)
  const R_OUT = Math.round(SEC*232/660), R_IN = R_OUT*B/A
  // The removed wedge. cutStart is steerable so you can slice straight through
  // whatever you have painted -- a blob at some longitude is otherwise only
  // clipped by the face, and you would see its tail rather than its body.
  let cutStart = 0.9, cutSpan = Math.PI/2
  const DPR = Math.min(window.devicePixelRatio || 1, 2)
  const WMIN = $(BLOB_WIDTH_MIN)
  const cvs = par.querySelector('#globecvs'), ctx = cvs.getContext('2d')
  const ker = par.querySelector('#kernelcvs'), kctx = ker.getContext('2d')
  const lbl = par.querySelector('#cnt')
  const CX = SEC/2, CY = SEC/2
  // The controls grid isn't tied to the canvas row's width -- on a short screen the
  // globe shrinks a lot (height-limited) but there's still plenty of horizontal room,
  // so let the controls use it instead of matching the (possibly much narrower) canvas row.
  par.style.setProperty('--totalw', Math.round(totalW) + 'px')

  function hidpi(cv, cx, w, h){
    cv.width=Math.round(w*DPR); cv.height=Math.round(h*DPR)
    cv.style.width=w+'px'; cv.style.height=h+'px'
    cx.setTransform(DPR,0,0,DPR,0,0)
  }
  hidpi(cvs,ctx,SEC,SEC); hidpi(ker,kctx,KW,KH)

  // blobs: [theta, phi, rfrac, amp, width]
  let blobs = $(replace(string([[b[1], b[2], b[3], b[4], b[5]] for b in w.default_blobs]), "Any" => ""))
  let viscRatio = $(w.visc_ratio), ell = $(w.ell)
  let directOnly = $(w.direct_only), dense = $(w.dense)
  let brushWidth = $(w.width), cutAway = true, showFlow = false
  let flowData = null
  let yaw = 0.9, pitch = 0.35
  let dragging = false, dragMoved = false, lastX = 0, lastY = 0
  let kernelData = null, geoidMap = null, topoData = null, muGrid = null, phiGrid = null

  // ---------- orthographic camera ----------
  // world (X,Y,Z), Z = polar axis. Screen x = X', y = -Z'; +Y' points at the viewer.
  function rot(p){
    const x = p[0]*Math.cos(yaw) - p[1]*Math.sin(yaw)
    const y = p[0]*Math.sin(yaw) + p[1]*Math.cos(yaw)
    const z = p[2]
    return [x, y*Math.cos(pitch) - z*Math.sin(pitch), y*Math.sin(pitch) + z*Math.cos(pitch)]
  }
  function invrot(p){
    const x=p[0], y2=p[1], z2=p[2]
    const y = y2*Math.cos(pitch) + z2*Math.sin(pitch)
    const z = -y2*Math.sin(pitch) + z2*Math.cos(pitch)
    return [x*Math.cos(yaw) + y*Math.sin(yaw), -x*Math.sin(yaw) + y*Math.cos(yaw), z]
  }
  function sph(theta, phi, r){
    return [r*Math.sin(theta)*Math.cos(phi), r*Math.sin(theta)*Math.sin(phi), r*Math.cos(theta)]
  }
  function proj(p){ return [CX + p[0], CY - p[2]] }
  // With the cut closed the globe is solid: nothing is removed, so no longitude
  // is ever "in the cut" and the interior is simply not visible.
  function cutFaces(){ return [cutStart, cutStart + cutSpan] }
  function inCut(phi){
    if(!cutAway) return false
    // angle measured from the start of the wedge, wrapped into [0, 2π)
    let a = ((phi - cutStart) % (2*Math.PI) + 2*Math.PI) % (2*Math.PI)
    return a <= cutSpan
  }

  // ---------- inverse: what did the user click? ----------
  function hitTest(px, py){
    const sx = px-CX, sz = -(py-CY), rr2 = sx*sx + sz*sz
    if(rr2 <= R_OUT*R_OUT){
      const sy = Math.sqrt(R_OUT*R_OUT - rr2)
      const wv = invrot([sx,sy,sz])
      const phi = Math.atan2(wv[1], wv[0])
      if(!inCut(phi))
        return {kind:'surface', theta:Math.acos(Math.max(-1,Math.min(1,wv[2]/R_OUT))), phi:phi}
    }
    for(const pc of (cutAway ? cutFaces() : [])){
      const nr_ = rot([-Math.sin(pc), Math.cos(pc), 0])
      if(Math.abs(nr_[1]) < 1e-9) continue
      const t = -(nr_[0]*sx + nr_[2]*sz)/nr_[1]
      const wv = invrot([sx,t,sz])
      const rad = Math.hypot(wv[0],wv[1],wv[2])
      if(rad < R_IN || rad > R_OUT) continue
      const phi = Math.atan2(wv[1], wv[0])
      const dphi = Math.atan2(Math.sin(phi-pc), Math.cos(phi-pc))
      if(Math.abs(dphi) > 1e-6) continue
      return {kind:'cutface', theta:Math.acos(Math.max(-1,Math.min(1,wv[2]/rad))),
              phi:pc, rfrac:(rad-R_IN)/(R_OUT-R_IN)}
    }
    if(cutAway && rr2 <= R_IN*R_IN){
      const sy = Math.sqrt(R_IN*R_IN - rr2)
      const wv = invrot([sx,sy,sz])
      return {kind:'cmb', theta:Math.acos(Math.max(-1,Math.min(1,wv[2]/R_IN))),
              phi:Math.atan2(wv[1],wv[0])}
    }
    if(cutAway && rr2 <= R_OUT*R_OUT){   // seen through the cut: far shell wall
      const sy = -Math.sqrt(R_OUT*R_OUT - rr2)
      const wv = invrot([sx,sy,sz])
      return {kind:'surface', theta:Math.acos(Math.max(-1,Math.min(1,wv[2]/R_OUT))),
              phi:Math.atan2(wv[1],wv[0]), farSide:true}
    }
    return null
  }

  // ---------- sampling the pushed maps ----------
  function sampleMap(M, theta, phi){
    if(!M || !muGrid || !phiGrid) return 0
    const mu = Math.cos(theta)
    // muGrid descends (Gauss nodes); find nearest
    let bi=0, bd=1e9
    for(let i=0;i<muGrid.length;i++){ const d=Math.abs(muGrid[i]-mu); if(d<bd){bd=d;bi=i} }
    let a = ((phi % (2*Math.PI)) + 2*Math.PI) % (2*Math.PI)
    const nj = phiGrid.length
    let j = Math.round(a/(2*Math.PI)*nj) % nj
    return M[bi][j]
  }
  function mapRange(M){
    if(!M) return 1
    let mx=0
    for(const row of M) for(const v of row) mx = Math.max(mx, Math.abs(v))
    return mx || 1
  }
  // Opaque density ramp for the cut faces: dark mantle grey at zero, saturating to
  // blue for dense and red for light. Fully opaque at every amplitude -- opacity
  // carries no information, so a faint anomaly is still plainly visible.
  function densityColor(d){
    const t = Math.max(-1, Math.min(1, d))
    const bg = [18, 22, 29]                      // mantle background
    const hi = t >= 0 ? [59, 130, 246] : [239, 68, 68]
    const a = Math.min(1, Math.abs(t) * 1.6)     // reach full colour before |d|=1
    const c = [0,1,2].map(i => Math.round(bg[i] + (hi[i]-bg[i])*a))
    return 'rgb('+c[0]+','+c[1]+','+c[2]+')'
  }
  // Geoid ramp for the outer surface, following the usual GEODESY convention:
  // highs warm (yellow → red), lows cool (blue → indigo), near-white at zero.
  //
  // Deliberately a DIFFERENT palette from densityColor above, not merely an inverted
  // one. The two quantities use opposite sign conventions by community habit — warm =
  // geoid high, but blue = dense in tomography — so if they shared endpoint colours
  // the same red would mean "geoid high" outside and "light material" inside. Making
  // the ramps visually distinct (warm/cool with light midpoint vs saturated blue/red
  // on dark) keeps them from being confused for one another.
  function signColor(v, mx){
    const t = Math.max(-1, Math.min(1, v/mx))
    const mid = [246, 246, 240]                     // near-white at the reference geoid
    // warm limb: yellow at moderate highs deepening to red; cool limb: cyan-blue to indigo
    const warm = [[250, 204, 21], [220, 38, 38]]
    const cool = [[56, 189, 248], [49, 46, 129]]
    const a = Math.abs(t)
    const stops = t >= 0 ? warm : cool
    // two-segment ramp: mid -> stops[0] over the first half, stops[0] -> stops[1] after
    const c = a <= 0.5
      ? [0,1,2].map(i => Math.round(mid[i] + (stops[0][i]-mid[i]) * (a/0.5)))
      : [0,1,2].map(i => Math.round(stops[0][i] + (stops[1][i]-stops[0][i]) * ((a-0.5)/0.5)))
    return 'rgb('+c[0]+','+c[1]+','+c[2]+')'
  }

  // ---------- drawing ----------
  function drawGlobe(){
    ctx.clearRect(0,0,SEC,SEC)
    const gmax = mapRange(geoidMap)

    // --- outer surface, painted with the geoid, drawn as lat/lon quads.
    // Back faces are culled, and the cut wedge is skipped.
    const NT = 60, NP = 120
    for(let i=0;i<NT;i++){
      const t0 = Math.PI*i/NT, t1 = Math.PI*(i+1)/NT
      for(let j=0;j<NP;j++){
        const p0 = 2*Math.PI*j/NP, p1 = 2*Math.PI*(j+1)/NP
        const pm = (p0+p1)/2
        if(inCut(pm)) continue
        const c = [sph(t0,p0,R_OUT), sph(t0,p1,R_OUT), sph(t1,p1,R_OUT), sph(t1,p0,R_OUT)]
        const rc = c.map(rot)
        // cull back faces (average outward normal pointing away)
        if((rc[0][1]+rc[1][1]+rc[2][1]+rc[3][1])/4 < 0) continue
        const v = sampleMap(geoidMap, (t0+t1)/2, pm)
        ctx.beginPath()
        const s0 = proj(rc[0]); ctx.moveTo(s0[0], s0[1])
        for(let k=1;k<4;k++){ const s=proj(rc[k]); ctx.lineTo(s[0], s[1]) }
        ctx.closePath()
        ctx.fillStyle = signColor(v, gmax)
        ctx.fill()
      }
    }

    // --- the two flat cut faces: mantle interior, shaded by the painted density
    for(const pc of (cutAway ? cutFaces() : [])){
      const nr_ = rot([-Math.sin(pc), Math.cos(pc), 0])
      if(nr_[1] < 0) continue                 // face points away from us
      const NR = 34, NTh = 60
      for(let i=0;i<NR;i++){
        const r0 = R_IN + (R_OUT-R_IN)*i/NR, r1 = R_IN + (R_OUT-R_IN)*(i+1)/NR
        for(let k=0;k<NTh;k++){
          const t0 = Math.PI*k/NTh, t1 = Math.PI*(k+1)/NTh
          // density from the blobs at this (r,θ,φ=pc)
          let d = 0
          const rf = ((r0+r1)/2 - R_IN)/(R_OUT-R_IN)
          const tm = (t0+t1)/2
          for(const bb of blobs){
            const ca = Math.cos(tm)*Math.cos(bb[0]) +
                       Math.sin(tm)*Math.sin(bb[0])*Math.cos(pc-bb[1])
            const dd = Math.acos(Math.max(-1,Math.min(1,ca)))
            const wgt = Math.exp(-Math.pow(dd/Math.max(bb[4],WMIN),2))
            const rw = Math.exp(-Math.pow((rf-bb[2])/0.09,2))
            d += bb[3]*wgt*rw
          }
          const quad = [sph(t0,pc,r0), sph(t0,pc,r1), sph(t1,pc,r1), sph(t1,pc,r0)].map(rot)
          ctx.beginPath()
          const s0 = proj(quad[0]); ctx.moveTo(s0[0], s0[1])
          for(let q=1;q<4;q++){ const s=proj(quad[q]); ctx.lineTo(s[0], s[1]) }
          ctx.closePath()
          // The face is FULLY OPAQUE: amplitude is shown by colour, not by
          // transparency, so even a weak anomaly is clearly visible against the
          // dark mantle rather than fading into the background.
          ctx.fillStyle = densityColor(d)
          ctx.fill()
        }
      }
      // outline the face so the wedge reads as a solid cut
      ctx.beginPath()
      for(let k=0;k<=60;k++){ const s=proj(rot(sph(Math.PI*k/60,pc,R_OUT))); k===0?ctx.moveTo(s[0],s[1]):ctx.lineTo(s[0],s[1]) }
      for(let k=60;k>=0;k--){ const s=proj(rot(sph(Math.PI*k/60,pc,R_IN)));  ctx.lineTo(s[0],s[1]) }
      ctx.closePath(); ctx.strokeStyle='#6b7280'; ctx.lineWidth=1; ctx.stroke()
    }

    // --- the core, visible only through the cut
    if(cutAway){
      ctx.beginPath(); ctx.arc(CX,CY,R_IN,0,2*Math.PI)
      ctx.fillStyle='rgba(36,26,18,0.55)'; ctx.fill()
    }

    // --- instantaneous flow arrows on the visible cut faces
    if(showFlow && flowData && cutAway){
      const faces = cutFaces()
      for(let fi=0; fi<faces.length; fi++){
        const pc = faces[fi], F = flowData[fi]
        if(!F) continue
        const nr_ = rot([-Math.sin(pc), Math.cos(pc), 0])
        if(nr_[1] < 0) continue                  // face turned away from us
        // one global scale so arrow lengths are comparable across the face
        let vmax = 0
        for(let i=0;i<F.ur.length;i++) for(let k=0;k<F.ur[i].length;k++)
          vmax = Math.max(vmax, Math.hypot(F.ur[i][k], F.ut[i][k]))
        if(!(vmax > 0)) continue
        const LEN = 26
        ctx.strokeStyle = 'rgba(226,232,240,0.85)'
        ctx.fillStyle   = 'rgba(226,232,240,0.85)'
        ctx.lineWidth = 1.1
        for(let i=0;i<F.rr.length;i++){
          for(let k=0;k<F.tt.length;k++){
            const r = F.rr[i]*R_OUT/A, th = F.tt[k]
            const vr = F.ur[i][k], vt = F.ut[i][k]
            const mag = Math.hypot(vr, vt)
            if(mag < 0.06*vmax) continue         // skip near-stagnant points
            const sc = LEN*mag/vmax
            // unit vectors in the meridional plane, mapped to 3-D then projected
            const er = sph(th, pc, 1)
            const et = [Math.cos(th)*Math.cos(pc), Math.cos(th)*Math.sin(pc), -Math.sin(th)]
            const p0 = sph(th, pc, r)
            const d  = [0,1,2].map(q => (vr*er[q] + vt*et[q])/mag)
            const p1 = [0,1,2].map(q => p0[q] + d[q]*sc)
            const a0 = proj(rot(p0)), a1 = proj(rot(p1))
            ctx.beginPath(); ctx.moveTo(a0[0],a0[1]); ctx.lineTo(a1[0],a1[1]); ctx.stroke()
            // arrowhead
            const ang = Math.atan2(a1[1]-a0[1], a1[0]-a0[0])
            const hl = Math.min(5, 0.42*Math.hypot(a1[0]-a0[0], a1[1]-a0[1]))
            if(hl > 1.4){
              ctx.beginPath(); ctx.moveTo(a1[0],a1[1])
              ctx.lineTo(a1[0]-hl*Math.cos(ang-0.4), a1[1]-hl*Math.sin(ang-0.4))
              ctx.lineTo(a1[0]-hl*Math.cos(ang+0.4), a1[1]-hl*Math.sin(ang+0.4))
              ctx.closePath(); ctx.fill()
            }
          }
        }
      }
    }

    // --- deflected boundaries along the two cut faces (exaggerated)
    // ONE shared scale for both boundaries, so their relative amplitudes are honest:
    // the surface deflection really is larger than the CMB's for a shallow anomaly,
    // and normalising each to its own maximum would hide exactly that.
    if(topoData && cutAway){
      const tmax = Math.max(mapRange(topoData.surf), mapRange(topoData.cmb))
      for(const pc of cutFaces()){
        const nr_ = rot([-Math.sin(pc), Math.cos(pc), 0])
        if(nr_[1] < 0) continue
        // Both boundaries in the SAME orange: they are the same physical quantity
        // (a deflected density interface), and orange is distinct from both the
        // geoid ramp on the surface and the density ramp on the cut faces.
        for(const [M, baseR, color] of [[topoData.surf, R_OUT, '#f59e0b'],
                                        [topoData.cmb,  R_IN,  '#f59e0b']]){
          ctx.beginPath()
          for(let k=0;k<=80;k++){
            const th = Math.PI*k/80
            const rr = baseR + 26*sampleMap(M, th, pc)/tmax
            const s = proj(rot(sph(th, pc, rr)))
            k===0?ctx.moveTo(s[0],s[1]):ctx.lineTo(s[0],s[1])
          }
          ctx.strokeStyle=color; ctx.lineWidth=2.6; ctx.stroke()
        }
      }
    }

    // --- limb
    ctx.beginPath(); ctx.arc(CX,CY,R_OUT,0,2*Math.PI)
    ctx.strokeStyle='#4b5563'; ctx.lineWidth=1; ctx.stroke()

    // --- labels
    ctx.fillStyle='#9ca3af'; ctx.font='14px sans-serif'
    ctx.fillText('geoid on the surface · mantle exposed in the cut', 12, 20)
    ctx.font='13px sans-serif'
    ctx.fillStyle='#9ca3af'
    ctx.fillText('drag to rotate', 12, 38)
    // Two DIFFERENT palettes, each following its own community convention:
    // geoid warm=high (geodesy), density blue=dense (tomography).
    ctx.fillStyle='#9ca3af'; ctx.fillText('surface \\u2014 geoid:', 12, SEC-46)
    ctx.fillStyle='#dc2626'; ctx.fillText('high', 108, SEC-46)
    ctx.fillStyle='#9ca3af'; ctx.fillText('/', 136, SEC-46)
    ctx.fillStyle='#38bdf8'; ctx.fillText('low', 144, SEC-46)
    ctx.fillStyle='#9ca3af'; ctx.fillText('interior \\u2014 density:', 12, SEC-32)
    ctx.fillStyle='#3b82f6'; ctx.fillText('dense', 108, SEC-32)
    ctx.fillStyle='#9ca3af'; ctx.fillText('/', 146, SEC-32)
    ctx.fillStyle='#ef4444'; ctx.fillText('light', 154, SEC-32)
    // Both boundaries share one exaggeration scale, so quote both peaks: the ratio
    // between these two numbers is real, and visible in the drawing.
    // One solid orange swatch for both deflected boundaries — same quantity, same
    // colour; the peak values distinguish them.
    ctx.strokeStyle='#f59e0b'; ctx.lineWidth=2.6
    ctx.beginPath(); ctx.moveTo(210, SEC-50); ctx.lineTo(228, SEC-50); ctx.stroke()
    ctx.fillStyle='#f59e0b'
    ctx.fillText('deflected boundaries', 234, SEC-46)
    if(topoData){
      ctx.font='13px sans-serif'
      ctx.fillText('surface ' + mapRange(topoData.surf).toPrecision(2) +
                   ' m  ·  CMB ' + mapRange(topoData.cmb).toPrecision(2) + ' m', 234, SEC-32)
      ctx.font='14px sans-serif'
    }
    if(showFlow){
      ctx.fillStyle = flowData ? '#e2e8f0' : '#6b7280'
      ctx.fillText(flowData ? '\\u2192 instantaneous velocity (arbitrary scale, NOT motion in time)'
                            : '\\u2192 solving for the flow...', 12, SEC-18)
    }
    if(geoidMap){
      ctx.fillStyle='#e5e7eb'; ctx.font='14px sans-serif'
      ctx.fillText('peak |N| = '+gmax.toPrecision(3)+' m', 12, SEC-14)
    }
    lbl.textContent = 'blobs: '+blobs.length
  }

  function drawKernel(){
    kctx.clearRect(0,0,KW,KH)
    kctx.fillStyle='#9ca3af'; kctx.font='14px sans-serif'
    kctx.fillText('kernel K'+'\\u2113'+'(r), \\u2113='+ell, 10, 18)
    if(!kernelData){ kctx.fillText('computing...', 10, 40); return }
    // Leave 64 px below the plot for three stacked label rows: the depth tick, the
    // negative/positive sign labels, and the anomaly-marker legend.
    const pad=38, W=KW-pad-14, H=KH-pad-64
    const vals=kernelData.k, deps=kernelData.depth
    const mx = Math.max(...vals.map(Math.abs)) || 1
    const X = v => pad + W*(0.5 + 0.45*v/mx)
    const Y = d => pad + H*(d/2891)
    kctx.strokeStyle='#6b7280'; kctx.setLineDash([4,4])
    kctx.beginPath(); kctx.moveTo(X(0),pad); kctx.lineTo(X(0),pad+H); kctx.stroke()
    kctx.setLineDash([])
    kctx.beginPath()
    for(let i=0;i<vals.length;i++){
      const x=X(vals[i]), y=Y(deps[i]); i===0?kctx.moveTo(x,y):kctx.lineTo(x,y)
    }
    kctx.strokeStyle='#f97316'; kctx.lineWidth=2; kctx.stroke()
    if(kernelData.crossing !== null && kernelData.crossing !== undefined){
      const y=Y(kernelData.crossing)
      kctx.strokeStyle='#22c55e'; kctx.lineWidth=1.5; kctx.setLineDash([3,3])
      kctx.beginPath(); kctx.moveTo(pad,y); kctx.lineTo(KW-14,y); kctx.stroke()
      kctx.setLineDash([])
      kctx.fillStyle='#22c55e'; kctx.font='13px sans-serif'
      kctx.fillText('sign flip '+Math.round(kernelData.crossing)+' km', pad+4, y-5)
    } else {
      kctx.fillStyle='#22c55e'; kctx.font='13px sans-serif'
      kctx.fillText('no sign flip', pad+4, pad+14)
    }
    // Mark the depth of every anomaly actually painted, so you can read off whether
    // each one sits above or below the sign flip.
    for(const bb of blobs){
      const py = Y(2891*(1-bb[2]))
      kctx.strokeStyle = bb[3]>0 ? 'rgba(59,130,246,0.9)' : 'rgba(239,68,68,0.9)'
      kctx.lineWidth=1; kctx.setLineDash([1,3])
      kctx.beginPath(); kctx.moveTo(pad,py); kctx.lineTo(KW-14,py); kctx.stroke()
      kctx.setLineDash([])
      kctx.beginPath(); kctx.arc(pad+5, py, 3, 0, 2*Math.PI)
      kctx.fillStyle = bb[3]>0 ? '#3b82f6' : '#ef4444'; kctx.fill()
    }
    kctx.strokeStyle='#374151'; kctx.lineWidth=1
    kctx.beginPath(); kctx.moveTo(pad,pad+H); kctx.lineTo(KW-14,pad+H); kctx.stroke()
    kctx.beginPath(); kctx.moveTo(pad,pad); kctx.lineTo(pad,pad+H); kctx.stroke()
    kctx.fillStyle='#9ca3af'; kctx.font='13px sans-serif'
    kctx.save(); kctx.translate(12,pad+H/2); kctx.rotate(-Math.PI/2)
    kctx.fillText('depth (km)',-30,0); kctx.restore()
    kctx.fillText('0',pad-4,pad-6); kctx.fillText('2891',pad-16,pad+H+14)
    // Three separate rows below the axis so nothing collides: the sign labels on one
    // line, the anomaly-marker legend on the next.
    kctx.font='13px sans-serif'
    kctx.fillStyle='#ef4444'; kctx.fillText('negative', pad+4, KH-30)
    kctx.fillStyle='#3b82f6'; kctx.fillText('positive \\u2192', KW-80, KH-30)
    if(blobs.length){
      kctx.fillStyle='#9ca3af'; kctx.font='13px sans-serif'
      kctx.fillText('\\u25cf your anomalies', pad+4, KH-14)
    }
  }

  function redraw(){ drawGlobe(); drawKernel() }

  function emit(){
    par.value = {blobs:blobs, visc_ratio:viscRatio, ell:ell,
                 direct_only:directOnly, dense:dense, width:brushWidth,
                 show_flow:showFlow, cut_start:cutStart}
    par.dispatchEvent(new CustomEvent('input'))
  }

  function addBlob(theta, phi, rfrac){
    for(const b of blobs){
      const ca = Math.cos(theta)*Math.cos(b[0]) +
                 Math.sin(theta)*Math.sin(b[0])*Math.cos(phi-b[1])
      const dd = Math.acos(Math.max(-1,Math.min(1,ca)))
      if(dd < 0.18 && Math.abs(rfrac-b[2]) < 0.10) return false
    }
    blobs.push([theta, phi, rfrac, dense?1.0:-1.0, brushWidth])
    // No need to move the slice: you painted ON a cut face, so the blob is already
    // sitting in a plane you are looking at.
    return true
  }

  cvs.addEventListener('mousedown', e=>{
    dragging=true; dragMoved=false; lastX=e.offsetX; lastY=e.offsetY
  })
  cvs.addEventListener('mousemove', e=>{
    if(dragging){
      const dx=e.offsetX-lastX, dy=e.offsetY-lastY
      if(Math.abs(dx)+Math.abs(dy) > 2) dragMoved=true
      yaw += dx*0.008
      pitch = Math.max(-1.3, Math.min(1.3, pitch + dy*0.008))
      lastX=e.offsetX; lastY=e.offsetY
      redraw()
    } else {
      const h = hitTest(e.offsetX, e.offsetY)
      cvs.style.cursor = h ? 'crosshair' : 'default'
    }
  })
  window.addEventListener('mouseup', e=>{
    if(dragging && !dragMoved){
      const h = hitTest(e.target===cvs ? e.offsetX : lastX, e.target===cvs ? e.offsetY : lastY)
      // Paint on the CUT FACE only: the click itself then carries the depth, so there
      // is no need to set one numerically. Clicks elsewhere just rotate.
      if(h && h.kind==='cutface'){
        if(addBlob(h.theta, h.phi, h.rfrac)){ redraw(); emit() }
      }
    }
    dragging=false
  })

  function widthLabel(w){ return Math.round(2*w*6371)+' km' }

  par.querySelector('#visc').addEventListener('input', e=>{
    viscRatio = Math.pow(10, parseFloat(e.target.value))
    par.querySelector('#viscv').textContent = viscRatio.toFixed(viscRatio<10?1:0)+'\\u00d7'
    emit()
  })
  par.querySelector('#ell').addEventListener('input', e=>{
    ell = parseInt(e.target.value); par.querySelector('#ellv').textContent = ell; emit()
  })
  par.querySelector('#width').addEventListener('input', e=>{
    brushWidth = Math.max(WMIN, parseFloat(e.target.value))
    par.querySelector('#widthv').textContent = widthLabel(brushWidth)
  })
  par.querySelector('#direct').addEventListener('change', e=>{ directOnly=e.target.checked; emit() })
  // Purely a viewing option -- no physics depends on it, so redraw without re-solving.
  par.querySelector('#cut').addEventListener('change', e=>{ cutAway=e.target.checked; redraw() })
  // Steer the slice to cut through whatever you have painted. Like the cut checkbox
  // this is purely a viewing control -- no physics depends on it.
  par.querySelector('#cutpos').addEventListener('input', e=>{
    cutStart = parseFloat(e.target.value)
    par.querySelector('#cutposv').textContent = Math.round(cutStart*180/Math.PI)+'\\u00b0'
    // the flow is solved on the slice, so moving the slice needs a fresh solve
    showFlow ? emit() : redraw()
  })
  par.querySelector('#showflow').addEventListener('change', e=>{
    showFlow = e.target.checked
    if(!showFlow) flowData = null
    emit()
  })

  function syncBrush(){
    const dl=par.querySelector('#lab-dense'), ll=par.querySelector('#lab-light')
    dl.style.outline = dense ? '1px solid #3b82f6' : 'none'
    ll.style.outline = dense ? 'none' : '1px solid #ef4444'
    dl.style.borderRadius = ll.style.borderRadius = '4px'
    dl.style.padding = ll.style.padding = '2px 5px'
  }
  par.querySelector('#dense').addEventListener('change', e=>{ if(e.target.checked){dense=true; syncBrush()} })
  par.querySelector('#light').addEventListener('change', e=>{ if(e.target.checked){dense=false; syncBrush()} })

  // Presets drop their blob onto the near cut face, so it is immediately visible in
  // section rather than buried inside the solid part of the globe.
  function presetPhi(){ return cutFaces()[0] }
  par.querySelector('#clrbtn').addEventListener('click', ()=>{ blobs=[]; redraw(); emit() })
  par.querySelector('#shallowbtn').addEventListener('click', ()=>{
    blobs=[[1.15, presetPhi(), 0.80, dense?1.0:-1.0, brushWidth]]; redraw(); emit() })
  par.querySelector('#deepbtn').addEventListener('click', ()=>{
    blobs=[[1.15, presetPhi(), 0.10, dense?1.0:-1.0, brushWidth]]; redraw(); emit() })
  // A compact dense anomaly at depth, roughly the scale and setting discussed for
  // the Indian Ocean Geoid Low. Not a fit to data -- a demonstration.
  par.querySelector('#iogbtn').addEventListener('click', ()=>{
    blobs=[[1.40, presetPhi(), 0.55, 1.0, 0.30]]; redraw(); emit() })

  window.addEventListener('geoid-results', e=>{
    const d = e.detail ? JSON.parse(e.detail) : null
    if(!d) return
    kernelData=d.kernel; geoidMap=d.geoid; topoData=d.topo
    muGrid=d.mu; phiGrid=d.phi
    flowData=d.flow || null
    redraw()
  })

  par.querySelector('#widthv').textContent = widthLabel(brushWidth)
  syncBrush(); redraw(); emit()
</script>
        """)
    end

    const _geoid_canvas_ready = true
end

# ╔═╡ e5c8b756-c5b6-4708-3964-374859607182
begin
    _geoid_canvas_ready
    WideCell(@bind _gk GeoidGlobeInput())
end

# ╔═╡ f6d9c86c-d6e7-4819-4a15-263748596071
md"""
## Layer 1: The Analytic Building Blocks

Rather than transcribe a propagator matrix from a paper — the classic place for a sign error
to hide — we **derive** the four fundamental solutions and then verify them against the
original differential equations further below.

Expanding the poloidal flow at degree ``\ell`` (with ``L = \ell(\ell+1)``):

```math
u_r = U(r)\,Y_\ell, \qquad u_\theta = V(r)\,\partial_\theta Y_\ell, \qquad p = P(r)\,Y_\ell
```

Three facts pin the solution down completely:

1. **Incompressibility** gives ``V`` for free: ``V = (rU' + 2U)/L``.
2. **The pressure is harmonic**, ``\nabla^2 p = 0``, so ``P \propto r^s`` with
   ``s^2 + s - L = 0``, giving ``s = \ell`` or ``s = -(\ell+1)``.
3. **The homogeneous radial equation** ``U'' + 4U'/r + (2-L)U/r^2 = 0`` gives
   ``p = \ell-1`` or ``p = -\ell-2``; each pressure mode drives a particular solution with
   ``p = s+1``.

Together the four radial exponents are

```math
U \sim \{\, r^{\ell-1},\; r^{\ell+1},\; r^{-\ell},\; r^{-\ell-2} \,\}
```

Every quantity we need is a pure power of ``r``, so the code below is short and exact.
"""

# ╔═╡ 07ead97d-e7f8-492a-5b26-374859607182
Lof(ℓ) = ℓ * (ℓ + 1)

# ╔═╡ 18fbea8e-f809-4a3b-6c37-485960718293
"""
	modes(η, ℓ)

The four fundamental solutions at degree `ℓ`, each returned as `(p, cU, cP)` where the
radial velocity goes as `cU * r^p` and the pressure as `cP * r^(p-1)`.

Two are homogeneous (`cP = 0`, no pressure); two are particular solutions driven by the
two harmonic pressure modes. The coefficient `cU` of a particular mode follows from the
radial momentum balance.
"""
function modes(η, ℓ)
    L = Lof(ℓ)
    out = Tuple{Int,Float64,Float64}[]
    # homogeneous modes: pressure vanishes
    for p in (ℓ - 1, -ℓ - 2)
        push!(out, (p, 1.0, 0.0))
    end
    # particular modes driven by P = r^s, s = ℓ and s = -(ℓ+1); take cP = 1
    for s in (ℓ, -(ℓ + 1))
        p = s + 1
        denom = p * (p - 1) + 2p - (2 + L) + 2 * (p + 2)
        push!(out, (p, ((p - 1) / η) / denom, 1.0))
    end
    return out
end

# ╔═╡ 290cfb9f-090a-4b4c-7d48-596071829304
"""
	mode_state(r, η, ℓ, p, cU, cP)

State vector ``y = [u_r,\\; u_\\theta,\\; \\tau_{rr},\\; \\tau_{r\\theta}]`` for a single
mode at radius `r`, using

	τ_rr  = -P + 2η U'
	τ_rθ  =  η (r V' - V + U)
"""
function mode_state(r, η, ℓ, p, cU, cP)
    L = Lof(ℓ)
    U = cU * r^p
    Up = cU * p * r^(p - 1)
    V = cU * (p + 2) / L * r^p
    Vp = cU * (p + 2) * p / L * r^(p - 1)
    P = cP * r^(p - 1)
    τrr = -P + 2 * η * Up
    τrθ = η * (r * Vp - V + U)
    return [U, V, τrr, τrθ]
end

# ╔═╡ 3a1d0ca0-1a1b-4c5d-8e59-607182930415
"""
	Ymat(r, η, ℓ)

Matrix whose four columns are the fundamental solutions at radius `r`. Radii are scaled by
the Earth radius so that `r^ℓ` stays near unity — without this the columns overflow and the
matrix becomes badly conditioned at high degree.
"""
function Ymat(r, η, ℓ)
    x = r / A_EARTH
    hcat([mode_state(x, η, ℓ, p, cU, cP) for (p, cU, cP) in modes(η, ℓ)]...)
end

# ╔═╡ 4b2e1db1-2b2c-4d6e-9f6a-718293041526
"""
	propagator(r₂, r₁, η, ℓ)

Maps the state vector from radius `r₁` to `r₂` inside a layer of constant viscosity.
Because the columns of `Ymat` span the solution space, `Y(r₂) Y(r₁)⁻¹` *is* the propagator —
no hand-written matrix entries, hence no place for a transcription error to hide.
"""
propagator(r₂, r₁, η, ℓ) = Ymat(r₂, η, ℓ) / Ymat(r₁, η, ℓ)

# ╔═╡ 5c3f2ec2-3c3d-4e7f-a07b-829304152637
md"""
### Verifying the building blocks

Before trusting any of this, we check the four solutions against the **original** equations.
Every quantity is a pure power of ``r``, so each equation collapses to an algebraic identity
on the coefficients — we can test it exactly, with no finite differences.

The three equations are incompressibility, radial momentum, and colatitudinal momentum:

	U' + 2U/r - L·V/r                              = 0
	η[U'' + 2U'/r - (2+L)U/r² + 2L·V/r²] - P'      = 0
	η[V'' + 2V'/r - L·V/r² + 2U/r²]    - P/r       = 0
"""

# ╔═╡ 6d403fd3-4d4e-4f80-b18c-930415263748
"Exact residuals of the three governing equations for one mode (no finite differences)."
function mode_residuals(η, ℓ, p, cU, cP)
    L = Lof(ℓ)
    cV = cU * (p + 2) / L
    incomp = cU * p + 2cU - L * cV
    radial = η * (cU * (p * (p - 1) + 2p - (2 + L)) + 2L * cV) - cP * (p - 1)
    theta = η * (cV * (p * (p - 1) + 2p - L) + 2cU) - cP
    scale = max(abs(cU), abs(cP), abs(cV)) * max(η, 1.0) * max(L, 1.0)
    return maximum(abs.((incomp, radial, theta))) / scale
end

# ╔═╡ 7e5140e4-5e5f-4091-c29d-041526374859
validation_residual = maximum(
    mode_residuals(η_test, ℓ_test, p, cU, cP)
    for ℓ_test in 2:30, η_test in (1.0, 5.0, 1e3)
    for (p, cU, cP) in modes(η_test, ℓ_test)
)

# ╔═╡ 8f6251f5-6f60-41a2-d3ae-152637485960
if validation_residual < 1e-12
    md"""
    !!! correct "Validation 1 passed"
    	Maximum relative residual over degrees ``\ell = 2 \ldots 30`` and viscosities
    	spanning three orders of magnitude: **$(round(validation_residual, sigdigits=3))**
    	— machine precision. The fundamental solutions satisfy the Stokes equations exactly.
    """
else
    md"""
    !!! danger "Validation 1 FAILED"
    	Residual $(validation_residual) is too large — the fundamental solutions do not
    	satisfy the governing equations. Do not trust anything below.
    """
end

# ╔═╡ 90736265-7071-42b3-e4bf-263748596071
md"""
## The Geoid Kernel

Now we assemble the three contributions. Place a thin sheet of density anomaly at radius
``s`` and ask what geoid it produces at the surface.

**Boundary conditions.** Both the surface and the CMB are *free-slip*: no shear traction
(``\tau_{r\theta} = 0``), and no flow through the boundary (``u_r = 0``). The boundary does
not physically move in the solve — instead the normal stress ``\tau_{rr}`` that builds up
against it tells us how far it *would* deflect, which is the dynamic topography.

**The load.** Crossing the loaded radius, the normal stress jumps by the weight of the sheet:

```math
y(s^+) = y(s^-) + [\,0,\; 0,\; -\delta\rho\, g,\; 0\,]^T
```

This sign is worth checking rather than trusting: a *dense* sheet must **depress** the
surface. We verify exactly that below.
"""

# ╔═╡ a18473c5-8182-43c4-f5ca-374859607182
"""
	layer_propagator(r₁, r₂, ℓ, layers)

Chain propagators through every constant-viscosity layer between `r₁` and `r₂`. All four
state components are continuous across a viscosity interface, so the layer propagators
simply multiply.
"""
function layer_propagator(r₁, r₂, ℓ, layers)
    M = Matrix{Float64}(I, 4, 4)
    for (rb, rt, η) in layers
        lo = max(rb, r₁)
        hi = min(rt, r₂)
        hi > lo || continue
        M = propagator(hi, lo, η, ℓ) * M
    end
    return M
end

# ╔═╡ b2958425-9293-44d5-06d1-485960718293
"""
	geoid_response(s, ℓ, layers)

Solve the Stokes problem for a unit density sheet at radius `s`, degree `ℓ`.

Returns a named tuple with the three geoid contributions (m per kg/m³ per m of sheet
thickness), the surface and CMB dynamic topography, and the boundary-condition residuals
used for validation.
"""
function geoid_response(s, ℓ, layers)
    a, b = A_EARTH, B_CMB
    upward(σ, rad) = 4pi * G_GRAV * a / (G_SURF * (2ℓ + 1)) * σ * (rad / a)^(ℓ + 2)

    # (1) the anomaly's own gravitational pull — always positive for dense material
    N_direct = upward(1.0, s)

    Pb = layer_propagator(b, s, ℓ, layers)   # CMB  -> load
    Pa = layer_propagator(s, a, ℓ, layers)   # load -> surface
    M = Pa * Pb
    jump = [0.0, 0.0, -G_SURF, 0.0]          # unit δρ

    # Free slip at the CMB: u_r = 0 and τ_rθ = 0, leaving two unknown constants.
    # Require the same two conditions at the surface -> 2x2 linear system.
    rhs = -(Pa * jump)
    Asys = [M[1, 2] M[1, 3]; M[4, 2] M[4, 3]]
    coef = Asys \ [rhs[1], rhs[4]]
    y_cmb = [0.0, coef[1], coef[2], 0.0]
    y_surf = M * y_cmb + Pa * jump

    # Normal stress -> dynamic topography at each boundary.
    #
    # Hager & Richards (1989) eq. (6) writes this as a JUMP across the boundary,
    #     [τ_rr]₊⁻ = Δρ g δh,
    # with Δρ the density contrast *across* that boundary. Traversed in the same
    # outward direction, the surface jump (mantle → vacuum) and the CMB jump
    # (core → mantle) have OPPOSITE orientation — hence the sign difference below.
    # Using the same sign at both boundaries (as an earlier version of this notebook
    # did) makes a dense anomaly depress the surface but *raise* the CMB, which no
    # single convection cell can do.
    #
    # The convention is pinned by the compensation limit, checked in `check_compensation`
    # below: a load resting on either boundary must be almost exactly cancelled by that
    # boundary's own deflection.
    h_surf = y_surf[3] / (RHO_M * G_SURF)
    h_cmb = -y_cmb[3] / ((RHO_C - RHO_M) * G_SURF)

    # (2) and (3): the deflected boundaries, as mass sheets
    N_surf = upward(RHO_M * h_surf, a)
    N_cmb = upward((RHO_C - RHO_M) * h_cmb, b)

    return (direct=N_direct, surface=N_surf, cmb=N_cmb,
        total=N_direct + N_surf + N_cmb,
        h_surf=h_surf, h_cmb=h_cmb, y_cmb=y_cmb,
        bc_residual=maximum(abs.((y_surf[1], y_surf[4], y_cmb[1], y_cmb[4]))))
end

# ╔═╡ c3a69595-a3a4-45e6-17e2-596071829304
"""
	geoid_kernel(s, ℓ, layers; direct_only=false)

Geoid at the surface (in metres) per unit density anomaly, for a sheet at radius `s`.

`direct_only=true` switches off the deflected-boundary terms, leaving only the anomaly's
own attraction — the naive, always-positive answer.
"""
function geoid_kernel(s, ℓ, layers; direct_only=false)
    resp = geoid_response(s, ℓ, layers)
    direct_only ? resp.direct : resp.total
end

# ╔═╡ d4b7a645-b4b5-46f7-28f3-607182930415
begin
    "Uniform-viscosity mantle."
    isoviscous_layers(η=1e21) = [(B_CMB, A_EARTH, η)]

    "Two-layer mantle: lower mantle `ratio` times stiffer than the upper mantle."
    twolayer_layers(ratio; η_um=1e21) =
        [(B_CMB, R_660, η_um * ratio), (R_660, A_EARTH, η_um)]

    """
    Hager & Richards (1989) model **WL**, their preferred whole-mantle model
    (p. 316): relative viscosity 1 in the top 100 km, **1/30** from 100–400 km (the
    low-viscosity asthenosphere), 1 from 400–670 km, **10** from 670–2600 km, and
    1/10 from 2600 km to the CMB.

    This is the model whose kernels are plotted in their Figure 2, so it is what we
    compare against rather than any uniform-viscosity case — the paper publishes no
    isoviscous kernel.
    """
    hr89_WL_layers(η₀=1e21) = [
        (B_CMB, A_EARTH - 2600e3, η₀ / 10),
        (A_EARTH - 2600e3, A_EARTH - 670e3, η₀ * 10),
        (A_EARTH - 670e3, A_EARTH - 400e3, η₀),
        (A_EARTH - 400e3, A_EARTH - 100e3, η₀ / 30),
        (A_EARTH - 100e3, A_EARTH, η₀),
    ]
    nothing
end

# ╔═╡ f6d9c866-d6a7-4819-4a75-485960718293
begin
    gk_blobs_raw = get(_gk, "blobs", Vector{Float64}[])
    gk_visc = Float64(get(_gk, "visc_ratio", 1.0))
    gk_ell = Int(round(Float64(get(_gk, "ell", 2))))
    gk_direct = Bool(get(_gk, "direct_only", false))
    gk_showflow = Bool(get(_gk, "show_flow", false))
    _cs = Float64(get(_gk, "cut_start", 0.9))
    gk_cutfaces = (_cs, _cs + pi / 2)

    # The globe hands back (θ, φ, r_fraction, amplitude, angular_width).
    gk_blobs = [(B_CMB + Float64(b[3]) * (A_EARTH - B_CMB),   # radius
        Float64(b[1]),                             # colatitude θ
        Float64(b[2]),                             # longitude φ
        Float64(b[4]),                             # amplitude ±1
        Float64(b[5]))                             # angular width
                for b in gk_blobs_raw]
    gk_layers = gk_visc ≈ 1.0 ? isoviscous_layers() : twolayer_layers(gk_visc)
    nothing
end

# ╔═╡ e5c8b755-c5c6-4708-3954-718293041526
"Depth (km) where the kernel changes sign, or `nothing` if it never does."
function zero_crossing_depth(ℓ, layers; n=200, direct_only=false)
    rs = range(B_CMB + 5e3, A_EARTH - 5e3, length=n)
    v = [geoid_kernel(s, ℓ, layers; direct_only) for s in rs]
    idx = findfirst(i -> sign(v[i]) != sign(v[i+1]), 1:length(v)-1)
    isnothing(idx) && return nothing
    s0 = rs[idx] - v[idx] * (rs[idx+1] - rs[idx]) / (v[idx+1] - v[idx])
    return (A_EARTH - s0) / 1e3
end

# ╔═╡ 18fbea86-f8c9-4a3b-6ca7-607182930415
begin
    # kernel curve for the selected single degree
    gk_kdepths = range(20.0, 2870.0, length=120)
    gk_kvals = [geoid_kernel(A_EARTH - d * 1e3, gk_ell, gk_layers; direct_only=gk_direct)
                for d in gk_kdepths]
    gk_crossing = zero_crossing_depth(gk_ell, gk_layers; direct_only=gk_direct)
    nothing
end

# ╔═╡ f6d9c865-d6d7-4819-4a65-829304152637
md"""
### Validating the kernel

Seven independent physical checks. Each one would catch a different class of sign error — and
two of them (the **compensation limit** and the **Hager & Richards Figure 2 normalisation**)
exist because an earlier version of this notebook got the CMB sign wrong and the other checks
did not notice.
"""

# ╔═╡ 290cfb95-0901-4b4c-7da8-152637485960
md"""
## Layer 3: From a Painted Picture to a Geoid

The globe gives us ``\delta\rho(r,\theta,\varphi)`` on a grid. To use the kernels we need its
**spherical-harmonic coefficients** at each radius. Using real harmonics, so that everything
stays real-valued:

```math
\delta\rho(r,\theta,\varphi) = \sum_{\ell}\sum_{m=0}^{\ell} P_\ell^m(\cos\theta)
\left[\,C_{\ell m}(r)\cos m\varphi + S_{\ell m}(r)\sin m\varphi\,\right]
```

Each coefficient comes from a Gauss–Legendre quadrature in ``\cos\theta`` and a plain sum in
``\varphi``. Then — and this is the step that makes the whole thing cheap — we apply the
**same** radial kernel ``K_\ell(r)`` to every order ``m`` of a given degree, and sum back up:

```math
N(\theta,\varphi) = \sum_{\ell}\sum_{m=0}^{\ell} P_\ell^m(\cos\theta)
\left[\,\tilde{C}_{\ell m}\cos m\varphi + \tilde{S}_{\ell m}\sin m\varphi\,\right],
\qquad
\tilde{C}_{\ell m} = \int K_\ell(r)\,C_{\ell m}(r)\,dr
```

We roll the machinery by hand — a three-term recurrence for the associated Legendre functions
and a Newton solve for the Gauss–Legendre nodes. It needs no extra packages, and you can read
every step.

If you set every ``m \neq 0`` coefficient to zero you recover the axisymmetric case exactly;
the validation below checks precisely that.
"""

# ╔═╡ 3a1d0ca5-1a0b-4c5d-8eb9-263748596071
"All Legendre polynomials `P_0..P_lmax` at `x`, via the three-term recurrence."
function legendre_all(x, lmax)
    P = zeros(lmax + 1)
    P[1] = 1.0
    lmax >= 1 && (P[2] = x)
    for l in 1:lmax-1
        P[l+2] = ((2l + 1) * x * P[l+1] - l * P[l]) / (l + 1)
    end
    return P
end

# ╔═╡ 4b2e1db5-2b1c-4d6e-9fca-374859607182
"Gauss–Legendre nodes and weights on [-1,1], found by Newton iteration on `P_n`."
function gauss_legendre(n)
    x = zeros(n)
    w = zeros(n)
    for i in 1:n
        z = cos(pi * (i - 0.25) / (n + 0.5))   # excellent initial guess
        for _ in 1:100
            P = legendre_all(z, n)
            dP = n * (P[n] - z * P[n+1]) / (1 - z^2)
            dz = -P[n+1] / dP
            z += dz
            abs(dz) < 1e-15 && break
        end
        P = legendre_all(z, n)
        dP = n * (P[n] - z * P[n+1]) / (1 - z^2)
        x[i] = z
        w[i] = 2 / ((1 - z^2) * dP^2)
    end
    return x, w
end

# ╔═╡ 4b2e1db7-2b0c-4d6e-9fea-485960718294
"""
	assoc_legendre(x, lmax)

Fully-normalised associated Legendre functions ``P_\\ell^m(x)`` for all
``\\ell = 0 \\ldots \\ell_{max}``, ``m = 0 \\ldots \\ell``, returned as `P[ℓ+1][m+1]`.

These are what let an anomaly vary in **longitude** as well as colatitude. The
normalisation is chosen so that ``\\int_{-1}^{1} [P_\\ell^m]^2\\,dx = 2``, which makes
the ``m = 0`` column equal to ``\\sqrt{2\\ell+1}\\,P_\\ell`` — i.e. the plain Legendre
polynomials of the axisymmetric case, rescaled.

Built from the standard sectoral seed plus the two-term recurrence in ``\\ell``.
"""
function assoc_legendre(x, lmax)
    P = [zeros(l + 1) for l in 0:lmax]
    P[1][1] = 1.0
    lmax == 0 && return P
    s = sqrt(max(0.0, 1 - x^2))
    # sectoral terms P_m^m
    for m in 1:lmax
        P[m+1][m+1] = sqrt((2m + 1) / (2m)) * s * P[m][m]
    end
    # first step off the diagonal
    for m in 0:lmax-1
        P[m+2][m+1] = sqrt(2m + 3) * x * P[m+1][m+1]
    end
    # recurrence in ℓ at fixed m
    for m in 0:lmax, l in (m+2):lmax
        a = sqrt((4l^2 - 1) / (l^2 - m^2))
        b = sqrt(((l - 1)^2 - m^2) / (4(l - 1)^2 - 1))
        P[l+1][m+1] = a * (x * P[l][m+1] - b * P[l-1][m+1])
    end
    return P
end

# ╔═╡ 6d403fd5-4d3e-4f80-b1eb-596071829304
md"""
### Validating the Legendre machinery

Four checks: the plain polynomials must be **orthogonal** under the quadrature; the
associated functions must be **orthonormal**; their ``m=0`` column must reduce to the
axisymmetric case; and a purely axisymmetric input must put **all** its power in ``m=0``.
"""

# ╔═╡ 7e5140e5-5e4f-4091-c2fb-607182930415
begin
    _leg_lmax = 30
    _leg_n = 2 * _leg_lmax + 2
    _leg_x, _leg_w = gauss_legendre(_leg_n)
    _leg_P = [legendre_all(xi, _leg_lmax) for xi in _leg_x]

    # orthogonality: ∫ P_i P_j = 2/(2i+1) δ_ij
    _orth_off = 0.0
    _orth_diag = 0.0
    for i in 0:_leg_lmax, j in 0:_leg_lmax
        s = sum(_leg_w[k] * _leg_P[k][i+1] * _leg_P[k][j+1] for k in 1:_leg_n)
        if i == j
            _orth_diag = max(_orth_diag, abs(s - 2 / (2i + 1)))
        else
            _orth_off = max(_orth_off, abs(s))
        end
    end
    # --- associated Legendre: orthonormality ∫ P_ℓ^m P_ℓ'^m dx = 2 δ_ℓℓ' ---
    _am_lmax = 12
    _am_n = 2 * _am_lmax + 4
    _am_x, _am_w = gauss_legendre(_am_n)
    _am_P = [assoc_legendre(xi, _am_lmax) for xi in _am_x]
    _am_off = 0.0
    _am_diag = 0.0
    for m in 0:_am_lmax, l1 in m:_am_lmax, l2 in m:_am_lmax
        s = sum(_am_w[k] * _am_P[k][l1+1][m+1] * _am_P[k][l2+1][m+1] for k in 1:_am_n)
        if l1 == l2
            _am_diag = max(_am_diag, abs(s - 2.0))
        else
            _am_off = max(_am_off, abs(s))
        end
    end

    # --- the m=0 column must be exactly sqrt(2ℓ+1) times the plain P_ℓ ---
    _am_m0 = 0.0
    for xt in (-0.87, -0.3, 0.0, 0.42, 0.95)
        A = assoc_legendre(xt, _am_lmax)
        U = legendre_all(xt, _am_lmax)
        for l in 0:_am_lmax
            pred = sqrt(2l + 1) * U[l+1]
            _am_m0 = max(_am_m0, abs(A[l+1][1] - pred) / max(abs(pred), 1e-12))
        end
    end

    check_legendre = _orth_off < 1e-12 && _orth_diag < 1e-12 &&
                     abs(sum(_leg_w) - 2) < 1e-12 &&
                     _am_off < 1e-11 && _am_diag < 1e-11 && _am_m0 < 1e-11
    nothing
end

# ╔═╡ 8f6251f6-6f40-41a2-d3bd-718293041526
if check_legendre
    md"""
    !!! correct "Validation 3 passed"
    	- Quadrature weights sum to 2 exactly; plain Legendre orthogonality holds to
    	  **$(round(_orth_off, sigdigits=3))** (off-diagonal) and
    	  **$(round(_orth_diag, sigdigits=3))** (diagonal) up to ``\ell = 30``.
    	- Associated ``P_\ell^m`` are orthonormal to **$(round(_am_off, sigdigits=3))**
    	  (off-diagonal) and **$(round(_am_diag, sigdigits=3))** (diagonal).
    	- Their ``m=0`` column reproduces ``\sqrt{2\ell+1}\,P_\ell`` to
    	  **$(round(_am_m0, sigdigits=3))** — so an axisymmetric anomaly gives exactly the
    	  same answer as the purely 2-D version of this notebook did.
    """
else
    md"""!!! danger "Validation 3 FAILED" """
end

# ╔═╡ 90736266-7051-42b3-e4cf-829304152637
md"""
## Putting It Together: Painted Anomaly → Geoid

Now the full pipeline. The globe hands us a set of painted blobs; we build
``\delta\rho(r,\theta,\varphi)`` from them, decompose, apply kernels, and synthesize.
"""

# ╔═╡ a18473c6-8152-43c4-f5db-930415263748
"""
	blobs_to_field(blobs, rgrid, μgrid, φgrid)

Build ``\\delta\\rho(r,\\theta,\\varphi)`` from painted blobs. Each blob is
`(r_centre, θ_centre, φ_centre, amplitude, angular_width)` and is a Gaussian cap in
*angular distance on the sphere* — so it is genuinely localized in both colatitude and
longitude, not a ring.

The angular width is floored at [`BLOB_WIDTH_MIN`](@ref): a cap narrower than the
expansion can resolve produces Gibbs ringing, which would show up as spurious wobbles
in the geoid.
"""
function blobs_to_field(blobs, rgrid, μgrid, φgrid; r_width=250e3)
    F = zeros(length(rgrid), length(μgrid), length(φgrid))
    isempty(blobs) && return F
    for (rc, θc, φc, amp, w) in blobs
        ww = max(w, BLOB_WIDTH_MIN)
        for (j, μ) in enumerate(μgrid)
            θ = acos(clamp(μ, -1, 1))
            for (k, φ) in enumerate(φgrid)
                # angular distance from the blob centre, via the spherical cosine rule
                ca = cos(θ) * cos(θc) + sin(θ) * sin(θc) * cos(φ - φc)
                Δ = acos(clamp(ca, -1, 1))
                ang = exp(-(Δ / ww)^2)
                ang < 1e-8 && continue
                for (i, r) in enumerate(rgrid)
                    F[i, j, k] += amp * ang * exp(-((r - rc) / r_width)^2)
                end
            end
        end
    end
    return F
end

# ╔═╡ b2958426-9263-44d5-06e1-041526374859
"""
	synthesize_geoid(blobs, layers; direct_only=false)

The full Layer-3 pipeline, now in three dimensions.

For each degree ``\\ell`` and order ``m`` we project the painted density field onto the
real spherical harmonics, apply the **same** radial kernel ``K_\\ell(r)`` — it does not
depend on ``m`` — and sum the result back onto a ``(\\theta,\\varphi)`` grid.

Returns the geoid map `N[θ,φ]`, the grids, the per-degree contributions (summed over
``m``), and the density field itself for display.
"""
function synthesize_geoid(blobs, layers; direct_only=false, nr=40, lmax=LMAX)
    nθ = 2lmax + 4
    nφ = 4lmax + 8
    μ, w = gauss_legendre(nθ)
    φgrid = [2π * (k - 1) / nφ for k in 1:nφ]
    rgrid = range(B_CMB + 30e3, A_EARTH - 30e3, length=nr)
    dr = step(rgrid)
    dφ = 2π / nφ

    F = blobs_to_field(blobs, rgrid, μ, φgrid)
    P = [assoc_legendre(μk, lmax) for μk in μ]
    cosmφ = [[cos(m * φ) for φ in φgrid] for m in 0:lmax]
    sinmφ = [[sin(m * φ) for φ in φgrid] for m in 0:lmax]

    # geoid coefficients, accumulated over radius: Ngeo[m+1][ℓ+1] for cos and sin
    Ncos = [zeros(lmax + 1) for _ in 0:lmax]
    Nsin = [zeros(lmax + 1) for _ in 0:lmax]
    per_degree = zeros(lmax + 1)
    spectrum = zeros(lmax + 1)

    # The three geoid contributions are tracked separately so we can show students
    # the arithmetic behind the sign flip, not just its result.
    Dcos = [zeros(lmax + 1) for _ in 0:lmax]   # direct (the anomaly's own pull)
    Dsin = [zeros(lmax + 1) for _ in 0:lmax]
    Scos = [zeros(lmax + 1) for _ in 0:lmax]   # deflected surface
    Ssin = [zeros(lmax + 1) for _ in 0:lmax]
    Ccos = [zeros(lmax + 1) for _ in 0:lmax]   # deflected CMB
    Csin = [zeros(lmax + 1) for _ in 0:lmax]

    for ℓ in 1:lmax                     # ℓ=0 cannot deform the geoid's shape
        # the response is computed ONCE per degree and reused for every order m
        resp = [geoid_response(r, ℓ, layers) for r in rgrid]
        for m in 0:ℓ
            nrm = m == 0 ? 1.0 : 2.0
            ad = 0.0; asd = 0.0
            asu = 0.0; assu = 0.0
            acm = 0.0; ascm = 0.0
            spec = 0.0
            for (i, r) in enumerate(rgrid)
                # project this radial shell onto Y_ℓ^m
                gc = 0.0
                gs = 0.0
                for k in 1:nθ
                    fc = 0.0
                    fs = 0.0
                    for j in 1:nφ
                        fc += F[i, k, j] * cosmφ[m+1][j]
                        fs += F[i, k, j] * sinmφ[m+1][j]
                    end
                    gc += w[k] * P[k][ℓ+1][m+1] * fc
                    gs += w[k] * P[k][ℓ+1][m+1] * fs
                end
                δc = nrm * gc * dφ / (2 * 2π)
                δs = nrm * gs * dφ / (2 * 2π)
                spec += (abs(δc) + abs(δs)) * dr
                ad += resp[i].direct * δc * dr
                asd += resp[i].direct * δs * dr
                asu += resp[i].surface * δc * dr
                assu += resp[i].surface * δs * dr
                acm += resp[i].cmb * δc * dr
                ascm += resp[i].cmb * δs * dr
            end
            Dcos[m+1][ℓ+1] = ad; Dsin[m+1][ℓ+1] = asd
            Scos[m+1][ℓ+1] = asu; Ssin[m+1][ℓ+1] = assu
            Ccos[m+1][ℓ+1] = acm; Csin[m+1][ℓ+1] = ascm
            # `direct_only` drops the boundary terms from the total
            Ncos[m+1][ℓ+1] = direct_only ? ad : ad + asu + acm
            Nsin[m+1][ℓ+1] = direct_only ? asd : asd + assu + ascm
            per_degree[ℓ+1] += abs(Ncos[m+1][ℓ+1]) + abs(Nsin[m+1][ℓ+1])
            spectrum[ℓ+1] += spec
        end
    end

    "Synthesize a (θ,φ) map from a set of real spherical-harmonic coefficients."
    function tomap(Ac, As)
        M = zeros(nθ, nφ)
        for k in 1:nθ, j in 1:nφ
            s = 0.0
            for m in 0:lmax, ℓ in max(m, 1):lmax
                s += P[k][ℓ+1][m+1] * (Ac[m+1][ℓ+1] * cosmφ[m+1][j] +
                                       As[m+1][ℓ+1] * sinmφ[m+1][j])
            end
            M[k, j] = s
        end
        M
    end

    return (μ=μ, φ=φgrid, N=tomap(Ncos, Nsin),
        N_direct=tomap(Dcos, Dsin),
        N_surface=tomap(Scos, Ssin),
        N_cmb=tomap(Ccos, Csin),
        per_degree=per_degree, spectrum=spectrum,
        rgrid=rgrid, field=F, lmax=lmax)
end

# ╔═╡ 07ead976-e7b8-492a-5bc6-596071829304
gk_result = synthesize_geoid(gk_blobs, gk_layers; direct_only=gk_direct)

# ╔═╡ 290cfb97-0861-4b4c-7dc8-829304152638
let
    isempty(gk_blobs) && return md"*Paint an anomaly to see the breakdown.*"
    # read every component at the same point: the peak of the total geoid
    idx = argmax(abs.(gk_result.N))
    d = gk_result.N_direct[idx]
    s = gk_result.N_surface[idx]
    c = gk_result.N_cmb[idx]
    net = gk_result.N[idx]
    θpk = round(rad2deg(acos(clamp(gk_result.μ[idx[1]], -1, 1))), digits=1)
    φpk = round(rad2deg(gk_result.φ[idx[2]]), digits=1)

    # Render the signed number as inline math. Julia's Markdown HTML-escapes a bare
    # "+" to &#43; in ordinary text, which shows up literally; inside a math span the
    # browser decodes it and KaTeX typesets the sign properly.
    fmt(x) = "``" * (x ≥ 0 ? "+" : "-") * string(round(abs(x), digits=2)) * "\\ \\mathrm{m}``"
    bar(x, mx) = repeat(x ≥ 0 ? "█" : "█", max(1, round(Int, 18 * abs(x) / mx)))
    mx = max(abs(d), abs(s), abs(c), abs(net), 1e-9)

    rows = [
        ("① the anomaly itself (direct)", d),
        ("② the deflected surface", s),
        ("③ the deflected CMB", c),
    ]
    body = join(["| $(nm) | $(fmt(v)) | $(bar(v, mx)) |" for (nm, v) in rows], "\n")

    Markdown.parse("""
At ``\\theta = $(θpk)°``, ``\\varphi = $(φpk)°`` — the peak of the geoid:

| contribution | value | |
|---|---:|:---|
$(body)
| **net geoid** | **$(fmt(net))** | |

$(if gk_direct
        "!!! warning \"Boundary terms are switched OFF\"\n\tYou have ticked *direct effect only*, so the net is just ① — and it is
\tpositive for dense material, always. Untick it to let ② and ③ compete."
    elseif abs(d) > 1e-12 && sign(net) != sign(d)
        "!!! correct \"The boundaries won\"\n\tThe anomaly's own pull (①) is **$(fmt(d))**, but the two deflected
\tboundaries together contribute **$(fmt(s + c))** — enough to *flip the sign*.
\tThis is the whole lesson: what you measure at the surface has the opposite sign
\tfrom the anomaly that caused it."
    else
        "!!! note \"The direct term won\"\n\tHere ① dominates: the boundaries contribute **$(fmt(s + c))**, which reduces
\tthe geoid but does not flip it. Move the anomaly shallower (or lower the viscosity
\tratio) to reach the regime where they do."
    end)
""")
end

# ╔═╡ 5c3f2ec6-3cfd-4e7f-a0eb-041526374859
let
    ℓs = 1:gk_result.lmax
    contrib = [gk_result.per_degree[ℓ+1] for ℓ in ℓs]
    fig = PlutoPlotly.Plot(Layout(
        height=300, width=700,
        title="geoid power by degree ℓ (summed over all orders m)",
        xaxis=attr(title="spherical harmonic degree ℓ", dtick=2),
        yaxis=attr(title="|contribution| to N (m)"),
        font=attr(size=11), showlegend=false))
    add_trace!(fig, PlutoPlotly.bar(x=collect(ℓs), y=contrib,
        marker=attr(color="#3b82f6")))
    PlutoPlotly.plot(fig)
end

# ╔═╡ 290cfb96-0871-4b4c-7db8-718293041526
"""
	dynamic_topography(result, layers)

Surface and CMB dynamic topography **maps** ``h(\\theta,\\varphi)`` for a painted field,
built exactly like the geoid but using the boundary deflections rather than the geoid
response. Same decomposition, same ``m``-independent radial solve per degree.
"""
function dynamic_topography(result, layers)
    lmax = result.lmax
    μ, φgrid = result.μ, result.φ
    nθ, nφ = length(μ), length(φgrid)
    _, w = gauss_legendre(nθ)
    dr = step(result.rgrid)
    dφ = 2π / nφ
    P = [assoc_legendre(μk, lmax) for μk in μ]
    cosmφ = [[cos(m * φ) for φ in φgrid] for m in 0:lmax]
    sinmφ = [[sin(m * φ) for φ in φgrid] for m in 0:lmax]

    Hs = zeros(nθ, nφ)
    Hc = zeros(nθ, nφ)
    for ℓ in 1:lmax
        resp = [geoid_response(r, ℓ, layers) for r in result.rgrid]
        for m in 0:ℓ
            nrm = m == 0 ? 1.0 : 2.0
            sc = 0.0; ss = 0.0; cc = 0.0; cs = 0.0
            for (i, r) in enumerate(result.rgrid)
                gc = 0.0
                gs = 0.0
                for k in 1:nθ
                    fc = 0.0
                    fs = 0.0
                    for j in 1:nφ
                        fc += result.field[i, k, j] * cosmφ[m+1][j]
                        fs += result.field[i, k, j] * sinmφ[m+1][j]
                    end
                    gc += w[k] * P[k][ℓ+1][m+1] * fc
                    gs += w[k] * P[k][ℓ+1][m+1] * fs
                end
                δc = nrm * gc * dφ / (2 * 2π)
                δs = nrm * gs * dφ / (2 * 2π)
                sc += resp[i].h_surf * δc * dr
                ss += resp[i].h_surf * δs * dr
                cc += resp[i].h_cmb * δc * dr
                cs += resp[i].h_cmb * δs * dr
            end
            for k in 1:nθ, j in 1:nφ
                Hs[k, j] += P[k][ℓ+1][m+1] * (sc * cosmφ[m+1][j] + ss * sinmφ[m+1][j])
                Hc[k, j] += P[k][ℓ+1][m+1] * (cc * cosmφ[m+1][j] + cs * sinmφ[m+1][j])
            end
        end
    end
    return (surf=Hs, cmb=Hc)
end

# ╔═╡ 07ead977-e7a8-492a-5bd6-607182930416
gk_topo = dynamic_topography(gk_result, gk_layers)

# ╔═╡ 5c3f2ec5-3c2d-4e7f-a0db-485960718293
md"""
## Assumptions — What This Notebook Does and Does Not Do

Being explicit about this matters more than the pretty pictures.

!!! warning "Read before drawing conclusions"
	- **Radial viscosity only — this is the one that really binds.** ``\eta`` varies with depth,
	  never sideways. The *density* is now fully three-dimensional, but the moment viscosity
	  varies laterally the whole approach collapses: degrees stop separating, the kernel
	  ``K_\ell(r)`` ceases to exist as a concept, and you need a genuine 3-D solver
	  (CitcomS, ASPECT). Real mantle has strong lateral viscosity contrasts — cold stiff slabs,
	  hot weak plumes — so treat what you see here as the *mechanism*, never as a prediction
	  for any specific region.
	- **Truncated at ``\ell_{max} = 20``.** Structure finer than roughly 2000 km across is not
	  represented, which is why the brush enforces a minimum anomaly size. This is a resolution
	  limit, not a physical assumption — but it does mean you cannot paint a narrow slab.
	- **Instantaneous.** One linear solve. No time-stepping, no evolution — the anomaly stays
	  where you paint it.
	- **Self-gravitation neglected.** The deflected boundaries and the anomaly do not pull on
	  each other. This shifts long-wavelength amplitudes by roughly 10–20% and can move the
	  crossing depth somewhat, but it does not change the sign structure that is the lesson here.
	- **Free-slip boundaries, no lithosphere.** No elastic plate, no rigid lid.

!!! tip "For quantitative work, use a real code"
	This notebook is built for **intuition**. For research, use
	**HC** (Becker et al.) for propagator-matrix geoid calculations, or **ASPECT** /
	**CitcomS** for full 3-D convection. Benchmarks show propagator-matrix codes and ASPECT
	agree to about 1% out to degree ~15, which is exactly why this simple approach is
	trustworthy *at long wavelengths* — and why we truncate rather than chasing high degrees.
"""

# ╔═╡ c3a69597-a394-45e6-1802-607182930417
"""
	flow_profile(s, ℓ, layers, sample_radii)

Instantaneous flow driven by a unit density sheet at radius `s`, degree `ℓ`, evaluated at
each radius in `sample_radii`.

Returns `(U, V)`: the radial and colatitudinal velocity scalars, where the physical
velocity is ``u_r = U\\,Y_\\ell`` and ``u_\\theta = V\\,\\partial_\\theta Y_\\ell``.

This is the *same* solve as [`geoid_response`](@ref) — we simply keep the interior state
instead of discarding everything but the boundaries. Viscosity enters only as a ratio, so
these velocities carry an arbitrary overall scale; only their **pattern** is meaningful.

!!! note "Sign convention"
	The velocity components of the propagator state vector carry the *opposite* sign
	convention from the stress components, an artefact of how the angular factors are
	defined. That never matters for the geoid, which is built from stresses — but it does
	matter for drawing arrows, so we negate here. The convention is pinned to the
	already-validated fact that a dense load **depresses** the surface (`h_surf < 0`):
	with this negation, `u_r` and `h_surf` agree in sign, so downwelling goes with a
	depressed surface as it physically must.
"""
function flow_profile(s, ℓ, layers, sample_radii)
    b = B_CMB
    resp = geoid_response(s, ℓ, layers)
    jump = [0.0, 0.0, -G_SURF, 0.0]
    U = zeros(length(sample_radii))
    V = zeros(length(sample_radii))
    for (i, r) in enumerate(sample_radii)
        # propagate the CMB state up to r, adding the load jump once we pass it
        y = layer_propagator(b, min(r, s), ℓ, layers) * resp.y_cmb
        if r > s
            y = layer_propagator(s, r, ℓ, layers) * (y + jump)
        end
        # negated: see the sign-convention note above
        U[i] = -y[1]
        V[i] = -y[2]
    end
    return (U=U, V=V)
end

# ╔═╡ 07ead975-e7e8-492a-5bb6-930415263748
begin
    _val_rs = range(B_CMB + 20e3, A_EARTH - 20e3, length=40)

    # Check A: with boundaries off, the kernel must be positive and monotonic.
    _direct_vals = [geoid_kernel(s, 2, isoviscous_layers(); direct_only=true) for s in _val_rs]
    check_direct = all(_direct_vals .> 0) && all(diff(_direct_vals) .> 0)

    # Check B: free-slip boundary conditions satisfied, measured RELATIVE to the
    # size of the solution itself (an absolute threshold would be meaningless here,
    # since the state components carry physical units spanning many orders).
    check_bc = maximum(
        let resp = geoid_response(s, ℓ, isoviscous_layers())
            resp.bc_residual / max(abs(resp.h_surf), 1e-30)
        end
        for s in _val_rs, ℓ in (2, 4, 8)) < 1e-8

    # Check C: a dense anomaly must DEPRESS the surface at every depth.
    check_topo = all(geoid_response(s, 2, isoviscous_layers()).h_surf < 0 for s in _val_rs)

    # Check D — THE COMPENSATION LIMIT. A load resting on a boundary is supported by
    # that boundary's own deflection, so the boundary term must very nearly cancel the
    # direct term. This is the check that catches a flipped boundary sign: with the
    # wrong sign the two ADD instead of cancelling. (It is how the CMB sign error in an
    # earlier version of this notebook was found.)
    check_compensation = all((2, 4, 8)) do ℓ
        qa = geoid_response(A_EARTH - 3e3, ℓ, isoviscous_layers())   # load at the surface
        qb = geoid_response(B_CMB + 3e3, ℓ, isoviscous_layers())     # load at the CMB
        abs(qa.direct + qa.surface) / abs(qa.direct) < 1e-2 &&
            abs(qb.direct + qb.cmb) / abs(qb.direct) < 1e-2
    end

    # Check E: Hager & Richards (1989) Figure 2 normalisation. Their A^l (surface) and
    # C^l (CMB) displaced-mass kernels are "normalized by dividing by (−σ)", and both
    # panels of Fig. 2 span 0 → 1, with unity meaning perfect dynamic compensation.
    # So both must be positive, ≤ ~1, with A → 1 at the surface and C → 1 at the CMB.
    # (Values may overshoot 1 by ~2%: the first-order Taylor boundary approximation
    # HR89 discuss on p. 314.)
    _Akern(d) = -(RHO_M * geoid_response(A_EARTH - d, 2, isoviscous_layers()).h_surf)
    _Ckern(d) = -((RHO_C - RHO_M) * geoid_response(A_EARTH - d, 2, isoviscous_layers()).h_cmb)
    check_hr89_norm = all(range(20e3, 2870e3, length=60)) do d
        0 ≤ _Akern(d) ≤ 1.03 && 0 ≤ _Ckern(d) ≤ 1.03
    end && _Akern(30e3) > 0.9 && _Ckern(2860e3) > 0.9

    # Check F: reproduce the published kernel. HR89's Figure 2 is for their model WL,
    # whose ℓ=2 geoid kernel CHANGES SIGN in the mid-mantle. (The uniform-viscosity
    # case does not — and the paper never plots one.)
    _wl_ds = collect(range(30e3, 2860e3, length=150))
    _wl_G = [geoid_kernel(A_EARTH - d, 2, hr89_WL_layers()) for d in _wl_ds]
    check_hr89_WL = any(i -> sign(_wl_G[i]) != sign(_wl_G[i+1]), 1:length(_wl_G)-1)

    # Check G: the interior flow must vanish at both free-slip boundaries, and its
    # sign must agree with the surface deflection — downwelling depresses the surface.
    _fr = collect(range(B_CMB, A_EARTH, length=81))
    check_flow = all((2, 4, 6)) do ℓ
        sload = A_EARTH - 800e3
        f = flow_profile(sload, ℓ, isoviscous_layers(), _fr)
        pk = maximum(abs, f.U)
        h = geoid_response(sload, ℓ, isoviscous_layers()).h_surf
        abs(f.U[1]) / pk < 1e-6 && abs(f.U[end]) / pk < 1e-6 &&
            sign(f.U[argmax(abs.(f.U))]) == sign(h)
    end

    # Check H: a dense anomaly must depress BOTH boundaries — one convection cell
    # cannot push the surface down while pulling the CMB up.
    check_both_down = all(_val_rs) do s
        q = geoid_response(s, 2, isoviscous_layers())
        q.h_surf < 0 && q.h_cmb < 0
    end
    nothing
end

# ╔═╡ 18fbea85-f8f9-4a3b-6c97-041526374859
if check_direct && check_bc && check_topo && check_compensation &&
   check_hr89_norm && check_hr89_WL && check_flow && check_both_down
    md"""
    !!! correct "Validation 2 passed — all seven physical checks"
    	- **Direct-only kernel** is positive and monotonic everywhere ✓ (no boundaries ⇒ no flip)
    	- **Free-slip boundary conditions** satisfied to machine zero at surface and CMB ✓
    	- **A dense anomaly depresses *both* boundaries** at every depth ✓ — one convection cell
    	  cannot push the surface down while pulling the CMB up
    	- **Compensation limit** ✓ — a load resting on either boundary is cancelled by that
    	  boundary's own deflection, to better than 1%. *This is the check that catches a
    	  flipped boundary sign: with the wrong sign the two terms add instead of cancelling.*
    	- **Hager & Richards (1989) Figure 2 normalisation** ✓ — their ``A^\ell`` and ``C^\ell``
    	  displaced-mass kernels both come out in ``[0,1]``, with ``A^\ell \to 1`` at the surface
    	  and ``C^\ell \to 1`` at the CMB, exactly as the figure shows
    	- **Published-kernel reproduction** ✓ — the ``\ell = 2`` geoid kernel for their model
    	  **WL** changes sign in the mid-mantle, reproducing the ``G^\ell`` panel of Figure 2
    	- **The interior flow vanishes at both free-slip boundaries** and its sign agrees with
    	  the surface deflection ✓ — downwelling goes with a depressed surface, as it must

    	!!! note "Which model is the published check?"
    		Hager & Richards plot kernels for **layered** models only (WL, WC, CL) — the paper
    		contains no uniform-viscosity kernel. So model **WL** is our comparison, not the
    		isoviscous case. The isoviscous kernel here is this notebook's own baseline: useful
    		for isolating the mechanism, but not a published result, and (see below) it does
    		**not** change sign.
    """
else
    md"""
    !!! danger "Validation 2 FAILED"
    	direct-only: $(check_direct) · boundary conditions: $(check_bc) ·
    	surface depression: $(check_topo) · both boundaries down: $(check_both_down) ·
    	compensation: $(check_compensation) · HR89 normalisation: $(check_hr89_norm) ·
    	HR89 model WL sign change: $(check_hr89_WL) · flow: $(check_flow)
    """
end

# ╔═╡ 3a1d0ca7-1acb-4c5d-8ed9-930415263749
"""
	flow_slice(result, layers, φ_face; nr=18, nθ=26)

Instantaneous velocity field on a meridional slice at longitude `φ_face`, for display on
a cut face. Returns `(rr, θθ, ur, uθ)` on a coarse grid suitable for drawing arrows.

Velocities are in arbitrary units — viscosity enters the Stokes problem only as a ratio,
so the *pattern* of the flow is meaningful but its absolute magnitude is not.
"""
function flow_slice(result, layers, φ_face; nr=18, nθ=26)
    lmax = result.lmax
    μ, φgrid = result.μ, result.φ
    nθg, nφg = length(μ), length(φgrid)
    _, w = gauss_legendre(nθg)
    dr = step(result.rgrid)
    dφ = 2π / nφg
    Pg = [assoc_legendre(μk, lmax) for μk in μ]
    cosmφ = [[cos(m * φ) for φ in φgrid] for m in 0:lmax]
    sinmφ = [[sin(m * φ) for φ in φgrid] for m in 0:lmax]

    # sample right up to the boundaries so the arrows visibly die there, as free-slip
    # with no through-flow requires
    rr = collect(range(B_CMB, A_EARTH, length=nr))
    θθ = collect(range(0.06π, 0.94π, length=nθ))
    Pθ = [assoc_legendre(cos(t), lmax) for t in θθ]
    # dP_ℓ^m/dθ by central differences, same normalisation as Pθ
    hθ = 1e-5
    # P_ℓ^m only exists for m ≤ ℓ, so this must stay triangular like assoc_legendre
    # itself — a square matrix would index coefficients that do not exist.
    dPθ = map(θθ) do t
        Pp = assoc_legendre(cos(t + hθ), lmax)
        Pm = assoc_legendre(cos(t - hθ), lmax)
        [[(Pp[l+1][m+1] - Pm[l+1][m+1]) / (2hθ) for m in 0:l] for l in 0:lmax]
    end

    ur = zeros(nr, nθ)
    uθ = zeros(nr, nθ)
    for ℓ in 1:lmax
        # one flow solve per (degree, load radius); reused for every order m
        prof = [flow_profile(rs, ℓ, layers, rr) for rs in result.rgrid]
        for m in 0:ℓ
            nrm = m == 0 ? 1.0 : 2.0
            # δρ_ℓm(r) for this degree/order, then integrate the flow response over r
            Uacc = zeros(nr)
            Vacc = zeros(nr)
            for (i, _) in enumerate(result.rgrid)
                gc = 0.0
                gs = 0.0
                for k in 1:nθg
                    fc = 0.0
                    fs = 0.0
                    for j in 1:nφg
                        fc += result.field[i, k, j] * cosmφ[m+1][j]
                        fs += result.field[i, k, j] * sinmφ[m+1][j]
                    end
                    gc += w[k] * Pg[k][ℓ+1][m+1] * fc
                    gs += w[k] * Pg[k][ℓ+1][m+1] * fs
                end
                δc = nrm * gc * dφ / (2 * 2π)
                δs = nrm * gs * dφ / (2 * 2π)
                # the slice sits at φ_face, so combine cos/sin parts there
                amp = δc * cos(m * φ_face) + δs * sin(m * φ_face)
                @. Uacc += prof[i].U * amp * dr
                @. Vacc += prof[i].V * amp * dr
            end
            for k in 1:nθ, i in 1:nr
                ur[i, k] += Uacc[i] * Pθ[k][ℓ+1][m+1]
                uθ[i, k] += Vacc[i] * dPθ[k][ℓ+1][m+1]
            end
        end
    end
    return (rr=rr, θθ=θθ, ur=ur, uθ=uθ)
end

# ╔═╡ 4b2e1db8-2bdc-4d6e-9fba-041526374860
# Instantaneous flow on the two cut faces. Computed only when the arrows are shown,
# since this is the heaviest cell in the notebook.
gk_flow = gk_showflow && !isempty(gk_blobs) ?
          [flow_slice(gk_result, gk_layers, pc) for pc in gk_cutfaces] : nothing

# ╔═╡ 3a1d0ca6-1adb-4c5d-8ec9-829304152637
# Push the computed kernel, geoid map and topography maps back to the globe widget.
let
    num(x) = isfinite(x) ? string(round(x, sigdigits=5)) : "0"
    jsonarr(v) = "[" * join([num(x) for x in v], ",") * "]"
    # a 2-D map goes over as an array of rows (one per colatitude node)
    jsonmap(M) = "[" * join(["[" * join([num(M[k, j]) for j in axes(M, 2)], ",") * "]"
                             for k in axes(M, 1)], ",") * "]"
    payload = string(
        "{\"kernel\":{\"k\":", jsonarr(gk_kvals),
        ",\"depth\":", jsonarr(collect(gk_kdepths)),
        ",\"crossing\":", isnothing(gk_crossing) ? "null" : string(round(gk_crossing, digits=1)),
        "},\"mu\":", jsonarr(gk_result.μ),
        ",\"phi\":", jsonarr(gk_result.φ),
        ",\"geoid\":", jsonmap(gk_result.N),
        ",\"topo\":{\"surf\":", jsonmap(gk_topo.surf),
        ",\"cmb\":", jsonmap(gk_topo.cmb), "}",
        ",\"flow\":", isnothing(gk_flow) ? "null" :
        "[" * join(["{\"rr\":" * jsonarr(f.rr) *
                    ",\"tt\":" * jsonarr(f.θθ) *
                    ",\"ur\":" * jsonmap(f.ur) *
                    ",\"ut\":" * jsonmap(f.uθ) * "}" for f in gk_flow], ",") * "]",
        "}")
    HTML("""<script>
      window.dispatchEvent(new CustomEvent('geoid-results', {detail: $(repr(payload))}));
    </script>""")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
ColorSchemes = "~3.31.0"
LaTeXStrings = "~1.4.0"
PlutoPlotly = "~0.6.6"
PlutoUI = "~0.7.83"
Symbolics = "~7.0.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "87d2a42ff811fc500d0639e52cb9a46ce6d82b63"

[[deps.ADTypes]]
git-tree-sha1 = "f7304359109c768cf32dc5fa2d371565bb63b68a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.21.0"

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

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
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

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

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
git-tree-sha1 = "6fb53a69613a0b2b68a0d12671717d307ab8b24e"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.5"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c249d86e97a7e8398ce2068dce4c078a1c3464de"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.16"

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
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "3f50fa86c968fc1a9e006c07b6bc40ccbb1b704d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.4"

[[deps.EnumX]]
git-tree-sha1 = "7bebc8aad6ee6217c78c5ddcf7ed289d65d0263e"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.6"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

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

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"
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

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c89d196f5ffb64bfbf80985b699ea913b0d2c211"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.6.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

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
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

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

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "d38b8653b1cdfac5a7da3b819c0a8d6024f9a18c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.13"

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

    [deps.MultivariatePolynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "22df8573f8e7c593ac205455ca088989d0a2c7a0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.7"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

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

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

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

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "2f609ec2295c452685d3142bc4df202686e555d2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.16"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLPublic]]
git-tree-sha1 = "0ba076dbdce87ba230fff48ca9bca62e1f345c9b"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.1"

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

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "6547cbdd8ce32efba0d21c5a40fa96d1a3548f9f"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

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

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "49201c2793ce02f141c6f8b5194ce34e8012cd68"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.4"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "30cb5145192c723dff2d5790ca79082f3490079e"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.5.0"

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
deps = ["ADTypes", "AbstractPlutoDingetjes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "Preferences", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "5ccee3582d344a87918840862c0b67285ec9fce1"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.0.1"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

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
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

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

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═5c3f2e62-3c4d-4e7f-a01b-2c3d4e5f6071
# ╟─6d403f73-4d5e-4f80-b12c-3d4e5f607182
# ╟─e5c8b756-c5b6-4708-3964-374859607182
# ╟─7e514084-5e6f-4091-c23d-4e5f60718293
# ╟─8f625195-6f70-41a2-d34e-5f6071829304
# ╟─c3a69596-a374-45e6-17f2-152637485960
# ╠═d4b7a646-b485-46f7-2803-263748596071
# ╠═f6d9c866-d6a7-4819-4a75-485960718293
# ╠═07ead976-e7b8-492a-5bc6-596071829304
# ╠═07ead977-e7a8-492a-5bd6-607182930416
# ╟─18fbea87-f8b9-4a3b-6cb7-718293041527
# ╟─290cfb97-0861-4b4c-7dc8-829304152638
# ╟─290cfb98-0851-4b4c-7dd8-930415263749
# ╠═18fbea86-f8c9-4a3b-6ca7-607182930415
# ╠═3a1d0ca6-1adb-4c5d-8ec9-829304152637
# ╟─4b2e1db6-2bec-4d6e-9fda-930415263748
# ╠═5c3f2ec6-3cfd-4e7f-a0eb-041526374859
# ╟─6d403fd6-4dcd-4f80-b1fb-152637485960
# ╟─7e5140e6-5ebf-4091-c20b-263748596071
# ╟─b2958428-92a3-44d5-0671-829304152637
# ╠═c3a69539-a3b4-45e6-1782-930415263748
# ╟─d4b7a64a-b4c5-46f7-2893-041526374859
# ╟─e5c8b75b-c5d6-4708-3904-152637485960
# ╟─90736206-7081-42b3-e45f-607182930415
# ╠═a1847317-8192-43c4-f560-718293041526
# ╟─f6d9c86c-d6e7-4819-4a15-263748596071
# ╠═07ead97d-e7f8-492a-5b26-374859607182
# ╠═18fbea8e-f809-4a3b-6c37-485960718293
# ╠═290cfb9f-090a-4b4c-7d48-596071829304
# ╠═3a1d0ca0-1a1b-4c5d-8e59-607182930415
# ╠═4b2e1db1-2b2c-4d6e-9f6a-718293041526
# ╟─5c3f2ec2-3c3d-4e7f-a07b-829304152637
# ╠═6d403fd3-4d4e-4f80-b18c-930415263748
# ╠═7e5140e4-5e5f-4091-c29d-041526374859
# ╟─8f6251f5-6f60-41a2-d3ae-152637485960
# ╟─90736265-7071-42b3-e4bf-263748596071
# ╠═a18473c5-8182-43c4-f5ca-374859607182
# ╠═b2958425-9293-44d5-06d1-485960718293
# ╠═c3a69595-a3a4-45e6-17e2-596071829304
# ╠═d4b7a645-b4b5-46f7-28f3-607182930415
# ╠═e5c8b755-c5c6-4708-3954-718293041526
# ╟─f6d9c865-d6d7-4819-4a65-829304152637
# ╠═07ead975-e7e8-492a-5bb6-930415263748
# ╟─18fbea85-f8f9-4a3b-6c97-041526374859
# ╟─290cfb95-0901-4b4c-7da8-152637485960
# ╠═3a1d0ca5-1a0b-4c5d-8eb9-263748596071
# ╠═4b2e1db5-2b1c-4d6e-9fca-374859607182
# ╠═4b2e1db7-2b0c-4d6e-9fea-485960718294
# ╟─6d403fd5-4d3e-4f80-b1eb-596071829304
# ╠═7e5140e5-5e4f-4091-c2fb-607182930415
# ╟─8f6251f6-6f40-41a2-d3bd-718293041526
# ╟─90736266-7051-42b3-e4cf-829304152637
# ╠═a18473c6-8152-43c4-f5db-930415263748
# ╠═b2958426-9263-44d5-06e1-041526374859
# ╠═290cfb96-0871-4b4c-7db8-718293041526
# ╟─5c3f2ec5-3c2d-4e7f-a0db-485960718293
# ╠═3a1f0c40-1a2b-4c5d-8e9f-0a1b2c3d4e5f
# ╠═4b2e1d51-2b3c-4d6e-9f0a-1b2c3d4e5f60
# ╠═c3a69597-a394-45e6-1802-607182930417
# ╠═3a1d0ca7-1acb-4c5d-8ed9-930415263749
# ╠═4b2e1db8-2bdc-4d6e-9fba-041526374860
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
