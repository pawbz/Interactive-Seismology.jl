### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Viscoelastic Rheology: Elastic, Viscoelastic, or Viscous?"
#> tags = ["elasticity", "rheology", "viscoelasticity"]
#> layout = "layout.jlhtml"
#> description = "Pull a block, release it, repeat -- then discover that whether a material acts elastic, viscoelastic, or viscous depends on the clock you're comparing it to, not just its recipe."

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

# ╔═╡ 46fa4dc1-b407-4584-91bf-d136b77622bc
begin
    using PlutoUI
end

# ╔═╡ 52597587-4c6b-4420-8242-6dbe166e280d
TableOfContents()

# ╔═╡ faf8d654-5c42-46a5-967d-93c3897273e1
md"""
# Why Does a Rod Stretch Proportionally -- and What Happens If You Pull Too Hard?

Hooke's law, ``F \propto \Delta \ell``, is usually the very first thing anyone learns
about elastic materials. It is also, on its own, an oddly incomplete statement: a real
rod has a length, a cross-section, and a material, and none of those show up in
``F = k\Delta\ell`` -- only a single number, ``k``, does. Is ``k`` really *the* constant
of the material, or is it secretly a statement about how big your particular rod happens
to be?

And Hooke's law only describes the *elastic* part of the story. Pull hard enough on
almost anything and it stops springing back -- it yields, deforms permanently, and
"releases" you along a different line than the one you loaded it on.

This notebook builds both ideas from a single hands-on widget: pull a block, release it,
repeat. Then it asks a harder question -- when does a material even *behave* elastically
in the first place? The same block, the same underlying physics, can look like a seismic
wave, a locked fault accumulating strain, a relaxing postseismic transient, or slowly
convecting mantle -- depending entirely on how the timescale you're asking about compares
to two clocks the material carries internally. This notebook is the physical intuition
behind that comparison.
"""

# ╔═╡ 6444aa15-39e9-4f0f-9e90-809899dd19ff
md"""
## Pull It, Release It, Repeat

Drag the handle in free-play mode -- pull it out, then ease it back. Small pulls spring
back exactly the way they came: the point on the graph below traces the same line up and
down. Pull harder, past some threshold, and the material *yields* -- the friction-slider
block visibly slips, and easing back off no longer retraces the same path. It comes back
down a *parallel* line, stopping short of zero strain: a permanent, plastic set.

Do this several times, at different peak pulls. The stress-strain graph keeps every past
cycle (faintly) alongside the current one, so the pattern you're mapping out -- where the
line stops being straight, how steep the recovery slope is, how much permanent strain
each cycle leaves behind -- accumulates as *your own* record, not something stated up
front.

The block is drawn as a lattice of small fibers, not a plain bar, on purpose: a real rod
is secretly made of many such springs, running both along its length and across it, and
it's worth asking what happens to the *slope* of this graph as you imagine subdividing the
rod more finely. See the aside below.

Watch the lattice narrow as you pull it, too -- the vertical fibers show the lateral
response, the same fibers as the horizontal ones, just felt as a 2-D continuum instead of
a 1-D chain. While the material is purely elastic, that narrowing is governed by ``\nu``
(Poisson's ratio), the second elastic constant alongside ``E``. But once the material
starts yielding or creeping, the narrowing visibly speeds up -- real inelastic flow is
close to volume-preserving regardless of ``\nu``, which is the actual reason a stretched
rod visibly *necks* just before it breaks.

Both ``E`` and ``\nu`` are set by scrubbing the material itself, not a slider -- drag the
spring (stiffness) or the dashpot (viscosity) between the wall and the block, or the
*bracing* strip below, and watch the lattice actually change shape to match. See the last
aside below for why bracing, specifically, is what ``\nu`` is really about.
"""

# ╔═╡ d3f8a2c1-7e5b-4c9a-8d16-2f4a6e8c0b31
md"""
## The Mirror Image: Fixing the Strain Instead of the Stress

Every pull above was **force-controlled**: the handle sets the applied stress, and
whatever strain results is left to creep further the longer you hold a steady pull. That
is a genuinely different experiment from squeezing something to a fixed shape and
clamping it there -- **displacement-controlled**. Switch the toggle above the block to
"Stretch & hold" to try it: now the handle sets the *strain* directly, and once you stop
moving it (hold the mouse steady, exactly the same "hold" gesture as before), simulated
time advances again -- except this time it's the *stress* that visibly relaxes, decaying
toward zero while your strain stays pinned exactly where you left it.

This is the textbook **stress relaxation** test, the mirror image of the **creep** test
above: same Maxwell dashpot, same material, but which quantity is held fixed and which is
left to respond has flipped. Held at a fixed strain, stress decays exponentially,
``\sigma(t)=\sigma_0e^{-t/\tau_{\rm relax}}`` (`stress_from_strain_with_relaxation` in the
Appendix), with its own relaxation time ``\tau_{\rm relax}=3\eta/E`` -- close to, but not
identical to, the shear-based ``\tau_M=\eta/\mu`` from the timescale strip above the block
(the two differ by exactly ``3/(2(1+\nu))``, the same Trouton-ratio bookkeeping from the
aside below, just applied to a different clock). Geophysically, this is the more literal
picture of what happens right after an earthquake locks a fault segment at a new, fixed
offset: it's the *stress* that offset left behind which relaxes over the postseismic years
that follow, not a stress someone keeps re-applying by hand.
"""

# ╔═╡ 3e1a5c9a-8b7f-4a2f-9d1e-6f5c4b8e2d3a
md"""
## Two Clocks, Compared to Whatever Clock You're Asking About

Pulling the block by hand is really the *interseismic* part of the seismic cycle in
disguise, per the first aside below: pull slowly, hold, release. But "slowly" compared to
*what*? Start from the one law that's always true regardless of material -- momentum
conservation, ``\rho\ddot{\mathbf u} = \nabla\cdot\boldsymbol\sigma``, which says nothing
yet about whether the material is elastic, viscous, or plastic. Rescale space by a length
``L``, time by a duration ``T``, and displacement by ``U`` (its 1-D shear form is
``\rho u_{tt} = \mu u_{xx}``) and a single dimensionless group falls out in front of the
inertia term:

```math
I = \frac{L}{V_sT} = \frac{\tau_{\rm wave}}{T}, \qquad
\tau_{\rm wave} = \frac{L}{V_s}, \qquad V_s = \sqrt{\mu/\rho}.
```

``\tau_{\rm wave}`` is how long a shear wave takes to cross a body of size ``L`` -- how
quickly the material can even *find out* that a force was applied somewhere else in it.
``I \sim 1`` means inertia and elastic restoring forces are comparable: a seismic wave.
``I \ll 1`` means elastic information crosses the region far faster than the loading
changes, so ``\nabla\cdot\sigma\approx0`` at every instant -- quasi-statics.

Do the same rescaling to the Maxwell rheology sitting in this block's dashpot -- but the
dashpot is fundamentally a *shear* element (it only ever relaxes the deviatoric,
shape-changing part of the stress), while the handle applies a plain uniaxial pull. Split
that pull into its deviatoric and volumetric parts (`` [`axial_maxwell_creep_rate`](@ref)
``, Appendix) and only ``2/3`` of it -- the deviatoric third -- ever reaches the dashpot;
the rest is a lossless volume change the dashpot never feels. Working through that split
turns the rheology into:

```math
\dot\varepsilon = \frac{\dot\sigma}{E} + \frac{\sigma}{3\eta}, \qquad
De = \frac{\tau_M}{T}, \qquad \tau_M = \frac{\eta}{\mu}.
```

``\eta`` here is the material's *shear* viscosity -- the same one that sets
``\tau_M=\eta/\mu`` -- so the ``3`` in front of it is not a fitted fudge factor. It is the
classic **Trouton ratio**: pulling a Newtonian material along one axis resists exactly
``3\times`` as hard, per unit strain rate, as shearing it does, purely from how much of a
uniaxial pull is actually shear. Skip that factor (an easy mistake -- the 1-D
spring-and-dashpot cartoon doesn't show it) and matching this block's creep to a real
mantle viscosity needs a fictitious, ``3\times``-too-small ``\eta``.

``De \gg 1``: the dashpot hasn't had time to move, ``\sigma\approx E\varepsilon`` --
elastic. ``De \ll 1``: plenty of time to relax, ``\sigma\approx 3\eta\dot\varepsilon`` --
viscous. ``De\sim1``: both terms matter -- viscoelastic. ``I`` and ``De`` are *independent*
questions -- one asks whether inertia matters, the other whether flow matters -- so
together they sort every regime this widget can reach:

| ``T`` vs. the two clocks | Regime | Preset |
|---|---|---|
| ``T\sim\tau_{\rm wave}`` (``I\sim1``) | elastic, dynamic | *(needs real inertia -- not modeled here, see the note)* |
| ``T\gg\tau_{\rm wave}``, ``De\gg1`` | elastic, quasi-static | Elastic Loading, Interseismic Loading |
| ``T\gg\tau_{\rm wave}``, ``De\sim1`` | viscoelastic | Postseismic Relaxation |
| ``T\gg\tau_{\rm wave}``, ``De\ll1`` | quasi-static, viscous | Mantle Convection |

The timescale strip above the block always plots ``\tau_{\rm wave}`` and ``\tau_M`` as
markers on one log-scale number line spanning a hundredth of a second to hundreds of
millions of years, alongside whichever ``T`` is currently in play: a preset's own realistic
duration, or -- while you hold a hand-pull -- *your own* pull's duration, live. The four
presets each fast-forward through a scripted, realistic forcing over a few real seconds so
you can actually watch ``De\gg1``, ``De\sim1``, and ``De\ll1`` play out (Elastic Loading and
Interseismic Loading both sit at ``De\gg1``, just at very different durations); a genuinely
elastodynamic ``I\sim1`` demonstration would need a real wave equation with mass in it,
which this quasi-static (massless) engine deliberately doesn't attempt.

!!! note "Your hand cannot honestly reach these timescales -- so it's stretched"
    A real hold can only ever last a few seconds, but ``\tau_M`` for real rock spans years
    to millions of years. Rather than pretending a few seconds of dragging *are* geological
    time (which would need an unrealistically tiny ``\eta`` to show anything, exactly the
    compromise this notebook used to make), holding the handle maps your hold's real
    duration, on a log scale, across the *entire* timescale axis -- a ~12 second hold
    sweeps from ~0.01 s of simulated time to ~300 million years. That simulated duration,
    not your real hold time, is what actually gets fed to the creep calculation, so holding
    long enough genuinely creeps, honestly, at whatever ``\eta`` you've dialed in.
"""

# ╔═╡ 64f12034-4161-4d2b-b1c2-1809ee2e6e2c
md"""
## Three Asides: What the Springs Are Really Doing, and What This Looks Like as a Fault

**Why the slope, not the spring, is the material property.** Imagine the rod really is
``n`` identical springs of stiffness ``k`` in series, spaced ``dx`` apart, so its total
length is ``L = n\,dx``. Each spring carries the same force ``F`` and stretches ``F/k``,
so the whole rod stretches ``nF/k``. Dividing by area to get stress and by length to get
strain:

```math
\frac{F/A}{(nF/k)/(n\,dx)} = \frac{k\,dx}{A}
```

the ``n`` cancels. Subdivide the same rod into twice as many, half as stiff springs and
you measure the *same* slope. The raw spring constant ``k`` depends on how big a chunk you
happened to cut; ``k\,dx/A`` -- the modulus -- does not. See `effective_modulus` in the
Appendix for the same statement checked numerically across several chain lengths.

**Elastic rebound and the seismic cycle.** Relabel the widget: the handle is a point on a
locked fault, slowly dragged by plate motion. Pulling it out is the *interseismic* period
-- elastic strain quietly accumulating for years to centuries. Reaching the yield point
and dropping back down the unloading slope is the *earthquake* -- a sudden release of the
stored elastic strain. That vertical drop, from peak stress to the post-slip stress, is
exactly the seismological quantity **stress drop**. Repeating the cycle is the seismic
cycle itself. This is Reid's 1910 elastic rebound theory, and the widget already produces
every piece of it. (The yield condition here is a 1-D stand-in for the full Mohr-Coulomb
failure criterion explored interactively in `stress-tensor.jl`.)

**Two materials, same springs, different Poisson's ratio.** Picture two lattices with the
*exact* same number of springs per unit area and the *exact* same spring constant --
arranged differently:

```
Material A          Material B

o---o---o            o---o---o
|   |   |            |\ /|\ /|
|   |   |            | X | X |
|   |   |            |/ \|/ \|
o---o---o            o---o---o
```

Material A is a plain grid -- horizontal and vertical springs only. Material B adds
diagonal braces in every cell. Pull Material A sideways: every horizontal spring only
ever senses motion *along itself*, so nothing anywhere is ever pushed or pulled in the
vertical direction -- not approximately, not a little, *at all*. Its Poisson's ratio is
exactly zero, for a structural reason, not because anyone assumed it. Material B's
diagonals are the only springs that feel *both* directions at once, so pulling
horizontally now genuinely drags the lattice inward vertically too -- a nonzero, derivable
``\nu``. `lattice_homogenized_nu` in the Appendix solves exactly this small mechanics
problem, and its self-check confirms Material A's ``\nu=0`` directly from the solver, not
by assumption. Material A itself isn't one of the three bracing levels above, though --
an unbraced pin-jointed grid is a genuine *mechanism* (nothing stops it from freely
shearing), not something you could hand someone as a material. All three scrubbable
levels are braced, differing only in how much, so every one of them is a real, stable
lattice with its own honestly-derived ``\nu``.
"""

# ╔═╡ 4a425525-dba3-4dcc-b9f6-7deaf4fb0588
md"""
## The Other Two Numbers: Yielding and Flow

``I`` and ``De`` both ask *"is there enough time?"* -- for elastic information to
propagate, or for viscous flow to relax stress. Two more questions round out the picture,
and neither one is controlled by time at all.

**Does it yield?** ``\sigma_y`` and the hardening modulus ``H``, already scrubbable above,
are exactly the material-strength side of a third, genuinely independent ratio:

```math
\Pi_Y = \frac{\sigma}{\sigma_y}.
```

``\Pi_Y \ll 1``: safely below strength, whatever ``I`` or ``De`` happen to be.
``\Pi_Y \sim 1``: plastic (or brittle) failure begins -- a slowly-loaded, elastic
(``De\gg1``) fault still ruptures the instant ``\Pi_Y`` reaches 1. Pull hard enough in
free-play and you cross this threshold directly; the friction-slider block visibly
slipping *is* ``\Pi_Y\sim1``.

**Once it's flowing, does inertia matter again?** Deep in the ``De\ll1`` corner the
material creeps like a fluid, and it's more natural to ask about the flow itself rather
than about waves crossing a solid:

```math
Re = \frac{\rho V L}{\eta}.
```

Mantle convection has flow speeds ``V\sim\text{cm/yr}`` and ``\eta\sim10^{21}\,\text{Pa·s}``
-- an almost inconceivably small ``Re``. That's creeping (Stokes) flow, not turbulence,
which is exactly why the Mantle Convection preset above never looks anything like a
turbulent fluid despite moving enormous volumes of rock over its own ``\tau_M``. This
notebook doesn't compute ``Re`` directly (there's no actual flow velocity in a pulled
block, only a strain history), but it's the reason that regime is safely creeping rather
than chaotic.
"""

# ╔═╡ 9b7df3b0-c0f0-4c96-8314-23842a2d97dd
md"""
## References

- Reid, H.F. (1910). *The mechanism of the earthquake*, the origin of elastic rebound
  theory.
- The Mohr-Coulomb failure criterion referenced above is explored interactively in
  `stress-tensor.jl`.
- The dilatation/deviatoric strain split behind the lateral (necking) strain formula
  comes from `strain-tensor.jl`.
- The spring-chain visual reuses the drawing routine from `coupled_oscillations.jl`.
"""

# ╔═╡ df3cccd7-5a20-452b-89e9-1f20b818808a
md"## Appendix"

# ╔═╡ 6ab7bf3b-19e1-477a-b046-88df986855cb
md"## Discovering Elastic and Plastic Response"

# ╔═╡ 33fb6caf-2564-41bc-bd5f-8fba3670d887
"""
	elastoplastic_response(strain_path, E, sigma_y0; H=0.0)

Integrate a 1-D elastic-(perfectly-)plastic rheology -- a spring in series with a
Coulomb-friction slider -- along a prescribed strain history `strain_path`, using a
standard return-mapping algorithm: a trial stress `E*(ε - εₚ)` is accepted if it stays
within the yield surface `|σ| ≤ σ_y0 + H|εₚ|`, otherwise the plastic strain `εₚ` is
advanced just enough to bring the stress back onto that surface.

Returns `(stress_path, plastic_strain_path)`, both the same length as `strain_path`.
Setting `H=0` gives perfectly-plastic behavior (yield stress never changes); `H>0` gives
linear kinematic hardening (the yield stress grows with accumulated plastic strain).
"""
function elastoplastic_response(strain_path, E, sigma_y0; H=0.0)
    n = length(strain_path)
    stress = zeros(Float64, n)
    eps_p_path = zeros(Float64, n)
    eps_p = 0.0
    for i in 1:n
        eps = strain_path[i]
        trial = E * (eps - eps_p)
        yield_now = sigma_y0 + H * abs(eps_p)
        if abs(trial) <= yield_now
            stress[i] = trial
        else
            s = sign(trial)
            eps_p = (E * eps - s * sigma_y0) / (E + H)
            stress[i] = E * (eps - eps_p)
        end
        eps_p_path[i] = eps_p
    end
    return stress, eps_p_path
end

# ╔═╡ 6cb2a1e4-9f61-4b8a-9c3e-7d5b6a2e4f10
"""
	strain_from_stress(stress_path, E, sigma_y0; H=0.0, plastic_cap=1.5)

The force-controlled dual of [`elastoplastic_response`](@ref): given a *prescribed stress*
history `stress_path` (the widget drags an applied force, not a displacement), return the
resulting `(strain_path, plastic_strain_path)`. For each target stress `σ`, if `|σ|` is
within the current yield surface `σ_y0 + H|εₚ|`, the response is purely elastic from the
current plastic strain: `ε = εₚ + σ/E`. If `|σ|` would exceed that surface, plastic strain
advances just enough that the yield surface catches up to `σ`.

!!! warning "Perfect plasticity has no equilibrium strain under load control"
    With `H=0`, no finite strain can push the stress past `σ_y0` -- physically, holding a
    load at yield causes the material to flow without bound. `plastic_cap` sets the
    largest plastic strain this function will ever report (default `1.5`, i.e. 150%),
    which implicitly clamps the *achievable* stress to `σ_y0 + H·plastic_cap` -- exactly
    `σ_y0` when `H=0`. Callers should clamp `stress_path` to that same bound before calling,
    so the clamp is visible in the input, not silently absorbed here.
"""
function strain_from_stress(stress_path, E, sigma_y0; H=0.0, plastic_cap=1.5)
    n = length(stress_path)
    strain = zeros(Float64, n)
    eps_p_path = zeros(Float64, n)
    eps_p = 0.0
    # H=0 ("perfectly plastic", this notebook's default) is a numerical singularity for
    # this algebraic update: with the OLD `max(H, 1e-9)` floor, sigma_max = sigma_y0 +
    # H*plastic_cap collapsed to exactly sigma_y0, so `sig` was clamped there BEFORE the
    # yield check ran -- `abs(sig) > yield_now` (both sigma_y0) could then never be true,
    # and eps_p stayed 0 forever regardless of how hard stress_path pushed past yield.
    # Regularizing H with a floor tied to E instead gives sigma_max real headroom above
    # sigma_y0, so the yield check can actually fire, and (since eps_p's formula below uses
    # the CLAMPED sig, which the same sigma_max bounds) eps_p is naturally capped at exactly
    # plastic_cap once sig saturates -- matching this function's own docstring.
    H_eff = max(H, 0.1 * E)
    sigma_max = sigma_y0 + H_eff * plastic_cap
    for i in 1:n
        sig_raw = stress_path[i]
        sig = clamp(sig_raw, -sigma_max, sigma_max)
        yield_now = sigma_y0 + H_eff * abs(eps_p)
        # check the RAW attempted stress, not the already-clamped one -- clamping first and
        # then checking the clamped value against the same bound it was clamped to is
        # exactly the bug this replaces
        if abs(sig_raw) > yield_now
            s = sign(sig_raw)
            eps_p = s * (abs(sig) - sigma_y0) / H_eff
        end
        strain[i] = eps_p + sig / E
        eps_p_path[i] = eps_p
    end
    return strain, eps_p_path
end

# ╔═╡ 5e2c9a3f-1b6d-4f2a-9c0e-7a3f8d1b6e42
md"### Why the Dashpot Only Feels a Third of the Pull"

# ╔═╡ 7f4a1d8e-2c5b-4e9f-a0d3-6b1c8e4f2a97
"""
	axial_maxwell_creep_rate(sigma, eta)

The creep-strain RATE a Maxwell dashpot contributes under a plain uniaxial stress
`sigma`, given the material's *shear* viscosity `eta` -- the same `eta` that sets
``\\tau_M = \\eta/\\mu`` in [`wave_relaxation_clocks`](@ref). The dashpot is fundamentally
a shear element: it only ever relaxes the deviatoric (shape-changing) part of the stress
tensor, ``\\dot e_{ij} = \\dot s_{ij}/(2\\mu) + s_{ij}/(2\\eta)``. A uniaxial stress
`sigma` (lateral stresses held at zero) splits into a purely elastic, lossless
volumetric part ``p=\\sigma/3`` and a deviatoric part ``s_{11}=\\sigma-p=2\\sigma/3``;
only that deviatoric third ever reaches the dashpot, giving
``\\dot\\varepsilon_{\\rm creep}=\\dot e_{11}=s_{11}/(2\\eta)=\\sigma/(3\\eta)``.

The `3` is not a fitted fudge factor -- it's the classic *Trouton ratio*: pulling a
Newtonian material along one axis resists exactly ``3\\times`` as hard, per unit strain
rate, as shearing it does. Used by [`strain_from_stress_with_creep`](@ref); multiply by
an elapsed time and accumulate to get creep strain.
"""
axial_maxwell_creep_rate(sigma, eta) = sigma / (3 * eta)

# ╔═╡ a0900587-c705-4fb9-99b6-002b4fe257c8
"""
	strain_from_stress_with_creep(stress_path, time_path, E, sigma_y0, eta; H=0.0, plastic_cap=1.5)

The "silly putty" extension of [`strain_from_stress`](@ref): the same yield-capped
elastic/plastic response, plus a Maxwell dashpot in series that creeps continuously
under *any* sustained nonzero stress, at rate [`axial_maxwell_creep_rate`](@ref)`(σ,
η) = σ/(3η)` -- the Trouton-ratio-corrected rate for a *uniaxial* stress acting on a
dashpot that is fundamentally a shear element (see the derivation there and in
[`wave_relaxation_clocks`](@ref)). Unlike yielding (which depends only on how *hard* you
pull), creep depends on how *long* you hold a given pull: integrated here as
`Δε_creep = σ·Δt/(3η)` between successive samples of `time_path` (seconds).

Returns `(strain_path, residual_path)`, where `residual_path` is the *combined*
plastic-plus-creep strain left behind once stress returns to zero -- the quantity the
widget's release animation snaps back to, since only the elastic part is recoverable.

Setting `η` very large makes the creep term negligible over any realistic drag duration,
recovering [`strain_from_stress`](@ref) exactly -- this is why the widget's default `η`
is large: creep is an *additional* effect layered on top of yielding, off by default.
"""
function strain_from_stress_with_creep(stress_path, time_path, E, sigma_y0, eta; H=0.0, plastic_cap=1.5)
    n = length(stress_path)
    strain = zeros(Float64, n)
    residual_path = zeros(Float64, n)
    eps_p = 0.0
    eps_c = 0.0
    # H=0 ("perfectly plastic", the default everywhere in this notebook) is a genuine
    # numerical singularity for this algebraic (not incremental) plasticity update: real
    # perfectly-plastic flow is UNBOUNDED at constant stress, so dividing by H directly
    # would blow up. Regularize with a floor tied to E, not a bare epsilon -- 1e-9 (the old
    # floor) makes the post-yield slope ~1e9x steeper than the elastic one, so `abs(sig) >
    # yield_now` below could never fire in the first place without ALSO clamping sig to
    # sigma_max = sigma_y0 (since plastic_cap*H=0 too), which is exactly why yielding never
    # triggered before this fix. 0.1E keeps the post-yield curve visibly bent but still
    # finite -- about a 10x-shallower rise than the elastic slope.
    H_eff = max(H, 0.1 * E)
    sigma_max = sigma_y0 + H_eff * plastic_cap
    for i in 1:n
        sig_raw = stress_path[i]
        sig = clamp(sig_raw, -sigma_max, sigma_max)
        yield_now = sigma_y0 + H_eff * abs(eps_p)
        # check the RAW attempted stress, not the already-clamped one -- clamping first
        # and then checking the clamped value against the same bound it was clamped to
        # is exactly the bug this replaces
        if abs(sig_raw) > yield_now
            s = sign(sig_raw)
            eps_p = s * (abs(sig) - sigma_y0) / H_eff
        end
        if i > 1
            dt = max(time_path[i] - time_path[i-1], 0.0)
            eps_c += axial_maxwell_creep_rate(sig, eta) * dt
        end
        strain[i] = eps_p + eps_c + sig / E
        residual_path[i] = eps_p + eps_c
    end
    return strain, residual_path
end

# ╔═╡ f3a1b2c4-6d5e-4a8b-9c1f-2e7d4b6a8f01
md"### Squeezing Instead of Pulling: Stress Relaxation Under Fixed Strain"

# ╔═╡ a7c3e9f1-4b2d-4e6a-8f19-3d5c7b9e1a02
"""
	axial_stress_relaxation_time(E, eta)

The axial Maxwell relaxation time ``\\tau_{\\rm relax}=3\\eta/E`` that
[`stress_from_strain_with_relaxation`](@ref) decays with once strain is held fixed --
the constant-*strain* dual of the constant-*stress* creep rate
[`axial_maxwell_creep_rate`](@ref), both carrying the same Trouton-ratio ``3``. Related to
the shear-based ``\\tau_M=\\eta/\\mu`` from [`wave_relaxation_clocks`](@ref) by
``\\tau_{\\rm relax}=\\tau_M\\cdot 3/(2(1+\\nu))`` -- the same dashpot, timed two different
ways (a uniaxial pull vs. a pure shear), so the two clocks agree in order of magnitude but
are not numerically identical. See the self-check below.
"""
axial_stress_relaxation_time(E, eta) = 3 * eta / E

# ╔═╡ b8d4f0a2-5c3e-4f7b-9a2d-4e6d8c0f2b13
"""
	stress_from_strain_with_relaxation(strain_path, time_path, E, eta; sigma_cap=Inf)

The displacement-controlled dual of [`strain_from_stress_with_creep`](@ref): given a
*prescribed total strain* history `strain_path` (the widget clamps the block at a fixed
elongation instead of pulling with a fixed force), return the stress the Maxwell element
carries as a result.

Differentiating the same decomposition used everywhere in this notebook,
``\\varepsilon=\\sigma/E+\\varepsilon_c``, ``\\dot\\varepsilon_c=\\sigma/(3\\eta)`` (the
Trouton-corrected rate from [`axial_maxwell_creep_rate`](@ref)), at *prescribed*
``\\varepsilon(t)`` gives a first-order linear ODE for ``\\sigma``:

```math
\\dot\\sigma + \\frac{\\sigma}{\\tau} = E\\dot\\varepsilon, \\qquad
\\tau = \\tau_{\\rm relax} = \\frac{3\\eta}{E}
```

(see [`axial_stress_relaxation_time`](@ref)). Rather than a forward-Euler step, which
would need a fine ``\\Delta t/\\tau`` to avoid drifting, each sample interval is integrated
*exactly*, treating the strain rate as constant (``r=\\Delta\\varepsilon/\\Delta t``) over
that one interval -- the standard variation-of-parameters solution:

```math
\\sigma_i = \\sigma_{i-1}e^{-\\Delta t/\\tau} + Er\\tau\\left(1-e^{-\\Delta t/\\tau}\\right)
```

When `strain_path` stops changing (``r=0``, the "fix it" moment), this collapses to the
textbook **stress relaxation** solution ``\\sigma(t)=\\sigma_0e^{-t/\\tau}`` exactly, to
machine precision, no matter how coarsely the hold is sampled -- see the self-check below.

Ignores plasticity -- a relaxation test is normally run below yield, so unlike
[`strain_from_stress_with_creep`](@ref) there is no return-mapping here, only an optional
hard `sigma_cap` (default `Inf`) that clips runaway stress if `strain_path` is pushed
unrealistically far. Returns `(stress_path, eps_c_path)`, where `eps_c_path =
strain_path .- stress_path./E` is whatever the elastic relation can't otherwise explain --
the creep residual, in the units [`lateral_strain_from_axial`](@ref) expects.
"""
function stress_from_strain_with_relaxation(strain_path, time_path, E, eta; sigma_cap=Inf)
    n = length(strain_path)
    stress = zeros(Float64, n)
    tau = axial_stress_relaxation_time(E, eta)
    stress[1] = clamp(E * strain_path[1], -sigma_cap, sigma_cap)
    for i in 2:n
        dt = max(time_path[i] - time_path[i-1], 0.0)
        rate = dt > 0 ? (strain_path[i] - strain_path[i-1]) / dt : 0.0
        x = dt / tau
        decay = exp(-x)
        # 1-exp(-x) loses precision catastrophically as x->0 (exactly the near-elastic,
        # huge-eta limit this function is supposed to reduce to exactly) -- expm1 is the
        # standard fix, accurate to machine precision for small x where a naive `1-decay`
        # was actually the dominant source of error, not the Euler-vs-exact integration
        # this function exists to avoid (confirmed directly: `1-decay` gave ~1e-4 relative
        # error per step here, `-expm1(-x)` brings it back to ~1e-13).
        one_minus_decay = -expm1(-x)
        stress[i] = clamp(stress[i-1] * decay + E * rate * tau * one_minus_decay, -sigma_cap, sigma_cap)
    end
    eps_c_path = strain_path .- stress ./ E
    return stress, eps_c_path
end

# ╔═╡ 1a52d2bb-c6ac-4e0d-8a00-39ec56e0cdff
"""
	lateral_strain_from_axial(strain_path, residual_path, nu)

The lateral (necking) strain that accompanies an axial strain history, given the combined
plastic+creep residual already computed by [`strain_from_stress_with_creep`](@ref). The
elastic part of the axial strain contracts laterally by the material's Poisson's ratio
`ν`; the *inelastic* part (`residual_path`, either plastic or viscous) is treated as
volume-preserving flow -- the standard reason real rods visibly neck once they yield or
creep, regardless of what their elastic `ν` happens to be.
"""
function lateral_strain_from_axial(strain_path, residual_path, nu)
    return [-nu * (strain_path[i] - residual_path[i]) - 0.5 * residual_path[i]
            for i in eachindex(strain_path)]
end

# ╔═╡ c1d10e8c-c3e6-4ea5-8dad-d83b87c397a8
"""
	deborah_regime(tau_wave_s, tau_M_s, T_s)

Classify a forcing period `T` against the two clocks from [`wave_relaxation_clocks`](@ref):
whether inertia matters (`T` comparable to `τ_wave`) and, if not, whether the material
responds elastically, transiently, or viscously over `T` -- the Deborah number
`De = τ_M/T`. Returns `(R_wave, De, regime)`, `regime` one of `"seismic wave"`, `"elastic
quasi-static"`, `"viscoelastic"`, `"viscous flow"`.

!!! note "Why two ratios, not one"
    `R_wave = τ_wave/T` alone decides whether inertia matters at all; only once it's small
    does `De = τ_M/T` get to decide elastic vs. viscoelastic vs. viscous. A material can
    have `De ≫ 1` (elastic) and still be in the wave regime if `T` is short enough that
    `R_wave` is not small -- that's exactly the seismic-wave case.
"""
function deborah_regime(tau_wave_s, tau_M_s, T_s)
    R_wave = tau_wave_s / T_s
    De = tau_M_s / T_s
    regime = if R_wave > 0.15
        "seismic wave"
    elseif De > 3
        "elastic quasi-static"
    elseif De > 0.3
        "viscoelastic"
    else
        "viscous flow"
    end
    return (R_wave=R_wave, De=De, regime=regime)
end

# ╔═╡ 9d2f94d6-6ff6-4f13-8a71-64eadf2c1a40
"""
	scripted_forcing_response(E_GPa, nu, sigma_y_MPa, H_MPa, eta_Pas, kind, amplitude_MPa, duration_s; n=300)

Drive [`strain_from_stress_with_creep`](@ref) with one of two realistic, non-periodic
forcing shapes instead of a sine wave: `kind=:ramp` linearly ramps stress from 0 to
`amplitude_MPa` over `duration_s` -- the interseismic-loading story, steady tectonic
loading with essentially no time for creep. `kind=:step` jumps instantly to
`amplitude_MPa` and holds it for `duration_s` -- the postseismic-relaxation /
mantle-convection story, an elastic jump followed by creep. Same GPa/MPa/Pa·s/s unit
convention as [`wave_relaxation_clocks`](@ref); this is the SAME underlying rheology as
the free-play pull-release, just driven by a scripted history instead of a mouse gesture,
so `duration_s` can be compared directly against the two clocks.

Returns `(t_s, stress_Pa, strain, lateral)`.
"""
function scripted_forcing_response(E_GPa, nu, sigma_y_MPa, H_MPa, eta_Pas, kind, amplitude_MPa, duration_s; n=300)
    t = collect(range(0, duration_s; length=n))
    stress = if kind === :ramp
        (amplitude_MPa * 1.0e6) .* (t ./ duration_s)
    else
        fill(amplitude_MPa * 1.0e6, n)
    end
    strain, residual = strain_from_stress_with_creep(
        stress, t, E_GPa * 1.0e9, sigma_y_MPa * 1.0e6, eta_Pas; H=H_MPa * 1.0e6)
    lateral = lateral_strain_from_axial(strain, residual, nu)
    return t, stress, strain, lateral
end

# ╔═╡ 89cd485a-c305-44c0-a54c-24941c7f4e28
"""
	lattice_homogenized_nu(k_edge, k_diag; nx=3, ny=3, k_anchor_ratio=1.0e-6)

The effective Poisson's ratio of a small 2-D spring lattice under uniaxial tension --
derived, not assumed. Builds an `nx`&times;`ny` grid of pin-jointed nodes connected by
orthogonal edge springs (stiffness `k_edge`) and, within each cell, a pair of diagonal
braces (stiffness `k_diag`); clamps the left edge, prescribes a small stretch on the
right edge, leaves every other degree of freedom completely free, and solves the linear
elastic equilibrium problem directly. `k_anchor_ratio·k_edge` is a tiny stabilizing spring
to ground added to every degree of freedom -- needed because a *purely* orthogonal
lattice (`k_diag=0`) has a genuine zero-energy mechanism (see the note) and is otherwise
singular; its effect vanishes as `k_anchor_ratio → 0`.

Returns `(nu_eff, u)`: the effective Poisson's ratio, and the full displacement field
(returned only so the self-check below can verify the solver directly).

!!! note "Why a plain grid has no mechanism to couple x and y at all"
    A spring aligned exactly along `x` only ever senses relative motion of its endpoints
    *along* `x` -- to linear order, sliding an endpoint sideways doesn't change the
    spring's length. So in a lattice built from purely horizontal and purely vertical
    springs, pulling in `x` cannot exert force in `y` anywhere -- not "a small amount",
    none. That's `nu_eff=0`; the diagonal braces are what create the coupling.
"""
function lattice_homogenized_nu(k_edge, k_diag; nx=3, ny=3, k_anchor_ratio=1.0e-6)
    nodeid(i, j) = (j - 1) * nx + i  # 1-based; i in 1:nx (columns), j in 1:ny (rows)
    pos(n) = (Float64((n - 1) % nx), Float64((n - 1) ÷ nx))
    ndof = 2 * nx * ny
    K = zeros(ndof, ndof)
    dofs_of(n) = (2n - 1, 2n)
    function add_spring!(n1, n2, k)
        k == 0 && return
        x1, y1 = pos(n1)
        x2, y2 = pos(n2)
        dx, dy = x2 - x1, y2 - y1
        L = hypot(dx, dy)
        c, s = dx / L, dy / L
        kloc = k .* [c*c c*s -c*c -c*s; c*s s*s -c*s -s*s; -c*c -c*s c*c c*s; -c*s -s*s c*s s*s]
        dofs = [dofs_of(n1)..., dofs_of(n2)...]
        K[dofs, dofs] .+= kloc
    end
    for j in 1:ny, i in 1:nx
        n = nodeid(i, j)
        i < nx && add_spring!(n, nodeid(i + 1, j), k_edge)
        j < ny && add_spring!(n, nodeid(i, j + 1), k_edge)
        if i < nx && j < ny
            add_spring!(n, nodeid(i + 1, j + 1), k_diag)
            add_spring!(nodeid(i + 1, j), nodeid(i, j + 1), k_diag)
        end
    end
    k_anchor = k_anchor_ratio * k_edge
    for d in 1:ndof
        K[d, d] += k_anchor
    end

    delta = 0.05
    fixed_dofs = Int[]
    fixed_vals = Float64[]
    for j in 1:ny
        nL = nodeid(1, j)
        push!(fixed_dofs, 2nL - 1, 2nL)
        push!(fixed_vals, 0.0, 0.0)
        nR = nodeid(nx, j)
        push!(fixed_dofs, 2nR - 1)
        push!(fixed_vals, delta)
    end
    free_dofs = setdiff(1:ndof, fixed_dofs)
    Kff = K[free_dofs, free_dofs]
    Kfc = K[free_dofs, fixed_dofs]
    uf = -(Kff \ (Kfc * fixed_vals))

    u = zeros(ndof)
    u[fixed_dofs] = fixed_vals
    u[free_dofs] = uf

    uy_top = u[2*nodeid(nx, ny)]
    uy_bot = u[2*nodeid(nx, 1)]
    eps_x = delta / (nx - 1)
    eps_y = (uy_top - uy_bot) / (ny - 1)
    nu_eff = -eps_y / eps_x
    return nu_eff, u
end

# ╔═╡ 95f0ad68-ee90-42b2-91a3-f2788271964f
let
    nx_t, ny_t = 3, 3
    nodeid_t(i, j) = (j - 1) * nx_t + i
    nu_plain_t, u_plain_t = lattice_homogenized_nu(1.0, 0.0; nx=nx_t, ny=ny_t)
    nu_light_t, u_light_t = lattice_homogenized_nu(1.0, 0.4; nx=nx_t, ny=ny_t)
    nu_full_t, u_full_t = lattice_homogenized_nu(1.0, 1.0; nx=nx_t, ny=ny_t)

    # for k_diag=0 there is no mechanism coupling x and y at all -- every uy should be
    # zero, not just nu_eff
    max_uy_plain = maximum(abs, u_plain_t[2:2:end])
    @assert max_uy_plain < 1.0e-4

    # antisymmetry about the middle row: an independent check on the solver itself,
    # unrelated to the k_diag=0 special case
    max_asym = maximum(abs(u_full_t[2*nodeid_t(i, 1)] + u_full_t[2*nodeid_t(i, ny_t)]) for i in 1:nx_t)
    mid_uy = maximum(abs(u_full_t[2*nodeid_t(i, 2)]) for i in 1:nx_t)
    @assert max_asym < 1.0e-8
    @assert mid_uy < 1.0e-8
    @assert nu_light_t > 0 && nu_full_t > 0
    @assert nu_full_t > nu_light_t  # more bracing, more coupling

    md"""
    !!! correct "Self-check"
        A plain orthogonal grid has genuinely zero lateral displacement everywhere
        (max ``|u_y|=$(round(max_uy_plain, sigdigits=2))``, not just a small number) when
        pulled --- ``\nu=0`` isn't assumed, it falls out because nothing in the lattice
        can move a node in ``y``. The solver also respects the top-bottom antisymmetry
        the geometry demands (asymmetry ``$(round(max_asym, sigdigits=2))``, middle-row
        motion ``$(round(mid_uy, sigdigits=2))``, both numerical noise). Adding diagonal
        braces creates real coupling: ``\nu=$(round(nu_light_t, digits=3))`` lightly
        braced, ``\nu=$(round(nu_full_t, digits=3))`` fully braced --- more bracing, more
        coupling, exactly the Material A vs. B story.
    """
end

# ╔═╡ 7d5499f8-dfeb-489b-8cfc-bd2d10886bf7
begin
    # three genuinely braced (nonzero-diagonal) presets for the widget itself -- the pure
    # k_diag=0 case is the clean textbook limit used in the self-check above to validate
    # the solver, but it isn't a material anyone could hand you, so it isn't offered as
    # one of the three scrubbable options
    nu_light, _ = lattice_homogenized_nu(1.0, 0.2)
    nu_medium, _ = lattice_homogenized_nu(1.0, 0.6)
    nu_full, _ = lattice_homogenized_nu(1.0, 1.2)
    (light=round(nu_light, digits=4), medium=round(nu_medium, digits=4), full=round(nu_full, digits=4))
end

# ╔═╡ 1f8e58a0-c7a3-44f5-99b0-d5bbb9ee89b5
"""
	effective_modulus(k, dx, A)

The intensive (material) elastic modulus implied by a 1-D chain of identical springs of
stiffness `k` spaced `dx` apart, forming a rod of cross-sectional area `A`. Each of `n`
springs in series carries the same force `F` and stretches `F/k`, so the total elongation
is `n·F/k` over a length `L = n·dx`, giving `stress/strain = (F/A) / (n·F/(k·dx·n)) =
k·dx/A` -- independent of `n`. This is why the raw spring constant `k` is not itself a
material property, but `k·dx/A` is.
"""
effective_modulus(k, dx, A) = k * dx / A

# ╔═╡ 0c92d271-1808-4109-b173-5dc4f52fa198
let
    E_test, sy_test = 100.0, 10.0

    strains_a = collect(range(0, 0.05; length=50))
    stress_a, _ = elastoplastic_response(strains_a, E_test, sy_test)
    err_a = maximum(abs.(stress_a .- E_test .* strains_a))

    path_b = collect(vcat(range(0, 0.3; length=200), range(0.3, -0.05; length=200)))
    stress_b, _ = elastoplastic_response(path_b, E_test, sy_test)
    slope_b = (stress_b[250] - stress_b[210]) / (path_b[250] - path_b[210])
    idx0 = 200 + findfirst(<=(0), stress_b[201:end])
    eps_resid_numeric = path_b[idx0]
    eps_resid_theory = 0.3 - sy_test / E_test

    k_test, dx_test, A_test = 5.0, 0.1, 2.0
    moduli_direct = [begin
        F = 3.0
        strain = (ns * F / k_test) / (ns * dx_test)
        (F / A_test) / strain
    end for ns in (5, 10, 50, 200)]

    @assert isapprox(err_a, 0.0; atol=1e-9)
    @assert isapprox(slope_b, E_test; atol=1e-6)
    @assert isapprox(eps_resid_numeric, eps_resid_theory; atol=5e-3)
    @assert all(isapprox.(moduli_direct, effective_modulus(k_test, dx_test, A_test); atol=1e-9))

    # strain_from_stress is the force-controlled dual: below yield it must invert
    # elastoplastic_response exactly (round-trip strain -> stress -> strain).
    strains_c = collect(range(-0.05, 0.05; length=40))
    stress_c, _ = elastoplastic_response(strains_c, E_test, sy_test)
    strain_back, _ = strain_from_stress(stress_c, E_test, sy_test)
    err_c = maximum(abs.(strain_back .- strains_c))
    @assert isapprox(err_c, 0.0; atol=1e-9)

    # with H=0, no amount of applied stress beyond sigma_y0 should be achievable --
    # strain_from_stress must clamp the resulting stress at sigma_y0, never exceed it.
    strain_runaway, _ = strain_from_stress([50 * sy_test], E_test, sy_test; H=0.0)
    stress_runaway, _ = elastoplastic_response(strain_runaway, E_test, sy_test; H=0.0)
    @assert isapprox(stress_runaway[1], sy_test; atol=1e-6)

    # the check above alone does NOT prove plastic strain actually accumulated -- a version
    # that silently clamped stress at sigma_y0 while eps_p stayed 0 forever passes it too
    # (that was the actual bug: sig got clamped to exactly sigma_y0 before the yield check
    # ever ran, so eps_p never moved). Confirm directly: eps_p must be nonzero once pushed
    # past yield, and it must grow as the overstress grows (not just saturate instantly).
    _, epsp_runaway = strain_from_stress([50 * sy_test], E_test, sy_test; H=0.0)
    _, epsp_mild = strain_from_stress([1.2 * sy_test], E_test, sy_test; H=0.0)
    @assert epsp_runaway[1] > 0
    @assert epsp_mild[1] > 0
    @assert epsp_runaway[1] > epsp_mild[1]

    # strain_from_stress_with_creep must reduce to strain_from_stress exactly when eta is
    # huge (creep negligible over any realistic drag duration).
    stress_d = fill(5.0, 20)
    time_d = collect(range(0, 2.0; length=20))
    strain_noc, _ = strain_from_stress(stress_d, E_test, sy_test)
    strain_hugeeta, _ = strain_from_stress_with_creep(stress_d, time_d, E_test, sy_test, 1.0e12)
    err_d = maximum(abs.(strain_hugeeta .- strain_noc))
    @assert isapprox(err_d, 0.0; atol=1e-8)

    # a held constant stress below yield should creep linearly in time: eps_creep =
    # sigma*t/(3*eta) -- the Trouton-ratio-corrected analytic Maxwell creep solution for
    # a uniaxial step load (see axial_maxwell_creep_rate).
    eta_test = 200.0
    sigma_held = 3.0  # below sy_test=10, so no plastic contribution -- isolates the creep term
    T_hold = 5.0
    stress_e = fill(sigma_held, 100)
    time_e = collect(range(0, T_hold; length=100))
    strain_e, residual_e = strain_from_stress_with_creep(stress_e, time_e, E_test, sy_test, eta_test)
    creep_theory = sigma_held * T_hold / (3 * eta_test)
    @assert isapprox(residual_e[end], creep_theory; atol=1e-3)
    @assert isapprox(strain_e[end], sigma_held / E_test + creep_theory; atol=1e-3)

    # lateral_strain_from_axial: purely elastic (residual=0) must give exactly -nu*strain;
    # fully inelastic (residual==strain, no elastic part left) must give exactly -0.5*strain,
    # the incompressible limit -- independent of nu, since nothing elastic remains.
    nu_test = 0.28
    strain_f = collect(range(-0.05, 0.05; length=30))
    zero_residual = zeros(length(strain_f))
    lateral_elastic = lateral_strain_from_axial(strain_f, zero_residual, nu_test)
    err_lat_elastic = maximum(abs.(lateral_elastic .- (-nu_test .* strain_f)))
    @assert isapprox(err_lat_elastic, 0.0; atol=1e-12)

    lateral_incompressible = lateral_strain_from_axial(strain_f, strain_f, nu_test)
    err_lat_incompr = maximum(abs.(lateral_incompressible .- (-0.5 .* strain_f)))
    @assert isapprox(err_lat_incompr, 0.0; atol=1e-12)

    md"""
    !!! correct "Self-check"
        Below yield the response is linear to within ``$(round(err_a, sigdigits=3))`` ---
        numerical noise, not curvature. Unloading after yield retraces at slope
        ``$(round(slope_b, digits=2))`` (theory ``E=$(E_test)``), leaving a residual
        strain of ``$(round(eps_resid_numeric, digits=4))`` at zero stress (theory
        ``$(round(eps_resid_theory, digits=4))``). And the modulus ``k \cdot dx / A``,
        computed directly from chains of 5, 10, 50, and 200 springs in series, comes out
        to ``$(round.(moduli_direct, digits=6))`` in every case --- unchanged, confirming
        it is independent of how finely the same rod is subdivided. The force-controlled
        dual, `strain_from_stress`, round-trips the displacement-controlled response to
        within ``$(round(err_c, sigdigits=3))`` below yield, and correctly refuses to let
        applied stress exceed ``\sigma_{y0}`` when ``H=0`` (it clamps at exactly
        ``$(round(stress_runaway[1], digits=6))``, not the ``50\sigma_{y0}`` that was
        asked for) -- perfect plasticity under load control has no equilibrium beyond
        yield, and the widget respects that instead of returning nonsense. Adding the
        creep term changes nothing (``$(round(err_d, sigdigits=3))`` difference) when
        ``\eta`` is huge, and holding a below-yield stress of ``$(sigma_held)`` for
        ``$(T_hold)`` s at ``\eta=$(eta_test)`` creeps by ``$(round(residual_e[end], digits=4))``
        --- matching the Trouton-corrected Maxwell solution ``\sigma t/(3\eta)`` =
        $(round(creep_theory, digits=4)) exactly (not the uncorrected ``\sigma t/\eta``
        a plain 1-D reading of the dashpot law would suggest): the same force, held longer,
        deforms more, with no yielding involved at all. And the lateral (necking) strain
        matches ``-\nu\varepsilon`` exactly while purely elastic, then switches to the
        incompressible ``-0.5\varepsilon`` exactly once nothing elastic is left ---
        necking is a property of the *flow*, not of ``\nu``.
    """
end

# ╔═╡ e76b20bd-28f7-418c-8fb3-42a40e710497
md"## Two Independent Constants"

# ╔═╡ 02d7acba-4408-4158-b67d-76516d239309
begin
    """
    	E_from_KG(K, G)

    Young's modulus implied by bulk modulus `K` and shear modulus `G` (isotropic
    elasticity): `E = 9KG / (3K + G)`.
    """
    E_from_KG(K, G) = 9 * K * G / (3 * K + G)

    """
    	nu_from_KG(K, G)

    Poisson's ratio implied by `K`, `G`: `\\nu = (3K - 2G) / (2(3K + G))`. Reduces to
    exactly `0.5` when `G=0` -- the shear-free, incompressible-shape limit of a fluid.
    """
    nu_from_KG(K, G) = (3 * K - 2 * G) / (2 * (3 * K + G))

    """
    	K_from_Enu(E, nu)

    Bulk modulus implied by Young's modulus `E` and Poisson's ratio `\\nu` -- the inverse
    of `E_from_KG`/`nu_from_KG`.
    """
    K_from_Enu(E, nu) = E / (3 * (1 - 2 * nu))

    """
    	G_from_Enu(E, nu)

    Shear modulus implied by `E`, `\\nu` -- the inverse of `E_from_KG`.
    """
    G_from_Enu(E, nu) = E / (2 * (1 + nu))
end

# ╔═╡ 3d8b6f21-9e4a-4c7d-b502-1f6a9d3e8c05
let
    sigma_t, eta_t = 12.0e6, 4.0e19

    # the deviatoric decomposition this function is derived from, spelled out directly:
    # uniaxial sigma with zero lateral stress -> mean stress p, deviatoric s_ij
    p_t = sigma_t / 3
    s11_t, s22_t, s33_t = sigma_t - p_t, -p_t, -p_t
    @assert isapprox(s11_t + s22_t + s33_t, 0.0; atol=1e-6)   # deviator is traceless, always
    @assert isapprox(s11_t, 2 * sigma_t / 3; atol=1e-6)

    # the dashpot law e_dot = s/(2*eta) applied to JUST the deviatoric part must reproduce
    # this function exactly
    edot_direct = s11_t / (2 * eta_t)
    @assert isapprox(edot_direct, axial_maxwell_creep_rate(sigma_t, eta_t); rtol=1e-12)

    # the elastic identity that makes E the right (unmodified) coefficient for the RATE
    # term while eta needs the factor of 3: 1/(3mu) + 1/(9K) must equal 1/E exactly, for
    # any K, mu -- this is what lets the widget keep a single E for elastic loading while
    # still needing 3eta, not eta, for creep
    K_t, mu_t = 50.0e9, 30.0e9
    E_t = E_from_KG(K_t, mu_t)
    lhs = 1 / (3 * mu_t) + 1 / (9 * K_t)
    @assert isapprox(lhs, 1 / E_t; rtol=1e-10)

    md"""
    !!! correct "Self-check"
        A uniaxial stress of $(sigma_t/1e6) MPa splits into a deviatoric part
        $(round(s11_t/1e6, digits=2)) MPa --- exactly ``2/3`` of the total, confirmed
        directly from the (traceless, by construction) stress deviator --- and a lossless
        volumetric third. Feeding only that deviatoric part into the dashpot's own law
        ``s/(2\eta)`` reproduces `axial_maxwell_creep_rate` exactly, and the same
        decomposition's elastic half confirms ``1/(3\mu)+1/(9K)=1/E`` to
        $(round(abs(lhs*E_t-1), sigdigits=2)) relative error --- exactly why the widget
        can keep a single, unmodified ``E`` for the elastic slope while still needing
        ``3\eta``, not ``\eta``, for creep.
    """
end

# ╔═╡ 3cf3b5c7-f651-4a0d-842e-ba54b5bb4c59
"""
	wave_relaxation_clocks(E_GPa, nu, eta_Pas, L_km, rho)

The two clocks that decide whether a material communicates stress as a wave or relaxes it
viscously: the wave-crossing time `τ_wave = L/Vs` (`Vs = sqrt(μ/ρ)`, reusing
[`G_from_Enu`](@ref) for `μ`) and the Maxwell relaxation time `τ_M = η/μ`. Takes the same
human-scale units (`E` in GPa, `L` in km) and converts to SI internally so the returned
clocks are in seconds -- real Earth numbers, not schematic ones.

Returns `(tau_wave_s, tau_M_s)`.
"""
function wave_relaxation_clocks(E_GPa, nu, eta_Pas, L_km, rho)
    mu_Pa = G_from_Enu(E_GPa, nu) * 1.0e9
    Vs = sqrt(mu_Pa / rho)
    return (L_km * 1.0e3) / Vs, eta_Pas / mu_Pa
end

# ╔═╡ c9e5a1b3-6d4f-4a8c-8b3e-5f7e9d1a3c24
let
    E_test2, eta_test2 = 300.0, 500.0  # tau_relax = 3*eta/E = 5 (arbitrary consistent units)
    tau_test2 = axial_stress_relaxation_time(E_test2, eta_test2)
    @assert isapprox(tau_test2, 5.0; atol=1e-12)

    # reduces to plain elasticity when eta is huge (no creep to speak of): a ramping strain
    # history should give stress = E*strain exactly at every sample, same style of check as
    # strain_from_stress_with_creep's own huge-eta reduction above.
    strain_ramp = collect(range(0, 0.01; length=50))
    time_ramp = collect(range(0, 1.0; length=50))
    stress_elastic, _ = stress_from_strain_with_relaxation(strain_ramp, time_ramp, E_test2, 1.0e12)
    err_elastic = maximum(abs.(stress_elastic .- E_test2 .* strain_ramp))
    @assert isapprox(err_elastic, 0.0; atol=1e-9)

    # the classic relaxation test: strain held flat at eps0 from t=0 onward -- exactly the
    # "fix it" moment the widget's hold gesture represents -- sampled at only a handful of
    # points across 6*tau_relax, deliberately coarse, since the whole point of integrating
    # each interval exactly (not with forward-Euler) is that it shouldn't matter.
    eps0 = 0.02
    n_hold = 6
    T_hold = 6 * tau_test2
    time_relax = collect(range(0, T_hold; length=n_hold))
    strain_relax = fill(eps0, n_hold)
    stress_relax, eps_c_relax = stress_from_strain_with_relaxation(strain_relax, time_relax, E_test2, eta_test2)
    sigma0 = E_test2 * eps0  # stress[1] == E*strain_path[1] exactly -- the function's own first line
    theory = sigma0 .* exp.(-time_relax ./ tau_test2)
    err_relax = maximum(abs.(stress_relax .- theory)) / sigma0
    eps_c_theory = eps0 * (1 - exp(-T_hold / tau_test2))  # everything not-yet-relaxed strain has become creep

    @assert isapprox(err_relax, 0.0; atol=1e-12)
    @assert stress_relax[end] < 0.01 * sigma0
    @assert isapprox(eps_c_relax[end], eps_c_theory; atol=1e-10)

    # independent validation: a crude explicit-Euler integration of the SAME ODE on a very
    # fine grid (not the function under test), as ground truth, versus our function
    # evaluated at only a handful of points on the identical strain-vs-time path -- if the
    # exact per-interval formula is right, evaluating it coarsely should still land on the
    # fine-grid answer, which an explicit-Euler integrator at that same coarse spacing
    # could not. The one thing the coarse grid must still do is place a sample exactly at
    # the ramp's own kink (t_ramp_b): "exact per interval" means exact GIVEN that the true
    # strain rate really is constant across each sampled interval, and a coarse interval
    # straddling an unresolved kink violates that premise by construction -- that would be
    # testing whether 8 points can resolve a kink they were never given the location of,
    # not whether the integrator is correct.
    eps0_b, t_ramp_b, t_hold_b = 0.015, 2.0, 25.0
    strain_of_t(t) = t <= t_ramp_b ? eps0_b * t / t_ramp_b : eps0_b
    n_fine = 200_000
    t_fine = collect(range(0, t_ramp_b + t_hold_b; length=n_fine))
    sigma_fine = zeros(n_fine)
    sigma_fine[1] = E_test2 * strain_of_t(t_fine[1])
    for i in 2:n_fine
        dtf = t_fine[i] - t_fine[i-1]
        rate = (strain_of_t(t_fine[i]) - strain_of_t(t_fine[i-1])) / dtf
        sigma_fine[i] = sigma_fine[i-1] + dtf * (E_test2 * rate - sigma_fine[i-1] / tau_test2)
    end
    n_coarse = 8
    t_coarse = vcat([0.0, t_ramp_b], collect(range(t_ramp_b, t_ramp_b + t_hold_b; length=n_coarse-1))[2:end])
    strain_coarse = strain_of_t.(t_coarse)
    stress_coarse, _ = stress_from_strain_with_relaxation(strain_coarse, t_coarse, E_test2, eta_test2)
    err_coarse_vs_fine = abs(stress_coarse[end] - sigma_fine[end]) / (E_test2 * eps0_b)
    @assert err_coarse_vs_fine < 1e-5

    # ties the two clocks together: the axial relaxation time and the shear-based tau_M
    # from wave_relaxation_clocks must differ by exactly 3/(2(1+nu)), the same
    # Trouton-ratio bookkeeping applied to a different clock
    nu_check = 0.25
    E_check_GPa, eta_check_Pas = 75.0, 1.0e19
    _, tauM_check = wave_relaxation_clocks(E_check_GPa, nu_check, eta_check_Pas, 15.0, 2700.0)
    tau_relax_check = axial_stress_relaxation_time(E_check_GPa * 1.0e9, eta_check_Pas)
    ratio_check = tau_relax_check / tauM_check
    ratio_theory = 3 / (2 * (1 + nu_check))
    @assert isapprox(ratio_check, ratio_theory; rtol=1e-9)

    md"""
    !!! correct "Self-check"
        With ``\eta`` huge, `stress_from_strain_with_relaxation` reduces to plain
        elasticity, ``\sigma=E\varepsilon``, exactly ($(round(err_elastic, sigdigits=3))
        error). The real test: hold a strain of $(eps0) fixed for ``6\tau_{\rm relax}``,
        sampled at only $(n_hold) points -- deliberately coarse, since integrating each
        interval exactly (rather than with forward-Euler) should make that not matter.
        Stress matches the textbook ``\sigma_0e^{-t/\tau_{\rm relax}}`` to
        $(round(err_relax, sigdigits=3)) relative error (machine precision, as it should
        be -- the strain never changes during the hold, so the "exact per interval"
        formula's own approximation, that the strain rate is constant over each interval,
        is exactly true here, not merely a good approximation), decaying to
        $(round(stress_relax[end]/sigma0*100, sigdigits=2))% of its starting value --
        essentially fully relaxed, with essentially all of the held strain now creep
        (``\varepsilon_c\to`` $(round(eps_c_relax[end], sigdigits=4)), matching
        ``\varepsilon_0(1-e^{-6})=`` $(round(eps_c_theory, sigdigits=4)) to
        $(round(abs(eps_c_relax[end]/eps_c_theory-1)*100, sigdigits=2))%): with nothing
        left to relax, there is almost nothing left to be elastic either. Independently,
        evaluating a *different*, genuinely ramp-then-hold strain path at only
        $(n_coarse) points lands within $(round(err_coarse_vs_fine*100, sigdigits=2))%
        of a $(n_fine)-step explicit-Euler integration of the identical ODE -- exactly the
        coarse-sampling robustness the exact-per-interval update was built for, confirmed
        against an independent (if crude) numerical reference rather than just against
        itself. And the two relaxation clocks agree exactly where theory says they should:
        ``\tau_{\rm relax}/\tau_M=`` $(round(ratio_check, digits=4)), matching
        ``3/(2(1+\nu))=`` $(round(ratio_theory, digits=4)) to machine precision -- the same
        dashpot, the same Trouton-ratio bookkeeping, just timing a uniaxial pull instead of
        a pure shear.
    """
end

# ╔═╡ 6a1af493-35f5-4813-8940-1992c72cb1e5
let
    # crust-like material, real numbers (K=50, G=30 GPa)
    E_t, nu_t, rho_t = E_from_KG(50.0, 30.0), nu_from_KG(50.0, 30.0), 2700.0
    eta_elastic_t = 1.0e24  # near-elastic -- eta this large should make tau_M enormous
    tw_t, tm_t = wave_relaxation_clocks(E_t, nu_t, eta_elastic_t, 15.0, rho_t)

    yr = 365.25 * 24 * 3600
    # four T values, each picked to sit squarely inside one of the four named regimes
    r_wave = deborah_regime(tw_t, tm_t, 0.5 * tw_t)
    r_elastic = deborah_regime(tw_t, tm_t, 30.0 * yr)
    r_visco = deborah_regime(tw_t, 1.0 * yr, 1.0 * yr)          # a much smaller tau_M for this check
    r_viscous = deborah_regime(tw_t, 270.0 * yr, 50.0e6 * yr)

    @assert r_wave.regime == "seismic wave"
    @assert r_elastic.regime == "elastic quasi-static"
    @assert r_visco.regime == "viscoelastic"
    @assert r_viscous.regime == "viscous flow"
    @assert tw_t > 0 && tw_t < 60  # a crustal S-wave crossing a ~15 km fault: seconds, not years
    @assert tm_t > 1.0e10  # near-elastic: tau_M should dwarf any human timescale

    # elastic limit of scripted_forcing_response: with eta huge, a :ramp must match the
    # pure elastic prediction sigma(t)/E almost exactly (a little real creep remains --
    # eta=1e30 here is very large but not literally infinite)
    t_r, stress_r, strain_r, _ = scripted_forcing_response(E_t, nu_t, 100.0, 0.0, 1.0e30, :ramp, 60.0, 100.0)
    strain_theory_r = stress_r ./ (E_t * 1.0e9)
    err_r = maximum(abs.(strain_r .- strain_theory_r)) / maximum(abs.(strain_theory_r))

    # a :step at finite eta must show BOTH pieces of the story: the very first sample is
    # the pure elastic jump (no creep has had time to accumulate yet), and the last sample
    # must match the Trouton-corrected analytic sum sigma/E + sigma/(3*mu) exactly, one
    # Maxwell time later (duration = tau_M for THIS material, see mu_asth below)
    mu_asth = G_from_Enu(120.0, nu_t) * 1.0e9
    _, tau_M_asth = wave_relaxation_clocks(120.0, nu_t, 5.0e18, 50.0, 3300.0)
    t_s, stress_s, strain_s, _ = scripted_forcing_response(120.0, nu_t, 200.0, 0.0, 5.0e18, :step, 20.0, tau_M_asth)
    jump_theory = 20.0e6 / (120.0e9)
    creep_theory_s = 20.0e6 / (3 * mu_asth)

    @assert isapprox(err_r, 0.0; atol=1e-3)
    @assert isapprox(strain_s[1], jump_theory; rtol=1e-6)
    @assert isapprox(strain_s[end], jump_theory + creep_theory_s; rtol=1e-3)
    @assert strain_s[end] > 1.5 * strain_s[1]  # creep has clearly overtaken the initial jump

    md"""
    !!! correct "Self-check"
        A ``15`` km crustal fault, real ``V_s`` and ``\eta`` plugged into
        `wave_relaxation_clocks`, gives ``\tau_{\rm wave}\approx`` $(round(tw_t, digits=2))
        s --- an actual seismic-wave crossing time, not a schematic number. Four forcing
        periods, picked relative to that same material's two clocks, land in exactly the
        four possible regimes: seismic wave, elastic quasi-static, viscoelastic, viscous
        flow --- nothing assumed, `deborah_regime` derives the label from
        ``\tau_{\rm wave}/T`` and ``\tau_M/T`` alone (only the last three appear as widget
        presets; the wave regime needs real inertia this quasi-static engine doesn't have).
        `scripted_forcing_response`'s `:ramp` matches the pure elastic prediction to within
        $(round(err_r*100, sigdigits=2))% once ``\eta`` is large enough that creep is
        negligible, and its `:step` correctly separates the instantaneous elastic jump
        ( $(round(strain_s[1], sigdigits=3)), matching ``\sigma/E`` exactly) from the
        creep that follows: $(round(strain_s[end], sigdigits=3)) after one Maxwell time,
        matching the Trouton-corrected sum ``\sigma/E+\sigma/(3\mu)`` =
        $(round(jump_theory+creep_theory_s, sigdigits=3)) to within
        $(round(abs(strain_s[end]/(jump_theory+creep_theory_s)-1)*100, sigdigits=2))% ---
        not the ``2\times`` overshoot a naive ``\sigma/\eta``
        dashpot would have given --- the same rheology as free-play's pull-release, just
        driven by a scripted history instead of a mouse gesture.
    """
end

# ╔═╡ 6cbc6bab-3e0a-415a-acd0-5ea78d7e0f36
begin
    # three realistic regime bundles for Widget A's preset buttons -- material + geometry +
    # a SCRIPTED forcing (a ramp or a step-and-hold, not a periodic wave) whose kind and
    # duration are chosen relative to THAT material's own two clocks, so each preset
    # genuinely lands in its named regime rather than assuming round numbers happen to
    # work. nu comes from the SAME structurally-derived bracing levels the scrub already
    # offers (nu_light/medium/full, from the "7d5499f8" cell above) -- not a free-floating
    # number. There is deliberately no "seismic wave" / inertial preset: this engine is
    # quasi-static throughout (no mass), so a wave-crossing regime would need a genuinely
    # different solver to demonstrate honestly -- tau_wave is still computed and shown for
    # comparison, it just never gets its own scripted scenario.
    local yr = 365.25 * 24 * 3600
    local E_crust = E_from_KG(50.0, 30.0)
    local E_asth = 120.0
    local _, tm_asth = wave_relaxation_clocks(E_asth, nu_medium, 5.0e18, 50.0, 3300.0)

    regime_presets = (
        elastic=(name="Elastic loading", E=E_crust, nu=nu_medium, bracingIdx=1, eta=1.0e24, rho=2700.0,
            L=15.0, sigma_y=500.0, H=0.0, kind=:ramp, amplitude=30.0, duration=1.0 * yr),
        interseismic=(name="Interseismic loading", E=E_crust, nu=nu_medium, bracingIdx=1, eta=1.0e24, rho=2700.0,
            L=15.0, sigma_y=100.0, H=0.0, kind=:ramp, amplitude=60.0, duration=100.0 * yr),
        postseismic=(name="Postseismic relaxation", E=E_asth, nu=nu_medium, bracingIdx=1, eta=5.0e18, rho=3300.0,
            L=50.0, sigma_y=200.0, H=0.0, kind=:step, amplitude=20.0, duration=tm_asth),
        convection=(name="Mantle convection", E=300.0, nu=nu_full, bracingIdx=2, eta=1.0e21, rho=4500.0,
            L=2900.0, sigma_y=50.0, H=0.0, kind=:step, amplitude=20.0, duration=50_000.0 * yr),
    )
    regime_presets
end

# ╔═╡ bdbf5dc7-ecd0-4a11-a3dd-a6d0a597cc72
begin
    """
    	dilatation_strain(p, K)

    The isotropic (dilatational) strain per axis produced by isotropic tension `p` in a
    material of bulk modulus `K`: `e = p / (3K)`. Pure volume change, no shape change.
    """
    dilatation_strain(p, K) = p / (3 * K)

    """
    	shear_strain(tau, G)

    Engineering shear strain produced by shear stress `tau` in a material of shear
    modulus `G`: `\\gamma = \\tau / G`. Undefined for a fluid, `G=0` -- callers must
    guard this case.
    """
    shear_strain(tau, G) = tau / G

    """
    	uniaxial_strains(sigma, E, nu)

    Axial and lateral strain produced by uniaxial stress `sigma` in a material of Young's
    modulus `E` and Poisson's ratio `nu`: `(e_axial, e_lateral) = (sigma/E, -nu·sigma/E)`.
    Undefined for a fluid, `E=0` -- callers must guard this case.
    """
    uniaxial_strains(sigma, E, nu) = (sigma / E, -nu * sigma / E)

    """
    	vp(K, G, rho)

    P-wave speed for an isotropic elastic medium: `Vp = sqrt((K + 4G/3) / rho)`.
    """
    vp(K, G, rho) = sqrt((K + 4G / 3) / rho)

    """
    	vs(K, G, rho)

    S-wave speed for an isotropic elastic medium: `Vs = sqrt(G / rho)`. Exactly zero for
    a fluid (`G=0`) -- a fluid cannot carry a shear wave.
    """
    vs(K, G, rho) = sqrt(G / rho)
end

# ╔═╡ c104eee8-b783-47bf-8a56-7769a77fe426
let
    K_test, G_test, rho_test = 50.0e9, 30.0e9, 2700.0
    E_test = E_from_KG(K_test, G_test)
    nu_test = nu_from_KG(K_test, G_test)

    p_ref, tau_ref, sigma_ref = 5.0e9, 5.0e9, 5.0e9
    e_dil = dilatation_strain(p_ref, K_test)
    K_meas = p_ref / (3 * e_dil)

    gamma = shear_strain(tau_ref, G_test)
    G_meas = tau_ref / gamma

    e_ax, e_lat = uniaxial_strains(sigma_ref, E_test, nu_test)
    E_meas = sigma_ref / e_ax
    nu_meas = -e_lat / e_ax

    @assert isapprox(K_meas, K_test; rtol=1e-8)
    @assert isapprox(G_meas, G_test; rtol=1e-8)
    @assert isapprox(E_meas, E_test; rtol=1e-8)
    @assert isapprox(nu_meas, nu_test; rtol=1e-8)

    vpv, vsv = vp(K_test, G_test, rho_test), vs(K_test, G_test, rho_test)
    ratio = vpv / vsv
    ratio_theory = sqrt(2 * (1 - nu_test) / (1 - 2 * nu_test))
    @assert isapprox(ratio, ratio_theory; rtol=1e-8)
    @assert vs(K_test, 0.0, rho_test) == 0.0

    md"""
    !!! correct "Self-check"
        Three independent experiments -- pure dilatation, pure shear, and a uniaxial pull
        -- were each asked to report back the modulus that drives them. All three agree
        with the driving ``K=$(round(K_test/1e9, digits=1))\,\text{GPa}``,
        ``G=$(round(G_test/1e9, digits=1))\,\text{GPa}`` to 8 significant figures: no
        fourth experiment could reveal a third independent number. And
        ``V_p/V_s=$(round(ratio, digits=4))`` matches the ``\nu``-only formula exactly --
        the two wave speeds encode nothing beyond the same two constants.
    """
end

# ╔═╡ 00961e05-3259-494e-8123-38675b04dfcc
md"## The Interactive Widget"

# ╔═╡ 58011730-0491-400a-83a1-ebaf42c6b574
begin
    """
    	ElastoplasticLoadingInput(; E=75.0, sigma_y=100.0, H=0.0, eta=1.0e24, nu=0.2412)

    Initial state for the pull-release-repeat widget: elastic modulus `E` (GPa), yield
    stress `sigma_y` and hardening modulus `H` (MPa), viscosity `eta` (Pa·s), Poisson's
    ratio `nu`, density `rho` (kg/m³) and length scale `L` (km) -- all real units, used
    both for manual dragging and for the preset-triggered scripted scenarios. There is
    deliberately only ONE parameter set now: manual pulling at a real `η` simply shows no
    visible creep over a few seconds of dragging (correctly -- creep needs `τ_M`, which is
    years to millions of years), and that is the whole point of the preset animations,
    which fast-forward through the real, relevant duration instead. `presetKey` tracks
    which of the Appendix's `regime_presets` the scripted-scenario panel should compute.
    The student's drag history is tracked entirely client-side and read back through this
    same bond.

    !!! note "Why literal defaults, not `nu_full`/`regime_presets` directly"
        Every default here is a plain number, deliberately not a reference to another
        notebook variable. A `@bind`-constructing cell's dependency on values used only
        inside ITS callee's default-argument expressions isn't reliably tracked by Pluto's
        static analysis -- confirmed directly: referencing `nu_full`/`regime_presets` here
        threw `UndefVarError` on a fresh client connecting to a cold-started kernel, even
        though `pluto-collab restart` reported zero errors (that check doesn't reproduce a
        new client's own bind-evaluation path). The fix is the same wherever this pattern
        shows up in a Pluto notebook -- literal defaults only.
    """
    struct ElastoplasticLoadingInput
        E::Float64
        sigma_y::Float64
        H::Float64
        eta::Float64
        nu::Float64
        rho::Float64
        L::Float64
        presetKey::String
    end
    ElastoplasticLoadingInput(; E=75.0, sigma_y=100.0, H=0.0, eta=1.0e24, nu=0.2412,
        rho=2700.0, L=15.0, presetKey="interseismic") =
        ElastoplasticLoadingInput(E, sigma_y, H, eta, nu, rho, L, presetKey)

    Base.get(w::ElastoplasticLoadingInput) = Dict{String,Any}(
        "E" => w.E, "sigma_y" => w.sigma_y, "H" => w.H, "eta" => w.eta, "nu" => w.nu,
        "rho" => w.rho, "L" => w.L, "presetKey" => w.presetKey,
        "stressHistory" => [0.0], "strainHistory" => [0.0], "timeHistory" => [0.0],
        "controlMode" => "stress", "epoch" => 0)

    """
    	Base.show(io, ::MIME"text/html", w::ElastoplasticLoadingInput)

    Render the pull-release scene (a block made of an internal 2-D spring lattice, pulled
    from a handle, with a scrubbable spring+dashpot Maxwell element between the wall and
    the block), a "two clocks" comparison panel, and a graph column that shows the live
    pull-release stress-strain trace during manual dragging or a preset's scripted
    time-series right after that preset fires.
    """
    function Base.show(io::IO, ::MIME"text/html", w::ElastoplasticLoadingInput)
        write(io, """
        <div id="eplwidget">
        <style>
        pluto-cell:has(#eplwidget) { width: min(85vw, 1150px) !important;
          margin-left: calc((100% - min(85vw, 1150px)) / 2) !important; }
        #eplwidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #eplwidget .epl-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #eplwidget .epl-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #eplwidget .epl-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #eplwidget .epl-main-row{display:flex;gap:14px;align-items:flex-start;flex-wrap:wrap}
        #eplwidget .epl-block-col{flex:3 1 420px;min-width:0}
        #eplwidget .epl-graph-col{flex:1 1 220px;min-width:200px}
        #eplwidget .epl-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px;margin-bottom:8px}
        #eplwidget .epl-panel-title{font-size:14px;font-weight:700;color:#e5e7eb;text-align:center;margin-bottom:6px}
        #eplwidget .epl-caption{font-size:12px;color:#9ca3af;text-align:center;margin-top:4px;margin-bottom:10px}
        #eplwidget .epl-clocks-panel{width:100%;box-sizing:border-box;text-align:center;font-size:13px;color:#e5e7eb;
          background:#0a0f18;border:1px solid #eab308;border-radius:6px;padding:8px;margin-bottom:8px}
        #eplwidget .epl-clocks-panel b{color:#facc15}
        #eplwidget #epl-clocks-text{margin-top:2px;line-height:1.5}
        #eplwidget .epl-readouts{width:100%;box-sizing:border-box;text-align:center;font-size:14px;color:#e5e7eb;
          background:#050505;border:1px solid #2f3744;border-radius:6px;padding:8px;margin-bottom:14px}
        #eplwidget canvas{display:block;width:100%}
        #eplwidget #epl-scene{cursor:default}
        #eplwidget .epl-presets-row{display:flex;justify-content:center;gap:8px;flex-wrap:wrap;margin-bottom:14px}
        #eplwidget .epl-preset-btn{border-radius:4px;border:1px solid #6b7280;background:#0b0b0b;color:#e5e7eb;padding:6px 10px;font-size:12px;cursor:pointer}
        #eplwidget .epl-preset-btn:hover{background:#1f2937}
        #eplwidget .epl-preset-btn.active{border-color:#38bdf8;color:#38bdf8}
        #eplwidget .epl-mode-row{display:flex;justify-content:center;gap:8px;flex-wrap:wrap;margin-bottom:4px}
        #eplwidget .epl-mode-btn{border-radius:4px;border:1px solid #6b7280;background:#0b0b0b;color:#e5e7eb;padding:6px 10px;font-size:12px;cursor:pointer}
        #eplwidget .epl-mode-btn:hover{background:#1f2937}
        #eplwidget .epl-mode-btn.active{border-color:#4ade80;color:#4ade80}
        #eplwidget .epl-mode-desc{font-size:12px;color:#9ca3af;text-align:center;margin-bottom:12px}
        #eplwidget .epl-control-group{background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 14px;margin-bottom:14px}
        #eplwidget .epl-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #eplwidget .epl-control-row{display:grid;grid-template-columns:110px minmax(0,1fr) 90px;gap:8px;align-items:center;margin:6px 0}
        #eplwidget .epl-control-row label{font-size:13px;color:#9ca3af}
        #eplwidget .epl-control-row input[type=range]{width:100%;min-width:0}
        #eplwidget .epl-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #eplwidget .epl-material-desc{font-size:12px;color:#9ca3af;text-align:center;margin-top:4px}
        #eplwidget canvas.epl-scrub{height:34px;cursor:grab;border:1px solid #2f3744;border-radius:4px;background:#000}
        #eplwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px;cursor:pointer}
        #eplwidget button:hover{background:#767676}
        @media (max-width: 640px) { #eplwidget .epl-block-col{flex-basis:100%} #eplwidget .epl-graph-col{flex-basis:100%} }
        </style>

        <div class="epl-title">
          <div class="epl-title-desc">One block, three regimes -- elastic, viscoelastic, or viscous flow, depending on the clock you compare it to.</div>
          <div class="epl-title-hint">drag the handle by hand, or click a preset to fast-forward through a realistic scenario -- hold a hand-pull steady and watch &sigma;_y/plasticity respond (real &eta; won't visibly creep in a few seconds, that's what the presets are for) &middot; drag the spring for stiffness, the dashpot for viscosity</div>
        </div>

        <div class="epl-main-row">
          <div class="epl-block-col">
            <div class="epl-mode-row" id="epl-mode-row">
              <button class="epl-mode-btn" data-mode="stress" type="button">Pull it (constant stress → creep)</button>
              <button class="epl-mode-btn" data-mode="strain" type="button">Stretch &amp; hold (constant strain → relax)</button>
            </div>
            <div class="epl-mode-desc" id="epl-mode-desc"></div>
            <div class="epl-panel-title">Pull the Block</div>
            <div class="epl-panel"><canvas id="epl-scene"></canvas></div>
            <div class="epl-caption">drag the spring to change stiffness, the dashpot to change viscosity &middot; scrub bracing below to change &nu; &middot; strain exaggerated for visibility</div>

            <div class="epl-clocks-panel">
              <canvas id="epl-timescale"></canvas>
              <div id="epl-clocks-text">&tau;_wave = -- &middot; &tau;_M = --</div>
            </div>
            <div class="epl-readouts" id="epl-readouts">stress = -- &middot; strain = -- &middot; lateral = --</div>

            <div class="epl-control-group">
              <div class="epl-control-title">Material</div>
              <div class="epl-control-row"><label>bracing</label><canvas id="epl-bracing-scrub" class="epl-scrub"></canvas><span class="epl-value" id="epl-bracing-v"></span></div>
              <div class="epl-control-row"><label id="epl-sy-label">&sigma;_y (yield, MPa)</label><input type="range" id="epl-sy" min="10" max="500" step="1" value="$(w.sigma_y)"><span class="epl-value" id="epl-sy-v"></span></div>
              <div class="epl-control-row"><label>H (hardening)</label><input type="range" id="epl-H" min="0" max="60" step="1" value="$(w.H)"><span class="epl-value" id="epl-H-v"></span></div>
              <div class="epl-material-desc" id="epl-material-desc"></div>
              <button id="epl-reset" type="button">Reset</button>
            </div>

            <div class="epl-presets-row" id="epl-presets-row">
              <button class="epl-preset-btn" data-preset="elastic" type="button">Elastic Loading</button>
              <button class="epl-preset-btn" data-preset="interseismic" type="button">Interseismic Loading</button>
              <button class="epl-preset-btn" data-preset="postseismic" type="button">Postseismic Relaxation</button>
              <button class="epl-preset-btn" data-preset="convection" type="button">Mantle Convection</button>
            </div>
          </div>

          <div class="epl-graph-col">
            <div class="epl-panel-title">Stress-Strain</div>
            <div class="epl-panel"><canvas id="epl-graph"></canvas></div>
            <div class="epl-caption">accumulates every hand-pull cycle</div>

            <div class="epl-panel-title" id="epl-ts-title">Time Series</div>
            <div class="epl-panel"><canvas id="epl-timeseries"></canvas></div>
            <div class="epl-caption" id="epl-ts-caption">stress (left axis, MPa) &middot; axial &amp; lateral strain (right axis) vs. time</div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        // Pluto can execute a cell's <script> tag more than once against the SAME
        // persistent DOM nodes -- e.g. a client's initial connection replays a few
        // intermediate cell outputs before the final one, and later edits rerun this cell
        // again. Every setup call below (event listeners, the rAF playback loop) targets
        // elements that survive across those reruns, so without this guard each rerun
        // stacks ANOTHER full set of duplicate listeners/loops racing each other on the
        // same canvas -- confirmed directly: a stray duplicate execution was overwriting
        // freshly-drawn labels with stale ones on every redraw. Only the FIRST execution
        // for this particular widget instance actually wires anything up.
        if(!par._eplInitialized){
        par._eplInitialized = true;
        // state.stressHistory is the exact input the student has applied (force control);
        // displayStrain is Julia's (or a local elastic prediction of) the resulting strain,
        // index-aligned with stressHistory -- never push to one without the other. Every
        // field below is a REAL unit (GPa/Pa·s/MPa/kg per m³/km) -- there is only one
        // parameter set now; manual dragging simply shows no visible creep at a realistic
        // eta (correctly -- creep needs tau_M, which presets fast-forward through instead).
        let state = { E: $(w.E), sigma_y: $(w.sigma_y), H: $(w.H), eta: $(w.eta), nu: $(w.nu),
          rho: $(w.rho), L: $(w.L), presetKey: "$(w.presetKey)",
          bracingIdx: 1, // matches the struct's default nu (nu_medium) above -- keep in sync
          // 'stress': the free-play pull-release-repeat test above (force is driven, strain
          // creeps). 'strain': its displacement-controlled dual -- drag sets a target strain
          // directly, and holding it fixed lets the derived stress relax instead. Both modes
          // share every other field below; only which quantity is driven vs. derived flips,
          // and Julia (see stress_from_strain_with_relaxation in the Appendix) is what
          // actually computes whichever side isn't driven, exactly as strain_from_stress_with_creep
          // already does for the 'stress' mode.
          controlMode: 'stress',
          stressHistory: [0], timeHistory: [0], stress: 0, strain: 0, epsP: 0, lateralStrain: 0,
          currentT: null, currentTLabel: '', graphMode: 'live' };
        let pushed = null;
        let pushedAnim = null; // {animTime, animStress(MPa), animStrain, animLateral} for state.presetKey
        let awaitingPresetAnim = null; // preset key we're waiting for FRESH anim data for
        let animPlaying = false, animStartMs = 0;
        const ANIM_REAL_MS = 4000; // every preset animates over the same ~4 real seconds, whatever its true duration
        let graphRevertTimer = null;
        let displayStrain = [0];
        let displayLateral = [0]; // parallel to displayStrain -- the transverse (necking) strain history
        // index into state.timeHistory/stressHistory/displayStrain/displayLateral where
        // each fresh pull begins -- each new mousedown restarts dragSimTime near zero (see
        // dragSimTime above), so connecting cycle N's tail directly to cycle N+1's head as
        // one continuous line jumps backward in (log) time and draws a zigzagging tangle.
        // drawTimeSeries uses this to draw each cycle as its own polyline instead.
        let cycleStarts = [0];
        // click-to-toggle visibility for the Time Series legend -- tsLegendBoxes is
        // recomputed (from actual measured text width) every drawTimeSeries call, so the
        // click handler below always hit-tests against the labels' real current position
        let tsVisible = { stress: true, axial: true, lateral: true };
        let tsLegendBoxes = {};
        let commitInFlight = false;
        let dragging = false;
        let releasing = false;
        let cursorPx = null;
        let holdTimer = null;
        const t0 = performance.now();
        function elapsed(){ return (performance.now()-t0)/1000; }

        // par.clientWidth is already forced synchronously by this widget's own
        // `pluto-cell:has(#eplwidget){width:min(85vw,1150px)!important}` CSS rule (applied
        // before this script ever runs, unlike WideCell's async JS resize elsewhere), so
        // it's trustworthy here without deferring via a ResizeObserver. Capping it further
        // against window.innerWidth*0.85 is still wrong whenever that fraction is smaller
        // than the real (CSS-bounded) clientWidth -- it silently shrinks every canvas
        // below its own wrapper's width, leaving a dead stripe of unused black background.
        const availW = Math.min(par.clientWidth || (window.innerWidth*0.85 || 900), 1150);
        const PW = Math.max(300, Math.floor(availW*0.72) - 14);
        const PH = Math.max(190, Math.round(PW*0.34));
        const PW_GRAPH = Math.max(200, Math.floor(availW*0.25));
        const PH_GRAPH = Math.max(200, Math.round(PW_GRAPH*0.85));
        const PH_TS = Math.max(190, Math.round(PW_GRAPH*0.8));
        const DPR = window.devicePixelRatio || 1;

        function hidpi(canvas, ctx, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          ctx.setTransform(DPR,0,0,DPR,0,0);
        }

        const sceneCv = par.querySelector('#epl-scene'), sceneCtx = sceneCv.getContext('2d');
        hidpi(sceneCv, sceneCtx, PW, PH);
        const graphCv = par.querySelector('#epl-graph'), graphCtx = graphCv.getContext('2d');
        hidpi(graphCv, graphCtx, PW_GRAPH, PH_GRAPH);
        const tsCv = par.querySelector('#epl-timeseries'), tsCtx = tsCv.getContext('2d');
        hidpi(tsCv, tsCtx, PW_GRAPH, PH_TS);
        const TS_H = 70;
        const timescaleCv = par.querySelector('#epl-timescale'), timescaleCtx = timescaleCv.getContext('2d');
        hidpi(timescaleCv, timescaleCtx, PW, TS_H);
        const readoutEl = par.querySelector('#epl-readouts');
        const clocksTextEl = par.querySelector('#epl-clocks-text');
        const tsTitleEl = par.querySelector('#epl-ts-title');
        const tsCaptionEl = par.querySelector('#epl-ts-caption');
        const materialDescEl = par.querySelector('#epl-material-desc');
        const syLabelEl = par.querySelector('#epl-sy-label');
        const modeDescEl = par.querySelector('#epl-mode-desc');

        // bracing sets nu to a value the Appendix actually DERIVES from a small
        // spring-lattice solve, not an assigned number
        const bracingPresets = [
          { name: 'light bracing',  nu: $(nu_light) },
          { name: 'medium bracing', nu: $(nu_medium) },
          { name: 'full bracing',   nu: $(nu_full) }
        ];
        // real-unit regime bundles from the Appendix's regime_presets -- each sets material
        // + geometry + a SCRIPTED forcing (kind/amplitude/duration, not a period) chosen
        // relative to THAT material's own two clocks, so clicking one always lands in its
        // named regime
        const REGIME_PRESETS = {
          elastic: { name: "$(regime_presets.elastic.name)", E: $(regime_presets.elastic.E), bracingIdx: $(regime_presets.elastic.bracingIdx), eta: $(regime_presets.elastic.eta), rho: $(regime_presets.elastic.rho), L: $(regime_presets.elastic.L), sigma_y: $(regime_presets.elastic.sigma_y), kind: "$(regime_presets.elastic.kind)", amplitude: $(regime_presets.elastic.amplitude), duration: $(regime_presets.elastic.duration) },
          interseismic: { name: "$(regime_presets.interseismic.name)", E: $(regime_presets.interseismic.E), bracingIdx: $(regime_presets.interseismic.bracingIdx), eta: $(regime_presets.interseismic.eta), rho: $(regime_presets.interseismic.rho), L: $(regime_presets.interseismic.L), sigma_y: $(regime_presets.interseismic.sigma_y), kind: "$(regime_presets.interseismic.kind)", amplitude: $(regime_presets.interseismic.amplitude), duration: $(regime_presets.interseismic.duration) },
          postseismic: { name: "$(regime_presets.postseismic.name)", E: $(regime_presets.postseismic.E), bracingIdx: $(regime_presets.postseismic.bracingIdx), eta: $(regime_presets.postseismic.eta), rho: $(regime_presets.postseismic.rho), L: $(regime_presets.postseismic.L), sigma_y: $(regime_presets.postseismic.sigma_y), kind: "$(regime_presets.postseismic.kind)", amplitude: $(regime_presets.postseismic.amplitude), duration: $(regime_presets.postseismic.duration) },
          convection: { name: "$(regime_presets.convection.name)", E: $(regime_presets.convection.E), bracingIdx: $(regime_presets.convection.bracingIdx), eta: $(regime_presets.convection.eta), rho: $(regime_presets.convection.rho), L: $(regime_presets.convection.L), sigma_y: $(regime_presets.convection.sigma_y), kind: "$(regime_presets.convection.kind)", amplitude: $(regime_presets.convection.amplitude), duration: $(regime_presets.convection.duration) }
        };
        const bracingCv = par.querySelector('#epl-bracing-scrub'), bracingCtx = bracingCv.getContext('2d');
        const SCRUB_H = 34;
        function sizeScrub(cv, ctx){
          const w = cv.getBoundingClientRect().width || 150;
          hidpi(cv, ctx, w, SCRUB_H);
          return w;
        }
        const bracingScrubW = sizeScrub(bracingCv, bracingCtx);

        function drawScrubStrip(ctx, w, presets, selectedIdx){
          ctx.clearRect(0, 0, w, SCRUB_H);
          const zoneW = w / 3;
          for(let z=0; z<3; z++){
            const zx = z*zoneW, cx = zx + zoneW/2, cy = SCRUB_H/2;
            if(z === selectedIdx){
              ctx.fillStyle = 'rgba(56,189,248,0.18)';
              ctx.fillRect(zx, 0, zoneW, SCRUB_H);
            }
            const col = z === selectedIdx ? '#38bdf8' : '#6b7280';
            ctx.strokeStyle = col; ctx.fillStyle = col; ctx.lineWidth = 1;
            // all three levels are genuinely braced -- no "zero" icon, see the aside
            // for why an unbraced grid isn't offered as one of the three options
            const s = 9;
            ctx.beginPath();
            ctx.moveTo(cx-s,cy-s); ctx.lineTo(cx+s,cy-s); ctx.lineTo(cx+s,cy+s); ctx.lineTo(cx-s,cy+s); ctx.closePath();
            ctx.stroke();
            ctx.lineWidth = z >= 2 ? 1.4 : 1;
            ctx.beginPath(); ctx.moveTo(cx-s,cy-s); ctx.lineTo(cx+s,cy+s); ctx.stroke();
            if(z >= 1){ ctx.beginPath(); ctx.moveTo(cx+s,cy-s); ctx.lineTo(cx-s,cy+s); ctx.stroke(); }
          }
          ctx.strokeStyle = '#2f3744'; ctx.lineWidth = 1;
          for(let z=1; z<3; z++){ ctx.beginPath(); ctx.moveTo(z*zoneW,0); ctx.lineTo(z*zoneW,SCRUB_H); ctx.stroke(); }
        }

        // geometry: wall -- [spring] -- [dashpot] -- lattice block -- handle. The spring
        // and dashpot ARE the Maxwell element (elastic + viscous, in series) that the
        // lattice's own bulk response approximates; drawing them explicitly is what makes
        // them directly draggable instead of buried in a slider.
        const WALL_X = 30;
        const SPRING_X0 = WALL_X + 12, SPRING_LEN = Math.round(PW*0.09);
        const SPRING_X1_ELEM = SPRING_X0 + SPRING_LEN;
        const DASHPOT_X0 = SPRING_X1_ELEM + Math.round(PW*0.02), DASHPOT_LEN = Math.round(PW*0.07);
        const DASHPOT_X1 = DASHPOT_X0 + DASHPOT_LEN;
        const BLOCK_X0 = DASHPOT_X1 + 16; // left edge of the deforming lattice
        const SPRING_X1 = BLOCK_X0 + Math.round(PW*0.24); // rest position of the lattice's right edge
        const BLOCK_TRAVEL = Math.round(PW*0.20); // max growth of the right edge, responsive to panel width
        const BLOCK_H = 78;
        const MESH_NX = 6, MESH_NY = 5; // fixed mesh -- no longer tied to a density scrub
        const MESH_SPRING_AMP = 2.5, MESH_SPRING_FREQ = 5;
        const LATERAL_AMP = 3; // visual exaggeration of necking, for visibility
        const BLOCK_H_MIN = BLOCK_H * 0.2;
        const PX_PER_STRESS_UNIT = 2.5; // px of drag per MPa -- tuned for the sigma_y range below
        // px of drag per unit strain -- calibrated against PX_PER_STRESS_UNIT at the
        // widget's own defaults (E=75 GPa, sigma_y=100 MPa) so reaching the elastic limit
        // (eps=sigma_y/(1000E)~=0.00133) takes ~200px, comparable to how far a stress-mode
        // pull needs to reach the same 100 MPa (250px at 2.5px/MPa). The OLD value (30000)
        // reached that same strain in only ~40px -- any normal drag instantly overshot into
        // the region localStressFromStrain and Julia disagreed about (see its own comment),
        // which was the dominant cause of the reported jitter, not just a UI-feel issue.
        const PX_PER_STRAIN_UNIT = 150000;
        const PLASTIC_CAP = 1.5; // must match strain_from_stress's default in the Appendix
        // floors for the pull-release graph's axes (see graphBounds below) -- just large
        // enough that a near-zero sigma_y/E material doesn't zoom the axis in on noise
        const GRAPH_STRAIN_MAX = 1e-4, GRAPH_STRESS_MAX = 20;
        const E_BOUNDS = [20, 400];       // GPa
        const ETA_BOUNDS = [1e17, 1e24];  // Pa·s, log-scale
        const SY_BOUNDS = [10, 500];      // MPa

        function emit(){
          commitInFlight = true;
          // displayStrain doubles as the driven strain history in 'strain' mode (nothing else
          // ever needs to correct it there -- see pushStrain) and as the Julia-corrected
          // creep-aware strain trace in 'stress' mode; either way it's always the right thing
          // to send as strainHistory. Symmetric with stressHistory below.
          par.value = { E: state.E, sigma_y: state.sigma_y, H: state.H, eta: state.eta, nu: state.nu,
            rho: state.rho, L: state.L, presetKey: state.presetKey, controlMode: state.controlMode,
            stressHistory: state.stressHistory, strainHistory: displayStrain,
            timeHistory: state.timeHistory, epoch: epoch };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledEmit(){ if(!commitInFlight) emit(); }

        // local mirror of the Appendix's strain_from_stress -- used only for smooth,
        // immediate visual feedback while dragging; Julia's response is authoritative
        // and overwrites this the moment it arrives. Stress/strain here are real units
        // (MPa in, dimensionless strain out via E in GPa) -- same formula, real numbers.
        function localStrainFromStress(sigmaRaw, epsPPrev){
          // H=0 (the default everywhere) is a numerical singularity for this algebraic
          // update -- see the matching comment on strain_from_stress_with_creep in the
          // Appendix, which this mirrors. 100*state.E is 0.1*E converted GPa->MPa to match
          // H's own units.
          const H_eff = Math.max(state.H, 100*state.E);
          const capMag = state.sigma_y + H_eff*PLASTIC_CAP;
          const sigma = Math.max(-capMag, Math.min(capMag, sigmaRaw));
          const yieldNow = state.sigma_y + H_eff*Math.abs(epsPPrev);
          let epsP = epsPPrev;
          // check the RAW attempted stress, not the already-clamped one -- clamping first
          // and checking the clamped value against the same bound it was clamped to is
          // exactly the bug this replaces (yielding could never trigger)
          if(Math.abs(sigmaRaw) > yieldNow){
            const s = Math.sign(sigmaRaw);
            epsP = s * (Math.abs(sigma) - state.sigma_y) / H_eff;
          }
          return [epsP + (sigma*1e6)/(state.E*1e9), epsP, sigma];
        }

        // local mirror of the elastic branch of the Appendix's stress_from_strain_with_relaxation
        // -- used only for smooth, immediate visual feedback while dragging in 'strain' mode;
        // Julia's relaxation-corrected response is authoritative and overwrites this the
        // moment it arrives. Ignores relaxation entirely (deliberately -- relaxation only
        // happens over elapsed time, which a single instantaneous drag sample doesn't have).
        // Clamping this at sigma_y (the YIELD stress) was a real bug: the real Appendix
        // function models no plasticity at all and only clamps at a much more generous
        // sigma_cap safety valve (see the Julia cell that calls it) -- with a tight sigma_y
        // clamp here instead, any drag past the elastic limit (trivially easy at the drag
        // sensitivity below) pinned this local preview flat while Julia's real, unclamped-
        // until-much-higher answer kept climbing, so every round trip snapped the displayed
        // stress between two very different values -- the actual cause of the "jittery"
        // feel reported live. Matching Julia's own sigma_cap formula exactly fixes it.
        function localStressFromStrain(epsRaw){
          const trial = 1000 * state.E * epsRaw; // GPa*strain -> MPa
          const H_eff = Math.max(state.H, 100 * state.E); // matches Julia's 0.1*E_val*1000 (GPa->MPa)
          const capMag = state.sigma_y + H_eff * 1.5;
          return Math.max(-capMag, Math.min(capMag, trial));
        }

        // local mirror of the Appendix's lateral_strain_from_axial -- same approximation
        // as localStrainFromStress (plastic-only, ignores creep, corrected by Julia's
        // authoritative response the moment it arrives).
        function localLateralStrain(strain, epsP){
          return -state.nu*(strain-epsP) - 0.5*epsP;
        }

        // local mirrors of the Appendix's wave_relaxation_clocks/deborah_regime -- pure
        // algebra on the current material, no history involved, so no need to wait on
        // Julia for these (unlike the strain response itself)
        function GFromEnu(E,nu){ return E/(2*(1+nu)); }
        function waveRelaxationClocks(E_GPa, nu, eta_Pas, L_km, rho){
          const mu_Pa = GFromEnu(E_GPa, nu) * 1e9;
          const Vs = Math.sqrt(mu_Pa/rho);
          return [ (L_km*1e3)/Vs, eta_Pas/mu_Pa ];
        }
        function deborahRegime(tauWave, tauM, T){
          const Rwave = tauWave/T, De = tauM/T;
          let regime;
          if(Rwave > 0.15) regime = 'seismic wave';
          else if(De > 3) regime = 'elastic quasi-static';
          else if(De > 0.3) regime = 'viscoelastic';
          else regime = 'viscous flow';
          return { Rwave, De, regime };
        }
        function formatDuration(s){
          const yr = 365.25*24*3600, a = Math.abs(s);
          if(a < 60) return s.toFixed(2)+' s';
          if(a < 3600) return (s/60).toFixed(2)+' min';
          if(a < 86400) return (s/3600).toFixed(2)+' hr';
          if(a < yr) return (s/86400).toFixed(2)+' day';
          if(a < 1e3*yr) return (s/yr).toFixed(2)+' yr';
          if(a < 1e6*yr) return (s/(1e3*yr)).toFixed(2)+' kyr';
          return (s/(1e6*yr)).toFixed(2)+' Myr';
        }
        // the SAME log-seconds range the timescale number line spans (see drawTimescale
        // below) -- shared here so a manual pull can be mapped onto exactly that range
        const TS_MIN_LOG = -2, TS_MAX_LOG = 16; // seconds: ~0.01s to ~300 Myr
        // nonlinear drag-to-simulated-time map: a real hand-hold can only ever last a few
        // seconds, but the whole point of this widget is comparing timescales that span
        // seconds to millions of years, so a HELD pull's real duration is mapped, on a log
        // scale, across that entire range -- roughly a 12-second hold sweeps from ~0.01s
        // to ~300 Myr of SIMULATED time, which is what actually gets fed to the creep
        // integrator (via state.timeHistory) so holding long enough genuinely creeps, at
        // whatever eta is currently dialed in, exactly like a preset would show.
        const DRAG_MAX_REAL_S = 12;
        function dragSimTime(realDt){
          if(realDt <= 0) return 0;
          const frac = Math.min(1, realDt/DRAG_MAX_REAL_S);
          return Math.pow(10, TS_MIN_LOG + frac*(TS_MAX_LOG-TS_MIN_LOG));
        }
        // plain-language qualifiers for where a value sits in a realistic Earth range --
        // bins are illustrative, not a precise rock-physics lookup
        function describeStiffness(E_GPa){
          if(E_GPa < 50) return 'weak crust';
          if(E_GPa < 100) return 'crust-like';
          if(E_GPa < 200) return 'upper mantle';
          return 'lower mantle';
        }
        function describeViscosity(eta_Pas){
          if(eta_Pas > 1e22) return 'lithosphere, effectively elastic';
          if(eta_Pas > 1e20) return 'lower mantle';
          if(eta_Pas > 1e18) return 'asthenosphere';
          return 'weak zone';
        }
        // shared color scale for BOTH the dashpot fluid and the lattice fibers, so the
        // whole block visibly agrees with the dashpot about how viscous the material is
        function viscosityColor(etaFrac){
          const r = Math.round(255*(1-0.55*etaFrac)), g = Math.round(230*(1-0.7*etaFrac)), b = Math.round(255-90*etaFrac);
          return 'rgb('+r+','+g+','+b+')';
        }
        function currentEtaFrac(){
          return (Math.log(state.eta)-Math.log(ETA_BOUNDS[0]))/(Math.log(ETA_BOUNDS[1])-Math.log(ETA_BOUNDS[0]));
        }

        function pushStress(sigmaRaw, force){
          const [strain, epsP, sigma] = localStrainFromStress(sigmaRaw, state.epsP);
          const lateral = localLateralStrain(strain, epsP);
          const last = state.stressHistory[state.stressHistory.length-1];
          state.stress = sigma; state.strain = strain; state.epsP = epsP;
          state.lateralStrain = lateral;
          if(!force && Math.abs(sigma-last) <= 0.3) return false;
          // simulated (not real wall-clock) time -- see dragSimTime above -- so a held
          // pull genuinely creeps at whatever eta is currently dialed in, exactly like a
          // preset would show, instead of needing an unrealistically tiny eta to see
          // anything move within a human-length hold
          const simT = dragSimTime(elapsed()-dragStartTime);
          state.stressHistory.push(sigma);
          state.timeHistory.push(simT);
          displayStrain.push(strain);
          displayLateral.push(lateral);
          if(state.stressHistory.length > 4000){
            state.stressHistory.shift(); state.timeHistory.shift(); displayStrain.shift(); displayLateral.shift();
          }
          // live "your pull" marker on the timescale panel -- crawls across the SAME log
          // axis as tau_wave/tau_M while you hold, not just once you release
          state.currentT = Math.max(0.01, 2*simT);
          state.currentTLabel = 'your pull (holding)';
          return true;
        }

        // the displacement-controlled dual of pushStress: epsRaw IS the driven quantity here
        // (no correction needed, unlike strain in 'stress' mode -- see localStressFromStrain),
        // and the local elastic stress guess pushed into stressHistory is what Julia's
        // relaxation-aware response (stress_from_strain_with_relaxation, Appendix) overwrites
        // once it arrives, mirroring exactly how displayStrain gets overwritten in 'stress' mode.
        function pushStrain(epsRaw, force){
          const sigma = localStressFromStrain(epsRaw);
          const lateral = localLateralStrain(epsRaw, 0); // relaxation residual unknown until Julia responds
          const last = displayStrain[displayStrain.length-1];
          state.strain = epsRaw; state.stress = sigma; state.lateralStrain = lateral;
          if(!force && Math.abs(epsRaw-last) <= 2e-6) return false;
          const simT = dragSimTime(elapsed()-dragStartTime);
          state.stressHistory.push(sigma);
          state.timeHistory.push(simT);
          displayStrain.push(epsRaw);
          displayLateral.push(lateral);
          if(state.stressHistory.length > 4000){
            state.stressHistory.shift(); state.timeHistory.shift(); displayStrain.shift(); displayLateral.shift();
          }
          state.currentT = Math.max(0.01, 2*simT);
          state.currentTLabel = 'your pull (holding)';
          return true;
        }

        // keeps sampling (time, currentStress/currentStrain) while the block is held still
        // under a fixed pull -- without this, holding the mouse steady would never advance
        // state.timeHistory, and creep/relaxation (which only happen over elapsed time) would
        // be invisible unless you kept wiggling the mouse. Same mechanism drives both modes;
        // only which quantity is re-pushed (held fixed) while time advances differs.
        function startHoldSampling(){
          if(holdTimer) return;
          holdTimer = setInterval(() => {
            if(!dragging) return;
            if(state.controlMode === 'strain') pushStrain(state.strain, true);
            else pushStress(state.stress, true);
            draw();
            throttledEmit();
          }, 150);
        }
        function stopHoldSampling(){
          if(holdTimer){ clearInterval(holdTimer); holdTimer = null; }
        }

        function drawSpring(x1,x2,y,amp,freq){
          sceneCtx.beginPath();
          const steps = 60;
          for(let j=0;j<=steps;j++){
            const f = j/steps, x = x1+(x2-x1)*f, yy = y + Math.sin(f*Math.PI*freq)*amp;
            j===0 ? sceneCtx.moveTo(x,yy) : sceneCtx.lineTo(x,yy);
          }
          sceneCtx.stroke();
        }
        // vertical variant -- same zigzag math, x/y roles swapped, for the springs that
        // show the lateral (necking) response
        function drawSpringV(x,y1,y2,amp,freq){
          sceneCtx.beginPath();
          const steps = 40;
          for(let j=0;j<=steps;j++){
            const f = j/steps, y = y1+(y2-y1)*f, xx = x + Math.sin(f*Math.PI*freq)*amp;
            j===0 ? sceneCtx.moveTo(xx,y) : sceneCtx.lineTo(xx,y);
          }
          sceneCtx.stroke();
        }
        // a real dashpot symbol -- a fluid-filled cylinder, closed on the spring side, open
        // on the block side, with a piston plate inside connected to a rod that exits the
        // open end -- rather than the plain rectangle-with-a-line placeholder this replaced.
        // `full` controls how large a version to draw (the macroscopic element vs. the tiny
        // copies embedded in the lattice fibers below).
        function drawDashpotIcon(ctx, x0, x1, midY, dh, color, full=true){
          const len = x1 - x0, cylLen = len * 0.62, cylX1 = x0 + cylLen;
          ctx.fillStyle = color;
          ctx.fillRect(x0, midY-dh, cylLen, 2*dh);
          ctx.strokeStyle = full ? '#6b7280' : color;
          ctx.lineWidth = full ? 1.4 : 1;
          ctx.beginPath();
          ctx.moveTo(cylX1, midY-dh); ctx.lineTo(x0, midY-dh); ctx.lineTo(x0, midY+dh); ctx.lineTo(cylX1, midY+dh);
          ctx.stroke();
          const pistonX = x0 + cylLen*0.58;
          const inset = full ? 2 : 1;
          ctx.fillStyle = full ? '#9ca3af' : color;
          ctx.fillRect(pistonX - (full?1.5:1), midY-dh+inset, full?3:1.4, 2*dh-2*inset);
          ctx.strokeStyle = full ? '#9ca3af' : color;
          ctx.lineWidth = full ? 2.2 : 1.2;
          ctx.beginPath(); ctx.moveTo(pistonX, midY); ctx.lineTo(x1, midY); ctx.stroke();
        }
        // the wall-to-block Maxwell element: a coil (stiffness) in series with a dashpot
        // (viscosity) -- literally what the lattice's own bulk response approximates, drawn
        // explicitly so both can be scrubbed directly instead of buried in a slider. More
        // turns = stiffer spring; darker, thicker-looking fluid = more viscous dashpot.
        function drawMaxwellElement(midY){
          const ctx = sceneCtx;
          const eFrac = (state.E-E_BOUNDS[0])/(E_BOUNDS[1]-E_BOUNDS[0]);
          const etaFrac = currentEtaFrac();
          ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 1.5;
          ctx.beginPath(); ctx.moveTo(WALL_X, midY); ctx.lineTo(SPRING_X0, midY); ctx.stroke();
          ctx.strokeStyle = '#e5e7eb'; ctx.lineWidth = 1.6;
          drawSpring(SPRING_X0, SPRING_X1_ELEM, midY, 10 - 4*eFrac, 3 + 9*eFrac);
          ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 1.5;
          ctx.beginPath(); ctx.moveTo(SPRING_X1_ELEM, midY); ctx.lineTo(DASHPOT_X0, midY); ctx.stroke();
          const dh = 15;
          drawDashpotIcon(ctx, DASHPOT_X0, DASHPOT_X1, midY, dh, viscosityColor(etaFrac));
          ctx.strokeStyle = '#9ca3af';
          ctx.beginPath(); ctx.moveTo(DASHPOT_X1, midY); ctx.lineTo(BLOCK_X0, midY); ctx.stroke();
          ctx.fillStyle = '#e5e7eb'; ctx.textAlign = 'center';
          // short numeric labels only -- the descriptive qualifier ("crust-like" etc.) has
          // more room in the Material control panel below than this narrow gap allows
          ctx.font = '10px sans-serif';
          ctx.fillText(state.E.toFixed(0)+' GPa', (SPRING_X0+SPRING_X1_ELEM)/2, midY-dh-14);
          ctx.fillText(state.eta.toExponential(1)+' Pa·s', (DASHPOT_X0+DASHPOT_X1)/2, midY+dh+18);
        }
        function pointerOnSpring(px,py,midY){ return px>=SPRING_X0-6 && px<=SPRING_X1_ELEM+6 && Math.abs(py-midY)<=18; }
        function pointerOnDashpot(px,py,midY){ return px>=DASHPOT_X0-6 && px<=DASHPOT_X1+6 && Math.abs(py-midY)<=18; }

        function strainMaxCurrent(){
          // realistic elastic strains are ~1e-5 to a few %, nothing like the old
          // arbitrary-unit widget's ~0.15 floor -- keep a small floor just so the block
          // doesn't visually flatten to a line before any pull has happened
          let m = 0.001;
          for(const e of displayStrain) m = Math.max(m, Math.abs(e));
          // during a preset's own scripted playback, the relevant range is THAT preset's
          // strain, not the (possibly still-empty) manual pull history -- without this, a
          // preset whose strain exceeds the tiny manual-pull floor (mantle convection
          // creeps to several percent) stretched the block's right edge far off-canvas
          if(state.graphMode === 'preset' && pushedAnim){
            for(const e of pushedAnim.animStrain) m = Math.max(m, Math.abs(e));
          }
          return m * 1.15;
        }
        function currentRightX(eMax){
          // defensive clamp -- the block's right edge should never be able to leave the
          // canvas regardless of how eMax and state.strain end up related
          const ratio = Math.max(-1, Math.min(1, state.strain/eMax));
          return SPRING_X1 + ratio * BLOCK_TRAVEL;
        }
        function currentBlockH(){
          return Math.max(BLOCK_H_MIN, BLOCK_H * (1 + LATERAL_AMP*state.lateralStrain));
        }
        function currentHandleX(){ return currentRightX(strainMaxCurrent()) + 18; }
        function pointerOnHandle(px, py, midY){
          // only the handle applies force -- the lattice is the medium responding, not
          // something you grab directly
          const hx = currentHandleX();
          const dx = px-hx, dy = py-midY;
          return Math.sqrt(dx*dx + dy*dy) <= 14;
        }

        function drawScene(){
          const ctx = sceneCtx, W = PW, H = PH;
          ctx.clearRect(0,0,W,H);
          const midY = H/2;
          const rightX = currentRightX(strainMaxCurrent());
          const blockH = currentBlockH();
          const leftX = BLOCK_X0;
          const topY = midY - blockH/2, botY = midY + blockH/2;
          // bright reference outline of the undeformed rest shape -- exactly what Reset
          // brings the material back to, so you can always see how far it has drifted
          ctx.setLineDash([3,3]);
          ctx.strokeStyle = '#facc15'; ctx.lineWidth = 1.5;
          ctx.strokeRect(BLOCK_X0, midY - BLOCK_H/2, SPRING_X1 - BLOCK_X0, BLOCK_H);
          ctx.setLineDash([]);
          const restHandleX = rightX + 18;
          const handleX = dragging && cursorPx !== null ? cursorPx : restHandleX;
          const wallHalf = BLOCK_H/2 + 8;
          ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 2;
          ctx.beginPath(); ctx.moveTo(WALL_X, midY-wallHalf); ctx.lineTo(WALL_X, midY+wallHalf); ctx.stroke();
          const nhatch = 9;
          for(let i=0;i<nhatch;i++){
            const y0 = midY - wallHalf + (2*wallHalf)*i/(nhatch-1);
            ctx.beginPath(); ctx.moveTo(WALL_X, y0); ctx.lineTo(WALL_X-8, y0+7); ctx.stroke();
          }
          drawMaxwellElement(midY);
          // the deforming element itself: a lattice of small "fibers" (springs) running
          // both directions, filling the rod's own rectangular volume -- horizontal fibers
          // stretch with axial strain, vertical fibers compress with the lateral response,
          // felt together as one continuum rather than as separate external springs. The
          // mesh itself is fixed now (stiffness lives in the spring/dashpot element above);
          // only the bracing pattern below is read live from the scrub selection.
          const xs = []; for(let i=0;i<MESH_NX;i++) xs.push(leftX + (rightX-leftX)*i/(MESH_NX-1));
          const ys = []; for(let j=0;j<MESH_NY;j++) ys.push(topY + (botY-topY)*j/(MESH_NY-1));
          // fiber color mirrors the dashpot's own viscosity color scale -- the whole
          // continuum visibly agrees with the dashpot about how fluid-like it currently is.
          // Coil amplitude/frequency echo the macroscopic spring's own E-dependence (see
          // drawMaxwellElement above) -- scrubbing stiffness retunes the whole lattice, not
          // just the one external spring.
          const fiberColor = viscosityColor(currentEtaFrac());
          const meshEFrac = (state.E-E_BOUNDS[0])/(E_BOUNDS[1]-E_BOUNDS[0]);
          const meshAmp = MESH_SPRING_AMP * (1.3 - 0.55*meshEFrac), meshFreq = MESH_SPRING_FREQ * (0.7 + 0.5*meshEFrac);
          // each horizontal fiber is its own tiny Maxwell element -- a short spring in
          // series with a miniature copy of the same dashpot icon used above -- since a
          // horizontal fiber IS the axial spring+dashpot path, just felt at fiber scale
          // instead of as one external element
          for(const y of ys){
            for(let i=0;i<MESH_NX-1;i++){
              const segX0 = xs[i], segX1 = xs[i+1], segLen = segX1-segX0;
              const springEnd = segX0 + segLen*0.6;
              ctx.strokeStyle = fiberColor; ctx.lineWidth = 1;
              drawSpring(segX0, springEnd, y, meshAmp, meshFreq*0.6);
              drawDashpotIcon(ctx, springEnd + segLen*0.05, segX1, y, Math.min(3.5, meshAmp*0.9+1), fiberColor, false);
            }
          }
          ctx.strokeStyle = fiberColor; ctx.lineWidth = 1;
          for(const x of xs){
            for(let j=0;j<MESH_NY-1;j++) drawSpringV(x, ys[j], ys[j+1], meshAmp, meshFreq);
          }
          // diagonal braces -- the whole point of the bracing scrub: these are what
          // couple horizontal stretching to vertical motion in the first place (see the
          // aside). All three levels are genuinely braced (no "zero" option -- see the
          // aside for why an unbraced grid isn't offered as a real material); they differ
          // in how many diagonals per cell and how boldly drawn.
          ctx.strokeStyle = state.bracingIdx >= 2 ? 'rgba(229,231,235,0.85)' : 'rgba(229,231,235,0.45)';
          ctx.lineWidth = state.bracingIdx >= 2 ? 1.1 : 0.7;
          for(let i=0;i<MESH_NX-1;i++){
            for(let j=0;j<MESH_NY-1;j++){
              ctx.beginPath(); ctx.moveTo(xs[i],ys[j]); ctx.lineTo(xs[i+1],ys[j+1]); ctx.stroke();
              if(state.bracingIdx >= 1){
                ctx.beginPath(); ctx.moveTo(xs[i+1],ys[j]); ctx.lineTo(xs[i],ys[j+1]); ctx.stroke();
              }
            }
          }
          const isYielding = !releasing && Math.abs(state.stress) >= state.sigma_y - 1e-6;
          const nodeColor = isYielding ? '#ef4444' : (releasing ? '#4ade80' : '#38bdf8');
          ctx.fillStyle = nodeColor;
          for(const x of xs){ for(const y of ys){ ctx.beginPath(); ctx.arc(x,y,2.5,0,2*Math.PI); ctx.fill(); } }
          ctx.strokeStyle = '#4b5563'; ctx.lineWidth = 1;
          ctx.strokeRect(leftX, topY, rightX-leftX, botY-topY);
          // cord from the element's right edge to the handle -- its length visibly shows
          // strain running ahead of (or staying with) the applied pull, since the two are
          // no longer the same thing
          ctx.strokeStyle = dragging ? '#f3f4f6' : '#9ca3af'; ctx.lineWidth = 2;
          ctx.beginPath(); ctx.moveTo(rightX, midY); ctx.lineTo(handleX, midY); ctx.stroke();
          ctx.beginPath(); ctx.arc(handleX, midY, 9, 0, 2*Math.PI);
          ctx.fillStyle = dragging ? '#f3f4f6' : '#9ca3af'; ctx.fill();
          if(dragging){
            // how long THIS pull has been held, in SIMULATED time (see dragSimTime) --
            // the whole point of eta/creep is that this number, not just the force,
            // decides how much the material remembers; the applied stress alongside it is
            // what's actually driving that response
            const realHeld = elapsed()-dragStartTime;
            ctx.fillStyle = '#facc15'; ctx.font = 'bold 12px sans-serif'; ctx.textAlign = 'center';
            ctx.fillText(formatDuration(dragSimTime(realHeld)) + '  (' + state.stress.toFixed(1) + ' MPa)', handleX, midY-18);
          }
          ctx.fillStyle = '#9ca3af'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
          ctx.fillText(isYielding ? 'slipping' : (releasing ? 'releasing…' : 'elastic'), W/2, H-10);
        }

        function graphX(eps, eMax, W){ return 30 + (eps+eMax)/(2*eMax) * (W-45); }
        function graphY(sig, sMax, H, padT, padB){ return padT + (1 - (sig+sMax)/(2*sMax)) * (H-padT-padB); }

        // smart axis bounds for the stress-strain plot: scaled to the CURRENT material's
        // own sigma_y/E (with a small floor -- GRAPH_STRESS_MAX/GRAPH_STRAIN_MAX -- so a
        // near-zero sigma_y doesn't zoom the axis in on noise) rather than one fixed
        // constant, so a soft material's modest pull doesn't collapse to a flat line
        // against an oversized axis. Recomputed from state, not per-drag data, so the axis
        // stays stable WHILE dragging -- but still widens if actual plastic/creep strain
        // ever outgrows the elastic-scale guess (the `for` loops below).
        function graphBounds(){
          let sMax = Math.max(GRAPH_STRESS_MAX, state.sigma_y * 1.6);
          const elasticAtYield = (state.sigma_y*1.0e6) / (state.E*1.0e9);
          let eMax = Math.max(GRAPH_STRAIN_MAX, elasticAtYield * 3.5);
          for(const v of state.stressHistory) sMax = Math.max(sMax, Math.abs(v)*1.1);
          for(const v of displayStrain) eMax = Math.max(eMax, Math.abs(v)*1.1);
          return { sMax, eMax };
        }

        // the manual pull-release stress-strain trace -- ALWAYS shows the hand-pull
        // history (unaffected by preset clicks, which don't populate stressHistory);
        // unchanged behavior from before, just resized to the narrower graph column
        function drawGraph(){
          const ctx = graphCtx, W = PW_GRAPH, H = PH_GRAPH;
          ctx.clearRect(0,0,W,H); // missing this left the LAST drawn trace on screen forever
          // after Reset (nothing left to draw over it, so the old pixels just sat there --
          // drawTimeSeries already clears, which is why only this panel looked "stuck")
          const padL=34, padT=12, padB=20;
          const { sMax, eMax } = graphBounds();
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
          ctx.beginPath(); ctx.moveTo(padL, padT); ctx.lineTo(padL, H-padB); ctx.lineTo(W-8, H-padB); ctx.stroke();
          const zeroY = graphY(0, sMax, H, padT, padB);
          ctx.strokeStyle = '#2f3744';
          ctx.beginPath(); ctx.moveTo(padL, zeroY); ctx.lineTo(W-8, zeroY); ctx.stroke();
          // axis tick markers: extremes + zero on each axis, reflecting the smart bounds.
          // Alignment/offsets below are deliberately per-position, not a uniform 'center' +
          // fixed offset -- a naive version centers the leftmost x-label directly under the
          // y-axis's own bottom label, and the two collide (same corner, near-identical y).
          ctx.fillStyle = '#6b7280'; ctx.font = '9px sans-serif';
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
          const xTickVals = [-eMax, 0, eMax];
          xTickVals.forEach((ev, i) => {
            const x = graphX(ev, eMax, W);
            ctx.beginPath(); ctx.moveTo(x, H-padB); ctx.lineTo(x, H-padB+3); ctx.stroke();
            // leftmost label reads rightward from its tick, rightmost reads leftward --
            // both point away from the corner they'd otherwise collide with
            ctx.textAlign = i === 0 ? 'left' : (i === xTickVals.length-1 ? 'right' : 'center');
            ctx.fillText(ev === 0 ? '0' : ev.toExponential(1), x, H-padB+13);
          });
          const yTickVals = [sMax, 0, -sMax];
          ctx.textAlign = 'right';
          yTickVals.forEach((sv, i) => {
            const y = graphY(sv, sMax, H, padT, padB);
            ctx.beginPath(); ctx.moveTo(padL-3, y); ctx.lineTo(padL, y); ctx.stroke();
            // top label nudged down (it starts right at the canvas edge), bottom label
            // nudged up (away from the x-axis tick labels drawn just below it)
            const yOff = i === 0 ? 8 : (i === yTickVals.length-1 ? -2 : 3);
            ctx.fillText(sv === 0 ? '0' : sv.toFixed(0), padL-5, y+yOff);
          });
          if(displayStrain.length >= 2){
            const n = Math.min(displayStrain.length, state.stressHistory.length);
            const OLD = [100,116,139], NEW = [239,68,68];
            const GRAD_WINDOW = 150;
            for(let i=1;i<n;i++){
              const age = (n-1) - i;
              const t = Math.max(0, 1 - age/GRAD_WINDOW);
              const r = Math.round(OLD[0]+(NEW[0]-OLD[0])*t);
              const g = Math.round(OLD[1]+(NEW[1]-OLD[1])*t);
              const b = Math.round(OLD[2]+(NEW[2]-OLD[2])*t);
              const a = (0.3 + 0.7*t).toFixed(2);
              ctx.strokeStyle = 'rgba('+r+','+g+','+b+','+a+')';
              // capped smaller than before -- a thick stroke over a near-zero-length
              // segment (many samples land almost on top of each other during a slow or
              // held pull) renders as a fat round blob, not a thin line, which is what
              // made the trace look like a chain of beads instead of a curve
              ctx.lineWidth = 0.8 + 0.7*t;
              ctx.beginPath();
              ctx.moveTo(graphX(displayStrain[i-1], eMax, W), graphY(state.stressHistory[i-1], sMax, H, padT, padB));
              ctx.lineTo(graphX(displayStrain[i], eMax, W), graphY(state.stressHistory[i], sMax, H, padT, padB));
              ctx.stroke();
            }
            const lastX = graphX(displayStrain[n-1], eMax, W), lastY = graphY(state.stressHistory[n-1], sMax, H, padT, padB);
            ctx.beginPath(); ctx.arc(lastX, lastY, 2.5, 0, 2*Math.PI);
            ctx.fillStyle = '#ef4444'; ctx.strokeStyle = '#f3f4f6'; ctx.lineWidth = 0.8;
            ctx.fill(); ctx.stroke();
          }
          ctx.fillStyle = '#9ca3af'; ctx.font = '10px sans-serif'; ctx.textAlign='center';
          ctx.fillText('strain (MPa/GPa)', W/2, H-4);
          ctx.save(); ctx.translate(9, H/2); ctx.rotate(-Math.PI/2);
          ctx.fillText('stress (MPa)', 0, 0); ctx.restore();
        }

        // dual-axis time series: stress(t) on an implicit left scale, axial AND lateral
        // strain(t) sharing an implicit right scale (they're the same units, directly
        // comparable to each other, unlike stress) -- no numeric tick labels on either
        // axis, matching this widget's existing "colored legend, self-normalized" style;
        // reads either the live pull-release history or a preset's scripted response,
        // whichever state.graphMode currently points at (see applyPreset/graphRevertTimer)
        function drawTimeSeries(){
          const ctx = tsCtx, W = PW_GRAPH, H = PH_TS;
          ctx.clearRect(0,0,W,H);
          const padL=8, padR=8, padT=26, padB=24;
          const plotW = W-padL-padR, plotH = H-padT-padB;

          // click any legend label to hide/show that trace (see tsCv 'click' listener
          // below) -- boxes measured fresh each draw so hit-testing always matches what's
          // actually on screen; hidden labels dim to gray instead of disappearing, so
          // there's still something to click to bring them back
          ctx.font='10px sans-serif'; ctx.textAlign='left';
          let lx = padL;
          for(const [key, label, color] of [['stress','stress','#3b82f6'], ['axial','axial strain','#ef4444'], ['lateral','lateral strain','#f59e0b']]){
            const w = ctx.measureText(label).width;
            ctx.fillStyle = tsVisible[key] ? color : '#4b5563';
            ctx.fillText(label, lx, 11);
            tsLegendBoxes[key] = { x0: lx-2, x1: lx+w+2, y0: 1, y1: 14 };
            lx += w + 10;
          }

          const zeroY = padT + plotH/2;
          function drawFrame(){
            ctx.strokeStyle='#374151'; ctx.lineWidth=1;
            ctx.beginPath(); ctx.moveTo(padL,padT); ctx.lineTo(padL,H-padB); ctx.lineTo(W-padR,H-padB); ctx.stroke();
            ctx.strokeStyle='#2f3744';
            ctx.beginPath(); ctx.moveTo(padL,zeroY); ctx.lineTo(W-padR,zeroY); ctx.stroke();
          }
          // round geological reference points (matching the tau_wave/tau_M timescale strip
          // above the block) rather than the exact data endpoints -- "16.3 day" reads as
          // an arbitrary number, "10 day"/"1 yr" reads as a timescale you can compare to
          const yrTS = 365.25*24*3600;
          const TS_REF_TICKS = [
            [1,'1 s'], [60,'1 min'], [3600,'1 hr'], [86400,'1 day'],
            [yrTS,'1 yr'], [10*yrTS,'10 yr'], [100*yrTS,'100 yr'], [1.0e3*yrTS,'1 kyr'],
            [1.0e4*yrTS,'10 kyr'], [1.0e5*yrTS,'100 kyr'], [1.0e6*yrTS,'1 Myr'],
            [1.0e7*yrTS,'10 Myr'], [1.0e8*yrTS,'100 Myr'],
          ];
          function drawTimeTicks(xOfRel, tMinPos, tMaxPos){
            ctx.fillStyle = '#6b7280'; ctx.font='9px sans-serif'; ctx.textAlign='center';
            ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
            let picks = TS_REF_TICKS.filter(([s]) => s >= tMinPos*0.9 && s <= tMaxPos*1.1);
            if(picks.length > 4){
              // too many round ticks would crowd this narrow panel -- thin to ~4, evenly
              // spaced through the filtered (already log-ordered) list
              const step = (picks.length-1)/3;
              picks = [0,1,2,3].map(i => picks[Math.round(i*step)]);
            } else if(picks.length === 0){
              // the pull never reached even one round reference point (e.g. a sub-second
              // hold) -- fall back to the trace's own endpoints so the axis isn't empty
              picks = [[tMinPos, formatDuration(tMinPos)], [tMaxPos, formatDuration(tMaxPos)]];
            }
            for(const [s,label] of picks){
              const x = xOfRel(s);
              ctx.beginPath(); ctx.moveTo(x, H-padB); ctx.lineTo(x, H-padB+3); ctx.stroke();
              ctx.fillText(label, x, H-padB+11);
            }
          }
          // stress and strain share one implicit vertical axis (independently normalized),
          // so a single set of numeric y-ticks would be ambiguous -- fold each trace's own
          // range into the bottom caption line instead of a separate row (there isn't
          // vertical room for one at this panel's height)
          function scaleLegendText(sMax, eMax){
            return '±' + sMax.toFixed(sMax<10?1:0) + ' MPa · ±' + eMax.toExponential(1) + ' strain';
          }

          if(state.graphMode === 'preset' && pushedAnim){
            const times = pushedAnim.animTime, stresses = pushedAnim.animStress;
            const axial = pushedAnim.animStrain, lateral = pushedAnim.animLateral;
            if(times.length < 2){
              ctx.fillStyle = '#6b7280'; ctx.font='11px sans-serif'; ctx.textAlign='center';
              ctx.fillText('pull the block to see its history here', W/2, H/2);
              return;
            }
            let tMinPos = Infinity;
            for(const t of times) if(t > 0 && t < tMinPos) tMinPos = t;
            if(!isFinite(tMinPos)) tMinPos = 1e-2;
            const tMaxPos = Math.max(tMinPos * 10, times[times.length-1]);
            const logMin = Math.log10(tMinPos), logMax = Math.log10(tMaxPos);
            const logSpan = Math.max(1e-9, logMax - logMin);
            const xOf = i => {
              const t = times[i];
              if(t <= 0) return padL;
              const frac = (Math.log10(t) - logMin) / logSpan;
              return padL + Math.max(0, Math.min(1, frac)) * plotW;
            };
            // preset times already start at 0, so "relative to cycle start" is just t itself
            const xOfRel = t => t <= 0 ? padL : padL + Math.max(0, Math.min(1, (Math.log10(t)-logMin)/logSpan)) * plotW;
            drawFrame();

            let sMax=1e-30; for(const v of stresses) sMax=Math.max(sMax, Math.abs(v));
            let eMax=1e-30; for(const v of axial) eMax=Math.max(eMax, Math.abs(v));
            for(const v of lateral) eMax=Math.max(eMax, Math.abs(v));

            function drawTrace(arr, maxAbs, color){
              ctx.strokeStyle=color; ctx.lineWidth=1.6;
              ctx.beginPath();
              for(let i=0;i<arr.length;i++){
                const x=xOf(i), y=zeroY - (arr[i]/maxAbs)*(plotH/2-2);
                i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
              }
              ctx.stroke();
            }
            if(tsVisible.stress) drawTrace(stresses, sMax, '#3b82f6');
            if(tsVisible.axial) drawTrace(axial, eMax, '#ef4444');
            if(tsVisible.lateral) drawTrace(lateral, eMax, '#f59e0b');

            // an instantaneous elastic jump only exists for a preset's step-and-hold
            // scenarios -- mark it distinctly from the creep that follows (only meaningful
            // while the axial-strain trace it points at is actually shown)
            if(tsVisible.axial){
              const jy = zeroY - (axial[0]/eMax)*(plotH/2-2);
              ctx.setLineDash([2,2]); ctx.strokeStyle='#6b7280'; ctx.lineWidth=1;
              ctx.beginPath(); ctx.moveTo(xOf(0), H-padB); ctx.lineTo(xOf(0), jy); ctx.stroke();
              ctx.setLineDash([]);
              ctx.fillStyle='#6b7280'; ctx.font='9px sans-serif'; ctx.textAlign='left';
              ctx.fillText('elastic jump', xOf(0)+3, jy+9);
            }

            drawTimeTicks(xOfRel, tMinPos, tMaxPos);
            const tMax = times[times.length-1] - times[0];
            ctx.fillStyle = '#9ca3af'; ctx.font='9px sans-serif'; ctx.textAlign='center';
            ctx.fillText('t: 0 to ' + formatDuration(tMax) + '  ·  ' + scaleLegendText(sMax, eMax), W/2, H-4);
            return;
          }

          // live mode: each pull is its own cycle (see cycleStarts above) -- draw each as
          // its own polyline instead of one line connecting cycle boundaries, which jumps
          // backward in (log) time and zigzags. Recent cycles opaque, older ones faded --
          // a slow gradient so several recent pulls stay comparable, not just the latest.
          const times = state.timeHistory, stresses = state.stressHistory;
          const axial = displayStrain, lateral = displayLateral;
          if(times.length < 2){
            ctx.fillStyle = '#6b7280'; ctx.font='11px sans-serif'; ctx.textAlign='center';
            ctx.fillText('pull the block to see its history here', W/2, H/2);
            return;
          }
          const cycles = [];
          for(let c=0;c<cycleStarts.length;c++){
            const s = cycleStarts[c], e = (c+1<cycleStarts.length) ? cycleStarts[c+1] : times.length-1;
            if(e > s) cycles.push([s,e]);
          }
          if(cycles.length === 0){
            ctx.fillStyle = '#6b7280'; ctx.font='11px sans-serif'; ctx.textAlign='center';
            ctx.fillText('pull the block to see its history here', W/2, H/2);
            return;
          }
          // shared log-time axis across cycles, each measured relative to ITS OWN start,
          // so cycles of very different real durations still overlay comparably
          let tMinPos = Infinity, tMaxPos = 0;
          for(const [s,e] of cycles){
            const t0 = times[s];
            for(let i=s;i<=e;i++){
              const rel = times[i]-t0;
              if(rel > 0 && rel < tMinPos) tMinPos = rel;
              if(rel > tMaxPos) tMaxPos = rel;
            }
          }
          if(!isFinite(tMinPos)) tMinPos = 1e-2;
          tMaxPos = Math.max(tMaxPos, tMinPos*10);
          const logMin = Math.log10(tMinPos), logMax = Math.log10(tMaxPos);
          const logSpan = Math.max(1e-9, logMax - logMin);
          const xOfRel = rel => rel <= 0 ? padL : padL + Math.max(0, Math.min(1, (Math.log10(rel)-logMin)/logSpan)) * plotW;
          drawFrame();

          let sMax=1e-30; for(const v of stresses) sMax=Math.max(sMax, Math.abs(v));
          let eMax=1e-30; for(const v of axial) eMax=Math.max(eMax, Math.abs(v));
          for(const v of lateral) eMax=Math.max(eMax, Math.abs(v));

          const nCycles = cycles.length;
          cycles.forEach(([s,e], idx) => {
            const t0 = times[s];
            ctx.globalAlpha = nCycles === 1 ? 1 : Math.max(0.12, 1 - 0.15*(nCycles-1-idx));
            function drawSeg(arr, maxAbs, color){
              ctx.strokeStyle=color; ctx.lineWidth=1.6;
              ctx.beginPath();
              for(let i=s;i<=e;i++){
                const x = xOfRel(times[i]-t0), y = zeroY - (arr[i]/maxAbs)*(plotH/2-2);
                i===s?ctx.moveTo(x,y):ctx.lineTo(x,y);
              }
              ctx.stroke();
            }
            if(tsVisible.stress) drawSeg(stresses, sMax, '#3b82f6');
            if(tsVisible.axial) drawSeg(axial, eMax, '#ef4444');
            if(tsVisible.lateral) drawSeg(lateral, eMax, '#f59e0b');
          });
          ctx.globalAlpha = 1;

          drawTimeTicks(xOfRel, tMinPos, tMaxPos);
          ctx.fillStyle = '#9ca3af'; ctx.font='9px sans-serif'; ctx.textAlign='center';
          ctx.fillText((nCycles>1 ? nCycles+' pulls, ' : '') + 'each 0 to ' + formatDuration(tMaxPos) + '  ·  ' + scaleLegendText(sMax, eMax), W/2, H-4);
        }

        // click a Time Series legend label to hide/show that trace; hovering one shows a
        // pointer cursor so the affordance is discoverable without a hint label eating
        // space in this already-tight panel
        function tsLegendKeyAt(px, py){
          for(const key of ['stress','axial','lateral']){
            const b = tsLegendBoxes[key];
            if(b && px>=b.x0 && px<=b.x1 && py>=b.y0 && py<=b.y1) return key;
          }
          return null;
        }
        tsCv.addEventListener('click', ev => {
          const r = tsCv.getBoundingClientRect();
          const key = tsLegendKeyAt(ev.clientX-r.left, ev.clientY-r.top);
          if(key){ tsVisible[key] = !tsVisible[key]; draw(); }
        });
        tsCv.addEventListener('mousemove', ev => {
          const r = tsCv.getBoundingClientRect();
          tsCv.style.cursor = tsLegendKeyAt(ev.clientX-r.left, ev.clientY-r.top) ? 'pointer' : 'default';
        });

        // log-scale timescale number line: a fixed axis from 0.01 s to ~300 Myr (the same
        // TS_MIN_LOG/TS_MAX_LOG range dragSimTime uses above), with reference ticks at
        // round timescales, and colored markers for tau_wave, tau_M, and (once there's one
        // to show) the current comparison period T -- a single glance at relative POSITION
        // says more than three separately-read numbers do
        function tsX(seconds, padL, padR){
          const lg = Math.log10(Math.max(seconds, 1e-2));
          const frac = (lg-TS_MIN_LOG)/(TS_MAX_LOG-TS_MIN_LOG);
          return padL + Math.max(0, Math.min(1, frac)) * (PW-padL-padR);
        }
        function drawTimescale(){
          const ctx = timescaleCtx, W = PW, H = TS_H;
          ctx.clearRect(0,0,W,H);
          const padL = 8, padR = 8, axisY = H*0.52;
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
          ctx.beginPath(); ctx.moveTo(padL, axisY); ctx.lineTo(W-padR, axisY); ctx.stroke();
          const yr = 365.25*24*3600;
          const refs = [[1,'s'], [3600,'hr'], [yr,'yr'], [1e3*yr,'kyr'], [1e6*yr,'Myr']];
          ctx.fillStyle = '#6b7280'; ctx.font = '9px sans-serif'; ctx.textAlign = 'center';
          for(const [s,label] of refs){
            const x = tsX(s, padL, padR);
            ctx.strokeStyle = '#374151';
            ctx.beginPath(); ctx.moveTo(x, axisY-3); ctx.lineTo(x, axisY+3); ctx.stroke();
            ctx.fillText(label, x, axisY+15);
          }
          const [tauWaveS, tauMS] = waveRelaxationClocks(state.E, state.nu, state.eta, state.L, state.rho);
          function marker(seconds, color, label, above){
            const x = tsX(seconds, padL, padR);
            const y2 = above ? axisY-16 : axisY+16;
            ctx.strokeStyle = color; ctx.lineWidth = 2;
            ctx.beginPath(); ctx.moveTo(x, axisY); ctx.lineTo(x, y2); ctx.stroke();
            ctx.beginPath(); ctx.arc(x, y2, 3.5, 0, 2*Math.PI); ctx.fillStyle = color; ctx.fill();
            ctx.font = 'bold 10px sans-serif'; ctx.textAlign = 'center';
            ctx.fillText(label, x, above ? y2-6 : y2+13);
          }
          marker(tauWaveS, '#38bdf8', 'τ_wave', true);
          marker(tauMS, '#f59e0b', 'τ_M', true);
          if(state.currentT != null) marker(state.currentT, '#facc15', 'T', false);
          return [tauWaveS, tauMS];
        }
        function drawClocksPanel(){
          const [tauWaveS, tauMS] = drawTimescale();
          if(state.currentT == null){
            clocksTextEl.innerHTML = '&tau;_wave = <b>' + formatDuration(tauWaveS) + '</b> &middot; &tau;_M = <b>'
              + formatDuration(tauMS) + '</b> &middot; pull the block or click a preset to compare your own timescale';
            return;
          }
          const {Rwave, De, regime} = deborahRegime(tauWaveS, tauMS, state.currentT);
          clocksTextEl.innerHTML = 'T (' + state.currentTLabel + ') = <b>' + formatDuration(state.currentT) + '</b>'
            + ' &middot; &tau;_wave/T = ' + Rwave.toExponential(2) + ' &middot; &tau;_M/T = ' + De.toExponential(2)
            + ' &middot; regime = <b>' + regime + '</b>';
        }

        function drawReadouts(){
          readoutEl.textContent = 'stress = ' + state.stress.toFixed(2) + ' MPa'
            + '  ·  strain = ' + state.strain.toExponential(3) + '  ·  lateral = ' + state.lateralStrain.toExponential(3);
        }

        // flags the sigma_y slider's own label so the user has a persistent, in-place cue
        // that yielding happened -- not just the transient "slipping" caption on the scene
        // canvas, which disappears the moment the pull eases off
        function updateYieldLabel(){
          const yieldingNow = !releasing && Math.abs(state.stress) >= state.sigma_y - 1e-6;
          const everYielded = Math.abs(state.epsP) > 1e-9;
          syLabelEl.textContent = 'σ_y (yield, MPa)' + (yieldingNow ? '  — yielding now' : (everYielded ? '  — yielded' : ''));
          syLabelEl.style.color = yieldingNow ? '#ef4444' : (everYielded ? '#f59e0b' : '');
        }
        function draw(){
          drawScene(); drawClocksPanel(); drawReadouts(); drawGraph(); drawTimeSeries();
          updateYieldLabel();
        }

        function tsTitleFor(mode){
          return mode === 'preset' ? 'Time Series -- Scripted Scenario' : 'Time Series -- Your Pull History';
        }
        function tsCaptionFor(mode){
          if(mode === 'preset') return 'stress (left scale) · axial & lateral strain (right scale) vs. time, for the preset that just fired';
          return 'stress (left scale, MPa) · axial & lateral strain (right scale) vs. time, from your own hand-pull';
        }
        function setPresetButtons(){
          par.querySelectorAll('.epl-preset-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.preset === state.presetKey);
          });
        }
        function modeDescFor(mode){
          if(mode === 'strain') return "you're setting the DISPLACEMENT -- hold the handle steady (or let go) and simulated time advances while your strain stays fixed, watching the derived stress relax";
          return "you're setting the FORCE -- hold the handle steady and simulated time advances while your stress stays fixed, watching the derived strain creep";
        }
        function setModeButtons(){
          par.querySelectorAll('.epl-mode-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.mode === state.controlMode);
          });
          modeDescEl.textContent = modeDescFor(state.controlMode);
        }
        function syncControls(){
          const syEl = par.querySelector('#epl-sy');
          syEl.value = state.sigma_y;
          par.querySelector('#epl-sy-v').textContent = state.sigma_y.toFixed(0);
          par.querySelector('#epl-H-v').textContent = state.H.toFixed(0);
          par.querySelector('#epl-bracing-v').textContent = 'ν=' + state.nu.toFixed(3);
          drawScrubStrip(bracingCtx, bracingScrubW, bracingPresets, state.bracingIdx);
          materialDescEl.textContent = 'E: ' + describeStiffness(state.E) + '  ·  η: ' + describeViscosity(state.eta);
          tsTitleEl.textContent = tsTitleFor(state.graphMode);
          tsCaptionEl.textContent = tsCaptionFor(state.graphMode);
          setPresetButtons();
          setModeButtons();
        }
        function resetMaterial(){
          // any parameter change invalidates the current deformation -- it was computed
          // under the OLD material, so a strain left over from before would no longer
          // mean anything for what's now selected. This is also exactly what the Reset
          // button itself does. Also reverts the graph/clocks panel to the plain
          // free-play baseline -- callers that just fired a preset set graphMode/currentT
          // again immediately afterward.
          releasing = false;
          dragging = false;
          cursorPx = null;
          animPlaying = false;
          if(graphRevertTimer){ clearTimeout(graphRevertTimer); graphRevertTimer = null; }
          state.stress = 0; state.strain = 0; state.epsP = 0; state.lateralStrain = 0;
          state.stressHistory = [0];
          state.timeHistory = [0]; // simulated seconds, not wall-clock -- see dragSimTime
          cycleStarts = [0];
          displayStrain = [0];
          displayLateral = [0];
          state.graphMode = 'live';
          state.currentT = null; state.currentTLabel = '';
          epoch++;
        }

        function onSlider(event){
          const id = event.target.id;
          if(id === 'epl-sy') state.sigma_y = Number(event.target.value);
          else if(id === 'epl-H') state.H = Number(event.target.value);
          else return;
          resetMaterial();
          syncControls();
          draw();
          throttledEmit();
        }
        par.querySelectorAll('input[type=range]').forEach(el => el.addEventListener('input', onSlider));

        // scrub bracing instead of dialing a slider -- click or drag across the strip
        // snaps to whichever of the 3 discrete zones the pointer is over
        let draggingBracing = false;
        function zoneFromEvent(cv, w, ev){
          const r = cv.getBoundingClientRect();
          const x = ev.clientX - r.left;
          return Math.max(0, Math.min(2, Math.floor(3*x/w)));
        }
        function selectBracing(idx){
          if(idx === state.bracingIdx) return;
          state.bracingIdx = idx;
          state.nu = bracingPresets[idx].nu;
          resetMaterial();
          syncControls(); draw(); emit();
        }
        bracingCv.addEventListener('mousedown', ev => {
          draggingBracing = true;
          selectBracing(zoneFromEvent(bracingCv, bracingScrubW, ev));
        });
        window.addEventListener('mousemove', ev => {
          if(draggingBracing) selectBracing(zoneFromEvent(bracingCv, bracingScrubW, ev));
        });
        window.addEventListener('mouseup', () => { draggingBracing = false; });

        par.querySelector('#epl-reset').addEventListener('click', () => {
          state.presetKey = null; // no preset is "active" after a plain reset
          resetMaterial();
          tsVisible = { stress: true, axial: true, lateral: true }; // un-hide any legend toggles
          syncControls();
          draw();
          emit();
        });

        // clicking a preset: reset to that scenario's material, then wait for Julia's
        // freshly-computed scripted response before starting the one-shot playback (see
        // the epl-results handler and playbackTick below) -- clicking the SAME preset
        // again simply replays it, no special-casing needed
        function applyPreset(key){
          const p = REGIME_PRESETS[key];
          state.E = p.E; state.eta = p.eta; state.rho = p.rho; state.L = p.L;
          state.sigma_y = p.sigma_y; state.H = 0;
          state.bracingIdx = p.bracingIdx; state.nu = bracingPresets[p.bracingIdx].nu;
          state.controlMode = 'stress'; // every scripted regime preset is a force history
          state.presetKey = key;
          resetMaterial();
          state.graphMode = 'preset';
          state.currentT = p.duration;
          state.currentTLabel = p.name.toLowerCase();
          awaitingPresetAnim = key;
          syncControls(); draw(); emit();
        }
        par.querySelectorAll('.epl-preset-btn').forEach(btn => {
          btn.addEventListener('click', () => applyPreset(btn.dataset.preset));
        });

        // switching control mode discards the current deformation the same way a material
        // change does -- a strain left over from a force-controlled pull (or vice versa)
        // doesn't mean anything once the driven/derived roles have flipped
        function selectMode(mode){
          if(mode === state.controlMode) return;
          state.controlMode = mode;
          state.presetKey = null; // presets are scripted, force-controlled scenarios -- switching modes exits that context
          resetMaterial();
          syncControls(); draw(); emit();
        }
        par.querySelectorAll('.epl-mode-btn').forEach(btn => {
          btn.addEventListener('click', () => selectMode(btn.dataset.mode));
        });

        function scenePointerXY(ev){
          const r = sceneCv.getBoundingClientRect();
          return [ev.clientX - r.left, ev.clientY - r.top];
        }

        let dragAnchorPx = 0, dragAnchorStress = 0, dragAnchorStrain = 0;
        let dragStartTime = 0; // elapsed() at mousedown -- also becomes "your pull period" on release
        // bumped on any hard state invalidation (release settling, reset) so a response
        // computed for a request sent *before* that point -- which can easily still be
        // in flight -- gets recognized as stale and discarded when it finally arrives,
        // rather than silently clobbering the settled/reset state
        let epoch = 0;

        function startRelease(){
          releasing = true;
          epoch++;
          const fromSigma = state.stress;
          // the authoritative residual (plastic + creep) from Julia's last response --
          // NOT state.epsP, which is only the local elastic/plastic prediction and knows
          // nothing about creep, so it would ignore anything accumulated while holding
          const epsPFixed = pushed && typeof pushed.plasticStrain === 'number' ? pushed.plasticStrain : state.epsP;
          // freeze simulated time at whatever the hold reached -- releasing is a ~350ms UI
          // animation, not meant to represent additional geological time, so every sample
          // pushed during it reuses this SAME simulated instant (dt=0, no extra creep)
          const releaseSimT = dragSimTime(elapsed()-dragStartTime);
          const dur = 350;
          const t0r = performance.now();
          function step(now){
            if(!releasing) return; // a new drag cancelled this animation
            const f = Math.min(1, (now-t0r)/dur);
            if(f < 1){
              const eased = 1 - Math.pow(1-f, 3);
              const sigma = fromSigma * (1-eased);
              const strain = epsPFixed + (sigma*1e6)/(state.E*1e9);
              const lateral = localLateralStrain(strain, epsPFixed);
              const last = state.stressHistory[state.stressHistory.length-1];
              state.stress = sigma; state.strain = strain;
              state.lateralStrain = lateral;
              if(Math.abs(sigma-last) > 0.3){
                state.stressHistory.push(sigma); state.timeHistory.push(releaseSimT);
                displayStrain.push(strain); displayLateral.push(lateral);
              }
              draw();
              requestAnimationFrame(step);
            } else {
              // force the exact final point so stressHistory/timeHistory/displayStrain stay
              // aligned and the trace visibly lands on (epsPFixed, 0), not a threshold-skipped neighbor
              const finalLateral = localLateralStrain(epsPFixed, epsPFixed);
              state.stress = 0; state.strain = epsPFixed; state.epsP = epsPFixed;
              state.lateralStrain = finalLateral;
              state.stressHistory.push(0); state.timeHistory.push(releaseSimT);
              displayStrain.push(epsPFixed); displayLateral.push(finalLateral);
              releasing = false;
              // "your own period": treat the pull-and-release as half a cycle, so a
              // student can directly compare their own hand-timescale against tau_wave/tau_M
              state.currentT = Math.max(0.02, 2 * releaseSimT);
              state.currentTLabel = 'your last pull';
              draw();
              emit();
            }
          }
          requestAnimationFrame(step);
        }

        // relative-pixel-delta dragging: every control below (spring, dashpot) updates
        // its parameter from the MOVEMENT of the pointer since the last frame, not its
        // absolute position -- the same trick works for both linear (add scaled delta)
        // and log-scale (multiply by exp(scaled delta)) controls without needing a
        // separate fixed pixel-range mapping per control.
        let activeDrag = null; // 'handle' | 'spring' | 'dashpot' | null
        let dragLastX = 0;

        sceneCv.addEventListener('mousedown', ev => {
          const [px, py] = scenePointerXY(ev);
          const midY = PH/2;
          if(pointerOnHandle(px, py, midY)){
            activeDrag = 'handle';
            releasing = false;
            dragging = true;
            cursorPx = px;
            sceneCv.style.cursor = 'grabbing';
            dragAnchorPx = px;
            dragAnchorStress = state.stress;
            dragAnchorStrain = state.strain;
            dragStartTime = elapsed();
            // start a fresh time-series cycle here, unless the previous one never
            // actually grew (e.g. a click that didn't turn into a drag)
            const lastIdx = state.timeHistory.length - 1;
            if(lastIdx > cycleStarts[cycleStarts.length-1]) cycleStarts.push(lastIdx);
            startHoldSampling();
            draw();
            return;
          }
          if(pointerOnSpring(px, py, midY)){ activeDrag = 'spring'; dragLastX = ev.clientX; return; }
          if(pointerOnDashpot(px, py, midY)){ activeDrag = 'dashpot'; dragLastX = ev.clientX; return; }
        });
        sceneCv.addEventListener('mousemove', ev => {
          if(activeDrag) return; // window listener below handles the active drag
          const [px, py] = scenePointerXY(ev);
          const midY = PH/2;
          if(pointerOnHandle(px, py, midY)) sceneCv.style.cursor = 'grab';
          else if(pointerOnSpring(px,py,midY) || pointerOnDashpot(px,py,midY)) sceneCv.style.cursor = 'ew-resize';
          else sceneCv.style.cursor = 'default';
        });

        window.addEventListener('mousemove', ev => {
          if(activeDrag === 'handle'){
            if(!dragging) return;
            const [px] = scenePointerXY(ev);
            cursorPx = px;
            if(state.controlMode === 'strain'){
              pushStrain(dragAnchorStrain + (px - dragAnchorPx) / PX_PER_STRAIN_UNIT);
            } else {
              pushStress(dragAnchorStress + (px - dragAnchorPx) / PX_PER_STRESS_UNIT);
            }
            draw();
            throttledEmit();
            return;
          }
          if(activeDrag === 'spring'){
            const dx = ev.clientX - dragLastX; dragLastX = ev.clientX;
            const sens = (E_BOUNDS[1]-E_BOUNDS[0])/250;
            state.E = Math.max(E_BOUNDS[0], Math.min(E_BOUNDS[1], state.E + dx*sens));
            resetMaterial(); syncControls(); draw(); throttledEmit();
            return;
          }
          if(activeDrag === 'dashpot'){
            const dx = ev.clientX - dragLastX; dragLastX = ev.clientX;
            const v = state.eta * Math.exp(dx*0.0645); // log-scale: ~250px sweeps the full 1e17-1e24 range
            state.eta = Math.max(ETA_BOUNDS[0], Math.min(ETA_BOUNDS[1], v));
            resetMaterial(); syncControls(); draw(); throttledEmit();
          }
        });
        window.addEventListener('mouseup', () => {
          if(activeDrag === 'handle'){
            dragging = false;
            cursorPx = null;
            sceneCv.style.cursor = 'grab';
            stopHoldSampling();
            startRelease();
          }
          activeDrag = null;
        });

        par.addEventListener('epl-results', event => {
          commitInFlight = false;
          const detail = event.detail || null;
          // a response can arrive after a release/reset already invalidated the request
          // it answers (a stale in-flight round trip) -- the epoch tag is what actually
          // identifies that, not the response's shape or length. Only assign it to
          // `pushed` at all when it's current: startRelease() later trusts whatever
          // `pushed.plasticStrain` last held, so a rejected response must never even
          // land there, or it silently corrupts the *next* release's residual target.
          // Julia always echoes back BOTH strain and stress, whichever one is the driven
          // input (an exact echo, no correction needed) and whichever is derived (the
          // authoritative creep- or relaxation-corrected value) -- see the Appendix cell
          // that builds EplPush. Overwriting both uniformly here means this handler needs
          // no controlMode branch of its own: in 'stress' mode strain is what gets
          // corrected, in 'strain' mode stress is, and the other one just round-trips.
          if(detail && detail.strain && detail.stress && detail.epoch === epoch){
            pushed = detail;
            displayStrain = pushed.strain.slice();
            state.strain = displayStrain[displayStrain.length-1];
            state.stressHistory = pushed.stress.slice();
            state.stress = state.stressHistory[state.stressHistory.length-1];
            state.epsP = typeof pushed.plasticStrain === 'number' ? pushed.plasticStrain : state.epsP;
            if(pushed.lateral && pushed.lateral.length){
              displayLateral = pushed.lateral.slice();
              state.lateralStrain = displayLateral[displayLateral.length-1];
            }
          }
          // the preset's scripted response has no epoch/staleness concern -- it's a pure
          // function of the current material + presetKey, always safe to adopt. Only
          // actually START the one-shot playback if this arrival is the fresh data a
          // just-clicked preset is waiting for -- otherwise a slider tweak's incidental
          // recompute would silently kick off an unwanted animation.
          if(detail && detail.animTime){
            pushedAnim = {
              animTime: detail.animTime,
              animStress: detail.animStress.map(v => v/1.0e6), // Pa -> MPa, matching sigma_y's units
              animStrain: detail.animStrain,
              animLateral: detail.animLateral
            };
            if(awaitingPresetAnim){
              animPlaying = true;
              animStartMs = performance.now();
              awaitingPresetAnim = null;
            }
          }
          draw();
        });

        // one-shot preset playback: advances a SIMULATED time (not wall-clock -- a
        // preset's duration can be anywhere from decades to tens of thousands of years)
        // over a fixed ~4 real seconds regardless of the true duration, reading the
        // block's strain/lateral/stress off the precomputed response by nearest-sample
        // lookup, then STOPS (leaves the final frame showing) -- unlike the sine-wave
        // version this replaced, there is no reason to loop forever. A few seconds after
        // it finishes, the graph quietly reverts to the plain pull-release trace.
        //
        // If this cell's script ever re-executes (Pluto reruns a cell's <script> on every
        // output update, e.g. while a live-collab edit is still settling) while REUSING the
        // same canvas element, a naive rAF loop from the PREVIOUS execution keeps running
        // forever -- nothing ever cancels it -- so two (or more) loops end up racing to
        // paint the same canvas, each with its own stale closure over old drawing code.
        // Stashing the handle on `par` (which persists across reruns, unlike this script's
        // own local scope) lets each new execution cancel whatever its predecessor left
        // running, so there's only ever one live loop.
        if(par._eplRafHandle) cancelAnimationFrame(par._eplRafHandle);
        function playbackTick(nowMs){
          if(!par.isConnected) return; // this script instance's widget is gone; stop for good
          if(animPlaying && pushedAnim){
            const frac = Math.min(1, (nowMs - animStartMs)/ANIM_REAL_MS);
            const n = pushedAnim.animTime.length;
            const idx = Math.min(n-1, Math.floor(frac*n));
            state.stress = pushedAnim.animStress[idx];
            state.strain = pushedAnim.animStrain[idx];
            state.lateralStrain = pushedAnim.animLateral[idx];
            draw();
            if(frac >= 1){
              animPlaying = false;
              graphRevertTimer = setTimeout(() => {
                state.graphMode = 'live';
                syncControls(); draw();
              }, 4000);
            }
          }
          par._eplRafHandle = requestAnimationFrame(playbackTick);
        }
        par._eplRafHandle = requestAnimationFrame(playbackTick);

        syncControls();
        draw();
        } // end of the `if(!par._eplInitialized)` idempotency guard
        }
        </script>
        """)
    end

    const _epl_ready = true
end

# ╔═╡ 199e5d68-391c-4c20-871e-bf1222f28b3f
begin
    _epl_ready
    WideCell(@bind epl ElastoplasticLoadingInput(); max_width=1150)
end

# ╔═╡ df2b9bd9-f0b0-48db-af39-d1d9fbfbc9f4
begin
    struct EplPush
        strain::String
        stress::String
        plasticStrain::Float64
        lateral::String
        epoch::Int
        animTime::String
        animStress::String
        animStrain::String
        animLateral::String
    end
    function Base.show(io::IO, ::MIME"text/html", p::EplPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('eplwidget');
        if(w){
          w.dispatchEvent(new CustomEvent('epl-results', { detail: {
            strain: [$(p.strain)],
            stress: [$(p.stress)],
            plasticStrain: $(p.plasticStrain),
            lateral: [$(p.lateral)],
            epoch: $(p.epoch),
            animTime: [$(p.animTime)],
            animStress: [$(p.animStress)],
            animStrain: [$(p.animStrain)],
            animLateral: [$(p.animLateral)],
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 84911c1a-b218-40d2-b6af-0cff4180eef8
begin
    local stress_history = epl isa AbstractDict ? Float64.(epl["stressHistory"]) : [0.0]
    local strain_history_in = epl isa AbstractDict ? Float64.(get(epl, "strainHistory", [0.0])) : [0.0]
    local time_history = epl isa AbstractDict ? Float64.(epl["timeHistory"]) : [0.0]
    local control_mode = epl isa AbstractDict ? get(epl, "controlMode", "stress") : "stress"
    local E_val = epl isa AbstractDict ? epl["E"] : 75.0
    local sy_val = epl isa AbstractDict ? epl["sigma_y"] : 100.0
    local H_val = epl isa AbstractDict ? epl["H"] : 0.0
    local eta_val = epl isa AbstractDict ? epl["eta"] : 1.0e24
    local nu_val = epl isa AbstractDict ? epl["nu"] : 0.2412
    local epoch_val = epl isa AbstractDict ? Int(epl["epoch"]) : 0
    # convert from the widget's real-but-scaled units (MPa/GPa) to SI (Pa) before calling
    # the creep/relaxation integrators -- both work in raw Pa, and skipping this conversion
    # (as a previous version of this cell did) makes the elastic term ~1e3x too large and the
    # creep term ~1e6x too small, since sig/E and sig*dt/eta stop cancelling correctly.
    # scripted_forcing_response (below) does this same conversion internally.
    local strain_path, residual_path, stress_path_MPa
    if control_mode == "strain"
        # displacement-controlled: strain is the DRIVEN input here (echoed straight back
        # below -- see the widget's epl-results handler, which trusts it unconditionally
        # since nothing needs to correct it), and stress is what Julia derives via the
        # relaxation ODE. sigma_cap mirrors strain_from_stress's own plastic_cap*H_eff
        # bound, just to stop a wildly overdriven strain from reporting an unbounded stress
        # -- the function itself does not model plastic flow (see its docstring). E_val is
        # GPa, sy_val/H_val are MPa -- E_val*1000 converts it to the same MPa scale before
        # comparing/summing (mixing raw GPa and MPa numbers here was a real bug: it made
        # this cap ~1000x too tight, clamping barely above sigma_y instead of the generous,
        # rarely-hit ceiling this is meant to be).
        local sigma_cap_pa = (sy_val + max(H_val, 0.1 * E_val * 1000.0) * 1.5) * 1.0e6
        local stress_path_pa, eps_c_path = stress_from_strain_with_relaxation(
            strain_history_in, time_history, E_val * 1.0e9, eta_val; sigma_cap=sigma_cap_pa)
        strain_path = strain_history_in
        residual_path = eps_c_path
        stress_path_MPa = stress_path_pa ./ 1.0e6
    else
        strain_path, residual_path = strain_from_stress_with_creep(
            stress_history .* 1.0e6, time_history, E_val * 1.0e9, sy_val * 1.0e6, eta_val; H=H_val * 1.0e6)
        stress_path_MPa = stress_history
    end
    local lateral_path = lateral_strain_from_axial(strain_path, residual_path, nu_val)

    # the CURRENTLY SELECTED preset's scripted forcing (kind/amplitude/duration) applied to
    # the CURRENT material (E_val/eta_val/nu_val above, which may have been scrubbed away
    # from that preset's own defaults) -- computed unconditionally (cheap) regardless of
    # whether a preset animation is actually playing right now, so replaying or tweaking
    # the material never needs an extra round trip before the data is ready
    local preset_key_raw = epl isa AbstractDict ? epl["presetKey"] : "interseismic"
    local preset_key = preset_key_raw === nothing ? "interseismic" : preset_key_raw
    local preset = regime_presets[Symbol(preset_key)]
    local anim_t, anim_stress, anim_strain, anim_lateral = scripted_forcing_response(
        E_val, nu_val, sy_val, H_val, eta_val, preset.kind, preset.amplitude, preset.duration)

    EplPush(join(strain_path, ","), join(stress_path_MPa, ","), residual_path[end], join(lateral_path, ","),
        epoch_val, join(anim_t, ","), join(anim_stress, ","), join(anim_strain, ","), join(anim_lateral, ","))
end

# ╔═╡ 3851d568-3eaf-4a4f-b3a2-213a03958982
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.83"
"""

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
# ╟─46fa4dc1-b407-4584-91bf-d136b77622bc
# ╠═52597587-4c6b-4420-8242-6dbe166e280d
# ╟─faf8d654-5c42-46a5-967d-93c3897273e1
# ╟─6444aa15-39e9-4f0f-9e90-809899dd19ff
# ╟─199e5d68-391c-4c20-871e-bf1222f28b3f
# ╠═d3f8a2c1-7e5b-4c9a-8d16-2f4a6e8c0b31
# ╟─3e1a5c9a-8b7f-4a2f-9d1e-6f5c4b8e2d3a
# ╟─64f12034-4161-4d2b-b1c2-1809ee2e6e2c
# ╟─4a425525-dba3-4dcc-b9f6-7deaf4fb0588
# ╟─9b7df3b0-c0f0-4c96-8314-23842a2d97dd
# ╟─df3cccd7-5a20-452b-89e9-1f20b818808a
# ╟─6ab7bf3b-19e1-477a-b046-88df986855cb
# ╠═33fb6caf-2564-41bc-bd5f-8fba3670d887
# ╠═6cb2a1e4-9f61-4b8a-9c3e-7d5b6a2e4f10
# ╟─5e2c9a3f-1b6d-4f2a-9c0e-7a3f8d1b6e42
# ╠═7f4a1d8e-2c5b-4e9f-a0d3-6b1c8e4f2a97
# ╠═3d8b6f21-9e4a-4c7d-b502-1f6a9d3e8c05
# ╠═a0900587-c705-4fb9-99b6-002b4fe257c8
# ╠═f3a1b2c4-6d5e-4a8b-9c1f-2e7d4b6a8f01
# ╠═a7c3e9f1-4b2d-4e6a-8f19-3d5c7b9e1a02
# ╠═b8d4f0a2-5c3e-4f7b-9a2d-4e6d8c0f2b13
# ╠═c9e5a1b3-6d4f-4a8c-8b3e-5f7e9d1a3c24
# ╠═1a52d2bb-c6ac-4e0d-8a00-39ec56e0cdff
# ╠═3cf3b5c7-f651-4a0d-842e-ba54b5bb4c59
# ╠═c1d10e8c-c3e6-4ea5-8dad-d83b87c397a8
# ╠═9d2f94d6-6ff6-4f13-8a71-64eadf2c1a40
# ╠═6a1af493-35f5-4813-8940-1992c72cb1e5
# ╠═6cbc6bab-3e0a-415a-acd0-5ea78d7e0f36
# ╠═89cd485a-c305-44c0-a54c-24941c7f4e28
# ╠═95f0ad68-ee90-42b2-91a3-f2788271964f
# ╠═7d5499f8-dfeb-489b-8cfc-bd2d10886bf7
# ╠═1f8e58a0-c7a3-44f5-99b0-d5bbb9ee89b5
# ╠═0c92d271-1808-4109-b173-5dc4f52fa198
# ╟─84911c1a-b218-40d2-b6af-0cff4180eef8
# ╟─e76b20bd-28f7-418c-8fb3-42a40e710497
# ╠═02d7acba-4408-4158-b67d-76516d239309
# ╠═bdbf5dc7-ecd0-4a11-a3dd-a6d0a597cc72
# ╠═c104eee8-b783-47bf-8a56-7769a77fe426
# ╟─00961e05-3259-494e-8123-38675b04dfcc
# ╠═58011730-0491-400a-83a1-ebaf42c6b574
# ╠═df2b9bd9-f0b0-48db-af39-d1d9fbfbc9f4
# ╟─3851d568-3eaf-4a4f-b3a2-213a03958982
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
