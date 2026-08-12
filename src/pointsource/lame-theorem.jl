### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Lamé's Theorem"
#> layout = "layout.jlhtml"
#> tags = ["pointsource"]
#> description = "A point force radiates P and S waves in one shot. Drag the force direction, drag a receiver, and watch the near-field, far-field-P, and far-field-S terms separate."

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

# ╔═╡ 021492fb-0d1a-4f6e-a656-bdf6a8e74a90
begin
    using PlutoUI
    using LinearAlgebra
    using SpecialFunctions
end

# ╔═╡ cbc8f361-986c-424d-a649-6e841c072e72
TableOfContents()

# ╔═╡ dfaa5e67-7d3d-492a-8b56-a66bd0cd3970
md"""
# Lamé's Theorem

Hit an elastic whole-space with a single point force and it doesn't radiate one kind of wave —
it radiates two, a longitudinal P wave and a transverse S wave, from the same source at the same
instant. Lamé's theorem is the reason: any displacement field solving the elastic wave equation
splits, exactly, into a curl-free part and a divergence-free part, and those two parts turn out
to obey their *own* simple scalar/vector wave equations at speeds `` \alpha `` and `` \beta ``.
This notebook builds the complete time-domain solution for a single point force and hands you a
widget to explore it: drag the force direction, drag a receiver, and watch the near-field,
far-field-P, and far-field-S terms trade off as you move.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ f80221b6-eec6-4044-896b-68cca2bfcb6d
md"""
!!! note "Lamé's Theorem — Helmholtz Potentials"
	If the displacement field `` \mathbf u = \mathbf u(\mathbf x, t) `` satisfies the elastic wave equation
	```math
	\rho\ddot{\mathbf u} = \mathbf f + (\lambda+2\mu)\nabla(\nabla\cdot\mathbf u) - \mu\,\nabla\times(\nabla\times\mathbf u),
	```
	then there exist potentials `` \phi `` and `` \boldsymbol\psi `` such that
	```math
	\mathbf u = \nabla\phi + \nabla\times\boldsymbol\psi, \qquad \nabla\cdot\boldsymbol\psi = 0,
	```
	each obeying its own wave equation, sourced by the matching split of the body force
	`` \mathbf f = \nabla f_1 + \nabla\times\mathbf f_2 ``:
	```math
	\rho\ddot\phi = f_1 + (\lambda+2\mu)\nabla^2\phi, \qquad \rho\ddot{\boldsymbol\psi} = \mathbf f_2 + \mu\,\nabla^2\boldsymbol\psi.
	```
	In words: **any** elastic wave, however complicated, is exactly a sum of a curl-free piece
	(P waves, speed `` \alpha=\sqrt{(\lambda+2\mu)/\rho} ``) and a divergence-free piece (S waves,
	speed `` \beta=\sqrt{\mu/\rho} ``) — the split this whole notebook is built around.
"""

# ╔═╡ baf5aa58-0d0f-4e2c-803e-b1a8d98580a4
md"""
## P Radiation

Drag the receiver dot around the widget above while watching the **P** field. Its displacement
always points exactly along the line from the source to the receiver — push in the direction the
force points, pull behind it:

```math
\boxed{u_i^{\rm P}(\mathbf x,t) = \frac{\gamma_i\gamma_1}{4\pi\rho\alpha^2 r}\,X_0\!\left(t-\frac r\alpha\right)}
```

where `` \gamma_i=\partial r/\partial x_i `` is the unit vector from source to receiver and
`` X_0(t) `` is the force's time history. The `` \gamma_i\gamma_1 `` factor is exactly
`` \cos\theta ``, `` \theta `` being the angle between the receiver direction and the force —
plotted against angle it is the classic **dumbbell** lobe: maximal along the force axis, zero
perpendicular to it, and sign-flipped behind the source.

!!! note "Longitudinal waves"
	The far-field P displacement at any receiver has a direction *parallel* to the direction from
	the source — this is what "longitudinal" means. It is only exactly true far from the source;
	see the curl check in the Appendix for how "far" that has to be.
"""

# ╔═╡ b954f277-9840-465e-926e-5c4963751487
md"""
## S Radiation

Now watch the **S** field instead. Its displacement is always *perpendicular* to the
source-receiver line:

```math
\boxed{u_i^{\rm S}(\mathbf x,t) = \frac{\delta_{i1}-\gamma_i\gamma_1}{4\pi\rho\beta^2 r}\,X_0\!\left(t-\frac r\beta\right)}
```

`` \delta_{i1}-\gamma_i\gamma_1 `` has magnitude `` \sin\theta `` — zero along the force axis,
maximal at 90° from it, the **donut** lobe that exactly complements P's dumbbell. Between the two,
every direction around the source is covered: wherever P is weak, S is strong, and vice versa.

!!! note "Transverse waves"
	The far-field S displacement at any receiver has a direction *perpendicular* to the direction
	from the source — the same "far enough" caveat as the P term above applies.
"""

# ╔═╡ fc71ba84-35f3-4269-b8c9-3308782ae9df
md"""
## Near Field vs. Far Field

The P and S terms above both decay as `` 1/r `` — pure geometric spreading. But they are not the
whole story: there is a third term, the **near field**, that decays as `` 1/r^3 `` (through a time
integral of the source rather than a simple time-shift) and carries *both* P- and S-type motion
at once. It falls off much faster with distance, so it only matters close to the source — drag the
receiver toward the source in the widget above and watch the near-field contribution (toggle it
off/on with the checkbox) grow relative to the far-field terms.

The crossover is set by comparing the receiver distance `` r `` to a wavelength: the **ratio**
readout in the widget is `` \omega r/\alpha ``, using the pulse's characteristic frequency
`` \omega = 2\pi/T ``, where `` T `` is the **period** slider -- i.e. `` 2\pi\times `` the
number of wavelengths between source and receiver. Drag it to see the crossover shift: a short
period (high frequency) pushes the far field to dominate close in; a long period keeps the near
field significant out to much larger distances.

```math
\boxed{\text{ratio}\gg1 \Rightarrow \text{far field dominates}\qquad \text{ratio}\ll1 \Rightarrow \text{near field dominates}}
```

!!! warning "The crossover is further out than `` \text{ratio}\sim1 `` suggests"
	Don't take "``\gg1``" too literally here: because of a `` 1/r^2 ``, not `` 1/r^3 ``, effective
	decay of this specific near-field kernel at large `` r `` (the integration window widens with
	`` r ``, partly canceling the formal `` 1/r^3 `` prefactor), the near field is still
	several times *stronger* than far-field P out to `` \text{ratio}\approx10 ``, and only clearly
	subdominant past `` \text{ratio}\approx20 \text{--} 30 ``. The inequality is real, just slower
	to kick in than the formula's exponents alone would suggest — drag the receiver far out and
	watch the near-field checkbox's effect on the seismogram shrink to confirm.

This is exactly the same near/far tradeoff that shows up in the Appendix's numerical check of the
"P is curl-free / S is divergence-free" claims above — both are only exactly true in the
`` \text{ratio}\to\infty `` limit.
"""

# ╔═╡ 7f488eb9-9372-41b4-95ec-69abb753ef96
md"""
## Beyond a Single Force: A General Moment Tensor

A single point force is a useful teaching example, but it isn't what an earthquake actually looks
like. A fault slipping past itself exerts no net force on the surrounding rock — for every push
there's an equal and opposite pull, so the source can't be represented by a vector (a force) at
all. It takes a symmetric second-rank tensor, the **moment tensor** `` M ``.

Remarkably, Lamé's theorem still holds exactly unchanged: the far-field response of *any* point
source in a homogeneous, isotropic medium splits into a purely radial P part and a purely
tangential S part. For a moment-tensor source (Aki & Richards eq. 4.29), the radiation pattern
generalizes from `` \cos\theta ``/`` \sin\theta `` to

```math
R^P(\hat\gamma) = \hat\gamma^{\!\top} M \hat\gamma, \qquad
\vec{R}^S(\hat\gamma) = M\hat\gamma - \left(\hat\gamma^{\!\top} M \hat\gamma\right)\hat\gamma,
```

where `` \hat\gamma `` is the unit vector from source to receiver. `` R^P `` is a signed scalar —
push or pull, exactly like `` \cos\theta `` was for the point force. `` \vec{R}^S `` is a genuine
vector, lying in the plane perpendicular to `` \hat\gamma `` (still purely tangential, still
divergence-free).

Edit the moment tensor matrix directly below, or drag the strike/dip/rake sliders to generate a
pure double-couple — the mechanism of an ordinary shear-slip earthquake, and exactly the
four-lobed "beachball" pattern seismologists read a fault's orientation from.
"""

# ╔═╡ e62a5517-0584-42fd-9372-07d6d0dc8119
_mt_ntheta, _mt_nphi = 18, 28

# ╔═╡ 672e6bb6-90a8-4a18-85a8-d854e797aa00
_mt_tgrid_wave = range(-30.0, 30.0; length=800)

# ╔═╡ ad81924a-96a4-486c-a24e-34fe50ca9fd3
md"""
## Appendix
"""

# ╔═╡ b32d43d5-7b7e-4c90-8660-a0328ef6bb60
md"""
### Geometry & Angular Patterns
"""

# ╔═╡ 9d59cfce-5ac4-46d9-bf66-4defb8a03546
begin
    """
    	radial_pattern(::Val{:P}, θrel)
    	radial_pattern(::Val{:S}, θrel)
    	radial_pattern(::Val{:Near}, θrel)

    Radiation-pattern factor for the component of displacement *along* the source-receiver
    direction, as a function of the angle `θrel` between the receiver and the force direction.
    Far-field P is purely radial (`cos θrel`); far-field S has no radial component at all; the
    near field has both, with a `2cos θrel` radial part -- derived from
    `` \\partial_i\\partial_1(1/r) = (3\\gamma_i\\gamma_1-\\delta_{i1})/r^3 `` projected onto the
    radial direction (verified against finite differences of `1/r` in the check below).
    """
    radial_pattern(::Val{:P}, θrel) = cos(θrel)
    radial_pattern(::Val{:S}, θrel) = 0.0
    radial_pattern(::Val{:Near}, θrel) = 2 * cos(θrel)

    """
    	transverse_pattern(::Val{:P}, θrel)
    	transverse_pattern(::Val{:S}, θrel)
    	transverse_pattern(::Val{:Near}, θrel)

    Radiation-pattern factor for the component of displacement *perpendicular* to the
    source-receiver direction. Far-field S is purely transverse (`sin θrel`); far-field P has
    none; the near field has both (see [`radial_pattern`](@ref)).
    """
    transverse_pattern(::Val{:P}, θrel) = 0.0
    transverse_pattern(::Val{:S}, θrel) = sin(θrel)
    transverse_pattern(::Val{:Near}, θrel) = sin(θrel)
end

# ╔═╡ 0b0ea15d-957e-467a-beea-961284f67c63
md"""
### Source Time Functions
"""

# ╔═╡ cbc35a5e-0cc7-4ece-aac0-8520759d70fa
begin
    """
    	sourcewave(::Val{:gaussian}, t, σ²)

    Unit-strength source time function `` X_0(t) `` : a Gaussian pulse of width `σ²`.
    """
    sourcewave(::Val{:gaussian}, t, σ²) = exp(-t^2 / σ²)

    """
    	sourcewave_deriv(::Val{:gaussian}, t, σ²)

    `` \\dot{X}_0(t) ``, the time derivative of [`sourcewave`](@ref). Aki & Richards eq. (4.23)
    uses `` X_0 `` itself in the near-field integral but `` \\dot{X}_0 `` in the far-field P and S
    terms -- the classic tell is a step-force (`` X_0=H(t) ``): its far field is then an impulsive
    `` \\delta(t) `` pulse riding on top of a near field that builds smoothly to a permanent static
    offset, which is exactly the textbook picture.
    """
    sourcewave_deriv(::Val{:gaussian}, t, σ²) = -2t / σ² * exp(-t^2 / σ²)

    """
    	nearfield_kernel(::Val{:gaussian}, t, r, α, β, σ²)

    The near-field time integral `` K(t,r)=\\int_{r/\\alpha}^{r/\\beta}\\tau\\,X_0(t-\\tau)\\,d\\tau ``
    in closed form. The integration window `[r/α, r/β]` is exactly the time it takes the P and S
    fronts to sweep past a point at distance `r` -- everything in between is the near field.

    !!! note "A bug this fixes"
    	An earlier version of this notebook wrote the Gaussian case with the same `erf(...)` term
    	twice (`tp` where the second occurrence should have read `ts`), which silently evaluated
    	to exactly zero for every input. The corrected closed form below was verified against
    	direct numerical quadrature of the defining integral before shipping.
    """
    function nearfield_kernel(::Val{:gaussian}, t, r, α, β, σ²)
        tp, ts = r / α, r / β
        t * 0.5 * sqrt(σ² * π) * (erf((t - tp) / sqrt(σ²)) - erf((t - ts) / sqrt(σ²))) +
        0.5 * σ² * (sourcewave(Val(:gaussian), t - tp, σ²) - sourcewave(Val(:gaussian), t - ts, σ²))
    end
end

# ╔═╡ 9b7a2e1c-4c1a-4a6e-8f3a-6a8c1d2f5e01
md"""
!!! correct "Checking `sourcewave_deriv` against a finite difference"
	`sourcewave_deriv` is compared below to a centered finite difference of `sourcewave` itself,
	at several sample times -- the exact check that would have caught the far-field terms
	silently using `` X_0 `` instead of `` \dot{X}_0 ``.
"""

# ╔═╡ 9b7a2e1c-4c1a-4a6e-8f3a-6a8c1d2f5e02
let
    h = 1e-6
    for t in (-3.0, -0.7, 0.0, 0.4, 2.1)
        fd_gauss = (sourcewave(Val(:gaussian), t + h, 5.0) - sourcewave(Val(:gaussian), t - h, 5.0)) / 2h
        @assert isapprox(sourcewave_deriv(Val(:gaussian), t, 5.0), fd_gauss; atol=1e-6)
    end
    "sourcewave_deriv matches a finite-difference derivative of sourcewave at every sample point ✓"
end

# ╔═╡ 18b1d5bd-de79-41a4-9324-963d17a145f7
md"""
### Displacement Green's Functions
"""

# ╔═╡ 153caeb6-34d2-4fa6-a0f5-7c76f101a04a
begin
    const RHO = 5e12 # kg/km^3 -- a representative crustal density, in this notebook's km/(km/s) units
    const SIGMA2_DEFAULT = 5.0 # fixed Gaussian pulse width for the moment-tensor widget (no period control there)

    """
    	period_to_sigma2(T)

    Convert a dominant period `T` (seconds) to the Gaussian pulse-width parameter `σ²` used by
    [`sourcewave`](@ref), via the same characteristic frequency `` \\omega=1/\\sqrt{\\sigma^2} ``
    already used for the near/far-field ratio readout -- i.e. `` T = 2\\pi\\sqrt{\\sigma^2} ``.
    Letting the user set a period directly is far more intuitive than exposing `σ²` itself, and
    it's exactly the knob that controls how far out the near field stays significant relative to
    the far field (see the "Near Field vs. Far Field" section above).
    """
    period_to_sigma2(T) = (T / (2π))^2

    """
    	pfar_displacement(θrel, r, t, α, sourcetype, params...)

    Far-field P displacement (a signed scalar along the purely radial direction), Aki & Richards
    eq. (4.23)'s middle term. Uses [`sourcewave_deriv`](@ref) (`` \\dot{X}_0 ``), not `sourcewave`
    itself -- see that function's docstring for why.
    """
    pfar_displacement(θrel, r, t, α, sourcetype, params...) =
        radial_pattern(Val(:P), θrel) / (4π * RHO * α^2 * r) * sourcewave_deriv(sourcetype, t - r / α, params...)

    """
    	sfar_displacement(θrel, r, t, β, sourcetype, params...)

    Far-field S displacement (a signed scalar along the purely transverse direction), Aki &
    Richards eq. (4.23)'s last term. Uses [`sourcewave_deriv`](@ref), not `sourcewave` itself.
    """
    sfar_displacement(θrel, r, t, β, sourcetype, params...) =
        transverse_pattern(Val(:S), θrel) / (4π * RHO * β^2 * r) * sourcewave_deriv(sourcetype, t - r / β, params...)

    """
    	near_displacement(θrel, r, t, α, β, sourcetype, params...)

    Near-field displacement, returned as `(radial, transverse)` -- unlike the far-field terms it
    has both components at once, scaled by [`nearfield_kernel`](@ref) and decaying as `1/r³`.
    """
    function near_displacement(θrel, r, t, α, β, sourcetype, params...)
        K = nearfield_kernel(sourcetype, t, r, α, β, params...)
        scale = K / (4π * RHO * r^3)
        return (radial_pattern(Val(:Near), θrel) * scale, transverse_pattern(Val(:Near), θrel) * scale)
    end
end

# ╔═╡ 92539fa5-60db-4e54-b02c-79da24b2d9b9
md"""
### Verifying the Radiation Patterns

The claims above -- P is exactly curl-free, S is exactly divergence-free -- are only exactly
true in the far-field *limit*. The full closed-form expressions (evaluated at real, finite `r`)
also pick up a small correction from differentiating the `1/r` geometric-spreading factor itself,
which is the *same* near/far tradeoff the rest of this notebook is about. So a fair numerical
check should show `curl(P)/|P|` and `div(S)/|S|` as small but **nonzero**, shrinking as
`` \omega r/\alpha `` grows -- not exactly zero.
"""

# ╔═╡ ed9930ba-ba7b-4295-a4c7-beb820fd2d66
let
    α, β = 4.0, 2.0
    # a Gaussian source (never exactly zero) rather than a monochromatic one, so this check
    # can't accidentally land on a sin(...)=0 zero-crossing and divide by a near-zero norm
    st, σ² = Val(:gaussian), 5.0
    force = [1.0, 0.0, 0.0] # unit force along x -- a plain Cartesian frame just for this check

    pfar_vec(x, y, z, t) = begin
        r = norm([x, y, z])
        γ = [x, y, z] ./ r
        (γ .* dot(γ, force)) ./ (4π * RHO * α^2 * r) .* sourcewave(st, t - r / α, σ²)
    end
    sfar_vec(x, y, z, t) = begin
        r = norm([x, y, z])
        γ = [x, y, z] ./ r
        (force .- γ .* dot(γ, force)) ./ (4π * RHO * β^2 * r) .* sourcewave(st, t - r / β, σ²)
    end
    function curl_fd(F, x, y, z, t; h=1e-3)
        [(F(x, y + h, z, t)[3] - F(x, y - h, z, t)[3]) / 2h - (F(x, y, z + h, t)[2] - F(x, y, z - h, t)[2]) / 2h,
            (F(x, y, z + h, t)[1] - F(x, y, z - h, t)[1]) / 2h - (F(x + h, y, z, t)[3] - F(x - h, y, z, t)[3]) / 2h,
            (F(x + h, y, z, t)[2] - F(x - h, y, z, t)[2]) / 2h - (F(x, y + h, z, t)[1] - F(x, y - h, z, t)[1]) / 2h]
    end
    div_fd(F, x, y, z, t; h=1e-3) =
        (F(x + h, y, z, t)[1] - F(x - h, y, z, t)[1]) / 2h + (F(x, y + h, z, t)[2] - F(x, y - h, z, t)[2]) / 2h +
        (F(x, y, z + h, t)[3] - F(x, y, z - h, t)[3]) / 2h

    map([15.0, 40.0, 100.0]) do r
        # evaluate right at the P pulse's own peak (retarded time zero, never a zero-crossing)
        x, y, z, t = r / sqrt(3), r / sqrt(3), r / sqrt(3), r / α
        curlP_rel = norm(curl_fd(pfar_vec, x, y, z, t)) / max(norm(pfar_vec(x, y, z, t)), 1e-300)
        divS_rel = abs(div_fd(sfar_vec, x, y, z, t)) / max(norm(sfar_vec(x, y, z, t)), 1e-300)
        effω = 1 / sqrt(σ²)
        (; r, ωr_over_α=round(effω * r / α; digits=2), curlP_rel=round(curlP_rel; sigdigits=3), divS_rel=round(divS_rel; sigdigits=3))
    end
end

# ╔═╡ 54313d17-0de1-4d94-a9ed-948eadff5cf2
md"""
### Field Sampling
"""

# ╔═╡ 6787dae7-467b-4bf6-a643-5f318b109a32
"""
	pointforce_seismogram(fhat, rhat, r, α, β, sourcetype, params, tgrid)

Evaluate the radial/transverse displacement time series at a receiver in direction `rhat`
(a 3D unit vector from the source) at distance `r`, for a point force along direction `fhat`
(also a unit vector) -- decomposed into near-field, far-field-P, and far-field-S contributions
*separately* -- so the widget can toggle each on/off without a further Julia call. Only the
angle between `fhat` and `rhat` matters (the radiation is axisymmetric about the force axis),
computed here as `θrel = acos(fhat·rhat)`.

Returns a named tuple `(near_radial, near_transverse, pfar_radial, sfar_transverse)`, each a
vector the same length as `tgrid`.
"""
function pointforce_seismogram(fhat, rhat, r, α, β, sourcetype, params, tgrid)
    θrel = acos(clamp(dot(fhat, rhat), -1.0, 1.0))
    near_radial = Vector{Float64}(undef, length(tgrid))
    near_transverse = Vector{Float64}(undef, length(tgrid))
    pfar_radial = Vector{Float64}(undef, length(tgrid))
    sfar_transverse = Vector{Float64}(undef, length(tgrid))
    for (i, t) in enumerate(tgrid)
        nr, nt = near_displacement(θrel, r, t, α, β, sourcetype, params...)
        near_radial[i] = nr
        near_transverse[i] = nt
        pfar_radial[i] = pfar_displacement(θrel, r, t, α, sourcetype, params...)
        sfar_transverse[i] = sfar_displacement(θrel, r, t, β, sourcetype, params...)
    end
    return (; near_radial, near_transverse, pfar_radial, sfar_transverse)
end

# ╔═╡ 27257ff9-d835-476a-923b-fa8dc57cbd69
"""
	sourcewave_deriv_table(sourcetype, params, tgrid)

Sample `` \\dot{X}_0(t) `` (see [`sourcewave_deriv`](@ref)) on `tgrid` -- a small 1D lookup table
the widget interpolates client-side to animate the far-field wavefronts. It's the derivative, not
`sourcewave` itself, because the far-field P/S terms it's driving are `` \\dot{X}_0 `` (see
[`pfar_displacement`](@ref)). The angular pattern and `1/r` geometric spreading are simple enough
to compute directly in JS; the source waveform is the one piece of real physics worth pushing
over from Julia.
"""
sourcewave_deriv_table(sourcetype, params, tgrid) = [sourcewave_deriv(sourcetype, t, params...) for t in tgrid]

# ╔═╡ 970be60f-3c20-4f57-909d-c163b388a9dc
_mt_wave = sourcewave_deriv_table(Val(:gaussian), (SIGMA2_DEFAULT,), _mt_tgrid_wave)

# ╔═╡ 757a080d-c03f-4d8f-9bb8-e0eb8aa4eca0
md"""
### The Interactive Widget
"""

# ╔═╡ a3154298-98dd-4625-b8b0-f60d1070c31a
begin
    struct PointForceRadiationInput
        fx::Float64
        fy::Float64
        fz::Float64
        rdx::Float64
        rdy::Float64
        rdz::Float64
        rdist::Float64
        alpha::Float64
        beta::Float64
        period::Float64
    end
    function PointForceRadiationInput(; fx=0.0, fy=0.0, fz=1.0, rdx=0.6, rdy=0.3, rdz=0.74,
        rdist=60.0, alpha=4.0, beta=2.0, period=10.0)
        fn = hypot(fx, fy, fz)
        fx, fy, fz = fx / fn, fy / fn, fz / fn
        rn = hypot(rdx, rdy, rdz)
        rdx, rdy, rdz = rdx / rn, rdy / rn, rdz / rn
        PointForceRadiationInput(fx, fy, fz, rdx, rdy, rdz, rdist, alpha, beta, period)
    end

    Base.get(w::PointForceRadiationInput) = Dict{String,Any}(
        "fx" => w.fx, "fy" => w.fy, "fz" => w.fz,
        "rdx" => w.rdx, "rdy" => w.rdy, "rdz" => w.rdz, "rdist" => w.rdist,
        "alpha" => w.alpha, "beta" => w.beta, "period" => w.period)

    function Base.show(io::IO, ::MIME"text/html", w::PointForceRadiationInput)
        write(io, """
        <div id="pfrwidget">
        <style>
        pluto-cell:has(#pfrwidget) { width: min(80vw, 1500px) !important;
          margin-left: calc((100% - min(80vw, 1500px)) / 2) !important; }
        #pfrwidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #pfrwidget .pfr-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #pfrwidget .pfr-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #pfrwidget .pfr-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #pfrwidget .pfr-workspace{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start}
        #pfrwidget .pfr-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #pfrwidget .pfr-caption{font-size:13px;color:#9ca3af;text-align:center;margin-top:4px}
        #pfrwidget canvas{display:block}
        #pfrwidget #pfr-radiation{cursor:grab}
        #pfrwidget .pfr-controls{flex:0 0 260px;width:260px;display:flex;flex-direction:column;
          gap:8px;font:14px sans-serif}
        #pfrwidget .pfr-control-group{background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
        #pfrwidget .pfr-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #pfrwidget .pfr-control-row{display:grid;grid-template-columns:70px minmax(0,1fr) 52px;gap:6px;align-items:center;margin:5px 0}
        #pfrwidget .pfr-control-row label{font-size:13px;color:#9ca3af}
        #pfrwidget .pfr-control-row input[type=range]{width:100%;min-width:0}
        #pfrwidget .pfr-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #pfrwidget .pfr-check-row{display:flex;align-items:center;gap:6px;margin:5px 0;font-size:13px;color:#d1d5db}
        #pfrwidget .pfr-actions{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
        #pfrwidget .pfr-readout{font-size:13px;color:#d1d5db;line-height:1.6}
        #pfrwidget .pfr-readout b{color:#e5e7eb}
        #pfrwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:13px;cursor:pointer}
        #pfrwidget button.active{background:#2563eb;border-color:#93c5fd}
        #pfrwidget .pfr-seis-panel{width:100%;box-sizing:border-box;margin-top:14px}
        </style>

        <div class="pfr-title">
          <div class="pfr-title-desc">A point force radiates P and S waves at once — watch the wavefronts expand, painted by their own particle motion.</div>
          <div class="pfr-title-hint">press Play to expand the P/S wavefronts, their surfaces arrowed by push/pull polarity &middot; drag empty space to rotate &middot; drag the yellow arrow to set the force direction &middot; drag the white dot to set the receiver direction</div>
        </div>

        <div class="pfr-workspace">
          <div>
            <div class="pfr-panel"><canvas id="pfr-radiation"></canvas></div>
            <div class="pfr-caption" id="pfr-caption"></div>
          </div>
          <div class="pfr-controls">
            <div class="pfr-control-group">
              <div class="pfr-control-title">View</div>
              <div class="pfr-actions">
                <button id="pfr-view-anim" type="button">Wavefront animation</button>
                <button id="pfr-view-pattern" type="button">Radiation pattern</button>
              </div>
            </div>
            <div class="pfr-control-group">
              <div class="pfr-control-title">Medium</div>
              <div class="pfr-control-row"><label>α (km/s)</label><input type="range" id="pfr-alpha" min="2" max="8" step="0.1" value="$(w.alpha)"><span class="pfr-value" id="pfr-alpha-v"></span></div>
              <div class="pfr-control-row"><label>β (km/s)</label><input type="range" id="pfr-beta" min="1" max="5" step="0.1" value="$(w.beta)"><span class="pfr-value" id="pfr-beta-v"></span></div>
              <div class="pfr-control-row"><label>period T</label><input type="range" id="pfr-period" min="2" max="30" step="0.5" value="$(w.period)"><span class="pfr-value" id="pfr-period-v"></span></div>
            </div>
            <div class="pfr-control-group">
              <div class="pfr-control-title">Wavefront</div>
              <div class="pfr-actions"><button id="pfr-play" type="button">Play</button><button id="pfr-reset" type="button">Reset</button></div>
            </div>
            <div class="pfr-control-group">
              <div class="pfr-control-title">Seismogram terms</div>
              <div class="pfr-check-row"><input type="checkbox" id="pfr-chk-near" checked><label for="pfr-chk-near">Near field</label></div>
              <div class="pfr-check-row"><input type="checkbox" id="pfr-chk-pfar" checked><label for="pfr-chk-pfar">Far-field P</label></div>
              <div class="pfr-check-row"><input type="checkbox" id="pfr-chk-sfar" checked><label for="pfr-chk-sfar">Far-field S</label></div>
            </div>
            <div class="pfr-control-group">
              <div class="pfr-control-title">Readouts</div>
              <div class="pfr-readout" id="pfr-readout"></div>
            </div>
          </div>
        </div>

        <div class="pfr-seis-panel">
          <div class="pfr-panel"><canvas id="pfr-seismogram"></canvas></div>
          <div class="pfr-caption">displacement at the receiver &middot; click a legend label to show/hide that series</div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        const WORLD_HALF = 100; // km -- must match pointforce_seismogram/sourcewave_deriv_table's implicit domain in the Appendix

        let fx=$(w.fx), fy=$(w.fy), fz=$(w.fz);
        let rdx=$(w.rdx), rdy=$(w.rdy), rdz=$(w.rdz), rdist=$(w.rdist);
        let alpha=$(w.alpha), beta=$(w.beta), period=$(w.period);
        let effOmega = 0; // set from the first 'pfr-update' push (Julia computes it from period/σ², not JS)

        const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1500);
        const CONTROLS_W = 260+16;
        const SEC = Math.max(400, Math.min(availW - CONTROLS_W, 575));
        const DPR = window.devicePixelRatio || 1;

        const cv = par.querySelector('#pfr-radiation');
        const ctx = cv.getContext('2d');
        function hidpi(canvas, context, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          context.setTransform(DPR,0,0,DPR,0,0);
        }
        hidpi(cv, ctx, SEC, SEC);

        const seisCv = par.querySelector('#pfr-seismogram');
        const sctx = seisCv.getContext('2d');
        const SEIS_W = SEC + CONTROLS_W - 16, SEIS_H = 190;
        hidpi(seisCv, sctx, SEIS_W, SEIS_H);

        // ---- 3D camera: orthographic drag-to-rotate, same rot/unrot/proj/screenToXYZ family
        // as StressTensorInput in stress-tensor.jl -- p here is in NORMALIZED units, where a
        // vector of length 1 is the outer reference sphere (radius WORLD_HALF km).
        let yaw = -0.6, pitch = 0.35;
        const CX = SEC/2, CY = SEC/2, RPIX = SEC*0.36;
        function rot(p){
          const x = p[0]*Math.cos(yaw) - p[1]*Math.sin(yaw);
          const y = p[0]*Math.sin(yaw) + p[1]*Math.cos(yaw);
          const z = p[2];
          return [x, y*Math.cos(pitch) - z*Math.sin(pitch), y*Math.sin(pitch) + z*Math.cos(pitch)];
        }
        function unrot(P){
          const y1 = P[1]*Math.cos(pitch) + P[2]*Math.sin(pitch);
          const z1 = -P[1]*Math.sin(pitch) + P[2]*Math.cos(pitch);
          const px = P[0]*Math.cos(yaw) + y1*Math.sin(yaw);
          const py = -P[0]*Math.sin(yaw) + y1*Math.cos(yaw);
          return [px, py, z1];
        }
        function projUnit(p){ const r = rot(p); return [CX + r[0]*RPIX, CY - r[2]*RPIX]; }
        function projKm(pkm){ return projUnit([pkm[0]/WORLD_HALF, pkm[1]/WORLD_HALF, pkm[2]/WORLD_HALF]); }
        // inverse of projUnit, restricted to the front hemisphere of the unit sphere -- used to
        // set a DIRECTION (force or receiver) by dragging, independent of display radius
        function screenToXYZ(mx, my){
          let px = (mx-CX)/RPIX, pz = (CY-my)/RPIX;
          let r2 = px*px + pz*pz;
          if(r2 > 1){ const s = 1/Math.sqrt(r2); px *= s; pz *= s; r2 = 1; }
          const py = Math.sqrt(Math.max(0, 1-r2));
          const v = unrot([px, py, pz]);
          const n = Math.hypot(v[0], v[1], v[2]) || 1;
          return [v[0]/n, v[1]/n, v[2]/n];
        }

        // ---- data pushed from Julia: source waveform lookup table + seismogram arrays ----
        let waveData = null; // {tgridWave, wave}
        let seisData = null; // {tgridSeis, nearRadial, nearTransverse, pfarRadial, sfarTransverse}

        function lookupWave(t){
          if(!waveData) return 0;
          const {tgridWave, wave} = waveData;
          const n = tgridWave.length;
          const t0 = tgridWave[0], t1 = tgridWave[n-1];
          if(t <= t0) return wave[0];
          if(t >= t1) return wave[n-1];
          const frac = (t-t0)/(t1-t0)*(n-1);
          const i0 = Math.floor(frac), i1 = Math.min(n-1, i0+1);
          const f = frac - i0;
          return wave[i0]*(1-f) + wave[i1]*f;
        }

        // the wavefront shells sample the source waveform at retarded time 0 (the front's own
        // arrival instant) -- but this is now Ẋ (see sourcewave_deriv's docstring), and Ẋ is an
        // ODD function of t for a symmetric pulse like a Gaussian, so it's exactly zero AT t=0.
        // Sampling literally at the front would make every arrow invisible. Instead take the
        // largest-magnitude sample anywhere in the already-pushed waveform array -- Julia builds
        // that array to span the physically relevant window at fixed resolution regardless of
        // the pulse's width, so scanning all of it (not a fixed-width/fixed-sample-count local
        // window) is correct for any period the user picks, not just the one it was tuned for.
        function peakWaveNear(){
          if(!waveData) return 0;
          let best = 0, bestAbs = -1;
          for(let i=0; i<waveData.wave.length; i++){
            const v = waveData.wave[i];
            if(Math.abs(v) > bestAbs){ bestAbs = Math.abs(v); best = v; }
          }
          return best;
        }

        let playing = false, rafId = null, tPhase = 0;
        let viewMode = 'anim'; // 'anim' | 'pattern'
        function tLoopMax(){ return 1.3*WORLD_HALF/Math.min(alpha,beta); }

        // ---- Fibonacci sphere: an evenly-spaced set of sample directions, reused both for the
        // static radiation-pattern glyph and for the animated per-point vector arrows ----
        function fibonacciSphere(n){
          const pts = [];
          const golden = Math.PI*(3-Math.sqrt(5));
          for(let i=0;i<n;i++){
            const yv = 1 - (n===1?0:(i/(n-1))*2);
            const rr = Math.sqrt(Math.max(0,1-yv*yv));
            const theta = golden*i;
            pts.push([Math.cos(theta)*rr, yv, Math.sin(theta)*rr]);
          }
          return pts;
        }
        const FIB_PTS = fibonacciSphere(84);

        // decompose the force direction fhat relative to a receiver direction gamma into the
        // radial (cosT = cosThetaRel) and transverse (that, a genuine 3D unit vector) parts --
        // this is the JS-side geometry counterpart of near_displacement/pfar/sfar_displacement's
        // (radial,transverse) split in the Appendix, just evaluated per sample point here
        function basisAt(gamma, fhat){
          let c = gamma[0]*fhat[0] + gamma[1]*fhat[1] + gamma[2]*fhat[2];
          c = Math.max(-1, Math.min(1, c));
          let tx = fhat[0]-c*gamma[0], ty = fhat[1]-c*gamma[1], tz = fhat[2]-c*gamma[2];
          let tn = Math.hypot(tx,ty,tz);
          if(tn < 1e-6){
            const arb = Math.abs(gamma[1])<0.9 ? [0,1,0] : [1,0,0];
            const d = arb[0]*gamma[0]+arb[1]*gamma[1]+arb[2]*gamma[2];
            tx = arb[0]-d*gamma[0]; ty = arb[1]-d*gamma[1]; tz = arb[2]-d*gamma[2];
            tn = Math.hypot(tx,ty,tz) || 1;
          }
          return {cosT: c, that: [tx/tn, ty/tn, tz/tn]};
        }

        // Marker conventions shared across this repo's widgets (see the pluto-widget-style
        // skill): a source is drawn as a star, a receiver as a downward-pointing triangle.
        function drawStarMarker(cx, cy, r, fill, stroke){
          const spikes = 5, rOuter = r, rInner = r * 0.45;
          ctx.beginPath();
          for(let i=0; i<spikes*2; i++){
            const rad = i % 2 === 0 ? rOuter : rInner;
            const ang = -Math.PI/2 + i*Math.PI/spikes;
            const x = cx + rad*Math.cos(ang), y = cy + rad*Math.sin(ang);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.closePath();
          ctx.fillStyle = fill; ctx.fill();
          ctx.strokeStyle = stroke; ctx.lineWidth = 1; ctx.stroke();
        }
        function drawTriangleDownMarker(cx, cy, r, fill, stroke){
          ctx.beginPath();
          for(let i=0; i<3; i++){
            const ang = Math.PI/2 + i*2*Math.PI/3; // first vertex points down (canvas y grows downward)
            const x = cx + r*Math.cos(ang), y = cy + r*Math.sin(ang);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.closePath();
          ctx.fillStyle = fill; ctx.fill();
          ctx.strokeStyle = stroke; ctx.lineWidth = 1.5; ctx.stroke();
        }

        // arrowhead-only helper for the (already screen-projected) force handle line
        function drawArrowHead2D(x0, y0, x1, y1, color, size){
          size = size || 10;
          const ang = Math.atan2(y1-y0, x1-x0);
          ctx.beginPath(); ctx.moveTo(x1,y1);
          ctx.lineTo(x1-size*Math.cos(ang-0.4), y1-size*Math.sin(ang-0.4));
          ctx.lineTo(x1-size*Math.cos(ang+0.4), y1-size*Math.sin(ang+0.4));
          ctx.closePath(); ctx.fillStyle = color; ctx.fill();
        }

        // ---- Radiation-pattern view: a smooth SHADED SURFACE (not a point cloud) for the P
        // dumbbell (radius |cosThetaRel|, red=push/theta<90, blue=pull/theta>90) and the S donut
        // (radius sinThetaRel, cyan, magnitude-shaded), built on a (theta,phi) grid whose pole is
        // aligned with the live force direction -- gives each lobe real connectivity so it can be
        // filled as a mesh of quads instead of a scatter of dots. Both lobes are collected into
        // one depth-sorted painter's-algorithm pass so they occlude each other correctly no matter
        // how the camera is rotated.
        function forceFrame(fhat){
          const arb = Math.abs(fhat[1]) < 0.9 ? [0,1,0] : [1,0,0];
          const d = arb[0]*fhat[0]+arb[1]*fhat[1]+arb[2]*fhat[2];
          let e1 = [arb[0]-d*fhat[0], arb[1]-d*fhat[1], arb[2]-d*fhat[2]];
          const n = Math.hypot(e1[0],e1[1],e1[2]) || 1;
          e1 = [e1[0]/n, e1[1]/n, e1[2]/n];
          const e2 = [fhat[1]*e1[2]-fhat[2]*e1[1], fhat[2]*e1[0]-fhat[0]*e1[2], fhat[0]*e1[1]-fhat[1]*e1[0]];
          return {e1, e2};
        }

        const LOBE_NTHETA = 26, LOBE_NPHI = 40;
        function lobeVertex(theta, phi, fhat, e1, e2, radiusFn){
          const st = Math.sin(theta), ct = Math.cos(theta);
          const g = [
            e1[0]*st*Math.cos(phi) + e2[0]*st*Math.sin(phi) + fhat[0]*ct,
            e1[1]*st*Math.cos(phi) + e2[1]*st*Math.sin(phi) + fhat[1]*ct,
            e1[2]*st*Math.cos(phi) + e2[2]*st*Math.sin(phi) + fhat[2]*ct,
          ];
          const r = radiusFn(ct);
          return [g[0]*r, g[1]*r, g[2]*r];
        }

        function collectLobeQuads(fhat, e1, e2, radiusFn, colorFn, quads){
          const grid = [];
          for(let i=0;i<=LOBE_NTHETA;i++){
            const theta = Math.PI*i/LOBE_NTHETA;
            const row = [];
            for(let j=0;j<=LOBE_NPHI;j++){
              const phi = 2*Math.PI*j/LOBE_NPHI;
              row.push(lobeVertex(theta, phi, fhat, e1, e2, radiusFn));
            }
            grid.push(row);
          }
          for(let i=0;i<LOBE_NTHETA;i++){
            for(let j=0;j<LOBE_NPHI;j++){
              const v00=grid[i][j], v01=grid[i][j+1], v11=grid[i+1][j+1], v10=grid[i+1][j];
              const ctr = [(v00[0]+v01[0]+v11[0]+v10[0])/4, (v00[1]+v01[1]+v11[1]+v10[1])/4, (v00[2]+v01[2]+v11[2]+v10[2])/4];
              quads.push({ verts: [v00,v01,v11,v10], depth: rot(ctr)[1], color: colorFn(ctr) });
            }
          }
        }

        function drawRadiationSurfaces(){
          const fhat = [fx,fy,fz];
          const {e1, e2} = forceFrame(fhat);
          const quads = [];
          collectLobeQuads(fhat, e1, e2,
            ct => 0.08 + Math.abs(ct)*0.62,
            ctr => {
              const n = Math.hypot(ctr[0],ctr[1],ctr[2]) || 1;
              const cosT = (ctr[0]*fhat[0]+ctr[1]*fhat[1]+ctr[2]*fhat[2])/n;
              const mag = Math.min(1, n/0.7);
              return cosT >= 0
                ? 'rgba(239,68,68,' + (0.25+0.65*mag) + ')'
                : 'rgba(59,130,246,' + (0.25+0.65*mag) + ')';
            }, quads);
          collectLobeQuads(fhat, e1, e2,
            ct => 0.08 + Math.sqrt(Math.max(0,1-ct*ct))*0.52,
            ctr => {
              const n = Math.hypot(ctr[0],ctr[1],ctr[2]) || 1;
              const mag = Math.min(1, n/0.6);
              return 'rgba(34,211,238,' + (0.20+0.55*mag) + ')';
            }, quads);
          quads.sort((a,b)=>a.depth-b.depth);
          for(const q of quads){
            ctx.beginPath();
            const [x0,y0] = projUnit(q.verts[0]); ctx.moveTo(x0,y0);
            for(let k=1;k<4;k++){ const [x,y] = projUnit(q.verts[k]); ctx.lineTo(x,y); }
            ctx.closePath();
            ctx.fillStyle = q.color; ctx.fill();
          }
        }

        // animated expanding P/S wavefronts: each front is drawn as an actual spherical SURFACE
        // (radius alpha*tPhase for P, beta*tPhase for S -- true geometric wavefronts, not a
        // stand-in silhouette) carrying a DENSE field of displacement-direction arrows sampled
        // on that surface -- radial for P (longitudinal), tangential for S (transverse) --
        // colored red/blue by the sign of the particle motion (push/outward vs. pull/inward for
        // P; the two shear senses for S). Because each shell's radius is defined to be exactly
        // v*tPhase, every point on it is at its own retarded time zero -- i.e. the front has
        // JUST arrived there -- so the arrows always show the wave's leading-edge motion, with
        // no artificial fade-in needed (unlike a fixed-position particle, a point exactly on the
        // expanding front is never "not yet reached").
        const SHELL_ARROW_DIRS = fibonacciSphere(140);

        // a FIXED (not per-frame) amplitude reference, set once whenever fresh Julia data
        // arrives -- normalizing arrow length/color against a value that itself changed every
        // frame is what made an earlier version of this animation flicker.
        let refAmp = 1e-12;
        const REF_R = 0.5*WORLD_HALF;

        // one expanding shell: wireframe surface (front hemisphere only, so the far side of the
        // sphere doesn't clutter the view) plus its dense arrow field.
        function drawWavefrontShell(R, fhat, wireColor, radial){
          if(!(R > 1 && R < WORLD_HALF*1.05)) return;
          const Rn = Math.min(1, R/WORLD_HALF);

          ctx.strokeStyle = wireColor; ctx.lineWidth = 1;
          const NLAT=4, NLON=6, NSEG=32;
          for(let i=1;i<NLAT;i++){
            const theta = Math.PI*i/NLAT;
            ctx.beginPath(); let started=false;
            for(let j=0;j<=NSEG;j++){
              const phi = 2*Math.PI*j/NSEG;
              const p = [Math.sin(theta)*Math.cos(phi)*Rn, Math.sin(theta)*Math.sin(phi)*Rn, Math.cos(theta)*Rn];
              if(rot(p)[1] < 0){ started=false; continue; }
              const s = projUnit(p);
              if(!started){ ctx.moveTo(s[0],s[1]); started=true; } else ctx.lineTo(s[0],s[1]);
            }
            ctx.stroke();
          }
          for(let j=0;j<NLON;j++){
            const phi = Math.PI*j/NLON;
            ctx.beginPath(); let started=false;
            for(let i=0;i<=NSEG;i++){
              const theta = Math.PI*i/NSEG;
              const p = [Math.sin(theta)*Math.cos(phi)*Rn, Math.sin(theta)*Math.sin(phi)*Rn, Math.cos(theta)*Rn];
              if(rot(p)[1] < 0){ started=false; continue; }
              const s = projUnit(p);
              if(!started){ ctx.moveTo(s[0],s[1]); started=true; } else ctx.lineTo(s[0],s[1]);
            }
            ctx.stroke();
          }

          const srcAtFront = peakWaveNear();
          const ARROW_LEN_MAX = WORLD_HALF*0.11;
          for(const g of SHELL_ARROW_DIRS){
            if(rot(g)[1] < 0.02) continue; // front hemisphere only -- avoids an occluded ball of arrows
            const {cosT, that} = basisAt(g, fhat);
            let amp, dir;
            if(radial){
              amp = cosT*srcAtFront/R;
              dir = g;
            } else {
              const sinT = Math.sqrt(Math.max(0, 1-cosT*cosT));
              amp = sinT*srcAtFront/R;
              dir = that;
            }
            const mag = Math.min(1, Math.abs(amp)/refAmp);
            if(mag < 0.05) continue;
            const len = mag*ARROW_LEN_MAX*Math.sign(amp);
            const base = [g[0]*R, g[1]*R, g[2]*R];
            const tip = [base[0]+dir[0]*len, base[1]+dir[1]*len, base[2]+dir[2]*len];
            const [x0,y0] = projKm(base), [x1,y1] = projKm(tip);
            const color = amp >= 0
              ? 'rgba(239,68,68,' + (0.4+0.55*mag) + ')'
              : 'rgba(59,130,246,' + (0.4+0.55*mag) + ')';
            ctx.strokeStyle = color; ctx.lineWidth = 1.6;
            ctx.beginPath(); ctx.moveTo(x0,y0); ctx.lineTo(x1,y1); ctx.stroke();
            drawArrowHead2D(x0,y0,x1,y1,color,5);
          }
        }

        function drawWavefronts(){
          if(!waveData) return;
          const fhat = [fx,fy,fz];
          drawWavefrontShell(alpha*tPhase, fhat, 'rgba(249,115,22,0.22)', true);
          drawWavefrontShell(beta*tPhase, fhat, 'rgba(34,211,238,0.22)', false);
        }

        // fixed x1/x2/x3 reference frame the force/receiver direction cosines are actually
        // defined against -- without this there's no visible cue for which screen direction is
        // "1" vs "2" vs "3" (same convention as StressTensorInput in stress-tensor.jl).
        function drawAxes(){
          const AXLEN = 1.18;
          const axes = [ [[1,0,0],'x₁'], [[0,1,0],'x₂'], [[0,0,1],'x₃'] ];
          ctx.lineWidth = 1;
          ctx.font = '13px sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
          for(const [v,label] of axes){
            ctx.strokeStyle = 'rgba(156,163,175,0.55)';
            const [x0,y0] = projUnit([-v[0]*AXLEN,-v[1]*AXLEN,-v[2]*AXLEN]);
            const [x1,y1] = projUnit([v[0]*AXLEN,v[1]*AXLEN,v[2]*AXLEN]);
            ctx.beginPath(); ctx.moveTo(x0,y0); ctx.lineTo(x1,y1); ctx.stroke();
            ctx.fillStyle = '#9ca3af'; ctx.fillText(label, x1, y1);
          }
          ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
        }

        function draw(){
          ctx.clearRect(0,0,SEC,SEC);
          drawAxes();

          // near-field zone indicator: a translucent disk (screen-facing billboard) sized by a
          // characteristic wavelength, centered on the source
          // effOmega pushed from Julia (computed from the period slider) -- see 'pfr-update' below
          const nearR = Math.min(alpha,beta)/effOmega;
          const [ncx,ncy] = projUnit([0,0,0]);
          ctx.beginPath(); ctx.arc(ncx, ncy, Math.max(4, (nearR/WORLD_HALF)*RPIX), 0, 2*Math.PI);
          ctx.fillStyle = 'rgba(250,204,21,0.10)'; ctx.fill();
          ctx.strokeStyle = 'rgba(250,204,21,0.35)'; ctx.lineWidth = 1; ctx.stroke();

          if(viewMode === 'pattern'){ drawRadiationSurfaces(); } else { drawWavefronts(); }

          // force arrow -- fixed display length, direction only, drag anywhere on the sphere
          const FORCE_DISP_R = 0.30;
          const [sx0,sy0] = projUnit([0,0,0]);
          const [sx1,sy1] = projUnit([fx*FORCE_DISP_R, fy*FORCE_DISP_R, fz*FORCE_DISP_R]);
          ctx.strokeStyle = '#111827'; ctx.lineWidth = 6;
          ctx.beginPath(); ctx.moveTo(sx0,sy0); ctx.lineTo(sx1,sy1); ctx.stroke();
          ctx.strokeStyle = '#facc15'; ctx.lineWidth = 3;
          ctx.beginPath(); ctx.moveTo(sx0,sy0); ctx.lineTo(sx1,sy1); ctx.stroke();
          drawArrowHead2D(sx0,sy0,sx1,sy1, '#facc15');

          // source marker
          drawStarMarker(sx0, sy0, 7, '#e5e7eb', '#000');
          ctx.fillStyle = '#9ca3af'; ctx.font = '12px sans-serif'; ctx.textAlign='left';
          ctx.fillText('source', sx0+10, sy0-8);

          // receiver -- direction (and, since a recent change, distance too) set by dragging
          const RECV_DISP_R = rdist / WORLD_HALF;
          const [rpx,rpy] = projUnit([rdx*RECV_DISP_R, rdy*RECV_DISP_R, rdz*RECV_DISP_R]);
          drawTriangleDownMarker(rpx, rpy, 7, '#f5f3ef', '#0a0f18');
          ctx.fillStyle = '#e5e7eb'; ctx.font = '12px sans-serif';
          ctx.fillText('receiver', rpx+9, rpy-8);

          par.querySelector('#pfr-caption').textContent = 'drag to rotate  ·  r = ' + rdist.toFixed(0) + ' km';

          updateReadout();
          drawSeismogram();
        }

        function updateReadout(){
          const r = rdist;
          const tp = r/alpha, ts = r/beta;
          // effOmega pushed from Julia (computed from the period slider) -- see 'pfr-update' below
          const ratio = effOmega*r/alpha;
          par.querySelector('#pfr-readout').innerHTML =
            'distance <b>' + r.toFixed(1) + '</b> km<br>' +
            'P arrival <b>' + tp.toFixed(2) + '</b> s<br>' +
            'S arrival <b>' + ts.toFixed(2) + '</b> s<br>' +
            'S&minus;P <b>' + (ts-tp).toFixed(2) + '</b> s<br>' +
            'ratio ωr/α <b>' + ratio.toFixed(2) + '</b> (' + (ratio>3?'far field':(ratio<0.5?'near field':'transition')) + ')<br>' +
            '<span style="color:#ef4444">&#9679;</span> push/forward &nbsp; <span style="color:#3b82f6">&#9679;</span> pull/back';
        }

        function niceSeisY(v, vmax, h, pad){
          return h/2 - (vmax>0 ? v/vmax : 0)*(h/2-pad);
        }

        // which seismogram series are currently plotted -- clicking a legend label toggles its
        // entry, same interaction as the comparison-panel legend in seismic-interferometry.jl
        let seisShow = { radial: true, transverse: true, u1: false, u2: false, u3: false };
        let seisLegendHits = [];
        const SEIS_SERIES = [
          { key: 'radial', label: 'radial', color: '#f97316', dash: null },
          { key: 'transverse', label: 'transverse', color: '#38bdf8', dash: [6,3] },
          { key: 'u1', label: 'u₁', color: '#facc15', dash: null },
          { key: 'u2', label: 'u₂', color: '#4ade80', dash: null },
          { key: 'u3', label: 'u₃', color: '#c084fc', dash: null },
        ];
        function seisLegendHit(px, py){
          return seisLegendHits.find(h => px>=h.x0 && px<=h.x1 && py>=h.y0 && py<=h.y1);
        }

        function drawSeismogram(){
          sctx.clearRect(0,0,SEIS_W,SEIS_H);
          sctx.strokeStyle = '#374151'; sctx.lineWidth = 1; sctx.strokeRect(0.5,0.5,SEIS_W-1,SEIS_H-1);
          seisLegendHits = [];
          if(!seisData){
            sctx.fillStyle = '#6b7280'; sctx.font = '12px sans-serif'; sctx.fillText('computing...', 10, 18);
            return;
          }
          const {tgridSeis, nearRadial, nearTransverse, pfarRadial, sfarTransverse} = seisData;
          const n = tgridSeis.length;
          const showNear = par.querySelector('#pfr-chk-near').checked;
          const showP = par.querySelector('#pfr-chk-pfar').checked;
          const showS = par.querySelector('#pfr-chk-sfar').checked;

          // radial/transverse are scalars along rhat/that; u1,u2,u3 are the same physical
          // displacement re-expressed in the fixed x1/x2/x3 lab frame -- pure recombination of
          // already-computed physics, no new Green's-function evaluation needed here.
          const {that} = basisAt([rdx,rdy,rdz], [fx,fy,fz]);
          const radial = new Array(n), transverse = new Array(n);
          const u1 = new Array(n), u2 = new Array(n), u3 = new Array(n);
          for(let i=0;i<n;i++){
            radial[i] = (showNear?nearRadial[i]:0) + (showP?pfarRadial[i]:0);
            transverse[i] = (showNear?nearTransverse[i]:0) + (showS?sfarTransverse[i]:0);
            u1[i] = radial[i]*rdx + transverse[i]*that[0];
            u2[i] = radial[i]*rdy + transverse[i]*that[1];
            u3[i] = radial[i]*rdz + transverse[i]*that[2];
          }
          const seriesValues = { radial, transverse, u1, u2, u3 };

          let vmax = 1e-300;
          for(const s of SEIS_SERIES){
            if(!seisShow[s.key]) continue;
            const arr = seriesValues[s.key];
            for(let i=0;i<n;i++) vmax = Math.max(vmax, Math.abs(arr[i]));
          }

          const t0 = tgridSeis[0], t1 = tgridSeis[n-1];
          const PAD = 30;
          function xOf(t){ return PAD + (t-t0)/(t1-t0)*(SEIS_W-PAD-10); }

          sctx.strokeStyle = '#2f3744'; sctx.beginPath();
          sctx.moveTo(PAD, SEIS_H/2); sctx.lineTo(SEIS_W-10, SEIS_H/2); sctx.stroke();

          for(const s of SEIS_SERIES){
            if(!seisShow[s.key]) continue;
            const arr = seriesValues[s.key];
            sctx.strokeStyle = s.color; sctx.lineWidth = 1.8; sctx.setLineDash(s.dash||[]);
            sctx.beginPath();
            for(let i=0;i<n;i++){
              const x = xOf(tgridSeis[i]), y = niceSeisY(arr[i], vmax, SEIS_H, 14);
              i===0 ? sctx.moveTo(x,y) : sctx.lineTo(x,y);
            }
            sctx.stroke(); sctx.setLineDash([]);
          }

          // time cursor synced to the wavefront animation
          const cxp = xOf(tPhase);
          if(cxp>=PAD && cxp<=SEIS_W-10){
            sctx.strokeStyle = '#e5e7eb'; sctx.lineWidth = 1;
            sctx.beginPath(); sctx.moveTo(cxp, 6); sctx.lineTo(cxp, SEIS_H-6); sctx.stroke();
          }

          sctx.fillStyle = '#6b7280'; sctx.font = '11px sans-serif'; sctx.textAlign='left';
          sctx.fillText('t=' + t0.toFixed(0) + 's', 2, SEIS_H-2);
          sctx.textAlign='right'; sctx.fillText('t=' + t1.toFixed(0) + 's', SEIS_W-2, SEIS_H-2);

          // clickable legend -- toggles seisShow[key] and redraws, exactly like the comparison
          // panel's legendItem/comparisonLegendHit pattern in seismic-interferometry.jl
          sctx.font = '12px sans-serif'; sctx.textAlign = 'left';
          let lx = PAD;
          for(const s of SEIS_SERIES){
            const active = seisShow[s.key];
            sctx.fillStyle = active ? s.color : 'rgba(156,163,175,0.45)';
            sctx.fillText(s.label, lx, 12);
            const w = sctx.measureText(s.label).width;
            seisLegendHits.push({ key: s.key, x0: lx-3, y0: 2, x1: lx+w+3, y1: 16 });
            lx += w + 14;
          }
        }

        function syncControls(){
          par.querySelector('#pfr-view-anim').classList.toggle('active', viewMode==='anim');
          par.querySelector('#pfr-view-pattern').classList.toggle('active', viewMode==='pattern');
          par.querySelector('#pfr-play').style.display = viewMode==='anim' ? '' : 'none';
          par.querySelector('#pfr-alpha').value = alpha; par.querySelector('#pfr-alpha-v').textContent = alpha.toFixed(1);
          par.querySelector('#pfr-beta').value = beta; par.querySelector('#pfr-beta-v').textContent = beta.toFixed(1);
          par.querySelector('#pfr-period').value = period; par.querySelector('#pfr-period-v').textContent = period.toFixed(1)+' s';
        }

        let commitInFlight = false;
        function commit(){
          commitInFlight = true;
          par.value = { fx, fy, fz, rdx, rdy, rdz, rdist, alpha, beta, period };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledCommit(){ if(!commitInFlight) commit(); }

        par.addEventListener('pfr-update', e=>{
          seisData = { tgridSeis: e.detail.tgridSeis, nearRadial: e.detail.nearRadial,
            nearTransverse: e.detail.nearTransverse, pfarRadial: e.detail.pfarRadial, sfarTransverse: e.detail.sfarTransverse };
          waveData = { tgridWave: e.detail.tgridWave, wave: e.detail.wave };
          effOmega = e.detail.effOmega;
          const maxAbsWave = Math.max(1e-12, ...waveData.wave.map(Math.abs));
          refAmp = maxAbsWave / REF_R;
          commitInFlight = false;
          draw();
        });

        // ---- dragging: force direction, receiver direction, or (empty space) rotate the
        // camera -- closest-handle-wins hit test, since the two handles can land close
        // together on screen depending on orientation ----
        let dragMode = null, lastX = 0, lastY = 0;
        cv.addEventListener('mousedown', e=>{
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left, py = e.clientY-rect.top;
          const FORCE_DISP_R = 0.30;
          const RECV_DISP_R = rdist / WORLD_HALF;
          const [fpx,fpy] = projUnit([fx*FORCE_DISP_R, fy*FORCE_DISP_R, fz*FORCE_DISP_R]);
          const [rpx,rpy] = projUnit([rdx*RECV_DISP_R, rdy*RECV_DISP_R, rdz*RECV_DISP_R]);
          const dArrow = Math.hypot(px-fpx, py-fpy), dRecv = Math.hypot(px-rpx, py-rpy);
          let best = null, bestD = 14;
          if(dArrow < bestD){ bestD = dArrow; best = 'force'; }
          if(dRecv < bestD){ bestD = dRecv; best = 'receiver'; }
          dragMode = best || 'rotate';
          lastX = px; lastY = py;
        });
        window.addEventListener('mousemove', e=>{
          if(!dragMode) return;
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left, py = e.clientY-rect.top;
          if(dragMode==='rotate'){
            const dx = px-lastX, dy = py-lastY;
            yaw += dx*0.008;
            pitch = Math.max(-1.4, Math.min(1.4, pitch - dy*0.008));
            lastX = px; lastY = py;
            draw();
            return;
          }
          const v = screenToXYZ(px, py);
          if(dragMode==='force'){ fx=v[0]; fy=v[1]; fz=v[2]; }
          else if(dragMode==='receiver'){
            rdx=v[0]; rdy=v[1]; rdz=v[2];
            // radial screen-space distance from the source sets rdist too, so one drag
            // controls both direction and distance -- no separate slider needed
            const [sx0,sy0] = projUnit([0,0,0]);
            const pixelDist = Math.hypot(px-sx0, py-sy0);
            rdist = Math.max(5, Math.min(95, (pixelDist/RPIX)*WORLD_HALF));
          }
          draw(); throttledCommit();
        });
        window.addEventListener('mouseup', ()=>{
          if(dragMode==='force' || dragMode==='receiver') commit();
          dragMode = null;
        });

        // ---- controls ----
        par.addEventListener('input', e=>{
          if(e.target===par) return;
          e.stopImmediatePropagation();
          const id = e.target.id, v = e.target.value;
          if(id==='pfr-alpha'){ alpha=+v; par.querySelector('#pfr-alpha-v').textContent=alpha.toFixed(1); }
          else if(id==='pfr-beta'){ beta=+v; par.querySelector('#pfr-beta-v').textContent=beta.toFixed(1); }
          else if(id==='pfr-period'){ period=+v; par.querySelector('#pfr-period-v').textContent=period.toFixed(1)+' s'; }
          else return;
          draw(); throttledCommit();
        }, true);

        par.addEventListener('change', e=>{
          if(e.target===par) return;
          e.stopImmediatePropagation();
          const id = e.target.id;
          if(id==='pfr-chk-near'||id==='pfr-chk-pfar'||id==='pfr-chk-sfar'){ drawSeismogram(); return; }
          if(id==='pfr-alpha'||id==='pfr-beta'||id==='pfr-period'){ draw(); commit(); return; }
        }, true);

        // simulated time-units advanced per REAL second of wall-clock playback -- tied to the
        // animation's own rAF timestamp rather than a fixed per-frame increment, so playback
        // speed no longer depends on (and stutters with) the actual frame rate
        const SIM_SPEED = tLoopMax() / 13;
        const playBtn = par.querySelector('#pfr-play');
        let lastTs = null;
        function stepAnim(ts){
          if(lastTs === null) lastTs = ts;
          const dt = Math.min(0.1, (ts-lastTs)/1000); // clamp so a stalled/backgrounded tab can't leap
          lastTs = ts;
          tPhase += dt*SIM_SPEED;
          if(tPhase > tLoopMax()) tPhase = 0;
          draw();
          rafId = requestAnimationFrame(stepAnim);
        }
        playBtn.addEventListener('click', ()=>{
          playing = !playing;
          playBtn.textContent = playing ? 'Pause' : 'Play';
          if(playing){ lastTs = null; rafId = requestAnimationFrame(stepAnim); }
          else if(rafId){ cancelAnimationFrame(rafId); rafId = null; }
        });

        function stopAnim(){
          playing = false; playBtn.textContent = 'Play';
          if(rafId){ cancelAnimationFrame(rafId); rafId = null; }
        }
        par.querySelector('#pfr-view-anim').addEventListener('click', ()=>{
          if(viewMode==='anim') return;
          viewMode = 'anim'; syncControls(); draw();
        });
        par.querySelector('#pfr-view-pattern').addEventListener('click', ()=>{
          if(viewMode==='pattern') return;
          stopAnim(); viewMode = 'pattern'; syncControls(); draw();
        });

        par.querySelector('#pfr-reset').addEventListener('click', ()=>{
          fx=$(w.fx); fy=$(w.fy); fz=$(w.fz);
          rdx=$(w.rdx); rdy=$(w.rdy); rdz=$(w.rdz); rdist=$(w.rdist);
          alpha=$(w.alpha); beta=$(w.beta); period=$(w.period); tPhase=0;
          syncControls(); draw(); commit();
        });

        // clickable seismogram legend -- toggle a series on/off, or hover for a pointer cursor
        seisCv.addEventListener('click', e=>{
          const rect = seisCv.getBoundingClientRect();
          const hit = seisLegendHit(e.clientX-rect.left, e.clientY-rect.top);
          if(hit){ seisShow[hit.key] = !seisShow[hit.key]; drawSeismogram(); }
        });
        seisCv.addEventListener('mousemove', e=>{
          const rect = seisCv.getBoundingClientRect();
          seisCv.style.cursor = seisLegendHit(e.clientX-rect.left, e.clientY-rect.top) ? 'pointer' : 'default';
        });

        syncControls(); draw();
        }
        </script>

        """)
    end

    const _pfr_ready = true
end

# ╔═╡ 40f2f5d6-1f06-4828-aaa7-ef7b7e42d90e
begin
    _pfr_ready
    WideCell(@bind pfr PointForceRadiationInput(); max_width=1500)
end

# ╔═╡ 1e9f4218-324c-411e-8394-6b79201dbde3
# The bond starts as `nothing` until the widget's first real interaction in a live browser
# reports back -- fall back to the same defaults the widget itself opens with.
pfr_safe = pfr isa AbstractDict ? pfr : Dict{String,Any}(
    "fx" => 0.0, "fy" => 0.0, "fz" => 1.0, "rdx" => 0.6, "rdy" => 0.3, "rdz" => 0.74,
    "rdist" => 60.0, "alpha" => 4.0, "beta" => 2.0, "period" => 10.0)

# ╔═╡ 0999e269-5dd7-4ab5-b091-b5f4178d6ba8
let
    fhat = (pfr_safe["fx"], pfr_safe["fy"], pfr_safe["fz"])
    rhat = (pfr_safe["rdx"], pfr_safe["rdy"], pfr_safe["rdz"])
    θdeg = round(acos(clamp(dot(fhat, rhat), -1.0, 1.0)) * 180 / π; digits=0)
    md"""
    Receiver at **$(θdeg)°** from the force axis, distance r = **$(round(pfr_safe["rdist"]; digits=1))** km ·
    α = **$(pfr_safe["alpha"])**, β = **$(pfr_safe["beta"])** km/s ·
    source: **Gaussian pulse, period T = $(pfr_safe["period"]) s**
    """
end

# ╔═╡ 7b330e11-d451-4992-bcdb-10c6dc770399
md"""
`PfrPush` does no physics -- it takes the already-computed seismogram + source-waveform table and
hands them to the *already-rendered* [`PointForceRadiationInput`](@ref) widget by dispatching a
browser `CustomEvent` at the widget's own `<div>`, the same pattern this repo's other widgets use
(see e.g. `FieldPush` in `fault-dislocation.jl`).
"""

# ╔═╡ 52162e1d-c250-4941-9a3b-ffed90ca972f
_pfr_flatten(v) = join(v, ",")

# ╔═╡ 9642334f-90f6-4b8b-a015-758ec6b5246b
begin
    struct PfrPush
        tgrid_seis::Any
        seis::Any
        tgrid_wave::Any
        wave::Any
        effomega::Float64
    end
    function Base.show(io::IO, ::MIME"text/html", p::PfrPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('pfrwidget');
        if(w){
          w.dispatchEvent(new CustomEvent('pfr-update', { detail: {
            tgridSeis: [$(_pfr_flatten(p.tgrid_seis))],
            nearRadial: [$(_pfr_flatten(p.seis.near_radial))],
            nearTransverse: [$(_pfr_flatten(p.seis.near_transverse))],
            pfarRadial: [$(_pfr_flatten(p.seis.pfar_radial))],
            sfarTransverse: [$(_pfr_flatten(p.seis.sfar_transverse))],
            tgridWave: [$(_pfr_flatten(p.tgrid_wave))],
            wave: [$(_pfr_flatten(p.wave))],
            effOmega: $(p.effomega),
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 48aa578a-e733-4384-a034-6824fef86bfe
_pfr_sourcetype = Val(:gaussian)

# ╔═╡ d5590b12-e38e-40c0-923a-7691c40fb559
_pfr_sigma2 = period_to_sigma2(pfr_safe["period"])

# ╔═╡ 6d3f9b8a-1a2e-4c9f-9b3a-2f6a7d8c5e10
_pfr_sourceparams = (_pfr_sigma2,)

# ╔═╡ e3d9ea0e-0098-46f1-87a1-113b62bcc12f
_pfr_tgrid_seis = range(-2.0, 1.3 * 100.0 / min(pfr_safe["alpha"], pfr_safe["beta"]); length=400)

# ╔═╡ 94a6b641-a3f8-414c-a840-164f4d25d4f3
_pfr_seis = pointforce_seismogram(
    (pfr_safe["fx"], pfr_safe["fy"], pfr_safe["fz"]),
    (pfr_safe["rdx"], pfr_safe["rdy"], pfr_safe["rdz"]),
    pfr_safe["rdist"], pfr_safe["alpha"], pfr_safe["beta"],
    _pfr_sourcetype, _pfr_sourceparams, _pfr_tgrid_seis)

# ╔═╡ 515b9d85-8a37-4791-89f5-442111bc83a6
_pfr_tgrid_wave = range(-100.0 / min(pfr_safe["alpha"], pfr_safe["beta"]) - 1,
    1.3 * 100.0 / min(pfr_safe["alpha"], pfr_safe["beta"]) + 1; length=800)

# ╔═╡ 6fe7c3f3-121a-44cc-a8ef-7fd887565535
_pfr_wave = sourcewave_deriv_table(_pfr_sourcetype, _pfr_sourceparams, _pfr_tgrid_wave)

# ╔═╡ 491afad3-1ea6-4942-a533-a0d559f80fd1
PfrPush(_pfr_tgrid_seis, _pfr_seis, _pfr_tgrid_wave, _pfr_wave, 1 / sqrt(_pfr_sigma2))

# ╔═╡ 321aa3aa-d29e-481a-bdcc-aa0655a72bda
md"""
### Moment Tensor Radiation
"""

# ╔═╡ c6616a5c-225d-4415-9155-c1f769c14eea
begin
    """
    	latlon_grid_directions(ntheta, nphi)

    A `(ntheta+1)*(nphi+1)` grid of unit directions covering the whole sphere -- `theta` (polar
    angle) the outer index, `phi` (azimuth) the inner. The widget's JS regenerates this exact
    grid independently (pure geometry, no physics dependence) so it can zip Julia's per-direction
    physics results back in by flat index, without ever sending the directions themselves over
    the wire.
    """
    function latlon_grid_directions(ntheta::Int, nphi::Int)
        dirs = Vector{NTuple{3,Float64}}(undef, (ntheta + 1) * (nphi + 1))
        k = 1
        for i in 0:ntheta, j in 0:nphi
            theta = π * i / ntheta
            phi = 2π * j / nphi
            dirs[k] = (sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
            k += 1
        end
        return dirs
    end

    """
    	mt_radial_pattern(M, γ)

    Far-field P radiation-pattern amplitude for a general symmetric moment tensor `M` in
    direction `γ` -- the signed scalar `` \\hat\\gamma^\\top M \\hat\\gamma ``, generalizing the
    point force's `` \\cos\\theta `` pattern (Aki & Richards eq. 4.29).
    """
    mt_radial_pattern(M::AbstractMatrix, γ::AbstractVector) = dot(γ, M * γ)

    """
    	mt_transverse_vector(M, γ)

    Far-field S radiation-pattern VECTOR for a general symmetric moment tensor `M` in direction
    `γ` -- `` M\\hat\\gamma - (\\hat\\gamma^\\top M \\hat\\gamma)\\hat\\gamma ``, the component of
    `` M\\gamma `` tangential to `γ`. Generalizes the point force's `` \\sin\\theta\\,\\hat{t} ``
    pattern; unlike the point-force case this is a genuine 3-vector, not a single scalar times a
    fixed tangential direction, since a general `M` has no single preferred axis.
    """
    function mt_transverse_vector(M::AbstractMatrix, γ::AbstractVector)
        Mγ = M * γ
        return Mγ .- dot(γ, Mγ) .* γ
    end

    """
    	mt_field_grid(M, ntheta, nphi)

    Evaluate [`mt_radial_pattern`](@ref) and [`mt_transverse_vector`](@ref) at every point of
    [`latlon_grid_directions`](@ref)`(ntheta, nphi)` -- this is the widget's entire field data,
    recomputed whenever the matrix changes and pushed to JS in one shot. Returns
    `(pradial, stangent)`, concretely-typed `Vector{Float64}`/`Vector{NTuple{3,Float64}}` indexed
    the same way as the direction grid.
    """
    function mt_field_grid(M::AbstractMatrix, ntheta::Int, nphi::Int)
        dirs = latlon_grid_directions(ntheta, nphi)
        pradial = Vector{Float64}(undef, length(dirs))
        stangent = Vector{NTuple{3,Float64}}(undef, length(dirs))
        for (k, d) in enumerate(dirs)
            γ = [d[1], d[2], d[3]]
            pradial[k] = mt_radial_pattern(M, γ)
            st = mt_transverse_vector(M, γ)
            stangent[k] = (st[1], st[2], st[3])
        end
        return pradial, stangent
    end

    """
    	double_couple_matrix(strike, dip, rake)

    Symmetric, traceless moment tensor for a pure double-couple (shear-slip) source -- angles in
    degrees, Aki & Richards eq. (4.85) convention (`x`=North, `y`=East, `z`=Down). Normalized to
    unit Frobenius norm so the widget's matrix values stay in a fixed, comparable range regardless
    of mechanism.
    """
    function double_couple_matrix(strike::Real, dip::Real, rake::Real)
        φ, δ, λ = deg2rad(strike), deg2rad(dip), deg2rad(rake)
        mxx = -(sin(δ) * cos(λ) * sin(2φ) + sin(2δ) * sin(λ) * sin(φ)^2)
        myy = sin(δ) * cos(λ) * sin(2φ) - sin(2δ) * sin(λ) * cos(φ)^2
        mzz = sin(2δ) * sin(λ)
        mxy = sin(δ) * cos(λ) * cos(2φ) + 0.5 * sin(2δ) * sin(λ) * sin(2φ)
        mxz = -(cos(δ) * cos(λ) * cos(φ) + cos(2δ) * sin(λ) * sin(φ))
        myz = -(cos(δ) * cos(λ) * sin(φ) - cos(2δ) * sin(λ) * cos(φ))
        M = Symmetric([mxx mxy mxz; mxy myy myz; mxz myz mzz])
        return M ./ norm(M)
    end
end

# ╔═╡ 3908f04c-36b5-4a27-b51a-e4affa9467fc
md"""
!!! correct "Checking that a double couple is really traceless"
	A pure shear dislocation conserves volume, so its moment tensor must be traceless
	(`` M_{xx}+M_{yy}+M_{zz}=0 ``) for *every* strike/dip/rake -- checked below, algebraically
	guaranteed by the formula but worth confirming numerically too.
"""

# ╔═╡ 5606eb83-7ee5-4589-885c-e35f04c3f944
let
    for (strike, dip, rake) in ((0.0, 90.0, 0.0), (30.0, 45.0, 60.0), (200.0, 70.0, -30.0), (355.0, 10.0, 170.0))
        M = double_couple_matrix(strike, dip, rake)
        @assert abs(tr(M)) < 1e-10
    end
    "double_couple_matrix produces a traceless (pure shear) moment tensor for every mechanism ✓"
end

# ╔═╡ a9a7cebc-ee46-4884-9a06-0b30f335a613
md"""
!!! correct "Lamé's theorem still holds for a moment tensor"
	The same finite-difference check used above for the point force, repeated for a moment-tensor
	far field: `curl(P)` and `div(S)` should again be small but nonzero, shrinking with distance.
"""

# ╔═╡ 420d095f-c215-4e3f-800f-3e2ea2692321
let
    α, β = 4.0, 2.0
    st, σ² = Val(:gaussian), 5.0
    M = double_couple_matrix(30.0, 60.0, 90.0)

    pfarM_vec(x, y, z, t) = begin
        r = norm([x, y, z])
        γ = [x, y, z] ./ r
        (γ .* mt_radial_pattern(M, γ)) ./ (4π * RHO * α^3 * r) .* sourcewave(st, t - r / α, σ²)
    end
    sfarM_vec(x, y, z, t) = begin
        r = norm([x, y, z])
        γ = [x, y, z] ./ r
        mt_transverse_vector(M, γ) ./ (4π * RHO * β^3 * r) .* sourcewave(st, t - r / β, σ²)
    end
    function curl_fd(F, x, y, z, t; h=1e-3)
        [(F(x, y + h, z, t)[3] - F(x, y - h, z, t)[3]) / 2h - (F(x, y, z + h, t)[2] - F(x, y, z - h, t)[2]) / 2h,
            (F(x, y, z + h, t)[1] - F(x, y, z - h, t)[1]) / 2h - (F(x + h, y, z, t)[3] - F(x - h, y, z, t)[3]) / 2h,
            (F(x, y, z + h, t)[2] - F(x, y, z - h, t)[2]) / 2h - (F(x, y + h, z, t)[1] - F(x, y - h, z, t)[1]) / 2h]
    end
    div_fd(F, x, y, z, t; h=1e-3) =
        (F(x + h, y, z, t)[1] - F(x - h, y, z, t)[1]) / 2h + (F(x, y + h, z, t)[2] - F(x, y - h, z, t)[2]) / 2h +
        (F(x, y, z + h, t)[3] - F(x, y, z - h, t)[3]) / 2h

    map([15.0, 40.0, 100.0]) do r
        x, y, z, t = r / sqrt(3), r / sqrt(3), r / sqrt(3), r / α
        curlP_rel = norm(curl_fd(pfarM_vec, x, y, z, t)) / max(norm(pfarM_vec(x, y, z, t)), 1e-300)
        divS_rel = abs(div_fd(sfarM_vec, x, y, z, t)) / max(norm(sfarM_vec(x, y, z, t)), 1e-300)
        (; r, curlP_rel=round(curlP_rel; sigdigits=3), divS_rel=round(divS_rel; sigdigits=3))
    end
end

# ╔═╡ 1f36f0f5-7692-434e-bc9b-bed27f6c5bf4
md"""
### The Moment Tensor Widget
"""

# ╔═╡ 81421235-9cd9-428a-8629-03000f640a3d
begin
    struct MomentTensorInput
        mxx::Float64
        myy::Float64
        mzz::Float64
        mxy::Float64
        mxz::Float64
        myz::Float64
        strike::Float64
        dip::Float64
        rake::Float64
    end
    function MomentTensorInput(; strike=30.0, dip=60.0, rake=90.0)
        M = double_couple_matrix(strike, dip, rake)
        MomentTensorInput(M[1, 1], M[2, 2], M[3, 3], M[1, 2], M[1, 3], M[2, 3], strike, dip, rake)
    end

    Base.get(w::MomentTensorInput) = Dict{String,Any}(
        "mxx" => w.mxx, "myy" => w.myy, "mzz" => w.mzz,
        "mxy" => w.mxy, "mxz" => w.mxz, "myz" => w.myz,
        "strike" => w.strike, "dip" => w.dip, "rake" => w.rake, "source" => "sdr")

    function Base.show(io::IO, ::MIME"text/html", w::MomentTensorInput)
        write(io, """
        <div id="mtwidget">
        <style>
        pluto-cell:has(#mtwidget) { width: min(80vw, 1500px) !important;
          margin-left: calc((100% - min(80vw, 1500px)) / 2) !important; }
        #mtwidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #mtwidget .mt-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #mtwidget .mt-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #mtwidget .mt-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #mtwidget .mt-workspace{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start}
        #mtwidget .mt-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #mtwidget .mt-caption{font-size:13px;color:#9ca3af;text-align:center;margin-top:4px}
        #mtwidget canvas{display:block;cursor:grab}
        #mtwidget .mt-controls{flex:0 0 280px;width:280px;display:flex;flex-direction:column;
          gap:8px;font:14px sans-serif}
        #mtwidget .mt-control-group{background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
        #mtwidget .mt-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #mtwidget .mt-control-row{display:grid;grid-template-columns:60px minmax(0,1fr) 44px;gap:6px;align-items:center;margin:5px 0}
        #mtwidget .mt-control-row label{font-size:13px;color:#9ca3af}
        #mtwidget .mt-control-row input[type=range]{width:100%;min-width:0}
        #mtwidget .mt-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #mtwidget .mt-actions{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
        #mtwidget .mt-readout{font-size:13px;color:#d1d5db;line-height:1.6}
        #mtwidget .mt-readout b{color:#e5e7eb}
        #mtwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:13px;cursor:pointer}
        #mtwidget button.active{background:#2563eb;border-color:#93c5fd}
        #mtwidget .mt-jmatrix{display:inline-grid;grid-template-columns:repeat(3,auto);gap:4px 6px;padding:4px 8px;position:relative;margin:4px auto}
        #mtwidget .mt-jmatrix::before, #mtwidget .mt-jmatrix::after{content:'';position:absolute;top:2px;bottom:2px;width:6px;border:2px solid #9ca3af}
        #mtwidget .mt-jmatrix::before{left:0;border-right:none}
        #mtwidget .mt-jmatrix::after{right:0;border-left:none}
        #mtwidget .mt-jmatrix input{width:52px;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:4px 3px;font-size:12px;text-align:center;font-variant-numeric:tabular-nums}
        #mtwidget .mt-jmatrix input:focus{outline:2px solid #38bdf8;border-color:#38bdf8}
        #mtwidget .mt-jhint{font-size:11px;color:#6b7280;text-align:center;margin-top:4px}
        </style>

        <div class="mt-title">
          <div class="mt-title-desc">A fault has no net force -- its source is a moment tensor, and Lamé's theorem still splits it into radial P and tangential S.</div>
          <div class="mt-title-hint">edit the matrix, or drag strike/dip/rake for a double couple &middot; drag empty space to rotate &middot; press Play to expand the wavefronts</div>
        </div>

        <div class="mt-workspace">
          <div>
            <div class="mt-panel"><canvas id="mt-canvas"></canvas></div>
            <div class="mt-caption" id="mt-caption"></div>
          </div>
          <div class="mt-controls">
            <div class="mt-control-group">
              <div class="mt-control-title">View</div>
              <div class="mt-actions">
                <button id="mt-view-anim" type="button">Wavefront animation</button>
                <button id="mt-view-pattern" type="button">Radiation pattern</button>
              </div>
            </div>
            <div class="mt-control-group">
              <div class="mt-control-title">Moment tensor M</div>
              <div class="mt-jmatrix">
                <input type="number" step="0.05" id="mt-mxx" value="$(w.mxx)">
                <input type="number" step="0.05" id="mt-mxy" value="$(w.mxy)">
                <input type="number" step="0.05" id="mt-mxz" value="$(w.mxz)">
                <input type="number" step="0.05" id="mt-myx" value="$(w.mxy)">
                <input type="number" step="0.05" id="mt-myy" value="$(w.myy)">
                <input type="number" step="0.05" id="mt-myz" value="$(w.myz)">
                <input type="number" step="0.05" id="mt-mzx" value="$(w.mxz)">
                <input type="number" step="0.05" id="mt-mzy" value="$(w.myz)">
                <input type="number" step="0.05" id="mt-mzz" value="$(w.mzz)">
              </div>
              <div class="mt-jhint">symmetric -- editing an off-diagonal cell updates its mirror</div>
            </div>
            <div class="mt-control-group">
              <div class="mt-control-title">Double-couple generator</div>
              <div class="mt-control-row"><label>strike</label><input type="range" id="mt-strike" min="0" max="360" step="1" value="$(w.strike)"><span class="mt-value" id="mt-strike-v"></span></div>
              <div class="mt-control-row"><label>dip</label><input type="range" id="mt-dip" min="0" max="90" step="1" value="$(w.dip)"><span class="mt-value" id="mt-dip-v"></span></div>
              <div class="mt-control-row"><label>rake</label><input type="range" id="mt-rake" min="-180" max="180" step="1" value="$(w.rake)"><span class="mt-value" id="mt-rake-v"></span></div>
            </div>
            <div class="mt-control-group">
              <div class="mt-control-title">Wavefront</div>
              <div class="mt-actions"><button id="mt-play" type="button">Play</button><button id="mt-reset" type="button">Reset</button></div>
            </div>
            <div class="mt-control-group">
              <div class="mt-control-title">Readouts</div>
              <div class="mt-readout" id="mt-readout"></div>
            </div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        const WORLD = 1.0;

        let mxx=$(w.mxx), myy=$(w.myy), mzz=$(w.mzz), mxy=$(w.mxy), mxz=$(w.mxz), myz=$(w.myz);
        let strike=$(w.strike), dip=$(w.dip), rake=$(w.rake);
        let pendingSource = 'sdr';

        const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1500);
        const CONTROLS_W = 280+16;
        const SEC = Math.max(400, Math.min(availW - CONTROLS_W, 575));
        const DPR = window.devicePixelRatio || 1;

        const cv = par.querySelector('#mt-canvas');
        const ctx = cv.getContext('2d');
        function hidpi(canvas, context, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          context.setTransform(DPR,0,0,DPR,0,0);
        }
        hidpi(cv, ctx, SEC, SEC);

        // ---- 3D camera: identical orthographic drag-to-rotate family used throughout this
        // repo's widgets (StressTensorInput, PointForceRadiationInput, ...) -- pure display
        // geometry, no physics content, safe to keep client-side.
        let yaw = -0.6, pitch = 0.35;
        const CX = SEC/2, CY = SEC/2, RPIX = SEC*0.36;
        function rot(p){
          const x = p[0]*Math.cos(yaw) - p[1]*Math.sin(yaw);
          const y = p[0]*Math.sin(yaw) + p[1]*Math.cos(yaw);
          const z = p[2];
          return [x, y*Math.cos(pitch) - z*Math.sin(pitch), y*Math.sin(pitch) + z*Math.cos(pitch)];
        }
        function projUnit(p){ const r = rot(p); return [CX + r[0]*RPIX, CY - r[2]*RPIX]; }

        function drawAxes(){
          const AXLEN = 1.18;
          const axes = [ [[1,0,0],'x₁'], [[0,1,0],'x₂'], [[0,0,1],'x₃'] ];
          ctx.lineWidth = 1;
          ctx.font = '13px sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
          for(const [v,label] of axes){
            ctx.strokeStyle = 'rgba(156,163,175,0.55)';
            const [x0,y0] = projUnit([-v[0]*AXLEN,-v[1]*AXLEN,-v[2]*AXLEN]);
            const [x1,y1] = projUnit([v[0]*AXLEN,v[1]*AXLEN,v[2]*AXLEN]);
            ctx.beginPath(); ctx.moveTo(x0,y0); ctx.lineTo(x1,y1); ctx.stroke();
            ctx.fillStyle = '#9ca3af'; ctx.fillText(label, x1, y1);
          }
          ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
        }

        function drawStarMarker(cx, cy, r, fill, stroke){
          const spikes = 5, rOuter = r, rInner = r * 0.45;
          ctx.beginPath();
          for(let i=0; i<spikes*2; i++){
            const rad = i % 2 === 0 ? rOuter : rInner;
            const ang = -Math.PI/2 + i*Math.PI/spikes;
            const x = cx + rad*Math.cos(ang), y = cy + rad*Math.sin(ang);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.closePath();
          ctx.fillStyle = fill; ctx.fill();
          ctx.strokeStyle = stroke; ctx.lineWidth = 1; ctx.stroke();
        }

        function drawArrowHead2D(x0, y0, x1, y1, color, size){
          size = size || 6;
          const ang = Math.atan2(y1-y0, x1-x0);
          ctx.beginPath(); ctx.moveTo(x1,y1);
          ctx.lineTo(x1-size*Math.cos(ang-0.4), y1-size*Math.sin(ang-0.4));
          ctx.lineTo(x1-size*Math.cos(ang+0.4), y1-size*Math.sin(ang+0.4));
          ctx.closePath(); ctx.fillStyle = color; ctx.fill();
        }

        // ---- data pushed from Julia: mt_field_grid's per-direction P/S amplitudes, and the
        // source-derivative lookup table -- this widget computes NO physics client-side, only
        // geometry (camera, direction-grid regeneration) and drawing.
        let fieldData = null; // {ntheta, nphi, pradial, stx, sty, stz}
        let waveData = null;  // {tgridWave, wave}

        function latlonDirs(ntheta, nphi){
          const dirs = [];
          for(let i=0;i<=ntheta;i++){
            const theta = Math.PI*i/ntheta;
            for(let j=0;j<=nphi;j++){
              const phi = 2*Math.PI*j/nphi;
              dirs.push([Math.sin(theta)*Math.cos(phi), Math.sin(theta)*Math.sin(phi), Math.cos(theta)]);
            }
          }
          return dirs;
        }

        function lookupWave(t){
          if(!waveData) return 0;
          const {tgridWave, wave} = waveData;
          const n = tgridWave.length;
          const t0 = tgridWave[0], t1 = tgridWave[n-1];
          if(t <= t0) return wave[0];
          if(t >= t1) return wave[n-1];
          const frac = (t-t0)/(t1-t0)*(n-1);
          const i0 = Math.floor(frac), i1 = Math.min(n-1, i0+1);
          const f = frac - i0;
          return wave[i0]*(1-f) + wave[i1]*f;
        }
        function peakWaveNear(){
          if(!waveData) return 0;
          let best = 0, bestAbs = -1;
          for(let i=0; i<waveData.wave.length; i++){
            const v = waveData.wave[i];
            if(Math.abs(v) > bestAbs){ bestAbs = Math.abs(v); best = v; }
          }
          return best;
        }

        let refAmp = 1e-12;
        let playing = false, rafId = null, tPhase = 0, lastTs = null;
        let viewMode = 'anim';
        const ALPHA_MT = 0.55, BETA_MT = 0.32; // fixed display wave speeds -- shell expansion rate only, not physics content

        function drawWavefrontShell(R, dirs, wireColor, radial){
          if(!(R > 0.02 && R < WORLD*1.05)) return;
          ctx.strokeStyle = wireColor; ctx.lineWidth = 1;
          const {ntheta, nphi} = fieldData;
          for(let i=0;i<=ntheta;i++){
            ctx.beginPath(); let started=false;
            for(let j=0;j<=nphi;j++){
              const k = i*(nphi+1)+j;
              const g = dirs[k];
              const p = [g[0]*R, g[1]*R, g[2]*R];
              if(rot(p)[1] < 0){ started=false; continue; }
              const s = projUnit(p);
              if(!started){ ctx.moveTo(s[0],s[1]); started=true; } else ctx.lineTo(s[0],s[1]);
            }
            ctx.stroke();
          }

          const srcAtFront = peakWaveNear();
          const ARROW_LEN_MAX = WORLD*0.09;
          for(let k=0;k<dirs.length;k++){
            const g = dirs[k];
            if(rot(g)[1] < 0.02) continue;
            let amp, dir;
            if(radial){
              amp = fieldData.pradial[k]*srcAtFront;
              dir = g;
            } else {
              const sx=fieldData.stx[k], sy=fieldData.sty[k], sz=fieldData.stz[k];
              const mag = Math.hypot(sx,sy,sz);
              amp = mag*srcAtFront;
              dir = mag>1e-9 ? [sx/mag, sy/mag, sz/mag] : [0,0,0];
            }
            const mrel = Math.min(1, Math.abs(amp)/refAmp);
            if(mrel < 0.05) continue;
            const len = mrel*ARROW_LEN_MAX*Math.sign(amp);
            const base = [g[0]*R, g[1]*R, g[2]*R];
            const tip = [base[0]+dir[0]*len, base[1]+dir[1]*len, base[2]+dir[2]*len];
            const [x0,y0] = projUnit(base), [x1,y1] = projUnit(tip);
            const color = amp >= 0 ? 'rgba(239,68,68,' + (0.4+0.55*mrel) + ')' : 'rgba(59,130,246,' + (0.4+0.55*mrel) + ')';
            ctx.strokeStyle = color; ctx.lineWidth = 1.6;
            ctx.beginPath(); ctx.moveTo(x0,y0); ctx.lineTo(x1,y1); ctx.stroke();
            drawArrowHead2D(x0,y0,x1,y1,color,4);
          }
        }

        function drawWavefronts(){
          if(!fieldData) return;
          const dirs = latlonDirs(fieldData.ntheta, fieldData.nphi);
          drawWavefrontShell(ALPHA_MT*tPhase, dirs, 'rgba(249,115,22,0.22)', true);
          drawWavefrontShell(BETA_MT*tPhase, dirs, 'rgba(34,211,238,0.22)', false);
        }

        function drawRadiationSurfaces(){
          if(!fieldData) return;
          const {ntheta, nphi, pradial, stx, sty, stz} = fieldData;
          const dirs = latlonDirs(ntheta, nphi);
          const quads = [];
          function vertex(k, r){ const g = dirs[k]; return [g[0]*r, g[1]*r, g[2]*r]; }
          for(let i=0;i<ntheta;i++){
            for(let j=0;j<nphi;j++){
              const k00=i*(nphi+1)+j, k01=i*(nphi+1)+j+1, k11=(i+1)*(nphi+1)+j+1, k10=(i+1)*(nphi+1)+j;
              // P lobe (radial amplitude)
              const pAvg = (pradial[k00]+pradial[k01]+pradial[k11]+pradial[k10])/4;
              const pr = 0.08 + Math.abs(pAvg)*0.6;
              const pverts = [vertex(k00,pr),vertex(k01,pr),vertex(k11,pr),vertex(k10,pr)];
              const pctr = [(pverts[0][0]+pverts[1][0]+pverts[2][0]+pverts[3][0])/4,(pverts[0][1]+pverts[1][1]+pverts[2][1]+pverts[3][1])/4,(pverts[0][2]+pverts[1][2]+pverts[2][2]+pverts[3][2])/4];
              const pcolor = pAvg >= 0 ? 'rgba(239,68,68,' + (0.3+0.6*Math.min(1,Math.abs(pAvg)*1.4)) + ')' : 'rgba(59,130,246,' + (0.3+0.6*Math.min(1,Math.abs(pAvg)*1.4)) + ')';
              quads.push({verts: pverts, depth: rot(pctr)[1], color: pcolor});
              // S lobe (tangential magnitude)
              const smag = (Math.hypot(stx[k00],sty[k00],stz[k00])+Math.hypot(stx[k01],sty[k01],stz[k01])+Math.hypot(stx[k11],sty[k11],stz[k11])+Math.hypot(stx[k10],sty[k10],stz[k10]))/4;
              const sr = 0.08 + smag*0.5;
              const sverts = [vertex(k00,sr),vertex(k01,sr),vertex(k11,sr),vertex(k10,sr)];
              const sctr = [(sverts[0][0]+sverts[1][0]+sverts[2][0]+sverts[3][0])/4,(sverts[0][1]+sverts[1][1]+sverts[2][1]+sverts[3][1])/4,(sverts[0][2]+sverts[1][2]+sverts[2][2]+sverts[3][2])/4];
              const scolor = 'rgba(34,211,238,' + (0.2+0.55*Math.min(1,smag*1.4)) + ')';
              quads.push({verts: sverts, depth: rot(sctr)[1], color: scolor});
            }
          }
          quads.sort((a,b)=>a.depth-b.depth);
          for(const q of quads){
            ctx.beginPath();
            const [x0,y0] = projUnit(q.verts[0]); ctx.moveTo(x0,y0);
            for(let k=1;k<4;k++){ const [x,y] = projUnit(q.verts[k]); ctx.lineTo(x,y); }
            ctx.closePath();
            ctx.fillStyle = q.color; ctx.fill();
          }
        }

        function draw(){
          ctx.clearRect(0,0,SEC,SEC);
          drawAxes();
          if(viewMode === 'pattern'){ drawRadiationSurfaces(); } else { drawWavefronts(); }
          const [sx0,sy0] = projUnit([0,0,0]);
          drawStarMarker(sx0, sy0, 7, '#e5e7eb', '#000');
          ctx.fillStyle = '#9ca3af'; ctx.font = '12px sans-serif'; ctx.textAlign='left';
          ctx.fillText('source', sx0+10, sy0-8);
          par.querySelector('#mt-caption').textContent = 'drag to rotate';
          updateReadout();
        }

        function updateReadout(){
          const trace = mxx+myy+mzz;
          par.querySelector('#mt-readout').innerHTML =
            'trace(M) <b>' + trace.toFixed(2) + '</b>' + (Math.abs(trace)<0.05 ? ' (pure shear)' : '') + '<br>' +
            '<span style="color:#ef4444">&#9679;</span> push/outward &nbsp; <span style="color:#3b82f6">&#9679;</span> pull/inward';
        }

        function syncControls(){
          par.querySelector('#mt-view-anim').classList.toggle('active', viewMode==='anim');
          par.querySelector('#mt-view-pattern').classList.toggle('active', viewMode==='pattern');
          par.querySelector('#mt-play').style.display = viewMode==='anim' ? '' : 'none';
          par.querySelector('#mt-mxx').value = mxx.toFixed(3);
          par.querySelector('#mt-myy').value = myy.toFixed(3);
          par.querySelector('#mt-mzz').value = mzz.toFixed(3);
          par.querySelector('#mt-mxy').value = mxy.toFixed(3);
          par.querySelector('#mt-myx').value = mxy.toFixed(3);
          par.querySelector('#mt-mxz').value = mxz.toFixed(3);
          par.querySelector('#mt-mzx').value = mxz.toFixed(3);
          par.querySelector('#mt-myz').value = myz.toFixed(3);
          par.querySelector('#mt-mzy').value = myz.toFixed(3);
          par.querySelector('#mt-strike').value = strike; par.querySelector('#mt-strike-v').textContent = strike.toFixed(0)+'°';
          par.querySelector('#mt-dip').value = dip; par.querySelector('#mt-dip-v').textContent = dip.toFixed(0)+'°';
          par.querySelector('#mt-rake').value = rake; par.querySelector('#mt-rake-v').textContent = rake.toFixed(0)+'°';
        }

        let commitInFlight = false;
        function commit(){
          commitInFlight = true;
          par.value = { mxx, myy, mzz, mxy, mxz, myz, strike, dip, rake, source: pendingSource };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledCommit(){ if(!commitInFlight) commit(); }

        par.addEventListener('mt-update', e=>{
          fieldData = { ntheta: e.detail.ntheta, nphi: e.detail.nphi, pradial: e.detail.pradial,
            stx: e.detail.stx, sty: e.detail.sty, stz: e.detail.stz };
          waveData = { tgridWave: e.detail.tgridWave, wave: e.detail.wave };
          mxx=e.detail.mxx; myy=e.detail.myy; mzz=e.detail.mzz;
          mxy=e.detail.mxy; mxz=e.detail.mxz; myz=e.detail.myz;
          const maxAbsWave = Math.max(1e-12, ...waveData.wave.map(Math.abs));
          const maxAbsField = Math.max(1e-12, ...fieldData.pradial.map(Math.abs));
          refAmp = maxAbsWave*maxAbsField;
          commitInFlight = false;
          syncControls(); draw();
        });

        // ---- camera drag (empty space only -- no draggable handles on this widget, the matrix
        // and sliders are the input) ----
        let dragging = false, lastX = 0, lastY = 0;
        cv.addEventListener('mousedown', e=>{
          const rect = cv.getBoundingClientRect();
          dragging = true; lastX = e.clientX-rect.left; lastY = e.clientY-rect.top;
        });
        window.addEventListener('mousemove', e=>{
          if(!dragging) return;
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left, py = e.clientY-rect.top;
          const dx = px-lastX, dy = py-lastY;
          yaw += dx*0.008;
          pitch = Math.max(-1.4, Math.min(1.4, pitch - dy*0.008));
          lastX = px; lastY = py;
          draw();
        });
        window.addEventListener('mouseup', ()=>{ dragging = false; });

        // ---- matrix editor: 9 inputs, off-diagonal pairs mirrored ----
        const MIRROR = { 'mt-mxy':'mt-myx', 'mt-myx':'mt-mxy', 'mt-mxz':'mt-mzx', 'mt-mzx':'mt-mxz', 'mt-myz':'mt-mzy', 'mt-mzy':'mt-myz' };
        const FIELD_OF = { 'mt-mxx':'mxx','mt-myy':'myy','mt-mzz':'mzz',
          'mt-mxy':'mxy','mt-myx':'mxy','mt-mxz':'mxz','mt-mzx':'mxz','mt-myz':'myz','mt-mzy':'myz' };

        par.addEventListener('input', e=>{
          if(e.target===par) return;
          e.stopImmediatePropagation();
          const id = e.target.id, v = +e.target.value;
          if(id in FIELD_OF){
            const field = FIELD_OF[id];
            if(field==='mxx') mxx=v; else if(field==='myy') myy=v; else if(field==='mzz') mzz=v;
            else if(field==='mxy') mxy=v; else if(field==='mxz') mxz=v; else if(field==='myz') myz=v;
            if(MIRROR[id]) par.querySelector('#'+MIRROR[id]).value = v;
            pendingSource = 'matrix';
            throttledCommit();
            return;
          }
          if(id==='mt-strike'){ strike=v; par.querySelector('#mt-strike-v').textContent=strike.toFixed(0)+'°'; }
          else if(id==='mt-dip'){ dip=v; par.querySelector('#mt-dip-v').textContent=dip.toFixed(0)+'°'; }
          else if(id==='mt-rake'){ rake=v; par.querySelector('#mt-rake-v').textContent=rake.toFixed(0)+'°'; }
          else return;
          pendingSource = 'sdr';
          throttledCommit();
        }, true);

        par.querySelector('#mt-view-anim').addEventListener('click', ()=>{
          if(viewMode==='anim') return;
          viewMode = 'anim'; syncControls(); draw();
        });
        function stopAnim(){
          playing = false; par.querySelector('#mt-play').textContent = 'Play';
          if(rafId){ cancelAnimationFrame(rafId); rafId = null; }
        }
        par.querySelector('#mt-view-pattern').addEventListener('click', ()=>{
          if(viewMode==='pattern') return;
          stopAnim(); viewMode = 'pattern'; syncControls(); draw();
        });

        function tLoopMax(){ return 1.3*WORLD/Math.min(ALPHA_MT,BETA_MT); }
        const SIM_SPEED = tLoopMax() / 13;
        const playBtn = par.querySelector('#mt-play');
        function stepAnim(ts){
          if(lastTs===null) lastTs = ts;
          const dt = Math.min(0.1, (ts-lastTs)/1000);
          lastTs = ts;
          tPhase += dt*SIM_SPEED;
          if(tPhase > tLoopMax()) tPhase = 0;
          draw();
          rafId = requestAnimationFrame(stepAnim);
        }
        playBtn.addEventListener('click', ()=>{
          playing = !playing;
          playBtn.textContent = playing ? 'Pause' : 'Play';
          if(playing){ lastTs=null; rafId = requestAnimationFrame(stepAnim); }
          else if(rafId){ cancelAnimationFrame(rafId); rafId=null; }
        });
        par.querySelector('#mt-reset').addEventListener('click', ()=>{
          strike=$(w.strike); dip=$(w.dip); rake=$(w.rake); tPhase=0;
          pendingSource = 'sdr';
          syncControls(); draw(); commit();
        });

        syncControls(); draw();
        }
        </script>

        """)
    end

    const _mt_ready = true
end

# ╔═╡ 73444fbb-3a43-442b-a85d-29d243644a92
begin
    _mt_ready
    WideCell(@bind mt MomentTensorInput(); max_width=1500)
end

# ╔═╡ 61dea2a2-6d32-4bca-ba7e-1cf39ce83e35
# The bond starts as `nothing` until the widget's first real interaction in a live browser
# reports back -- fall back to a concrete double-couple example so the very first render (and
# the physics cells below) have real numbers to work with.
mt_safe = mt isa AbstractDict ? mt : let
    M0 = double_couple_matrix(30.0, 60.0, 90.0)
    Dict{String,Any}(
        "mxx" => M0[1, 1], "myy" => M0[2, 2], "mzz" => M0[3, 3],
        "mxy" => M0[1, 2], "mxz" => M0[1, 3], "myz" => M0[2, 3],
        "strike" => 30.0, "dip" => 60.0, "rake" => 90.0, "source" => "sdr")
end

# ╔═╡ ae12217e-0d9c-4156-8d36-3c7d8274cc34
# "source" tells us whether the strike/dip/rake sliders or the matrix cells themselves were the
# last thing the user touched -- the matrix is the actual source of truth for the physics either
# way, but only one of the two inputs can be authoritative for any given commit.
_mt_M = get(mt_safe, "source", "sdr") == "matrix" ?
        Symmetric([mt_safe["mxx"] mt_safe["mxy"] mt_safe["mxz"]
            mt_safe["mxy"] mt_safe["myy"] mt_safe["myz"]
            mt_safe["mxz"] mt_safe["myz"] mt_safe["mzz"]]) :
        double_couple_matrix(mt_safe["strike"], mt_safe["dip"], mt_safe["rake"])

# ╔═╡ 12685459-3bc9-44c7-b056-a8b913bd974c
_mt_pradial, _mt_stangent = mt_field_grid(_mt_M, _mt_ntheta, _mt_nphi)

# ╔═╡ 1e745b82-3ff6-49cb-a52e-a71fca6434ca
md"""
`MtPush` mirrors `PfrPush`: it carries no physics of its own, just the already-computed field
grid and source-waveform table, handed to the *already-rendered* [`MomentTensorInput`](@ref)
widget via a `CustomEvent`.
"""

# ╔═╡ 0bffa9c8-e05d-47f9-9bd1-4a096cd99747
begin
    struct MtPush
        ntheta::Int
        nphi::Int
        mxx::Float64
        myy::Float64
        mzz::Float64
        mxy::Float64
        mxz::Float64
        myz::Float64
        pradial::Vector{Float64}
        stangent::Vector{NTuple{3,Float64}}
        tgrid_wave::Any
        wave::Any
    end
    function Base.show(io::IO, ::MIME"text/html", p::MtPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('mtwidget');
        if(w){
          w.dispatchEvent(new CustomEvent('mt-update', { detail: {
            ntheta: $(p.ntheta), nphi: $(p.nphi),
            mxx: $(p.mxx), myy: $(p.myy), mzz: $(p.mzz), mxy: $(p.mxy), mxz: $(p.mxz), myz: $(p.myz),
            pradial: [$(_pfr_flatten(p.pradial))],
            stx: [$(_pfr_flatten([s[1] for s in p.stangent]))],
            sty: [$(_pfr_flatten([s[2] for s in p.stangent]))],
            stz: [$(_pfr_flatten([s[3] for s in p.stangent]))],
            tgridWave: [$(_pfr_flatten(p.tgrid_wave))],
            wave: [$(_pfr_flatten(p.wave))],
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 05b19046-6389-41e8-97f9-c13084b65872
MtPush(_mt_ntheta, _mt_nphi, _mt_M[1, 1], _mt_M[2, 2], _mt_M[3, 3], _mt_M[1, 2], _mt_M[1, 3], _mt_M[2, 3],
    _mt_pradial, _mt_stangent, _mt_tgrid_wave, _mt_wave)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[compat]
PlutoUI = "~0.7.83"
SpecialFunctions = "~2.8.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "0eb42bd4276548eb4fd54b3732ee3f5404667bf5"

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

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "cd4b115137894ced9830a92bcdb95a6bd8f38880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.2"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

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

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

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
# ╠═021492fb-0d1a-4f6e-a656-bdf6a8e74a90
# ╠═cbc8f361-986c-424d-a649-6e841c072e72
# ╟─dfaa5e67-7d3d-492a-8b56-a66bd0cd3970
# ╟─f80221b6-eec6-4044-896b-68cca2bfcb6d
# ╟─40f2f5d6-1f06-4828-aaa7-ef7b7e42d90e
# ╟─1e9f4218-324c-411e-8394-6b79201dbde3
# ╟─0999e269-5dd7-4ab5-b091-b5f4178d6ba8
# ╟─baf5aa58-0d0f-4e2c-803e-b1a8d98580a4
# ╟─b954f277-9840-465e-926e-5c4963751487
# ╟─fc71ba84-35f3-4269-b8c9-3308782ae9df
# ╟─7f488eb9-9372-41b4-95ec-69abb753ef96
# ╟─73444fbb-3a43-442b-a85d-29d243644a92
# ╟─61dea2a2-6d32-4bca-ba7e-1cf39ce83e35
# ╟─ae12217e-0d9c-4156-8d36-3c7d8274cc34
# ╠═e62a5517-0584-42fd-9372-07d6d0dc8119
# ╠═12685459-3bc9-44c7-b056-a8b913bd974c
# ╠═672e6bb6-90a8-4a18-85a8-d854e797aa00
# ╠═970be60f-3c20-4f57-909d-c163b388a9dc
# ╠═05b19046-6389-41e8-97f9-c13084b65872
# ╟─ad81924a-96a4-486c-a24e-34fe50ca9fd3
# ╟─b32d43d5-7b7e-4c90-8660-a0328ef6bb60
# ╠═9d59cfce-5ac4-46d9-bf66-4defb8a03546
# ╟─0b0ea15d-957e-467a-beea-961284f67c63
# ╠═cbc35a5e-0cc7-4ece-aac0-8520759d70fa
# ╠═9b7a2e1c-4c1a-4a6e-8f3a-6a8c1d2f5e01
# ╠═9b7a2e1c-4c1a-4a6e-8f3a-6a8c1d2f5e02
# ╟─18b1d5bd-de79-41a4-9324-963d17a145f7
# ╠═153caeb6-34d2-4fa6-a0f5-7c76f101a04a
# ╟─92539fa5-60db-4e54-b02c-79da24b2d9b9
# ╠═ed9930ba-ba7b-4295-a4c7-beb820fd2d66
# ╟─54313d17-0de1-4d94-a9ed-948eadff5cf2
# ╠═6787dae7-467b-4bf6-a643-5f318b109a32
# ╠═27257ff9-d835-476a-923b-fa8dc57cbd69
# ╟─757a080d-c03f-4d8f-9bb8-e0eb8aa4eca0
# ╠═a3154298-98dd-4625-b8b0-f60d1070c31a
# ╟─7b330e11-d451-4992-bcdb-10c6dc770399
# ╠═52162e1d-c250-4941-9a3b-ffed90ca972f
# ╠═9642334f-90f6-4b8b-a015-758ec6b5246b
# ╠═48aa578a-e733-4384-a034-6824fef86bfe
# ╠═d5590b12-e38e-40c0-923a-7691c40fb559
# ╠═6d3f9b8a-1a2e-4c9f-9b3a-2f6a7d8c5e10
# ╠═e3d9ea0e-0098-46f1-87a1-113b62bcc12f
# ╠═94a6b641-a3f8-414c-a840-164f4d25d4f3
# ╠═515b9d85-8a37-4791-89f5-442111bc83a6
# ╠═6fe7c3f3-121a-44cc-a8ef-7fd887565535
# ╠═491afad3-1ea6-4942-a533-a0d559f80fd1
# ╟─321aa3aa-d29e-481a-bdcc-aa0655a72bda
# ╠═c6616a5c-225d-4415-9155-c1f769c14eea
# ╟─3908f04c-36b5-4a27-b51a-e4affa9467fc
# ╠═5606eb83-7ee5-4589-885c-e35f04c3f944
# ╟─a9a7cebc-ee46-4884-9a06-0b30f335a613
# ╠═420d095f-c215-4e3f-800f-3e2ea2692321
# ╟─1f36f0f5-7692-434e-bc9b-bed27f6c5bf4
# ╠═81421235-9cd9-428a-8629-03000f640a3d
# ╟─1e745b82-3ff6-49cb-a52e-a71fca6434ca
# ╠═0bffa9c8-e05d-47f9-9bd1-4a096cd99747
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
