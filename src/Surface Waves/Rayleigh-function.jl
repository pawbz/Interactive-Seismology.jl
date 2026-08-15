### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Rayleigh Function"
#> tags = ["surfacewaves"]
#> layout = "layout.jlhtml"
#> description = "Particles may go round and round, but don't worry, let the Rayleigh wave pass by! This notebook will help you understand the physics behind it."

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

# ╔═╡ 5b807994-416d-11ed-22ff-77b37ff3cfac
begin
    using Symbolics
    using SymbolicUtils
    using PlutoUI
    using PlutoTeachingTools
end

# ╔═╡ d921cf1c-ca8b-44d9-8d00-e1ed55267647
ChooseDisplayMode()

# ╔═╡ 915243cb-9619-4eab-8ad4-1d9b821e8d01
PlutoUI.TableOfContents()

# ╔═╡ 07ee82c2-d17b-4b85-9ef6-c107b76d0cfe
md"""
# Rayleigh Function
Rayleigh waves are made from coupled, evanescent P and SV waves that satisfy a stress-free surface. Start with the experiment below: play the wavefield to watch the coupled P and SV displacement roll along the free surface, find the physical zero of the Rayleigh function, and follow one particle through its depth-dependent elliptical orbit. Then use the appendix to derive the same surface-wave condition symbolically.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 318ac5d5-bd8b-417b-a919-b46ee555ca1a
md"""
### What to look for

- The non-zero root of the secular function is always below the shear speed `\beta`; the density cancels from this normalized homogeneous half-space equation.
- At the free surface, the particle orbit is retrograde. Increase depth to see its size decay and its sense of rotation change.
- Increasing angular frequency makes the same depth look deeper in wavelengths, so both the orbit and the wavefield attenuate more quickly with depth.
- In the wavefield panel, toggle P and SV off one at a time: neither component alone satisfies the free-surface condition — only their sum does.
"""

# ╔═╡ c8253470-7d3c-436f-bced-141c06786ff9
md"""
In this document, we denote 2D spatial coordinates, angular frequency, time and necessary differential operators.
"""

# ╔═╡ 8afdf9d3-951d-4e53-8478-d5a075acbbe3
@syms x::Real z::Real ω::Real t::Real

# ╔═╡ aafebaf4-9091-42f0-8e5c-cdae92c21203
begin
    Dx = Differential(x)
    Dz = Differential(z)
end

# ╔═╡ 91506b9c-af2f-44ab-a200-04ba88bc53df
@syms ı::Complex{Real} # imaginary unit, going to substitute with im later

# ╔═╡ 6eb3e9d6-2afa-450d-915c-81b0bfa92e9d
md"""
The medium is assumed to be homogeneous with uniform P- and S-wave velocities
`\alpha` and `\beta`. We use `p` to denote the horizontal slowness. As both the P- and SV-waves are inhomogeneous, a necessary condition is
```math
\frac{1}{\alpha} < \frac{1}{\beta} < p.
```
The Lame parameters are denoted using `\lambda` and `\mu`, and `\rho` is mass density.
"""

# ╔═╡ e19ff895-cee0-45ad-994d-26b91e7c4d3c
@syms α::Real p::Real β::Real ρ::Real λ::Real μ::Real

# ╔═╡ 3e9930db-b54c-45af-8bb7-fc86568ae504
tip(md"The inhomegeneous P- and SV-waves are coupled due to the fact that they share the same horizonal slowness `p`. If they are not coupled, one can show that the stress-free conditions cannot be satisfied.")

# ╔═╡ 9409762a-7803-4001-afbd-62c5d40dea2e
md"""
The vertical components of the slowness vectors (corresponding to `\alpha` and `\beta`) are imaginary. Here, we use `\eta` to denote the imaginary part of the vertical components of the slowness vector.
"""

# ╔═╡ 296dbce3-8e24-49b2-a8c9-5f59c6586160
@syms ηₚ::Real ηₛ::Real

# ╔═╡ a38858bd-d1ea-4918-8445-9b1306840979
md"The dispersion relation gives."

# ╔═╡ ca8f7333-8bb0-42f8-82d4-b99f1b03c2c9
ηₚ ~ sqrt(p^2 - 1 / α^2)

# ╔═╡ ce12c0d3-1601-431c-9283-a88a44b725b4
ηₛ ~ sqrt(p^2 - 1 / β^2)

# ╔═╡ 1213a9c2-0194-4f57-8b64-5f810b46ee67
@syms Aₚ A Aₛ

# ╔═╡ ef389316-a822-4a31-865e-5751c8bc6fe0
md"""
## Inhomogeneous Plane Waves
We now derive the expression for an inhomogeneous plane wave with amplitude `A` and frequency `ω`, propagating in the `x` direction with slowness `p`. These waves exhibit exponential decay with increasing distance from the surface at `z=0`.
"""

# ╔═╡ e88d9f22-d833-4d4c-b651-a7d2e318323e
inhomo_plane(p, η, A) = A * exp(-ω * η * z) * exp(ı * ω * (t - (p * x)))

# ╔═╡ 59c61d01-e0e3-405b-b048-1ceebee0d29b
md"""
We now consider two harmonic plane waves that satisfy the scalar wave equations for P and SV potentials `\phi` and `\psi`. `A_p` and `A_s` denote the amplitudes of these potentials.
"""

# ╔═╡ 0ffac267-6380-4288-8080-5f1ae30646d8
ϕ = inhomo_plane(p, ηₚ, Aₚ)

# ╔═╡ ff414220-762b-4e82-950c-833dd0282772
ψ = inhomo_plane(p, ηₛ, Aₛ)

# ╔═╡ 596dccc2-7853-4562-8f25-110108ff97c2
md"""
## Total Displacement
The particle displacement due to the P wave.
"""

# ╔═╡ a61d687c-5e13-423a-98b0-9b73371a8044
uₚ = expand_derivatives.([Dx(ϕ), Dz(ϕ)])

# ╔═╡ a254e7a7-d5ef-467a-98eb-f4bb1ceed97c
md"Similarly for SV wave, we have."

# ╔═╡ 019ac47c-a41f-4da2-adf5-fd1d6bac7f5d
uₛ = expand_derivatives.([-Dz(ψ), Dx(ψ)])

# ╔═╡ 28291690-152c-4bca-b194-f1aac4b102c7
md"The total displacement field is given by summing the P and SV components. We denote the x and z components of the total displacement field using `ux` and `uz`."

# ╔═╡ ff965bbb-61bc-4595-8082-fb7a8b4640ba
ux, uz = uₚ .+ uₛ

# ╔═╡ 127d393b-1cf5-4b4a-bb98-b271d0167b97
md"""
## Stress-free Condition
We now evaluate the stress components `σxz` and `σzz` on the free surface and apply the boundary conditions to estimate the ratio `Aₛ/Aₚ` using both boundary conditions. Subsequently, we derive a condition that must be satisfied for this ratio to be unique. This condition provides the Rayleigh function, and the roots of this function give the phase velocities of the Rayleigh waves.
"""

# ╔═╡ 0aa629dd-0201-463e-9765-3bc002cdc252
σxz = expand_derivatives(μ * (Dx(uz) + Dz(ux))) |> simplify

# ╔═╡ 83366f56-03a2-4134-8d94-8d54943ce3ff
σxz_z0 = substitute(σxz, [z => 0, ı * ı => -1]) |> simplify

# ╔═╡ 6f6af471-2d5f-44ee-bb4f-ccd441dd8248
σxz_z0 ~ 0

# ╔═╡ 6308d52d-809c-4f1f-b508-61388b1d6557
Aratio1 = Symbolics.symbolic_linear_solve(σxz_z0 ~ 0, Aₛ) / Aₚ |> simplify

# ╔═╡ a963b758-1ea5-4248-8544-88a612815d88
σzz = expand_derivatives(λ * (Dx(ux) + Dz(uz)) + 2 * μ * (Dz(uz))) |> simplify

# ╔═╡ 60409882-81c5-4645-85ac-4c40ac92327b
σzz_z0 = simplify(substitute(σzz, [z => 0, ı * ı => -1]))

# ╔═╡ d2a84a1e-6529-467b-905a-fdcc21135245
Aratio2 = Symbolics.symbolic_linear_solve(σzz_z0 ~ 0, Aₛ) / Aₚ |> simplify

# ╔═╡ 3daee10c-54d5-4071-b577-b4a8abd90976
md"""
## Rayleigh Function
For a non-trivial solution for `Aₚ` and `Aₛ`, the following equation must be satisfied.
"""

# ╔═╡ a22543fe-63cd-4cfb-8791-dfc2eb929359
R1 = diff(arguments(Aratio1) .* reverse(arguments(Aratio2))) |> first |> simplify

# ╔═╡ 288c375f-ff18-45a8-bc12-f5b3cc4fa377
md"""
We will reparameterize the Rayleigh function using `\alpha`, `\beta`, and `\rho`, and substitute `ı^2 = -1`. This reparameterized function must be satisfied for a given medium. The root of this function provides the horizontal slowness `p` at which both the inhomogeneous P and SV plane waves propagate.
"""

# ╔═╡ d8369538-ab73-4575-bf86-31b9111523cc
md"""
Towards plotting, we will now substitute the UI values of `\alpha`, `\beta` and `\rho`.
"""

# ╔═╡ a7c6f76d-fbcd-4fa6-a12f-868428cedef8
md"""
We will now find the root of the Rayleigh function using the bisection method.
"""

# ╔═╡ 8ad847c8-3d83-4c15-97f1-9b4490cbf2b0
md"""
## Rayleigh Waves
We shall now substitute the zero of the Rayleigh function into the expressions of inhomogeneous P and SV planewaves that are channeled along the free-surface.
"""

# ╔═╡ 8d54de78-38f7-4fd0-87e8-cb3ed6009b87
simplify(substitute(ux, [Aₚ => 1, Aₛ => Aratio1]))

# ╔═╡ 6bdcc2ac-3319-432d-a286-f706bf68234e
simplify(substitute(uz, [Aₚ => 1, Aₛ => Aratio1]))

# ╔═╡ 379feb55-32e5-4867-b108-b10c2c0782c3
md"""
## Appendix
"""

# ╔═╡ aa111199-0001-4001-8001-100000000001
md"""
### Symbolic Reparametrization
"""

# ╔═╡ 996bd1d1-09f9-4774-b73f-0f4e58f52455
"""
	subs_αβρ(x)

Reparametrize a symbolic expression `x` using the elastic speeds `α`, `β` and
density `ρ`, instead of the Lame parameters `λ`, `μ` and the vertical
slownesses `ηₚ`, `ηₛ`. Also substitutes `ı^2 => -1` so the result is a real
symbolic expression in `α`, `β`, `ρ` and `p`.
"""
function subs_αβρ(x)
    simplify(substitute(x, [ı * ı => -1, λ => α * α * ρ - 2 * β * β * ρ, μ => β * β * ρ, ηₚ * ηₚ => p * p - 1 / α^2, ηₛ * ηₛ => p * p - 1 / β^2, ηₚ => sqrt(p * p - 1 / α^2), ηₛ => sqrt(p * p - 1 / β^2), μ => β * β * ρ]))
end

# ╔═╡ 1a0819fb-b998-4fac-91d3-b3b375cdfc42
R = subs_αβρ(R1)

# ╔═╡ 18fa43a1-4e0d-4910-9125-d41d6d174ef1
"""
	subs_αβρ_plot(x, αp_val, βp_val, ρp_val)

Substitute the widget's current `α`, `β` and `ρ` slider values into an
expression already reparametrized by [`subs_αβρ`](@ref).
"""
function subs_αβρ_plot(x, αp_val, βp_val, ρp_val)
    simplify(substitute(subs_αβρ(x), [α => αp_val, β => βp_val, ρ => ρp_val]))
end

# ╔═╡ 3053571a-096c-448a-9b48-60f9a1b2fce6
"""
	subs_all_plot(x, αp_val, βp_val, ρp_val, cᵣ_val, ωp_val)

Substitute the widget's current elastic speeds, density, Rayleigh phase
speed `cᵣ_val` (as `p = 1/cᵣ_val`) and angular frequency into an expression
already reparametrized by [`subs_αβρ_plot`](@ref), also fixing `Aₚ = 1`.
"""
function subs_all_plot(x, αp_val, βp_val, ρp_val, cᵣ_val, ωp_val)
    simplify(substitute(subs_αβρ_plot(x, αp_val, βp_val, ρp_val), [p => inv(cᵣ_val), ı => im, ω => ωp_val, Aₚ => 1,]))
end

# ╔═╡ aa111199-0002-4002-8002-100000000002
md"""
### Numerical Rayleigh Function and Wavefield

The symbolic derivation above yields the Rayleigh function as an expression
in `α`, `β`, `ρ`, `p` and `ω`; density and the overall amplitude cancel out
of the root-finding condition, so the functions below reimplement the
dimensionless secular function directly and solve it numerically, then use
the same closed-form potentials to build the particle orbit and the 2-D
wavefield snapshots that drive the widget.
"""

# ╔═╡ bb111111-1111-4111-8111-111111111111
begin
    """The dimensionless Rayleigh secular function for an isotropic half-space."""
    function rayleigh_function(alpha::Real, beta::Real, c::Real)
        q = Float64(c) / Float64(beta)
        r = Float64(c) / Float64(alpha)
        (2 - q^2)^2 - 4 * sqrt(max(0.0, 1 - q^2)) * sqrt(max(0.0, 1 - r^2))
    end

    """Find the non-zero, sub-shear root by a bracketed bisection."""
    function rayleigh_speed(alpha::Real, beta::Real)
        a, b = Float64(alpha), Float64(beta)
        @assert a > b > 0 "Rayleigh waves require alpha > beta > 0"
        f(c) = rayleigh_function(a, b, c)
        left, fleft = 0.02 * b, f(0.02 * b) # exclude the trivial c = 0 root
        for right in range(0.03 * b, 0.999 * b; length = 400)
            fright = f(right)
            if fleft * fright <= 0
                lo, hi = left, right
                for _ in 1:64
                    mid = (lo + hi) / 2
                    f(lo) * f(mid) <= 0 ? (hi = mid) : (lo = mid)
                end
                return (lo + hi) / 2
            end
            left, fleft = right, fright
        end
        error("Could not bracket the non-zero Rayleigh root")
    end

    """Numerical P--SV particle orbit at one position and depth; physics stays in Julia."""
    function rayleigh_orbit(alpha::Real, beta::Real, omega::Real, depth::Real; samples::Int = 181)
        a, b, ω, z = Float64(alpha), Float64(beta), Float64(omega), Float64(depth)
        c = rayleigh_speed(a, b)
        p = inv(c)
        ηp = sqrt(p^2 - inv(a)^2)
        ηs = sqrt(p^2 - inv(b)^2)
        # The shear-traction condition sets this complex SV/P potential ratio.
        As_over_Ap = 2im * p * ηp / (ηs^2 + p^2)
        phases = range(0.0, 2π; length = samples)
        displacement(phase, zvalue) = begin
            ϕ = exp(-ω * ηp * zvalue) * cis(phase)
            ψ = As_over_Ap * exp(-ω * ηs * zvalue) * cis(phase)
            ux = real(-im * ω * p * ϕ + ω * ηs * ψ)
            uz = real(-ω * ηp * ϕ - im * ω * p * ψ)
            (ux, uz)
        end
        orbit = displacement.(phases, z)
        surface = displacement.(phases, 0.0)
        orbit_x = first.(orbit)
        orbit_z = last.(orbit)
        scale = max(maximum(hypot.(orbit_x, orbit_z)), eps(Float64))
        surface_scale = max(maximum(hypot.(first.(surface), last.(surface))), eps(Float64))
        return (
            c_r = c, phases = collect(phases), orbit_x = orbit_x ./ scale, orbit_z = orbit_z ./ scale,
            attenuation = scale / surface_scale,
        )
    end

    """
    	rayleigh_wavefield_frames(alpha, beta, omega; nx=26, nz=14, nframes=16, xmax=500.0, zmax=250.0)

    Evaluate the coupled inhomogeneous P- and SV-wave displacement field of a
    Rayleigh wave on a `(x, z)` grid, at `nframes` equally-spaced phases over
    one full temporal period. Returns the P-only and SV-only displacement
    components separately (both normalized by the same overall scale) so a
    viewer can toggle either contribution on or off without recomputation.
    """
    function rayleigh_wavefield_frames(alpha::Real, beta::Real, omega::Real;
        nx::Int = 26, nz::Int = 14, nframes::Int = 16, xmax::Real = 500.0, zmax::Real = 250.0)
        a, b, ω = Float64(alpha), Float64(beta), Float64(omega)
        c = rayleigh_speed(a, b)
        p = inv(c)
        ηp = sqrt(p^2 - inv(a)^2)
        ηs = sqrt(p^2 - inv(b)^2)
        As_over_Ap = 2im * p * ηp / (ηs^2 + p^2)

        xs = collect(range(0.0, xmax; length = nx))
        zs = collect(range(0.0, zmax; length = nz))
        ts = collect(range(0.0, 2π / ω; length = nframes + 1))[1:nframes]

        uxP = Array{Float64}(undef, nframes, nx, nz)
        uzP = similar(uxP)
        uxS = similar(uxP)
        uzS = similar(uxP)

        for (it, t) in enumerate(ts), (ix, x) in enumerate(xs), (iz, z) in enumerate(zs)
            ϕ = exp(-ω * ηp * z) * cis(ω * (t - p * x))
            ψ = As_over_Ap * exp(-ω * ηs * z) * cis(ω * (t - p * x))
            uxP[it, ix, iz] = real(-im * ω * p * ϕ)
            uzP[it, ix, iz] = real(-ω * ηp * ϕ)
            uxS[it, ix, iz] = real(ω * ηs * ψ)
            uzS[it, ix, iz] = real(-im * ω * p * ψ)
        end

        norm = max(maximum(hypot.(uxP .+ uxS, uzP .+ uzS)), eps())
        return (
            xs = xs, zs = zs, nx = nx, nz = nz, nframes = nframes,
            uxP = uxP ./ norm, uzP = uzP ./ norm, uxS = uxS ./ norm, uzS = uzS ./ norm,
        )
    end

    """
    	rayleigh_widget_payload(alpha, beta, omega, time_fraction, depth)

    Bundle everything the "Rayleigh secular function" and "Particle orbit"
    panels need for one widget update: the secular-function curve sampled
    below `\\beta`, the normalized particle orbit at `depth`, and the single
    snapshot point selected by `time_fraction` (a fraction of one cycle).
    """
    function rayleigh_widget_payload(alpha::Real, beta::Real, omega::Real, time_fraction::Real, depth::Real)
        c_r = rayleigh_speed(alpha, beta)
        curve_c = collect(range(0.0, 0.999 * Float64(beta); length = 240))
        curve_r = rayleigh_function.(Float64(alpha), Float64(beta), curve_c)
        orbit = rayleigh_orbit(alpha, beta, omega, depth)
        i = clamp(round(Int, Float64(time_fraction) * (length(orbit.phases) - 1)) + 1, 1, length(orbit.phases))
        return (
            curve_c = curve_c, curve_r = curve_r, orbit_x = orbit.orbit_x, orbit_z = orbit.orbit_z,
            snapshot_x = orbit.orbit_x[i], snapshot_z = orbit.orbit_z[i], c_r = c_r,
            surface_ratio = orbit.attenuation,
        )
    end
end

# ╔═╡ aa111199-0003-4003-8003-100000000003
md"""
### The Interactive Widget
"""

# ╔═╡ bb222222-2222-4222-8222-222222222222
begin
    struct RayleighFunctionInput
        alpha::Float64
        beta::Float64
        density::Float64
        omega::Float64
        time::Float64
        depth::Float64
    end

    RayleighFunctionInput(; alpha=6.0, beta=4.49, density=2.9, omega=0.15, time=0.25, depth=0.0) =
        RayleighFunctionInput(Float64(alpha), Float64(beta), Float64(density), Float64(omega), Float64(time), Float64(depth))

    Base.get(w::RayleighFunctionInput) = Dict{String,Any}(
        "alpha" => w.alpha, "beta" => w.beta, "density" => w.density,
        "omega" => w.omega, "time" => w.time, "depth" => w.depth,
    )

    function Base.show(io::IO, ::MIME"text/html", w::RayleighFunctionInput)
        write(io, """
        <div id="rf-widget"><style>
        #rf-widget{width:100%;max-width:1400px;margin:auto;color:#e5e7eb;font:14px system-ui,sans-serif}#rf-widget *{box-sizing:border-box}
        #rf-widget .rf-title{padding:10px 14px;margin-bottom:10px;text-align:center;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px}#rf-widget .rf-title b{font-size:17px}#rf-widget .rf-hint{margin-top:3px;color:#9ca3af;font-size:13px}
        #rf-widget .rf-workspace{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}#rf-widget .rf-panel,#rf-widget .rf-group{min-width:0;padding:10px;background:#050505;border:1px solid #2f3744;border-radius:6px}#rf-widget .rf-panel-title{font-size:15px;font-weight:700;color:#f3f4f6}.rf-caption{min-height:32px;margin:4px 0 7px;color:#9ca3af;font-size:13px;line-height:1.3}
        #rf-widget .rf-wide{grid-column:1/-1}#rf-widget canvas{display:block;width:100%;height:290px;background:#000;border:1px solid #374151;border-radius:4px}#rf-wave{height:300px!important}
        #rf-widget .rf-legend{display:flex;gap:12px;flex-wrap:wrap;align-items:center;margin-top:7px;font-size:12px}.rf-key{display:inline-flex;align-items:center;gap:5px;cursor:pointer;user-select:none}.rf-swatch{width:18px;height:3px;border-radius:2px;background:currentColor}.rf-dash{background:repeating-linear-gradient(90deg,currentColor 0 6px,transparent 6px 10px)}.rf-key input{accent-color:currentColor;margin:0}
        #rf-widget .rf-controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:8px;margin-top:10px}.rf-group-title{margin-bottom:7px;font-size:16px;font-weight:700;color:#f3f4f6}.rf-row{display:grid;grid-template-columns:minmax(76px,120px) minmax(70px,1fr) minmax(44px,72px);gap:7px;align-items:center;margin:7px 0;color:#d1d5db}#rf-widget input[type=range]{width:100%;min-width:0;accent-color:#f59e0b}.rf-value{color:#fbbf24;text-align:right;font-variant-numeric:tabular-nums}#rf-widget button.rf-btn{border:1px solid #9ca3af;border-radius:4px;padding:6px 12px;background:#606060;color:#f3f4f6;font-size:14px;cursor:pointer}#rf-widget button.rf-btn:hover{background:#737373}.rf-readout{margin-top:8px;color:#d1d5db;font-size:13px}.rf-readout b{color:#fbbf24}@media(max-width:760px){#rf-widget .rf-workspace{grid-template-columns:1fr}#rf-widget canvas{height:250px}}
        </style><div class="rf-title"><b>Find the surface-wave speed, then watch the coupled P&ndash;SV field roll by.</b><div class="rf-hint">Change the elastic speeds &middot; play the wavefield &middot; compare the secular zero with particle motion</div></div>
        <div class="rf-workspace">
        <section class="rf-panel rf-wide"><div class="rf-panel-title">Rayleigh wavefield</div><div class="rf-caption" id="rf-wave-caption">Snapshot of the coupled evanescent P and SV displacement field near the free surface.</div><canvas id="rf-wave"></canvas>
        <div class="rf-legend"><button id="rf-play" type="button" class="rf-btn">Play</button><label class="rf-key" style="color:#f59e0b"><input type="checkbox" id="rf-show-p" checked><i class="rf-swatch"></i>P</label><label class="rf-key" style="color:#60a5fa"><input type="checkbox" id="rf-show-s" checked><i class="rf-swatch"></i>SV</label><span class="rf-key" style="color:#ef4444"><i class="rf-swatch"></i>free surface</span></div>
        </section>
        <section class="rf-panel"><div class="rf-panel-title">Rayleigh secular function</div><div class="rf-caption">The physical, non-zero crossing determines the sub-shear Rayleigh phase speed.</div><canvas id="rf-secular"></canvas><div class="rf-legend"><span class="rf-key" style="color:#f59e0b"><i class="rf-swatch"></i>R(c)</span><span class="rf-key" style="color:#22d3ee"><i class="rf-swatch rf-dash"></i>cᵣ</span><span class="rf-key" style="color:#a78bfa"><i class="rf-swatch rf-dash"></i>β</span></div></section>
        <section class="rf-panel"><div class="rf-panel-title">Particle orbit at selected depth</div><div class="rf-caption">One normalized cycle at x = 0. The marker is the selected phase; the readout reports the physical attenuation.</div><canvas id="rf-orbit"></canvas><div class="rf-legend"><span class="rf-key" style="color:#60a5fa"><i class="rf-swatch"></i>particle path</span><span class="rf-key" style="color:#fbbf24"><i class="rf-swatch"></i>snapshot</span></div></section>
        </div>
        <div class="rf-controls"><section class="rf-group"><div class="rf-group-title">Isotropic half-space</div><label class="rf-row"><span>P speed α</span><input id="rf-alpha" type="range" min="4" max="10" step="0.1" value="$(w.alpha)"><span id="rf-alpha-v" class="rf-value"></span></label><label class="rf-row"><span>S speed β</span><input id="rf-beta" type="range" min="2" max="6" step="0.01" value="$(w.beta)"><span id="rf-beta-v" class="rf-value"></span></label><label class="rf-row"><span>Density ρ</span><input id="rf-density" type="range" min="2" max="4" step="0.1" value="$(w.density)"><span id="rf-density-v" class="rf-value"></span></label></section><section class="rf-group"><div class="rf-group-title">Observe the orbit</div><label class="rf-row"><span>Angular ω</span><input id="rf-omega" type="range" min="0.05" max="0.35" step="0.01" value="$(w.omega)"><span id="rf-omega-v" class="rf-value"></span></label><label class="rf-row"><span>Phase</span><input id="rf-time" type="range" min="0" max="1" step="0.01" value="$(w.time)"><span id="rf-time-v" class="rf-value"></span></label><label class="rf-row"><span>Depth z</span><input id="rf-depth" type="range" min="0" max="25" step="0.5" value="$(w.depth)"><span id="rf-depth-v" class="rf-value"></span></label></section><section class="rf-group"><div class="rf-group-title">Readout</div><div id="rf-readout" class="rf-readout">Computing Rayleigh orbit…</div><button id="rf-reset" type="button" class="rf-btn">Reset experiment</button></section></div></div>
        <script>
        (function(){
        if (window.__rfCleanup) { window.__rfCleanup() }
        const rfInstance = {};
        window.__rfCurrentInstance = rfInstance;
        const rfController = new AbortController();
        window.__rfCleanup = () => rfController.abort();
        const rfSignal = { signal: rfController.signal };

        const par=document.currentScript.previousElementSibling,by=id=>par.querySelector('#'+id);
        const state={alpha:$(w.alpha),beta:$(w.beta),density:$(w.density),omega:$(w.omega),time:$(w.time),depth:$(w.depth)};
        const ids={'rf-alpha':'alpha','rf-beta':'beta','rf-density':'density','rf-omega':'omega','rf-time':'time','rf-depth':'depth'};
        let result=null;
        let showP=true, showS=true, playing=false, frameIdx=0, lastFrameTime=0, rafId=null;
        function labels(){by('rf-alpha-v').textContent=state.alpha.toFixed(1)+' km/s';by('rf-beta-v').textContent=state.beta.toFixed(2)+' km/s';by('rf-density-v').textContent=state.density.toFixed(1)+' g/cm³';by('rf-omega-v').textContent=state.omega.toFixed(2)+' rad/s';by('rf-time-v').textContent=(360*state.time).toFixed(0)+'°';by('rf-depth-v').textContent=state.depth.toFixed(1)+' km'}
        function emit(){if(state.beta>=.98*state.alpha){state.beta=Math.max(2,.98*state.alpha);by('rf-beta').value=state.beta}par.value={...state};par.dispatchEvent(new CustomEvent('input'))}
        function setup(canvas){const rect=canvas.getBoundingClientRect(),d=devicePixelRatio||1,W=Math.max(1,rect.width),H=Math.max(1,rect.height);canvas.width=Math.round(W*d);canvas.height=Math.round(H*d);const x=canvas.getContext('2d');x.setTransform(d,0,0,d,0,0);x.fillStyle='#000';x.fillRect(0,0,W,H);return[x,W,H]}
        function axes(x,W,H,title,xleft,xright,ylow,yhigh){const m={l:42,r:14,t:28,b:30},px=v=>m.l+(v-xleft)/(xright-xleft)*(W-m.l-m.r),py=v=>H-m.b-(v-ylow)/(yhigh-ylow)*(H-m.t-m.b);x.strokeStyle='#1f2937';x.lineWidth=1;for(let i=0;i<5;i++){let y=m.t+i*(H-m.t-m.b)/4;x.beginPath();x.moveTo(m.l,y);x.lineTo(W-m.r,y);x.stroke()}if(ylow<0&&yhigh>0){x.strokeStyle='#4b5563';x.beginPath();x.moveTo(m.l,py(0));x.lineTo(W-m.r,py(0));x.stroke()}x.fillStyle='#e5e7eb';x.font='600 13px sans-serif';x.fillText(title,m.l,18);x.fillStyle='#9ca3af';x.font='12px sans-serif';x.fillText(xleft.toFixed(1),m.l,H-8);x.textAlign='right';x.fillText(xright.toFixed(1),W-m.r,H-8);x.textAlign='left';return{m,px,py}}
        function curve(x,xs,ys,px,py,color,dash=[]){x.beginPath();for(let i=0;i<xs.length;i++){const X=px(xs[i]),Y=py(ys[i]);i?x.lineTo(X,Y):x.moveTo(X,Y)}x.strokeStyle=color;x.lineWidth=2.3;x.setLineDash(dash);x.stroke();x.setLineDash([])}
        function secular(){if(!result)return;const [x,W,H]=setup(by('rf-secular')),ys=result.curve_r,limit=Math.max(...ys.map(Math.abs),.05)*1.12,a=axes(x,W,H,'phase velocity c (km/s)',0,state.beta,-limit,limit);curve(x,result.curve_c,ys,a.px,a.py,'#f59e0b');for(const item of [[result.c_r,'#22d3ee'],[state.beta,'#a78bfa']]){x.strokeStyle=item[1];x.lineWidth=1.7;x.setLineDash([5,4]);x.beginPath();x.moveTo(a.px(item[0]),a.m.t);x.lineTo(a.px(item[0]),H-a.m.b);x.stroke();x.setLineDash([])}}
        function orbit(){if(!result)return;const [x,W,H]=setup(by('rf-orbit')),R=1.25,cx=W/2,cy=H/2,S=Math.min(W,H)*.32,px=v=>cx+v/R*S,py=v=>cy-v/R*S;x.strokeStyle='#374151';x.lineWidth=1;x.beginPath();x.moveTo(cx-S,cy);x.lineTo(cx+S,cy);x.moveTo(cx,cy-S);x.lineTo(cx,cy+S);x.stroke();curve(x,result.orbit_x,result.orbit_z,px,py,'#60a5fa');x.fillStyle='#fbbf24';x.beginPath();x.arc(px(result.snapshot_x),py(result.snapshot_z),5,0,2*Math.PI);x.fill();x.fillStyle='#9ca3af';x.font='12px sans-serif';x.fillText('horizontal',W-74,cy-7);x.fillText('vertical',cx+6,24);x.fillText('normalised displacement',12,H-10)}
        function drawArrow(x,sx,sy,dx,dy,color){const len=Math.hypot(dx,dy);if(len<0.6)return;const ux=dx/len,uy=dy/len,ah=4,aw=2.4;x.strokeStyle=color;x.lineWidth=1.6;x.beginPath();x.moveTo(sx,sy);x.lineTo(sx+dx,sy+dy);x.stroke();const tx=sx+dx,ty=sy+dy,bx=tx-ux*ah,by2=ty-uy*ah;x.beginPath();x.moveTo(tx,ty);x.lineTo(bx-uy*aw,by2+ux*aw);x.lineTo(bx+uy*aw,by2-ux*aw);x.closePath();x.fillStyle=color;x.fill()}
        function idxWF(fr,ix,iz,nx,nz){return (fr*nx+ix)*nz+iz}
        function wavefield(){
          if(!result) return;
          const [x,W,H]=setup(by('rf-wave'));
          const m={l:46,r:14,t:10,b:26};
          const xmax=result.wf_xs[result.wf_xs.length-1], zmax=result.wf_zs[result.wf_zs.length-1];
          const px=v=>m.l+v/xmax*(W-m.l-m.r), py=v=>m.t+v/zmax*(H-m.t-m.b);
          const nx=result.wf_nx, nz=result.wf_nz, fr=frameIdx;
          const arrowScale=26;
          for(const ix of result.wf_xs.keys()){
            for(const iz of result.wf_zs.keys()){
              let ux=0, uz=0;
              if(showP){ const b=idxWF(fr,ix,iz,nx,nz); ux+=result.wf_uxP[b]; uz+=result.wf_uzP[b]; }
              if(showS){ const b=idxWF(fr,ix,iz,nx,nz); ux+=result.wf_uxS[b]; uz+=result.wf_uzS[b]; }
              const sx=px(result.wf_xs[ix]), sy=py(result.wf_zs[iz]);
              drawArrow(x, sx, sy, ux*arrowScale, uz*arrowScale, '#60a5fa');
            }
          }
          x.strokeStyle='#ef4444'; x.lineWidth=2; x.beginPath(); x.moveTo(m.l,py(0)); x.lineTo(W-m.r,py(0)); x.stroke();
          x.fillStyle='#9ca3af'; x.font='12px sans-serif'; x.fillText('distance (km)', m.l, H-8);
          x.save(); x.translate(14,H/2); x.rotate(-Math.PI/2); x.textAlign='center'; x.fillText('depth (km)',0,0); x.restore();
        }
        function draw(){secular();orbit();wavefield();if(result)by('rf-readout').innerHTML='c<sub>R</sub> = <b>'+result.c_r.toFixed(3)+' km/s</b> = '+(100*result.c_r/state.beta).toFixed(1)+'% of β<br>Amplitude at '+state.depth.toFixed(1)+' km: <b>'+(100*result.surface_ratio).toFixed(1)+'%</b> of surface';}
        function tick(now){
          if (window.__rfCurrentInstance !== rfInstance) return;
          if(playing && result && now-lastFrameTime>85){ frameIdx=(frameIdx+1)%result.wf_nframes; lastFrameTime=now; wavefield(); }
          rafId=requestAnimationFrame(tick);
        }
        function onControl(event){const key=ids[event.target.id];if(!key)return;event.stopImmediatePropagation();state[key]=Number(event.target.value);labels();emit()}
        par.addEventListener('input',onControl,true);
        par.addEventListener('change',onControl,true);
        by('rf-reset').onclick=()=>{Object.assign(state,{alpha:6,beta:4.49,density:2.9,omega:.15,time:.25,depth:0});for(const [id,key]of Object.entries(ids))by(id).value=state[key];labels();emit()};
        by('rf-show-p').addEventListener('change',e=>{showP=e.target.checked;wavefield()});
        by('rf-show-s').addEventListener('change',e=>{showS=e.target.checked;wavefield()});
        by('rf-play').addEventListener('click',()=>{playing=!playing;by('rf-play').textContent=playing?'Pause':'Play';if(playing){lastFrameTime=0}});
        window.addEventListener('rayleigh-function-results',event=>{result=event.detail?JSON.parse(event.detail):null;draw()},rfSignal);
        new ResizeObserver(draw).observe(par);
        rafId=requestAnimationFrame(tick);
        labels(); draw();
        })();
        </script></div>
        """)
    end

    const _rf_ready = true
end

# ╔═╡ 544d41b1-b6ce-4844-8269-a0251b391962
begin
    # The widget implementation lives in the Appendix. This explicit reference
    # makes a cold Pluto load wait for its definition.
    _rf_ready
    PlutoUI.WideCell(@bind _rf RayleighFunctionInput(); max_width=1400)
end

# ╔═╡ 31a317e2-0864-4552-83ee-ba4138ab5c0b
begin
    # Convert the browser value at one boundary; the remainder of the notebook
    # only sees concrete, physically admissible numerical values.
    rf_safe = _rf isa AbstractDict ? _rf : Dict{String,Any}()
    αp = clamp(Float64(get(rf_safe, "alpha", 6.0)), 4.0, 10.0)
    βp = clamp(Float64(get(rf_safe, "beta", 4.49)), 2.0, min(6.0, 0.98 * αp))
    ρp = clamp(Float64(get(rf_safe, "density", 2.9)), 2.0, 4.0)
    ωp = clamp(Float64(get(rf_safe, "omega", 0.15)), 0.05, 0.35)
    rf_time = clamp(Float64(get(rf_safe, "time", 0.25)), 0.0, 1.0)
    rf_depth = clamp(Float64(get(rf_safe, "depth", 0.0)), 0.0, 25.0)
end

# ╔═╡ bd7562f0-cb14-4ddb-b8c7-ffc2a37b079e
let
    payload = rayleigh_widget_payload(αp, βp, ωp, rf_time, rf_depth)
    wf = rayleigh_wavefield_frames(αp, βp, ωp)
    number(x) = isfinite(x) ? string(round(Float64(x), digits = 7)) : "0"
    array(values) = "[" * join(number.(values), ",") * "]"
    # (frame,ix,iz) -> flatten with iz fastest, then ix, then frame slowest
    flat(A) = vec(permutedims(A, (3, 2, 1)))
    message = string(
        "{\"curve_c\":", array(payload.curve_c), ",\"curve_r\":", array(payload.curve_r),
        ",\"orbit_x\":", array(payload.orbit_x), ",\"orbit_z\":", array(payload.orbit_z),
        ",\"snapshot_x\":", number(payload.snapshot_x), ",\"snapshot_z\":", number(payload.snapshot_z),
        ",\"c_r\":", number(payload.c_r), ",\"surface_ratio\":", number(payload.surface_ratio),
        ",\"wf_xs\":", array(wf.xs), ",\"wf_zs\":", array(wf.zs),
        ",\"wf_nx\":", wf.nx, ",\"wf_nz\":", wf.nz, ",\"wf_nframes\":", wf.nframes,
        ",\"wf_uxP\":", array(flat(wf.uxP)), ",\"wf_uzP\":", array(flat(wf.uzP)),
        ",\"wf_uxS\":", array(flat(wf.uxS)), ",\"wf_uzS\":", array(flat(wf.uzS)), "}",
    )
    HTML("""<script>
      window.dispatchEvent(new CustomEvent('rayleigh-function-results', {detail: $(repr(message))}));
    </script>""")
end

# ╔═╡ c830a6e9-61bf-49b2-95d1-8dea1b89c154
RUI = build_function(subs_αβρ_plot(R, αp, βp, ρp), p, expression=Val{false})

# ╔═╡ a53444b0-cb7f-4729-a262-38b33766bb74
begin
    cᵣ = rayleigh_speed(αp, βp)
    md"The root is cᵣ=$(cᵣ) (km/s)"
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
PlutoTeachingTools = "~0.3.1"
PlutoUI = "~0.7.72"
SymbolicUtils = "~3.7.2"
Symbolics = "~6.22.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "a8f6a03aa7eac564a5fb447e969ae86e7d7eaa02"

[[deps.ADTypes]]
git-tree-sha1 = "27cecae79e5cc9935255f90c53bb831cc3c870d7"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.18.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

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
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

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
git-tree-sha1 = "6aaafea90a56dc1fc8cbc15e3cf26d6bc81eb0a3"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.10"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "980f01d6d3283b3dbdfd7ed89405f96b7256ad57"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "2.0.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.Combinatorics]]
git-tree-sha1 = "8010b6bb3388abe68d95743dcbea77650bb2eddf"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.3"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

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

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3bc002af51045ca3b47d2e1787d6ce02e68b943a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.122"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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
git-tree-sha1 = "ca693f8707a77a0e365d49fe4622203b72b6cf1d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.3"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

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

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "173e4d8f14230a7523ae11b9a3fa9edb3e0efd78"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.14.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

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
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
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

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "277779adfedf4a30d66b64edc75dc6bb6d52a16e"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.6"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

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

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "e24491cb83551e44a69b9106c50666dea9d953ab"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.4.4"

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

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

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
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "fade91fe9bee7b142d332fc6ab3f0deea29f637b"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.9"

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
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoLinks", "PlutoUI"]
git-tree-sha1 = "8252b5de1f81dc103eb0293523ddf917695adea1"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.3.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "f53232a27a8c1c836d3998ae1e17d898d4df2a46"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "c05b4c6325262152483a1ecb6c69846d2e01727b"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.34"

    [deps.PreallocationTools.extensions]
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "51bdb23afaaa551f923a0e990f7c44a4451a26f1"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.39.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "276237d98b5bfc2541d92e952bfa99d362011ac7"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.11.0"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "86a8a8b783481e1ea6b9c91dd949cb32191f8ab4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.15"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "7680fbbc8a4fdf9837b4cae5e3fbebe53ec8e4ff"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.122.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseEnzymeExt = "Enzyme"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBaseMeasurementsExt = "Measurements"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReverseDiffExt = "ReverseDiff"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "c1053ba68ede9e4005fc925dd4e8723fcd96eef8"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.9.0"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

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

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

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

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "a136f98cefaf3e2924a66bd75173d1c891ab7453"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.7"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

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
deps = ["SymbolicUtils"]
git-tree-sha1 = "f75c7deb7e11eea72d2c1ea31b24070b713ba061"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.3"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "04e9157537ba51dad58336976f8d04b9ab7122f0"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.7.2"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "IfElse", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "493c70c2a94a96dff97523afd2fb4ea7fd346b93"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.22.1"

    [deps.Symbolics.extensions]
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

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
# ╠═d921cf1c-ca8b-44d9-8d00-e1ed55267647
# ╠═915243cb-9619-4eab-8ad4-1d9b821e8d01
# ╟─07ee82c2-d17b-4b85-9ef6-c107b76d0cfe
# ╟─544d41b1-b6ce-4844-8269-a0251b391962
# ╟─31a317e2-0864-4552-83ee-ba4138ab5c0b
# ╟─318ac5d5-bd8b-417b-a919-b46ee555ca1a
# ╟─bd7562f0-cb14-4ddb-b8c7-ffc2a37b079e
# ╟─c8253470-7d3c-436f-bced-141c06786ff9
# ╠═8afdf9d3-951d-4e53-8478-d5a075acbbe3
# ╠═aafebaf4-9091-42f0-8e5c-cdae92c21203
# ╠═91506b9c-af2f-44ab-a200-04ba88bc53df
# ╟─6eb3e9d6-2afa-450d-915c-81b0bfa92e9d
# ╠═e19ff895-cee0-45ad-994d-26b91e7c4d3c
# ╟─3e9930db-b54c-45af-8bb7-fc86568ae504
# ╟─9409762a-7803-4001-afbd-62c5d40dea2e
# ╠═296dbce3-8e24-49b2-a8c9-5f59c6586160
# ╟─a38858bd-d1ea-4918-8445-9b1306840979
# ╠═ca8f7333-8bb0-42f8-82d4-b99f1b03c2c9
# ╠═ce12c0d3-1601-431c-9283-a88a44b725b4
# ╠═1213a9c2-0194-4f57-8b64-5f810b46ee67
# ╟─ef389316-a822-4a31-865e-5751c8bc6fe0
# ╠═e88d9f22-d833-4d4c-b651-a7d2e318323e
# ╟─59c61d01-e0e3-405b-b048-1ceebee0d29b
# ╠═0ffac267-6380-4288-8080-5f1ae30646d8
# ╠═ff414220-762b-4e82-950c-833dd0282772
# ╟─596dccc2-7853-4562-8f25-110108ff97c2
# ╠═a61d687c-5e13-423a-98b0-9b73371a8044
# ╟─a254e7a7-d5ef-467a-98eb-f4bb1ceed97c
# ╠═019ac47c-a41f-4da2-adf5-fd1d6bac7f5d
# ╟─28291690-152c-4bca-b194-f1aac4b102c7
# ╠═ff965bbb-61bc-4595-8082-fb7a8b4640ba
# ╟─127d393b-1cf5-4b4a-bb98-b271d0167b97
# ╠═0aa629dd-0201-463e-9765-3bc002cdc252
# ╠═83366f56-03a2-4134-8d94-8d54943ce3ff
# ╠═6f6af471-2d5f-44ee-bb4f-ccd441dd8248
# ╠═6308d52d-809c-4f1f-b508-61388b1d6557
# ╠═a963b758-1ea5-4248-8544-88a612815d88
# ╠═60409882-81c5-4645-85ac-4c40ac92327b
# ╠═d2a84a1e-6529-467b-905a-fdcc21135245
# ╟─3daee10c-54d5-4071-b577-b4a8abd90976
# ╠═a22543fe-63cd-4cfb-8791-dfc2eb929359
# ╟─288c375f-ff18-45a8-bc12-f5b3cc4fa377
# ╠═1a0819fb-b998-4fac-91d3-b3b375cdfc42
# ╟─d8369538-ab73-4575-bf86-31b9111523cc
# ╠═c830a6e9-61bf-49b2-95d1-8dea1b89c154
# ╟─a7c6f76d-fbcd-4fa6-a12f-868428cedef8
# ╠═a53444b0-cb7f-4729-a262-38b33766bb74
# ╟─8ad847c8-3d83-4c15-97f1-9b4490cbf2b0
# ╠═8d54de78-38f7-4fd0-87e8-cb3ed6009b87
# ╠═6bdcc2ac-3319-432d-a286-f706bf68234e
# ╟─379feb55-32e5-4867-b108-b10c2c0782c3
# ╠═5b807994-416d-11ed-22ff-77b37ff3cfac
# ╟─aa111199-0001-4001-8001-100000000001
# ╠═996bd1d1-09f9-4774-b73f-0f4e58f52455
# ╠═18fa43a1-4e0d-4910-9125-d41d6d174ef1
# ╠═3053571a-096c-448a-9b48-60f9a1b2fce6
# ╟─aa111199-0002-4002-8002-100000000002
# ╠═bb111111-1111-4111-8111-111111111111
# ╟─aa111199-0003-4003-8003-100000000003
# ╠═bb222222-2222-4222-8222-222222222222
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
