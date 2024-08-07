### A Pluto.jl notebook ###
# v0.19.36

#> [frontmatter]
#> title = "Rayleigh Function"
#> description = "Particles may go round and round, but don't worry, let the Rayleigh wave pass by! This notebook will help you understand the physics behind it."

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5b807994-416d-11ed-22ff-77b37ff3cfac
begin
    using Symbolics
    using SymbolicUtils
    using LinearAlgebra
    using Latexify
    using LaTeXStrings
    using PlutoUI
    using PlutoTeachingTools
	using PlutoPlotly
    using Measures
    using Roots
	using Parameters
end

# ╔═╡ d921cf1c-ca8b-44d9-8d00-e1ed55267647
ChooseDisplayMode()

# ╔═╡ 915243cb-9619-4eab-8ad4-1d9b821e8d01
PlutoUI.TableOfContents()

# ╔═╡ 07ee82c2-d17b-4b85-9ef6-c107b76d0cfe
md"""
# Rayleigh Function
This notebook shall help us thoroughly grasp the theory for Rayleigh waves. While discussing reflection coefficients, we have seen examples in which the incident wave is homogeneous with at least one of the scattered waves being homogeneous. Here, all the plane waves interacting with the boundary are inhomogeneous. Specifically, we consider inhomogeneous P- and SV-waves interacting with the free surface of the medium.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 544d41b1-b6ce-4844-8269-a0251b391962
TwoColumnWideLeft(md"""
Choose the medium parameters of the half-space using these sliders. Accordingly, we 
shall plot the Rayleigh function, where its root gives the phase velocity of the Rayleigh wave $cᵣ$. It can be noticed that this velocity is a few percent less than the shear-wave speed.
""",
md"""
α (km/sec) $(@bind αp PlutoUI.Slider(range(4.0, stop=10.0, step=0.1), default=6.0, show_value=true))

β (km/sec) $(@bind βp PlutoUI.Slider(range(3.0, stop=6.0, step=0.1), default=4.49, show_value=true))

ρ (gm/cm³) $(@bind ρp PlutoUI.Slider(range(2.0, stop=6.0, step=0.1), default=2.9, show_value=true))
""")

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
$\alpha$ and $\beta$. We use $p$ to denote the horizontal slowness. As both the P- and SV-waves are inhomogeneous, a necessary condition is
```math
\frac{1}{\alpha} < \frac{1}{\beta} < p.
```
The Lame parameters are denoted using $\lambda$ and $\mu$, and $\rho$ is mass density.
"""

# ╔═╡ e19ff895-cee0-45ad-994d-26b91e7c4d3c
@syms α::Real p::Real β::Real ρ::Real λ::Real μ::Real

# ╔═╡ 3e9930db-b54c-45af-8bb7-fc86568ae504
tip(md"The inhomegeneous P- and SV-waves are coupled due to the fact that they share the same horizonal slowness $p$. If they are not coupled, one can show that the stress-free conditions cannot be satisfied.")

# ╔═╡ 9409762a-7803-4001-afbd-62c5d40dea2e
md"""
The vertical components of the slowness vectors (corresponding to $\alpha$ and $\beta$) are imaginary. Here, we use $\eta$ to denote the imaginary part of the vertical components of the slowness vector.
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
Let's now write an expression for the inhomogeneous plane wave, with amplitude $A$ and frequency $ω$, traveling along the $x$ direction with slowness $p$. We consider wave types that exponentially decay with distance from the surface $z=0$.
"""

# ╔═╡ e88d9f22-d833-4d4c-b651-a7d2e318323e
inhomo_plane(p, η, A) = A * exp(-ω * η * z) * exp(ı * ω * (t - (p * x)))

# ╔═╡ 59c61d01-e0e3-405b-b048-1ceebee0d29b
md"We now consider two harmonic plane waves that satisfy the scalar wave equations for P and SV potentials $\phi$ and $\psi$. $A_p$ and $A_s$ denote the amplitudes of these potentials."

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

# ╔═╡ 318ac5d5-bd8b-417b-a919-b46ee555ca1a
TwoColumnWideLeft(md"""
For this homogeneous medium, $cᵣ$ is independent of the frequency $ω$ i.e., the wave does not disperse. In order to visualize the particle motion, let's first choose the frequency and adjust the time to plot the snapshot below. It can be noticed that 
- the particles move in an elliptical trajectory that is retrograde at the free surface;
- depending on the frequency, the magnitude of amplitude exponentially decreases with depth;
- below a certain depth, depending on the frequency, the particle motion is dominated by the SV component and becomes prograde.
""", md"""
angular frequency: $(@bind ωp PlutoUI.Slider(range(0.05, stop=0.3, length=100), show_value=true)) $(@bind tp Clock(0.5))
(experimental) plot individual components $(@bind uplot MultiCheckBox([uₚ=>"P", uₛ=>"SV"], default=[uₚ, uₛ]))
""")

# ╔═╡ 28291690-152c-4bca-b194-f1aac4b102c7
md"The total displacement field is given by summing the P and SV components. We denote the x and z components of the total displacement field using `ux` and `uz`."

# ╔═╡ ff965bbb-61bc-4595-8082-fb7a8b4640ba
ux, uz = uₚ .+ uₛ

# ╔═╡ 127d393b-1cf5-4b4a-bb98-b271d0167b97
md"""
## Stress-free Condition
We shall now evaluate the stress components `σxz` and `σzz` on the free surface and apply the boundary condition to estimate the ratio `Aₛ/Aₚ`.
"""

# ╔═╡ 0aa629dd-0201-463e-9765-3bc002cdc252
σxz = expand_derivatives(μ * (Dx(uz) + Dz(ux))) |> simplify

# ╔═╡ 83366f56-03a2-4134-8d94-8d54943ce3ff
σxz_z0 = substitute(σxz, [z => 0, ı * ı => -1]) |> simplify

# ╔═╡ 6f6af471-2d5f-44ee-bb4f-ccd441dd8248
σxz_z0 ~ 0

# ╔═╡ 6308d52d-809c-4f1f-b508-61388b1d6557
Aratio1 = Symbolics.solve_for(σxz_z0 ~ 0, Aₛ) / Aₚ |> simplify

# ╔═╡ a963b758-1ea5-4248-8544-88a612815d88
σzz = expand_derivatives(λ * (Dx(ux) + Dz(uz)) + 2 * μ * (Dz(uz))) |> simplify

# ╔═╡ 60409882-81c5-4645-85ac-4c40ac92327b
σzz_z0 = simplify(substitute(σzz, [z => 0, ı * ı => -1]))

# ╔═╡ d2a84a1e-6529-467b-905a-fdcc21135245
Aratio2 = Symbolics.solve_for(σzz_z0 ~ 0, Aₛ) / Aₚ |> simplify

# ╔═╡ 3daee10c-54d5-4071-b577-b4a8abd90976
md"""
## Rayleigh Function
In order that there is a non-trivial solution for `Aₚ` and `Aₛ`, we need the following equation to be satisfied.
"""

# ╔═╡ a22543fe-63cd-4cfb-8791-dfc2eb929359
R1 = diff(arguments(Aratio1) .* reverse(arguments(Aratio2))) |> first |> simplify

# ╔═╡ 288c375f-ff18-45a8-bc12-f5b3cc4fa377
md"""
Let's reparameterize using α, β, and ρ, and substitute ı^2=>-1 to obtain the Rayleigh
function that must be satisfied. For a given medium, the root of this function gives the horizontal slowness $p$ that both the inhomogeneous P and SV plane waves travel with.
"""

# ╔═╡ d8369538-ab73-4575-bf86-31b9111523cc
md"""
Towards plotting, we will now substitute the UI values of $\alpha$, $\beta$ and $\rho$.
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
simplify(substitute(ux, [Aₚ=>1, Aₛ=>Aratio1]))

# ╔═╡ 6bdcc2ac-3319-432d-a286-f706bf68234e
simplify(substitute(uz, [Aₚ=>1, Aₛ=>Aratio1]))

# ╔═╡ 379feb55-32e5-4867-b108-b10c2c0782c3
md"""
## Appendix
"""

# ╔═╡ 996bd1d1-09f9-4774-b73f-0f4e58f52455
# Reparametrize an expression `x` using α, β and ρ, instead of λ, μ, η.
# Also substitutes ı^2=>-1.
function subs_αβρ(x)
    simplify(substitute(x, [ı * ı => -1, λ => α * α * ρ - 2 * β * β * ρ, μ => β * β * ρ, ηₚ * ηₚ => p * p - 1 / α^2, ηₛ * ηₛ => p * p - 1 / β^2, ηₚ => sqrt(p * p - 1 / α^2), ηₛ => sqrt(p * p - 1 / β^2), μ => β * β * ρ]))
end

# ╔═╡ 1a0819fb-b998-4fac-91d3-b3b375cdfc42
R = subs_αβρ(R1)

# ╔═╡ 18fa43a1-4e0d-4910-9125-d41d6d174ef1
# Substitutes UI values of α, β, and ρ.
function subs_αβρ_plot(x)
    simplify(substitute(subs_αβρ(x), [α => αp, β => βp, ρ => ρp]))
end

# ╔═╡ c830a6e9-61bf-49b2-95d1-8dea1b89c154
RUI = build_function(subs_αβρ_plot(R), p, expression=Val{false})

# ╔═╡ a53444b0-cb7f-4729-a262-38b33766bb74
begin
    cmin = 1.0 # minimum possible speed
    cmax = βp - eps(Float64) # maximum speed, otherwise the SV wave will no longer be evanascent
    # root
    cᵣ = find_zero(c -> RUI(inv(c)), (cmin, cmax), Bisection())
    md"The root is cᵣ=$(cᵣ) (km/s)"
end

# ╔═╡ 3053571a-096c-448a-9b48-60f9a1b2fce6
function subs_all_plot(x)
    simplify(substitute(subs_αβρ_plot(x), [p => inv(cᵣ), ı => im, ω => ωp, Aₚ => 1,]))
end

# ╔═╡ 872a5a67-b217-424a-bb88-627b7b2db7a9
uxp = build_function(substitute(subs_all_plot(first(sum(uplot))), [Aₚ => 1, Aₛ => subs_all_plot(Aratio1)]), x, z, t, expression=Val{false})

# ╔═╡ b3faea59-7b99-4070-a9fa-29fc9dc9f48b
uzp = build_function(substitute(subs_all_plot(last(sum(uplot))), [Aₚ => 1, Aₛ => subs_all_plot(Aratio1)]), x, z, t, expression=Val{false})

# ╔═╡ 79d729c0-57cb-4100-812d-8907a4b3910d
md"""
### Plots
"""

# ╔═╡ 4a4dc2cb-46ce-4fed-b67b-708943753f90
function plot_displacement_field()
	xgrid = range(0, stop=500, length=10)
	zgrid = range(0, stop=250, length=10)
	# evaluate ux and uz on the grids
    ux = [real(uxp(x, z, tp)) for x in xgrid, z in zgrid]
    uz = [real(uzp(x, z, tp)) for x in xgrid, z in zgrid]

	   
    # normalize ux and uz
    ux ./= abs(uxp(0, 0, 0))
    uz ./= abs(uxp(0, 0, 0))
    strength = vec(sqrt.(ux .^ 2 .+ uz .^ 2))

	xq = [x1 for x1 in xgrid, z1 in zgrid]
	yq = [z1 for x1 in xgrid, z1 in zgrid]


scale=25
	
vels = [attr(
        x = xq[i] + scale * ux[i], y = yq[i] + scale * uz[i],
        showarrow = true, 
        axref = "x",
        ayref = "y",
	text="∘",
		arrowcolor="blue",
        arrowsize=strength[i],
        arrowwidth=2,
        arrowhead=3,
        ax =  xq[i],
        ay = yq[i])
for i in 1:length(xq)]

layout = Layout(; uirevision=1, template="white", yaxis=attr(title="Depth", range=(400, -100)),xaxis=attr(title="Distance", range=(-100, 600)), title_text="Rayleigh Displacement Field",
            autosize = false, width = 700, height = 500, 
            margin = attr(l = 50, r = 50, b = 50, t = 65),
            showlegend = false,
            annotations = vels
            )

	fig = Plot(layout)
	
add_trace!(fig, scatter(y=[0,0], x=extrema(xgrid),
	mode="lines+markers+text",
    text=["", "Free Surface"],
	textfont=attr(
        family="sans serif",
        size=15,
        color="Red"
	),
	textposition="top center",
    line=attr(
        color="Red",
        width=2,
    )
))
	
    return plot(fig)
end

# ╔═╡ bd7562f0-cb14-4ddb-b8c7-ffc2a37b079e
plot_displacement_field()

# ╔═╡ 38f7ae72-edae-4ebf-9600-1104e4fd92cf
function plot_Rayleigh_function()
    fig = Plot(Layout(title="Rayleigh Function", xaxis_title="Phase Velocity (km/s)",width=700, height=350))
    # declare a range of speeds for plotting
    crange = range(cmin, stop=cmax, length=100)
	add_trace!(fig, PlutoPlotly.scatter(x=crange, y=[RUI(inv(c)) for c in crange], line=attr(width=3, color=:black), name="Rayleigh Function"))
	add_hline!(fig, 0.0, line=attr(dash="dot"))
	add_trace!(fig, PlutoPlotly.scatter(x=[cᵣ, cᵣ], y=[-1, 1],
		name="Rayleigh Phase Velocity (cᵣ)",
    line=attr(
        color="Orange",
        width=4,
    )
))
		add_trace!(fig, PlutoPlotly.scatter(x=[αp, αp], y=[-1, 1],
			name="P-wave Velocity (α)",
    line=attr(
        color="Red",
        width=4,
    )
))
		add_trace!(fig, PlutoPlotly.scatter(x=[βp, βp], y=[-1, 1],
			name="S-wave Velocity (β)",
    line=attr(
        color="Blue",
        width=4,
    )
))
	relayout!(fig, template = :plotly_dark)
	PlutoPlotly.plot(fig)
end

# ╔═╡ 31a317e2-0864-4552-83ee-ba4138ab5c0b
plot_Rayleigh_function()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Measures = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
LaTeXStrings = "~1.3.0"
Latexify = "~0.16.1"
Measures = "~0.3.2"
Parameters = "~0.12.3"
PlutoPlotly = "~0.4.1"
PlutoTeachingTools = "~0.2.13"
PlutoUI = "~0.7.52"
Roots = "~2.0.20"
SymbolicUtils = "~1.4.0"
Symbolics = "~5.8.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "498af720fba0accf3427fc25524b459899bece99"

[[deps.ADTypes]]
git-tree-sha1 = "5d2e21d7b0d8c22f67483ef95ebdc39c0e6b6003"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.4"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "c3c29bf6363b3ac3e421dc8b2ba8e33bdacbd245"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.32.5"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "1c9b6f39f40dba0ef22244a175e2d4e42c8f6ee7"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.0"

[[deps.Bijections]]
git-tree-sha1 = "71281c0c28f97e0adeed24fdaa6bf7d37177f297"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.5"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "c0216e792f518b39b22212127d4a84dc31e4e386"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.5"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "3d5873f811f582873bb9871fc9c451784d5dc8c7"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.102"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "51b4b84d33ec5e0955b55ff4b748b99ce2c3faa9"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.7"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "fea68c84ba262b121754539e6ea0546146515d4f"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.3"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "a20eaa3ad64254c61eeb5f230d9306e937405434"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.6.1"
weakdeps = ["SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

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

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "ExprTools", "Logging", "MultivariatePolynomials", "Primes", "Random", "SIMD", "SnoopPrecompile"]
git-tree-sha1 = "44f595de4f6485ab5ba71fe257b5eadaa3cf161e"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.4.4"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "8e59ea773deee525c99a8018409f64f19fb719e6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.7"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "81dc6aefcbe7421bd62cb6ca0e700779330acff8"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.25"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "cd04158424635efd05ff38d5f55843397b7416a9"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.14.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

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

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "6c2e970692b6f4fed2508865c43a0f67f3820cf4"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.2"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "6985021d02ab8c509c841bb8b2becd3145a7b490"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.3.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "fcf8fd477bd7f33cb8dbb1243653fb0d415c256c"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.25"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

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

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "BaseDirs", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "Reexport", "TOML"]
git-tree-sha1 = "9fefc3bfea24f08474e86e86743ee7f8f1bf12a0"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.4.1"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "542de5acb35585afcf202a6d3361b430bc1c3fbd"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.13"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "4c9f306e5d6603ae203c2000dd460d81a5251489"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9ebcd48c498668c7fa0e97a9cae873fbee7bfee1"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "d7087c013e8a496ff396bae843b1e16d9a30ede8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.10"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "609c26951d80551620241c3d7090c71a73da75ab"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.6"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "06b5ac80ff1b88bd82df92c1c1875eea3954cd6e"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.20"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "0e270732477b9e551d884e6b07e23bb2ec947790"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.5"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "ChainRulesCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "cc82f4c0ee2493cb2c1b3df66a2ef7fc45b6eacb"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.2.0"

    [deps.SciMLBase.extensions]
    PyCallExt = "PyCall"
    ZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "65c2e6ced6f62ea796af251eb292a0e131a3613b"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.6"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "0adf069a2a490c47273727e029371b31d44b72b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.5"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "2f3fa844bcd33e40d8c29de5ee8dded7a0a70422"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.4.0"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "4e703925f8d0954974f466ee08e866b1513a6d08"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.8.0"

    [deps.Symbolics.extensions]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a1f34829d5ac0ef499f6d84428bd6b4c71f02ead"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.0"

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

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "21c8fc7cd598ef49f11bc9e94871f5d7740e34b9"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "977aed5d006b840e2e40c0b48984f7463109046d"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
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
# ╠═996bd1d1-09f9-4774-b73f-0f4e58f52455
# ╠═18fa43a1-4e0d-4910-9125-d41d6d174ef1
# ╠═3053571a-096c-448a-9b48-60f9a1b2fce6
# ╠═872a5a67-b217-424a-bb88-627b7b2db7a9
# ╠═b3faea59-7b99-4070-a9fa-29fc9dc9f48b
# ╟─79d729c0-57cb-4100-812d-8907a4b3910d
# ╠═4a4dc2cb-46ce-4fed-b67b-708943753f90
# ╠═38f7ae72-edae-4ebf-9600-1104e4fd92cf
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
