### A Pluto.jl notebook ###
# v0.19.30

#> [frontmatter]
#> title = "Seismic Source Theory"
#> description = "Representation of the displacement field due to slip across fault planes in the source region."

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

# ╔═╡ 1df6b4f0-4485-11ee-1284-4967e4928ca8
begin
    using PlutoUI
    using PlutoTeachingTools
    using Symbolics
    using SymbolicUtils
    using LinearAlgebra
    using Einsum
    using PlutoHooks
    using PlutoPlotly
    using TikzPictures
end

# ╔═╡ 5f6e65e0-48c8-46d8-bc2f-460b93eb64fe
ChooseDisplayMode()

# ╔═╡ 3376d2e5-9eb7-4d5d-8682-3961550ed9b8
TableOfContents()

# ╔═╡ e646e728-823a-4d9a-b51d-36e544e18426
md"""
# Seismic Source Theory
This notebook provides a simple demonstration of the representation theorem and its use in seismic source theory. For internal sources, it derives expressions for _seismic moment-density tensor_ in the case of an isotropic elastic media, where a displacement discontinuity is assumed across two adjacent surfaces
internal to a volume.


| Fault Orientation | Slip | Measured Displacement|
|:---------- | :------------:|:------------:|
| $(@bind fault_input Radio(["1", "2", "3"], default="1")) |   $(@bind slip_input Radio(["1", "2", "3"], default="2"))         | $(@bind measurement_input Radio(["1", "2", "3"], default="1"))   
"""

# ╔═╡ 9aaa0c16-9322-4873-bf3b-363cc2261381
md"The displacement field in the second state is given by integrating the moment tensor density over the source region, which constitutes of geologic faults or underground volume explosions."

# ╔═╡ e2a789a9-0c20-4475-a4c7-4e6de420167f
md"""
## Setting Up Two States
"""

# ╔═╡ c504c3ba-f37c-408d-92c8-3e12ce9f7e97
@variables x[1:3]

# ╔═╡ 410e2695-7744-4b3e-8ff2-6601bdcd33dd
md"### Displacement Fields
Declaring symbols for vector displacement fields in the temporal frequency domain.
"

# ╔═╡ 47d69634-9452-48b3-b2d8-8ad629218381
@variables u[1:3] v[1:3]

# ╔═╡ 074f62ae-3886-4582-9afe-fb7df0f5d46c
md"### Body Force Densities
Body-force densities in both the states. For internal seismic sources, we assume that the body forces are absent in the case of the second state."

# ╔═╡ b07db20f-d693-4aae-8a9a-f26a7cccacd1
@variables f[1:3] h[1:3]

# ╔═╡ 9ec615ae-0490-44d1-99b8-519ac1e5a7da
md"### c_ijkl
The elastic constants for an isotropic medium are $\lambda$ and $\mu$. These constants are identical for both the states, and their values at the source region affect the resultant displacement field."

# ╔═╡ e8bca0a0-4807-455a-8f7d-357b250724fb
@variables λ μ

# ╔═╡ e5ae04e6-f9af-447b-90ae-0a594addfb9d
Ciso = [[λ + 2μ, λ, λ, 0, 0, 0];; [λ, λ + 2μ, λ, 0, 0, 0];; [λ, λ, λ + 2μ, 0, 0, 0];; [0, 0, 0, μ, 0, 0];; [0, 0, 0, 0, μ, 0];; [0, 0, 0, 0, 0, μ]]

# ╔═╡ 352fccff-dd16-4cbe-8b41-d93654ea99c2
function get_cijkl(C)
    [C[i*(isequal(i, j))+(1-isequal(i, j))*(9-i-j), k*isequal(k, l)+(1-isequal(k, l))*(9-k-l)] for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
end

# ╔═╡ fb5fd18b-f21d-42d3-abbd-86a66549469d
ciso = get_cijkl(Ciso);

# ╔═╡ e515a71a-82d4-41fd-a271-21066821285a
md"### Strain Tensors"

# ╔═╡ 07c8663f-0383-4a2f-bf61-b1557bd1ab52
@einsum eᵤ[i, j] := 0.5 * (∇[i](u[j]) + ∇[j](u[i]))

# ╔═╡ 222cf38a-3357-44a2-9b53-d4d783b14d52
@einsum eᵥ[i, j] := 0.5 * (∇[i](v[j]) + ∇[j](v[i]))

# ╔═╡ 9f723fc0-d8d1-4188-a4b1-3808432ee8e2
md"### Stress Tensors"

# ╔═╡ 45b7bb47-e898-4bc6-8238-237804df9933
@einsum σᵤ[i, j] := ciso[i, j, k, l] * eᵤ[k, l]

# ╔═╡ 0a1f32ec-3294-4c1b-ba26-b97edaa32926
@einsum σᵥ[i, j] := ciso[i, j, k, l] * eᵥ[k, l]

# ╔═╡ cf624a40-d9ea-4fe5-bd86-f90802abb6d9
md"""
## Betti's Theorem
"""

# ╔═╡ 7edf2698-b97d-4685-887b-26413fc849b3
@syms ∭(f) ∬(f) # volume and surface integrals

# ╔═╡ 691683af-5e51-43a8-9d22-a57f2eb739fd
Bv1 = dot(u, h)

# ╔═╡ 26b671cd-fbc6-4e9e-895b-dc3a3bcfcdd9
Bv2 = dot(v, f)

# ╔═╡ 68581620-11a0-458b-8288-1ffea41b5d6c
md"""
## Elastodynamic Green's Tensor
We shall denote the spatial coordinate in the source region using $\xi$. The ith component of the displacement field in the frequency domain, due to an impulsive force acting in the jth direction, is given by  
$G_{ji}(\xi, x)$.
"""

# ╔═╡ b4a4d08f-e75c-4a1e-b222-2eba54c18347
@variables ξ[1:3]

# ╔═╡ 9c98532d-650f-48cb-8432-940a0bced104
@syms G₁₁(x₁, x₂, x₃, ξ₁, ξ₂, ξ₃) G₁₂(x₁, x₂, x₃, ξ₁, ξ₂, ξ₃) G₁₃(x₁, x₂, x₃, ξ₁, ξ₂, ξ₃) G₂₁(x₁, x₂, x₃, ξ₁, ξ₂, ξ₃) G₂₂(x₁, x₂, x₃, ξ₁, ξ₂, ξ₃) G₂₃(x₁, x₂, x₃, ξ₁, ξ₂, ξ₃) G₃₁(x₁, x₂, x₃, ξ₁, ξ₂, ξ₃) G₃₂(x₁, x₂, x₃, ξ₁, ξ₂, ξ₃) G₃₃(x₁, x₂, x₃, ξ₁, ξ₂, ξ₃)

# ╔═╡ 1b41f7ad-ae5d-4859-a2bb-ebbea96d2580
md"These symbols look a bit clumsy, anyway, let us form the Green's tensor"

# ╔═╡ 8ecd011b-222b-4c8d-9a30-1eda078d35f2
G = [G₁₁(x..., ξ...) G₁₂(x..., ξ...) G₁₃(x..., ξ...); G₂₁(x..., ξ...) G₂₂(x..., ξ...) G₂₃(x..., ξ...); G₃₁(x..., ξ...) G₃₂(x..., ξ...) G₃₃(x..., ξ...)]

# ╔═╡ 2aec2d90-b0b2-4f68-a11d-ebdfa8d906ff
md"## Spatial Reciprocity
Assuming homogeneous boundary conditions everywhere on $S+\Sigma^++\Sigma^-$, we can show that"

# ╔═╡ a22f2772-77be-4b50-a67f-e64a9cc06e52
G1 = [G₁₁(ξ..., x...) G₂₁(ξ..., x...) G₃₁(ξ..., x...); G₁₂(ξ..., x...) G₂₂(ξ..., x...) G₃₂(ξ..., x...); G₁₃(ξ..., x...) G₂₃(ξ..., x...) G₃₃(ξ..., x...)]

# ╔═╡ 39841a19-6f0b-4b83-8c1f-c07bb5e29239
G1 ~ G

# ╔═╡ 06537e8b-d975-4f18-a5bf-3e2a0f3ae569
md"""
## Displacement Discontinuity
"""

# ╔═╡ 99e4c5c0-d957-44eb-a8a1-5dd5b6616c82
@variables v⁺[1:3] v⁻[1:3]

# ╔═╡ 3f9657b3-9b1c-43b5-89ac-30374c861c70
md"## Appendix"

# ╔═╡ 9c448d50-8528-4293-a31d-183e9fcc5e15
tikz = @ingredients("tikz.jl")

# ╔═╡ c658a4ce-4b0c-404d-87ec-3c7861395e11
∇ = [Differential(ξ[1]), Differential(x[2]), Differential(ξ[3])]

# ╔═╡ d6ce17ee-455a-4c1f-acff-4eb8f03c938b
md"### UI"

# ╔═╡ ceb932f5-b378-41cc-8b1e-71b270294224
fault_normal = map(x -> isequal(x, fault_input), ["1", "2", "3"])

# ╔═╡ a0d9be72-dfea-4081-b0e6-c1625d22c587
Bs1 = dot(v, σᵤ * fault_normal) |> simplify

# ╔═╡ cc6c1397-862a-4b0a-a3c1-ed244530c267
Bs2 = dot(u, σᵥ * fault_normal) |> simplify

# ╔═╡ 4db54c1f-d60c-4583-ae70-f0486dbc96b3
∭(Bv1 - Bv2) ~ ∬(Bs1 - Bs2)

# ╔═╡ 3ad231b9-f397-489f-bdaa-add6ef22b58c
slip_component = map(x -> isequal(x, slip_input), ["1", "2", "3"])

# ╔═╡ 7c1fa279-9c2a-4a50-88ac-6e74a55ad2b6
measured_component = map(x -> isequal(x, measurement_input), ["1", "2", "3"])

# ╔═╡ bbf7d3fa-d752-455f-823c-dbe74085927e
iv_comp = findfirst(measured_component)

# ╔═╡ e3e9e883-c0d4-48cc-ac60-5943d87109cc
Bv1_sub = substitute(Bv1, [[u[k] => G[k, iv_comp] for k in 1:3]...]) |> expand_derivatives

# ╔═╡ baeca205-2658-4060-9f21-c7a4206a8839
Bs1_sub = substitute(Bs1, [[v[k] => slip_component[k] * (v⁺[k] - v⁻[k]) for k in 1:3]..., [u[k] => G[k, iv_comp] for k in 1:3]...]) |> expand_derivatives

# ╔═╡ c12bd922-7acc-4501-b7eb-191b8eed344b
v[iv_comp] ~ ∬(Bs1_sub)

# ╔═╡ 3784a77f-b718-4122-8cdd-eb03cb32313f
md"""
### Plots
"""

# ╔═╡ e322f25e-e798-4af4-ac79-f2ed382a0fc9
xyz = Iterators.product([:X, :Y, :Z], [:X, :Y, :Z]); # need an iterator for elements of moment tensor

# ╔═╡ 5cbbefc7-4287-46cf-9987-0666cacc14b2
begin
    struct X end
    struct Y end
    struct Z end
    struct XY end
    struct YZ end
    struct XZ end
end

# ╔═╡ 95717de5-2634-4b93-b25a-06337aa17a49
begin
    dcouple = 0.5 # half distance b/w the elements of a force couple
    body_force_locations(::X) = [[dcouple, 0.0, 0.0], [-dcouple, 0.0, 0.0]]
    body_force_locations(::Y) = [[0.0, dcouple, 0.0], [0.0, -dcouple, 0.0]]
    body_force_locations(::Z) = [[0.0, 0.0, dcouple], [0.0, 0.0, -dcouple]]
end

# ╔═╡ a6131176-8851-4504-a9ff-884543438109
# collect forces (just for plotting)
srclocs = reshape(cat([body_force_locations(eval(x[2])()) for x in xyz]..., dims=1), 1, 1, 9 * 2);

# ╔═╡ 9887499c-0ec0-4086-838b-a4c5c1542e73
function plot_moment_tensor(normal)
    ps = Plot(
        Layout(
            scene=attr(
                xaxis_range=(-5, 5),
                yaxis_range=(-5, 5),
                zaxis_range=(-5, 5),
                xaxis_title="x (m)",
                yaxis_title="y (m)",
                zaxis_title="z (m)",
                title="Equivalent Force Distribution"
            ),
            margin=attr(l=0, r=0, b=0, t=0),
            showlegend=false,
        )
    )

    # limits!(ps, x=(-5, 5), y=(-5, 5), z=(-5, 5))

    # arrows!(
    #     ps,
    #     x=vec(getindex.(srclocs, 1)),
    #     y=vec(getindex.(srclocs, 2)),
    #     z=vec(getindex.(srclocs, 3)),
    #     u=vec(getindex.(srcstrength[1, 1, :], 2)),
    #     v=vec(getindex.(srcstrength[1, 2, :], 2)),
    #     w=vec(getindex.(srcstrength[1, 3, :], 2)),
    #     line = attr(color=:black, width=0.05),
    #     sizemode = "absolute",
    #     sizeref = 0.5,
    #     color=:black,
    # )

    # Create points on the plane to form a grid
    x = range(-5, stop=5, length=50)
    y = range(-5, stop=5, length=50)
    xgrid = first(Iterators.product(x, y))
    ygrid = last(Iterators.product(x, y))
    # Calculate z values for the plane using the plane equation
    z = [[(-normal[1] .* x1 .- normal[2] .* y1) ./ normal[3] for x1 in x] for y1 in y]

    # Plot the plane
    # add_trace!(ps, surface(x=x, y=y, z=z))
    add_trace!(ps, cone(x=[1], y=[1], z=[1], u=[1], v=[1], w=[1], sizeref=1, sizemode="absolute",),)

    return plot(ps)
end

# ╔═╡ 8164e70a-c3af-4deb-9816-6fea1996403e
plot_moment_tensor([1, 1, 1])

# ╔═╡ 74617983-46c6-4558-b0c4-282a1622821c
function plot_volume(u=" ", G=false)
    L1 = L"""

    \tikzset{
      Pacman/.pic={
    \shadedraw[inner color=gray,outer color=gray!80!white,draw=black,thick] 
    (1,1) -- node[left, xshift=-0.5cm, yshift=0.25cm] (1) {$\Sigma^+$} (5:2) arc(5:360:2) node[above, xshift=-0.5cm, yshift=0.7cm] (4) {$\Sigma^-$}-- cycle;
      }
    }

    \pic[] at (1,0) {Pacman};

    \node at ($(current bounding box)+(-0.5cm,-1cm)$) {V};
    \node at ($(current bounding box)+(-2cm,1cm)$) {S};
    \node at ($(current bounding box)+(0cm,2.5cm)$) {%$u};

    """
    L2 = L"""
     \node[label={[yellow]90:$\xi$}] (xi) at (2,1) {};
     \draw [->, yellow, thick] node[midway, above left] {$x$} (2,1) to [out=150,in=0] node [midway, below, yshift=-0.5cm] {$G(x, \xi)$} (0,0);
     	"""
    if (G)
        return tikz.plot(L1 * L2)
    else
        return tikz.plot(L1)
    end
end

# ╔═╡ 789bea01-e014-43fb-8539-25ff5f4e6878
function two_state_tikz()
    TwoColumn(md"""
    $(plot_volume(L"\mathrm{State\,I}: u, \sigma_u, f", true))
    """, md"""
     $(plot_volume(L"\mathrm{State\,II}: v, \sigma_v, h", false))
     """)
end

# ╔═╡ e762644d-807a-4a31-9b07-cf59b13e0302
two_state_tikz()

# ╔═╡ Cell order:
# ╠═5f6e65e0-48c8-46d8-bc2f-460b93eb64fe
# ╠═3376d2e5-9eb7-4d5d-8682-3961550ed9b8
# ╟─e646e728-823a-4d9a-b51d-36e544e18426
# ╟─9aaa0c16-9322-4873-bf3b-363cc2261381
# ╟─c12bd922-7acc-4501-b7eb-191b8eed344b
# ╟─e762644d-807a-4a31-9b07-cf59b13e0302
# ╟─e2a789a9-0c20-4475-a4c7-4e6de420167f
# ╠═c504c3ba-f37c-408d-92c8-3e12ce9f7e97
# ╟─410e2695-7744-4b3e-8ff2-6601bdcd33dd
# ╠═47d69634-9452-48b3-b2d8-8ad629218381
# ╟─074f62ae-3886-4582-9afe-fb7df0f5d46c
# ╠═b07db20f-d693-4aae-8a9a-f26a7cccacd1
# ╟─9ec615ae-0490-44d1-99b8-519ac1e5a7da
# ╠═e8bca0a0-4807-455a-8f7d-357b250724fb
# ╠═e5ae04e6-f9af-447b-90ae-0a594addfb9d
# ╠═fb5fd18b-f21d-42d3-abbd-86a66549469d
# ╠═352fccff-dd16-4cbe-8b41-d93654ea99c2
# ╟─e515a71a-82d4-41fd-a271-21066821285a
# ╠═07c8663f-0383-4a2f-bf61-b1557bd1ab52
# ╠═222cf38a-3357-44a2-9b53-d4d783b14d52
# ╟─9f723fc0-d8d1-4188-a4b1-3808432ee8e2
# ╠═45b7bb47-e898-4bc6-8238-237804df9933
# ╠═0a1f32ec-3294-4c1b-ba26-b97edaa32926
# ╟─cf624a40-d9ea-4fe5-bd86-f90802abb6d9
# ╠═7edf2698-b97d-4685-887b-26413fc849b3
# ╠═691683af-5e51-43a8-9d22-a57f2eb739fd
# ╠═26b671cd-fbc6-4e9e-895b-dc3a3bcfcdd9
# ╠═a0d9be72-dfea-4081-b0e6-c1625d22c587
# ╠═cc6c1397-862a-4b0a-a3c1-ed244530c267
# ╠═4db54c1f-d60c-4583-ae70-f0486dbc96b3
# ╟─68581620-11a0-458b-8288-1ffea41b5d6c
# ╠═b4a4d08f-e75c-4a1e-b222-2eba54c18347
# ╠═9c98532d-650f-48cb-8432-940a0bced104
# ╟─1b41f7ad-ae5d-4859-a2bb-ebbea96d2580
# ╠═8ecd011b-222b-4c8d-9a30-1eda078d35f2
# ╟─2aec2d90-b0b2-4f68-a11d-ebdfa8d906ff
# ╠═a22f2772-77be-4b50-a67f-e64a9cc06e52
# ╠═39841a19-6f0b-4b83-8c1f-c07bb5e29239
# ╟─06537e8b-d975-4f18-a5bf-3e2a0f3ae569
# ╠═99e4c5c0-d957-44eb-a8a1-5dd5b6616c82
# ╠═bbf7d3fa-d752-455f-823c-dbe74085927e
# ╠═e3e9e883-c0d4-48cc-ac60-5943d87109cc
# ╠═baeca205-2658-4060-9f21-c7a4206a8839
# ╟─3f9657b3-9b1c-43b5-89ac-30374c861c70
# ╠═1df6b4f0-4485-11ee-1284-4967e4928ca8
# ╠═9c448d50-8528-4293-a31d-183e9fcc5e15
# ╠═c658a4ce-4b0c-404d-87ec-3c7861395e11
# ╟─d6ce17ee-455a-4c1f-acff-4eb8f03c938b
# ╠═ceb932f5-b378-41cc-8b1e-71b270294224
# ╠═3ad231b9-f397-489f-bdaa-add6ef22b58c
# ╠═7c1fa279-9c2a-4a50-88ac-6e74a55ad2b6
# ╟─3784a77f-b718-4122-8cdd-eb03cb32313f
# ╠═8164e70a-c3af-4deb-9816-6fea1996403e
# ╠═e322f25e-e798-4af4-ac79-f2ed382a0fc9
# ╠═a6131176-8851-4504-a9ff-884543438109
# ╠═5cbbefc7-4287-46cf-9987-0666cacc14b2
# ╠═95717de5-2634-4b93-b25a-06337aa17a49
# ╠═9887499c-0ec0-4086-838b-a4c5c1542e73
# ╠═74617983-46c6-4558-b0c4-282a1622821c
# ╠═789bea01-e014-43fb-8539-25ff5f4e6878
