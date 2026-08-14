### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ 41275829-5ad8-470d-b168-6d6c261ceb19
ChooseDisplayMode()

# ╔═╡ 69e5bc3b-494e-4b7d-9fea-ded035d544cc
TableOfContents()

# ╔═╡ 33a3705c-1660-4df6-bfae-23225a55bdc6
md"""# Adjoint State Method
In seismic imaging, the adjoint state method is a numerical method for efficiently computing the gradient of an objective function. 
This notebook presents a discrete version of the adjoint state formulation involving the seismic wave equation so that it can act as a reference while implementing
methods for full waveform inversion.

* We focus on the velocity-stress formulation, which is a widely used numerical method for solving wave propagation problems in geophysics. 
This method involves discretizing the governing equations of motion, such as the wave equation, using finite-difference techniques.

[Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: Pawan Bharadwaj, Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 35909f2f-3787-4c6a-9e74-f284d5eea635
md"## Medium
Created two symbolic variables ρ and μ⁻¹ which can be used to represent the density and inverse of the shear modulus of the medium as functions of position x."

# ╔═╡ d6188561-19ae-4006-b03e-9cf8e4f5e081
@variables ρ[1:3] K⁻¹[1:3]

# ╔═╡ d38130e3-611d-42b2-96b6-a6ca6210308b
md"## State Variables"

# ╔═╡ b4e9aa30-a16c-4b3a-a7be-c2fcf28f0904
md"In the velocity-stress formulation, the seismic wavefield is described in terms of two variables: the particle velocity and the stress. Two symbolic variables `v(x, t)` and `σ(x, t)` representing velocity and stress, respectively, as functions of space (`x`) and time (`t`)."

# ╔═╡ 857a5e06-2960-4ed4-826c-6bf3291387f0
@syms v(x, t) p(x, t)

# ╔═╡ 76fd526b-450a-4db1-bf9a-bc7ae59b895c
md"As the velocity and stress fields are staggered in time, we create corresponding vectors at each time value in t and t̂, respectively."

# ╔═╡ 4aec990b-338c-418b-a02d-39cec0cfc2c9
V = hcat(collect(broadcast(x -> collect(map(t -> v(x, t), vcat(t₀, t, T))), x))...) # velocity in the discrete world

# ╔═╡ c0f9c9af-aced-4066-a9d1-7af8e26c8a27
P = hcat(collect(broadcast(x -> collect(map(t -> p(x, t), vcat(t₀, t̂))), x̂))...) # stress in the discrete world

# ╔═╡ c0bb57b5-9d62-40b1-9256-648a984dbbae
md"## Forcing"

# ╔═╡ 61505a55-fc9e-4077-8365-8e166871b482
md"The symbol `f` represents the body force density, which is a quantity that describes the force per unit volume acting on a material."

# ╔═╡ 93ce7686-a297-438f-bfab-9ad1f9e9be10
@syms f(x, t)

# ╔═╡ ebefecd3-9dd2-457e-bf9b-97d4db9e983e
# F = collect(map(t -> f(0, t), t̂))

# ╔═╡ 2a2a5249-5293-4736-9a72-441f32775021
F = hcat(collect(broadcast(x -> collect(map(t -> f(x, t), vcat(t₀, t̂))), x))...)

# ╔═╡ bc20a677-5d50-4a80-b466-2b7edd5feac4
md"## Momentum Equation"

# ╔═╡ 33871626-5b9b-4718-a974-d94d7b048582
md"""
This line defines the initial momentum of the system, which is equal to the difference between the initial velocity of the particle v(x, t̂[1]) and its velocity at time t=0 (which is zero in this case) multiplied by the density of the medium ρ(x) and then subtracting the force acting on the particle at time t=0 f(0, x). It represents the balance of forces and momentum of the system at the initial time.
"""

# ╔═╡ 3cb20104-2872-4bfe-95a7-86d2b1bb6b0f
@syms δₓ⁻¹ δₜ⁻¹

# ╔═╡ 0d82f0b1-24a3-4adf-b645-4fea2fab6273
dPdx = diff(P, dims=2) * δₓ⁻¹

# ╔═╡ 84d8f063-2aad-4d57-9bdf-d44e37f173c7
md"""
We will now construct a vector that represents the time derivative of the velocity field. The first few elements can be computed using the `diff()` function applied to the `V` array. The last element is the time derivative of `v(x, t)` evaluated at the final time `T`, which is computed as the difference between `v(x, T)` and `v(x, last(t̂))`.
"""

# ╔═╡ 201ebc60-90d6-40a1-8d1f-e371057af060
dVdt = diff(V, dims=1) * δₜ⁻¹

# ╔═╡ de0ef73d-0df1-47bb-aec0-f1bdef92af83
Meqs = dVdt[:, 2:end-1] .* avx(ρt) - dPdx - F[:, 2:end-1];

# ╔═╡ 4ccfa601-b3b3-47b0-98f6-4f4ec6c1b794
ρt = repeat(reshape(ρ, 1, :), 4, 1)

# ╔═╡ a72d5528-f9a4-4b82-b4c8-644ecd0fab55
avx(ρt)

# ╔═╡ e8b1a70e-f63e-45c5-a525-2c55d99ddad4
# compute value of ρ on the velocity grid
function avx(ρt)
  hcat([0.5 * (ρt[:, i+1] .+ ρt[:, i]) for i in 1:size(ρt, 2)-1]...)
end

# ╔═╡ f54c2aa7-2afb-4806-b932-417e3b4a41e5
Meqs ~ 0

# ╔═╡ c5e41f1c-a734-4502-a60f-4ec0f9819e89
md"## Constitutive Relations
In linear elasticity, for example, the constitutive relation is given by Hooke's law, which states that stress is proportional to strain. "

# ╔═╡ 39570793-d23b-4c88-9c3b-c58690eb4ae8
∂V∂x = diff(V, dims=2) * δₓ⁻¹

# ╔═╡ fc20c71d-a23d-44be-8bf2-099a2400921b
md"Similar to velocity, we will construct a vector that represents the time derivative of the stress field."

# ╔═╡ 7343a50f-9835-4e10-97ed-5b213069044a
∂P∂t = diff(P, dims=1) * δₜ⁻¹

# ╔═╡ 0f894463-6513-49d9-98e4-12c3ca07ce51
K⁻¹t = repeat(reshape(K⁻¹, 1, :), 3, 1)

# ╔═╡ b2562463-09d0-40fc-9802-b49e105f946c
md"This equation states that the time derivative of stress minus the gradient of the stress scaled by the shear modulus is equal to zero. This is a simplified form of the constitutive relation."

# ╔═╡ a0886a56-c027-4bae-99ff-e7be53ba4a1f
Ceqs = ∂P∂t .* K⁻¹t - ∂V∂x[2:end-1, :];

# ╔═╡ 9678bb33-8e15-49ad-95b8-c8d4ab75f14f
Ceqs ~ 0

# ╔═╡ ebcbbf5e-6205-41a8-b31b-a6bed52cda5b
md"## Leap-frog Scheme"

# ╔═╡ 4d71efc9-6bf5-4aa6-8e0f-132814350351
plot_state_tikz()

# ╔═╡ 0b2ed009-e636-40d9-bf9c-8d600f3bae75
md"The leap-frog scheme updates the values of the velocity and stress components at half-integer time steps, staggered in time and space. The first two steps (stress and velocity updates) of this scheme are marked here."

# ╔═╡ 2b9eae49-a32d-4462-b5ff-8282d05a0a57
md"## Objective Function"

# ╔═╡ e301bd11-efb6-431f-a1fd-45b9a2922350
md"Variable to define an observed velocity field, which is a function of space and time."

# ╔═╡ 6332d014-79e2-451f-b65c-6a5179b7f85c
@syms v₀(x, t);

# ╔═╡ 110c5485-51c1-4837-89da-18854267ed7c
md"This defines a simple quadratic objective function that measures the difference between the model-predicted velocity field v and the observed velocity field v₀. Squaring the difference ensures that all differences are positive and also gives more weight to larger differences. The goal is to minimize this objective function to find the velocity field that best matches the observed data."

# ╔═╡ d5963693-e1c5-49f6-9c07-9187e893cd97
V₀ = hcat(collect(broadcast(x -> collect(map(t -> v₀(x, t), vcat(t))), x))...)

# ╔═╡ a8bd5aae-552f-42c6-8fc8-92ef5bade2f1
obj(v, v₀) = (v - v₀)^2

# ╔═╡ 65cc4d98-7cc8-400e-a41f-0c2cef5bde9e
md"The final value of the objective function."

# ╔═╡ e7536b06-2454-496b-8872-aec29edb37f6
J = sum(map(obj, V[2:end-1, :], V₀))

# ╔═╡ 80c1040c-88d2-45b5-89e1-b3cb8d6abaf7
md"## Adjoint State Variables
Adjoint state variables are introduced in the context of optimization problems. In this case, these variables are used to find gradients of an objective function with respect to the parameters mass density and shear modulus."

# ╔═╡ 28a345b9-d911-486c-a715-7d32a4ea41e8
@syms u(x, t) τ(x, t)

# ╔═╡ aee62f47-0bd7-492f-bc99-2392a16fc1dd
size(Meqs)

# ╔═╡ d447a8bf-2d42-4f93-a27e-9eed76050348
md"One multiplier for each of Meqs⋯"

# ╔═╡ 59d380f3-1bc0-41b4-a0e5-2cfbb0eaa333
U = hcat(collect(broadcast(x -> collect(map(t -> u(x, t), vcat(t̂, T̂))), x[2:end-1]))...)

# ╔═╡ 3b219c43-3056-49db-951e-d0e97a4f39db
size(Ceqs)

# ╔═╡ 723e81cd-36e1-44ea-99bf-b4cae2be03a2
md"One multiplier for each Ceqs⋯"

# ╔═╡ 1bf0a983-9592-4a72-a10e-af61316ce6e4
𝛕 = hcat(collect(broadcast(x -> collect(map(t -> τ(x, t), vcat([t[end-1], t[end]], T))), x̂))...)

# ╔═╡ b19344d2-d872-4bf5-a75f-1ca481779835
md"## Lagrangian"

# ╔═╡ e1048e24-e389-4e95-8157-00702f1b41a3
md"The Lagrangian component, for momentum equations `Meqs`, is given by"

# ╔═╡ fee19d3e-2235-4382-82d6-5a3cb452d317
L₁ = sum((U .* Meqs))

# ╔═╡ c2967e07-5d8b-4205-bb8b-add323494d86
md"The Lagrangian component, for constitutive equations `Ceqs`, is given by"

# ╔═╡ 55c08fb9-e29d-4500-84e7-46de7979639e
L₂ = sum(𝛕 .* Ceqs)

# ╔═╡ 911c9966-e983-46fd-8ca0-8ca06ab42bf0
md"Adding all the components of the Lagrangian together gives"

# ╔═╡ 3b8a80a7-4ab5-4b74-9aaa-3603e0ecccb5
L = L₁ + L₂ + J

# ╔═╡ 66ed7bf9-c5b7-4c27-a021-28bb5675c85f
md"## Adjoint State Equations"

# ╔═╡ 19dd77e9-4a88-4c06-9129-0c1391068900
∂vL ~ 0

# ╔═╡ 3bb85c65-89ed-4113-a498-3a0da01be0b1
∂σL ~ 0

# ╔═╡ bfd6f8e4-ea3e-4c87-b7f7-eba889e751fd
md"## Final Condition (Time Reversal)
To obtain the final condition, we shall compute the gradient of the Lagrangian with respect to `v(x, T)`."

# ╔═╡ fad3240f-9fd1-48c5-b1f3-515d9ea038bd
∂L∂VT ~ 0

# ╔═╡ 09c11c95-a82d-48e0-85ba-6db4a8c03f29
md"The solution of these adjoint equations is often obtained by using a time-reversed version of the original numerical solver, known as the backpropagation."

# ╔═╡ 57ea4b25-0f40-4816-85ba-05669770885b
plot_adjstate_tikz()

# ╔═╡ be1c590d-d70b-40f1-8370-bc294fb29c09
md"## Parameter Gradients
Lets compute the gradient of L with respect to ρ and μ."

# ╔═╡ bdc608ee-6a52-4130-b4a0-c27d513a0009
@bind iρd Select(1:3)

# ╔═╡ 7e56d621-a572-4077-83be-d3b002f4e808
# gradient w.r.t. mass density
∇ρ = Differential(ρ[iρd])(L) |> expand_derivatives#, [u(x, T̂) => 0, v(x, 0) => 0])

# ╔═╡ c150327f-b950-4a70-af44-722beee2069c
# gradient w.r.t. shear modulus
∇K⁻¹ = Differential(K⁻¹[2])(L) |> expand_derivatives #, [σ(x, 0) => 0, τ(x, T) => 0])

# ╔═╡ 3b2d5624-d365-4595-b95b-52825bc980d0
md"## Appendix"

# ╔═╡ c0d3e1c8-77d9-4f69-8f1a-97b4bec409e4
using Symbolics, SymbolicUtils, LinearAlgebra, ChainRules, PlutoUI, PlutoTeachingTools, TikzPictures

# ╔═╡ 5f6a1560-b22d-4e47-b6ee-7a306b0bc2d1
md"### Variables"

# ╔═╡ ec3bf735-5fff-44e6-8a79-a03d6eb4538f
md"These commands create three symbolic time variables `t`, and staggered time variables `t̂` and space variable `x`. We will consider three time steps for the leap-frog scheme. The `@syms` macro is used in the `Symbolics.jl` package to declare symbolic variables."

# ╔═╡ c9324a48-8cb8-45e1-909a-950333048d28
@variables t[1:3] t̂[1:3] # t̂ is staggered

# ╔═╡ 34e13768-7861-4324-8eb1-626275a3b2d1
collect(t), collect(t̂)

# ╔═╡ 3d57d1f0-4595-48a0-934c-7bc1d5dc15fc
@syms T T̂ # final time

# ╔═╡ 4876fb4e-4b00-4984-8932-a6ff576c5be6
@syms t₀

# ╔═╡ 24e8eafe-e040-4173-b705-4ade869791ba
@variables x[1:4] x̂[1:3] # space 

# ╔═╡ 6a1c66ec-6f3d-48a9-9809-2752e3818d18
md"### Lagrangian"

# ╔═╡ 69e15594-1d8e-4ddc-a1da-f968e8b2ee91
md"This code computes the gradient of the Lagrangian with respect to the  velocity field to obtain adjoint equations."

# ╔═╡ 6a2bafca-e4cb-4ce9-b322-783c8e10fbfd
∂vL = broadcast(V[2:end, 2:end-1]) do v
  Differential(v)(L) |> expand_derivatives
end;

# ╔═╡ 8d09bb0b-40d9-4d79-89a3-ad0d8679b08c
md"This code computes the gradient of the Lagrangian with respect to the stress field to obtain adjoint equations."

# ╔═╡ 2613e6c6-ec09-4574-aecf-bf2a2266bc55
∂σL = broadcast(P[2:end, :]) do σ
  Differential(σ)(L) |> expand_derivatives
end;

# ╔═╡ 269f6922-7f14-45b9-91f2-554a03a92b57
md"### Final Condition"

# ╔═╡ 7d7db129-4524-48ec-a9a9-61d7b1ee967d
∂L∂VT = broadcast(V[end, :]) do v
  Differential(v)(L) |> expand_derivatives
end;

# ╔═╡ 8d67dc24-a8b1-4a4f-9811-f46e7b09fc08
md"### Gradient Expressions"

# ╔═╡ a5d63b00-4520-4f0c-a2f2-bf6476c4a1a2
nt = 3

# ╔═╡ 653a7d32-d0cc-447f-95eb-3b8ee866b6dd
Vs1 = V[:, :]

# ╔═╡ 218136d4-bd69-4d6b-92b7-98d8f320e17f
Vs2 = reverse(Vs1, dims=1)

# ╔═╡ fb82077d-43b8-4950-96ec-ffc72d400f26
Us1 = hcat(fill(0, 4), U[:, :], fill(0, 4))

# ╔═╡ e11a771a-18b9-42b2-bb4e-c3eaadbbc3c1
Us2 = reverse(Us1, dims=1)

# ╔═╡ b3d653be-450b-49a9-bc3f-118cd5ec9921
begin
  ∇ρ1 = 0.0
  for it in 2:5
    ∇ρ1 += 0.5 * δₜ * (Us2[it-1, 2] * (Vs2[it-1, 2] - Vs2[it, 2])) + 0.5 * δₜ * (Us2[it-1, 3] * (Vs2[it-1, 3] - Vs2[it, 3]))
  end
end

# ╔═╡ 24346ac7-f089-4278-9364-d21071352500
∇ρ

# ╔═╡ f64849ef-b48f-4783-b213-c108454dcb9a
∇ρ1

# ╔═╡ e5c727d8-752c-4a3c-b971-292f4187416e
∇ρ - ∇ρ1

# ╔═╡ 4011a061-68a8-4cc5-875b-7ae27331ad4f
Ps1 = P

# ╔═╡ 4fa98b1c-0777-4f66-a029-d11488aecb2a
𝛕s1 =  𝛕[:, :]

# ╔═╡ 86321b6a-3d2d-487a-be68-5838e0ff6abf
Ps2 = reverse(Ps1, dims=1)

# ╔═╡ d52cbdaf-f13c-4ecb-8e45-537631de013b
𝛕s2  = reverse(𝛕s1, dims=1)

# ╔═╡ e3a2c0b4-5a67-4431-a9db-8116a8fefdbd
begin
  ∇K⁻¹1 = 0.0
  for it in 2:4
    global ∇K⁻¹1 += δₜ * (𝛕s2[it-1, 2] * (+Ps2[it-1, 2] - Ps2[it, 2]))
  end
end

# ╔═╡ 503a3912-31cc-45f2-b009-dca4700359f7
∇K⁻¹ - ∇K⁻¹1

# ╔═╡ ffb7b8c0-00c6-45a5-8786-8abc6918e997
∇K⁻¹

# ╔═╡ df5ba723-4ea9-4b52-a4ae-120bf37d56df
∇K⁻¹1

# ╔═╡ 75563ce3-9228-4ca5-89e1-5e704024c209
let
  ∇K⁻¹2 = 0.0
	pp = 0
	pap = 0
	p = 0
	pa = 0
  for it in 1:3
	  pp = p
	  pap = pa
	  p = Ps2[it, 2]
	  pa = 𝛕s2[it, 2]
	  @show p, pa
    # global ∇K⁻¹1 += δₜ * (𝛕s2[it-1, 2] * (+Ps2[it-1, 2] - Ps2[it, 2]))
  end
end

# ╔═╡ 38257847-5ae7-4a5a-a937-22a6729a3640
md"### Tikz"

# ╔═╡ 5a9e17d9-2552-48fd-b3ad-0a1e50279953
# t1, t2, t3, tv are texts
# s1 and s2 are labels with sizes
plot_state_tikz() = TikzPicture(L"""
   
   \tikzstyle{vertex}=[circle,minimum size=20pt,inner sep=0pt]

  \foreach \name/\x in {v(0)/1, v(t_1)/5, v(t_2)/9, v(t_3)/13, v(T)/17}
    \node[vertex,fill=black!25,] (v-\x) at (\x,0) {$\name$};

 \foreach \name/\x in {\sigma(0)/3, \sigma(\hat{t}_1)/7, \sigma(\hat{t}_2)/11, \sigma(\hat{t}_3)/15}
    \node[vertex,fill=red!25] (s-\x) at (\x,-1) {$\name$};

 \foreach \name/\x in {f(0)/3, f(\hat{t}_1)/7, f(\hat{t}_2)/11, f(\hat{t}_3)/15}
    \node[vertex,fill=blue!25] (f-\x) at (\x,1) {$\name$};

 \draw[->,thick, red] (s-3) -- (s-7) node[above, midway, red] {1};
\draw[->,thick, red] (v-5) -- (s-7) node[above, midway, red] {1};
\draw[->,thick, blue] (s-7) -- (v-9) node[below, midway, blue] {2};
	\draw[->,thick, blue] (v-5) -- (v-9) node[above, midway, blue] {2};

\node[above right=-2mm of s-3] {=0};
\node[above right=-2mm of v-1] {=0};
  """, options=tikz_default_options, preamble=tikz_preamble, width="20cm")

# ╔═╡ 21af98b7-712d-4b25-a9fa-41d008f97962
# t1, t2, t3, tv are texts
# s1 and s2 are labels with sizes
plot_adjstate_tikz() = TikzPicture(L"""


   \tikzstyle{vertex}=[circle,minimum size=20pt,inner sep=0pt]

  \foreach \name/\x in {u(0)/1, u(\hat{t}_1)/5, u(\hat{t}_2)/9, u(\hat{t}_3)/13, u(\hat{T})/17}
    \node[vertex,fill=black!25,] (v-\x) at (\x,0) {$\name$};

 \foreach \name/\x in {\tau(t_1)/3, \tau(t_2)/7, \tau(t_3)/11, \tau(T)/15}
    \node[vertex,fill=red!25] (s-\x) at (\x,-1) {$\name$};

 \foreach \name/\x in { g(t_1)/7, g(t_2)/11, g(t_3)/15}
    \node[vertex,fill=blue!25] (f-\x) at (\x,1) {$\name$};


 \draw[->,thick, red] (s-15) -- (s-11) node[above, midway, red] {1};
\draw[->,thick, red] (v-13) -- (s-11) node[above, midway, red] {1};
\draw[->,thick, blue] (s-11) -- (v-9) node[below, midway, blue] {2};
	\draw[->,thick, blue] (v-13) -- (v-9) node[above, midway, blue] {2};

\node[above right=-2mm of s-15] {=0};
\node[above right=-2mm of v-17] {=0};


  """, options=tikz_default_options, preamble=tikz_preamble, width="20cm")

# ╔═╡ 8e514e64-0172-4b47-974d-efaa8e1f4990
tikz_default_options = raw"""
  background rectangle/.style={fill=white}, show background rectangle,
  """

# ╔═╡ 00199670-40fd-4daa-979f-fb414b116bed
tikz_preamble = raw"""
  \usepackage{tikz}
  \usepackage{tikz}
  \usetikzlibrary{fit, matrix, shapes.geometric}
  \tikzset{% use tikzset, not tikzstyle
      cell/.style={
          rectangle, rounded corners=5pt, draw,
      }
  }
  \tikzset{% use tikzset, not tikzstyle
      cellv/.style={
          rectangle, rounded corners=5pt, draw, rotate=90,
      }
  }
  \usepackage{xifthen}
  \usetikzlibrary{hobby}
  \usepackage{pgfplots}
  \usepackage{fontawesome}
  \usepackage{bm,amsfonts,amsmath}
  \usetikzlibrary{backgrounds,pgfplots.groupplots,snakes}
  \usepgfplotslibrary{patchplots}
  \pgfplotsset{try min ticks=2}
  \usepackage{pgfplotstable} 
  \usetikzlibrary{plotmarks,positioning,spy}
  \usetikzlibrary{shapes.geometric, arrows, fadings}
  \usepgfplotslibrary{groupplots, polar}
  \usepackage[space]{grffile}

  \usetikzlibrary{%
              decorations.pathreplacing,%
                  decorations.pathmorphing%
                  }
                  \usetikzlibrary{positioning,fit,backgrounds}



  \usetikzlibrary{shapes,arrows}
  \usetikzlibrary{decorations.markings}
  \usetikzlibrary{patterns}
  \usetikzlibrary{plotmarks}
  \usetikzlibrary{fit}
  \usetikzlibrary{intersections}
  \usepgfplotslibrary{fillbetween}

    \pgfplotsset{
                axis line style={black!10},
                    every axis label/.append style ={black!10},
                    every axis title/.append style ={black!10},
                        every tick label/.append style={black!10}  
                          }

  % need for pgfplots
  \newcommand{\axisz}{0cm}
  \newcommand{\axisx}{0cm}

  \usetikzlibrary{positioning}
  \usetikzlibrary{shapes.geometric}
  \usetikzlibrary{backgrounds}
  """

# ╔═╡ f9f5cb83-ab6b-4de2-9ea4-4a4d984f0489
md"""## References
[^plessix]: R.-E. Plessix, A review of the adjoint-state method for computing the gradient of a functional with geophysical applications, Geophysical Journal International, Volume 167, Issue 2, November 2006, Pages 495–503, https://doi.org/10.1111/j.1365-246X.2006.02978.x
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"

[compat]
ChainRules = "~1.72.6"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.72"
SymbolicUtils = "~3.32.0"
Symbolics = "~6.56.0"
TikzPictures = "~3.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.7"
manifest_format = "2.0"
project_hash = "ea5ae4d7b53363bb09936b99272f439910e87e73"

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
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "SparseInverseSubset", "Statistics", "StructArrays", "SuiteSparse"]
git-tree-sha1 = "3b704353e517a957323bd3ac70fa7b669b5f48d4"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.72.6"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

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

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

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
deps = ["OrderedCollections"]
git-tree-sha1 = "6c72198e6a101cccdd4c9731d3985e904ba26037"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.1"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

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
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "3f50fa86c968fc1a9e006c07b6bc40ccbb1b704d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.4"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

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

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

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

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.HarfBuzz_ICU_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "6ccbc4fdf65c8197738c2d68cc55b74b19c97ac2"
uuid = "655565e8-fb53-5cb3-b0cd-aec1ca0647ea"
version = "2.8.1+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

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

[[deps.ICU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "20b6765a3016e1fca0c9c93c80d50061b94218b7"
uuid = "a51ab1cf-af8e-5615-a023-bc2c838bba6b"
version = "69.1.0+0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

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

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd10d2cc78d34c0e2a3a36420ab607b611debfbb"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.7"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

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
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

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
version = "2023.12.12"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "d38b8653b1cdfac5a7da3b819c0a8d6024f9a18c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.13"
weakdeps = ["ChainRulesCore"]

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

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
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ad31332567b189f508a3ea8957a2640b1147ab00"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

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

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "dacc8be63916b078b592806acd13bb5e5137d7e9"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "f53232a27a8c1c836d3998ae1e17d898d4df2a46"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "02148a0cb2532f22c0589ceb75c110e168fb3d1f"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "21.9.0+0"

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

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

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
version = "1.11.0"

[[deps.SparseInverseSubset]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "52962839426b75b3021296f7df242e40ecfc0852"
uuid = "dc90abb0-5640-4711-901d-7e5b23a2fada"
version = "0.1.2"

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

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

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
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "ExproniconLite", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "a85b4262a55dbd1af39bb6facf621d79ca6a322d"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.32.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "1b09f5faec5284f505c40e68ba565115e7d48718"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.56.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

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

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "tectonic_jll"]
git-tree-sha1 = "79e2d29b216ef24a0f4f905532b900dcf529aa06"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.5.0"

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

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.tectonic_jll]]
deps = ["Artifacts", "Fontconfig_jll", "FreeType2_jll", "Graphite2_jll", "HarfBuzz_ICU_jll", "HarfBuzz_jll", "ICU_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "54867b00af20c70b52a1f9c00043864d8b926a21"
uuid = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"
version = "0.13.1+0"
"""

# ╔═╡ Cell order:
# ╠═41275829-5ad8-470d-b168-6d6c261ceb19
# ╠═69e5bc3b-494e-4b7d-9fea-ded035d544cc
# ╟─33a3705c-1660-4df6-bfae-23225a55bdc6
# ╟─35909f2f-3787-4c6a-9e74-f284d5eea635
# ╠═d6188561-19ae-4006-b03e-9cf8e4f5e081
# ╟─d38130e3-611d-42b2-96b6-a6ca6210308b
# ╟─b4e9aa30-a16c-4b3a-a7be-c2fcf28f0904
# ╠═857a5e06-2960-4ed4-826c-6bf3291387f0
# ╟─76fd526b-450a-4db1-bf9a-bc7ae59b895c
# ╠═4aec990b-338c-418b-a02d-39cec0cfc2c9
# ╠═c0f9c9af-aced-4066-a9d1-7af8e26c8a27
# ╟─c0bb57b5-9d62-40b1-9256-648a984dbbae
# ╟─61505a55-fc9e-4077-8365-8e166871b482
# ╠═93ce7686-a297-438f-bfab-9ad1f9e9be10
# ╠═ebefecd3-9dd2-457e-bf9b-97d4db9e983e
# ╠═2a2a5249-5293-4736-9a72-441f32775021
# ╟─bc20a677-5d50-4a80-b466-2b7edd5feac4
# ╟─33871626-5b9b-4718-a974-d94d7b048582
# ╠═3cb20104-2872-4bfe-95a7-86d2b1bb6b0f
# ╠═0d82f0b1-24a3-4adf-b645-4fea2fab6273
# ╟─84d8f063-2aad-4d57-9bdf-d44e37f173c7
# ╠═201ebc60-90d6-40a1-8d1f-e371057af060
# ╠═de0ef73d-0df1-47bb-aec0-f1bdef92af83
# ╠═4ccfa601-b3b3-47b0-98f6-4f4ec6c1b794
# ╠═a72d5528-f9a4-4b82-b4c8-644ecd0fab55
# ╠═e8b1a70e-f63e-45c5-a525-2c55d99ddad4
# ╟─f54c2aa7-2afb-4806-b932-417e3b4a41e5
# ╟─c5e41f1c-a734-4502-a60f-4ec0f9819e89
# ╠═39570793-d23b-4c88-9c3b-c58690eb4ae8
# ╟─fc20c71d-a23d-44be-8bf2-099a2400921b
# ╠═7343a50f-9835-4e10-97ed-5b213069044a
# ╠═0f894463-6513-49d9-98e4-12c3ca07ce51
# ╟─b2562463-09d0-40fc-9802-b49e105f946c
# ╠═a0886a56-c027-4bae-99ff-e7be53ba4a1f
# ╟─9678bb33-8e15-49ad-95b8-c8d4ab75f14f
# ╟─ebcbbf5e-6205-41a8-b31b-a6bed52cda5b
# ╟─4d71efc9-6bf5-4aa6-8e0f-132814350351
# ╟─0b2ed009-e636-40d9-bf9c-8d600f3bae75
# ╟─2b9eae49-a32d-4462-b5ff-8282d05a0a57
# ╟─e301bd11-efb6-431f-a1fd-45b9a2922350
# ╠═6332d014-79e2-451f-b65c-6a5179b7f85c
# ╟─110c5485-51c1-4837-89da-18854267ed7c
# ╠═d5963693-e1c5-49f6-9c07-9187e893cd97
# ╠═a8bd5aae-552f-42c6-8fc8-92ef5bade2f1
# ╟─65cc4d98-7cc8-400e-a41f-0c2cef5bde9e
# ╠═e7536b06-2454-496b-8872-aec29edb37f6
# ╟─80c1040c-88d2-45b5-89e1-b3cb8d6abaf7
# ╠═28a345b9-d911-486c-a715-7d32a4ea41e8
# ╠═aee62f47-0bd7-492f-bc99-2392a16fc1dd
# ╟─d447a8bf-2d42-4f93-a27e-9eed76050348
# ╠═59d380f3-1bc0-41b4-a0e5-2cfbb0eaa333
# ╠═3b219c43-3056-49db-951e-d0e97a4f39db
# ╟─723e81cd-36e1-44ea-99bf-b4cae2be03a2
# ╠═1bf0a983-9592-4a72-a10e-af61316ce6e4
# ╟─b19344d2-d872-4bf5-a75f-1ca481779835
# ╟─e1048e24-e389-4e95-8157-00702f1b41a3
# ╠═fee19d3e-2235-4382-82d6-5a3cb452d317
# ╟─c2967e07-5d8b-4205-bb8b-add323494d86
# ╠═55c08fb9-e29d-4500-84e7-46de7979639e
# ╟─911c9966-e983-46fd-8ca0-8ca06ab42bf0
# ╠═3b8a80a7-4ab5-4b74-9aaa-3603e0ecccb5
# ╟─66ed7bf9-c5b7-4c27-a021-28bb5675c85f
# ╟─19dd77e9-4a88-4c06-9129-0c1391068900
# ╟─3bb85c65-89ed-4113-a498-3a0da01be0b1
# ╟─bfd6f8e4-ea3e-4c87-b7f7-eba889e751fd
# ╠═fad3240f-9fd1-48c5-b1f3-515d9ea038bd
# ╟─09c11c95-a82d-48e0-85ba-6db4a8c03f29
# ╟─57ea4b25-0f40-4816-85ba-05669770885b
# ╟─be1c590d-d70b-40f1-8370-bc294fb29c09
# ╠═bdc608ee-6a52-4130-b4a0-c27d513a0009
# ╠═7e56d621-a572-4077-83be-d3b002f4e808
# ╟─c150327f-b950-4a70-af44-722beee2069c
# ╟─3b2d5624-d365-4595-b95b-52825bc980d0
# ╠═c0d3e1c8-77d9-4f69-8f1a-97b4bec409e4
# ╟─5f6a1560-b22d-4e47-b6ee-7a306b0bc2d1
# ╟─ec3bf735-5fff-44e6-8a79-a03d6eb4538f
# ╠═c9324a48-8cb8-45e1-909a-950333048d28
# ╠═34e13768-7861-4324-8eb1-626275a3b2d1
# ╠═3d57d1f0-4595-48a0-934c-7bc1d5dc15fc
# ╠═4876fb4e-4b00-4984-8932-a6ff576c5be6
# ╠═24e8eafe-e040-4173-b705-4ade869791ba
# ╟─6a1c66ec-6f3d-48a9-9809-2752e3818d18
# ╟─69e15594-1d8e-4ddc-a1da-f968e8b2ee91
# ╠═6a2bafca-e4cb-4ce9-b322-783c8e10fbfd
# ╟─8d09bb0b-40d9-4d79-89a3-ad0d8679b08c
# ╠═2613e6c6-ec09-4574-aecf-bf2a2266bc55
# ╟─269f6922-7f14-45b9-91f2-554a03a92b57
# ╠═7d7db129-4524-48ec-a9a9-61d7b1ee967d
# ╟─8d67dc24-a8b1-4a4f-9811-f46e7b09fc08
# ╠═a5d63b00-4520-4f0c-a2f2-bf6476c4a1a2
# ╠═653a7d32-d0cc-447f-95eb-3b8ee866b6dd
# ╠═218136d4-bd69-4d6b-92b7-98d8f320e17f
# ╠═fb82077d-43b8-4950-96ec-ffc72d400f26
# ╠═e11a771a-18b9-42b2-bb4e-c3eaadbbc3c1
# ╠═b3d653be-450b-49a9-bc3f-118cd5ec9921
# ╠═24346ac7-f089-4278-9364-d21071352500
# ╠═f64849ef-b48f-4783-b213-c108454dcb9a
# ╠═e5c727d8-752c-4a3c-b971-292f4187416e
# ╠═4011a061-68a8-4cc5-875b-7ae27331ad4f
# ╠═4fa98b1c-0777-4f66-a029-d11488aecb2a
# ╠═86321b6a-3d2d-487a-be68-5838e0ff6abf
# ╠═d52cbdaf-f13c-4ecb-8e45-537631de013b
# ╠═e3a2c0b4-5a67-4431-a9db-8116a8fefdbd
# ╠═503a3912-31cc-45f2-b009-dca4700359f7
# ╠═ffb7b8c0-00c6-45a5-8786-8abc6918e997
# ╠═df5ba723-4ea9-4b52-a4ae-120bf37d56df
# ╠═75563ce3-9228-4ca5-89e1-5e704024c209
# ╟─38257847-5ae7-4a5a-a937-22a6729a3640
# ╟─5a9e17d9-2552-48fd-b3ad-0a1e50279953
# ╟─21af98b7-712d-4b25-a9fa-41d008f97962
# ╟─8e514e64-0172-4b47-974d-efaa8e1f4990
# ╟─00199670-40fd-4daa-979f-fb414b116bed
# ╟─f9f5cb83-ab6b-4de2-9ea4-4a4d984f0489
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
