### A Pluto.jl notebook ###
# v0.19.22

#> [frontmatter]
#> title = "Adjoint State Method"
#> description = "This notebook presents a discrete symbolic version of the adjoint state formulation involving the seismic wave equation so that it can act as a reference while implementing methods for full waveform inversion."

using Markdown
using InteractiveUtils

# ╔═╡ c0d3e1c8-77d9-4f69-8f1a-97b4bec409e4
using Symbolics, SymbolicUtils, LinearAlgebra, ChainRules, PlutoUI, PlutoTeachingTools, TikzPictures

# ╔═╡ 41275829-5ad8-470d-b168-6d6c261ceb19
ChooseDisplayMode()

# ╔═╡ 69e5bc3b-494e-4b7d-9fea-ded035d544cc
TableOfContents()

# ╔═╡ 33a3705c-1660-4df6-bfae-23225a55bdc6
md"""# Adjoint State Method
In seismic imaging, the adjoint state method is a numerical method for efficiently computing the gradient of an objective function. 
This notebook presents a discrete version of the adjoint state formulation involving the seismic wave equation so that it can act as a reference while implementing
methods for full waveform inversion.
We focus on the velocity-stress formulation, which is a widely used numerical method for solving wave propagation problems in geophysics. This method involves discretizing the governing equations of motion, such as the wave equation, using finite-difference techniques.

[Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: Pawan Bharadwaj, Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 24e8eafe-e040-4173-b705-4ade869791ba
@variables x z

# ╔═╡ ec3bf735-5fff-44e6-8a79-a03d6eb4538f
md"These commands create three symbolic time variables `t`, and staggered time variables `𝚝`. In this notebook, we will consider three time steps for the leap-frog scheme."

# ╔═╡ c9324a48-8cb8-45e1-909a-950333048d28
@variables t[1:3]

# ╔═╡ f479ca51-fd2a-4261-98b2-e281daf49c5b
@variables 𝚝[1:3]

# ╔═╡ 35909f2f-3787-4c6a-9e74-f284d5eea635
md"## Medium
Created two symbolic variables ρ and μ which can be used to represent the density and shear modulus of the medium as functions of position x."

# ╔═╡ 764b4395-1a02-43b0-bbb9-0189d0d2a9b2
@syms ρ(x) μ(x) μ⁻¹(x)

# ╔═╡ d38130e3-611d-42b2-96b6-a6ca6210308b
md"## State Variables"

# ╔═╡ b4e9aa30-a16c-4b3a-a7be-c2fcf28f0904
md"In the velocity-stress formulation, the seismic wavefield is described in terms of two variables: the particle velocity and the stress. The `@syms` macro is used in the `Symbolics.jl` package to declare symbolic variables. In the given code snippet, the @syms macro is used to declare two symbolic variables `v(x, t)` and `σ(x, t)` representing velocity and stress, respectively, as functions of space (`x`) and time (`t`)."

# ╔═╡ 857a5e06-2960-4ed4-826c-6bf3291387f0
@syms v(x, t) σ(x, t)

# ╔═╡ 61505a55-fc9e-4077-8365-8e166871b482
md"The symbol `f` represents the body force density, which is a quantity that describes the force per unit volume acting on a material."

# ╔═╡ 93ce7686-a297-438f-bfab-9ad1f9e9be10
@syms f(x, t)

# ╔═╡ 76fd526b-450a-4db1-bf9a-bc7ae59b895c
md"Create vectors with velocity and stress fields at each time value in 𝚝 and t, respectively. "

# ╔═╡ c0624ec9-1d80-461d-a10a-5ec669814fb5
V = collect(map(t -> v(x, t), 𝚝)) # velocity

# ╔═╡ 32009437-f2b8-4520-851b-a5929b52d272
Σ = collect(map(t -> σ(x, t), t)) # stress

# ╔═╡ c0bb57b5-9d62-40b1-9256-648a984dbbae
md"## Forcing"

# ╔═╡ ebefecd3-9dd2-457e-bf9b-97d4db9e983e
F = collect(map(t -> f(x, t), t))

# ╔═╡ 5803bfab-354b-4f48-8864-7f3cb3ac3c08
md"## State Equations"

# ╔═╡ bc20a677-5d50-4a80-b466-2b7edd5feac4
md"### Momentum Equation"

# ╔═╡ 33871626-5b9b-4718-a974-d94d7b048582
md"""
This line defines the initial momentum of the system, which is equal to the difference between the initial velocity of the particle v(x, 𝚝[1]) and its velocity at time t=0 (which is zero in this case) multiplied by the density of the medium ρ(x) and then subtracting the force acting on the particle at time t=0 f(0, x). It represents the balance of forces and momentum of the system at the initial time.
"""

# ╔═╡ e5a37530-a0f2-470d-aada-39eaf3cd6054
M1 = (v(x, 𝚝[1]) - v(x, 0)) * ρ(x) - f(x, 0)

# ╔═╡ 84d8f063-2aad-4d57-9bdf-d44e37f173c7
md"""
We will now construct a vector that represents the time derivative of the velocity field. The first few elements can be computed using the `diff()` function applied to the `V` array. The last element is the time derivative of `v(x, t)` evaluated at the final time `𝚃`, which is computed as the difference between `v(x, 𝚃)` and `v(x, last(𝚝))`.
"""

# ╔═╡ 3d57d1f0-4595-48a0-934c-7bc1d5dc15fc
@syms 𝚃 # final time

# ╔═╡ 91f58348-d476-4799-ab01-a2b7b89a344c
DₜV = [diff(V)..., v(x, 𝚃) - v(x, last(𝚝))] # 

# ╔═╡ c959b590-aacf-4504-bb80-bd915acb0b3f
@syms Dₓ(v) # a function to return spatial derivative (with respect to x)

# ╔═╡ 9242d8be-a259-4a72-bbb7-4958884937c9
Meqs = vcat(M1, (DₜV .* ρ(x) - Dₓ.(Σ) - F));

# ╔═╡ f54c2aa7-2afb-4806-b932-417e3b4a41e5
Meqs ~ 0

# ╔═╡ c5e41f1c-a734-4502-a60f-4ec0f9819e89
md"### Constitutive Relation
In linear elasticity, for example, the constitutive relation is given by Hooke's law, which states that stress is proportional to strain. First, similar to velocity, we will construct a vector that represents the time derivative of the stress field."

# ╔═╡ d42c3324-4517-4ac1-b9d9-1af0656bf70e
𝚃

# ╔═╡ de66df57-c6d1-4dba-8e6f-d1d8fd6b59ff
DₜΣ = [σ(x, t[1]) - σ(x, 0), diff(Σ)..., σ(x, 𝚃) - σ(x, t[3])] # lives on t

# ╔═╡ b2562463-09d0-40fc-9802-b49e105f946c
md"This equation states that the time derivative of stress minus the gradient of the stress scaled by the shear modulus is equal to zero. This is a simplified form of the constitutive relation."

# ╔═╡ e73aa232-e0a9-4c7e-b399-6b0dee065a25
Ceqs = (DₜΣ * μ⁻¹(x) - Dₓ.(vcat(V, v(x, 𝚃))));

# ╔═╡ 9678bb33-8e15-49ad-95b8-c8d4ab75f14f
Ceqs ~ 0

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

# ╔═╡ a8bd5aae-552f-42c6-8fc8-92ef5bade2f1
obj(v, v₀) = (v - v₀)^2

# ╔═╡ 65cc4d98-7cc8-400e-a41f-0c2cef5bde9e
md"The final value of the objective function."

# ╔═╡ a6f25519-efdb-4430-a014-6a7ee0360df6
J = sum(collect(map(t -> obj(v(x, t), v₀(x, t)), 𝚝)))

# ╔═╡ 80c1040c-88d2-45b5-89e1-b3cb8d6abaf7
md"## Adjoint State Variables
Adjoint state variables are introduced in the context of optimization problems. In this case, these variables are used to find gradients of an objective function with respect to the parameters mass density and shear modulus."

# ╔═╡ 28a345b9-d911-486c-a715-7d32a4ea41e8
@syms u(x, t) τ(x, t)

# ╔═╡ d5095bc5-6a00-4e0e-9669-68e21deac773
U = vcat(collect(map(t -> u(x, t), 𝚝)), u(x, 𝚃))  # one u for each of Meqs

# ╔═╡ cc3a1c12-25b9-450a-a23d-8557fe351f83
T = vcat(collect(map(t -> τ(x, t), t)), τ(x, 𝚃)) # one for each Ceqs

# ╔═╡ b19344d2-d872-4bf5-a75f-1ca481779835
md"## Lagrangian"

# ╔═╡ e1048e24-e389-4e95-8157-00702f1b41a3
md"The Lagrangian component, for momentum equations `Meqs`, is given by"

# ╔═╡ fee19d3e-2235-4382-82d6-5a3cb452d317
L₁ = sum(U .* Meqs)

# ╔═╡ 0eb7d977-29fc-4be1-9e2d-16285431dcfb
md"The function `r_Dₓ_transpose` applies the transpose differential operator rule and returns the sum of the resulting terms after simplification."

# ╔═╡ c2967e07-5d8b-4205-bb8b-add323494d86
md"The Lagrangian component, for constitutive equations `Ceqs`, is given by"

# ╔═╡ 55c08fb9-e29d-4500-84e7-46de7979639e
L₂ = simplify(sum(T .* Ceqs), expand=true)

# ╔═╡ 911c9966-e983-46fd-8ca0-8ca06ab42bf0
md"Adding all the components of the Lagrangian together gives"

# ╔═╡ 66ed7bf9-c5b7-4c27-a021-28bb5675c85f
md"## Adjoint State Equations"

# ╔═╡ 69e15594-1d8e-4ddc-a1da-f968e8b2ee91
md"This code computes the gradient of the Lagrangian with respect to the  velocity field to obtain adjoint equations."

# ╔═╡ 8d09bb0b-40d9-4d79-89a3-ad0d8679b08c
md"This code computes the gradient of the Lagrangian with respect to the stress field to obtain adjoint equations."

# ╔═╡ bfd6f8e4-ea3e-4c87-b7f7-eba889e751fd
md"### Final Condition
To obtain the final condition, we shall compute the gradient of the Lagrangian with respect to `v(x, 𝚃)`."

# ╔═╡ 09c11c95-a82d-48e0-85ba-6db4a8c03f29
md"The solution of these adjoint equations is often obtained by using a time-reversed version of the original numerical solver, known as the backpropagation."

# ╔═╡ be1c590d-d70b-40f1-8370-bc294fb29c09
md"## Parameter Gradients
Lets compute the gradient of L with respect to ρ and μ."

# ╔═╡ a5d63b00-4520-4f0c-a2f2-bf6476c4a1a2
nt = 3

# ╔═╡ 653a7d32-d0cc-447f-95eb-3b8ee866b6dd
Vs = collect(map(t -> v(x, t), 𝚝)) # Simulated velocity fields

# ╔═╡ fb82077d-43b8-4950-96ec-ffc72d400f26
Us = reverse(collect(map(t -> u(x, t), 𝚝))) # Simulated adjoint velocity fields

# ╔═╡ f64849ef-b48f-4783-b213-c108454dcb9a
simplify(sum([Us[it] * (-Vs[nt-it] + Vs[nt-it+1]) for it in 1:nt-1]), expand=true) + Us[nt] * Vs[1]

# ╔═╡ 8e46496e-d44b-4c1e-96c1-e8b29d792e32
Ss = collect(map(t -> σ(x, t), t)) # Simulated stress fields

# ╔═╡ 3e2cb699-c73e-40c5-9bf7-83d0486eeed7
Ts = reverse(collect(map(t -> τ(x, t), t))) # Simulated adjoint stress fields

# ╔═╡ 0cbe023d-ba95-43f5-b06b-55e0a5333404
simplify(sum([Ts[it] * (-Ss[nt-it] + Ss[nt-it+1]) for it in 1:nt-1]), expand=true) + Ts[nt] * Ss[1]

# ╔═╡ 3b2d5624-d365-4595-b95b-52825bc980d0
md"## Appendix"

# ╔═╡ 73cec834-0b81-4b93-a00f-4e953d93b5de
md"This function implements the rule for transposing the spatial differential operator Dₓ in an expression using the adjoint relationship Dₓᵀ = -Dₓ. It takes as input a symbolic expression L that may contain Dₓ operators, and returns the expression with all Dₓ operators transposed using the Dₓᵀ operator."

# ╔═╡ 3f7e7312-236c-4710-a901-fc88e4c5afca
function r_Dₓ_transpose(L)
  r1 = @rule (~~z) * Dₓ(~x) * (~~y) => -Dₓ(prod(~~z) * prod(~~y)) * ~x
  return mapreduce(+, arguments(simplify(L, expand=true))) do l
    return (r1(l) === nothing) ? l : r1(l)
  end
end

# ╔═╡ 934ae935-70ad-439c-933c-2ef8a78ff5b4
L1 = r_Dₓ_transpose(L₁)

# ╔═╡ 5d69a093-d57d-499c-a5a4-a1248d1811c6
L2 = r_Dₓ_transpose(L₂)

# ╔═╡ 3b8a80a7-4ab5-4b74-9aaa-3603e0ecccb5
L = L1 + L2 + J

# ╔═╡ 6a2bafca-e4cb-4ce9-b322-783c8e10fbfd
∂vL = broadcast(V[1:end]) do v
  Differential(v)(L) |> expand_derivatives
end;

# ╔═╡ 19dd77e9-4a88-4c06-9129-0c1391068900
∂vL ~ 0

# ╔═╡ 2613e6c6-ec09-4574-aecf-bf2a2266bc55
∂σL = broadcast(Σ[1:end]) do σ
  Differential(σ)(L) |> expand_derivatives
end;

# ╔═╡ 3bb85c65-89ed-4113-a498-3a0da01be0b1
∂σL ~ 0

# ╔═╡ 7d7db129-4524-48ec-a9a9-61d7b1ee967d
∂vTL = Differential(v(x, 𝚃))(L) |> expand_derivatives;

# ╔═╡ 90b6fc69-0056-4971-9395-10c601ce4ef4
∂vTL ~ 0

# ╔═╡ a1f58735-5fb7-4920-ba2c-1821c013d4d9
∂σTL = Differential(σ(x, 𝚃))(L) |> expand_derivatives;

# ╔═╡ 1ede1ceb-f286-4a71-bea1-79e0fc3f74ff
∂σTL ~ 0

# ╔═╡ 7e56d621-a572-4077-83be-d3b002f4e808
# gradient w.r.t. mass density
∇ρ = substitute(Differential(ρ(x))(L) |> expand_derivatives, [u(x, 𝚃) => 0, v(x, 0) => 0])


# ╔═╡ c150327f-b950-4a70-af44-722beee2069c
# gradient w.r.t. shear modulus
∇μ⁻¹ = substitute(Differential(μ⁻¹(x))(L) |> expand_derivatives, [τ(x, 𝚃) => 0, σ(x, 0) => 0])

# ╔═╡ 38257847-5ae7-4a5a-a937-22a6729a3640
md"### Tikz"

# ╔═╡ 4483493d-cc6b-4e7a-b244-5409f1da013e
get_field_tikz(v, s, f="f") = L"""
    \tikzstyle{vertex}=[circle,minimum size=20pt,inner sep=0pt]

  \foreach \name/\x in {%$(v)(0)/1, %$(v)(\mathtt{t}_1)/5, %$(v)(\mathtt{t}_2)/9, %$(v)(\mathtt{t}_3)/13, %$(v)(\mathtt{T})/17}
    \node[vertex,fill=black!25,] (v-\x) at (\x,0) {$\name$};

 \foreach \name/\x in {%$(s)(0)/3, %$(s)(t_1)/7, %$(s)(t_2)/11, %$(s)(t_3)/15}
    \node[vertex,fill=red!25] (s-\x) at (\x,-1) {$\name$};

 \foreach \name/\x in {%$(f)(0)/3, %$(f)(\mathtt{t}_1)/7, %$(f)(\mathtt{t}_2)/11, %$(f)(\mathtt{t}_3)/15}
    \node[vertex,fill=blue!25] (f-\x) at (\x,1) {$\name$};
  """

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

# ╔═╡ 5a9e17d9-2552-48fd-b3ad-0a1e50279953
# t1, t2, t3, tv are texts
# s1 and s2 are labels with sizes
plot_state_tikz() = TikzPicture(L"""
   %$(get_field_tikz("v", L"\sigma", "f"))


 \draw[->,thick, red] (s-3) -- (s-7) node[above, midway, red] {1};
\draw[->,thick, red] (v-5) -- (s-7) node[above, midway, red] {1};
\draw[->,thick, blue] (s-7) -- (v-9) node[below, midway, blue] {2};
	\draw[->,thick, blue] (v-5) -- (v-9) node[above, midway, blue] {2};

\node[above right=-2mm of s-3] {=0};
\node[above right=-2mm of v-1] {=0};
  """, options=tikz_default_options, preamble=tikz_preamble, width="20cm")

# ╔═╡ 4d71efc9-6bf5-4aa6-8e0f-132814350351
plot_state_tikz()

# ╔═╡ 21af98b7-712d-4b25-a9fa-41d008f97962
# t1, t2, t3, tv are texts
# s1 and s2 are labels with sizes
plot_adjstate_tikz() = TikzPicture(L"""
   %$(get_field_tikz("u", L"\tau", "g"))




 \draw[->,thick, red] (s-15) -- (s-11) node[above, midway, red] {1};
\draw[->,thick, red] (v-13) -- (s-11) node[above, midway, red] {1};
\draw[->,thick, blue] (s-11) -- (v-9) node[below, midway, blue] {2};
	\draw[->,thick, blue] (v-13) -- (v-9) node[above, midway, blue] {2};

\node[above right=-2mm of s-15] {=0};
\node[above right=-2mm of v-17] {=0};


  """, options=tikz_default_options, preamble=tikz_preamble, width="20cm")

# ╔═╡ 57ea4b25-0f40-4816-85ba-05669770885b
plot_adjstate_tikz()

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
ChainRules = "~1.48.0"
PlutoTeachingTools = "~0.2.8"
PlutoUI = "~0.7.50"
SymbolicUtils = "~1.0.4"
Symbolics = "~5.1.0"
TikzPictures = "~3.4.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "29e65c331f97db9189ef00a4c7aed8127c2fd2d4"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.27.10"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "38911c7737e123b28182d89027f4216cfc8a9da7"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.3"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "Statistics", "StructArrays"]
git-tree-sha1 = "7d20c2fb8ab838e41069398685e7b6b5f89ed85b"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.48.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "d730914ef30a06732bdd9f763f6cc32e92ffbff1"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "89a9db8d28102b094992472d333674bd1a83ce2a"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "13027f188d26206b9e7b863036f87d2f2e7d013a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.87"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "988e2db482abeb69efc76ae8b6eba2e93805ee70"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.15"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "8b84876e31fa39479050e2d3395c4b3b210db8b0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.6"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExprTools]]
git-tree-sha1 = "c1d06d129da9f55715c6c212866f5b1bddc5fa00"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.9"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "fc86b4fd3eff76c3ce4f5e96e2fdfa6282722885"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.0.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

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
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "47f0f03eddecd7ad59c42b1dd46d5f42916aff63"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.11"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "432b5b03176f8182bd6841fbfc42c718506a2d5f"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.15"

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
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "3f91cd3f56ea48d4d2a75c2a65455c5fc74fa347"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "6a125e6a4cb391e0b9adbd1afa9e771c2179f8ef"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.23"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

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
git-tree-sha1 = "98dc144f1e0b299d49e8d23e56ad68d3e4f340a5"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.20"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

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
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "eaa98afe2033ffc0629f9d0d83961d66a021dfcc"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3295d296288ab1a0a2528feb424b854418acff57"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.2.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

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
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "8c8b07296990c12ac3a9eb9f74cd80f7e81c16b7"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.9"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "e11443687ac151ac6ef6699eb75f964bed8e1faa"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "0.87.0+2"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

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

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "140cddd2c457e4ebb0cdc7c2fd14a7fbfbdf206e"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.3"

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
git-tree-sha1 = "feafdc70b2e6684314e188d95fe66d116de834a7"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.2"

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

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "f139e81a81e6c29c40f1971c9e5309b09c03f2c3"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.6"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SnoopPrecompile", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "392d3e28b05984496af37100ded94dc46fa6c8de"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.91.7"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "e61e48ef909375203092a6e83508c8416df55a83"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.2.0"

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
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "63e84b7fdf5021026d0f17f76af7c57772313d99"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.21"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "5cb1f963f82e7b81305102dd69472fcd3e0e1483"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.5"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "fce1fd0b13f860128c8b8aab0bab475eeeeb7994"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.1.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Tectonic]]
deps = ["Pkg"]
git-tree-sha1 = "0b3881685ddb3ab066159b2ce294dc54fcf3b9ee"
uuid = "9ac5f52a-99c6-489f-af81-462ef484790f"
version = "0.8.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "Tectonic"]
git-tree-sha1 = "4e75374d207fefb21105074100034236fceed7cb"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.4.2"

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
git-tree-sha1 = "7bc1632a4eafbe9bd94cf1a784a9a4eb5e040a91"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "d5f4ec8c22db63bd3ccb239f640e895cfde145aa"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═41275829-5ad8-470d-b168-6d6c261ceb19
# ╠═69e5bc3b-494e-4b7d-9fea-ded035d544cc
# ╟─33a3705c-1660-4df6-bfae-23225a55bdc6
# ╠═24e8eafe-e040-4173-b705-4ade869791ba
# ╟─ec3bf735-5fff-44e6-8a79-a03d6eb4538f
# ╠═c9324a48-8cb8-45e1-909a-950333048d28
# ╠═f479ca51-fd2a-4261-98b2-e281daf49c5b
# ╟─35909f2f-3787-4c6a-9e74-f284d5eea635
# ╠═764b4395-1a02-43b0-bbb9-0189d0d2a9b2
# ╟─d38130e3-611d-42b2-96b6-a6ca6210308b
# ╟─b4e9aa30-a16c-4b3a-a7be-c2fcf28f0904
# ╠═857a5e06-2960-4ed4-826c-6bf3291387f0
# ╟─61505a55-fc9e-4077-8365-8e166871b482
# ╠═93ce7686-a297-438f-bfab-9ad1f9e9be10
# ╟─76fd526b-450a-4db1-bf9a-bc7ae59b895c
# ╠═c0624ec9-1d80-461d-a10a-5ec669814fb5
# ╠═32009437-f2b8-4520-851b-a5929b52d272
# ╟─c0bb57b5-9d62-40b1-9256-648a984dbbae
# ╠═ebefecd3-9dd2-457e-bf9b-97d4db9e983e
# ╟─5803bfab-354b-4f48-8864-7f3cb3ac3c08
# ╟─bc20a677-5d50-4a80-b466-2b7edd5feac4
# ╟─33871626-5b9b-4718-a974-d94d7b048582
# ╠═e5a37530-a0f2-470d-aada-39eaf3cd6054
# ╟─84d8f063-2aad-4d57-9bdf-d44e37f173c7
# ╠═91f58348-d476-4799-ab01-a2b7b89a344c
# ╠═3d57d1f0-4595-48a0-934c-7bc1d5dc15fc
# ╠═c959b590-aacf-4504-bb80-bd915acb0b3f
# ╠═9242d8be-a259-4a72-bbb7-4958884937c9
# ╠═f54c2aa7-2afb-4806-b932-417e3b4a41e5
# ╟─c5e41f1c-a734-4502-a60f-4ec0f9819e89
# ╠═d42c3324-4517-4ac1-b9d9-1af0656bf70e
# ╠═de66df57-c6d1-4dba-8e6f-d1d8fd6b59ff
# ╟─b2562463-09d0-40fc-9802-b49e105f946c
# ╠═e73aa232-e0a9-4c7e-b399-6b0dee065a25
# ╠═9678bb33-8e15-49ad-95b8-c8d4ab75f14f
# ╟─0b2ed009-e636-40d9-bf9c-8d600f3bae75
# ╟─4d71efc9-6bf5-4aa6-8e0f-132814350351
# ╟─2b9eae49-a32d-4462-b5ff-8282d05a0a57
# ╟─e301bd11-efb6-431f-a1fd-45b9a2922350
# ╠═6332d014-79e2-451f-b65c-6a5179b7f85c
# ╟─110c5485-51c1-4837-89da-18854267ed7c
# ╠═a8bd5aae-552f-42c6-8fc8-92ef5bade2f1
# ╟─65cc4d98-7cc8-400e-a41f-0c2cef5bde9e
# ╠═a6f25519-efdb-4430-a014-6a7ee0360df6
# ╟─80c1040c-88d2-45b5-89e1-b3cb8d6abaf7
# ╠═28a345b9-d911-486c-a715-7d32a4ea41e8
# ╠═d5095bc5-6a00-4e0e-9669-68e21deac773
# ╠═cc3a1c12-25b9-450a-a23d-8557fe351f83
# ╟─b19344d2-d872-4bf5-a75f-1ca481779835
# ╟─e1048e24-e389-4e95-8157-00702f1b41a3
# ╠═fee19d3e-2235-4382-82d6-5a3cb452d317
# ╟─0eb7d977-29fc-4be1-9e2d-16285431dcfb
# ╠═934ae935-70ad-439c-933c-2ef8a78ff5b4
# ╟─c2967e07-5d8b-4205-bb8b-add323494d86
# ╠═55c08fb9-e29d-4500-84e7-46de7979639e
# ╠═5d69a093-d57d-499c-a5a4-a1248d1811c6
# ╟─911c9966-e983-46fd-8ca0-8ca06ab42bf0
# ╠═3b8a80a7-4ab5-4b74-9aaa-3603e0ecccb5
# ╟─66ed7bf9-c5b7-4c27-a021-28bb5675c85f
# ╟─69e15594-1d8e-4ddc-a1da-f968e8b2ee91
# ╠═6a2bafca-e4cb-4ce9-b322-783c8e10fbfd
# ╠═19dd77e9-4a88-4c06-9129-0c1391068900
# ╟─8d09bb0b-40d9-4d79-89a3-ad0d8679b08c
# ╠═2613e6c6-ec09-4574-aecf-bf2a2266bc55
# ╠═3bb85c65-89ed-4113-a498-3a0da01be0b1
# ╟─bfd6f8e4-ea3e-4c87-b7f7-eba889e751fd
# ╠═7d7db129-4524-48ec-a9a9-61d7b1ee967d
# ╠═a1f58735-5fb7-4920-ba2c-1821c013d4d9
# ╠═90b6fc69-0056-4971-9395-10c601ce4ef4
# ╠═1ede1ceb-f286-4a71-bea1-79e0fc3f74ff
# ╟─09c11c95-a82d-48e0-85ba-6db4a8c03f29
# ╠═57ea4b25-0f40-4816-85ba-05669770885b
# ╟─be1c590d-d70b-40f1-8370-bc294fb29c09
# ╠═a5d63b00-4520-4f0c-a2f2-bf6476c4a1a2
# ╠═653a7d32-d0cc-447f-95eb-3b8ee866b6dd
# ╠═fb82077d-43b8-4950-96ec-ffc72d400f26
# ╠═7e56d621-a572-4077-83be-d3b002f4e808
# ╠═f64849ef-b48f-4783-b213-c108454dcb9a
# ╠═8e46496e-d44b-4c1e-96c1-e8b29d792e32
# ╠═3e2cb699-c73e-40c5-9bf7-83d0486eeed7
# ╠═c150327f-b950-4a70-af44-722beee2069c
# ╠═0cbe023d-ba95-43f5-b06b-55e0a5333404
# ╟─3b2d5624-d365-4595-b95b-52825bc980d0
# ╠═c0d3e1c8-77d9-4f69-8f1a-97b4bec409e4
# ╟─73cec834-0b81-4b93-a00f-4e953d93b5de
# ╠═3f7e7312-236c-4710-a901-fc88e4c5afca
# ╟─38257847-5ae7-4a5a-a937-22a6729a3640
# ╠═4483493d-cc6b-4e7a-b244-5409f1da013e
# ╠═5a9e17d9-2552-48fd-b3ad-0a1e50279953
# ╠═21af98b7-712d-4b25-a9fa-41d008f97962
# ╟─8e514e64-0172-4b47-974d-efaa8e1f4990
# ╟─00199670-40fd-4daa-979f-fb414b116bed
# ╟─f9f5cb83-ab6b-4de2-9ea4-4a4d984f0489
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
