### A Pluto.jl notebook ###
# v0.19.38

#> [frontmatter]
#> layout = "layout.jlhtml"
#> title = "Probabilistic Inverse Problems"
#> tags = ["module1"]
#> description = "Posterior Information = Conjunction(Prior Information, Theoretical Information)."

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

# ╔═╡ 04d19c81-08f8-4f07-a0c7-f51c1b39d271
using PlutoUI, Plots, Distributions, Measures, PlutoTeachingTools, StatsPlots, Symbolics, LaTeXStrings, TikzPictures

# ╔═╡ 02cd8bb4-d403-47f9-9de9-35d7b2b82bb8
ChooseDisplayMode()

# ╔═╡ 4c8ba5dc-a019-46ab-b7cc-b03cf4cee5f8
TableOfContents()

# ╔═╡ 91c82c81-8cfc-4b54-bed7-0592f5b39351
md"""
# Probabilistic Inverse Problems
- The fundamental idea of probabilistic inverse theory is to express all knowledge concerning model parameters, data, and the physical theory in terms of probabilities.

- Frequentist interpretation of probability is simple: it considers a large number of repeatable experiments. 

- When we don't consider repeatable experiments, we adopt the Bayesian interpretation of probabilities, where they represent state of information.
"""

# ╔═╡ 1b14a442-211e-40a6-b8d1-3854a69ef682
md"""## Repeating Experiments
The data measured at the receivers are no longer certain. They change from one experiment to another. For example, we can measure the traveltime at a specific receiver 10 times after repeating the experiment.
"""

# ╔═╡ c57d72de-920d-48b5-b768-661cd2435bcd
md"""
As we think of data as a set of random variables, we can plot the (normalized) histogram to roughly understand its probability density function.
"""

# ╔═╡ 65ed497b-18ed-4379-b220-439b51f3cb86
md"## Joint Probability
Usually, our data are composed of multiple random variables.

Independent if the conditional distribution is in.
"

# ╔═╡ c3dc0533-0cb5-4c2f-a1e6-6c18dde7674a
md"Are the elements of the data vector independent? Most likely not."

# ╔═╡ 29bf0e42-8dd0-41e9-a646-a7714440d09b
A = randn(2, 2) # mixing matrix, lets call this op 

# ╔═╡ e5d5ac74-bbfd-42e1-aeaf-5fbd17d59625
md"## Coordinate Transformation"

# ╔═╡ a2461a13-76f4-4a18-bcf4-fc9e7a18e5bb
plot(histogram(rand(Uniform(0.5, 4), 10000), normalize=true, label=nothing, xlabel="bulk modulus", c=:gray, xlim=(0, 5)), histogram(inv.(rand(Uniform(0.5, 4), 10000)), normalize=true, label=nothing, xlabel="compressibility", c=:gray), size=(700, 300), margin=1cm)

# ╔═╡ 54085478-12df-4d02-b2ed-ac9a373f1b91
mrange = range(0, 5, length=100)

# ╔═╡ 4465ed36-d7f8-4407-9270-0ec20dd7e08c
md"""
## Multi-dimensional Gaussian Distribution
```math
\frac{1}{\sqrt((2\pi)\mathrm{det}\,C_m)}\,\exp\left(-\frac{1}{2}\,(m-m_0)^\top\,C_m^{-1}\,(m-m_0)\right)
```
"""

# ╔═╡ b3212b47-a114-4a25-88f7-5a3cd3cb6542
@bind cov12 Slider(range(-0.7, stop=0.7, length=100), show_value=true)

# ╔═╡ 88060b9a-3c70-4bed-8b0c-ce6c5b9a3061
Cₘ = [1 cov12
    cov12 1]

# ╔═╡ 47c3f46f-158e-4eb6-aa8d-49c018168826
pm2D = MvNormal([3, 2], Cₘ)

# ╔═╡ 43363cb2-98d5-41d3-b276-97478829391a
Z = [pdf(pm2D, [x, y]) for y in mrange, x in mrange];

# ╔═╡ b2c008d8-48aa-4c32-b988-3552404a25b6
md"""
## Marginal densities
Given joint information about $m_1$ and $m_2$, if we are only interested in one of the parameters, say $m_1$, then we can marginalize over $m_2$.

```math
\rho(m_1) = \int_{m_2}\rho(m_1, m_2)\,\mathrm{d}m_2
```
is the marginal probability density that contains all the information on $m_1$ that is independent of the information on $m_2$.
"""

# ╔═╡ ab6453f1-9746-4825-84bc-32b552538e23
begin
    x, y = rand(Normal(), 1000), rand(Uniform(-10, 15), 1000)

    layout = @layout [a _
        b{0.8w,0.8h} c]

    default()
    plot(layout=layout, link=:both, size=(500, 500), margin=0.3cm)
    scatter!(x, y, subplot=2, framestyle=:box, fillcolor=:lightgrey, markercolor=:gray, grid=false, legend=false, xlabel="m₁", ylabel="m₂")
    histogram!([x y], subplot=[1 3], orientation=[:v :h], framestyle=:none, label=nothing, c=:gray)
end

# ╔═╡ d4123859-1228-444f-ac49-ec861aae3bdb
md"""
## Mean of measurements
Consider the mean of measurements that were collected by repeating the experiment. What would be the distribution of the mean? You may slide to increase the number of independent experiments.
"""

# ╔═╡ c7b84588-8e3a-4206-8ca5-0a7a85d395d5
@bind d_distribution Select([Normal(1, 1), Uniform(-2, 2)])

# ╔═╡ e8967bc9-daff-48b0-9846-f839848e3142
@bind nm Slider(range(1, stop=25, step=1), show_value=true)

# ╔═╡ 9e0ff7a6-3984-4be0-858c-b08e95bd46b1
dobs_mean = mean([rand(d_distribution, 1000) for i in 1:nm]);

# ╔═╡ 7b66bc33-4d0c-4f21-9a9c-7a208708d95d
begin
	histogram(dobs_mean, normalize=true, size=(500, 250), label=nothing, xlabel="Measurement", c=:gray, xlim=(-5, 5));
	vline!([1.0], w=3, c=:blue, label="True Measurement", margin=1cm, legend=:topleft);
end

# ╔═╡ 263c76d0-051b-4349-9f2f-c6ee82ad02f6
md"""## Theoretical Information
Θ(d, m) is the joint probability density  describing the correlations that correspond to our physical theory, together with inherent uncertainites of the theory (in this case, knowledge on the the outgoing wavenumber vector is uncertain, which is attributed to the raytracing)
"""

# ╔═╡ 3add2ee1-1de6-490d-a9b7-1f01d968cf29
md"""
```math
\Theta(d, m) = \theta(d\,|\,m)\,μ_m(m)
```
"""

# ╔═╡ 70df73b9-e418-47da-81e1-735fda13f131
md"""
```math
\theta(d\,|\,m) = 
\text{const.}\,\exp\left(-\frac{1}{2}\,(d-g(m))^\top\,C_T^{-1}\,(d-g(m))\right)
```
"""

# ╔═╡ a7c12637-c990-49a8-b60e-ddae71fda3e3
md"## Measurements"

# ╔═╡ a305940b-6699-4e04-a1d8-8182266b14a4
md"""
$\rho_D(d)$
"""

# ╔═╡ bedc65e0-e4bb-47c6-82cd-10f1cac4d9db
md"""## Prior Information
```math
\rho_M(m)
```
"""

# ╔═╡ 3cd564e1-cd99-4c25-8ae4-85904f997a4c
md"## Joint Prior Information
```math
\rho(d, m) = \rho_D(d)\,\rho_M(m)
```"

# ╔═╡ ac9bbed8-9d8f-461e-bea5-9337cec3ec80
md"""
## Combining states of information
Different experts or different datasets often provide different pieces of information that may be encoded by the
probability densities ρ1 (m) and ρ2 (m), respectively. 

```math
(\rho_1\land\rho_2)(m) = k\,\frac{\rho_1(m)\rho_2(m)}{\rho_h(m)}
```
The
homogeneous (uninformative) probability density in equation (3.47) ensures that the conjunction of information
ρ1 with the no-information ρh does not yield any new information, that is
```math
\rho_1 \land \rho_h = \rho_1
```
"""

# ╔═╡ da70f650-c306-461d-a607-8d2384946088
pm2D_1 = MvNormal([1, 1], [1 0; 0 1])

# ╔═╡ 8dca7b1e-7142-4bdc-81bf-645cc8c264c3
pm2D_2 = MvNormal([3, 3], [1 0.5; 0.5 1])

# ╔═╡ 6f7955f7-99fa-482c-8544-9f5bb9e75452
Z1 = [pdf(pm2D_1, [x, y]) for y in mrange, x in mrange];

# ╔═╡ 7610159c-f519-4dd6-bf49-e9015dfda7d5
Z2 = [pdf(pm2D_2, [x, y]) for y in mrange, x in mrange];

# ╔═╡ 83f55514-4c1e-4587-9dc9-576351dc0379
md"## Shannon's Measure"

# ╔═╡ 849ecb20-7d5b-48da-9fbc-3a8bfcb4fa9b
md"""
## Appendix
"""

# ╔═╡ 90a4f7a0-7a94-4522-a4e9-99f8777e8f53
# make a list of interesting distributions
dist_choices1 = Select([Normal(1, 0.1) => "Normal distribution", Normal(1.1, 0.1) => "Normal with bias", Uniform(0.9, 1.1) => "Uniform distribution (ignorance)", Dirac(1.0) => "δ distribution (perfect knowledge)"]);

# ╔═╡ 56c7e80d-068e-4c4a-a1fa-64e26daf147f
md"""
The distribution of each of the random variables in the data is usually unknown, and may not be Normal as per the default choice below. For example, you may choose another one here: $(@bind pdata1 dist_choices1)
"""

# ╔═╡ 760ad648-a0cf-4ec7-96bf-362346e34ba6
scatter(rand(pdata1, 10), ylim=(0.5, 1.5), xlabel="# Experiment", ylabel="Travel time (s)", label=nothing, size=(500, 250), c=:gray)

# ╔═╡ 9b3ac9a2-5bc0-43e3-a8f9-688490e19674
d1obs = rand(pdata1, 10000) # sample some data from a given distribution

# ╔═╡ f460a7f2-f812-4bdf-ba13-8a905f51ea55
begin
	pd1 = histogram(d1obs, normalize=true, size=(500, 250), label=nothing, xlabel="Measurement", c=:gray);
	plot!(range(0.6, 1.4, step=0.01), pdf.(pdata1, range(0.6, 1.4, step=0.01)), w=2, label="pdf");
	vline!([1.0], w=3, c=:blue, label="True Measurement", margin=1cm);
end

# ╔═╡ 58e5869f-03a7-4698-863b-c6eab02e927f
@bind pdata2 dist_choices1

# ╔═╡ 9e6e08cb-d26f-4f17-a662-1e5097cf7625
d2obs = rand(pdata2, 10000)

# ╔═╡ 859865b7-a2f8-4890-a127-7b96ec0513c7
dmix = mapslices(x -> A * x, hcat(d1obs, d2obs), dims=2);

# ╔═╡ 0eac2ee9-ce0c-46e2-9796-d3bb1b0eb40f
begin
    X = range(-8, 8, length=100)
    Y = range(-8, 8, length=100)
end

# ╔═╡ 22c3fe50-0ace-4a6b-a4b4-7e6f63a45941
md"### Locate Source Experiment"

# ╔═╡ 0e768c08-c806-4de9-b459-f52065b6ce0a
sxgrid = range(-2.0, stop=2.0, length=250)

# ╔═╡ 72ddf27e-de41-41f4-aeae-12f2ead3d89c
Tgrid = range(0.95, 1.05, length=250)

# ╔═╡ 7a0bf6ba-0353-42b5-b86b-1d43e22a22b9
cex = 5 # velocity in (km/sec)

# ╔═╡ 3332d7c6-8971-404b-86ae-7eaf03d0768d
@syms T 𝐱 𝐳

# ╔═╡ 1c2abe11-ac0f-413c-93da-e521ebf5ae36
T_expression = sqrt(abs2(𝐱) + abs2(𝐳)) * inv(cex)

# ╔═╡ b447013f-4e0c-44b1-ba17-167a34cd3f7b
get_traveltime = build_function(T_expression, 𝐱, 𝐳; expression=Val{false})

# ╔═╡ 459e9902-05c5-4c53-ad99-7f68143e205a
sz_expression = sqrt(abs2(cex * T) - abs2(𝐱))

# ╔═╡ 925c71dd-530f-4bd8-8685-d66747c277ec
# return sz for a given T and sx
get_sz = build_function(sz_expression, T, 𝐱; expression=Val{false})

# ╔═╡ 646a9806-dc7a-4814-b1f4-877e804055ce
# returns the derivative of sz w.r.t. T
get_dsz_dT = build_function(Symbolics.derivative(sz_expression, T), 𝐱, T; expression=Val{false})

# ╔═╡ 40797264-1787-459a-b554-79bf2b9cd8a6
sz_true = 5.0

# ╔═╡ ea8bfd4d-5cbc-42b2-a027-7eb568b40391
md"### UI"

# ╔═╡ a4435d5f-b1ce-492a-8836-1d5e847719ac
function source_loc_ex_input()
    nsx = length(sxgrid)
    nT = length(Tgrid)
    return PlutoUI.combine() do Child
        p = [md"""
True x location of source = $(Child("x", Slider(range(-1.5, stop=1.5, length=100), default=1.2, show_value=true))) and known mean of 
z location = $(Child("μz", Slider(range(1, stop=5, length=100), default=5, show_value=true)))
             """,
        ]
        σ = [md"""
mean in x = $(Child("μx", Slider(sxgrid, default=sxgrid[div(nsx,2)], show_value=true)))

and standard deviation in x = $(Child("σx", Slider(range(0.05, stop=10, length=100), default=10, show_value=true)))

and standard deviation in z = $(Child("σz", Slider(range(0.001, stop=.5, length=10), default=0.05, show_value=true)))
""",]
        priorT = [
            md"""
            variance in observed traveltime = $(Child("σT", Slider(range(0.01, stop=0.5, length=100), default=0.01, show_value=true)))
            """,]

        md"""### Locate The Source!
Consider a hypothetical experiment, where the x-coordinate of the source has to be estimated given the (P or S arrival) travel time measured at a single receiver. Here, luckily, the z-coordinate of the source has a known mean, and it has a given standard deviation  

True source location is $(p) Prior information is assumed to be of Gaussian type with $(σ)

		
The uncertainty in the observed data has a standard deviation 
$(priorT)
 """
    end
end

# ╔═╡ 8ea4d826-1aa4-4756-817c-b0ddb1fc29e4
@bind sex source_loc_ex_input()

# ╔═╡ 1c78616a-d5ef-49b3-954b-7b8fcefdb1b9
# uncertainity in sz, with a fixed mean sz_true
pz = Normal(sex.μz, sex.σz)

# ╔═╡ 93d2b4ec-516d-4cea-a1c2-49ff84dfd6d0
# theoretical information 
Θex = broadcast(Iterators.product(sxgrid, Tgrid)) do (sx, T)
    pdf(pz, get_sz(T, sx)) * get_dsz_dT(sx, T)
end;

# ╔═╡ 7bb8a2e6-6e67-4fc2-84fc-15e2e5b52e8a
# prior on sx
ρsx = Normal(sex.μx, sex.σx)

# ╔═╡ 2e8c23bc-3281-414d-b1db-113e6c2ee5a9
# generate observed traveltime
Tobs = get_traveltime(sex.x, sz_true)

# ╔═╡ 48653288-ebab-47cf-93d5-c7c44615c49a
# Gaussian uncertainty for observed traveltime
ρT = Normal(Tobs, sex.σT)

# ╔═╡ 834918a4-e2ee-4445-8068-d85e04a17580
# prior information
ρex = broadcast(Iterators.product(sxgrid, Tgrid)) do (sx, T)
    pdf(ρsx, sx) * pdf(ρT, T)
end;

# ╔═╡ 94125b78-3ffc-4c9d-bd47-eba34bdee94d
plot(sxgrid, sum(ρex .* Θex, dims=2), w=2, c=:black, title="Marginalized Posterior Model Information", xlabel="Source x location", size=(500, 250), label=nothing)

# ╔═╡ 1dc77f6a-daaa-44bf-9eee-c94b807eac9f
# posterior information
σex = ρex .* Θex;

# ╔═╡ 5ad29bb0-1002-4779-ba2d-7f2debcb86fa
md"### Plots"

# ╔═╡ f8fa0021-eef0-4c13-9178-33a86fd8ba4c
begin
    @userplot pheat2d

    @recipe function f(h::pheat2d)
        grid := true
        size --> (500, 400)
        margin := 2cm
        color := :amp
        seriestype := :contourf
        colorbar := nothing
        @series begin
            h.args
        end
    end
end

# ╔═╡ 06bcd777-c41c-4a51-97c4-e1771b72a2f6
pheat2d(mrange, mrange, Z, aspect_ratio = 1, title="2-D Gaussian Density Function")

# ╔═╡ f3eb3392-787b-4cb0-b177-580e2a1456cc
scatter2D(d1, d2, lim=:auto, title="") = scatter(collect(zip(d1, d2)), xlabel="Measurement 1", ylabel="Measurement 2", label=nothing, size=(500, 500), c=:gray, title=title, lim=lim, margin=1cm, framestyle=:box)

# ╔═╡ ee622609-7a34-4943-a64d-4e632c327158
scatter2D(d1obs, d2obs, (0.5, 1.5), "Independent Measurements")

# ╔═╡ 99f08cc3-6bd9-40ee-a61e-86320a7c8a3e
scatter2D(dmix[:, 1], dmix[:, 2], :auto, "Dependent Measurements")

# ╔═╡ 7fd8f37a-a26a-42ce-a3a4-9466093052d2
Θ_plot = pheat2d(Tgrid, sxgrid, Θex, xlabel="Travel time", title="Theoretical Information: Θ(d, m)");

# ╔═╡ cef06dca-9cad-48cb-bca6-1277c8da221a
σ_plot = pheat2d(Tgrid, sxgrid, σex, xlabel="Travel time", title="Posterior Information: σ(d, m)");

# ╔═╡ 36201994-8321-4b90-a7ae-931341543b6c
ρ_plot = pheat2d(Tgrid, sxgrid, ρex, xlabel="Travel time", ylabel="Source x location", title="Prior Information: ρ(d, m)");

# ╔═╡ 9b64082d-b04e-4eff-bd0e-f5897544b2d0
plot(ρ_plot, Θ_plot, σ_plot, layout=(1, 3), size=(1000, 350), margin=4mm)

# ╔═╡ d68c2f8b-7389-4cf3-827b-2edf79ae0dc7
md"### Tikz"

# ╔═╡ 1872d922-94f0-4168-b440-cf6965200aef
tikz_default_options = raw"""
  background rectangle/.style={fill=white}, show background rectangle,
  """

# ╔═╡ 2313405a-c587-478c-bac4-de4cbfe2e99f
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

# ╔═╡ e201a21a-d684-4359-bd8b-5c94fe56cd66
md"""
## References
- Mosegaard, Klaus, and Albert Tarantola. "Probabilistic approach to inverse problems." International Geophysics Series 81.A (2002): 237-268.
- Fichtner, Andreas. "Lecture Notes on Inverse Theory." (2021).
"""

# ╔═╡ Cell order:
# ╠═02cd8bb4-d403-47f9-9de9-35d7b2b82bb8
# ╠═4c8ba5dc-a019-46ab-b7cc-b03cf4cee5f8
# ╟─91c82c81-8cfc-4b54-bed7-0592f5b39351
# ╟─8ea4d826-1aa4-4756-817c-b0ddb1fc29e4
# ╟─9b64082d-b04e-4eff-bd0e-f5897544b2d0
# ╟─94125b78-3ffc-4c9d-bd47-eba34bdee94d
# ╟─1b14a442-211e-40a6-b8d1-3854a69ef682
# ╠═760ad648-a0cf-4ec7-96bf-362346e34ba6
# ╟─c57d72de-920d-48b5-b768-661cd2435bcd
# ╠═9b3ac9a2-5bc0-43e3-a8f9-688490e19674
# ╟─f460a7f2-f812-4bdf-ba13-8a905f51ea55
# ╟─56c7e80d-068e-4c4a-a1fa-64e26daf147f
# ╟─65ed497b-18ed-4379-b220-439b51f3cb86
# ╠═58e5869f-03a7-4698-863b-c6eab02e927f
# ╠═9e6e08cb-d26f-4f17-a662-1e5097cf7625
# ╠═ee622609-7a34-4943-a64d-4e632c327158
# ╟─c3dc0533-0cb5-4c2f-a1e6-6c18dde7674a
# ╠═29bf0e42-8dd0-41e9-a646-a7714440d09b
# ╠═859865b7-a2f8-4890-a127-7b96ec0513c7
# ╠═99f08cc3-6bd9-40ee-a61e-86320a7c8a3e
# ╟─e5d5ac74-bbfd-42e1-aeaf-5fbd17d59625
# ╠═a2461a13-76f4-4a18-bcf4-fc9e7a18e5bb
# ╠═54085478-12df-4d02-b2ed-ac9a373f1b91
# ╟─4465ed36-d7f8-4407-9270-0ec20dd7e08c
# ╠═47c3f46f-158e-4eb6-aa8d-49c018168826
# ╠═88060b9a-3c70-4bed-8b0c-ce6c5b9a3061
# ╠═b3212b47-a114-4a25-88f7-5a3cd3cb6542
# ╠═06bcd777-c41c-4a51-97c4-e1771b72a2f6
# ╠═43363cb2-98d5-41d3-b276-97478829391a
# ╟─b2c008d8-48aa-4c32-b988-3552404a25b6
# ╠═ab6453f1-9746-4825-84bc-32b552538e23
# ╟─d4123859-1228-444f-ac49-ec861aae3bdb
# ╟─c7b84588-8e3a-4206-8ca5-0a7a85d395d5
# ╠═e8967bc9-daff-48b0-9846-f839848e3142
# ╠═9e0ff7a6-3984-4be0-858c-b08e95bd46b1
# ╟─7b66bc33-4d0c-4f21-9a9c-7a208708d95d
# ╟─263c76d0-051b-4349-9f2f-c6ee82ad02f6
# ╟─3add2ee1-1de6-490d-a9b7-1f01d968cf29
# ╟─70df73b9-e418-47da-81e1-735fda13f131
# ╟─a7c12637-c990-49a8-b60e-ddae71fda3e3
# ╠═a305940b-6699-4e04-a1d8-8182266b14a4
# ╠═bedc65e0-e4bb-47c6-82cd-10f1cac4d9db
# ╠═3cd564e1-cd99-4c25-8ae4-85904f997a4c
# ╟─ac9bbed8-9d8f-461e-bea5-9337cec3ec80
# ╠═da70f650-c306-461d-a607-8d2384946088
# ╠═8dca7b1e-7142-4bdc-81bf-645cc8c264c3
# ╠═6f7955f7-99fa-482c-8544-9f5bb9e75452
# ╠═7610159c-f519-4dd6-bf49-e9015dfda7d5
# ╟─83f55514-4c1e-4587-9dc9-576351dc0379
# ╟─849ecb20-7d5b-48da-9fbc-3a8bfcb4fa9b
# ╠═04d19c81-08f8-4f07-a0c7-f51c1b39d271
# ╠═90a4f7a0-7a94-4522-a4e9-99f8777e8f53
# ╠═0eac2ee9-ce0c-46e2-9796-d3bb1b0eb40f
# ╟─22c3fe50-0ace-4a6b-a4b4-7e6f63a45941
# ╠═0e768c08-c806-4de9-b459-f52065b6ce0a
# ╠═72ddf27e-de41-41f4-aeae-12f2ead3d89c
# ╠═7a0bf6ba-0353-42b5-b86b-1d43e22a22b9
# ╠═3332d7c6-8971-404b-86ae-7eaf03d0768d
# ╠═1c2abe11-ac0f-413c-93da-e521ebf5ae36
# ╠═b447013f-4e0c-44b1-ba17-167a34cd3f7b
# ╠═459e9902-05c5-4c53-ad99-7f68143e205a
# ╠═925c71dd-530f-4bd8-8685-d66747c277ec
# ╠═646a9806-dc7a-4814-b1f4-877e804055ce
# ╠═40797264-1787-459a-b554-79bf2b9cd8a6
# ╠═1c78616a-d5ef-49b3-954b-7b8fcefdb1b9
# ╠═93d2b4ec-516d-4cea-a1c2-49ff84dfd6d0
# ╠═834918a4-e2ee-4445-8068-d85e04a17580
# ╠═1dc77f6a-daaa-44bf-9eee-c94b807eac9f
# ╠═7bb8a2e6-6e67-4fc2-84fc-15e2e5b52e8a
# ╠═2e8c23bc-3281-414d-b1db-113e6c2ee5a9
# ╠═48653288-ebab-47cf-93d5-c7c44615c49a
# ╟─ea8bfd4d-5cbc-42b2-a027-7eb568b40391
# ╠═a4435d5f-b1ce-492a-8836-1d5e847719ac
# ╟─5ad29bb0-1002-4779-ba2d-7f2debcb86fa
# ╠═f8fa0021-eef0-4c13-9178-33a86fd8ba4c
# ╠═f3eb3392-787b-4cb0-b177-580e2a1456cc
# ╠═7fd8f37a-a26a-42ce-a3a4-9466093052d2
# ╠═cef06dca-9cad-48cb-bca6-1277c8da221a
# ╠═36201994-8321-4b90-a7ae-931341543b6c
# ╟─d68c2f8b-7389-4cf3-827b-2edf79ae0dc7
# ╟─1872d922-94f0-4168-b440-cf6965200aef
# ╟─2313405a-c587-478c-bac4-de4cbfe2e99f
# ╟─e201a21a-d684-4359-bd8b-5c94fe56cd66
