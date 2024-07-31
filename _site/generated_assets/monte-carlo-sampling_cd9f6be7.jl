### A Pluto.jl notebook ###
# v0.19.40

#> [frontmatter]
#> chapter = "1"
#> title = "Monte Carlo Methods"
#> layout = "layout.jlhtml"
#> description = "How can we balance the trade-off between taking large steps with many rejections versus taking small steps with few rejections?"

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

# ╔═╡ e796296e-9941-11ec-34fa-cbc6c2c3eaa0
begin
    using Distributions, Plots, PlutoUI, PlutoTeachingTools
    using StatsPlots, MLUtils, Interpolations, Measures, ImageFiltering
end

# ╔═╡ 28b5ab19-2bce-4f37-985a-af550be109ef
TableOfContents()

# ╔═╡ 6bc469ef-3543-4fe7-8158-573981fe29f0
ChooseDisplayMode()

# ╔═╡ aa698b0c-aa70-477b-b6df-71c1e8aba009
md"# Monte Carlo Sampling"

# ╔═╡ 8b8e36b7-f10c-4d10-ac2e-711acad3d745
@bind sample_true_density Button("Change True Density Model")

# ╔═╡ 3ba056a9-86d4-45c6-a50b-1bf2d363205c
md"""#### Prior Movie
$(@bind sample_prior_density Button("Sample Prior Density Models"))
Smooth? $(@bind smooth_prior_flag CheckBox(default=false))
"""

# ╔═╡ b1928ff1-b74e-4c05-8778-c9c57fc06b2a
md"""#### Posterior Movie
$(@bind sample_posterior_density Button("Sample Posterior Density Models"))
Smooth? $(@bind smooth_posterior_flag CheckBox(default=false))
"""

# ╔═╡ 76558020-50a2-4de7-8ce3-823aaa35fde9
md"### Simpler Demo"

# ╔═╡ b61b56bc-1ed7-4603-9369-737f5f5c5c73
md"## Striking the balance"

# ╔═╡ 0eb49529-fd09-4faa-8577-0841db631f1b
Markdown.MD(Markdown.Admonition("wrg", "(Large Steps & Many Rejects) Vs. (Small Steps & Few Rejects)", [md"""
The acceptance rate of the Metropolis criterion has to be around $30-50\%$. If the acceptance_ex1 rate is, we are not moving fast enough in the model space. On the other hand, if the acceptance_ex1 rate is much smaller, then we are wasting computer resources as most of the samples are not accepted.
"""]))

# ╔═╡ d534b263-5ec5-41b3-af9f-fe3f0cb4b7c6
md"## Metropolis Criterion"

# ╔═╡ ab5a5cc2-1844-4e0a-b721-707d00464c72
Markdown.MD(Markdown.Admonition("info", "Metropolis Criterion", [md"""
Consider some random rules that define a random walk. This random walk, if unmodified, samples some initial probability distribution.
The metropolis algorithm modifies the random walk, by accepting or rejecting a proposed transition $m_i\to m_j$ using the following rules:
- if $L(m_j) ≥ L(m_i)$, then accept the proposed transition;
- else if $L(m_j) < L(m_i)$, then decide randomly whether to accept or reject, with the following probability of accepting the move to $m_j$
```math
P_{i\to j}=\frac{L(m_j)}{L(m_i)}.

```
It can be shown that the random walker samples the conjunction of the probability densities $\rho(m)$ and 
$L(m)$, given by [^1]
```math
σ(m) = \rho(m)\wedge\,L(m)=\mathrm{const}. \rho(m)L(m)
```

[^1]: We have assumed the homogeneous distribution $\mu(m)$ to be a constant for simplicity.
"""
]))

# ╔═╡ 6e6eb0a2-5e21-46cf-b768-ef0827728022
md"""
This function takes three mandatory inputs: `m0` which is the initial model, prior which is a function that samples from the prior distribution, and `L` which is a function that evaluates the likelihood of a model. The function also takes an optional input `N` which is the number of samples to generate and is set to 100 by default.

The function generates samples using the Metropolis-Hastings algorithm, which is a Markov Chain Monte Carlo (MCMC) method. It initializes the first sample to `m0`, and then repeatedly generates a new candidate sample using the prior function. It then evaluates the likelihood of the candidate sample using the `L` function and computes the Metropolis-Hastings ratio. If the ratio is greater than a randomly generated number from a Bernoulli distribution, the candidate sample is accepted as the next sample, and the process is repeated until N samples are generated or the loop breaks.

The function returns two objects: `msamples_ex1`, which is an array of the generated samples, and `acceptance_ex1_history`, which is an array that keeps track of whether a sample was accepted or rejected at each step of the Metropolis-Hastings algorithm.
"""

# ╔═╡ 0b20290a-200b-44aa-8b1e-65a63044be92
function generate_samples(m0, walk, L; N=10)

    msamples_ex1 = []
    push!(msamples_ex1, m0)
    acceptance_ex1_history = []
    metropolis_ratios = []

    # create an array to store previous model
    mp = copy(m0)

    # evaluate the initial likelihood
    p_current = (L === nothing) ? 1.0 : L(m0)

    i = 0
    while i ≤ N

        # sample from prior to get destination
        m = walk(mp)

        # evaluate prob of destination
        p_prop = (L === nothing) ? 1.0 : L(m)

        # metropolis ratio
        metropolis_ratio = p_prop / p_current
        # push!(metropolis_ratios, metropolis_ratio)

        if (rand(Bernoulli(min(1, metropolis_ratio))))
            mp = copy(m)
            p_current = p_prop
            push!(msamples_ex1, m)
            push!(acceptance_ex1_history, true)
            i += 1
        else
            # simply revert back
            push!(acceptance_ex1_history, false)
        end

    end
    return msamples_ex1, acceptance_ex1_history

end

# ╔═╡ ae40fbdc-614c-432c-8503-ca670bdebe8b
md"## Simple 2D Example"

# ╔═╡ 750fdbb0-695f-486c-a4e4-351e321ff853
# method to generate 2-D example target probability density
function Lex1(m)
    s1 = 1.0
    s2 = 0.5
    return exp(-0.5 * ((m[1] + 1.7)^2 + m[2]^2) / s1^2) / (
        2.0 * pi * s2
    ) + exp(-0.5 * ((m[1] - 2.0)^2 + m[2]^2) / s2^2) / (2.0 * pi * s2)
end

# ╔═╡ cd5826ab-8c86-4a8a-9eff-f6fbd241bdce
md"""
###### Use Metropolis criterion? $(@bind metro Select([Lex1=>"Yes", nothing=>"No"]))
"""

# ╔═╡ c75c0510-fbca-48db-94b9-2c1f62096a2d
range_ex1 = range(-10, stop=10, length=256)

# ╔═╡ 02ef9d57-a983-4756-9171-494f622f3db9
Lex1_grid = [Lex1([x, y]) for x in range_ex1, y in range_ex1];

# ╔═╡ 34c3d53a-4b2e-4cf4-b9d1-a6e389fca682
md"""
This function generates a 2D random walk step by drawing two independent random numbers from a normal distribution with mean 0 and standard deviation 1, scaling them by a factor of 10, and returning them as a 2D vector. The resulting step is uncorrelated with any previous steps and is independent of the current position. This is an example of an independent random walk.
"""

# ╔═╡ e0a88b43-2dab-411e-a45d-d758c6905359
md"""
This function generates a proposal for a 2D random walk by adding random perturbations to the current model m based on a Gaussian distribution with standard deviation sigma of 1.0. This is a local proposal, where each step is a small perturbation around the current point.
"""

# ╔═╡ 820f3079-f528-4496-b384-4cf3624468d3
heatmap(range_ex1, range_ex1, Lex1_grid, c=:grays, title="Likelihood Distribution", aspect_ratio=1)

# ╔═╡ 48a5055f-a239-4d8e-9369-74e1fc9a9bd9
md"## Gravity Example[^mosegaard1995]"

# ╔═╡ 55f94d7b-b083-49ed-8ad9-57083720dd45
md"""
Note that the complex prior knowledge used in this example renders the posterior distribution to be non-Gaussian.
"""

# ╔═╡ 9ee18955-a949-4e85-8157-5832c423f3e3
xrange_ex2 = range(1, stop=40, length=20)

# ╔═╡ 8e96a8d0-aca8-4018-8b3d-0706d9aa1cec
ρrange = (100, 6000)

# ╔═╡ 634db275-b82b-4439-9f1e-28f89bcfe2c1
md"### Samplers"

# ╔═╡ e023d351-5560-4f08-8f5b-5252844b3e0a
md"""
This sampler output neighbors of an input model
by creating or destroying a new interface
in the model.
"""

# ╔═╡ 3b31458a-3c5d-46e1-a5e8-c08cb334e4cb
md"""
### Physics
We measure the horizontal derivative of the vertical component of the acceleration due to gravity.
Contribution from a homogeneous half layer is given by
```math
G\Delta\rho\,\log\left(\frac{z_2^2+x^2}{z_1^2+x^2}\right),
```
where
- depth of the top pf the layer is $z_1$;
- the depth of the bottom is  $z_2$; 
- horizontal density contrast is $\Delta\rho$.
"""

# ╔═╡ 7997ef7f-47ba-419a-98af-0a5cdddb2915
function ggravity(m, xrange=xrange_ex2)
    Gconst = 6.6743 * 1e-11 # m3 kg-1 s-2
    ρ_left_half = 2570 # kg/m3
    l, ρ = splitobs(m, at=0.5)

    l = l .* 1e3
    return broadcast(xrange) do x
        xm = x * 1e3 # km to m
        data = 0.0
        d = 0 # first layer at zero depth
        for i in 1:length(l)
            data += Gconst * (ρ[i] - ρ_left_half) * log((abs2(d + l[i]) + abs2(xm)) * inv(abs2(d) + abs2(xm)))
            d += l[i] # depth of top of next layer
        end
        return data
    end
end

# ╔═╡ 02714b1b-7f31-405d-b47c-3a12b688db4b
md"### Prior"

# ╔═╡ 9d8b772d-a5c7-421b-b9aa-f851fd26896c
pρ = Gumbel(3000, 400)

# ╔═╡ 90a0357c-e86e-4c63-a419-d83231654f96
plot(pρ, w=2, label=nothing, title="Prior Mass Density Distribution", c=:black, size=(500, 250))

# ╔═╡ d707553d-c917-42cc-b47f-eb6043aee0c3
pl = Exponential(4)

# ╔═╡ fde40e79-70c4-464d-9f2c-b6eeb1b3446d
function randomwalk_density_independent(m, pl=pl, pρ=pρ)

    l = Vector{Float64}(undef, 0)
    d = 0 # cumulative sum of widths
    while d <= 100 # less than 100 km
        l1 = rand(pl)
        push!(l, l1)
        d += l1
    end
    # change last layer to obey 100 km constraint
    l[end] = 100 - sum(l[1:end-1])
    # pick density values
    ρ = rand(pρ, length(l))
    # model vector is vcat of l and d
    return cat(l, ρ, dims=1)
end

# ╔═╡ 5af96d1b-e9c0-4cca-86b1-6f40c7497a73
begin
    sample_true_density
    mex2_true = randomwalk_density_independent(0)
end

# ╔═╡ 8dad04b1-0f7f-4e52-87b4-caaf0f120bdb
dobs_ex2 = ggravity(mex2_true, xrange_ex2)

# ╔═╡ 28a067e1-f814-4e77-bb8c-9285b6e7cb11
msamples_ex2_prior = [randomwalk_density_independent(0) for i in 1:100]

# ╔═╡ 460442fe-dc55-4f88-895b-3e1bddabdc9b
function randomwalk_density_dependent(m, pl=pl, pρ=pρ)
    mout = copy(m)
    l, ρ = splitobs(mout, at=0.5)
    nlayers = length(l)
    if (rand(Bernoulli(0.5)))
        # change density of randomly chosen layer
        k = randobs(1:nlayers)
        ρ[k] = rand(pρ)
    else
        for wait in 1:5
            ib = randobs(2:nlayers)
            L = l[ib] + l[ib-1]
            l[ib] = rand(truncated(pl, 0, L))
            l[ib-1] = L - l[ib]
            ρ[ib] = rand(pρ)
            ρ[ib-1] = rand(pρ)
        end
    end
    return cat(l, ρ, dims=1)
end

# ╔═╡ 894be65f-61d9-4f9e-9cac-361ee4bf3e72
plot(pl, w=2, label=nothing, title="Prior Layer Thickness Distribution", c=:black, size=(500, 250))

# ╔═╡ bec76e39-df37-4aa4-89eb-c59abc52d3bb
md"### Likelihood"

# ╔═╡ eac5e9ca-1146-44bc-bd0d-39a11df3a215
md"""## References
[^mosegaard1995]: Monte Carlo sampling of solutions to inverse problems, Mosegaard and Tarantola, JGR, 1995.
"""

# ╔═╡ 00a66835-c1f0-40cd-b1c2-d0e12f9396c8
md"## Appendix"

# ╔═╡ 83b6cba5-a46d-4a88-809f-fdd36841ad84
md"### UI"

# ╔═╡ cea5e440-2194-4cc2-89e3-fa459d48509c
function sampler_input_ex2()

    return PlutoUI.combine() do Child
        p = [md"""
Initial model m₁ = $(Child("m₁0", Slider(range_ex1, default=-5, show_value=true))) and m₂ =
 $(Child("m₂0", Slider(range_ex1, default=-5, show_value=true)))
             """,
        ]
        σ = [md"""

Standard deviation for prior distribution for drawing independent samples = $(Child("σglobal", Slider(range(0.05, stop=10, length=100), default=10, show_value=true)))

Standard deviation for local proposal = $(Child("σlocal", Slider(range(0.001, stop=10, length=10), default=1, show_value=true)))
""",]


        md"""###### Design the random walk
$(p) 
$(σ)
"""
    end
end

# ╔═╡ 93cadb8e-a158-4553-b546-0fa384a6df01
function sampler_input_ex1()

    return PlutoUI.combine() do Child
        p = [md"""
Initial model m₁ = $(Child("m₁0", Slider(range_ex1, default=-5, show_value=true))) and m₂ =
 $(Child("m₂0", Slider(range_ex1, default=-5, show_value=true)))
             """,
        ]
        σ = [md"""

Standard deviation for prior distribution for drawing independent samples = $(Child("σglobal", Slider(range(0.05, stop=10, length=100), default=10, show_value=true)))

Standard deviation for local proposal = $(Child("σlocal", Slider(range(0.001, stop=10, length=10), default=1, show_value=true)))
""",]


        md"""###### Design the random walk
$(p) 
$(σ)
"""
    end
end

# ╔═╡ d2427a6c-1014-4fea-a185-543a5f2e7b02
@bind ex1 confirm(sampler_input_ex1())

# ╔═╡ f9c4bfe0-6fd6-4e72-a3f1-76757957373f
function randomwalk_ex1_independent(m)
    sigma = ex1[:σglobal]
    return sigma .* (randn(2) .- 0.5)
end

# ╔═╡ 6d2b530a-0376-47a0-9e15-c22bae6cac32
function randomwalk_ex1_local(m)
    sigma = ex1[:σlocal]
    return m .+ sigma .* randn(2)
end

# ╔═╡ aa3e9da8-8053-47e8-b859-519bb87be58a
md"""
Generate $(@bind Nex1 Slider(1:256, default=20, show_value=true)) samples which are 
$(@bind walk Select([randomwalk_ex1_independent=>"Independent", randomwalk_ex1_local=>"Markov Chain"]))
"""

# ╔═╡ 9d52551f-1384-487b-9a81-ddd42a98639f
msamples_ex1, acceptance_ex1 = generate_samples([ex1[:m₁0], ex1[:m₂0]], walk, metro, N=Nex1)

# ╔═╡ aa1c78b9-2cdb-49c8-8e5f-a34ff1a9c059
begin
    p1 = heatmap(range_ex1, range_ex1, Lex1_grid, c=:grays, lim=(-10, 10), xlabel="m₁", aspect_ratio=1, ylabel="m₂", title="Accepted $(Nex1)/$(length(acceptance_ex1)) with ratio $(round(count(acceptance_ex1)/length(acceptance_ex1), digits=3))")
    scatter!(p1, [msamples_ex1[1][1:1]], [msamples_ex1[1][2:2]], c=:yellow, m=:x, label="Initial Sample")
    @gif for isamp in 2:length(msamples_ex1)
        # println(m, acc)
        plot!(p1, [msamples_ex1[isamp-1][2], msamples_ex1[isamp][2]], [msamples_ex1[isamp-1][1], msamples_ex1[isamp][1]], c=:blue, label=nothing)
    end

end

# ╔═╡ b306ce65-792e-46e3-b522-a8e77d759a1d
md"### Plots"

# ╔═╡ ab042213-8885-4715-8ddd-ba32b8289f24
begin
    @userplot pdensity_model

    @recipe function f(h::pdensity_model; smooth=false)
        grid := true
        size --> (200, 300)
        margin --> 1cm
        color --> :black
        seriestype := :steppre
        colorbar := nothing
        yflip := true
        xlim --> (100, 6000)
        ylim --> (0, 100)
        w --> 2
        xlabel --> "Density (kg/m3)"
        ylabel --> "Depth (km)"
        xticks --> [1000, 5000]
        titlefontsize --> 10
        guidefontsize --> 7
        tickfontsize --> 7
        label --> nothing
        @series begin
            l, ρ = splitobs(h.args[1], at=0.5)
            l1 = vcat([0.0], cumsum(l))
            ρ1 = vcat(ρ, ρ[end:end])
            if h.args[2]
                itp = interpolate((l1,), ρ1, Gridded(Linear()))
                zrange = range(0.1, stop=99.9, length=256)
                ρ2, l2 = vec(imfilter(cat(itp.(zrange), dims=2), Kernel.gaussian(3))), zrange
            else
                ρ2, l2 = ρ1, l1
            end
            ρ2, l2
        end
    end
end

# ╔═╡ 343699e4-bc2f-41e0-99b1-c8dc3f15cc40
plot(pdensity_model(msamples_ex2_prior[1], false, title="Initial Sample"), pdensity_model(randomwalk_density_dependent(msamples_ex2_prior[1]), false, title="Perturbed Sample"), size=(500, 300))

# ╔═╡ 5cb023fd-1835-45fb-bc26-4dc34fc5216a
function plot_true_density_model()
    plot(pdensity_model(mex2_true, false, title="True Earth Model"), pdensity_model(mex2_true, true, title="Smoothed True Model"), size=(450, 300), margin=5mm)
end

# ╔═╡ 71f399d5-4d38-4b9b-9a3e-c9ca6c1f5cfe
TwoColumnWideRight(md"""
### Gravity Inversion
$(PlutoUI.LocalResource("gravity_inversion_example.png", :width=>300))
Standard deviation of normally distributed gravity measurements = 
$(@bind σex2 Slider(range(0.1e-8, stop=3e-8, length=10), default=0.1e-6, show_value=true))
""",
    md"$(plot_true_density_model())"
)

# ╔═╡ 0acabf48-b136-4dbb-9393-a7260e01c75d
pd_ex2 = MvNormal(dobs_ex2, σex2)

# ╔═╡ 69f0102a-0107-4c94-9423-6b5a63f88c2e
function Lex2(m, pd=pd_ex2, xrange=xrange_ex2)
    d = ggravity(m, xrange)
    return pdf(pd_ex2, d)
end

# ╔═╡ 6f9a7f23-91d1-429e-a532-3f02bc8ea431
msamples_ex2, acceptance_ex2 = generate_samples(randomwalk_density_independent(0), randomwalk_density_dependent, Lex2, N=100)

# ╔═╡ 077cbd07-1b6f-41ac-9c50-23132484f5b1
count(acceptance_ex2) / length(acceptance_ex2)

# ╔═╡ 4ba948e5-5cbd-42f7-a31a-9d467485d4a8
function plot_density_model(m)
    p1 = plot(playered_model(m, false))
    p2 = plot(playered_model(m, true))

    plot(p1, p2, size=(400, 300), margin=1cm)
end

# ╔═╡ 95492189-d9a7-4531-8f97-36148e26e438
function pdensity_models(mvec, smooth)
    plot([pdensity_model(m, smooth, xlabel="", ylabel="", margin=1mm, left_margin=-2.0 * Plots.mm, right_margin=-2.0 * Plots.mm, showaxis=false, axis=nothing) for m in mvec]..., layout=(1, length(mvec)), size=(650, 300), axis=nothing)
end

# ╔═╡ f0742a9f-3f65-440b-8124-fd087e14d214
begin
    sample_prior_density
    pdensity_models([randomwalk_density_independent(0) for i in 1:6], smooth_prior_flag)
end

# ╔═╡ f8343641-b00c-4116-8aad-ab94af42b0c0
begin
    sample_posterior_density
    pdensity_models(randobs(msamples_ex2[end-20:end], 6), smooth_posterior_flag)
end

# ╔═╡ 5af97c29-6768-43f6-8adc-9e8e88d1f966
pdata_ex2 = plot(xrange_ex2, dobs_ex2, w=2, size=(500, 250), xlim=extrema(xrange_ex2), label="Observed Data", xlabel="sensor horizontal position", ylabel="", titlefontsize=10, guidefontsize=7, tickfontsize=7,);

# ╔═╡ 8c0fc32a-415a-4912-8453-c156ed8d8ba3
function plot_prior_data_ex2()
    pdata_ex2_prior = deepcopy(pdata_ex2)
    for i in 1:20
        plot!(pdata_ex2_prior, xrange_ex2, ggravity(randomwalk_density_independent(0), xrange_ex2), w=2, c=:black, alpha=0.3, label=nothing, title="Prior")
    end
    pdata_ex2_prior
end

# ╔═╡ 32ebb870-4a43-432c-83c8-24c935fbecd6
function plot_post_data_ex2()
    pdata_ex2_post = deepcopy(pdata_ex2)
    for m in msamples_ex2[end-20:end]
        plot!(pdata_ex2_post, xrange_ex2, ggravity(m, xrange_ex2), w=2, c=:black, alpha=0.3, label=nothing, title="Posterior")
    end
    pdata_ex2_post
end

# ╔═╡ 49a1286f-218f-4070-b708-b14c60238839
md"""
#### Data Movie
$(plot(plot_prior_data_ex2(), plot_post_data_ex2(), size=(650, 250), link=:y))
"""

# ╔═╡ Cell order:
# ╠═28b5ab19-2bce-4f37-985a-af550be109ef
# ╠═6bc469ef-3543-4fe7-8158-573981fe29f0
# ╟─aa698b0c-aa70-477b-b6df-71c1e8aba009
# ╟─71f399d5-4d38-4b9b-9a3e-c9ca6c1f5cfe
# ╟─8b8e36b7-f10c-4d10-ac2e-711acad3d745
# ╟─3ba056a9-86d4-45c6-a50b-1bf2d363205c
# ╟─f0742a9f-3f65-440b-8124-fd087e14d214
# ╟─49a1286f-218f-4070-b708-b14c60238839
# ╟─b1928ff1-b74e-4c05-8778-c9c57fc06b2a
# ╟─f8343641-b00c-4116-8aad-ab94af42b0c0
# ╟─76558020-50a2-4de7-8ce3-823aaa35fde9
# ╟─d2427a6c-1014-4fea-a185-543a5f2e7b02
# ╟─aa3e9da8-8053-47e8-b859-519bb87be58a
# ╟─cd5826ab-8c86-4a8a-9eff-f6fbd241bdce
# ╟─aa1c78b9-2cdb-49c8-8e5f-a34ff1a9c059
# ╟─b61b56bc-1ed7-4603-9369-737f5f5c5c73
# ╟─0eb49529-fd09-4faa-8577-0841db631f1b
# ╟─d534b263-5ec5-41b3-af9f-fe3f0cb4b7c6
# ╟─ab5a5cc2-1844-4e0a-b721-707d00464c72
# ╟─6e6eb0a2-5e21-46cf-b768-ef0827728022
# ╠═0b20290a-200b-44aa-8b1e-65a63044be92
# ╟─ae40fbdc-614c-432c-8503-ca670bdebe8b
# ╠═750fdbb0-695f-486c-a4e4-351e321ff853
# ╠═c75c0510-fbca-48db-94b9-2c1f62096a2d
# ╠═02ef9d57-a983-4756-9171-494f622f3db9
# ╟─34c3d53a-4b2e-4cf4-b9d1-a6e389fca682
# ╠═f9c4bfe0-6fd6-4e72-a3f1-76757957373f
# ╟─e0a88b43-2dab-411e-a45d-d758c6905359
# ╠═6d2b530a-0376-47a0-9e15-c22bae6cac32
# ╟─820f3079-f528-4496-b384-4cf3624468d3
# ╠═9d52551f-1384-487b-9a81-ddd42a98639f
# ╟─48a5055f-a239-4d8e-9369-74e1fc9a9bd9
# ╟─55f94d7b-b083-49ed-8ad9-57083720dd45
# ╠═9ee18955-a949-4e85-8157-5832c423f3e3
# ╠═8e96a8d0-aca8-4018-8b3d-0706d9aa1cec
# ╟─634db275-b82b-4439-9f1e-28f89bcfe2c1
# ╠═fde40e79-70c4-464d-9f2c-b6eeb1b3446d
# ╟─e023d351-5560-4f08-8f5b-5252844b3e0a
# ╠═460442fe-dc55-4f88-895b-3e1bddabdc9b
# ╠═343699e4-bc2f-41e0-99b1-c8dc3f15cc40
# ╟─3b31458a-3c5d-46e1-a5e8-c08cb334e4cb
# ╠═7997ef7f-47ba-419a-98af-0a5cdddb2915
# ╠═5af96d1b-e9c0-4cca-86b1-6f40c7497a73
# ╠═8dad04b1-0f7f-4e52-87b4-caaf0f120bdb
# ╟─02714b1b-7f31-405d-b47c-3a12b688db4b
# ╠═28a067e1-f814-4e77-bb8c-9285b6e7cb11
# ╠═0acabf48-b136-4dbb-9393-a7260e01c75d
# ╠═9d8b772d-a5c7-421b-b9aa-f851fd26896c
# ╠═90a0357c-e86e-4c63-a419-d83231654f96
# ╠═d707553d-c917-42cc-b47f-eb6043aee0c3
# ╠═894be65f-61d9-4f9e-9cac-361ee4bf3e72
# ╟─bec76e39-df37-4aa4-89eb-c59abc52d3bb
# ╠═69f0102a-0107-4c94-9423-6b5a63f88c2e
# ╠═6f9a7f23-91d1-429e-a532-3f02bc8ea431
# ╠═077cbd07-1b6f-41ac-9c50-23132484f5b1
# ╟─eac5e9ca-1146-44bc-bd0d-39a11df3a215
# ╟─00a66835-c1f0-40cd-b1c2-d0e12f9396c8
# ╠═e796296e-9941-11ec-34fa-cbc6c2c3eaa0
# ╟─83b6cba5-a46d-4a88-809f-fdd36841ad84
# ╠═cea5e440-2194-4cc2-89e3-fa459d48509c
# ╠═93cadb8e-a158-4553-b546-0fa384a6df01
# ╟─b306ce65-792e-46e3-b522-a8e77d759a1d
# ╠═5cb023fd-1835-45fb-bc26-4dc34fc5216a
# ╠═ab042213-8885-4715-8ddd-ba32b8289f24
# ╠═4ba948e5-5cbd-42f7-a31a-9d467485d4a8
# ╠═95492189-d9a7-4531-8f97-36148e26e438
# ╠═5af97c29-6768-43f6-8adc-9e8e88d1f966
# ╠═8c0fc32a-415a-4912-8453-c156ed8d8ba3
# ╠═32ebb870-4a43-432c-83c8-24c935fbecd6
